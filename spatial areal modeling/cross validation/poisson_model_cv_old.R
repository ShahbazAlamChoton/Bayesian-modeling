poisson_spatial_model <- function (covariate_data, count_var, test_locator, 
                                   # X_train, X_test, 
                                   betas, theta_init, sigma2_init, tau2_init, 
                                   a_tau_pr, b_tau_pr, a_sigma_pr, b_sigma_pr, num_iter, burn_in, thin_in) {
  set.seed(729)
  #-----------------------
  # data setup
  #-----------------------
  # Train test setup
  X <- cbind(X_intercept = rep(1, nrow(covariate_data)), covariate_data); n <- nrow(covariate_data); p <- ncol(X) - 1
  train_locator <-  seq_along(1:n)[-test_locator]
  X_train <- X[train_locator, ]; n_train <- nrow(X_train); tX_train = t(X_train)  
  X_test <- X[test_locator, ]; n_test <- nrow(X_test)
  
  ####### prior info:
  # prior information for all beta, rest of the parameters prior coming from function argument
  # Beta prior setup
  mu_beta = matrix(c(2, 2, 1, -1, 0, 2), ncol = 1);  sigma_beta = 10^4 * diag(p + 1)
  
  ###### Assigning initial values 
  beta = matrix(betas, ncol = 1); theta = theta_init; tau2 = tau2_init; sigma2 = sigma2_init
  rm(betas)
  # rm(betas, theta_init, tau2_init, sigma2_init)
  # iteration has no effect on estimates of tau2, A_sigma2 and sigma_beta_inv, 
  # so calculating outside of the loop
  A_tau2 = ((n - eigen0_count)/2) + a_tau_pr          # posterior scale parameter for tau2 
  A_sigma2 = (n_train/2) + a_sigma_pr                       # posterior scale parameter for sigma2
  sigma_beta_inv = chol2inv(chol(sigma_beta))
  
  # Define overall Response (train + test) vector using train_response + initial test response
  log_lambda = matrix(log_lamda_train, ncol = 1)
  current_lambda = exp(log_lambda)
  #output matrix
  ncol_tempM = length(c(1, 1, beta, sigma2, tau2, theta, log_lambda, log_lambda)) 
  tempM <- matrix(0, nrow = num_iter, ncol = ncol_tempM)
  
  pb <- txtProgressBar(min = 1, max = num_iter, style = 3)
  for (j in 1:num_iter) {
    # tic = Sys.time()
    
    ###### update log_lambda: Metropolis-Hastings
    # propose a log_lambda value from Normal
    proposed_log_lambda <- rnorm(n_train, mean = log_lambda, sd = .02) # sd can be changed
    proposed_lambda = exp(proposed_log_lambda)
    mean_log_lambda = as.vector(X_train %*% beta + theta[train_locator])
    # current_lambda <- exp(log_lambda)
    # ratio
    l_posterior_proposal =  -(2 * sigma2)^-1 * (proposed_log_lambda - mean_log_lambda)^2 - proposed_lambda + count_var[train_locator] * log(proposed_lambda)
    l_posterior_current = - (2 * sigma2)^-1 * (log_lambda - mean_log_lambda)^2 - current_lambda + count_var[train_locator] * log(current_lambda)
    l_target_by_proposal = l_posterior_proposal - l_posterior_current # because of taking log
    
    # accepting proposal
    l_uniform_vector = log(runif(n_train, min = 0, max = 1))
    gain <- l_uniform_vector < l_target_by_proposal
    log_lambda[gain] <- proposed_log_lambda[gain] # updated log_lambda 
    current_lambda[gain] <- proposed_lambda[gain] # required for next iteration
    lambda_acc <- as.numeric(gain) # vector of log_lambda acceptance (binary)
    acc_rate_per_iteration <- mean(lambda_acc)  # lambda_acc is binary
    
    
    ###### Update theta: 
    res_nonspatial <-  log_lambda - X_train %*% beta # (log_lambda- Xb) part in theta posterior
    sum_sigma_w <- sigma2 * D_w_vector
    theta_post_var = (sigma2 * tau2)/(tau2 + sum_sigma_w)
    
    # theta= matrix(0, nrow = n_train + n_test, ncol = 1)
    
    for(i in 1: n_train) {  # update depends on prior and likelihood
      theta_position <- train_locator[i]
      theta_post_mean <-  (tau2 * res_nonspatial[i] + sigma2 * sum(theta[W_list[[theta_position]]]))/(tau2 + sum_sigma_w[theta_position])
      theta_post <- rnorm(1, theta_post_mean, sqrt(theta_post_var[theta_position]))
      theta[theta_position] <- theta_post
      
      rm(theta_position, theta_post_mean, theta_post)
    } 
    
    for(i in 1:n_test) { # update depends on prior only
      theta_position <- test_locator[i]
      theta_pred_mean <- sum(theta[W_list[[theta_position]]])/D_w_vector[theta_position]
      theta_pred_var <- tau2/D_w_vector[theta_position]
      theta_post_pred <- rnorm(1, theta_pred_mean, sqrt(theta_pred_var))
      theta[theta_position] <- theta_post_pred
      
      rm(theta_position, theta_pred_mean, theta_pred_var, theta_post_pred)
    }
    theta[island_locator] <- theta[island_locator] - mean(theta[island_locator])
    
    
    ###### Update tau2:
    # what we get from full conditional of tau2
    # prior on tau2 ~ IG( 2.5, 2.5)
    # posterior shape of tau2 is fixed and calculated before iteration loop
    # A_tau2 = ((n - eigen0_count)/2) + a_tau_pr
    # B_tau2 = (t(theta) %*% (D_w_tt - W_tt) %*% theta)/2 + b_tau_pr
    B_tau2 = .5 * (sum(D_w_vector * theta^2) - sum(theta * sapply(W_list, function(x) sum(theta[x])))) + b_tau_pr
    tau2 <- 1/rgamma(1, shape = A_tau2, rate = B_tau2)
    
    
    ###### Update sigma2:
    # what we get from full conditional of sigma2
    # prior on sigma2 ~ IG( 2.5, 2.5)
    # A_sigma2 = (n/2) + a_sigma_pr  # calculated before loop
    ss <- sum((res_nonspatial - theta[train_locator])^2) # res_nonspatial = log_lambda - X_train*beta
    B_sigma2 = ss/2 + b_sigma_pr
    sigma2 <- 1/rgamma(1, shape = A_sigma2, rate = B_sigma2)
    
    ###### Update beta:
    # what we get from full conditional of beta
    # sigma_beta_inv = chol2inv(chol(sigma_beta)) calculated before loop
    # x_transpose_sigma_inv = tX_train %*% ((sigma2^-1) * diag(n))
    x_transpose_sigma_inv = (sigma2^-1) * tX_train
    A_beta = x_transpose_sigma_inv %*% X_train + sigma_beta_inv; A_beta_inv = chol2inv(chol(A_beta))
    B_beta = x_transpose_sigma_inv %*% (log_lambda - theta[train_locator]) + sigma_beta_inv %*% mu_beta
    beta <- mvrnorm(1, A_beta_inv %*% B_beta, A_beta_inv)
    
    
    ####### Log-likelihood for all and training data
    fitted.mean_train = exp(log_lambda)
    # fitted.mean_train <- exp(X_train %*% beta + theta[train_locator])
    lg_lkhd_train <- sum(dpois(count_var[train_locator], lambda = fitted.mean_train, log = T))
    
    # Accumulates all parameters for a single iteration
    tempM[j, ] <-  c(lg_lkhd_train, acc_rate_per_iteration, beta, sigma2, tau2, theta, log_lambda, lambda_acc)
    
    rm(B_tau2, B_sigma2, x_transpose_sigma_inv,
       A_beta, B_beta, A_beta_inv, 
       acc_rate_per_iteration, lambda_acc, proposed_log_lambda, fitted.mean_train, lg_lkhd_train)
    
    setTxtProgressBar(pb, j) # to see the progress 
  } 
  
  output <- as.data.frame(tempM) %>% 
    .[(burn_in+1):num_iter, ] %>%
    .[seq(from = 1, to = num_iter - burn_in, by = thin_in), ] %>%
    set_names(c("lg_lkhd_train", "acc_rate_iter",  paste0("beta", 0:p),  "sigma2", "tau2", paste0("theta", 1:n), 
                paste0("log_lambda", 1:n_train), paste0("lambda", 1:n_train, "_acc")))
  
  ####### Count Prediction
  matrix(data = rnorm(6, mean = 0, sd = c(1, 100)), nrow = 3, ncol = 2, byrow = T)
  
  epsilon_test <- matrix(data = rnorm(n_test * nrow(output), mean = 0, sd = sqrt(output$sigma2)), nrow = n_test, ncol = nrow(output), byrow = T)
  lambda_test <- exp(X_test %*% t(output[, paste0("beta", 0:p)]) + t(output[, paste0("theta", test_locator)]) + epsilon_test)
  predicted_count <- apply(lambda_test, c(1, 2), function(x) rpois(1, x))
  
  prediction_summary <- data.frame(actual_count = count_var[test_locator],
                                   median_prediction = apply(predicted_count, 1, median),
                                   mean_prediction = apply(predicted_count, 1, mean),
                                   q_05 = apply(predicted_count, 1, quantile, .05),
                                   q_95 = apply(predicted_count, 1, quantile, .95)) %>% 
    mutate(include = if_else(actual_count >= q_05 & actual_count <= q_95, "Yes", "No"),
           abs_bias = abs(median_prediction - actual_count),
           interval_width = q_95 - q_05)
  
  ####### Parameter Suumary
  parameter_name <- c( paste0("beta", 0:p),  "sigma2", "tau2", paste0("theta", 1:n))
  parameter_summary <- data.frame(mean = apply(output[, parameter_name], 2, mean),
                                  median = apply(output[, parameter_name], 2, median),
                                  q_05 = apply(output[, parameter_name], 2, quantile, .05),
                                  q_95 = apply(output[, parameter_name], 2, quantile, .95)) %>%
    tibble:: rownames_to_column("parameter") %>%
    dplyr:: select(parameter, everything())
  
  ###### Morans I
  morani_iter <- apply(output[, grepl("theta", names(output))], 1, function(theta) ape::Moran.I(theta, W, scaled = FALSE)$observed)
  morans_I_summary <- data.frame(MI_thbar = ape:: Moran.I(apply(output[, grepl("theta", names(output))], 2, mean), W, scaled = F)$observed,
                                 MI_mean = mean(morani_iter),
                                 MI_median = median(morani_iter),
                                 MI_ql = quantile(morani_iter, .05),
                                 MI_qh = quantile(morani_iter, .95)) %>%
    tibble::remove_rownames()
  
  ###### DIC
  lg_lkhd_thbar_train <- sum(dpois(count_var[train_locator], lambda = exp(apply(output[, paste0("log_lambda", 1:n_train)], 2, mean)), log = T))
  D_bar_train <- -2 * mean(output$lg_lkhd_train)
  D_thbar_train <- -2 * lg_lkhd_thbar_train
  p.d <- D_bar_train - D_thbar_train
  DIC <- D_bar_train + p.d
  
  ###### Log Likelihood of all data
  beta_bar = matrix(parameter_summary$mean[grepl("beta", parameter_summary$parameter)], ncol = 1)
  theta_bar = matrix(parameter_summary$mean[grepl("theta", parameter_summary$parameter)], ncol = 1)
  fitted_bar <- exp(X %*% beta_bar + theta_bar)
  lg_lkhd_alldata <- sum(dpois(count_var, lambda = t(fitted_bar), log = T))
  
  performance_metric <- data.frame(lglike_prameter_bar = lg_lkhd_thbar_train,
                                   D_bar = D_bar_train,
                                   D_theta_bar = D_thbar_train,
                                   p.d = p.d,
                                   DIC = DIC)
  
  # count_var_name <- colnames(rawData$count_data)[current_response]
  file_name <- paste0(substr(count_var_name, start = 3, stop = nchar(count_var_name)), "_fold_", fold_number)
  all_output_list <- list(test_locator = test_locator, tested_count = count_var[test_locator], sample_parameter = output,
                         parameter_summary = parameter_summary,
                         prediction_summary = prediction_summary,
                         morans_I_summary = morans_I_summary,
                         performance_metric = performance_metric,
                         loglike_alldata_mparam = lg_lkhd_alldata)
  
  # assign(list_name, all_output_list)
  save(all_output_list, file = paste0("./output_file/", file_name, ".RData"))
  # return(all_output_list)
  
  
  rm(output, epsilon_test, lambda_test, predicted_count, prediction_summary, parameter_name,
  parameter_summary, morans_I_summary, morani_iter, lg_lkhd_thbar_train, D_bar_train,D_thbar_train, p.d, DIC,
  beta_bar, theta_bar, fitted_bar, lg_lkhd_alldata, performance_metric, file_name, all_output_list)

  # rm(tempM)
  
}

# #=========== prediction:
# output <- poisson_spatial_draws
# # matrix(data = rnorm(m * 2, mean = means, sd = sds[c(1, 3)]), nrow = m, ncol = 2, byrow = T)
# 
# # X_test needs to include in function argument
# # n_test = 6
# # X_test = cbind(X_intercept = rep(1, n_test), rawData$covariate_data[test_locator, ])
# epsilon_test <- matrix(data = rnorm(n_test * nrow(output), mean = 0, sd = t(sqrt(output$sigma2))), nrow = n_test, ncol = nrow(output), byrow = T)
# lambda_test <- exp(X_test %*% t(output[, paste0("beta", 0:p)]) + t(output[, paste0("theta", test_locator)]) + epsilon_test)
# predicted_count <- apply(lambda_test, c(1, 2), function(x) rpois(1, x))
# # or,
# # check <- rpois(length(fitted.test), fitted.test)
# # dim(check) <- dim(fitted.test)
# prediction_summary <- data.frame(actual_count = count_var[test_locator],
#                                  median_prediction = apply(predicted_count, 1, median),
#                                  mean_prediction = apply(predicted_count, 1, mean),
#                                  q_05 = apply(predicted_count, 1, quantile, .05),
#                                  q_95 = apply(predicted_count, 1, quantile, .95))
# 
# # rm(epsilon_test, lambda_test, predicted_count, prediction_summary)
# #=========== prediction ends
# 
# ###### Parameter Suumary
# parameter_name <- c( paste0("beta", 0:p),  "sigma2", "tau2", paste0("theta", 1:n))
# parameter_summary <- data.frame(mean = apply(output[, parameter_name], 2, mean),
#                                 median = apply(output[, parameter_name], 2, median),
#                                 q_05 = apply(output[, parameter_name], 2, quantile, .05),
#                                 q_95 = apply(output[, parameter_name], 2, quantile, .95)) %>%
#   tibble:: rownames_to_column("parameter") %>%
#   dplyr:: select(parameter, everything())
# 
# # =========Parameter Suumary ends
# 
# ###### Performance metric
# 
# # Morans I
# morani_iter <- apply(output[, grepl("theta", names(output))], 1, function(theta) ape::Moran.I(theta, W, scaled = FALSE)$observed)
# morans_I_summary <- data.frame(MI_thbar = ape:: Moran.I(apply(output[, grepl("theta", names(output))], 2, mean), W, scaled = F)$observed,
#                                MI_mean = mean(morani_iter),
#                                MI_ql = quantile(morani_iter, .05),
#                                MI_qh = quantile(morani_iter, .95)) %>%
#                     tibble::remove_rownames()
# 
# #DIC:
# # # from carbayes package:
# # fitted <- exp(as.numeric(X.standardised %*% beta) + offset + phi)
# # loglike <- dpois(x=as.numeric(Y), lambda=fitted, log=TRUE)
# # fitted.mean <- exp(X.standardised %*% mean.beta + mean.phi + offset)
# # deviance.fitted <- -2 * sum(dpois(x=Y, lambda=fitted.mean, log=TRUE), na.rm=TRUE)
# #
# # mean.deviance <- -2 * sum(samples.loglike, na.rm=TRUE) /   nrow(samples.loglike)
# # p.d <- mean.deviance - deviance.fitted
# # DIC <- deviance.fitted + 2 * p.d
# 
# # cbind(X_intercept = rep(1, n_train), rawData$covariate_data[train_locator, ])
# # prediction_lambda <- exp( %*% t(output[, paste0("beta", 0:p)]) + t(output[, grepl("theta", names(output))]) +
# #                            matrix(data = rnorm(n_test * nrow(output), mean = 0, sd = t(output$sigma2)), nrow = n_test, ncol = nrow(output), byrow = T))
# # X <- cbind(X_intercept = rep(1, n_train + n_test), rawData$covariate_data)
# 
# # thbar <- apply(output[, parameter_name], 2, mean)
# # fitted.thbar_train <- exp( X_train %*% as.matrix(thbar[grepl("beta", names(thbar))]) + as.matrix(thbar[grepl("theta", names(thbar))][train_locator]) )
# # fitted.thbar_train <- exp(apply(output[, paste0("log_lambda", 1:n_train)], 2, mean))
# 
# lg_lkhd_thbar_train <- sum(dpois(count_var[train_locator], lambda = exp(apply(output[, paste0("log_lambda", 1:n_train)], 2, mean)), log = T))
# D_bar_train <- -2 * mean(output$lg_lkhd_train)
# D_thbar_train <- -2 * lg_lkhd_thbar_train
# p.d <- D_bar_train - D_thbar_train
# DIC <- D_bar_train + p.d
# # rm(morani_iter)
# 
# 
# dpois(count)







# poisson_non_spatial_model <- function (X_train, betas, sigma2_init, a_sigma_pr, 
#                                        b_sigma_pr, num_iter, burn_in, thin_in) {
#   set.seed(729)
#   #-----------------------
#   # data setup
#   #-----------------------
#   # n, p, X_train set up
#   n_train <- nrow(X_train)        # number of training obs.
#   tX_train = t(X_train)           # training covariate matrix
#   p <- ncol(X_train)-1
#   n <- n_train + n_test
#   ####### prior info:
#   # prior information for all beta, rest of the parameters prior coming from function argument
#   # Beta prior setup
#   mu_beta = matrix(c(2, 2, 1, -1, 0, 2), ncol = 1);  sigma_beta = 10^4 * diag(p + 1)
#   
#   ###### Assigning initial values 
#   beta = matrix(betas, ncol = 1); sigma2 = sigma2_init
#   rm(betas, sigma2_init)
#   
#   # iteration has no effect on estimates of A_sigma2 and sigma_beta_inv, 
#   # so calculating outside of the loop
#   A_sigma2 = (n_train/2) + a_sigma_pr                       # posterior scale parameter for sigma2
#   sigma_beta_inv = chol2inv(chol(sigma_beta))
#   
#   # Define overall Response (train + test) vector using train_response + initial test response
#   log_lambda = matrix(log_lamda_train, ncol = 1)
#   current_lambda = exp(log_lambda)
#   #output matrix
#   ncol_tempM = length(c(1, 1, beta, sigma2, log_lambda, log_lambda)) 
#   tempM <- matrix(0, nrow = num_iter, ncol = ncol_tempM)
#   
#   pb <- txtProgressBar(min = 1, max = num_iter, style = 3)
#   
#   for (j in 1:num_iter) {
#     
#     # tic = Sys.time()
#     
#     ###### update log_lambda: Metropolis-Hastings
#     # propose a log_lambda value from Normal
#     proposed_log_lambda <- rnorm(n_train, mean = log_lambda, sd = .02) # sd can be changed
#     proposed_lambda = exp(proposed_log_lambda)
#     mean_log_lambda = as.vector(X_train %*% beta)
#     # current_lambda <- exp(log_lambda)
#     # ratio
#     l_posterior_proposal =  -(2 * sigma2)^-1 * (proposed_log_lambda - mean_log_lambda)^2 - proposed_lambda + count_var[train_locator] * log(proposed_lambda)
#     l_posterior_current = - (2 * sigma2)^-1 * (log_lambda - mean_log_lambda)^2 - current_lambda + count_var[train_locator] * log(current_lambda)
#     l_target_by_proposal = l_posterior_proposal - l_posterior_current # because of taking log
#     
#     # accepting proposal
#     l_uniform_vector = log(runif(n_train, min = 0, max = 1))
#     gain <- l_uniform_vector < l_target_by_proposal
#     log_lambda[gain] <- proposed_log_lambda[gain] # updated log_lambda 
#     current_lambda[gain] <- proposed_lambda[gain] # required for next iteration
#     lambda_acc <- as.numeric(gain) # vector of log_lambda acceptance (binary)
#     acc_rate_per_iteration <- mean(lambda_acc)  # lambda_acc is binary
#     
#     
#     ###### Update theta: 
#     res_nonspatial <-  log_lambda - X_train %*% beta # (log_lambda- Xb) part in theta posterior
#     sum_sigma_w <- sigma2 * D_w_vector
#     theta_post_var = (sigma2 * tau2)/(tau2 + sum_sigma_w)
#     
#     ###### Update sigma2:
#     # what we get from full conditional of sigma2
#     # prior on sigma2 ~ IG( 2.5, 2.5)
#     # A_sigma2 = (n/2) + a_sigma_pr  # calculated before loop
#     ss <- sum((res_nonspatial)^2) # res_nonspatial = log_lambda - X_train*beta
#     B_sigma2 = ss/2 + b_sigma_pr
#     sigma2 <- 1/rgamma(1, shape = A_sigma2, rate = B_sigma2)
#     
#     
#     ###### Update beta:
#     # what we get from full conditional of beta
#     # sigma_beta_inv = chol2inv(chol(sigma_beta)) calculated before loop
#     # x_transpose_sigma_inv = tX_train %*% ((sigma2^-1) * diag(n))
#     x_transpose_sigma_inv = (sigma2^-1) * tX_train
#     A_beta = x_transpose_sigma_inv %*% X_train + sigma_beta_inv; A_beta_inv = chol2inv(chol(A_beta))
#     B_beta = x_transpose_sigma_inv %*% log_lambda + sigma_beta_inv %*% mu_beta
#     beta <- mvrnorm(1, A_beta_inv %*% B_beta, A_beta_inv)
#     
#     
#     ####### Log-likelihood
#     lambda = exp(log_lambda)
#     lg_lkhd_train <- sum(dpois(count_var[train_locator], lambda = lambda[train_locator], log = T))
#     # lg_lkhd_all <- sum(dpois(count_data_raw, lambda = lambda, log = T))
#     
#     # Accumulates all parameters for a single iteration
#     tempM[j, ] <-  c(lg_lkhd_train, acc_rate_per_iteration, beta, sigma2, log_lambda, lambda_acc)
#     
#     rm(B_sigma2, x_transpose_sigma_inv,
#        A_beta, B_beta, A_beta_inv, 
#        acc_rate_per_iteration, lambda_acc, proposed_log_lambda, lambda, lg_lkhd_train)
#     
#     setTxtProgressBar(pb, j) # to see the progress 
#   } 
#   
#   output <- as.data.frame(tempM) %>% 
#     .[(burn_in+1):num_iter, ] %>%
#     .[seq(from = 1, to = num_iter - burn_in, by = thin_in), ] %>%
#     set_names(c("lg_lkhd_train", "acc_rate_iter",  paste0("beta", 0:p),  "sigma2", 
#                 paste0("log_lambda", 1:n_train), paste0("lambda", 1:n_train, "_acc")))
#   
#   rm(tempM)
#   return(output) 
# }
