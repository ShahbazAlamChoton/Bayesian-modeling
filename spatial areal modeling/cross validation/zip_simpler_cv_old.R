rm(list = ls())
graphics.off()
current_path <- rstudioapi::getSourceEditorContext()$path
setwd(dirname(current_path))
library(MCMCpack)
library(tidyverse)
library(dplyr)
library(DAAG)
library(igraph)
library(CARBayes)
set.seed(729)
options(scipen=0)

###### data setup
get_row_wise_adjacency_location <- function(row) {
  if (any(row == 1)) {
    return(which(row == 1)) 
  }
}

rawData <- get(load("DATA PATH"))
covariate_data <- rawData$covariate_data # excluding intercept
n <- nrow(covariate_data)
X <- matrix(cbind(1, covariate_data), nrow = nrow(covariate_data), byrow = F)
p <- ncol(X) - 1
###### Adjacency matrix
W <- rawData$W_mat
W_list <- apply(W, 1, get_row_wise_adjacency_location)
D_w <- diag(rowSums(W))
D_w_vector <- as.vector(rowSums(W))

###### Finding largest island
# igraph reference manual: Page 50
cluster_check <- igraph::graph_from_adjacency_matrix(W, mode = "undirected")
cluster_label <- unname(igraph::components(cluster_check)$membership)
table(cluster_label)
island_locator <- which(cluster_label == as.numeric(names(sort(table(cluster_label), decreasing = T)[1])))
eigen0_count <- length(unique(cluster_label)) #sum(eigen0 == 0)
rm(cluster_check, cluster_label)

# ===== command line for bird choice and MH sd choice
print(colSums(rawData$count_data!= 0, na.rm = T))
response <- as.numeric(readline("select response column number (1 - 10): "))
mh_sd <- as.numeric(readline("Type preferred sd value for log_lambda M-H: "))
# response <- 5 # un-comment if need to run loop for a single bird choice
#=====================================================

# Fold number for CV
cv_count <- 6 # hyper parameter
test_for_cv_matrix <- replicate(cv_count, sample(x = n, size = ceiling(n/cv_count), replace = FALSE))
# test_for_cv_matrix

param_test <- data.frame()
prediction_summary <- data.frame()

# CV loop starts
pb <- txtProgressBar(min = 1, max = cv_count, style = 3)
fold = 1
# for(fold in 1:cv_count) {
    
  fold_number <- fold # or fold_number
  test_loc <- test_for_cv_matrix[, fold_number]; len_test <- length(test_loc)
  train_loc <- seq_along(1:n)[-test_loc]
  
  X_train <- X[train_loc, ]; X_test <- X[test_loc, ]
  n_train <- nrow(X_train); n_test <- nrow(X_test)
  tX_train <- t(X_train)
  
  ###### initials for zip simpler model
  
  # parameters: 1) beta1, theta1; 2) beta2, theta2, sigma2; 3) omega; 4) z; 5) log_lambda
  
  num_iter = 20000; burn_in = 5000; thin_in = 30
  
  current_response <- response
  count_col <- rawData$count_data[, current_response] # hyper-parameter
  count_var <- count_col[train_loc]; count_test <- count_col[test_loc]
  
  count0_locator <- which(count_var == 0)
  count_var_name <- colnames(rawData$count_data)[current_response]
    
  ##### setting up initial thea1, theta2
  theta_init <- matrix(0, nrow = n, ncol = 1)
  
  ##### setting up initial z
  z_init <- rep(1, length(count_var))
  z_init[-count0_locator] <- 0
  
  ###### setting up initial beta2
  beta2_init <- unname(lm(log(count_var[z_init == 0] -1 +.25) ~ X_train[, -1][z_init == 0, ])$coef) # BETA2 initial value
  # is.na(beta2_init) <- 0
  
  ###### setting up initial log_lambda
  log_lambda_init <- numeric(length(count_var))
  log_lambda_init[z_init == 0] <- log(count_var[z_init == 0] -1 + 0.25)
  log_lambda_init[z_init == 1] <- X_train[z_init == 1, ] %*% beta2_init
  
  ###### setting up initial beta1
  beta1_0_init <- qnorm(mean(z_init))
  beta1_init <- c(beta1_0_init, rep(0, p))
  
  ###### setting up initial omega
  a <- if_else(z_init == 1, 0, -Inf)
  b <- if_else(z_init == 1, Inf, 0)
  mu_omega <- X_train %*% beta1_init + theta_init[train_loc, ]; sd_omega <- 1
  phi <- pnorm(mu_omega)
  A <- pnorm(a, mean = mu_omega, sd = sd_omega) + runif(n_train) * (pnorm(b, mean = mu_omega, sd = sd_omega) - pnorm(a, mean = mu_omega, sd = sd_omega))
  omega_init <- qnorm(A, mean = mu_omega, sd = sd_omega)
  
  # initial for sigma2, tau2_theta1, tau2_theta2
  sigma2_init = 1; tau2_init = 1
  
  ###### prior on sigma2 and tau2_theta1, tau2_theta2
  a_tau_pr = 2.01 ; b_tau_pr = 1.01
  a_sigma_pr = 2.01; b_sigma_pr = 1.01
  
  #-------------------- Start: ZIP simpler model
  
  set.seed(729)
  
  ####### prior info:
  # prior information for all beta1 and beta2 vector
  mu_beta1 = mu_beta2 = matrix(rep(0, p+1), ncol = 1); sigma_beta1 = sigma_beta2 = 10^4 * diag(p + 1)
  
  ###### Assigning initial values 
  beta1 = matrix(beta1_init, ncol = 1); beta2 = matrix(beta2_init, ncol = 1)
  theta1 = theta2 = theta_init; tau2_theta1 = tau2_theta2 = tau2_init; sigma2 = sigma2_init
  
  zee = z_init; zee0_locator <- which(zee==0); zee1_locator <- which(zee==1)
  omega = omega_init
  log_lambda = matrix(log_lambda_init, ncol = 1); current_lambda = exp(log_lambda)
  
  # iteration has no effect on estimates of tau2, A_sigma2 and sigma_beta_inv, 
  # so calculating outside of the loop
  A_tau2 = ((n_train - eigen0_count)/2) + a_tau_pr    # posterior scale parameter for tau2_theta1 and tau2_theta2
  A_sigma2 = length(zee0_locator)/2 + a_sigma_pr                 # posterior scale parameter for sigma2
  sigma_beta1_inv = chol2inv(chol(sigma_beta1)); sigma_beta2_inv = chol2inv(chol(sigma_beta2))
  
  #output matrix
  # ncol_tempM = length(c(test_loc, lambda_test, mu_test, y_tilda_test, omega_test, z_test, predicted_count,  beta1, theta1, beta2, theta2, acc_rate_per_iteration)
  ncol_tempM = 128
  tempM <- matrix(0, nrow = num_iter, ncol = ncol_tempM)
  
  for (j in 1:num_iter) {
    # tic = Sys.time()
    
    ###### update log_lambda: Metropolis-Hastings update for z = 0 and normal prior update for z = 1
    # propose a log_lambda value from Normal
    
    # change the sd value below to mh_sd when searching suitable sd==
    proposed_log_lambda <- rnorm(length(zee0_locator), mean = log_lambda[zee0_locator], sd = mh_sd) # sd can be changed mh_sd
    proposed_lambda = exp(proposed_log_lambda)
    mean_log_lambda = as.vector(X_train[zee0_locator, ] %*% beta2 + theta2[train_loc][zee0_locator])
    # current_lambda <- exp(log_lambda)
    # ratio
    l_posterior_proposal =  -(2 * sigma2)^-1 * (proposed_log_lambda - mean_log_lambda)^2 - proposed_lambda + (count_var[zee0_locator]-1) * log(proposed_lambda)
    l_posterior_current = - (2 * sigma2)^-1 * (log_lambda[zee0_locator] - mean_log_lambda)^2 - current_lambda[zee0_locator] + (count_var[zee0_locator] - 1) * log(current_lambda[zee0_locator])
    l_target_by_proposal = l_posterior_proposal - l_posterior_current # because of taking log
    
    # accepting proposal
    l_uniform_vector = log(runif(length(zee0_locator), min = 0, max = 1))
    gain = l_uniform_vector < l_target_by_proposal
    log_lambda[zee0_locator][gain] <- proposed_log_lambda[gain] # updated log_lambda 
    log_lambda[zee1_locator] <- rnorm(length(zee1_locator), mean = X_train[zee1_locator, ] %*% beta2 + theta2[train_loc][zee1_locator], sd = sqrt(sigma2))
    
    current_lambda <- exp(log_lambda) # required for next iteration
    lambda_acc <- as.numeric(gain) # vector of log_lambda acceptance (binary)
    acc_rate_per_iteration <- mean(lambda_acc)  # acceptance rate
    
    
    ###### Update theta2: 
    res_nonspatial_poisson <-  log_lambda[zee0_locator] - X_train[zee0_locator, ] %*% beta2 # (log_lambda- Xb) part in theta2 posterior
    sum_sigma_w_theta2 <- sigma2 * D_w_vector
    theta2_post_var = (sigma2 * tau2_theta2)/(tau2_theta2 + sum_sigma_w_theta2)
    
    # theta= matrix(0, nrow = n_train + n_test, ncol = 1)
    
    for(i in 1: length(zee0_locator)) {  # update depends on prior and likelihood
      theta2_position <- train_loc[zee0_locator][i]
      theta2_post_mean <-  (tau2_theta2 * res_nonspatial_poisson[i] + sigma2 * sum(theta2[W_list[[theta2_position]]]))/(tau2_theta2 + sum_sigma_w_theta2[theta2_position])
      theta2_post <- rnorm(1, theta2_post_mean, sqrt(theta2_post_var[theta2_position]))
      theta2[theta2_position] <- theta2_post
      
      rm(theta2_position, theta2_post_mean, theta2_post)
    } 
    
    for(i in 1: length(zee1_locator)) { # update depends on prior only
      theta2_position <- train_loc[zee1_locator][i]
      theta2_pred_mean <- sum(theta2[W_list[[theta2_position]]])/D_w_vector[theta2_position]
      theta2_pred_var <- tau2_theta2/D_w_vector[theta2_position]
      theta2_post_pred <- rnorm(1, theta2_pred_mean, sqrt(theta2_pred_var))
      theta2[theta2_position] <- theta2_post_pred
      
      rm(theta2_position, theta2_pred_mean, theta2_pred_var, theta2_post_pred)
    }
    
    for(i in 1:n_test) { # update depends on prior only
      theta2_position <- test_loc[i]
      theta2_pred_mean <- sum(theta2[W_list[[theta2_position]]])/D_w_vector[theta2_position]
      theta2_pred_var <- tau2_theta2/D_w_vector[theta2_position]
      theta2_post_pred <- rnorm(1, theta2_pred_mean, sqrt(theta2_pred_var))
      theta2[theta2_position] <- theta2_post_pred
      
      rm(theta2_position, theta2_pred_mean, theta2_pred_var, theta2_post_pred)
    }
    
    theta2[island_locator] <- theta2[island_locator] - mean(theta2[island_locator])
    theta2 <- matrix(0, ncol = 1, nrow = n)
    
    ###### Update tau2_theta2:
    # what we get from full conditional of tau2
    # prior on tau2 ~ IG( 2.5, 2.5)
    # posterior shape of tau2 is fixed and calculated before iteration loop
    # A_tau2 = ((n - eigen0_count)/2) + a_tau_pr
    # B_tau2 = (t(theta) %*% (D_w_tt - W_tt) %*% theta)/2 + b_tau_pr
    B_tau2_theta2 = .5 * (sum(D_w_vector * theta2^2) - sum(theta2 * sapply(W_list, function(x) sum(theta2[x])))) + b_tau_pr
    tau2_theta2 <- 1/rgamma(1, shape = A_tau2, rate = B_tau2_theta2) # shape parameter A_tau2 is common for both theta2 and theta1 vector
    
    
    ###### Update sigma2:
    # what we get from full conditional of sigma2
    # prior on sigma2 ~ IG( 2.5, 2.5)
    # A_sigma2 = (n/2) + a_sigma_pr  # calculated before loop
    ss_poisson <- sum((res_nonspatial_poisson - theta2[train_loc][zee0_locator])^2) # res_nonspatial = log_lambda - X_train*beta
    B_sigma2 = ss_poisson/2 + b_sigma_pr
    sigma2 <- 1/rgamma(1, shape = A_sigma2, rate = B_sigma2)
    
    
    ###### Update beta2:
    # what we get from full conditional of beta
    # sigma_beta_inv = chol2inv(chol(sigma_beta)) calculated before loop
    # x_transpose_sigma_inv = tX_train %*% ((sigma2^-1) * diag(n))
    x_transpose_sigma_inv_beta2 = (sigma2^-1) * tX_train[, zee0_locator]
    A_beta2 = x_transpose_sigma_inv_beta2 %*% X_train[zee0_locator, ] + sigma_beta2_inv; A_beta2_inv = chol2inv(chol(A_beta2))
    B_beta2 = x_transpose_sigma_inv_beta2 %*% (log_lambda[zee0_locator] - theta2[train_loc][zee0_locator]) + sigma_beta2_inv %*% mu_beta2
    beta2 <- mvrnorm(1, A_beta2_inv %*% B_beta2, A_beta2_inv)
    
    ###################################### 2nd part: bernoulli
    ###### Update omega:
    A = pnorm(a, mean = mu_omega, sd = sd_omega) + runif(n_train, 0, 1) * (pnorm(b, mean = mu_omega, sd = sd_omega) - pnorm(a, mean = mu_omega, sd = sd_omega))
    omega <- qnorm(A, mean = mu_omega, sd = sd_omega)
    
    ###### Update theta1: 
    # mu_omega <- X %*% beta1
    res_nonspatial_omega <-  omega - mu_omega # (omega- X * beta1) part in theta1 posterior
    sum_sigma_w_theta1 <- 1 * D_w_vector  #sigma2 * D_w_vector => here sigma2 = 1 
    theta1_post_var = (1 * tau2_theta1)/(tau2_theta1 + sum_sigma_w_theta1)
    
    # theta= matrix(0, nrow = n_train + n_test, ncol = 1)
    
    for(i in 1: n_train) {  # update depends on prior and likelihood
      theta1_position <- train_loc[i]
      theta1_post_mean <-  (tau2_theta1 * res_nonspatial_omega[i] + 1 * sum(theta1[W_list[[theta1_position]]]))/(tau2_theta1 + sum_sigma_w_theta1[theta1_position])
      theta1_post <- rnorm(1, theta1_post_mean, sqrt(theta1_post_var[theta1_position]))
      theta1[theta1_position] <- theta1_post
      # print(theta1_post)
      rm(theta1_position, theta1_post_mean, theta1_post)
    } 
    
    for(i in 1:n_test) { # update depends on prior only
      theta1_position <- test_loc[i]
      theta1_pred_mean <- sum(theta1[W_list[[theta1_position]]])/D_w_vector[theta1_position]
      theta1_pred_var <- tau2_theta1/D_w_vector[theta1_position]
      theta1_post_pred <- rnorm(1, theta1_pred_mean, sqrt(theta1_pred_var))
      theta1[theta1_position] <- theta1_post_pred
      
      rm(theta1_position, theta1_pred_mean, theta1_pred_var, theta1_post_pred)
    }
    theta1[island_locator] <- theta1[island_locator] - mean(theta1[island_locator])
    
    
    ###### Update tau2_theta1: same as tau2_theta2
    B_tau2_theta1 = .5 * (sum(D_w_vector * theta1^2) - sum(theta1 * sapply(W_list, function(x) sum(theta1[x])))) + b_tau_pr
    tau2_theta1 <- 1/rgamma(1, shape = A_tau2, rate = B_tau2_theta1) # shape parameter A_tau2 is common for both theta2 and theta1 vector
    
    
    ###### Update beta1: same as beta2
    x_transpose_sigma_inv_beta1 = 1 * tX_train #  (sigma2^-1) * tX => here, sigma2 = 1
    A_beta1 = x_transpose_sigma_inv_beta1 %*% X_train + sigma_beta1_inv; A_beta1_inv = chol2inv(chol(A_beta1))
    B_beta1 = x_transpose_sigma_inv_beta1 %*% (omega - theta1[train_loc]) + sigma_beta1_inv %*% mu_beta1
    beta1 <- mvrnorm(1, A_beta1_inv %*% B_beta1, A_beta1_inv)
    
    mu_omega <- X_train %*% beta1 + theta1[train_loc]; sd_omega = 1 # required for next iteration
    phi <- pnorm(mu_omega)                        #required for next iteration
    
    # #============================== Prediction
    # # ####### first part: test_lambda
    # lambda_test <- exp(rnorm(length(test_loc), X_test %*% beta2 + theta2[test_loc], sd = sqrt(sigma2)))
    # 
    # # cat("log_lambda from non-zero count:", round(exp(log_lambda[zee0_locator]), 2), "\n\n")
    # # cat("actual count:", count_col[test_loc])
    # # cat("lambda_test", lambda_test, "\n\n")
    # # cat("X_test * beta2", X_test %*% beta2, "\n\n")
    # 
    # 
    # y_tilda_test <- rpois(length(test_loc), lambda = lambda_test)
    # 
    # # ###### 2nd part:
    # # omega_test
    # mu_test <- X_test %*% beta1 + theta1[test_loc]
    # omega_test <- rnorm(length(test_loc), mean = mu_test, sd = 1)
    # z_test <- if_else(omega_test > 0, 1, 0)
    # 
    # # Prediction for test data using ZIP: count = (1 - Z_i) * Y_tilda_i;
    # # if Z = 0 -> count comes from Poisson part; if Z = 1 -> count = 0
    # predicted_count <- if_else(z_test == 0, y_tilda_test - 1, 0)
    # #==================================== Prediction End ======================
  
    # Accumulates all parameters for a single iteration
    # ncol_tempM = length(c(1, 1, beta1, beta2, sigma2, tau2_theta1, tau2_theta2, theta1, theta2, log_lambda, zee)) 
    tempM[j, ] <-  c(test_loc, lambda_test, mu_test, y_tilda_test, omega_test, z_test, predicted_count, beta1, theta1, beta2, theta2, sigma2, acc_rate_per_iteration)
    
    # if(j == 10) break 
    
    rm(acc_rate_per_iteration, lambda_acc, proposed_log_lambda, 
       res_nonspatial_poisson, sum_sigma_w_theta2,  theta2_post_var, B_tau2_theta2, 
       ss_poisson, B_sigma2,
       x_transpose_sigma_inv_beta2, A_beta2, A_beta2_inv, B_beta2, 
       A,
       res_nonspatial_omega, sum_sigma_w_theta1, theta1_post_var, B_tau2_theta1, 
       x_transpose_sigma_inv_beta1, A_beta1, A_beta1_inv, B_beta1)
  } 

  output <- as.data.frame(tempM) %>%
    # .[(burn_in+1):num_iter, ] %>%
    # .[seq(from = 1, to = num_iter - burn_in, by = thin_in), ] %>%
    set_names(c(paste0("test_loc_", 1:n_test), paste0("lambda_test_", 1:n_test), paste0("mu_test_", 1:n_test), paste0("ytilda_test_", 1:n_test), 
                paste0("omega_test_", 1:n_test), paste0("zee_test_", 1:n_test),   paste0("predicted_count_", 1:n_test),
                paste0("beta1_", 0:p), paste0("theta1_", 1:n), paste0("beta2_", 0:p), paste0("theta2_", 1:n), "sigma2", "mean_acc_rate" )) %>% 
    mutate(fold_number = fold)
  
#   
#   param_test <- rbind(param_test, output)
#   
#   prediction_temp <- data.frame(fold_number = fold, 
#                                 test_location = test_loc, 
#                                 actual_count = count_col[test_loc],
#                                 median_prediction = apply(output[, paste0("predicted_count_", 1:n_test)], 2, median, na.rm = T),
#                                 zee1_prop = apply(output[, paste0("zee_test_", 1:n_test)], 2, mean, na.rm = T),
#                                 y_tilda_mean = apply(output[, paste0("ytilda_test_", 1:n_test)], 2, mean, na.rm = T),
#                                 y_tilda_median = apply(output[, paste0("ytilda_test_", 1:n_test)], 2, median, na.rm = T),
#                                 mean_prediction = apply(output[, paste0("predicted_count_", 1:n_test)], 2, mean, na.rm = T),
#                                 pred_q05 = apply(output[, paste0("predicted_count_", 1:n_test)], 2, quantile, .05, na.rm = T),
#                                 pred_q95 = apply(output[, paste0("predicted_count_", 1:n_test)], 2, quantile, .95, na.rm = T)
#   ) %>% 
#     
#     mutate(include = if_else(actual_count >= pred_q05 & actual_count <= pred_q95, "Yes", "No"),
#            abs_bias = abs(median_prediction - actual_count),
#            interval_width = pred_q95 - pred_q05) %>% 
#     select(actual_count, median_prediction, pred_q05, pred_q95, zee1_prop, y_tilda_median, include, abs_bias, interval_width, fold_number, test_location, y_tilda_mean, everything()) %>% 
#     remove_rownames()
#   
#   prediction_summary <- rbind(prediction_summary, prediction_temp)
#   
#   rm(tempM, output, prediction_temp)
#   
#   setTxtProgressBar(pb, fold) # to see the progress 
# 
# }
# 
# 
# # ================= Test Summary: count 0, non 0 splitter
# # parameter_summary <- param_test
# test_summary_row = (num_iter - burn_in)/thin_in
# mu_test_summary <- lambda_test_summary <-  matrix(0, nrow = test_summary_row, ncol = n)
# 
# for(i in 1:cv_count) {
#   mu_test_summary[, ((i - 1) * n_test + 1) : (i * n_test)] = as.matrix(param_test[((i - 1) * test_summary_row + 1) : (i * test_summary_row),  paste0("mu_test_", 1:n_test)])
# }
# 
# mu_summary <- data.frame(mu_median = apply(mu_test_summary, 2, median),
#                          mu_rangeL = apply(mu_test_summary, 2, min), 
#                          mu_rangeU = apply(mu_test_summary, 2, max) )
# 
# for(i in 1:cv_count) {
#   lambda_test_summary[, ((i - 1) * n_test + 1) : (i * n_test)] = as.matrix(param_test[((i - 1) * test_summary_row + 1) : (i * test_summary_row),  paste0("lambda_test_", 1:n_test)])
#   # lambdaSummary[, ((i - 1) * 6 + 1) : (i * 6)] = as.matrix(parameter_summary[((i - 1) * 500 + 1) : (i * 500),  paste0("lambda_test_", 1:n_test)])
# }
# elambda_test_matrix <- exp(-lambda_test_summary)
# elambda_summary <- data.frame(elambda_median = apply(elambda_test_matrix, 2, median), 
#                               elambda_rangeL = apply(elambda_test_matrix, 2, min), 
#                               elambda_rangeU = apply(elambda_test_matrix, 2, max) )
# 
# # # Probability of having zero count: P(count = 0) = P(Y_tilda = 0|Z_i = 0) * P(Z_i = 0) + P(Y_tilda = 0|Z_i = 1) * P(Z_i = 1)
# # #                                                  Poisoon(0) * P(Z_i = 0) + P(Y_tilda = 0|Z_i = 1) * P(Z_i = 1)
# prob_being_count0 <- pnorm(mu_test_summary)
# count_splitter <- data.frame(actual_non0_count_binary = as.numeric(prediction_summary$actual_count > 0),
#                              prob_count0 = apply(prob_being_count0, 2, median)) %>% 
#   mutate(prob_count0_new = if_else(prob_count0 < .01, 0, if_else(prob_count0 > .99, 1, prob_count0)) ) %>% 
#   select(actual_non0_count_binary, prob_count0_new, prob_count0)
# 
# splitter_inacc <- mean(count_splitter$actual_non0_count_binary - count_splitter$prob_count0_new == 0)
# # ========================================= Splitter End ==========================================
# 
# all_output_list <- list(test_parameter = param_test,
#                         prediction_summary = prediction_summary,
#                         mu_summary = mu_summary,
#                         elambda_summary = elambda_summary,
#                         count_splitter = count_splitter,
#                         splitter_acc = 1 - splitter_inacc)
# 
# rm(test_summary_row, mu_test_summary, mu_summary,
#    lambda_test_summary, elambda_test_matrix, elambda_summary, 
#    prob_being_count0, count_splitter, splitter_acc)
# 
# # comment first line and uncomment 2nd line for file name when run usual==
# # file_name <- paste0(substr(count_var_name, start = 3, stop = nchar(count_var_name)), "_zip_sd_", mh_sd, "acc_", acc)
# file_name <- paste0(substr(count_var_name, start = 3, stop = nchar(count_var_name)), "_zipSimpler_cv")
# save(all_output_list, file = paste0("./output/", file_name, ".RData"))
# #==========================================================================
# 
