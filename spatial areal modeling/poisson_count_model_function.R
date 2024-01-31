poisson_count_model <- function (X, betas, theta_init, 
                                 sigma2_init, tau2_init, a_tau_pr, b_tau_pr, a_sigma_pr, 
                                 b_sigma_pr, num_iter, burn_in, thin_in) {
  set.seed(729)
  #-----------------------
  # data setup
  #-----------------------
  # n, p, X set up
  n <- nrow(X)        # n: number of training obs.
  tX = t(X)           # training covariate matrix
  p <- ncol(X)
  
  ####### prior info:
  # prior information for all beta, rest of the parameters prior coming from function argument
  # Beta prior setup
  mu_beta = matrix(c(2, 2, 1, -1, 0, 2), ncol = 1);  sigma_beta = 10^4 * diag(p + 1)
  
  ###### Assigning initial values 
  beta = matrix(betas, ncol = 1); theta = theta_init; tau2 = tau2_init; sigma2 = sigma2_init
  rm(betas, theta_init, tau2_init, sigma2_init)
  
  A_tau2 = ((n - eigen0_count)/2) + a_tau_pr # posterior scale parameter for tau2 
  A_sigma2 = (n/2) + a_sigma_pr # posterior scale parameter for sigma2
  sigma_beta_inv = chol2inv(chol(sigma_beta))
  
  # Define overall Response (train + test) vector using train_response + initial test response
  log_lambda = matrix(c(log_lamda_train, log_lamda_test_init), ncol = 1)
  current_lambda <- exp(log_lambda)
  #output matrix
  ncol_tempM = length(c(1, 1, 1, beta, sigma2, tau2, theta, log_lambda, log_lambda)) 
  tempM <- matrix(0, nrow = num_iter, ncol = ncol_tempM)
  
  pb <- txtProgressBar(min = 1, max = num_iter, style = 3)
  
  for (j in 1:num_iter) {
    
    # tic = Sys.time()
    
    # update log_lambda: Metropolis-Hastings
    # propose a log_lambda value from Normal
    mean_log_lambda <- as.vector(X %*% beta + theta)
    proposed_log_lambda <- rnorm(n, mean = log_lambda, sd = .5) # sd can be changed
    proposed_lambda <- exp(proposed_log_lambda)
    # current_lambda <- exp(log_lambda)
    # ratio
    posterior_proposal <- exp(- (2 * sigma2)^-1 * ((proposed_log_lambda - mean_log_lambda)^2) - proposed_lambda) * proposed_lambda^count_data_raw
    posterior_current <- exp(- (2 * sigma2)^-1 * ((log_lambda - mean_log_lambda)^2) - current_lambda) * current_lambda^count_data_raw
    
    target_by_proposal = posterior_proposal/posterior_current
    acceptance_ratio = if_else(target_by_proposal > 1, 1, target_by_proposal) # markov kernel cancel out
    
    # accepting proposal
    uniform_vector = runif(n, min = 0, max = 1)
    gain <- uniform_vector < acceptance_ratio
    log_lambda[gain] <- proposed_log_lambda[gain]
    current_lambda[gain] <- proposed_lambda[gain]
    lambda_acc <- as.numeric(gain) # vector of log_lambda acceptance (binary)
    acc_rate_per_iteration <- sum(lambda_acc)/n
    # gain <- which(uniform_vector < acceptance_ratio)
    # acc_rate_per_iteration <- length(gain)/n
    
    # Update theta: using functional form of posterior distribution
    # ** Note: in the case of combining train and test data altogether, 
    #         full conditional of theta has both likelihood and prior part.
    res_nonspatial <-  log_lambda - X %*% beta # (log_lambda- Xb) part in theta posterior mean 
    sum_sigma_w <- sigma2 * D_w_vector
    theta_post_var = (sigma2 * tau2)/(tau2 + sum_sigma_w)
    
    for(i in 1: n) {
      
      theta_post_mean_i =  (tau2 * res_nonspatial[i] + sigma2 * sum(theta[ W_list[[i]] ]))/(tau2 + sum_sigma_w[i])
      theta_post <- rnorm(1, theta_post_mean_i, sqrt(theta_post_var[i]))
      theta[i] <- theta_post
      
      rm(theta_post_mean_i, theta_post)
    }
    # theta <- theta - mean(theta)
    theta[island_locator] <- theta[island_locator] - mean(theta[island_locator])
    
    # Update tau2:
    # what we get from full conditional of tau2
    # prior on tau2 ~ IG( 2.5, 2.5)
    # posterior shape of tau2 is fixed and calculated before iteration loop
    # A_tau2 = ((n - eigen0_count)/2) + a_tau_pr
    # B_tau2 = (t(theta) %*% (D_w_tt - W_tt) %*% theta)/2 + b_tau_pr
    
    B_tau2 = .5 * (sum(D_w_vector * theta^2) - sum(theta * sapply(W_list, function(x) sum(theta[x])))) + b_tau_pr
    tau2 <- 1/rgamma(1, shape = A_tau2, rate = B_tau2)
    
    # Update sigma2:
    # what we get from full conditional of sigma2
    # prior on sigma2 ~ IG( 2.5, 2.5)
    # A_sigma2 = (n/2) + a_sigma_pr  # calculated before loop
    ss <- sum((res_nonspatial - theta)^2) # res_nonspatial = log_lambda - X*beta
    B_sigma2 = ss/2 + b_sigma_pr
    sigma2 <- 1/rgamma(1, shape = A_sigma2, rate = B_sigma2)
    
    # Update beta:
    # what we get from full conditional of beta
    # sigma_beta_inv = chol2inv(chol(sigma_beta)) calculated before loop
    # x_transpose_sigma_inv = tX %*% ((sigma2^-1) * diag(n))
    x_transpose_sigma_inv = (sigma2^-1) * tX
    A_beta = x_transpose_sigma_inv %*% X + sigma_beta_inv; A_beta_inv = chol2inv(chol(A_beta))
    B_beta = x_transpose_sigma_inv %*% (log_lambda - theta) + sigma_beta_inv %*% mu_beta
    beta <- mvrnorm(1, A_beta_inv %*% B_beta, A_beta_inv)
    
    
    # Log-likelihood
    lambda = exp(log_lambda)
    lg_lkhd_train <- sum(dpois(count_data_raw[-test_locator], lambda = lambda[-test_locator], log = T))
    lg_lkhd_all <- sum(dpois(count_data_raw, lambda = lambda, log = T))
    
    # Accumulates all parameters for a single iteration
    tempM[j, ] <-  c(lg_lkhd_train, lg_lkhd_all, acc_rate_per_iteration, beta, sigma2, tau2, theta, log_lambda, lambda_acc)
    
    rm(B_tau2, B_sigma2, x_transpose_sigma_inv,
       A_beta, B_beta, A_beta_inv, 
       acc_rate_per_iteration, lambda_acc, lambda, lg_lkhd_train, lg_lkhd_all)
    
    setTxtProgressBar(pb, j) # to see the progress 
    } 
  
  output <- as.data.frame(tempM) %>% 
    .[(burn_in+1):num_iter, ] %>%
    .[seq(from = 1, to = num_iter - burn_in, by = thin_in), ] %>%
    set_names(c("lg_lkhd_train", "lg_lkhd_all", "acc_rate_iter",  paste0("beta", 0:p),  "sigma2", "tau2", paste0("theta", 1:n), 
                paste0("l_lambda", 1:n), paste0("lambda", 1:n, "_acc")))
  
  rm(tempM)
  return(output) 
}

