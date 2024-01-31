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
options(scipen=999)

###### data setup
get_row_wise_adjacency_location <- function(row) {
  if (any(row == 1)) {
    return(which(row == 1)) 
  }
}
# source("poisson_spatial_model.R")

rawData <- get(load("DATA PATH"))
covariate_data <- rawData$covariate_data # excluding intercept
X <- matrix(cbind(1, covariate_data), nrow = nrow(covariate_data), byrow = F)
tX <- t(X)
n <- nrow(covariate_data); p <- ncol(covariate_data)
# cv_count <- 6 # hyper parameter
# test_for_cv_matrix <- replicate(cv_count, sample(x = n, size = ceiling(n/cv_count), replace = FALSE))

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

###### initials for zip model

# parameters: 1) beta1, theta1; 2) beta2, theta2, sigma2; 3) omega; 4) z; 5) log_lambda

num_iter = 20000; burn_in = 5000; thin_in = 5
print(colSums(rawData$count_data!= 0, na.rm = T))
# response <- 5
response <- as.numeric(readline("select response column number (1 - 10): "))
# fold_number <- 1
current_response <- response
count_var <- rawData$count_data[, current_response] # hyper-parameter
count0_locator <- which(count_var == 0)
count_var_name <- colnames(rawData$count_data)[current_response]
# column wise test locations for CV
# test_locator <- test_for_cv_matrix[, fold_number]

##### setting up initial thea1, theta2
theta_init <- matrix(0, nrow = n, ncol = 1)

##### setting up initial z
zero_location <- c(which(count_var > 0), sample(x = which(count_var == 0), size = ceiling(length(which(count_var==0))/2), replace = F))
z_init <- numeric(length(count_var))
z_init[-zero_location] <- 1

###### setting up initial beta2
beta2_init <- unname(lm(log(count_var[z_init == 0] + 0.25) ~ covariate_data[z_init == 0, ])$coef) # β(2) initial value

###### setting up initial log_lambda
log_lambda_init <- numeric(length(count_var))
log_lambda_init[z_init == 0] <- log(count_var[z_init == 0] + 0.25)
log_lambda_init[z_init == 1] <- X[z_init == 1, ] %*% beta2_init

###### setting up initial beta1
beta1_0_init <- qnorm(mean(z_init))
beta1_init <- c(beta1_0_init, rep(0, p))

###### setting up initial omega
a <- if_else(z_init == 1, 0, -Inf)
b <- if_else(z_init == 1, Inf, 0)
mu_omega <- X %*% beta1_init + theta_init; sd_omega <- 1
phi <- pnorm(mu_omega)
A <- pnorm(a, mean = mu_omega, sd = sd_omega) + runif(n) * (pnorm(b, mean = mu_omega, sd = sd_omega) - pnorm(a, mean = mu_omega, sd = sd_omega))
omega_init <- qnorm(A, mean = mu_omega, sd = sd_omega)

# initial for sigma2, tau2_theta1, tau2_theta2
sigma2_init = 1; tau2_init = 1

###### prior on sigma2 and tau2_theta1, tau2_theta2
a_tau_pr = 2.01 ; b_tau_pr = 1.01
a_sigma_pr = 2.01; b_sigma_pr = 1.01

# (tau_var <- b_tau_pr^2/((a_tau_pr -1)^2 * (a_tau_pr - 2))) # Variance formula for IG

#-------------------- Start: ZIP model

set.seed(729)

####### prior info:
# prior information for all beta1 and beta2 vector
mu_beta1 = mu_beta2 = matrix(c(2, 2, 1, -1, 0, 2), ncol = 1); sigma_beta1 = sigma_beta2 = 10^4 * diag(p + 1)

###### Assigning initial values 
beta1 = matrix(beta1_init, ncol = 1); beta2 = matrix(beta2_init, ncol = 1)
theta1 = theta2 = theta_init; tau2_theta1 = tau2_theta2 = tau2_init; sigma2 = sigma2_init

zee = z_init; omega = omega_init; log_lambda = matrix(log_lambda_init, ncol = 1)
current_lambda = exp(log_lambda)

# iteration has no effect on estimates of tau2, A_sigma2 and sigma_beta_inv, 
# so calculating outside of the loop
A_tau2 = ((n - eigen0_count)/2) + a_tau_pr    # posterior scale parameter for tau2_theta1 and tau2_theta2
A_sigma2 = (n/2) + a_sigma_pr                 # posterior scale parameter for sigma2
sigma_beta1_inv = chol2inv(chol(sigma_beta1)); sigma_beta2_inv = chol2inv(chol(sigma_beta2))

#output matrix
# c(lg_lkhd, acc_rate_per_iteration, beta1, beta2, sigma2, tau2_theta1, tau2_theta2, theta1, theta2, log_lambda, zee)
ncol_tempM = length(c(1, 1, beta1, beta2, sigma2, tau2_theta1, tau2_theta2, theta1, theta2, log_lambda, which(count_var == 0))) 
tempM <- matrix(0, nrow = num_iter, ncol = ncol_tempM)

pb <- txtProgressBar(min = 1, max = num_iter, style = 3)
for (j in 1:num_iter) {
  # tic = Sys.time()
  
  ###### update log_lambda: Metropolis-Hastings update for z = 0 and normal prior update for z = 1
  zee0_locator <- which(zee==0)
  zee1_locator <- which(zee==1)
  
  # propose a log_lambda value from Normal
  proposed_log_lambda <- rnorm(length(zee0_locator), mean = log_lambda[zee0_locator], sd = .615) # sd can be changed
  proposed_lambda = exp(proposed_log_lambda)
  mean_log_lambda = as.vector(X[zee0_locator, ] %*% beta2 + theta2[zee0_locator])
  # current_lambda <- exp(log_lambda)
  # ratio
  l_posterior_proposal =  -(2 * sigma2)^-1 * (proposed_log_lambda - mean_log_lambda)^2 - proposed_lambda + count_var[zee0_locator] * log(proposed_lambda)
  l_posterior_current = - (2 * sigma2)^-1 * (log_lambda[zee0_locator] - mean_log_lambda)^2 - current_lambda[zee0_locator] + count_var[zee0_locator] * log(current_lambda[zee0_locator])
  l_target_by_proposal = l_posterior_proposal - l_posterior_current # because of taking log
  
  # accepting proposal
  l_uniform_vector = log(runif(length(zee0_locator), min = 0, max = 1))
  gain = l_uniform_vector < l_target_by_proposal
  log_lambda[zee0_locator][gain] <- proposed_log_lambda[gain] # updated log_lambda 
  log_lambda[zee1_locator] <- rnorm(length(zee1_locator), mean = X[zee1_locator, ] %*% beta2 + theta2[zee1_locator], sd = sqrt(sigma2))
  
  current_lambda <- exp(log_lambda) # required for next iteration
  lambda_acc <- as.numeric(gain) # vector of log_lambda acceptance (binary)
  acc_rate_per_iteration <- mean(lambda_acc)  # acceptance rate
  
  
  ###### Update theta2: 
  res_nonspatial_poisson <-  log_lambda - X %*% beta2 # (log_lambda- Xb) part in theta2 posterior
  sum_sigma_w_theta2 <- sigma2 * D_w_vector
  theta2_post_var = (sigma2 * tau2_theta2)/(tau2_theta2 + sum_sigma_w_theta2)
  
  # theta= matrix(0, nrow = n_train + n_test, ncol = 1)
  
  for(i in 1: n) {  # update depends on prior and likelihood
    theta2_position <- i
    theta2_post_mean <-  (tau2_theta2 * res_nonspatial_poisson[i] + sigma2 * sum(theta2[W_list[[theta2_position]]]))/(tau2_theta2 + sum_sigma_w_theta2[theta2_position])
    theta2_post <- rnorm(1, theta2_post_mean, sqrt(theta2_post_var[theta2_position]))
    theta2[theta2_position] <- theta2_post
    
    rm(theta2_position, theta2_post_mean, theta2_post)
  } 
  
  # for(i in 1:n_test) { # update depends on prior only
  #   theta_position <- test_locator[i]
  #   theta_pred_mean <- sum(theta[W_list[[theta_position]]])/D_w_vector[theta_position]
  #   theta_pred_var <- tau2/D_w_vector[theta_position]
  #   theta_post_pred <- rnorm(1, theta_pred_mean, sqrt(theta_pred_var))
  #   theta[theta_position] <- theta_post_pred
  #   
  #   rm(theta_position, theta_pred_mean, theta_pred_var, theta_post_pred)
  # }
  theta2[island_locator] <- theta2[island_locator] - mean(theta2[island_locator])
  
  
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
  ss_poisson <- sum((res_nonspatial_poisson - theta2)^2) # res_nonspatial = log_lambda - X_train*beta
  B_sigma2 = ss_poisson/2 + b_sigma_pr
  sigma2 <- 1/rgamma(1, shape = A_sigma2, rate = B_sigma2)
  
  
  ###### Update beta2:
  # what we get from full conditional of beta
  # sigma_beta_inv = chol2inv(chol(sigma_beta)) calculated before loop
  # x_transpose_sigma_inv = tX_train %*% ((sigma2^-1) * diag(n))
  x_transpose_sigma_inv_beta2 = (sigma2^-1) * tX
  A_beta2 = x_transpose_sigma_inv_beta2 %*% X + sigma_beta2_inv; A_beta2_inv = chol2inv(chol(A_beta2))
  B_beta2 = x_transpose_sigma_inv_beta2 %*% (log_lambda - theta2) + sigma_beta2_inv %*% mu_beta2
  beta2 <- mvrnorm(1, A_beta2_inv %*% B_beta2, A_beta2_inv)
  
  ###################################### 2nd part: bernoulli
  ###### Update zee:
  zee[-count0_locator] = 0
  phi_count0_location <- phi[count0_locator]
  dpois0 <- dpois(0, lambda = current_lambda[count0_locator])
  post_prob_zee1 <- (1 * phi_count0_location)/(1 * phi_count0_location + dpois0 * (1 - phi_count0_location))
  runif_count0 <- runif(length(count0_locator))
  zee[count0_locator] <- if_else(runif_count0 <= post_prob_zee1, 1, 0)
  
  ###### Update omega:
  a = if_else(zee == 1, 0, -Inf)
  b = if_else(zee == 1, Inf, 0)
  A = pnorm(a, mean = mu_omega, sd = sd_omega) + runif(n, 0, 1) * (pnorm(b, mean = mu_omega, sd = sd_omega) - pnorm(a, mean = mu_omega, sd = sd_omega))
  omega <- qnorm(A, mean = mu_omega, sd = sd_omega)
  
  ###### Update theta1: 
  # mu_omega <- X %*% beta1
  res_nonspatial_omega <-  omega - mu_omega # (omega- X * beta1) part in theta1 posterior
  sum_sigma_w_theta1 <- 1 * D_w_vector  #sigma2 * D_w_vector => here sigma2 = 1 
  theta1_post_var = (1 * tau2_theta1)/(tau2_theta1 + sum_sigma_w_theta1)
  
  # theta= matrix(0, nrow = n_train + n_test, ncol = 1)
  
  for(i in 1: n) {  # update depends on prior and likelihood
    theta1_position <- i
    theta1_post_mean <-  (tau2_theta1 * res_nonspatial_omega[i] + 1 * sum(theta1[W_list[[theta1_position]]]))/(tau2_theta1 + sum_sigma_w_theta1[theta1_position])
    theta1_post <- rnorm(1, theta1_post_mean, sqrt(theta1_post_var[theta1_position]))
    theta1[theta1_position] <- theta1_post
    # print(theta1_post)
    rm(theta1_position, theta1_post_mean, theta1_post)
  } 

  theta1[island_locator] <- theta1[island_locator] - mean(theta1[island_locator])
  
  
  ###### Update tau2_theta1: same as tau2_theta2
  B_tau2_theta1 = .5 * (sum(D_w_vector * theta1^2) - sum(theta1 * sapply(W_list, function(x) sum(theta1[x])))) + b_tau_pr
  tau2_theta1 <- 1/rgamma(1, shape = A_tau2, rate = B_tau2_theta1) # shape parameter A_tau2 is common for both theta2 and theta1 vector
  
  
  ###### Update beta1: same as beta2
  x_transpose_sigma_inv_beta1 = 1 * tX #  (sigma2^-1) * tX => here, sigma2 = 1
  A_beta1 = x_transpose_sigma_inv_beta1 %*% X + sigma_beta1_inv; A_beta1_inv = chol2inv(chol(A_beta1))
  B_beta1 = x_transpose_sigma_inv_beta1 %*% (omega - theta1) + sigma_beta1_inv %*% mu_beta1
  beta1 <- mvrnorm(1, A_beta1_inv %*% B_beta1, A_beta1_inv)
  
  mu_omega <- X %*% beta1 + theta1; sd_omega = 1 # required for next iteration
  phi <- pnorm(mu_omega)                        #required for next iteration
  
  ####### Log-likelihood for all and training data
  # fitted.mean_train = exp(log_lambda)
  # # fitted.mean_train <- exp(X_train %*% beta + theta[train_locator])
  # lg_lkhd_train <- sum(dpois(count_var[train_locator], lambda = fitted.mean_train, log = T))
  
  count_prob = dpois(count_var, lambda = current_lambda) * (1- phi)
  count_prob[count0_locator] = count_prob[count0_locator] + phi[count0_locator]
  lg_lkhd <- sum(log(count_prob))
 

  # Accumulates all parameters for a single iteration
  # ncol_tempM = length(c(1, 1, beta1, beta2, sigma2, tau2_theta1, tau2_theta2, theta1, theta2, log_lambda, zee)) 
  # check <- c(lg_lkhd, acc_rate_per_iteration, beta1, beta2, sigma2, tau2_theta1, tau2_theta2, theta1, theta2, log_lambda)
  tempM[j, ] <-  c(lg_lkhd, acc_rate_per_iteration, beta1, beta2, sigma2, tau2_theta1, tau2_theta2, theta1, theta2, log_lambda, zee[count_var == 0])
  
  rm(acc_rate_per_iteration, lambda_acc, proposed_log_lambda, 
     res_nonspatial_poisson, sum_sigma_w_theta2,  theta2_post_var, B_tau2_theta2, 
     ss_poisson, B_sigma2,
     x_transpose_sigma_inv_beta2, A_beta2, A_beta2_inv, B_beta2,
     phi_count0_location, dpois0, runif_count0,
     a, b, A,
     res_nonspatial_omega, sum_sigma_w_theta1, theta1_post_var, B_tau2_theta1, 
     x_transpose_sigma_inv_beta1, A_beta1, A_beta1_inv, B_beta1,
     count_prob, lg_lkhd)
  
  # mean(tempM[, 2] )
  setTxtProgressBar(pb, j) # to see the progress 
} 

output <- as.data.frame(tempM) %>%
  .[(burn_in+1):num_iter, ] %>%
  .[seq(from = 1, to = num_iter - burn_in, by = thin_in), ] %>%
  set_names(c("lg_lkhd", "acc_rate_iter",  paste0("beta1_", 0:p), paste0("beta2_", 0:p),  "sigma2", "tau2_theta1", "tau2_theta2", paste0("theta1_", 1:n), paste0("theta2_", 1:n),
              paste0("log_lambda", 1:n), paste0("zee_count0_", 1:length(which(count_var == 0))) ))

# cat("\n", mean(output$acc_rate_iter))

# ####### Count Prediction
# matrix(data = rnorm(6, mean = 0, sd = c(1, 100)), nrow = 3, ncol = 2, byrow = T)
# 
# epsilon_test <- matrix(data = rnorm(n_test * nrow(output), mean = 0, sd = sqrt(output$sigma2)), nrow = n_test, ncol = nrow(output), byrow = T)
# lambda_test <- exp(X_test %*% t(output[, paste0("beta", 0:p)]) + t(output[, paste0("theta", test_locator)]) + epsilon_test)
# predicted_count <- apply(lambda_test, c(1, 2), function(x) rpois(1, x))
# 
# prediction_summary <- data.frame(actual_count = count_var[test_locator],
#                                  median_prediction = apply(predicted_count, 1, median),
#                                  mean_prediction = apply(predicted_count, 1, mean),
#                                  q_05 = apply(predicted_count, 1, quantile, .05),
#                                  q_95 = apply(predicted_count, 1, quantile, .95)) %>% 
#   mutate(include = if_else(actual_count >= q_05 & actual_count <= q_95, "Yes", "No"),
#          abs_bias = abs(median_prediction - actual_count),
#          interval_width = q_95 - q_05)
# 
# ####### Parameter Suumary
# parameter_name <- c( paste0("beta", 0:p),  "sigma2", "tau2", paste0("theta", 1:n))
# parameter_summary <- data.frame(mean = apply(output[, parameter_name], 2, mean),
#                                 median = apply(output[, parameter_name], 2, median),
#                                 q_05 = apply(output[, parameter_name], 2, quantile, .05),
#                                 q_95 = apply(output[, parameter_name], 2, quantile, .95)) %>%
#   tibble:: rownames_to_column("parameter") %>%
#   dplyr:: select(parameter, everything())

# paste0("beta1_", 0:p), paste0("beta2_", 0:p),  "sigma2", "tau2_theta1", "tau2_theta2", paste0("theta1_", 1:n), paste0("theta2_", 1:n),
# paste0("log_lambda", 1:n), paste0("zee_count0_", 1:length(which(count_var == 0)))

# coda::geweke.diag(output$beta1_0, frac1 = .1, frac2 = .5)$z
# 2*(1-pnorm(abs(coda::geweke.diag(output$beta1_0, frac1 = .1, frac2 = .5)$z)))

geweke_stat <- function(x){
  geweke_z = unname(coda::geweke.diag(x)$z)
  p_val = 2*(1-pnorm(abs(geweke_z)))
  # converge_y = if_else(p_val < .05, "y", "n")
  return(c(geweke_z = geweke_z, p_val = p_val))
}
# 2*pnorm(q=1.24, lower.tail=FALSE)
# test <- as.data.frame(t((apply(output[, c("beta1_0", "beta1_1")], 2, geweke_stat))))

parameter_name <- c( paste0("beta1_", 0:p), paste0("beta2_", 0:p),"sigma2", "tau2_theta1", "tau2_theta2", paste0("theta1_", 1:n), paste0("theta2_", 1:n), paste0("log_lambda", 1:n))
parameter_summary <- data.frame(mean = apply(output[, parameter_name], 2, mean),
                                median = apply(output[, parameter_name], 2, median),
                                q_05 = apply(output[, parameter_name], 2, quantile, .05),
                                q_95 = apply(output[, parameter_name], 2, quantile, .95)) %>% 
  bind_cols(as.data.frame(t((apply(output[, parameter_name], 2, geweke_stat))))) %>% 
  mutate(converge_y = if_else(p_val <.05, "y", "n")) %>% 
  tibble:: rownames_to_column("parameter") %>%
  # bind_cols(apply(output[, parameter_name], 2, geweke_stat)) %>% 
  dplyr:: select(parameter, everything())

# ###### Morans I
# morani_iter <- apply(output[, grepl("theta", names(output))], 1, function(theta) ape::Moran.I(theta, W, scaled = FALSE)$observed)
# morans_I_summary <- data.frame(MI_thbar = ape:: Moran.I(apply(output[, grepl("theta", names(output))], 2, mean), W, scaled = F)$observed,
#                                MI_mean = mean(morani_iter),
#                                MI_median = median(morani_iter),
#                                MI_ql = quantile(morani_iter, .05),
#                                MI_qh = quantile(morani_iter, .95)) %>%
#   tibble::remove_rownames()

### For theta 1 vector
morans_I_summary <- data.frame()
for( i in 1:2) { 
  theta_morani_iter <- apply(output[, grepl(paste0("theta",i,"_"), names(output))], 1, function(theta) ape::Moran.I(theta, W, scaled = FALSE)$observed)
  temp <- data.frame(theta = paste0("theta",i), 
                     MI_thbar = ape:: Moran.I(apply(output[, grepl(paste0("theta",i,"_"), names(output))], 2, mean), W, scaled = F)$observed, 
                     MI_mean = mean(theta_morani_iter), 
                     MI_median = median(theta_morani_iter), 
                     MI_ql = quantile(theta_morani_iter, .05), 
                     MI_qh = quantile(theta_morani_iter, .95)) %>% 
    tibble::remove_rownames()
  morans_I_summary <- bind_rows(morans_I_summary, temp)
  }
# print(morans_I_summary)

# ###### DIC
# lg_lkhd_thbar_train <- sum(dpois(count_var[train_locator], lambda = exp(apply(output[, paste0("log_lambda", 1:n_train)], 2, mean)), log = T))
# D_bar_train <- -2 * mean(output$lg_lkhd_train)
# D_thbar_train <- -2 * lg_lkhd_thbar_train
# p.d <- D_bar_train - D_thbar_train
# DIC <- D_bar_train + p.d

mu_omega_bar <- X %*% parameter_summary[parameter_summary$parameter %in% paste0("beta1_", 0:p), "mean"] + parameter_summary[parameter_summary$parameter %in% paste0("theta1_", 1:n), "mean"]
phi_bar <- pnorm(mu_omega)  
count_prob_bar = dpois(count_var, lambda = exp(parameter_summary[parameter_summary$parameter %in% paste0("log_lambda", 1:n), "mean"])) * (1- phi_bar)
count_prob_bar[count0_locator] = count_prob_bar[count0_locator] + phi_bar[count0_locator]
lg_lkhd_thbar <- sum(log(count_prob_bar))

D_bar <- -2 * mean(output$lg_lkhd)
D_thbar <- -2 * lg_lkhd_thbar
p.d <- D_bar - D_thbar
DIC <- D_bar + p.d

performance_metric <- data.frame(mean_lg_lkhd = mean(output$lg_lkhd), 
                                 lglike_prameter_bar = lg_lkhd_thbar,
                                 D_bar = D_bar,
                                 D_thbar = D_thbar,
                                 p.d = p.d,
                                 DIC = DIC)

file_name <- paste0(substr(count_var_name, start = 3, stop = nchar(count_var_name)))

all_output_list <- list(sample_parameter = output,
                        parameter_summary = parameter_summary,
                        morans_I_summary = morans_I_summary,
                        performance_metric = performance_metric)

save(all_output_list, file = paste0("./output/", file_name, ".RData"))

rm(output, parameter_name,parameter_summary, morans_I_summary, D_bar,D_thbar, p.d, DIC,
   performance_metric, file_name, all_output_list)

# ###### Log Likelihood of all data
# beta_bar = matrix(parameter_summary$mean[grepl("beta", parameter_summary$parameter)], ncol = 1)
# theta_bar = matrix(parameter_summary$mean[grepl("theta", parameter_summary$parameter)], ncol = 1)
# fitted_bar <- exp(X %*% beta_bar + theta_bar)
# lg_lkhd_alldata <- sum(dpois(count_var, lambda = t(fitted_bar), log = T))
# 
# performance_metric <- data.frame(lglike_prameter_bar = lg_lkhd_thbar_train,
#                                  D_bar = D_bar_train,
#                                  D_theta_bar = D_thbar_train,
#                                  p.d = p.d,
#                                  DIC = DIC)
# 
# # count_var_name <- colnames(rawData$count_data)[current_response]
# file_name <- paste0(substr(count_var_name, start = 3, stop = nchar(count_var_name)), "_fold_", fold_number)
# all_output_list <- list(test_locator = test_locator, tested_count = count_var[test_locator], sample_parameter = output,
#                         parameter_summary = parameter_summary,
#                         prediction_summary = prediction_summary,
#                         morans_I_summary = morans_I_summary,
#                         performance_metric = performance_metric,
#                         loglike_alldata_mparam = lg_lkhd_alldata)
# 
# # assign(list_name, all_output_list)
# save(all_output_list, file = paste0("./output_file/", file_name, ".RData"))
# # return(all_output_list)
# 
# 
# rm(output, epsilon_test, lambda_test, predicted_count, prediction_summary, parameter_name,
#    parameter_summary, morans_I_summary, morani_iter, lg_lkhd_thbar_train, D_bar_train,D_thbar_train, p.d, DIC,
#    beta_bar, theta_bar, fitted_bar, lg_lkhd_alldata, performance_metric, file_name, all_output_list)
# 
# # rm(tempM)
