rm(list = ls())

# # un-comment the rm line below when look for suitable sd choices
# rm(list=ls()[! ls() %in% c("mh_sd_choice","mh_sd", "sddf", "count_col", "response", "k")])
# #===========================================================================================

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

num_iter = 20000; burn_in = 8000; thin_in = 30

# # to look for suitable sd choice in M-H algorithm, comment below two line===
print(colSums(rawData$count_data!= 0, na.rm = T))
response <- as.numeric(readline("select response column number (1 - 10): "))
# # response <- 5 # un-comment if need to run loop for a single bird choice
# #===========================================================================

current_response <- response
count_var <- rawData$count_data[, current_response] # hyper-parameter
count0_locator <- which(count_var == 0)
count_var_name <- colnames(rawData$count_data)[current_response]

##### setting up initial theta
theta_init <- matrix(0, nrow = n, ncol = 1)

###### setting up initial log_lambda
log_lambda_init <- log(count_var +.25)

###### setting up initial beta
beta_init <- unname(lm(log_lambda_init ~ covariate_data)$coef) # BETA2 initial value

# initial for sigma2, tau2_
sigma2_init = 1; tau2_init = 1

###### prior on sigma2 and tau2_theta1, tau2_theta2
a_tau_pr = 2.01 ; b_tau_pr = 1.01
a_sigma_pr = 2.01; b_sigma_pr = 1.01

#-------------------- Start: poisson count model

set.seed(729)

####### prior info:
# prior information for all beta1 and beta2 vector
mu_beta = matrix(c(2, 2, 1, -1, 0, 2), ncol = 1); sigma_beta = 10^4 * diag(p + 1)

###### Assigning initial values 
beta = matrix(beta_init, ncol = 1)
theta = theta_init; tau2 = tau2_init; sigma2 = sigma2_init
rm(beta_init, theta_init, tau2_init, sigma2_init)

log_lambda = matrix(log_lambda_init, ncol = 1); current_lambda = exp(log_lambda)

# iteration has no effect on estimates of tau2, A_sigma2 and sigma_beta_inv, 
# so calculating outside of the loop
A_tau2 = ((n - eigen0_count)/2) + a_tau_pr    # posterior scale parameter for tau2
A_sigma2 = (n/2) + a_sigma_pr                 # posterior scale parameter for sigma2
sigma_beta_inv = chol2inv(chol(sigma_beta))

#output matrix
# c(lg_lkhd, acc_rate_per_iteration, beta1, beta2, sigma2, tau2_theta1, tau2_theta2, theta1, theta2, log_lambda, zee)
ncol_tempM = length(c(1, 1, beta, sigma2, tau2, theta, log_lambda, log_lambda, log_lambda)) 
tempM <- matrix(0, nrow = num_iter, ncol = ncol_tempM)

mh_sd <- as.numeric(readline("Type preferred sd value for log_lambda M-H: "))

pb <- txtProgressBar(min = 1, max = num_iter, style = 3)

for (j in 1:num_iter) {
  
  # tic = Sys.time()
  
  ###### update log_lambda: Metropolis-Hastings
  # propose a log_lambda value from Normal
  proposed_log_lambda <- rnorm(n, mean = log_lambda, sd = mh_sd) # sd can be changed mh_sd
  proposed_lambda = exp(proposed_log_lambda)
  mean_log_lambda = as.vector(X %*% beta + theta)
  
  # ratio
  l_posterior_proposal =  -(2 * sigma2)^-1 * (proposed_log_lambda - mean_log_lambda)^2 - proposed_lambda + count_var* log(proposed_lambda)
  l_posterior_current = - (2 * sigma2)^-1 * (log_lambda - mean_log_lambda)^2 - current_lambda + count_var * log(current_lambda)
  l_target_by_proposal = l_posterior_proposal - l_posterior_current # because of taking log
  
  # accepting proposal
  l_uniform_vector = log(runif(n, min = 0, max = 1))
  gain <- l_uniform_vector < l_target_by_proposal
  log_lambda[gain] <- proposed_log_lambda[gain] # updated log_lambda 
  current_lambda[gain] <- proposed_lambda[gain] # required for next iteration
  lambda_acc <- as.numeric(gain) # vector of log_lambda acceptance (binary)
  acc_rate_per_iteration <- mean(lambda_acc)  # lambda_acc is binary
  
  ###### Update theta: 
  res_nonspatial <-  log_lambda - X%*% beta # (log_lambda- Xb) part in theta posterior
  sum_sigma_w <- sigma2 * D_w_vector
  theta_post_var = (sigma2 * tau2)/(tau2 + sum_sigma_w)
  
  for(i in 1: n) {
    theta_post_mean =  (tau2 * res_nonspatial[i] + sigma2 * sum(theta[ W_list[[i]] ]))/(tau2 + sum_sigma_w[i])
    theta_post <- rnorm(1, theta_post_mean, sqrt(theta_post_var[i]))
    theta[i] <- theta_post

    rm(theta_post_mean, theta_post)
  }
  # theta <- theta - mean(theta)
  theta[island_locator] <- theta[island_locator] - mean(theta[island_locator])
  
  # Update tau2:
  # prior on tau2 ~ IG( 2.5, 2.5)
  # posterior shape of tau2 is fixed and calculated before iteration loop
  # A_tau2 = ((n - eigen0_count)/2) + a_tau_pr
  # B_tau2 = (t(theta) %*% (D_w_tt - W_tt) %*% theta)/2 + b_tau_pr
  
  B_tau2 = .5 * (sum(D_w_vector * theta^2) - sum(theta * sapply(W_list, function(x) sum(theta[x])))) + b_tau_pr
  tau2 <- 1/rgamma(1, shape = A_tau2, rate = B_tau2)
  
  ##### Update sigma2:
  # prior on sigma2 ~ IG( 2.5, 2.5)
  # A_sigma2 = (n/2) + a_sigma_pr  # calculated before loop
  ss <- sum((res_nonspatial - theta)^2) # res_nonspatial = log_lambda - X*beta
  B_sigma2 = ss/2 + b_sigma_pr
  sigma2 <- 1/rgamma(1, shape = A_sigma2, rate = B_sigma2)
  
  ##### Update beta:
  # sigma_beta_inv = chol2inv(chol(sigma_beta)) calculated before loop
  # x_transpose_sigma_inv = tX %*% ((sigma2^-1) * diag(n))
  x_transpose_sigma_inv = (sigma2^-1) * tX
  A_beta = x_transpose_sigma_inv %*% X + sigma_beta_inv; A_beta_inv = chol2inv(chol(A_beta))
  B_beta = x_transpose_sigma_inv %*% (log_lambda - theta) + sigma_beta_inv %*% mu_beta
  beta <- mvrnorm(1, A_beta_inv %*% B_beta, A_beta_inv)
  
  
  ##### Log-likelihood
  lg_lkhd_vector <- dpois(count_var, lambda = exp(log_lambda), log = T)
  lg_lkhd <- sum(dpois(count_var, lambda = exp(log_lambda), log = T))
  
  # Accumulates all parameters for a single iteration
  tempM[j, ] <-  c(lg_lkhd, acc_rate_per_iteration, beta, sigma2, tau2, theta, log_lambda, lambda_acc, lg_lkhd_vector)
  
  rm(B_tau2, B_sigma2, x_transpose_sigma_inv,
     A_beta, B_beta, A_beta_inv, 
     acc_rate_per_iteration, lambda_acc, lg_lkhd)
  
  setTxtProgressBar(pb, j) # to see the progress 
} 

output <- as.data.frame(tempM) %>%
  .[(burn_in+1):num_iter, ] %>%
  .[seq(from = 1, to = num_iter - burn_in, by = thin_in), ] %>%
  set_names(c("lg_lkhd", "acc_rate_iter",  paste0("beta_", 0:p), "sigma2", "tau2",  paste0("theta_", 1:n), paste0("lg_lambda", 1:n), paste0("lambda", 1:n, "_acc"), paste0("log_like_", 1:n)))

# # Comment below two lines when run usual====
# acc <- round(mean(output$acc_rate_iter) *100)
# if(acc < 45) next
# #===========================================

geweke_stat <- function(x){
  geweke_z = unname(coda::geweke.diag(x)$z)
  p_val = 2*(1-pnorm(abs(geweke_z)))
  # converge_y = if_else(p_val < .05, "y", "n")
  return(c(geweke_z = geweke_z, p_val = p_val))
}

parameter_name <- c( "lg_lkhd",  paste0("beta_", 0:p), "sigma2", "tau2",  paste0("theta_", 1:n), paste0("lg_lambda", 1:n) )
parameter_summary <- data.frame(mean = apply(output[, parameter_name], 2, mean),
                                median = apply(output[, parameter_name], 2, median),
                                q_05 = apply(output[, parameter_name], 2, quantile, .05),
                                q_95 = apply(output[, parameter_name], 2, quantile, .95)) %>% 
  bind_cols(as.data.frame(t((apply(output[, parameter_name], 2, geweke_stat))))) %>% 
  mutate(converge_y = if_else(p_val >.05, "y", "n")) %>% 
  tibble:: rownames_to_column("parameter") %>%
  # bind_cols(apply(output[, parameter_name], 2, geweke_stat)) %>% 
  dplyr:: select(parameter, everything())

# ###### Morans I

### For theta  vector
  theta_morani_iter <- apply(output[, grepl(paste0("theta_"), names(output))], 1, function(theta) ape::Moran.I(theta, W, scaled = FALSE)$observed)
  morans_I_summary <- data.frame(model_name = "Poisson", 
                                 theta = paste0("theta"), 
                                 MI_thbar = ape:: Moran.I(apply(output[, grepl(paste0("theta_"), names(output))], 2, mean), W, scaled = F)$observed, 
                                 MI_mean = mean(theta_morani_iter), 
                                 MI_median = median(theta_morani_iter), 
                                 MI_ql = quantile(theta_morani_iter, .05), 
                                 MI_qh = quantile(theta_morani_iter, .95)) %>% 
    tibble::remove_rownames()


# print(morans_I_summary)

# ###### DIC

lg_lkhd_thbar <-sum(dpois(count_var, lambda = exp(apply(output[, paste0("lg_lambda", 1:n)], 2, mean)), log = T))

D_bar <- -2 * mean(output$lg_lkhd)
D_thbar <- -2 * lg_lkhd_thbar
p.d <- D_bar - D_thbar
DIC <- D_bar + p.d

performance_metric <- data.frame(model_name = "Poisson", 
                                 mean_acc_rate_lambda = mean(output$acc_rate_iter),
                                 mean_lg_lkhd = mean(output$lg_lkhd), 
                                 lglike_prameter_bar = lg_lkhd_thbar,
                                 D_bar = D_bar,
                                 D_thbar = D_thbar,
                                 p.d = p.d,
                                 DIC = DIC)


# comment first line and uncomment 2nd line for file name when run usual==
# file_name <- paste0(substr(count_var_name, start = 3, stop = nchar(count_var_name)), "_poisson_sd_", mh_sd, "acc_", acc)
file_name <- paste0(substr(count_var_name, start = 3, stop = nchar(count_var_name)), "_poisson")
#=========================================================================

all_output_list <- list(sample_parameter = output,
                        parameter_summary = parameter_summary,
                        morans_I_summary = morans_I_summary,
                        performance_metric = performance_metric)

save(all_output_list, file = paste0("./output/", file_name, ".RData"))



# un-comment below two line when not looking for suitable sd in M-H ======
rm(output, parameter_name,parameter_summary, morans_I_summary, D_bar,D_thbar, p.d, DIC,
   performance_metric, file_name, all_output_list)
#=========================================================================
