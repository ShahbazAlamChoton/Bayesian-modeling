rm(list = ls())
graphics.off()
# setwd
library(MCMCpack)
library(tidyverse)

traindf <- get(load("./data/dataset_regression_train.Rdata"))



mv_gibbs_p1 <- function(traindf, num_iter, burn_in, thin_in, betas, alpha0, beta0) {    
  
  set.seed(729)
  # parameters
  y = matrix(traindf[, 1], ncol = 1)
  n = nrow(traindf)
  X = matrix(c(rep(1, n), traindf[, 2:4]), ncol = 4, byrow = F)
  
  beta_init <- matrix(betas, ncol = 1)
  sst0 = t(y - X %*% beta_init) %*% (y - X %*% beta_init)
  
  mu_beta = matrix(c(20, -10, 10, 15), ncol = 1)
  sigma_beta = 10^4 * matrix(c(0.90, 0.59, 0.93, 0.97,
                        0.59, 1.59, 1.12, 0.75,
                        0.93, 1.12, 1.26, 1.00,
                        0.97, 0.75, 1.00, 1.10), nrow = 4, byrow = TRUE)
  

  tempDf <- data.frame()
  # tempDf <- data.frame(m = numeric(), v = numeric())
  
  for (j in 1:num_iter) {
    
    #Update sigma2
    # what we get from full conditional of sigma2
    # prior on sigma2 ~ IG(2.1, 2.5)
    post_alpha = n/2 + alpha0
    post_beta = beta0 + 0.5 * sst0
    sigma2 <- 1/rgamma(1, shape = post_alpha, rate = post_beta)
    
    
    # Update mu
    # what we get from full conditional of mu
    # prior on beta ~ MVN( mu_beta, sigma_beta)
    
    A = sigma2^-1 * t(X) %*% X + solve(sigma_beta)
    B = sigma2^-1 * t(X) %*% y + solve(sigma_beta) %*% mu_beta
    beta <- mvrnorm(1, solve(A) %*% B, solve(A))
    
    # updating sst0 for next draw of sigma2 
    sst0 = t(y - X %*% beta) %*% (y - X %*% beta)
    
    tempDf <-  rbind(tempDf, c(beta, sigma2))
  }
  
  # thin <- seq(from = 1, to = num_iter - burn_in, by = thin_in)
  output <- tempDf %>% 
    .[(burn_in+1):num_iter, ] %>%
    .[seq(from = 1, to = num_iter - burn_in, by = thin_in), ] %>% 
    set_names(c("beta0", "beta1", "beta2", "beta3", "sigma2"))
  
  return(output)
}

mv_gibbs_p2 <- function(traindf, num_iter, burn_in, thin_in, betas, alpha0, beta0) {    
  
  set.seed(729)
  # parameters
  y = matrix(traindf[, 1], ncol = 1)
  n = nrow(traindf)
  X = matrix(c(rep(1, n), traindf[, 2:4]), ncol = 4, byrow = F)
  
  beta_init <- matrix(betas, ncol = 1)
  
  mu_beta = matrix(c(20, -10, 10, 15), ncol = 1)
  sigma_beta = 10^4 * matrix(c(0.90, 0.59, 0.93, 0.97,
                               0.59, 1.59, 1.12, 0.75,
                               0.93, 1.12, 1.26, 1.00,
                               0.97, 0.75, 1.00, 1.10), nrow = 4, byrow = TRUE)
  
  sst0 = t(y - X %*% beta_init) %*% (y - X %*% beta_init) + t(beta_init - betas) %*% solve(sigma_beta) %*% (beta_init - betas)
  
  tempDf <- data.frame()
  # tempDf <- data.frame(m = numeric(), v = numeric())
  
  for (j in 1:num_iter) {
    
    #Update sigma2
    # what we get from full conditional of sigma2
    # prior on sigma2 ~ IG(2.1, 2.5)
    post_alpha = n/2 + 2 + alpha0
    post_beta = beta0 + 0.5 * sst0
    sigma2 <- 1/rgamma(1, shape = post_alpha, rate = post_beta)
    
    
    # Update mu
    # what we get from full conditional of mu
    # prior on beta ~ MVN( mu_beta, sigma2 * sigma_beta)
    A = sigma2^-1 * t(X) %*% X + solve(sigma2 * sigma_beta)
    B = sigma2^-1 * t(X) %*% y + solve(sigma2 * sigma_beta) %*% mu_beta
    beta <- mvrnorm(1, solve(A) %*% B, solve(A))
    
    # updating sst0 for next draw of sigma2 
    sst0 = t(y - X %*% beta) %*% (y - X %*% beta)
    
    tempDf <-  rbind(tempDf, c(beta, sigma2))
  }
  
  # thin <- seq(from = 1, to = num_iter - burn_in, by = thin_in)
  output <- tempDf %>% 
    .[(burn_in+1):num_iter, ] %>%
    .[seq(from = 1, to = num_iter - burn_in, by = thin_in), ] %>% 
    set_names(c("beta0", "beta1", "beta2", "beta3", "sigma2"))
  
  return(output)
}

mv_gibbs_p3 <- function(traindf, num_iter, burn_in, thin_in, betas) {    
  
  set.seed(729)
  # parameters
  y = matrix(traindf[, 1], ncol = 1)
  n = nrow(traindf)
  X = matrix(c(rep(1, n), traindf[, 2:4]), ncol = 4, byrow = F)
  
  beta_init <- matrix(betas, ncol = 1)
  sst0 = t(y - X %*% beta_init) %*% (y - X %*% beta_init)
  
  mu_beta = matrix(c(20, -10, 10, 15), ncol = 1)
  sigma_beta = 10^4 * matrix(c(0.90, 0.59, 0.93, 0.97,
                               0.59, 1.59, 1.12, 0.75,
                               0.93, 1.12, 1.26, 1.00,
                               0.97, 0.75, 1.00, 1.10), nrow = 4, byrow = TRUE)
  
  
  tempDf <- data.frame()
  # tempDf <- data.frame(m = numeric(), v = numeric())
  
  for (j in 1:num_iter) {
    
    #Update sigma2
    # what we get from full conditional of sigma2
    # prior on sigma2 ~ IG(2.1, 2.5)
    post_alpha = n/2
    post_beta = 0.5 * sst0
    sigma2 <- 1/rgamma(1, shape = post_alpha, rate = post_beta)
    
    
    # Update mu
    # what we get from full conditional of mu
    # prior on beta ~ MVN( mu_beta, sigma_beta)
    
    A_inverse_b = solve(t(X) %*% X) %*% t(X) %*% y
    A_inverse = sigma2 * solve(t(X) %*% X) 
    beta <- mvrnorm(1, A_inverse_b, A_inverse)
    
    # updating sst0 for next draw of sigma2 
    sst0 = t(y - X %*% beta) %*% (y - X %*% beta)
    
    tempDf <-  rbind(tempDf, c(beta, sigma2))
  }
  
  # thin <- seq(from = 1, to = num_iter - burn_in, by = thin_in)
  output <- tempDf %>% 
    .[(burn_in+1):num_iter, ] %>%
    .[seq(from = 1, to = num_iter - burn_in, by = thin_in), ] %>% 
    set_names(c("beta0", "beta1", "beta2", "beta3", "sigma2"))
  
  return(output)
}

mv_gibbs_p4 <- function(traindf, num_iter, burn_in, thin_in, betas, initial_sigma2, a0, b0) {    
  
  set.seed(729)
  # parameters
  y = matrix(traindf[, 1], ncol = 1)
  n = nrow(traindf)
  X = matrix(c(rep(1, n), traindf[, 2:4]), ncol = 4, byrow = F)
  
  # Assigning initial values of betas and its predefined prior
  beta_init <- matrix(betas, ncol = 1)
  sst0 = t(y - X %*% beta_init) %*% (y - X %*% beta_init)
  
  mu_beta = matrix(c(20, -10, 10, 15), ncol = 1) #given value
  sigma_beta = 10^4 * matrix(c(0.90, 0.59, 0.93, 0.97,
                               0.59, 1.59, 1.12, 0.75,
                               0.93, 1.12, 1.26, 1.00,
                               0.97, 0.75, 1.00, 1.10), nrow = 4, byrow = TRUE)
  
  sigma2 <- initial_sigma2
  accepted_proposal <- 0
  tempDf <- data.frame()
  # tempDf <- data.frame(m = numeric(), v = numeric())
  
  for (j in 1:num_iter) {
    
    # update sigma2 using Metropolis-Hastings
    
    # propose a sigma^2 value from Inverse Gamma density proposal
    proposed_sigma2 <- 1/rgamma(1, shape = n/2, rate = sst0/2)
    
    # for inverse gamma as a proposal
    target_by_proposal <- (2 * b0)^-1 * ((log(proposed_sigma2 - a0))^2 - (log(sigma2 - a0))^2)
    
    acceptance_ratio <- min(1, target_by_proposal) # markov kernel cancel out
    # print(paste0("Acc ratio: ", acceptance_ratio," sigma2 ",  sigma2, " Proposal: ", proposed_sigma2))
    # pause()
    
      if(log(runif(1)) < acceptance_ratio){
        
        # accept proposed theta
        sigma2 <- proposed_sigma2
        accepted_proposal <- accepted_proposal + 1
      } 
    
    # Update mu
    # what we get from full conditional of mu
    # prior on beta ~ MVN( mu_beta, sigma_beta)
    
    A = sigma2^-1 * t(X) %*% X + solve(sigma_beta)
    B = sigma2^-1 * t(X) %*% y + solve(sigma_beta) %*% mu_beta
    beta <- mvrnorm(1, solve(A) %*% B, solve(A))
    
    # updating sst0 for next draw of sigma2 
    sst0 = t(y - X %*% beta) %*% (y - X %*% beta)
    
    tempDf <-  rbind(tempDf, c(beta, sigma2))
  }
  
  # thin <- seq(from = 1, to = num_iter - burn_in, by = thin_in)
  output <- tempDf %>% 
    .[(burn_in+1):num_iter, ] %>%
    .[seq(from = 1, to = num_iter - burn_in, by = thin_in), ] %>% 
    set_names(c("beta0", "beta1", "beta2", "beta3", "sigma2"))
  
  acceptance_rate <- accepted_proposal/num_iter
  cat(" Acceptance rate using Gamma proposal: ", round(acceptance_rate, 3))
  return(output)
}

parameter_sample_1 <- mv_gibbs_p1(traindf = traindf, num_iter = 50000, burn_in = 10000, thin_in = 5, betas = c(5.51,0, 0, 0), alpha0 = 2.1, beta0 = 2.5)
parameter_sample_2 <- mv_gibbs_p2(traindf = traindf, num_iter = 50000, burn_in = 10000, thin_in = 5, betas = c(1,1,1,1), alpha0 = 2.1, beta0 = 2.5)
parameter_sample_3 <- mv_gibbs_p3(traindf = traindf, num_iter = 50000, burn_in = 10000, thin_in = 5, betas = c(1,1,1,1))
parameter_sample_4 <- mv_gibbs_p4(traindf = traindf, num_iter = 50000, burn_in = 10000, thin_in = 5, betas = c(1,1,1,1), initial_sigma2 = 5, a0 = 0, b0 = 4)


# Summary and visualisation
sample_summary <- parameter_sample_4


apply(sample_summary, 2, median)
apply(sample_summary, 2, quantile, probs = c(.025, 0.975))


graphics.off()
par(mfrow = c(2, 3))
plot(1:8000, sample_summary$beta0 , "l", main = expression(paste("Posterior of ", beta[0])))
abline(h = 2.1, col="red3")
plot(1:8000, sample_summary$beta1, "l", main = expression(paste("Posterior of ", beta[1])))
abline(h = 1.25, col="red3")
plot(1:8000, sample_summary$beta2, "l", main = expression(paste("Posterior of ", beta[2])))
abline(h = 3.7, col="red3")
plot(1:8000, sample_summary$beta3, "l", main = expression(paste("Posterior of ", beta[3])))
abline(h = -0.5, col="red3")
plot(1:8000, sample_summary$sigma2, "l", main = expression(paste("Posterior of ", sigma^2)))
abline(h = 7.84, col="red3")

testdf <- get(load("./data/dataset_regression_test.Rdata"))
pred_y1 <- sample_summary$beta0 + sample_summary$beta1 * testdf[1, 2] + sample_summary$beta2 * testdf[1, 3] + sample_summary$beta3 * testdf[1, 4] + sqrt(sample_summary$sigma2) * rnorm(8000,0,1)
pred_y2 <- sample_summary$beta0 + sample_summary$beta1 * testdf[2, 2] + sample_summary$beta2 * testdf[2, 3] + sample_summary$beta3 * testdf[2, 4] + sqrt(sample_summary$sigma2) * rnorm(8000,0,1)

par(mfrow = c(1, 2))
hist(pred_y1, freq = FALSE, main = expression(paste("Predictive distribution of ", y[1])))
abline(v = c(1.137377, median(pred_y1)),  col = c("red3", "darkgreen"), lty = c(1,2), lwd = c(3, 3))
hist(pred_y2, freq = FALSE, main = expression(paste("Predictive distribution of ", y[2])))
abline(v = c(3.424860, median(pred_y2)),  col = c("red3", "darkgreen"), lty = c(1,2), lwd = c(3, 3))

y = matrix(traindf[, 1], ncol = 1)
n = nrow(traindf)
X = matrix(c(rep(1, n), traindf[, 2:4]), ncol = 4, byrow = F)

b_hat <- solve(t(X) %*% X) %*% (t(X) %*% y)
y1_new <- b_hat[1] + b_hat[2] * testdf[1, 2] + b_hat[3] * testdf[1, 3] + b_hat[4] * testdf[1, 4]
y2_new <- b_hat[1] + b_hat[2] * testdf[2, 2] + b_hat[3] * testdf[2, 3] + b_hat[4] * testdf[2, 4]

(sigma_estimate <- (t(y - X %*% b_hat) %*% (y - X %*% b_hat))/(n-3-1))
