rm(list = ls())
# setwd
# Load the required libraries
library(MCMCpack)
library(tidyverse)
library(ggplot2)
# install.packages("DAAG")
library(DAAG)
source("functions_gibbs_ex2.R")

# univariate data: for gibbs_p1 to gibbs_p4
grouped_dataset <- get(load("./data/grouped_dataset.Rdata"))

#Regression data: mv_gibbs_p1 to mv_gibbs_p4
traindf <- get(load("./data/ataset_regression_train.Rdata"))

#---------------------------------------------
# user defined parameter value for sampling:
#---------------------------------------------

#-------------------------------------------------------------------------------
#                       for univariate normal (grouped)
#-------------------------------------------------------------------------------
num_iter = 10000; burn_in = 5000
# v0 is a vector for initial mu which has length 3
v0 = rep(0, 3)
initial_sigma2 = 5
# alpha0, beta0 or a0, b0 are given parameter for sigma^2 prior 
alpha0 = 2.1; beta0 = 2.5
a0 = 0; b0 = 4

# for univariate grouped normal
normal_posterior_sample_1 <- gibbs_p1(y = grouped_dataset, num_iter = num_iter, burn_in = burn_in, v0 = v0, alpha0 = alpha0, beta0 = beta0)
normal_posterior_sample_2 <- gibbs_p2(y = grouped_dataset, num_iter = num_iter, burn_in = burn_in, v0 = v0, alpha0 = alpha0, beta0 = beta0)
normal_posterior_sample_3 <- gibbs_p3(y = grouped_dataset, num_iter = num_iter, burn_in = burn_in, v0 = v0)
normal_posterior_sample_4 <- gibbs_p4(y = grouped_dataset, num_iter = num_iter, burn_in = burn_in, v0 = v0, initial_sigma2 = initial_sigma2, a0 = a0, b0 = b0)

#-----------------------------
# Summary and visualization
#-----------------------------
sample_summary <- normal_posterior_sample_4

medians <- apply(sample_summary, 2, median)
quantiles <- apply(sample_summary, 2, quantile, probs = c(.025, 0.975))

# True values
Mu1 = 2.54; Mu2 =  6.2; Mu3 = 1.58; Sigma2 = 1.42 

graphics.off()
par(mfrow = c(2, 4))
plot(1:nrow(sample_summary), sample_summary$mu1 , "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste("Posterior of ", mu[1])))
abline(h = Mu1, col="red3")
plot(1:nrow(sample_summary), sample_summary$mu2, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste("Posterior of ", mu[2])))
abline(h = Mu2, col="red3")
plot(1:nrow(sample_summary), sample_summary$mu3, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste("Posterior of ", mu[3])))
abline(h = Mu3, col="red3")
plot(1:nrow(sample_summary), sample_summary$sigma2, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste("Posterior of ", sigma^2)))
abline(h = Sigma2, col="red3")


hist(sample_summary$mu1, freq = FALSE, xlab = expression(mu[1]), main = expression(paste("Histogram of ", mu[1])))
abline(v = c(Mu1, medians["mu1"], quantiles[, "mu1"]),  col = c("red3", "darkgreen", "navy", "navy"), lty = c(1,2, 4, 4), lwd = c(3, 3))
hist(sample_summary$mu2, freq = FALSE, xlab = expression(mu[2]), main = expression(paste("Histogram of ", mu[2])))
abline(v = c(Mu2, medians["mu2"], quantiles[, "mu2"]),  col = c("red3", "darkgreen", "navy", "navy"), lty = c(1,2, 4, 4), lwd = c(3, 3))
hist(sample_summary$mu3, freq = FALSE, xlab = expression(mu[3]), main = expression(paste("Histogram of ", mu[3])))
abline(v = c(Mu3, medians["mu3"], quantiles[, "mu3"]),  col = c("red3", "darkgreen", "navy", "navy"), lty = c(1,2, 4, 4), lwd = c(3, 3))
hist(sample_summary$sigma2, freq = FALSE, xlab = expression(sigma^2), main = expression(paste("Histogram of ", sigma^2)))
abline(v = c(Sigma2, medians["sigma2"], quantiles[, "sigma2"]),  col = c("red3", "darkgreen", "navy", "navy"), lty = c(1,2, 4, 4), lwd = c(3, 3))


#--------------------------------------------------------------------------
#                               for regression
#--------------------------------------------------------------------------
num_iter = 50000; burn_in = 10000; thin_in = 5
# for univariate data betas is a vector for initial beta which has length 4
betas = c(5.51, 0, 0, 0)
initial_sigma2 = 5
alpha0 = 2.1; beta0 = 2.5
a0 = 0; b0 = 4

reg_parameter_sample_1 <- mv_gibbs_p1(traindf = traindf, num_iter = num_iter, burn_in = burn_in, thin_in = thin_in, betas = betas, alpha0 = alpha0, beta0 = beta0)
reg_parameter_sample_2 <- mv_gibbs_p2(traindf = traindf, num_iter = num_iter, burn_in = burn_in, thin_in = thin_in, betas = betas, alpha0 = alpha0, beta0 = beta0)
reg_parameter_sample_3 <- mv_gibbs_p3(traindf = traindf, num_iter = num_iter, burn_in = burn_in, thin_in = thin_in, betas = betas)
reg_parameter_sample_4 <- mv_gibbs_p4(traindf = traindf, num_iter = num_iter, burn_in = burn_in, thin_in = thin_in, betas = betas, initial_sigma2 = initial_sigma2, a0 = a0, b0 = b0)

#-----------------------------
# Summary and visualization
#-----------------------------
sample_summary <- reg_parameter_sample_4

medians <- apply(sample_summary, 2, median)
quantiles <- apply(sample_summary, 2, quantile, probs = c(.025, 0.975))

# True values
Beta0 = 2.1; Beta1 = 1.25; Beta2 = 3.7; Beta3 = -0.5; Sigma2 =  2.8^2.0

graphics.off()
par(mfrow = c(2, 5))
plot(1:nrow(sample_summary), sample_summary$beta0 , "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste("Posterior of ", beta[0])))
abline(h = Beta0, col="red3")
plot(1:nrow(sample_summary), sample_summary$beta1, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste("Posterior of ", beta[1])))
abline(h = Beta1, col="red3")
plot(1:nrow(sample_summary), sample_summary$beta2, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste("Posterior of ", beta[2])))
abline(h = Beta2, col="red3")
plot(1:nrow(sample_summary), sample_summary$beta3, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste("Posterior of ", beta[3])))
abline(h = Beta3, col="red3")
plot(1:nrow(sample_summary), sample_summary$sigma2, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste("Posterior of ", sigma^2)))
abline(h = Sigma2, col="red3")

hist(sample_summary$beta0, freq = FALSE, xlab = expression(beta[0]), main = expression(paste("Histogram of ", beta[0])))
abline(v = c(Beta0, medians["beta0"], quantiles[, "beta0"]),  col = c("red3", "darkgreen", "navy", "navy"), lty = c(1,2, 4, 4), lwd = c(3, 3))
hist(sample_summary$beta1, freq = FALSE, xlab = expression(beta[1]), main = expression(paste("Histogram of ", beta[1])))
abline(v = c(Beta1, medians["beta1"], quantiles[, "beta1"]),  col = c("red3", "darkgreen", "navy", "navy"), lty = c(1,2, 4, 4), lwd = c(3, 3))
hist(sample_summary$beta2, freq = FALSE, xlab = expression(beta[2]), main = expression(paste("Histogram of ", beta[2])))
abline(v = c(Beta2, medians["beta2"], quantiles[, "beta2"]),  col = c("red3", "darkgreen", "navy", "navy"), lty = c(1,2, 4, 4), lwd = c(3, 3))
hist(sample_summary$beta3, freq = FALSE, xlab = expression(beta[3]), main = expression(paste("Histogram of ", beta[3])))
abline(v = c(Beta3, medians["beta3"], quantiles[, "beta3"]),  col = c("red3", "darkgreen", "navy", "navy"), lty = c(1,2, 4, 4), lwd = c(3, 3))
hist(sample_summary$sigma2, freq = FALSE, xlab = expression(sigma^2), main = expression(paste("Histogram of ", sigma^2)))
abline(v = c(Sigma2, medians["sigma2"], quantiles[, "sigma2"]),  col = c("red3", "darkgreen", "navy", "navy"), lty = c(1,2, 4, 4), lwd = c(3, 3))


#----------------------------------
# likelihood Calculation
#----------------------------------
normal_log_likelihood <- function(data, mu1, mu2, mu3, sigma2){
  # data = grouped_dataset; mu1 = sample_summary$mu1; mu2 = sample_summary$mu2
  # mu3 = sample_summary$mu3; sigma2 = sample_summary$sigma2
  data1 <- data[[1]]; data2 <- data[[2]]; data3 <- data[[3]]
  n1 <- length(data1); n2 <- length(data2); n3 <- length(data3)
  parapmeter_dim <- length(mu1)
  tempdf <- data.frame()
  
  for(i in 1: parapmeter_dim) {
    i <- 1
    L1 <- -(n/2) * log(2*pi*sigma2[i]) - (2 * sigma2[i])^-1 * sum((data1 - mu1[i])^2)
    L2 <- -(n/2) * log(2*pi*sigma2[i]) - (2 * sigma2[i])^-1 * sum((data2 - mu2[i])^2)
    L3 <- -(n/2) * log(2*pi*sigma2[i]) - (2 * sigma2[i])^-1 * sum((data3 - mu3[i])^2)
    # likelihood <- rbind(likelihood, L)
    tempdf <- rbind(tempdf, c(mu1[i], mu2[i], mu3[i], sigma2[i], L1, L2, L3))
  } 
  output <- tempdf %>% 
    set_names(c("mu1", "mu2", "mu3", "sigma2", "L1", "L2", "L3"))
  
  return(output)
  
}

joint_log_likelihood <- function(data, mu1, mu2, mu3, sigma2){
  # data = grouped_dataset; mu1 = sample_summary$mu1; mu2 = sample_summary$mu2
  # mu3 = sample_summary$mu3; sigma2 = sample_summary$sigma2
  
  data1 <- data[[1]]; data2 <- data[[2]]; data3 <- data[[3]]
  n1 <- length(data1); n2 <- length(data2); n3 <- length(data3)
  parameter_dim <- length(mu1)
  tempdf <- data.frame()
  
  for(i in 1: parameter_dim) { # parameter_dim
    # i <- 1
    L1 <- -(n1/2) * log(2*pi*sigma2[i]) - (2 * sigma2[i])^-1 * sum((data1 - mu1[i])^2)
    L2 <- -(n2/2) * log(2*pi*sigma2[i]) - (2 * sigma2[i])^-1 * sum((data2 - mu2[i])^2)
    L3 <- -(n3/2) * log(2*pi*sigma2[i]) - (2 * sigma2[i])^-1 * sum((data3 - mu3[i])^2)
    JL <- L1 + L2 + L3
    # likelihood <- rbind(likelihood, L)
    tempdf <- rbind(tempdf, c(mu1[i], mu2[i], mu3[i], sigma2[i], L1, L2, L3, JL))
  } 
  output <- tempdf %>% 
    set_names(c("mu1", "mu2", "mu3", "sigma2", "L1", "L2", "L3", "JL"))
  
  return(output)
  
}

gn_likelihood_p1 <- joint_log_likelihood(data = grouped_dataset, mu1 = normal_posterior_sample_1$mu1, mu2 = normal_posterior_sample_1$mu2,
                              mu3 = normal_posterior_sample_1$mu3, sigma2 = normal_posterior_sample_1$sigma2)
gn_likelihood_p2 <- joint_log_likelihood(data = grouped_dataset, mu1 = normal_posterior_sample_2$mu1, mu2 = normal_posterior_sample_2$mu2,
                                          mu3 = normal_posterior_sample_2$mu3, sigma2 = normal_posterior_sample_2$sigma2)
gn_likelihood_p3 <- joint_log_likelihood(data = grouped_dataset, mu1 = normal_posterior_sample_3$mu1, mu2 = normal_posterior_sample_3$mu2,
                                          mu3 = normal_posterior_sample_3$mu3, sigma2 = normal_posterior_sample_3$sigma2)
gn_likelihood_p4 <- joint_log_likelihood(data = grouped_dataset, mu1 = normal_posterior_sample_4$mu1, mu2 = normal_posterior_sample_4$mu2,
                                          mu3 = normal_posterior_sample_4$mu3, sigma2 = normal_posterior_sample_4$sigma2)

plot(1:nrow(gn_likelihood_p1), gn_likelihood_p1$JL, "l", xlab = "Iteration No", ylab = "Joint Log-likelihood", main = "Log likelihood plot")


gn_likelihood_p1[which(gn_likelihood_p1$JL == min(gn_likelihood_p1$JL)),]
gn_likelihood_p1[which(gn_likelihood_p1$L1 == min(gn_likelihood_p1$L1)),]
apply(normal_posterior_sample_1, 2, median)
