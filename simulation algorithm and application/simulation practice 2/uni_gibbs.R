rm(list = ls())
# setwd
# Load the required libraries
library(MCMCpack)
library(tidyverse)
library(ggplot2)
# install.packages("DAAG")
library(DAAG)

# Load the dataset
grouped_dataset <- get(load("./data/grouped_dataset.Rdata"))

gibbs_p1 <- function(y, num_iter, burn_in, v0, alpha0, beta0){    
  
  set.seed(729)
  
  # parameters
  n1 = length(y[[1]]); n2 = length(y[[2]]); n3 = length(y[[3]])
  mu0_1 = v0[1]; mu0_2 = v0[2]; mu0_3 = v0[3] # initial values of mu_i to update sigma2
  y1bar = mean(y[[1]]); y2bar = mean(y[[2]]); y3bar = mean(y[[3]])
  sst0 = sum((y[[1]] - mu0_1)^2, (y[[2]] - mu0_2)^2, (y[[3]] - mu0_3)^2)
  
  # n = length(y)
  # nybar = sum(y)
  # sigma2 = var(y)
  
  # hyperparameters
  # mu0 = v0[1]
  # sigma20 = v0[2]
  
  tempDf <- data.frame()
  # tempDf <- data.frame(m = numeric(), v = numeric())
  
  for (j in 1:num_iter) {
    
    #Update sigma2
    # what we get from full conditional of sigma2
    # prior on sigma2 ~ IG(2.1, 2.5)
    post_alpha = sum(n1, n2, n3)/2 + alpha0
    post_beta = beta0 + 0.5 * sst0
    sigma2 <- 1/rgamma(1, shape = post_alpha, rate = post_beta)
    
    
    # Update mu
    # what we get from full conditional of mu
    # prior on mu ~ N(mu0 = 20, tau0^2 = 10^4)
    a11 = 1/10^4 + n1/sigma2; b11 = 20/10^4 + (n1*y1bar)/sigma2
    mu1 <- rnorm(1, mean = b11/a11, sd = sqrt(a11)^-1)
    
    a12 = 1/10^4 + n2/sigma2; b12 = 20/10^4 + (n2*y2bar)/sigma2
    mu2 <- rnorm(1, mean = b12/a12, sd = sqrt(a12)^-1)
    
    a13 = 1/10^4 + n3/sigma2; b13 = 20/10^4 + (n3*y3bar)/sigma2
    mu3 <- rnorm(1, mean = b13/a13, sd = sqrt(a13)^-1)
    
    # updating sst0 for next draw of sigma2 
    sst0 = sum((y[[1]] - mu1)^2, (y[[2]] - mu2)^2, (y[[3]] - mu3)^2)
    
    tempDf <-  rbind(tempDf, c(mu1, mu2, mu3, sigma2))
  }
  
  output <- tempDf %>% 
    .[(burn_in+1):num_iter, ] %>% 
    set_names(c("mu1", "mu2", "mu3", "sigma2"))
  
  return(output)
}

gibbs_p2 <- function(y, num_iter, burn_in, v0, alpha0, beta0){    
  
  set.seed(729)
  
  # parameters
  n1 = length(y[[1]]); n2 = length(y[[2]]); n3 = length(y[[3]])
  mu0_1 = v0[1]; mu0_2 = v0[2]; mu0_3 = v0[3] # initial values of mu_i to update sigma2
  y1bar = mean(y[[1]]); y2bar = mean(y[[2]]); y3bar = mean(y[[3]])
  sst0 = sum((y[[1]] - mu0_1)^2, (y[[2]] - mu0_2)^2, (y[[3]] - mu0_3)^2)
  ssmu0 <- sum((mu0_1 - 20)^2, (mu0_2 - 20)^2, (mu0_3 - 20)^2)
  
  tempDf <- data.frame()
  
  for (j in 1:num_iter) {
    
    #Update sigma2
    # what we get from full conditional of sigma2
    # prior on sigma2 ~ IG(2.1, 2.5)
    post_alpha = sum(n1, n2, n3, 3)/2 + alpha0
    post_beta = beta0 + 2^-1 * sst0 + 200^-1 * ssmu0
    sigma2 <- 1/rgamma(1, shape = post_alpha, rate = post_beta)
    
    
    # Update mu
    # what we get from full conditional of mu
    # prior on mu ~ N(mu0 = 20, tau0^2 = 10^4)
    a11 = (100 * sigma2)^-1 + n1/sigma2; b11 = 20 * (100 * sigma2)^-1 + (n1*y1bar)/sigma2
    mu1 <- rnorm(1, mean = b11/a11, sd = sqrt(a11)^-1)
    
    a12 = (100 * sigma2)^-1 + n2/sigma2; b12 = 20 * (100 * sigma2)^-1 + (n2*y2bar)/sigma2
    mu2 <- rnorm(1, mean = b12/a12, sd = sqrt(a12)^-1)
    
    a13 = (100 * sigma2)^-1 + n3/sigma2; b13 = 20 * (100 * sigma2)^-1 + (n3*y3bar)/sigma2
    mu3 <- rnorm(1, mean = b13/a13, sd = sqrt(a13)^-1)
    
    # updating sst0 for next draw of sigma2 
    sst0 = sum((y[[1]] - mu1)^2, (y[[2]] - mu2)^2, (y[[3]] - mu3)^2)
    ssmu0 = sum((mu1 - 20)^2, (mu2 - 20)^2, (mu3 - 20)^2)
    
    tempDf <-  rbind(tempDf, c(mu1, mu2, mu3, sigma2))
  }
  
  output <- tempDf %>% 
    .[(burn_in+1):num_iter, ] %>% 
    set_names(c("mu1", "mu2", "mu3", "sigma2"))
  
  return(output)
}

gibbs_p3 <- function(y, num_iter, burn_in, v0){    
  
  set.seed(729)
  
  # parameters
  n1 = length(y[[1]]); n2 = length(y[[2]]); n3 = length(y[[3]])
  mu0_1 = v0[1]; mu0_2 = v0[2]; mu0_3 = v0[3] # initial values of mu_i to update sigma2
  y1bar = mean(y[[1]]); y2bar = mean(y[[2]]); y3bar = mean(y[[3]])
  sst0 = sum((y[[1]] - mu0_1)^2, (y[[2]] - mu0_2)^2, (y[[3]] - mu0_3)^2)
  
  tempDf <- data.frame()
  # tempDf <- data.frame(m = numeric(), v = numeric())
  
  for (j in 1:num_iter) {
    
    #Update sigma2
    # what we get from full conditional of sigma2
    # prior on sigma2 ~ IG(2.1, 2.5)
    post_alpha = sum(n1, n2, n3)/2 
    post_beta = 0.5 * sst0
    sigma2 <- 1/rgamma(1, shape = post_alpha, rate = post_beta)
    
    
    # Update mu
    # what we get from full conditional of mu
    # prior on mu ~ N(mu0 = 20, tau0^2 = 10^4)
    
    mu1 <- rnorm(1, mean = y1bar, sd = sqrt(sigma2/n1))
    mu2 <- rnorm(1, mean = y2bar, sd = sqrt(sigma2/n2))
    mu3 <- rnorm(1, mean = y3bar, sd = sqrt(sigma2/n3))
    
    # updating sst0 for next draw of sigma2 
    sst0 = sum((y[[1]] - mu1)^2, (y[[2]] - mu2)^2, (y[[3]] - mu3)^2)
    
    tempDf <-  rbind(tempDf, c(mu1, mu2, mu3, sigma2))
  }
  
  output <- tempDf %>% 
    .[(burn_in+1):num_iter, ] %>% 
    set_names(c("mu1", "mu2", "mu3", "sigma2"))
  
  return(output)
}

gibbs_p4 <- function(y, num_iter, burn_in, v0, initial_sigma2, a0, b0) {    
  
  # testvalues
  # y <- grouped_dataset; v0 = c(0,0,0); initial_sigma2 = 1
  set.seed(729)

  # parameters
  n1 = length(y[[1]]); n2 = length(y[[2]]); n3 = length(y[[3]])
  mu0_1 = v0[1]; mu0_2 = v0[2]; mu0_3 = v0[3] # initial values of mu_i to update sigma2
  y1bar = mean(y[[1]]); y2bar = mean(y[[2]]); y3bar = mean(y[[3]])
  sst0 = sum((y[[1]] - mu0_1)^2, (y[[2]] - mu0_2)^2, (y[[3]] - mu0_3)^2)
  
  sigma2 <- initial_sigma2
  accepted_proposal <- 0
  tempDf <- data.frame()
  # tempDf <- data.frame(m = numeric(), v = numeric())
  
  for (j in 1:num_iter) { 
    
    # update sigma2 using Metropolis-Hastings
    
      # propose a sigma^2 value from Gamma density proposal
      proposed_sigma2 <- 1/rgamma(1, shape = sum(n1, n2, n3)/2, rate = sst0/2)
    
      # Target ratio: decides whether we accept the proposed sigma^2 as current
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
    # prior on mu ~ N(mu0 = 20, tau0^2 = 10^4)
    a11 = 1/10^4 + n1/sigma2; b11 = 20/10^4 + (n1*y1bar)/sigma2
    mu1 <- rnorm(1, mean = b11/a11, sd = sqrt(a11)^-1)
    
    a12 = 1/10^4 + n2/sigma2; b12 = 20/10^4 + (n2*y2bar)/sigma2
    mu2 <- rnorm(1, mean = b12/a12, sd = sqrt(a12)^-1)
    
    a13 = 1/10^4 + n3/sigma2; b13 = 20/10^4 + (n3*y3bar)/sigma2
    mu3 <- rnorm(1, mean = b13/a13, sd = sqrt(a13)^-1)
    
    # updating sst0 for next draw of sigma2 
    sst0 = sum((y[[1]] - mu1)^2, (y[[2]] - mu2)^2, (y[[3]] - mu3)^2)
    
    tempDf <-  rbind(tempDf, c(mu1, mu2, mu3, sigma2))
  }
  
  output <- tempDf %>% 
    .[(burn_in+1):num_iter, ] %>% 
    set_names(c("mu1", "mu2", "mu3", "sigma2"))
  
  acceptance_rate <- accepted_proposal/num_iter
  cat(" Acceptance rate using Gamma proposal: ", round(acceptance_rate, 3))
  return(output)
}

# for prior 1 and prior 2
posterior_sample_1 <- gibbs_p1(y = grouped_dataset, num_iter = 10000, burn_in = 5000, v0 = c(0, 0, 0), alpha0 = 2.1, beta0 = 2.5)
posterior_sample_2 <- gibbs_p2(y = grouped_dataset, num_iter = 10000, burn_in = 5000, v0 = c(0, 0, 0), alpha0 = 2.1, beta0 = 2.5)

# for prior 3
posterior_sample_3 <- gibbs_p3(y = grouped_dataset, num_iter = 10000, burn_in = 5000, v0 = c(0, 0, 0))

# for prior 4
posterior_sample_4 <- gibbs_p4(y = grouped_dataset, num_iter = 10000, burn_in = 5000, v0 = c(2, 5, 2), initial_sigma2 = 2)

post_summary <- posterior_sample_4

apply(post_summary, 2, median)
apply(post_summary, 2, quantile, probs = c(.025, 0.975))



graphics.off()
par(mfrow = c(2, 2))
plot(1:5000, post_summary$mu1[1:5000], "l", main = expression(paste("Posterior of ", mu[1])))
abline(h = 2.54, col="red3", lty = 2, lwd = 3)
plot(1:5000, post_summary$mu2[1:5000], "l", main = expression(paste("Posterior of ", mu[2])))
abline(h = 6.2, col="red3", lty = 2, lwd = 3)
plot(1:5000, post_summary$mu3[1:5000], "l", main = expression(paste("Posterior of ", mu[3])))
abline(h = 1.58, col="red3", lty = 2, lwd = 3)
plot(1:5000, post_summary$sigma2[1:5000], "l", main = expression(paste("Posterior of ", sigma^2)))
abline(h = 1.42, col="red3", lty = 2, lwd = 3)


#---------------------
# Likelihood plot
#---------------------

normal_likelihood <- function(data, mu, sigma2){
  likelihood <- NULL
  n <- length(data)
  if(length(mu) == length(sigma2)){
    
    for(i in 1: length(mu)){
      
      L <- (2*pi*sigma2[i])^-(n/2) * exp(-(1/2*sigma2[i])*sum((data - mu[i])^2))
      likelihood <- rbind(likelihood, L)
    }
  } else {
    cat("Check out the length of the parameters!! \n")
  }
   return(likelihood)
}


test <- normal_likelihood(data = grouped_dataset[[1]], mu = sample_g1_p1$mu, sigma2 = sample_g1_p1$sigma2)

ldf <- sample_g1_p1 %>% 
  mutate(likelihood = normal_likelihood(data = grouped_dataset[[1]], mu = sample_g1_p1$mu, sigma2 = sample_g1_p1$sigma2))

ggplot(data = ldf) +
  # geom_tile(aes(x = mu, y = sigma2)) +
  stat_contour(aes(x = mu, y = sigma2, z = likelihood))
  # geom_point(data = opt_value_rev, aes(x = K, y = b), color = "red")