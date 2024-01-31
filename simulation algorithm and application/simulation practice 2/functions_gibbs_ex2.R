# Gibbs sampling for univariate distribution

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


# Gibbs sampling for univariate distribution

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
