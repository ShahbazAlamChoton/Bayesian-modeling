rm(list = ls())
graphics.off()
# setwd
library(MCMCpack)
library(tidyverse)

#-----------------------
# Custom function
#-----------------------
prob_matrix <- function(grid_length, m) {
  
  prob_mat <- matrix(0, nrow = grid_length, ncol = grid_length)
  
  for(i in 1 : m) {
    prob_mat[i, setdiff(1:(m+i), i)] = round((m + (i-1))^-1, 3)
  }
  
  for(i in (m+1) : (grid_length - m)) {
    prob_mat[i, setdiff((i-m):(i+m), i)] = round((2 * m)^-1, 3)
  }
  
  
  for( i in (grid_length - m + 1) : grid_length) {
    prob_mat[i, setdiff((i-m):grid_length, i)] = round(length(setdiff((i-m):grid_length, i))^-1, 3)
  }
  
  return(prob_mat)
}

# #-----------------------
# # data setup
# #-----------------------
traindf <- read_delim("./data/Point_referenced_GP_simulated_training_data_Jul_20_2023.txt", delim = " ") %>% 
  distinct(`horizontal_coordinate`, `vertical_coordinate`, .keep_all = T)
coorddf <- subset(traindf, select = c("horizontal_coordinate", "vertical_coordinate")) 
# dupdf <- traindf %>% 
#   filter(duplicated(cbind(`horizontal_coordinate`, `vertical_coordinate`)))
# 
# dupdf <- traindf %>% 
#   group_by(`horizontal_coordinate`, `vertical_coordinate`) %>% 
#   summarise(n = n(), .groups = 'drop') %>% 
#   filter(n > 1)

# dcord <- traindf %>% 
#   filter(`horizontal_coordinate` %in% dupdf$horizontal_coordinate & 
#            `vertical_coordinate` %in% dupdf$vertical_coordinate)
# horizontal: 97 unique value
#vertical: 167 unique value
# 4 duplicated value found in coordinates

n <- nrow(traindf)

#-----------------------
# Distance, Correlation R, lambda grid, lambda probability Calculation
#-----------------------
# distance matrix
dismat <- as.matrix(dist(coorddf, method = "euclidean", diag = T, upper = T, p = 2))
dmax <- max(dismat) # unique(sort(dismat, decreasing = T))[1]
dmin <- unique(sort(dismat, decreasing = F))[2]

# lambda grid
lambda_gridlength <- 150
point_consider <- 2
lambda_grid <- seq(from = 3/dmax, to = 3/dmin, length.out = lambda_gridlength) #lambda_range = (1.184, 300)
# probability matrix for lambda
lambda_prob_matrix <- prob_matrix(grid_length = lambda_gridlength, m = point_consider)

# # Validation
# first5 <- lambda_prob_matrix[1:5,1:20]
# last5 <- lambda_prob_matrix[96:100,91:100]

# R initial
R <- exp(-dismat)

#---------------------
# initials
#---------------------
set.seed(729)
num_iter = 50000; burn_in = 10000; thin_in = 5
# parameters: Beta, theta, sigma2, tau2, lambda
betas = c(2, 0, 0); theta_init = matrix(0, nrow = n, ncol = 1); 
sigma2_init = 0.4; tau2_init = 0.4; lambda_init = lambda_grid[14]

# prior on sigma2 and tau02 
a_tau_pr = 2.5 ; b_tau_pr = 2.5
a_sigma_pr = 2.5; b_sigma_pr = 2.5


sppr_draws <- function(traindf, num_iter, burn_in, thin_in, 
                                        betas, theta_init, sigma2_init, tau2_init, lambda_init,
                                        a_tau_pr, b_tau_pr, a_sigma_pr, b_sigma_pr) {    
  
  set.seed(729)
  
  #-----------------------
  # data setup
  #-----------------------
  # traindf <- read_delim("./data/Point_referenced_GP_simulated_training_data_Jul_20_2023.txt", delim = " ")
  # coorddf <- subset(traindf, select = c("horizontal_coordinate", "vertical_coordinate"))
  Y = matrix(traindf$response, ncol = 1)
  n = nrow(traindf)
  X = matrix(c(rep(1, n), traindf$covariate_01, traindf$covariate_02), ncol = 3, byrow = F)
  tX = t(X)
  p = ncol(X) - 1
  
  #-----------------------
  # prior info
  #-----------------------
  # prior information for all except beta are coming from function parameters
  # Beta prior setup
  mu_beta = matrix(c(2, 2, 1), ncol = 1);  sigma_beta = 10^4 * diag(p + 1)

  #---------------------------
  # Assigning initial values 
  #---------------------------
  # for beta: need to create a matrix from 'betas' (beta_initial)
  beta = matrix(betas, ncol = 1)
  # Assigning initial values to parameters
  theta = theta_init
  tau2 = tau2_init
  # sigma2 = sigma2_init
  # To draw lambdas using Metropolis-Hastings: 
  # _c for 'current' value
  #  _p for proposed value
  lambda_c = lambda_init
  R_c = R^lambda_c; det_R_c = det(R_c); inv_R_c = chol2inv(chol(R_c))
  
  ss = sum((Y - X %*% beta - theta)^2)
  # ss = t(Y - X %*% beta - theta) %*% (Y - X %*% beta - theta)
  
  accepted_proposal <- 0
  tempDf <- data.frame()
  pb <- txtProgressBar(min = 1, max = num_iter, style = 3)
  
  for (j in 1:num_iter) {
    
    # tic = Sys.time()

    # propose a lambda value from predefined probability matrix
    lambda_p <- sample(lambda_grid, 1, prob = lambda_prob_matrix[which(lambda_grid == lambda_c), ])
    R_p <- R^lambda_p; det_R_p <- det(R_p); inv_R_p <- chol2inv(chol(R_p))
    
    proposal <- lambda_prob_matrix[which(lambda_grid == lambda_p), which(lambda_grid == lambda_c)]/ lambda_prob_matrix[which(lambda_grid == lambda_c), which(lambda_grid == lambda_p)]
    target_by_proposal <- -.5*( (log(det_R_p) - log(det_R_c)) +  (tau2)^-1 * (t(theta) %*% (inv_R_p - inv_R_c) %*% theta)) + log(proposal)
    
    acceptance_ratio <- min(1, target_by_proposal) # markov kernel cancel out
    # print(paste0("Acc ratio: ", acceptance_ratio," sigma2 ",  sigma2, " Proposal: ", proposed_sigma2))
    # pause()
    
    if(log(runif(1)) < acceptance_ratio) {
      
      # if accepts, change current lambda to the proposed value
      lambda <- lambda_p; R_lambda <- R_p; R_lambda_inv = inv_R_p
      R_c = R_p; det_R_c = det_R_p; inv_R_c = inv_R_p
      
      accepted_proposal <- accepted_proposal + 1
      
    } else {
      
      lambda <- lambda_c; R_lambda = R_c; R_lambda_inv = inv_R_c
    }
    
    # R_lambda <- R^lambda
    # Update tau2
    # what we get from full conditional of tau2
    # prior on sigma2 ~ IG( 2.5, 2.5)
    A_tau2 = (n/2) + a_tau_pr
    B_tau2 = (t(theta) %*% R_lambda_inv %*% theta)/2 + b_tau_pr
    tau2 <- 1/rgamma(1, shape = A_tau2, rate = B_tau2)
    
    # Update sigma2
    # what we get from full conditional of sigma2
    # prior on sigma2 ~ IG( 2.5, 2.5)
    A_sigma2 = (n/2) + a_sigma_pr; B_sigma2 = ss/2 + b_sigma_pr
    sigma2 <- 1/rgamma(1, shape = A_sigma2, rate = B_sigma2)
    
    # Update theta
    # what we get from full conditional of theta
    sigma_inv = (sigma2^-1) * diag(n); sigma_theta_inv = (tau2^-1) * R_lambda_inv
    A_theta = sigma_inv + sigma_theta_inv
    B_theta = sigma_inv %*% Y - sigma_inv %*% X %*% beta
    A_theta_inv = chol2inv(chol(A_theta))
    theta <- mvrnorm(1, A_theta_inv %*% B_theta, A_theta_inv)
    
    # Update beta
    # what we get from full conditional of beta
    sigma_beta_inv = chol2inv(chol(sigma_beta))
    x_transpose_sigma_inv <- tX %*% sigma_inv
    A_beta = x_transpose_sigma_inv %*% X + sigma_beta_inv; A_beta_inv = chol2inv(chol(A_beta))
    B_beta = x_transpose_sigma_inv %*% Y - x_transpose_sigma_inv %*% theta + sigma_beta_inv %*% mu_beta
    beta <- mvrnorm(1, A_beta_inv %*% B_beta, A_beta_inv)

    
    # updating ss for next draw of sigma2 
    ss = sum((Y - X %*% beta - theta)^2)
    # ss = t(Y - X %*% beta - theta) %*% (Y - X %*% beta - theta) # re: sum of squares
    lg_lkhd <- -0.5 * (n * log(sigma2) + sigma2^(-1) * ss)
    
    tempDf <-  rbind(tempDf, c(lg_lkhd, beta, sigma2, tau2, lambda, theta))
    
    setTxtProgressBar(pb, j) # to see the progress
    # toc = Sys.time()
    # print(paste0("Iteration # ", j, " | Time taken: ", round(as.numeric(toc - tic, units = "secs"), 2)))
  }
  
  output <- tempDf %>% 
    .[(burn_in+1):num_iter, ] %>%
    .[seq(from = 1, to = num_iter - burn_in, by = thin_in), ] %>% 
    set_names(c("lg_lkhd", paste0("beta", 0:p),  "sigma2", "tau2", "lambda", paste0("theta_", 1:n)))
  
  acceptance_rate <- (accepted_proposal/num_iter)*100
  cat(" Acceptance rate of lambda: ", acceptance_rate, "%")
  return(output)
}

#-------------------- Start: Spatial Regression draws
tic = Sys.time()
spprDraws <- sppr_draws(traindf = traindf, num_iter = num_iter, burn_in = burn_in, thin_in = thin_in,
                                       betas = betas, theta_init = theta_init, 
                                       sigma2_init = sigma2_init, tau2_init = tau2_init, lambda_init = lambda_init,
                                       a_tau_pr = a_tau_pr, b_tau_pr = b_tau_pr,
                                       a_sigma_pr = a_sigma_pr, b_sigma_pr = b_sigma_pr)
# close(pb)
toc <- Sys.time()
paste0("Completion Time : ", round(as.numeric(toc - tic, units = "mins"), 2), " Minutes")

spprDraws <- spprDraws %>% 
  mutate(index = 1:nrow(spprDraws),
         th_mean = rowMeans(select(., starts_with("theta")), na.rm = T)) %>% 
  select(index, lg_lkhd, th_mean, everything())

#-------------------- End: Spatial Regression draws

#---------------------------
# Prediction
#---------------------------

# test data preperation
testdf <- read_delim("./data/Point_referenced_GP_simulated_test_data_Jul_20_2023.txt", delim = " ")
n_train <- nrow(traindf); n_test <- nrow(testdf)
testcord <- subset(testdf, select = c("horizontal_coordinate", "vertical_coordinate"))
X_test = matrix(c(rep(1, n_test), testdf$covariate_01, testdf$covariate_02), ncol = 3, byrow = F)
p_test <- ncol(X_test) - 1
# combined distance matrix: required for calculating conditional MVN
ttcord <- bind_rows(coorddf, testcord)
dismat_tt <- as.matrix(dist(ttcord, method = "euclidean", diag = T, upper = T, p = 2))

prdf <- data.frame()
for( i in 1:nrow(spprDraws)) {  
  
  # Parameters
  b <- matrix(as.matrix(spprDraws[i, c(paste0("beta", 0:p_test))]), ncol = 1)
  th <- matrix(as.matrix(spprDraws[i, c(paste0("theta_", 1:n_train))]), ncol = 1)
  tu2 <- spprDraws$tau2[i]; s2 <- spprDraws$sigma2[i]
  l <- spprDraws$lambda[i]
  r_tt <- exp(-dismat_tt)^l
  
  #  Matrix partition
  r_1 <- r_tt[1:n_train, 1:n_train]; r_10 <- r_tt[1:n_train, (n_train + 1):(n_train + n_test)]
  r_01 <- r_tt[(n_train + 1):(n_train + n_test), 1:n_train]; r_0 <- r_tt[(n_train + 1):(n_train + n_test), (n_train + 1):(n_train + n_test)]
  
  # Conditional MVN parameters for predicted response
  r_1_inv <- chol2inv(chol(r_1)); common_part <- r_01 %*% r_1_inv
  mu_s0_given_s1 <- common_part %*% th
  sigma_s0_given_s1 <- tu2 * (r_0 - common_part %*% r_10)
  
  mu_mvn <- matrix(0, nrow = n_test, ncol = 1)
  ptemp <- X_test %*% b + mvrnorm(1, mu = mu_s0_given_s1, Sigma = sigma_s0_given_s1) + mvrnorm(1, mu = mu_mvn, Sigma = s2*diag(n_test))
  
  prdf <- rbind(prdf, t(ptemp))
  rm(b, th, tu2, s2, l, r_tt, r_1, r_10, r_01, r_0, r_1_inv, common_part, mu_s0_given_s1, sigma_s0_given_s1, ptemp)

}

sp_predictiondf <- prdf %>% 
  set_names(c(paste0("y_", 1:n_test)))


sp_pred_medians <- apply(sp_predictiondf, 2, median)
sp_pred_quantiles <- apply(sp_predictiondf, 2, quantile, probs = c(.05, 0.95))
true_response <- testdf$response

#---------------------------------------------------------------------------------------------
#                               Multivariate Regression
#---------------------------------------------------------------------------------------------

# Custom function to draw from a non-spatial regression
mvr_draws <- function(traindf, num_iter, burn_in, thin_in, 
                      betas,
                      a_sigma_pr, b_sigma_pr) {    
  
  set.seed(729)
  #-----------------------
  # data setup
  #-----------------------
  # traindf <- read_delim("./data/Point_referenced_GP_simulated_training_data_Jul_20_2023.txt", delim = " ")
  # coorddf <- subset(traindf, select = c("horizontal_coordinate", "vertical_coordinate"))
  Y = matrix(traindf$response, ncol = 1)
  n = nrow(traindf)
  X = matrix(c(rep(1, n), traindf$covariate_01, traindf$covariate_02), ncol = 3, byrow = F)
  tX = t(X)
  p = ncol(X) - 1
  
  #-----------------------
  # prior info
  #-----------------------
  # prior information for all except beta are coming from function parameters
  # Beta prior setup
  mu_beta = matrix(c(2, 2, 1), ncol = 1);  sigma_beta = 10^4 * diag(p + 1)
  
  #---------------------------
  # Assigning initial values 
  #---------------------------
  # for beta: need to create a matrix from 'betas' (beta_initial)
  beta = matrix(betas, ncol = 1)
  # sigma2 = sigma2_init
  ss = sum((Y -  X %*% beta)^2)
  # ss = t(Y - X %*% beta) %*% (Y - X %*% beta)
  tempDf <- data.frame()
  pb <- txtProgressBar(min = 1, max = num_iter, style = 3)
  for (j in 1:num_iter) {
    
    # Update sigma2
    # what we get from full conditional of sigma2: updated inverse gamma
    # prior on sigma2 ~ IG( 2.5, 2.5)
    A_sigma2 = (n/2) + a_sigma_pr; B_sigma2 = ss/2 + b_sigma_pr
    sigma2 <- 1/rgamma(1, shape = A_sigma2, rate = B_sigma2)
    
    # Update beta
    # what we get from full conditional of beta; MVN
    sigma_inv = (sigma2^-1) * diag(n); sigma_beta_inv = chol2inv(chol(sigma_beta))
    x_transpose_sigma_inv <- tX %*% sigma_inv
    A_beta = x_transpose_sigma_inv %*% X + sigma_beta_inv; A_beta_inv = chol2inv(chol(A_beta))
    B_beta = x_transpose_sigma_inv %*% Y + sigma_beta_inv %*% mu_beta
    beta <- mvrnorm(1, A_beta_inv %*% B_beta, A_beta_inv)
    
    # updating ss for next draw of sigma2 
    # ss = t(Y - X %*% beta) %*% (Y - X %*% beta)
    ss = sum((Y -  X %*% beta)^2)
    lg_lkhd <- -0.5 * (n * log(sigma2) + sigma2^(-1) * ss)
    tempDf <- rbind(tempDf, c(lg_lkhd, beta, sigma2))
    
    setTxtProgressBar(pb, j) # to see the progress
  }
  
  output <- tempDf %>% 
    .[(burn_in+1):num_iter, ] %>%
    .[seq(from = 1, to = num_iter - burn_in, by = thin_in), ] %>% 
    set_names(c("lg_lkhd", paste0("beta", 0:p), "sigma2"))
  
  return(output)
}

#-------------------- Start: Non-spatial Regression draws
tic_mv <- Sys.time()
mvrDraws <- mvr_draws(traindf = traindf, num_iter = num_iter, burn_in = burn_in, thin_in = thin_in,
                        betas = betas,
                        a_sigma_pr = a_sigma_pr, b_sigma_pr = b_sigma_pr)

mvrDraws <- mvrDraws %>% 
  mutate(index = 1: nrow(mvrDraws)) %>% 
  select(index, everything())

# close(pb)
toc_mv <- Sys.time()
paste0("Completion Time : ", round(as.numeric(toc_mv - tic_mv, units = "secs"), 2), " seconds")
#-------------------- End: Non-spatial Regression draws


# Multivariate regression prediction

pr_mv_df <- data.frame()
for (i in 1: nrow(mvrDraws)) { 
  
  b_mvr <- matrix(as.matrix(mvrDraws[i, c(paste0("beta", 0:p_test))]), ncol = 1)
  s2_mvr <- mvrDraws$sigma2[i]
  
  mu_mvn <- matrix(0, nrow = n_test, ncol = 1)
  ptemp <- X_test %*% b_mvr + mvrnorm(1, mu = mu_mvn, Sigma = s2_mvr*diag(n_test))
  pr_mv_df <- rbind(pr_mv_df, t(ptemp))
  rm(b_mvr,  s2_mvr, ptemp)
}

mv_predictiondf <- pr_mv_df %>% 
  set_names(c(paste0("y_", 1:n_test)))

mv_pred_medians <- apply(mv_predictiondf, 2, median)
mv_pred_quantiles <- apply(mv_predictiondf, 2, quantile, probs = c(.05, 0.95))

# Log-likelihood quantiles
sp_lg_q <- quantile(spprDraws$lg_lkhd, probs = c(.05, 0.95))
mv_lg_q <- quantile(mvrDraws$lg_lkhd, probs = c(.05, 0.95))

# Length of credible intervals:
print(paste0("for spatial model: ", sp_pred_quantiles[2, ] - sp_pred_quantiles[1, ]))
print(paste0("for non-spatial model: ", mv_pred_quantiles[2, ] - mv_pred_quantiles[1, ]))

# # "#69b3a2"
# graphics.off()

# # Plot 1
# par(mfrow = c(2,5))
# plot(1:nrow(spprDraws), spprDraws$beta0, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste(beta[0])))
# plot(1:nrow(spprDraws), spprDraws$beta1, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste(beta[1])))
# plot(1:nrow(spprDraws), spprDraws$beta2, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste(beta[2])))
# plot(1:nrow(spprDraws), spprDraws$sigma2, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste(sigma^2)))
# plot(1:nrow(spprDraws), spprDraws$tau2, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste(tau^2)))
# 
# hist(spprDraws$beta0, freq = FALSE, xlab = expression(beta[0]), main = expression(paste(beta[0])))
# hist(spprDraws$beta1, freq = FALSE, xlab = expression(beta[1]), main = expression(paste(beta[1])))
# hist(spprDraws$beta2, freq = FALSE, xlab = expression(beta[2]), main = expression(paste(beta[2])))
# hist(spprDraws$sigma2, freq = FALSE, xlab = expression(sigma^2), main = expression(paste(sigma^2)))
# hist(spprDraws$tau2, freq = FALSE, xlab = expression(tau^2), main = expression(paste(tau^2)))

# plot 2

# hist(sp_predictiondf$y_1, freq = FALSE, xlab = expression(y[1]), main = expression(paste(y[1])))
# abline(v = c(true_response[1], sp_pred_medians["y_1"], sp_pred_quantiles[, "y_1"]), col = c("green4", "red3", "tomato3", "tomato3"), lty = c(1, 2, 4, 4), lwd = c(3, 3))

# par(mfrow = c(2, 3), xpd = F)

# hist(mv_predictiondf$y_1, freq = FALSE,breaks = 15 , col = "steelblue2", xlab = expression(y[1]), main = expression(paste(y[1])))
# hist(sp_predictiondf$y_1, freq = FALSE, breaks = 15,  col = "lightseagreen", xlab = expression(y[1]), main = expression(paste(y[1])),  add = T)
# abline(v = c(true_response[1], sp_pred_medians["y_1"], mv_pred_medians["y_1"], sp_pred_quantiles[, "y_1"], mv_pred_quantiles[, "y_1"]), col = c("red3", "lightseagreen", "steelblue2", "lightseagreen", "lightseagreen", "steelblue2", "steelblue2"), lty = c(1, 2, 2, 4, 4, 6, 6), lwd = c(3, 3))
# 
# hist(mv_predictiondf$y_2, freq = FALSE,breaks = 15 , col = "steelblue2", xlab = expression(y[2]), main = expression(paste(y[2])))
# hist(sp_predictiondf$y_2, freq = FALSE, breaks = 15,  col = "lightseagreen", xlab = expression(y[2]), main = expression(paste(y[2])),  add = T)
# abline(v = c(true_response[2], sp_pred_medians["y_2"], mv_pred_medians["y_2"], sp_pred_quantiles[, "y_2"], mv_pred_quantiles[, "y_2"]), col = c("red3", "lightseagreen", "steelblue2", "lightseagreen", "lightseagreen", "steelblue2", "steelblue2"), lty = c(1, 2, 2, 4, 4, 6, 6), lwd = c(3, 3))
# 
# hist(mv_predictiondf$y_3, freq = FALSE,breaks = 15 , col = "steelblue2", xlab = expression(y[3]), main = expression(paste(y[3])))
# hist(sp_predictiondf$y_3, freq = FALSE, breaks = 15,  col = "lightseagreen", xlab = expression(y[3]), main = expression(paste(y[3])),  add = T)
# abline(v = c(true_response[3], sp_pred_medians["y_3"], mv_pred_medians["y_3"], sp_pred_quantiles[, "y_3"], mv_pred_quantiles[, "y_3"]), col = c("red3", "lightseagreen", "steelblue2", "lightseagreen", "lightseagreen", "steelblue2", "steelblue2"), lty = c(1, 2, 2, 4, 4, 6, 6), lwd = c(3, 3))
# 
# hist(mv_predictiondf$y_4, freq = FALSE,breaks = 15 , col = "steelblue2", xlab = expression(y[4]), main = expression(paste(y[4])))
# hist(sp_predictiondf$y_4, freq = FALSE, breaks = 15,  col = "lightseagreen", xlab = expression(y[4]), main = expression(paste(y[4])),  add = T)
# abline(v = c(true_response[4], sp_pred_medians["y_4"], mv_pred_medians["y_4"], sp_pred_quantiles[, "y_4"], mv_pred_quantiles[, "y_4"]), col = c("red3", "lightseagreen", "steelblue2", "lightseagreen", "lightseagreen", "steelblue2", "steelblue2"), lty = c(1, 2, 2, 4, 4, 6, 6), lwd = c(3, 3))
# 
# hist(mv_predictiondf$y_5, freq = FALSE,breaks = 15 , col = "steelblue2", xlab = expression(y[5]), main = expression(paste(y[5])))
# hist(sp_predictiondf$y_5, freq = FALSE, breaks = 15,  col = "lightseagreen", xlab = expression(y[5]), main = expression(paste(y[5])),  add = T)
# abline(v = c(true_response[5], sp_pred_medians["y_5"], mv_pred_medians["y_5"], sp_pred_quantiles[, "y_5"], mv_pred_quantiles[, "y_5"]), col = c("red3", "lightseagreen", "steelblue2", "lightseagreen", "lightseagreen", "steelblue2", "steelblue2"), lty = c(1, 2, 2, 4, 4, 6, 6), lwd = c(3, 3))
# # x = "bottom"
# # legend("topright", legend=c("MVR","SPR"), col=c('steelblue2', "lightseagreen"), pt.cex=2, pch=15 )
# legend(-4,-.5, legend = c("Non-spatial Regression","Spatial Regression", "True Value"),
#        col= c("steelblue2", "lightseagreen", "red3"), lwd=1, cex=.5, horiz = TRUE, 
#        xpd = NA, bty = 6)
# 

# # Trace plot
# graphics.off()
# par(mfrow = c(2, 2))
# plot(1:nrow(mvrDraws), mvrDraws$beta0, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste("Full conditional of ", beta[0])))
# # abline(h = Mu1, col="red3")
# plot(1:nrow(mvrDraws), mvrDraws$beta1, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste("Full conditional of ", beta[1])))
# # abline(h = Mu2, col="red3")
# plot(1:nrow(mvrDraws), mvrDraws$beta2, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste("Full conditional of ", beta[2])))
# # abline(h = Mu3, col="red3")
# plot(1:nrow(mvrDraws), mvrDraws$sigma2, "l", xlab = "Iteration No", ylab = "Paramater Value", main = expression(paste("Full conditional of ", sigma^2)))
# # abline(h = Sigma2, col="red3")
# 
# # Likelihood plot
# 
# # Histograms and credible intervals
# medians <- apply(mvrDraws, 2, median)
# quantiles <- apply(mvrDraws, 2, quantile, probs = c(.05, 0.95))
# 
# par(mfrow = c(2, 2))
# hist(mvrDraws$beta0, freq = FALSE, xlab = expression(beta[0]), main = expression(paste("Histogram of ", beta[0])))
# abline(v = c(medians["beta0"], quantiles[, "beta0"]),  col = c("red3", "darkgreen", "darkgreen"), lty = c(1, 4, 4), lwd = c(3, 3))
# hist(mvrDraws$beta1, freq = FALSE, xlab = expression(beta[1]), main = expression(paste("Histogram of ", beta[1])))
# abline(v = c(medians["beta1"], quantiles[, "beta1"]),  col = c("red3", "darkgreen", "darkgreen"), lty = c(1, 4, 4), lwd = c(3, 3))
# hist(mvrDraws$beta2, freq = FALSE, xlab = expression(beta[2]), main = expression(paste("Histogram of ", beta[2])))
# abline(v = c(medians["beta2"], quantiles[, "beta2"]),  col = c("red3", "darkgreen", "darkgreen"), lty = c(1, 4, 4), lwd = c(3, 3))
# hist(mvrDraws$sigma2, freq = FALSE, xlab = expression(sigma^2), main = expression(paste("Histogram of ", sigma^2)))
# abline(v = c(medians["sigma2"], quantiles[, "sigma2"]),  col = c("red3", "darkgreen", "darkgreen"), lty = c(1, 4, 4), lwd = c(3, 3))
# 

