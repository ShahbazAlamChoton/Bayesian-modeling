set.seed(279)

pdf <- function(x){
  ifelse(0 < x & x < 2 * pi, exp(sin(x)), 0)
}

num_iter <- 10e4
X <- numeric(length = num_iter)
X[1] <- runif(1, 0, 2 * pi) # initial X

for(j in 2:num_iter) {
  
  candidate = rnorm(1, mean = X[j-1], sd = 1) # normal proposal
  prob_acc <- min(1, pdf(candidate)/pdf(X[j-1]))
  
  X[j] <- ifelse(prob_acc > runif(1), candidate, X[j-1])
}

#b: importance sampling
pdf <- function(x){
  ifelse(0 < x & x < 2 * pi, exp(sin(x)), 0)
}
# Self-normalized importance sampling to estimate E(X^2)
num_iter <- 10e4
samples <- numeric(num_iter)
weights <- numeric(num_iter)
  
for (i in 1:num_iter) {
  # Generate samples from the proposal distribution (uniformly between 0 and 2*pi)
  sample <- runif(1, min = 0, max = 2 * pi)
  samples[i] <- sample
  
  # Calculate importance weights
  weights[i] <- pdf(sample) / dunif(sample, min = 0, max = 2 * pi)
}

# Compute the estimated expectation E(X^2)
normalized_weights <- weights / sum(weights)
estimated_expectation <- sum((samples^2) * normalized_weights)
  
# Estimate E(X^2) using self-normalized importance sampling
print(paste("Estimated E(X^2):", estimated_expectation))

# to match with what we got using M-H
cat("E(X^2) using sample of X from Metropolis-Hastings: ",mean(X^2),
    "\nEstimated E(X^2) using self-normalized importance sampling:", estimated_expectation)

# Visualisation
# MH

trace_mh <-  ggplot() + 
  geom_line(aes(x = 1:length(X), y = X), col = "black", linewidth = .5) +
  labs(title = expression(paste("Traceplot of samples of X drawn using Metropolis-Hastings")),
       x = "Iteration Number",
       y = "Value") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size = 9, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(),
    legend.text = element_blank(), legend.position = c(.95, .9),
    legend.background = element_rect(color = "black", fill= NA, linewidth = 1, linetype= "solid"),
    panel.border = element_rect(color = "black", fill = NA, size = 1))

mh_hist <- ggplot() +
  geom_histogram(aes(x = X, y = after_stat(density)), bins = 30, fill = "white", col = "black") + 
  geom_density(aes( x = X), col = "#003153", lwd = 1) +
  # geom_line(aes( x = runif(10000, 0, 2 * pi), y = pdf(runif(10000, 0, 2 * pi))), col = "seagreen", lwd = 1) +
  labs(title = "Histogram with Density Curve of X",
       x = "",
       y = "Density", color = "Blue") +
  theme_classic() + 
  theme(
    # plot.title = element_markdown(lineheight = .75, size= 9, face="bold", hjust = 0.5),
    plot.title = element_text(size= 9, face = "bold", color = 'black', hjust = 0.5),
    axis.title = element_text(size = 9, face="bold"),
    legend.title = element_blank(), legend.position = "top",
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1)
  )

library(ggpubr)
library(ggrepel)
plot_X<- ggpubr::annotate_figure(ggarrange(trace_mh, mh_hist,  
                                          common.legend = T, nrow = 2, ncol = 1), top = text_grob(""
                                                                                                  , face = "bold", size = 1)) 
print(plot_X)
