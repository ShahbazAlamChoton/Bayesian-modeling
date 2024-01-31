#------------- Problem 2 (2): Inverse Transform Sampling

# Defining the PDF: f(x) = (1/4) * x^3
pdf <- function(x) {
  if (0 <= x && x < 2) {
    return((1/4) * x^3)
  } else {
    return(0)
  }
}

# Calculating the CDF: F(x) = (1/16) * x^4
cdf <- function(x) {
  if (0 <= x && x < 2) {
    return((1/16) * x^4)
  } else if (x >= 2) {
    return(1)
  } else {
    return(0)
  }
}

# Defining the inverse CDF: 2 * u^(1/4)
inverse_cdf <- function(u) {
  return(2 * u^(1/4))
}

# Generate 10,000 random samples using ITS
set.seed(729) 
n_samples <- 10000
u <- runif(n = n_samples, min = 0, max = 1)
random_samples <- inverse_cdf(u)

# Histogram of the generated samples
hist(random_samples, breaks = 50, main = "Histogram of Samples using Inverse Transform Sampling",
     xlab = "Random Samples", ylab = "Frequency", col = "white")
curve(1/4 * x^3, from = 0, to = 2, add = T, col = "dodgerblue4")
text(4,0.3, expression(paste("The line shows the actual density", f(X) == 1/4 * x^3)))


# problem 3

x <- seq(0, 2, length.out = 100)
y <- .25 * (x + 3/(2 * sqrt(2)) * sqrt(x))
plot(x, y)
