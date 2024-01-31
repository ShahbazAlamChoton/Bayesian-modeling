rm(list = ls())
graphics.off()
current_path <- rstudioapi::getSourceEditorContext()$path
setwd(dirname(current_path))
# library(MCMCpack)
library(tidyverse)
# library(dplyr)
# library(DAAG)
# library(igraph)
# library(CARBayes)
set.seed(729)
options(scipen = 999)
source("./common.function.R")

#=========== data Setup
rawData <- get(load("DATA PATH"))
rm(cross_sectional_data_01)

# Covariate and covariate matrix
covariate.data <- rawData$covariate_data # raw covariates info
covariate.scaled <- apply(covariate.data, 2, unitization)

# X <- matrix(cbind(1, covariate.scaled), nrow = nrow(covariate.scaled), byrow = F)
# n <- nrow(X); p <- ncol(X) - 1

# Response variable
count.df <- rawData$count_data[, -c(6, 9, 10)] # exclude column 6, 9, 10

# replace na with median:
count.df[which(is.na(count.df)) %% nrow(count.df), ceiling(which(is.na(count.df))/nrow(count.df))] <- apply(as.matrix(count.df[, ceiling(which(is.na(count.df))/nrow(count.df))]), 2, median, na.rm = T)
bird.names <- colnames(count.df) # contains only bird columns
#===============

#=============== Variables for poisson model (also for ZIP poisson part)
pois.covariates <- list()

for (i in 1:ncol(count.df)) {
  
  count.col = i
  bird.name <- bird.names[count.col]
  response.var <- count.df[, count.col]
  
  pois.glm.df <- dplyr:: bind_cols(response = response.var, covariate.scaled)
  
  pois.null.model = glm(response ~ 1, data = pois.glm.df, family = poisson(link = "log"))
  pois.full.model = glm(response ~ ., data = pois.glm.df, family = poisson(link = "log"))
  
  full.forward = step(pois.null.model,
                      scope = list(lower = formula(pois.null.model),
                                   upper = formula(pois.full.model)),
                      direction = "forward")

  pois.covariates[[bird.name]] <- names(full.forward$coef)[-1] # the first one is intercept
  rm(pois.glm.df, full.forward, pois.null.model, pois.full.model)
}

save(pois.covariates, file = paste0("./output/poisson_zip_poisson_part_vs.RData"))

#=============== Variables for zip_simpler poisson part
zip_simpler.pois.covariates <- list()

for (i in 1:ncol(count.df)) {
  
  count.col = i
  bird.name <- bird.names[count.col]
  response.var <- count.df[, count.col] - 1 #'[exception due to ZIP_simpler model]
  
  pois.glm.df <- dplyr:: bind_cols(response = response.var, covariate.scaled) %>% 
    filter(response >= 0) #'[ Because all those zero counts comes from Bernoulli part]
  
  pois.null.model = glm(response ~ 1, data = pois.glm.df, family = poisson(link = "log"))
  pois.full.model = glm(response ~ ., data = pois.glm.df, family = poisson(link = "log"))
  
  full.forward = step(pois.null.model,
                      scope = list(lower = formula(pois.null.model),
                                   upper = formula(pois.full.model)),
                      direction = "forward")
  
  zip_simpler.pois.covariates[[ bird.name ]] <- names(full.forward$coef)[-1] # the first one is intercept
  rm(pois.glm.df, full.forward, pois.null.model, pois.full.model)
}

save(zip_simpler.pois.covariates, file = paste0("./output/zip_simpler_poisson_part_vs.RData"))

#=============== Variable choice for Bernoulli part (for both ZIP & ZIP_simpler)
bin.covariates <- list()

for (i in 1:ncol(count.df)) {
  
  count.col = i
  bird.name <- bird.names[count.col]
  response.var <- count.df[, count.col]
  
  bin.glm.df <- dplyr:: bind_cols(response = if_else(response.var > 0, 1, 0), covariate.scaled)
  
  bin.null.model = glm(response ~ 1, data = bin.glm.df, family = binomial(link = "probit"))
  bin.full.model = glm(response ~ ., data = bin.glm.df, family = binomial(link = "probit"))
  
  full.forward = step(bin.null.model,
                      scope = list(lower = formula(bin.null.model),
                                   upper = formula(bin.full.model)),
                      direction = "forward")
  
  bin.covariates[[ bird.name ]] <- names(full.forward$coef)[-1] # the first one is intercept
  rm(bin.glm.df, full.forward, bin.null.model, bin.full.model)
  
  #'[ Warning: glm.fit: fitted probabilities numerically 0 or 1 occurred]
  #'[ Warning answer: https://stats.stackexchange.com/questions/336424/issue-with-complete-separation-in-logistic-regression-in-r]
  #'[Another sourse: https://www.statology.org/glm-fit-fitted-probabilities-numerically-0-or-1-occurred/]
}

save(bin.covariates, file = paste0("./output/zip_zip_simpler_bernoulli_part_vs.RData"))
