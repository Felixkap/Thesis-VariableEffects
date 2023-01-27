rm(list=ls())
library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ranger)
library(MixMatrix)
library(mvtnorm)

source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/RangerPredictFunction.R')


### Model Design
n <- 200
x <- matrix(rnorm(n, 0, 1), nrow = 200)
e <- rnorm(n, 0, 1)
y <- 2*-x + e
data <- data.frame(x, y)

### Fit Linear Regression Model 
linreg <- lm(formula = y ~ x, 
             data = data)
linreg.predictions <- linreg$fitted.values
linreg.se <- predict(linreg, newdata = data, se.fit = T)$se.fit

### Fit Random Forest Model from 'ranger'

rf <- ranger( formula = y ~ x, replace = F,
              data = data, 
              num.trees = 10, 
              keep.inbag = T, 
              min.node.size = 5) # min.node.size to smooth forest fit



### Predict RF Responses (Training Data)
rf.predict <- RangerForestPredict(rf$forest,
                                  data, 
                                  predict.all = T,
                                  inbag.counts = rf$inbag.counts)

# For each observation we get all predictions from all trees
rf.all.predictions <- rf.predict$predictions

# For each observation average predctions over all trees
rf.predictions <- rowMeans(rf.all.predictions)

#### Standard Error of Random Forest Predictions
rf.se <- RangerForestPredict(rf$forest, 
                             data, 
                             type = 'se', 
                             se.method = 'jack', 
                             inbag.counts = rf$inbag.counts)$se

rf.se2 <- predict(rf, type = 'se', se.method = 'infjack', data = data[,1, drop=F])



data <- cbind(data, linreg.predictions, rf.predictions, linreg.se, rf.se)


# Plot Linear Regression and Random Forest fit
# Random Forest Predictions are plotted including standard errors
gg_fit <- ggplot(data = data) + 
  geom_point(aes(x, y, colour='Data'), alpha=0.3) + 
  geom_line(aes(x=x, y=linreg.predictions, colour='LinReg')) +
  geom_point(aes(x=x, y=rf.predictions, colour='Predictions.RF'), alpha=0.2) +
  geom_line(aes(x=x, y=rf.predictions, colour='RF')) +
  scale_color_manual(name='Methods', 
                     values = c('Data'='black',
                                'LinReg'='red', 
                                'Predictions.RF'='green', 
                                'RF'='green'),
                     guide=guide_legend(override.aes = list(
                       linetype = c('blank', 'solid', 'blank', 'solid'),
                       shape = c(16, NA, 16, NA)
                     ))) +
  theme_bw()


gg_se <- ggplot(data = data, aes(x=y,y=rf.predictions)) + 
  geom_point(alpha=0.3) +
  geom_errorbar(aes(ymin=rf.predictions - 0.5*rf.se,
                    ymax=rf.predictions + 0.5*rf.se), width=.1)+
  geom_abline(slope = 1, linetype = 2 , color = 'red')+
  labs(x = 'True Y', y = 'Predicted Y', 
       title = 'Random Forest Prediction with Standard Errors 
       vs. True Outcome Variable')+
  theme_bw()


gg_fit
gg_se
#ggarrange(gg_fit, gg_se, nrow = 2)









###### Variable Effect and Standard Error #####
x_bar <- mean(x)
x_a <- mean(x) - sd(x)
x_b <- mean(x) + sd(x)
new_data <- data.frame(x = c(x_a, x_b))


#### CoVariance Matrix of two Random Forest Predictions
rf.predict <- RangerForestPredict(rf$forest, 
                                  new_data, 
                                  type='se', 
                                  se.method = 'jack_cov',
                                  predict.all = T,
                                  inbag.counts = rf$inbag.counts)

# Variance - Covariance Matrix between two predictions
rf.predict$cov






##### Simulation --- Distribution of Variable Effect and Standard Error
# Fit Random Forest n times
# Compute Variable Effect (as we defined it) n times
# Estimate Standard Error of Variable Effect n times directly
# ---> Based on estimated Variance-Covariance Matrix

# Compare Standard Error of simulated Variable Effects with mean of estimated Standard Errors of Variable Effects
sim_once <- function(N, num.trees, formula){
  
  ### Simulate Data
  n_vars <- length(unique(unlist(str_extract_all(formula,"x.[0-9]")))) # Number of variables
  if (n_vars == 1) {
    # Normal Distribution with Mean 0 and Variance 1
    x <- data.frame(x.1 = rnorm(N, 0, 1))
  } else{
    # Multivariate Normal Distribution with Mean Vector 0 and diagonal Variance Covariance Matrix 
    x <- data.frame(x = rmvnorm(N, mean = rep(0, n_vars), sigma = CSgenerate(n_vars, cor)))
  }
  
  
  # Get interaction terms
  all_terms <- unlist(strsplit(formulas, "(?<=.)(?=[+])|(?<=.)(?=[-])",perl = TRUE))
  n_coefs <- length(all_terms)
  for (term in all_terms) {
    if (length(unique(unlist(str_extract_all(term,"x.[0-9]")))) > 1) {
      interaction_term <- unique(unlist(str_extract_all(term,"x.[0-9]")))
    }
  }
  
  
  formula <- parse(text = formula)
  y <- eval(formula, x) + rnorm(N, 0, 1)
  data <- data.frame(cbind(x, y))
  
  
  
  rf <- ranger( formula = y ~ ., 
                data = data, 
                num.trees = num.trees, 
                keep.inbag = T)
  
  print('OK')
  ### Variable Effect
  x_bar <- mean(x)
  x_a <- mean(x) - sd(x)
  x_b <- mean(x) + sd(x)
  new_data <- data.frame(x = c(x_a, x_b))
  
  rf.predict <- RangerForestPredict(rf$forest, 
                                    new_data, 
                                    type='se', 
                                    se.method = 'jack_cov',
                                    predict.all = T,
                                    inbag.counts = rf$inbag.counts)
  
  
  
  ab_predictions <- rowMeans(rf.predict$predictions)
  effect <- (ab_predictions[2] - ab_predictions[1]) / 2*sd(x) 
  ### VAR[(f(B) - f(A)) / 2*sd(x)] 
  #   =( 1 / 4*sd(x)^2 ) * ( VAR[f(B)] + VAR[f(B)] - 2*COV[f(A), f(B)] )
  smaller_null <- sum((rf.predict$cov[1,1] + rf.predict$cov[2,2] 
                       - 2*rf.predict$cov[1,2]) /  (2*sd(x))^2 < 0)
  
  effect.var <- pmax((rf.predict$cov[1,1] + rf.predict$cov[2,2] 
                      - 2*rf.predict$cov[1,2]) /  (2*sd(x))^2 , 0)
  
  effect.se <- sqrt(effect.var)
  
  
  return(list(effect = effect, 
              effect.se = effect.se, 
              smaller_null = smaller_null))
  
}






sim_multi <- function(scenario){
  
  
  N <- scenario[["N"]]
  num.trees <- scenario[["N_Trees"]]
  formula <- scenario[["Formula"]]
  repeats <- scenario[["Repeats"]]
  
  res <- replicate(n = repeats,
                   expr = sim_once(N = N, num.trees = num.trees, formula = formula))
  
  
  ### Distribution of Test Statistics
  estimates <- data.frame(effect = unlist(res[1,]), 
                          effect.se = unlist(res[2,]),
                          smaller_null = unlist(res[3,]))
  
  
  effect_plot <- ggplot(estimates, aes(x = estimates[,1])) + 
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.1, color = 'black') +
    geom_density(alpha = .2, fill = "#FF6666") + 
    labs(x = 'Effect-Size Estimates', 
         title = 'Distribution of Variable Effect Estimates')+
    theme_bw()
  
  effect_se_plot <- ggplot(estimates, aes(x = estimates[,2])) + 
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.05, color = 'black') +
    geom_density(alpha = .2, fill = "#FF6666") + 
    labs(x = 'Standard Error Estimates', 
         title = 'Distribution of Standard Error Estimates of Variable Effect')+
    theme_bw()
  
  effect_plot_box <- ggplot(estimates, aes(x = estimates[,1])) + 
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.1, color = 'black') +
    geom_density(alpha = .2, fill = "#FF6666") + 
    labs(x = 'Effect-Size Estimates', 
         title = 'Distribution of Variable Effect Estimates')+
    theme_bw()
  
  effect_se_plot_box <- ggplot(estimates, aes(x = estimates[,2])) + 
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.05, color = 'black') +
    geom_density(alpha = .2, fill = "#FF6666") + 
    labs(x = 'Standard Error Estimates', 
         title = 'Distribution of Standard Error Estimates of Variable Effect')+
    theme_bw()
  
  plot_list <- list(effect_plot, effect_se_plot,
                    effect_plot_box, effect_se_plot_box)
  
  
  
  cat('Setting (N, num.trees, formula, repeats): ', 
      c(N, num.trees, formula, repeats),
      '\nMean of simulated Variable Effects: ', 
      mean(estimates$effect),
      '\nStandard Error of simulated Variable Effects: ', 
      sd(estimates$effect), 
      '.\nMean of stimulated estimates of Standard Errors of Variable Effects: ', 
      mean(estimates$effect.se), 
      '.\nSmaller Nulls: ', 
      sum(estimates$smaller_null), '\n')
  
  
  res_scenario <- list(formula, N, num.trees, mean(estimates$effect), 
                       sd(estimates$effect), mean(estimates$effect.se), sum(estimates$smaller_null),
                       plot_list)
  
  return(res_scenario)
}



###### Simulation Setup
formulas <- c("2*x.1")
n <- c(400)
num.trees <- 2000
repeats <- 1e3
scenarios <- data.frame(expand.grid(n, num.trees, formulas, repeats))
colnames(scenarios) = c("N", "N_Trees", "Formula", "Repeats")
scenarios[,"Formula"] <- as.character(scenarios[,"Formula"]) ### Formula became Factor
scenarios <- split(scenarios, seq(nrow(scenarios)))


result <- lapply(X = scenarios, FUN = sim_multi)

#test <- do.call(rbind, result)
#colnames(test) <- c("Formula", "N", "N_trees",  
#                    "Mean Effect", "SD Effect", "SE Estimate Effect", 
#                    "SE2 Estimate Effect", "N_Nulls")





#### How to get a smooth random forest fit
#### Lots of trees --> No negative standard error estimates of variable effects
#### Check Bias Correted Version of Covariance Estimation
#### Check Direkt SE Estimate without Bias corrected version
#### Compute non-parametric Bootstrap estimate for standard error of varibale effect 