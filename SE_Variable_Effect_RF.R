library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ranger)
library(randomForest)
library(parallel)

source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/RangerPredictFunction.R')


##### Simulation --- Distribution of Variable Effect and Standard Error
# Fit Random Forest n times
# Compute Variable Effect (as we defined it) n times
# Estimate Standard Error of Variable Effect n times directly
# ---> Based on estimated Variance-Covariance Matrix

# Compare Standard Error of simulated Variable Effects with mean of estimated Standard Errors of Variable Effects
sim_once <- function(N, num.trees, formula){
  
  
  x <- rnorm(N, 0, 1)
  e <- rnorm(N, 0, 1)
  xe <- data.frame(cbind(x,e))
  formula <- parse(text = formula)
  y <- eval(formula, xe)
  data <- data.frame(x, y)
  
  rf <- ranger( formula = y ~ x, 
                data = data, 
                num.trees = num.trees, 
                keep.inbag = T)
  
  
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
  
  rf.predict2 <- RangerForestPredict(rf$forest, 
                                    new_data, 
                                    type='se', 
                                    se.method = 'jack_effect',
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
  
  effect.se2 <- rf.predict2$se_effect
  
  return(list(effect = effect, 
              effect.se = effect.se, 
              effect.se2 = effect.se2, 
              smaller_null = smaller_null))
  
}






sim_multi <- function(scenario){
  
  
  N <- scenario[["N"]]
  num.trees <- scenario[["N_Trees"]]
  formula <- scenario[["Formula"]]
  repeats <- scenario[["Repeats"]]
  print(N)
  
  res <- replicate(n = repeats,
                   expr = sim_once(N = N, num.trees = num.trees, formula = formula))
  
  
  ### Distribution of Test Statistics
  estimates <- data.frame(effect = unlist(res[1,]), 
                          effect.se = unlist(res[2,]),
                          effect.se2 = unlist(res[3,]),
                          smaller_null = unlist(res[4,]))
  
  print(estimates)
  
  
  res_scenario <- c(formula, N, num.trees, round(mean(estimates$effect), 4), 
                    round(sd(estimates$effect), 4), round(mean(estimates$effect.se), 4), 
                    round(mean(estimates$effect.se2), 4), sum(estimates$smaller_null))
  
  return(res_scenario)
}




#### How to get a smooth random forest fit
#### Lots of trees --> No negative standard error estimates of variable effects
#### Check Bias Correted Version of Covariance Estimation
#### Check Direkt SE Estimate without Bias corrected version
#### Compute non-parametric Bootstrap estimate for standard error of varibale effect 





################################################################################
##########################     test      #######################################
################################################################################

#  formulas <- c("2*x+e", "2*-x^2", "3*sqrt(abs(x))+3*x+e")
#  n <- c(100,200,500)
#  n_trees <- c(50, 100, 200, 500)
#  repeats <- 1e3




################################################################################
##########################     test2    #######################################
################################################################################

#  formulas <- c("2*x+e", "2*-x^2", "3*sqrt(abs(x))+3*x+e")
#  n <- c(100,200,500)
#  n_trees <- c(50, 100, 200, 500)
#  repeats <- 1e4




