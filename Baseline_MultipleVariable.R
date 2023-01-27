
library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ranger)
library(MixMatrix)
library(mvtnorm)
library(stringr)
library(ggh4x)


source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/RangerPredictFunction.R')
source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/Helper.R')


##### Simulation --- Distribution of Variable Effect and Standard Error
# Fit Random Forest & Linear Regression Model N times
# Compute Variable Effects (as we defined it) N times
# Estimate Standard Error of Variable Effects N times directly
# ---> Based on estimated Variance-Covariance Matrix

# Compare Standard Error of simulated Variable Effects with mean of directly estimated Standard Errors of Variable Effects



sim_once <- function(k, N, num.trees, cor, formula){
  

  ### Simulate Data
  n_vars <- length(unique(unlist(str_extract_all(formula,"x.[0-9]")))) # Number of variables
  if (n_vars == 1) {
    # Normal Distribution with Mean 0 and Variance 1
    x <- data.frame(x.1 = rnorm(N, 0, 1))
  } else{
    # Multivariate Normal Distribution with Mean Vector 0 and diagonal Variance Covariance Matrix 
    x <- data.frame(x = rmvnorm(N, mean = rep(0, n_vars), sigma = CSgenerate(n_vars, cor)))
  }
  
  
  
  formula <- parse(text = formula)
  y <- eval(formula, x) + rnorm(N, 0, 1)
  data <- data.frame(cbind(x, y))
  
  
  rf <- ranger( formula = y ~ .,
                data = data,
                num.trees = num.trees,
                keep.inbag = T)
  
  linreg <- lm( formula = y ~ .,
                data = data)
  
  
  ### Variable Effect & Standard Error
  
  # True Variable Means and Standard Deviations
  true_mean <- rep(x = 0, times = n_vars)
  names(true_mean) <- unique(unlist(str_extract_all(formula,"x.[0-9]")))
  true_s <- rep(x = 1, times = n_vars)
  names(true_s) <- unique(unlist(str_extract_all(formula,"x.[0-9]")))
  
  
  # Estimated Variable Means
  x_bar <- colMeans(x)
  
  
  list_names <- c('effects.rf', 'effects.lm', 
                  #'effects.f.sample', 
                  #'effects.rf.true', 
                  'effects.f.true', 
                  'effects.se.rf', 'effects.se.lm'
                  #, 'smaller_nulls'
  )
  
  result_list <- vector(mode = 'list', length = length(list_names))
  result_list <- lapply(X = 1:length(list_names), FUN = function(X){result_list[[X]] = numeric(length = n_vars)})
  names(result_list) <- list_names
  

  for (v in 1:n_vars) {
    
    
    # New data points based on estimated variable means and standard deviations
    
    x_a <- replace(x_bar, v, x_bar[v] - k*sd(x[,v]))
    x_b <- replace(x_bar, v, x_bar[v] + k*sd(x[,v]))
    
    
    
    # New data points based on true variable means and standard deviations
    x_a.true <- replace(true_mean, v, true_mean[v] - k*true_s[v])
    x_b.true <- replace(true_mean, v, true_mean[v] + k*true_s[v])
    
    new_data <- data.frame(rbind(x_a, x_b))
    new_data.true <- data.frame(rbind(x_a.true, x_b.true))
    
    
    rf.predict <- RangerForestPredict(rf$forest,
                                      new_data,
                                      type='se',
                                      se.method = 'jack_cov',
                                      predict.all = T,
                                      inbag.counts = rf$inbag.counts)
    
    # rf.predict.true <- RangerForestPredict(rf$forest,
    #                                        new_data.true,
    #                                        type='se',
    #                                        se.method = 'jack_cov',
    #                                        predict.all = T,
    #                                        inbag.counts = rf$inbag.counts)
    # 
    lm.predict <- predict(linreg,
                          newdata = new_data,
                          se.fit = T)
    
    
    
    ### Computing Predictions
    
    # Using f_hat and sample mean&variance of X
    predictions.rf <- rowMeans(rf.predict$predictions)

    #predictions.rf.true <- rowMeans(rf.predict.true$predictions)
    predictions.lm <- lm.predict$fit
    
    # Using f and sample mean&variance of X
    #predictions.f.sample <- eval(formula, new_data)
    
    # Using f and true mean&variance of X
    predictions.f.true <- eval(formula, new_data.true)
    
    
    ### Computing Variable Effects
    
    # Using f_hat predictions based on sample mean&variance of X 
    result_list[['effects.rf']][v] <- (predictions.rf[2] - predictions.rf[1]) / (2*k*sd(x[,v]))
    #result_list[['effects.rf.true']][v] <- (predictions.rf[2] - predictions.rf[1]) / (2*k*true_s[v])
    result_list[['effects.lm']][v] <- (predictions.lm[2] - predictions.lm[1]) / (2*k*sd(x[,v]))
    
    
    # Using f predictions based on sample mean&variance of X 
    #result_list[['effects.f.sample']][v] <- (predictions.f.sample[2] - predictions.f.sample[1]) / (2*k*sd(x[,v]))
    
    # Using f predictions based on true mean&variance of X 
    result_list[['effects.f.true']][v] <- (predictions.f.true[2] - predictions.f.true[1]) / (2*k*true_s[v])
    
    
    result_list[['smaller_nulls']][v] <- sum((rf.predict$cov[1,1] + rf.predict$cov[2,2]
                                               - 2*rf.predict$cov[1,2]) /  (2*k*sd(x[,v]))^2 < 0)
    

    
    # Standard Error of Variable Effect (Coviariance Estimation)
    # VAR[(f(B) - f(A)) / 2*sd(x)]
    # =( 1 / 4*sd(x)^2 ) * ( VAR[f(B)] + VAR[f(B)] - 2*COV[f(A), f(B)] )
    effect.var <- pmax((rf.predict$cov[1,1] + rf.predict$cov[2,2]
                        - 2*rf.predict$cov[1,2]) /  (2*k*sd(x[,v]))^2, 0)
    
    
    
    
    # if (v==3) {
    #   print(paste('Variable', v))
    #   print('DataPoints')
    #   print(x_a)
    #   print(x_b)
    #   print('PREDICTIONS')
    #   print(predictions.rf)
    #   print('Variance Covariance')
    #   print(rf.predict$cov)
    #   print('SMALLER NULL')
    #   print((rf.predict$cov[1,1] + rf.predict$cov[2,2]
    #          - 2*rf.predict$cov[1,2]) /  (2*k*sd(x[,v]))^2 < 0)
    # }

    
    
    result_list[['effects.se.rf']][v] <- sqrt(effect.var)
    result_list[['effects.se.lm']][v] <- summary(linreg)$coef[,"Std. Error"][v+1]
    
  }
  
  return(result_list)
  
}






sim_multi <- function(scenario){
  

  k <- scenario[["k"]]
  N <- scenario[["N"]]
  num.trees <- scenario[["N_Trees"]]
  formula <- scenario[["Formula"]]
  repeats <- scenario[["Repeats"]]
  cor <- scenario[["Correlation"]]
  
  res <- replicate(n = repeats,
                   expr = sim_once(k=k, N = N, num.trees = num.trees,
                                   cor = cor, formula = formula))
  
  
  ### Distribution of Test Statistics
  res_scenario <- list(formula=formula, N=N, k=k, cor=cor, num.trees=num.trees,
                    effects.rf = do.call('rbind', res[1,]),
                    effects.lm = do.call('rbind', res[2,]),
                    #effects.f.sample = do.call('rbind', res[3,]),
                    #effects.rf.true = do.call('rbind', res[3,]),
                    effects.f.true = do.call('rbind', res[3,]),
                    effects.se.rf = do.call('rbind', res[4,]),
                    effects.se.lm = do.call('rbind', res[5,]),
                    n.nulls = do.call('rbind', res[6,])
  )
  
  
  return(res_scenario)
}



# ## Simulation Setup
# n <- c(50, 500) ; num.trees <- c(2000) ; repeats <- 1000; cor <- c(0); k <- c(1)
# formulas <- c("2*x.1+4*x.2-3*x.3+4*x.4-2.2*x.5")  #"-0.5*x.1^3+3*x.2+0.5*sqrt(abs(x.3))"
# scenarios <- data.frame(expand.grid(n, num.trees, formulas, repeats, cor, k))
# colnames(scenarios) = c("N", "N_Trees", "Formula", "Repeats", "Correlation", "k")
# scenarios[,"Formula"] <- as.character(scenarios[,"Formula"]) ### Formula became Factor
# scenarios <- split(scenarios, seq(nrow(scenarios)))
# 
# system.time(result <- lapply(X = scenarios, FUN = sim_multi))
# print_results(result)
# 
# effect_plots <- plot_effects(result)
# se_plot <- plot_se(result)
# effect_plots
# se_plot

n <- c(10, 20) ; num.trees <- c(40) ; repeats <- 20; cor <- c(0, 0.8); k <- c(0.2, 1)
formulas <- c(" 2*x.1+4*x.2-0.5*x.3")
scenarios <- data.frame(expand.grid(n, num.trees, formulas, repeats, cor, k))
colnames(scenarios) = c("N", "N_Trees", "Formula", "Repeats", "Correlation", "k")
scenarios[,"Formula"] <- as.character(scenarios[,"Formula"]) ### Formula became Factor
scenarios <- split(scenarios, seq(nrow(scenarios)))

system.time(result <- lapply(X = scenarios, FUN = sim_multi))

print_results(result)
effect_plots <- plot_effects(result)
se_plot <- plot_se(result)
effect_plots
se_plot

