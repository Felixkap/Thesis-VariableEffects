library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ranger)
library(MixMatrix)
library(mvtnorm)
library(stringr)
library(ggh4x)
library(randomForest)

source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/RangerPredictFunction.R')
source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/Helper.R')



sim_once <- function(k, N, num.trees, cor, formula){
  
  ### Simulate Data
  n_vars <- length(unique(unlist(str_extract_all(formula,"x.[0-9]+")))) # Number of variables
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
    if (length(unique(unlist(str_extract_all(term,"x.[0-9]+")))) > 1) {
      interaction_term <- unique(unlist(str_extract_all(term,"x.[0-9]")))
    }
  }
  
  for (rep in 1:5) {
    
    cat('REPITITION', rep)
    formula <- parse(text = formula)
    y <- eval(formula, x) + rnorm(N, 0, 1)
    data <- data.frame(cbind(x, y))
    
    rf <- ranger( formula = y ~ .,
                  data = data,
                  num.trees = num.trees,
                  keep.inbag = T)
    
    linreg <- lm( formula = y ~ .^2,
                  data = data)
    
    ### Variable Effect & Standard Error
    
    # True Variable Means and Standard Deviations
    true_mean <- rep(x = 0, times = n_vars)
    names(true_mean) <- unique(unlist(str_extract_all(formula,"x.[0-9]+")))
    true_s <- rep(x = 1, times = n_vars)
    names(true_s) <- unique(unlist(str_extract_all(formula,"x.[0-9]+")))
    
    # Estimated Variable Means
    x_bar <- colMeans(x)
    
    list_names <- c('effects.rf', 'effects.lm', 
                    'effects.f.true',  
                    'smaller_nulls')
    
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
      
      print(rf.predict$cov)
      print(rf.predict$cov[1,1]+rf.predict$cov[2,2] > 2*rf.predict$cov[1,2])
      
      lm.predict <- predict(linreg,
                            newdata = new_data,
                            se.fit = T)
      
      ### Computing Predictions
      
      # Using f_hat and sample mean&variance of X
      predictions.rf <- rowMeans(rf.predict$predictions)
      
      predictions.lm <- lm.predict$fit
      
      # Using f and true mean&variance of X
      predictions.f.true <- eval(formula, new_data.true)
      
      ### Computing Variable Effects
      
      # Using f_hat predictions based on sample mean&variance of X 
      result_list[['effects.rf']][v] <- (predictions.rf[2] - predictions.rf[1]) / (2*k*sd(x[,v]))
      result_list[['effects.lm']][v] <- (predictions.lm[2] - predictions.lm[1]) / (2*k*sd(x[,v]))
      
      # if (v==1) {
      #   print(c('RF', 'LM'))
      #   print(c(result_list[['effects.rf']][v], result_list[['effects.lm']][v]))
      #   cat("\n")
      # }
      
      # Using f predictions based on true mean&variance of X 
      result_list[['effects.f.true']][v] <- (predictions.f.true[2] - predictions.f.true[1]) / (2*k*true_s[v])
      
      
      result_list[['smaller_nulls']][v] <- sum((rf.predict$cov[1,1] + rf.predict$cov[2,2]
                                                - 2*rf.predict$cov[1,2]) /  (2*k*sd(x[,v]))^2 < 0)
      print(result_list[['smaller_nulls']][v])
    }
  
    
  }
  return(result_list)
}








### Repeat Simulation multiple times
sim_multi <- function(scenario){
  
  k <- scenario[["k"]]
  N <- scenario[["N"]]
  num.trees <- scenario[["N_Trees"]]
  formula <- scenario[["Formula"]]
  repeats <- scenario[["Repeats"]]
  cor <- scenario[["Correlation"]]
  
  all_terms <- extract_terms(formula)
  print(all_terms)
  
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
                       n.nulls = do.call('rbind', res[4,])
  )
  
  colnames(res_scenario$effects.rf) <- all_terms
  colnames(res_scenario$effects.lm) <- all_terms
  colnames(res_scenario$effects.f.true) <- all_terms
  
  return(res_scenario)
}




###### Simulation Setup
n <- c(2000) ; num.trees <- c(200) ; repeats <- 1; cor <- c(0); k=c(1)
#formulas <- c("2*x.1-3*x.2+4*x.3-0.5*x.4+2*x.5-3*x.6+4*x.7-0.5*x.8-2*x.9-3*x.10+4*x.11-0.5*x.12+2*x.13-3*x.14+4*x.15-0.5*x.16+2*x.17-3*x.18+4*x.19-0.5*x.20+2*x.21-3*x.22+4*x.23-0.5*x.24-2*x.25-3*x.26+4*x.27-0.5*x.28+2*x.29-3*x.30+4*x.31-0.5*x.32")
formulas <- c("2*x.1-3*x.2+4*x.3-0.5*x.4")
scenarios <- data.frame(expand.grid(n, num.trees, formulas, repeats, cor, k))
colnames(scenarios) = c("N", "N_Trees", "Formula", "Repeats", "Correlation", "k")
scenarios[,"Formula"] <- as.character(scenarios[,"Formula"]) ### Formula became Factor
scenarios <- split(scenarios, seq(nrow(scenarios)))

system.time(result <- lapply(X = scenarios, FUN = sim_multi))
colMeans(result$`1`$effects.rf)
colMeans(result$`1`$effects.lm)






