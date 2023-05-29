library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ranger)
library(MixMatrix)
library(mvtnorm)
library(stringr)
library(ggh4x)
library(iml)


source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/RangerPredictFunction.R')
source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/Simulation - Variable Effects & SE/Helper.R')

##### Simulation --- Distribution of Variable Effect and Jackknife-after-Bootstrap Standard Errors
# Generate REPEATS data sets
# Fit Random Forest & Linear Regression Model on each data set
# Compute all Variable Main and Low-Order Interaction Effects
# Estimate Standard Error of Variable Effects directly by means of Jackknife-after-Bootstrap Estimation

# Compare Standard Error of simulated Variable Effects with mean of Jackknife Standard Errors






sim_once <- function(k, N, num.trees, cor, formula, node_size, k_idx, sigma_e){
  
  ## @param k: numeric. Range over which effects of numerical variables are computed
  ## @param N: Integer. Number of observations
  ## @param cor: Numeric Value. Correlation between independent variables
  ## @param formula: String. Functional Relationship between predictors and outcome
  ## @param node_size: Integer: Minimum Number of observations in terminal leaves
  ## @param sigma_e: Numeric Value. Noise in the data (standard deviation of errors)
  
  
  ### Simulate Data
  n_vars <- length(unique(unlist(str_extract_all(formula,"x.[0-9]+")))) # Number of variables
  
  # Get all individual terms
  all_terms <- extract_terms(formula)
  
  # Number of effect terms
  n_coefs <- length(all_terms)
  
  # Get interaction terms
  all_interactions <- get_all_interactions(formula)
  
  # Generating Data Set
  train_data <- generate_data(N, n_vars, cor, formula, sigma_e)
  x <- train_data$x
  
  # Variable Names
  x_names <- names(x)
  
  # Predictor Names
  predictor_names <- c(x_names, all_interactions)
  
  # Estimated Variable Means 
  x_bar <- train_data$x_bar
  
  # True Variable Means and Standard Deviations
  true_mean <- rep(x = 0, times = n_vars)
  names(true_mean) <- unique(unlist(str_extract_all(formula,"x.[0-9]+")))
  true_s <- rep(x = 1, times = n_vars)
  names(true_s) <- unique(unlist(str_extract_all(formula,"x.[0-9]+")))
  
  # Random Forest
  rf <- ranger( formula = y ~ .,
                data = train_data$data,
                num.trees = num.trees,
                keep.inbag = T, 
                oob.error = T, # Save computational time
                min.bucket = node_size)
  
  # Linear Model
  linreg <- lm( formula = y ~ .^2,
                data = train_data$data)
  

  
  
  ### Setup result list containing components defined in 'list_names'
  
  list_names <- c('effects.rf', 'effects.lm', 
                  'effects.f.true', 
                  'effects.se.rf', 
                  'smaller_nulls')
  
  result_list <- vector(mode = 'list', length = length(list_names))
  result_list <- lapply(X = 1:length(list_names), FUN = function(X){
    tmp_vec <- numeric(length = n_coefs)
    names(tmp_vec)<- predictor_names
    result_list[[X]] = tmp_vec
    })
  names(result_list) <- list_names
  

  for (v in 1:n_vars) {
    # New data points based on estimated variable means and standard deviations
    
    # New data points based on estimated variable means and standard deviations
    x_a <- replace(x_bar, v, x_bar[v] - k*sd(x[,v]))
    x_b <- replace(x_bar, v, x_bar[v] + k*sd(x[,v]))
    
    # New data points based on true variable means and standard deviations
    x_a.true <- replace(true_mean, v, true_mean[v] - k*true_s[v])
    x_b.true <- replace(true_mean, v, true_mean[v] + k*true_s[v])
    
    new_data <- data.frame(rbind(x_a, x_b))
    new_data.true <- data.frame(rbind(x_a.true, x_b.true))
      
    # Predictions of Random Forest
    rf.predict <- RangerForestPredict(rf$forest,
                                      new_data,
                                      type='se',
                                      se.method = 'se_direct',
                                      num.trees = rf$num.trees,
                                      k_ratio = k, 
                                      save = FALSE,
                                      predict.all = T,
                                      inbag.counts = rf$inbag.counts)
    
    lm.predict <- predict(linreg,
                          newdata = new_data,
                          se.fit = T)
    
    ### Computing Variable Effects
    predictions.rf <- rowMeans(rf.predict$predictions)
    
    predictions.lm <- lm.predict$fit
    
    predictions.f.true <- eval(parse(text = formula), new_data.true)

    
    # Using f_hat predictions based on sample mean&variance of X 
    result_list[['effects.rf']][v] <- (predictions.rf[2] - predictions.rf[1]) / (2*k*sd(x[,v]))
    result_list[['effects.lm']][v] <- (predictions.lm[2] - predictions.lm[1]) / (2*k*sd(x[,v]))
    
    # Using f predictions based on true mean&variance of X 
    result_list[['effects.f.true']][v] <- (predictions.f.true[2] - predictions.f.true[1]) / (2*k*true_s[v])
    
    
    effect.se.direct <- rf.predict$se_effect
    result_list[['smaller_nulls']][v] <- sum(rf.predict$se_effect == 0)
    

    result_list[['effects.se.rf']][v] <- effect.se.direct
    
  
    
    if (exists('all_interactions')) {
      
      v_cnt <- n_vars + 1
      
      for (interaction_term in all_interactions) {
        
        v_idx <-  which(x_names %in% str_split(interaction_term, ':')[[1]])
        
        # New data points needed for low-order interaction effect
        x_a <- replace(x_bar, v_idx[1], x_bar[v_idx[1]] - k*sd(x[,v_idx[1]]))
        x_b <- replace(x_bar, v_idx[1], x_bar[v_idx[1]] + k*sd(x[,v_idx[1]]))
        
        # New data points based on true variable means and standard deviations
        x_a.true <- replace(true_mean, v_idx[1], true_mean[v_idx[1]] - k*true_s[v_idx[1]])
        x_b.true <- replace(true_mean, v_idx[1], true_mean[v_idx[1]] + k*true_s[v_idx[1]])
        
        x_a1 <- replace(x_a, v_idx[2], x_bar[v_idx[2]] - k*sd(x[,v_idx[2]]))
        x_b1 <- replace(x_b, v_idx[2], x_bar[v_idx[2]] - k*sd(x[,v_idx[2]]))
        x_a2 <- replace(x_a, v_idx[2], x_bar[v_idx[2]] + k*sd(x[,v_idx[2]]))
        x_b2 <- replace(x_b, v_idx[2], x_bar[v_idx[2]] + k*sd(x[,v_idx[2]]))
        
        x_a1.true <- replace(x_a.true, v_idx[2], true_mean[v_idx[2]] - k*true_s[v_idx[2]])
        x_b1.true <- replace(x_b.true, v_idx[2], true_mean[v_idx[2]] - k*true_s[v_idx[2]])
        x_a2.true <- replace(x_a.true, v_idx[2], true_mean[v_idx[2]] + k*true_s[v_idx[2]])
        x_b2.true <- replace(x_b.true, v_idx[2], true_mean[v_idx[2]] + k*true_s[v_idx[2]])
        
        new_data <- data.frame(rbind(x_a1, x_b1, 
                                    x_a2, x_b2))
        
        new_data.true <- data.frame(rbind(x_a1.true, x_b1.true, 
                                          x_a2.true, x_b2.true))
        
        # Predictions of Random Forest
        rf.predict <- RangerForestPredict(rf$forest,
                                          new_data,
                                          type='se',
                                          k_ratio = k, 
                                          num.trees = rf$num.trees,
                                          save = FALSE,
                                          se.method = 'se_direct',
                                          predict.all = T,
                                          inbag.counts = rf$inbag.counts)
        lm.predict <- predict(linreg,
                              newdata = new_data,
                              se.fit = T)
        
        ### Computing low-order interaction Effect
        predictions.rf <- rowMeans(rf.predict$predictions)
        predictions.lm <- lm.predict$fit
        predictions.f.true <- eval(parse(text = formula), new_data.true)
        
        result_list[['effects.rf']][v_cnt] <- ((predictions.rf[3] - predictions.rf[1]) - (predictions.rf[4] - predictions.rf[2])) / (-4*k^2*sd(x[,v_idx[1]])*sd(x[,v_idx[2]]))
        result_list[['effects.lm']][v_cnt] <- ((predictions.lm[3] - predictions.lm[1]) - (predictions.lm[4] - predictions.lm[2])) / (-4*k^2*sd(x[,v_idx[1]])*sd(x[,v_idx[2]]))
        result_list[['effects.f.true']][v_cnt] <- ((predictions.f.true[3] - predictions.f.true[1]) - (predictions.f.true[4] - predictions.f.true[2])) / (-4*k^2*true_s[v_idx[1]]*true_s[v_idx[2]])
        
        effect.se.direct <- rf.predict$se_effect
        result_list[['smaller_nulls']][v_cnt] <- sum(rf.predict$se_effect == 0)
        
        result_list[['effects.se.rf']][v_cnt] <- effect.se.direct 
        
        v_cnt <- v_cnt + 1 
        
      }
    }
  }
  
  return(result_list)
}







### Repeat Simulation multiple times
sim_multi <- function(scenario){
  
  ## @param scenario: DataFrame containing one row which defines settings of scenario
  
  k <- scenario[["k"]]
  N <- scenario[["N"]]
  num.trees <- scenario[["N_Trees"]]
  formula <- scenario[["Formula"]]
  longest_latex_formula <- scenario[["Longest_Latex_formula"]]
  repeats <- scenario[["Repeats"]]
  cor <- scenario[["Correlation"]]
  node_size <- scenario[["Node_Size"]]
  k_idx <- scenario[["k_idx"]]
  sigma_e <- scenario[["sigma_e"]]
  
  all_terms <- extract_terms(formula)
  
  res <- replicate(n = repeats,
                   expr = sim_once(k=k, N = N, num.trees = num.trees,
                                   cor = cor, formula = formula, 
                                   node_size = node_size,
                                   k_idx = k_idx, sigma_e = sigma_e))
  
  
  
  ### Distribution of Test Statistics
  res_scenario <- list(formula=formula, longest_latex_formula=longest_latex_formula, 
                       N=N, k=k, cor=cor, sigma_e = sigma_e,
                       num.trees=num.trees, node_size=node_size,
                       effects.rf = do.call('rbind', res[1,]),
                       effects.lm = do.call('rbind', res[2,]),
                       effects.f.true = do.call('rbind', res[3,]),
                       effects.se.rf = do.call('rbind', res[4,]),
                       n.nulls = do.call('rbind', res[5,])
  )
  
  return(res_scenario)
}



