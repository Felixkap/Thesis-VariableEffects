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
    
    
    # result_list[['smaller_nulls']][v] <- sum((rf.predict$cov[1,1] + rf.predict$cov[2,2]
    #                                            - 2*rf.predict$cov[1,2]) /  (2*k*sd(x[,v]))^2 < 0)
    
    
    
    # Standard Error of Variable Effect (Coviariance Estimation)
    # VAR[(f(B) - f(A)) / 2*sd(x)]
    # =( 1 / 4*sd(x)^2 ) * ( VAR[f(B)] + VAR[f(B)] - 2*COV[f(A), f(B)] )
    effect.var <- pmax((rf.predict$cov[1,1] + rf.predict$cov[2,2]
                        - 2*rf.predict$cov[1,2]) /  (2*k*sd(x[,v]))^2 , 0)
    
    
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
  estimates <- list(effects.rf = do.call('rbind', res[1,]),
                    effects.lm = do.call('rbind', res[2,]),
                    #effects.f.sample = do.call('rbind', res[3,]),
                    #effects.rf.true = do.call('rbind', res[3,]),
                    effects.f.true = do.call('rbind', res[3,]),
                    effects.se.rf = do.call('rbind', res[4,]),
                    effects.se.lm = do.call('rbind', res[5,])
                    #smaller_nulls = do.call('rbind', res[7,])
  )
  
  
  res_scenario <- list(formula=formula, N=N, k=k, cor=cor, num.trees=num.trees,
                       effects.rf=estimates$effects.rf,
                       effects.lm=estimates$effects.lm,
                       #effects.f.sample=estimates$effects.f.sample,
                       #effects.rf.true=estimates$effects.rf.true,
                       effects.f.true=estimates$effects.f.true,
                       SE_Estimates.rf=estimates$effects.se.rf,
                       SE_Estimates.lm=estimates$effects.se.lm
                       #n.nulls=estimates$smaller_nulls
  )
  
  return(res_scenario)
}





print_results <- function(result){
  
  
  for (i in 1:length(result)) {
    cat('Setting: N =', result[[i]]$N, '; k =', result[[i]]$k, '; Correlation =', result[[i]]$cor,
        sprintf(paste0(';\n%', nchar('Setting: ')*2, 's'), 'Formula ='), result[[i]]$formula, '; N_Trees =', result[[i]]$num.trees,
        '\nMean(s) of simulated RF Variable Effect(s):\n ',
        colMeans(result[[i]]$effects.rf),
        '\nMean(s) of simulated LM Variable Effect(s):\n ',
        colMeans(result[[i]]$effects.lm),
        # '\nMean(s) of simulated RF Variable Effect(s) with true Means/Variances:\n ',
        # colMeans(result[[i]]$effects.rf.true),
        '\nMean(s) of True Variable Effect(s):\n ',
        colMeans(result[[i]]$effects.f.true),
        '\nStandard Error of simulated Variable Effects (RF):\n ',
        apply(result[[i]]$effects.rf, 2, sd),
        '.\nMean of Standard Errors Estimates of Variable Effects (RF):\n ',
        colMeans(result[[i]]$SE_Estimates.rf), '\n\n')
  }
  
}




plot_results <- function(result){
  
  effect_vec <- c('effects.rf', 
                  'effects.lm' 
                  #'effects.rf.true'
  )
  
  mylist <- vector(mode = 'list', length = length(result)*3)
  cnt <- 1 ; max_vars <- 0
  for (res_i in 1:length(result)) {
    for (j in effect_vec) {
      
      current_var_names <- unique(unlist(str_extract_all(result[[res_i]][['formula']],"x.[0-9]")))
      if (length(current_var_names) > max_vars) {
        var_names <- current_var_names
        max_vars <- length(var_names) 
      }
      

      dat_tmp <- data.frame(x=result[[res_i]][[j]])
      dat <- cbind(dat_tmp, 
                   N = paste('N=', result[[res_i]]$N),
                   cor = paste('Cor=', result[[res_i]]$cor),
                   k = paste('k=', result[[res_i]]$k),
                   names(result[[res_i]][j]))
      names(dat)[length(names(dat))] <- 'effect.type'
  

      mylist[[cnt]] <-  dat %>%
        pivot_longer(cols = starts_with('x.'), names_to = 'variable')
      cnt <- cnt + 1
      
    }
  }
  
  mydat <- do.call('rbind', mylist) 
  
  mydat$effect.type <- as.factor(mydat$effect.type)
  
  levels(mydat$effect.type) <- c(effects.lm='LM', 
                                 effects.rf='RF'
                                 #effects.rf.true='RF.true'
  )
  
  
  df2 <- data.frame(variable = var_names, value = result[[1]]$effects.f.true[1,])
  
  
  plot_result <- ggplot(mydat, aes(x = value, y=effect.type)) +
    ggh4x::facet_grid2(variable~N+cor+k, scales="free_x", independent = "x")+
    geom_boxplot(aes(fill=effect.type)) + 
    geom_vline(data=df2, aes(xintercept=value, color='TrueEffect')) + 
    labs(y = 'Model', x = 'Effect Estimates', fill = "Model") + 
    theme(axis.text.x = element_text(angle = 45)) +
    scale_color_manual(name = "", values = c(TrueEffect = "orange"))
  
  return(plot_result)
}

# ## Simulation Setup
# n <- c(100) ; num.trees <- 20 ; repeats <- 100; cor <- c(0); k <- c(0.5)
# formulas <- c("2*x.1+4*x.2-0.5*x.3") #"-0.5*x.1^3+3*x.2+0.5*sqrt(abs(x.3))"
# scenarios <- data.frame(expand.grid(n, num.trees, formulas, repeats, cor, k))
# colnames(scenarios) = c("N", "N_Trees", "Formula", "Repeats", "Correlation", "k")
# scenarios[,"Formula"] <- as.character(scenarios[,"Formula"]) ### Formula became Factor
# scenarios <- split(scenarios, seq(nrow(scenarios)))
# 
# system.time(result <- lapply(X = scenarios, FUN = sim_multi))
# print_results(result)
# 
# result_plots <- plot_results(result)
# result_plots





