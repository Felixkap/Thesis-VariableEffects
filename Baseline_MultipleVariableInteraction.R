rm(list=ls())
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
# Fit Random Forest n times
# Compute Variable Effect (as we defined it) n times
# Estimate Standard Error of Variable Effect n times directly
# ---> Based on estimated Variance-Covariance Matrix

# Compare Standard Error of simulated Variable Effects with mean of estimated Standard Errors of Variable Effects



sim_once <- function(N, num.trees, cor, k, formula){
  
  
  ### Simulate Data
  n_vars <- length(unique(unlist(str_extract_all(formula,"x.[0-9]")))) # Number of variables
  n_coefs <- length(unlist(strsplit(formulas, "(?<=.)(?=[+])|(?<=.)(?=[-])",perl = TRUE)))
  if (n_vars == 1) {
    # Normal Distribution with Mean 0 and Variance 1
    x <- data.frame(x.1 = rnorm(N, 0, 1))
  } else{
    # Multivariate Normal Distribution with Mean Vector 0 and diagonal Variance Covariance Matrix 
    x <- data.frame(x = rmvnorm(N, mean = rep(0, n_vars), sigma = CSgenerate(n_vars, cor)))
  }
  
  
  # Compute interaction terms
  all_terms <- unlist(strsplit(formulas, "(?<=.)(?=[+])|(?<=.)(?=[-])",perl = TRUE))
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
  
  linreg <- lm( formula = y ~ (.)^2,
                data = data)

  
  
  ### Variable Effect & Standard Error
  
  # True Variable Means and Standard Deviations
  true_mean <- rep(x = 0, times = n_vars)
  names(true_mean) <- unique(unlist(str_extract_all(formula,"x.[0-9]")))
  true_s <- rep(x = 1, times = n_vars)
  names(true_s) <- unique(unlist(str_extract_all(formula,"x.[0-9]")))
  
  
  # Estimated Variable Means
  x_bar <- colMeans(x)

  
  list_names <- c('effects.rf', 'effects.lm', 'effects.f.true', 
                  'effects.se.rf', 'effects.se.lm', 'smaller_nulls')
  result_list <- vector(mode = 'list', length = length(list_names))

  result_list <- lapply(X = 1:length(list_names), FUN = function(X){result_list[[X]] = numeric(length = n_coefs)})
  names(result_list) <- list_names

  
  
  for (v in 1:n_vars) {
    
    
    if (v == as.numeric(gsub(".*?([0-9]+).*", "\\1", interaction_term[1]))) {
      
      # New data points based on estimated variable means and standard deviations
      x_a <- replace(x_bar, v, x_bar[v] - k*sd(x[,v]))
      x_b <- replace(x_bar, v, x_bar[v] + k*sd(x[,v]))
      
      
      # New data points based on true variable means and standard deviations
      x_a.true <- replace(true_mean, v, true_mean[v] - k*true_s[v])
      x_b.true <- replace(true_mean, v, true_mean[v] + k*true_s[v])

      x_a1 <- replace(x_a, n_vars, x_bar[v] - k*sd(x[,v]))
      x_b1 <- replace(x_b, n_vars, x_bar[v] - k*sd(x[,v]))
      x_a2 <- replace(x_a, n_vars, x_bar[v] + k*sd(x[,v]))
      x_b2 <- replace(x_b, n_vars, x_bar[v] + k*sd(x[,v]))

    
      x_a1.true <- replace(x_a.true, n_vars, true_mean[v] - k*true_s[v])
      x_b1.true <- replace(x_b.true, n_vars, true_mean[v] - k*true_s[v])
      x_a2.true <- replace(x_a.true, n_vars, true_mean[v] + k*true_s[v])
      x_b2.true <- replace(x_b.true, n_vars, true_mean[v] + k*true_s[v])
     
      new_data <- data.frame(rbind(x_a, x_b, 
                                   x_a1, x_b1, 
                                   x_a2, x_b2))
      new_data.true <- data.frame(rbind(x_a.true, x_b.true, 
                                        x_a1.true, x_b1.true, 
                                        x_a2.true, x_b2.true))


    } else {
      # New data points based on estimated variable means and standard deviations
      x_a <- replace(x_bar, v, x_bar[v] - k*sd(x[,v]))
      x_b <- replace(x_bar, v, x_bar[v] + k*sd(x[,v]))
      
      
      # New data points based on true variable means and standard deviations
      x_a.true <- replace(true_mean, v, true_mean[v] - k*true_s[v])
      x_b.true <- replace(true_mean, v, true_mean[v] + k*true_s[v])
      
      new_data <- data.frame(rbind(x_a, x_b))
      new_data.true <- data.frame(rbind(x_a.true, x_b.true))
    }
    
    
    rf.predict <- RangerForestPredict(rf$forest,
                                      new_data,
                                      type='se',
                                      se.method = 'jack_cov',
                                      predict.all = T,
                                      inbag.counts = rf$inbag.counts)
    
    lm.predict <- predict(linreg,
                          newdata = new_data,
                          se.fit = T)
    
    
    
    ### Computing Predictions
    
    # Using f_hat and sample mean&variance of X
    predictions.rf <- rowMeans(rf.predict$predictions)
    predictions.lm <- lm.predict$fit
    
    # Using f and sample mean&variance of X
    #predictions.f.sample <- eval(formula, new_data)
    
    # Using f and true mean&variance of X
    predictions.f.true <- eval(formula, new_data.true)

    
    
    ### Computing Variable Effects
    
    # Using f_hat predictions based on sample mean&variance of X 
    

    result_list[['effects.rf']][v] <- (predictions.rf[2] - predictions.rf[1]) / (2*k*sd(x[,v]))
    result_list[['effects.lm']][v] <- (predictions.lm[2] - predictions.lm[1]) / (2*k*sd(x[,v]))
    # Using f predictions based on sample mean&variance of X 
    result_list[['effects.f.sample']][v] <- (predictions.f.sample[2] - predictions.f.sample[1]) / (2*k*sd(x[,v]))
    
    # Using f predictions based on true mean&variance of X 
    result_list[['effects.f.true']][v] <- (predictions.f.true[2] - predictions.f.true[1]) / (2*k*true_s[v])
    
    
    # result_list[['smaller_nulls']][v] <- sum((rf.predict$cov[1,1] + rf.predict$cov[2,2]
    #                                           - 2*rf.predict$cov[1,2]) /  (2*k*sd(x[,v]))^2 < 0)
    
    # Standard Error of Variable Effect (Coviariance Estimation)
    # VAR[(f(B) - f(A)) / 2*sd(x)]
    # =( 1 / 4*sd(x)^2 ) * ( VAR[f(B)] + VAR[f(B)] - 2*COV[f(A), f(B)] )
    effect.var <- pmax((rf.predict$cov[1,1] + rf.predict$cov[2,2]
                        - 2*rf.predict$cov[1,2]) /  (2*k*sd(x[,v]))^2 , 0)
    
    
    result_list[['effects.var.rf']][v] <- effect.var
    result_list[['effects.var.lm']][v] <- (summary(linreg)$coef[,"Std. Error"][v+1])^2
    
    # Interaction Effect
    if (v == as.numeric(gsub(".*?([0-9]+).*", "\\1", interaction_term[1]))) {
      result_list[['effects.rf']][n_coefs] <- ((predictions.rf[6] - predictions.rf[4]) - (predictions.rf[5] - predictions.rf[3])) / (4*k^2*sd(x[,interaction_term[1]])*sd(x[,interaction_term[2]]))
      result_list[['effects.lm']][n_coefs] <- ((predictions.lm[6] - predictions.lm[4]) - (predictions.lm[5] - predictions.lm[3])) / (4*k^2*sd(x[,interaction_term[1]])*sd(x[,interaction_term[2]]))
      # Using f predictions based on sample mean&variance of X 
      #result_list[['effects.f.sample']][n_coefs] <- ((predictions.f.sample[6] - predictions.f.sample[4]) - (predictions.f.sample[5] - predictions.f.sample[3])) / (2*k*sd(x[,v]))
      
      # Using f predictions based on true mean&variance of X 
      result_list[['effects.f.true']][n_coefs] <- ((predictions.f.true[6] - predictions.f.true[4]) - (predictions.f.true[5] - predictions.f.true[3])) / (4*k^2*true_s[interaction_term[1]]*true_s[interaction_term[2]])
      
      # result_list[['smaller_nulls']][n_coefs] <- sum((sum(diag(rf.predict$cov))
      #                                                 - 2*rf.predict$cov[2,1]
      #                                                 - 2*rf.predict$cov[3,1] + 2*rf.predict$cov[3,2]
      #                                                 + 2*rf.predict$cov[4,1] - 2*rf.predict$cov[4,2] - 2*rf.predict$cov[4,3]) /  (2*sd(x[,v]))^2 < 0)
      
      effect.var <- pmax((sum(diag(rf.predict$cov))
                          - 2*rf.predict$cov[2,1]
                          - 2*rf.predict$cov[3,1] + 2*rf.predict$cov[3,2]
                          + 2*rf.predict$cov[4,1] - 2*rf.predict$cov[4,2] - 2*rf.predict$cov[4,3]) /  ((4*k^2*sd(x[,interaction_term[1]])*sd(x[,interaction_term[2]])))^2, 0)
      
      result_list[['effects.var.rf']][n_coefs] <- effect.var
      result_list[['effects.var.lm']][n_coefs] <- (summary(linreg)$coef[,"Std. Error"][v+1])^2
      
    }
    
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
        sqrt(colMeans(result[[i]]$SE_Estimates.rf)), '\n\n')
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


##### Simulation Setup
n <- c(10000) ; num.trees <- c(20) ; repeats <- 1e3; cor <- 0; k=0.3
formulas <- c("2*x.1+4*x.2-0.5*x.3+2*x.1*x.3") #"-0.5*x.1^3+3*x.2+0.5*sqrt(abs(x.3))"
scenarios <- data.frame(expand.grid(n, num.trees, formulas, repeats, cor))
colnames(scenarios) = c("N", "N_Trees", "Formula", "Repeats", "Correlation")
scenarios[,"Formula"] <- as.character(scenarios[,"Formula"]) ### Formula became Factor
scenarios <- split(scenarios, seq(nrow(scenarios)))


# Run Simulation
system.time(result <- lapply(X = scenarios, FUN = sim_multi))

print_results(result)
plot_results(result)






