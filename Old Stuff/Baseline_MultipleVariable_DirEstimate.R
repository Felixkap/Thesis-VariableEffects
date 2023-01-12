library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ranger)
library(MixMatrix)
library(mvtnorm)
library(stringr)

source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/RangerPredictFunction.R')


##### Simulation --- Distribution of Variable Effect and Standard Error
# Fit Random Forest n times
# Compute Variable Effect (as we defined it) n times
# Estimate Standard Error of Variable Effect n times directly
# ---> Based on estimated Variance-Covariance Matrix

# Compare Standard Error of simulated Variable Effects with mean of estimated Standard Errors of Variable Effects



### Helper function to create Plots
create_longformat <- function(result){
  
  n_scenarios <- length(result)
  effect_estimates_long_list <- vector(mode = 'list', length = n_scenarios)
  for (s in 1:n_scenarios) {
    
    if (ncol(result[[s]]$Effect_Estimates)==1) {
      
      effect_estimates_list <- data.frame(effect.x.1=result[[s]]$Effect_Estimates,
                                          se.cov.x.1=result[[s]]$SE_Estimates,
                                          se.direct.x.1=result[[s]]$SE_Estimates2,
                                          scenario=as.character(s))
      
    } else {
      
      effect_estimates_list <- data.frame(effect.x=result[[s]]$Effect_Estimates,
                                          se.cov.x=result[[s]]$SE_Estimates,
                                          se.direct.x=result[[s]]$SE_Estimates2,
                                          scenario=as.character(s))
    }
    
    
    
    effect_estimates_long_list[[s]] <- effect_estimates_list %>% 
      pivot_longer(cols = -one_of('scenario'), names_to = 'metric')
    
  }
  
  effect_estimates <- do.call('rbind', effect_estimates_long_list)
  
  
  return(effect_estimates)
  
}



split_formula <- function(form){
  
  form1 <- substr(form, start = 1, stop = ceiling(nchar(form)/2))
  form2 <- substr(form, start = ceiling(nchar(form)/2), stop = nchar(form))
  splitted_form <- paste0(form1, '\n', form2)
  return(splitted_form)
  
}




sim_once <- function(N, num.trees, cor, formula){
  
  
  n_vars <- length(unique(unlist(str_extract_all(formula,"x.[0-9]")))) # Number of variables
  if (n_vars == 1) {
    x <- data.frame(x.1 = rnorm(N, 0, 1))
  } else{
    print(CSgenerate(n_vars, cor))
    x <- data.frame(x = rmvnorm(N, mean = rep(0, n_vars), sigma = CSgenerate(n_vars, cor)))
  }
  
  e <- rnorm(N, 0, 1)
  xe <- data.frame(cbind(x,e))
  formula <- parse(text = formula)
  y <- eval(formula, xe)
  data <- data.frame(cbind(x, y))
  
  
  rf <- ranger( formula = y ~ ., 
                data = data, 
                num.trees = num.trees, 
                keep.inbag = T)
  
  
  
  
  ### Variable Effect & Standard Error
  x_bar <- colMeans(x)
  effects <- numeric(length = n_vars)
  effects.se <- numeric(length = n_vars)
  effects.se2 <- numeric(length = n_vars)
  smaller_nulls <- numeric(length = n_vars)
  
  for (v in 1:n_vars) {
    
    x_a <- replace(x_bar, v, x_bar[v] - sd(x[,v]))
    x_b <- replace(x_bar, v, x_bar[v] + sd(x[,v]))
    
    new_data <- data.frame(rbind(x_a, x_b))
    
    rf.predict <- RangerForestPredict(rf$forest, 
                                      new_data, 
                                      type='se', 
                                      se.method = 'jack_cov',
                                      predict.all = T,
                                      inbag.counts = rf$inbag.counts)
    
    
    
    # Variable Effect
    ab_predictions <- rowMeans(rf.predict$predictions)
    effects[v] <- (ab_predictions[2] - ab_predictions[1]) / 2*sd(x[,v]) 
    
    
    smaller_nulls[v] <- sum((rf.predict$cov[1,1] + rf.predict$cov[2,2] 
                             - 2*rf.predict$cov[1,2]) /  (2*sd(x[,v]))^2 < 0)
    
    
    # Standard Error of Variable Effect (Coviariance Estimation)
    # VAR[(f(B) - f(A)) / 2*sd(x)] 
    # =( 1 / 4*sd(x)^2 ) * ( VAR[f(B)] + VAR[f(B)] - 2*COV[f(A), f(B)] )
    effect.var <- pmax((rf.predict$cov[1,1] + rf.predict$cov[2,2] 
                        - 2*rf.predict$cov[1,2]) /  (2*sd(x[,v]))^2 , 0)
    
    
    
    effects.se[v] <- sqrt(effect.var)
    
    # Direct Standard Error Estimate of Variable Effect 
    effects.se2[v] <- rf.predict$se_effect
    
  }
  
  
  return(list(effect = effects, 
              effect.se = effects.se,
              effect.se2 = effects.se2,
              smaller_null = smaller_nulls))
  
}






sim_multi <- function(scenario){
  
  
  N <- scenario[["N"]]
  num.trees <- scenario[["N_Trees"]]
  formula <- scenario[["Formula"]]
  repeats <- scenario[["Repeats"]]
  cor <- scenario[["Correlation"]]
  
  res <- replicate(n = repeats,
                   expr = sim_once(N = N, num.trees = num.trees, 
                                   cor = cor, formula = formula))
  
  ### Distribution of Test Statistics
  estimates <- list(effects = do.call('rbind', res[1,]), 
                    effects.se = do.call('rbind', res[2,]),
                    effects.se2 = do.call('rbind', res[3,]),
                    smaller_nulls = do.call('rbind', res[4,]))
  
  
  
  res_scenario <- list(formula=formula, N=N, cor=cor, num.trees=num.trees, 
                       Effect_Estimates=estimates$effects, 
                       SE_Estimates=estimates$effects.se, 
                       n.nulls=estimates$smaller_nulls,
                       SE_Estimates2=estimates$effects.se2)
  
  return(res_scenario)
}







plot_results <- function(result){
  
  sim_results_long <- create_longformat(result)
  
  
  n_scenarios <- length(result)
  scenario_names <- character(n_scenarios)
  for (i in 1:n_scenarios) {
    if (nchar(result[[i]]$formula)> 25) {
      scenario_names[i] <- split_formula(result[[i]]$formula)
    } else {
      scenario_names[i] <- result[[i]]$formula
    }
  }
  
  names(scenario_names) <- 1:n_scenarios
  
  
  max_vars <- 0
  for (res in result) {
    n_vars <- length(unique(unlist(str_extract_all(res$formula,"x.[0-9]"))))
    if (n_vars > max_vars) {
      max_vars = n_vars
    }
  }
  
  effect_names <- character(length = max_vars)
  se_cov_names <- character(length = max_vars)
  se_direct_names <- character(length = max_vars)
  for (i in 1:max_vars) {
    effect_names[i] <- paste0('X', i)
    se_cov_names[i] <- paste0('X', i)
    se_direct_names[i] <- paste0('X', i)
  }
  
  names(effect_names) <- paste0('effect.x.', 1:max_vars)
  names(se_cov_names) <- paste0('se.cov.x.', 1:max_vars)
  names(se_direct_names) <- paste0('se.direct.x.', 1:max_vars)
  
  
  effect_effect_plot <- subset(sim_results_long, grepl("^effect", metric), drop = TRUE) %>%
    ggplot(aes(x = value)) +
    facet_grid(scenario~metric, scale="free", 
               labeller = labeller(scenario = as_labeller(scenario_names),
                                   metric = as_labeller(effect_names)))+
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.05, color = 'black') +
    geom_density(alpha = .2, fill = "#FF6666") +
    labs(x = 'Variable Effect Estimates',
         title = 'Distribution of Variable Effect Estimates')+
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text = element_text(size = 8))
  
  
  effect_se_plot <- subset(sim_results_long, grepl("^se.cov", metric), drop = TRUE) %>%
    ggplot(aes(x = value)) +
    facet_grid(scenario~metric, scale="free", 
               labeller = labeller(scenario = as_labeller(scenario_names),
                                   metric = as_labeller(se_cov_names)))+
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.05, color = 'black') +
    geom_density(alpha = .2, fill = "#FF6666") +
    labs(x = 'Standard Error Estimates',
         title = 'Distribution of (Cov) Standard Error Estimates of Variable Effects')+
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text = element_text(size = 8))
  
  effect_se2_plot <- subset(sim_results_long, grepl("^se.direct", metric), drop = TRUE) %>%
    ggplot(aes(x = value)) +
    facet_grid(scenario~metric, scale="free", 
               labeller = labeller(scenario = as_labeller(scenario_names),
                                   metric = as_labeller(se_direct_names)))+
    geom_histogram(aes(y = after_stat(density)), binwidth = 0.05, color = 'black') +
    geom_density(alpha = .2, fill = "#FF6666") +
    labs(x = 'Standard Error Estimates',
         title = 'Distribution of (direct) Standard Error Estimates of Variable Effects')+
    theme(strip.background=element_rect(fill=NA, color=NA),
          strip.text = element_text(size = 8))
  
  
  effect_boxplot <- subset(sim_results_long, grepl("^effect", metric), drop = TRUE) %>%
    ggplot(aes(x = scenario, y = value, color = scenario)) +
    facet_wrap(~metric, scale="free", 
               labeller = labeller(metric = as_labeller(effect_names)))+
    geom_boxplot() +
    labs(x = 'Scenarios',
         title = 'Distribution of Variable Effect Estimates')+
    theme(strip.background=element_rect(fill=NA, color=NA))
  
  
  effect_se_boxplot <- subset(sim_results_long, !grepl("^effect", metric), drop = TRUE) %>%
    mutate(se_cov=if_any(metric, function(x){grepl('cov',x)})) %>%
    mutate(variable=str_extract(metric, "x.[0-9]")) %>%
    ggplot(aes(x = scenario, y = value, color = scenario)) +
    facet_grid(se_cov~variable, labeller = labeller(se_cov = as_labeller(c(`TRUE` = 'SE Cov',
                                                                           `FALSE` = 'SE Direct'))))+
    geom_boxplot() +
    labs(x = 'Scenarios',
         y = 'Standard Error Estimate',
         title = 'Distribution of Standard Error Estimates of Variable Effects')+
    theme(strip.background=element_rect(fill=NA, color=NA))
  
  
  
  return(list(effect = effect_effect_plot,
              se = effect_se_plot,
              se2 = effect_se2_plot,
              effect_box = effect_boxplot,
              se_box = effect_se_boxplot))
  
}




print_results <- function(result){
  
  
  for (i in 1:length(result)) {
    cat('Setting: N=', result[[i]]$N, 
        '; Correlation=', result[[i]]$cor,
        '; Formula=', result[[i]]$formula, 
        '; N_Trees=', result[[i]]$num.trees,
        '\nMean(s) of simulated Variable Effect(s):\n ',  
        colMeans(result[[i]]$Effect_Estimates),
        '\nStandard Error of simulated Variable Effects:\n ', 
        apply(result[[i]]$Effect_Estimates, 2, sd), 
        '.\nMean of (Cov) Standard Errors Estimates of Variable Effects:\n ', 
        colMeans(result[[i]]$SE_Estimates), 
        '.\nMean of (direct) Standard Errors Estimates of Variable Effects:\n ', 
        colMeans(result[[i]]$SE_Estimates2),
        '.\nNumber of smaller Nulls\n', colSums(result[[i]]$n.nulls), '\n\n')
  }
  
}






###### Simulation Setup
#n <- 200 ; num.trees <- 2000 ; repeats <- 1e3; cor <- 0.8
#formulas <- c("2*x.1^3+2*x.2^2-0.5*x.3^3+0.2*x.4^4+e",
#              "-0.5*x.1^3+3*x.2+x.3^2+0.5*sqrt(abs(x.4))+e")
# scenarios <- data.frame(expand.grid(n, num.trees, formulas, repeats, cor))
# colnames(scenarios) = c("N", "N_Trees", "Formula", "Repeats", "Correlation")
# scenarios[,"Formula"] <- as.character(scenarios[,"Formula"]) ### Formula became Factor
# scenarios <- split(scenarios, seq(nrow(scenarios)))
# 
# 
# # Run Simulation
# system.time(result <- lapply(X = scenarios, FUN = sim_multi))
# 
# print_results(result)
# plot_results(result)

