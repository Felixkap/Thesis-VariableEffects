library(Matrix)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ranger)
library(MixMatrix)
library(mvtnorm)
library(stringr)
library(ggh4x)
library(latex2exp)
library(iml)
library(ALEPlot)

source('C:/Users/kpl/Thesis Project/Thesis-VariableEffects/RangerPredictFunction.R')


################################################################################
################################################################################
########################       HELPER FUNCTIONS      ###########################
################################################################################
################################################################################


plot_effects <- function(result){
  
  scenario_list <- data_per_scenario(var_name = 'x.', 
                                     effect_types = c('effects.rf'), 
                                     result = result)
  
  effects_df <- get_true_effects(result)
  alldat <- do.call('rbind', scenario_list) 
  
  alldat$effect.type <- as.factor(alldat$effect.type)
  
  if (sum(sapply(alldat[['variable']], function(x){str_count(x)!=3})) == 0) {
    title = 'Estimating Variable Main Effects (for a given Data Set)'
  } else {
    title = 'Estimating Variable Main and Interaction Effects (for a given Data Set)'
  }
  
  levels(alldat$effect.type) <- c(effects.rf='RF')
  
  
  plot_result <- ggplot(alldat, aes(y=value)) +
    ggh4x::facet_grid2(get_label_formula(alldat), axes= "all", scales="free")+
    geom_boxplot(aes(fill=variable)) + 
    geom_hline(data=effects_df, aes(yintercept=mean_val, color='TrueEffect')) + 
    labs(y = 'Effect Estimates', x = 'Scenarios', fill = "Variable", 
         title = title,
         caption = TeX(paste0('Remaining Settings: ', get_caption(result, alldat), '$'))) + 
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_line(),
          axis.title = element_text(size=18),
          panel.background = element_blank(),
          strip.text.x = element_text(size = 9),
          strip.text.y = element_text(size = 12),
          legend.title = element_text(size=15), #change legend title font size
          legend.text = element_text(size=15),
          plot.title = element_text(hjust = 0.5, lineheight = 0.9, size = 26),
          plot.caption = element_text(hjust = 0, size = 15)) +
    scale_color_manual(name = "", values = c(TrueEffect = "orange"))
  
  return(plot_result)
}



### Plot distribution of Variable Effect Standard Errors
plot_se_box <- function(result){
  
  scenario_list <- data_per_scenario(var_name = 'x.', 
                                     effect_types = c('effects.se.rf'), 
                                     result = result)
  

  alldat <- do.call('rbind', scenario_list)

  
  
  plot_result <- ggplot(alldat, aes(y = value)) +
    ggh4x::facet_grid2(get_label_formula(alldat), scales="free", axes="all")+
    geom_boxplot(aes(fill=variable)) + 
    labs(y = 'SE Estimates of Variable Effects', 
         x = 'Scenarios', fill = "Variable",
         title = 'Jackknife-after Bootstrap: Estimating Standard Errors of Variable Effects \n(for a given Data Set)',
         caption = TeX(paste0('Remaining Settings: ', get_caption(result, alldat), '$'))) +
    theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.line = element_line(),
          axis.title = element_text(size=18),
          panel.background = element_blank(),
          strip.text.x = element_text(size = 9),
          strip.text.y = element_text(size = 12),
          legend.title = element_text(size=15), #change legend title font size
          legend.text = element_text(size=15),
          plot.title = element_text(hjust = 0.5, lineheight = 0.9, size = 26),
          plot.caption = element_text(hjust = 0, size = 15))
  
  return(plot_result)
}


plot_se_dense <- function(result){
  
  scenario_list <- data_per_scenario(var_name = 'x.', 
                                     effect_types = c('effects.se.rf'), 
                                     result = result)
  
  alldat <- do.call('rbind', scenario_list)
  
  
  plot_result <- ggplot(alldat, aes(x = value)) +
    ggh4x::facet_grid2(get_label_formula(alldat), scales="free", axes="all", independent = "all")+
    geom_histogram(aes(y = after_stat(density), fill = variable), binwidth = 0.05) +
    geom_density() + 
    labs(x = 'SE Estimates of Variable Effects', fill = "Variable",
         title = 'Jackknife-after Bootstrap: Estimating Standard Errors of Variable Effects',
         caption = TeX(paste0('Remaining Settings: ', get_caption(result, alldat), '$'))) +
    theme(axis.text.x = element_text(angle = 45),
          axis.title = element_text(size=18),
          strip.text.x = element_text(size = 9),
          strip.text.y = element_text(size = 12),
          legend.title = element_text(size=15), #change legend title font size
          legend.text = element_text(size=15),
          plot.title = element_text(hjust = 0.5, lineheight = 0.9, size = 26),
          plot.caption = element_text(hjust = 0, size = 15)) 

  
  return(plot_result)
}




### Data per Scenario ###
data_per_scenario <- function(var_name = 'x.', effect_types, result){
  
  n_scenarios <- length(result)
  scenario_list <- vector(mode = 'list', length = n_scenarios*length(effect_types))
  cnt <- 1 ; max_vars <- 0
  for (scenario in 1:n_scenarios) {
    n_vars <- length(unique(unlist(str_extract_all(result[[scenario]]$formula,"x.[0-9]+"))))
    for (effect_type in effect_types) {
      all_terms_extracted <- extract_terms(result[[scenario]]$formula)
      if (length(all_terms_extracted) > max_vars) {
        var_names <- all_terms_extracted
        max_vars <- length(var_names) 
      }
      dat_tmp <- data.frame(result[[scenario]][[effect_type]])
      dat <- cbind(dat_tmp, 
                   N = paste('N=', result[[scenario]]$N),
                   cor = paste('Cor=', result[[scenario]]$cor),
                   k = paste('k=', result[[scenario]]$k),
                   num.trees = paste('Trees=', result[[scenario]]$num.trees),
                   node_size = paste('Node Size=', result[[scenario]]$node_size),
                   n_vars = paste('#Variables=', n_vars),
                   effect.type = names(result[[scenario]][effect_type]))
      scenario_list[[cnt]] <-  dat %>%
        pivot_longer(cols = starts_with('x.'), names_to = 'variable')
      cnt <- cnt + 1
    }
  }
  out <- scenario_list 
  return(out)
}




### Get true Variable Effects ###
get_true_effects <- function(result){
  scenario_list <- data_per_scenario(var_name = 'x.', 
                                     effect_types = c('effects.f.true'), 
                                     result = result)
  alldat <- do.call('rbind', scenario_list)
  mean_df <- alldat %>% 
    group_by(N, cor, k, num.trees, node_size, variable, n_vars) %>% 
    summarise(mean_val = mean(value))
  return(mean_df)
}




### Extract all terms of Formula ###
extract_terms <- function(formula){
  
  all_terms_extracted <- c()
  cnt <- 1
  all_terms <- unlist(strsplit(formula, "(?<=.)(?=[+])|(?<=.)(?=[-])",perl = TRUE))
  for (term in all_terms) {
    if (length(unique(unlist(str_extract_all(term,"x.[0-9]+")))) == 1) {
      all_terms_extracted[cnt] <- unique(unlist(str_extract_all(term,"x.[0-9]+")))
    } else if (length(unique(unlist(str_extract_all(term,"x.[0-9]+")))) == 0) {
      cnt <- cnt - 1
    } else {
      all_terms_extracted[cnt] <- paste(unlist(str_extract_all(term,"x.[0-9]+")), collapse = "")
    }
    cnt <- cnt + 1
  }
  return(all_terms_extracted)
}


### Label Formula in Plot
get_label_formula <- function(alldat){
  tmp <- alldat[,!names(alldat) %in% c("variable", "effect.type", "value")]
  col_idx <- sapply(colnames(tmp), function(x){length(unique(tmp[[x]]))!=1})
  formula_labels <- as.formula(paste0('variable~',paste0( colnames(tmp)[col_idx], collapse = '+')))
  #formula_labels <- as.formula(paste0(paste0(colnames(tmp)[col_idx], collapse = '+'), '~variable'))
  return(formula_labels)
}



### Create Caption in Plot
get_caption <- function(result, alldat){
  tmp <- alldat[,!names(alldat) %in% c("variable", "effect.type", "value")]
  col_idx <- sapply(colnames(tmp), function(x){length(unique(tmp[[x]]))!=1})
  longest_latex_formula <- result[[1]]$longest_latex_formula
  return(paste0(paste(tmp[1,][!col_idx], collapse = ";  "), ';  Formula= $', longest_latex_formula))
}


################################################################################
################################################################################
################################################################################
################################################################################
################################################################################








# Simulation Function
# Given a formula fit random forest N times and check effect and se_effect estimates
sim <- function(scenario){
  
  k <- scenario[['k']] 
  N <- scenario[['N']]
  num.trees <- scenario[['N_Trees']]
  cor <- scenario[['Correlation']] 
  formula <- scenario[['Formula']]
  node_size <- scenario[['Node_Size']]
  reps <- scenario[['Reps']]
  
  
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
  all_terms_t <- unlist(strsplit(formula, "(?<=.)(?=[+])|(?<=.)(?=[-])",perl = TRUE))
  all_terms <- character()
  idx <- 1
  for (term in 1:length(all_terms_t)) {
    if ((str_count(all_terms_t[term],"\\(") >= 1) & (str_count(all_terms_t[term],"\\)") == 0)) {
      all_terms[idx] <- paste0(all_terms_t[term], all_terms_t[term+1])
      idx = idx + 1
    } else if ((str_count(all_terms_t[term],"\\(") == 0) & (str_count(all_terms_t[term],"\\)") >= 1)){
      
    } else {
      all_terms[idx] <- all_terms_t[term]
      idx <- idx + 1
    }
  }
  
  n_coefs <- length(all_terms)
  for (term in all_terms) {
    if (length(unique(unlist(str_extract_all(term,"x.[0-9]+")))) > 1) {
      interaction_term <- unique(unlist(str_extract_all(term,"x.[0-9]+")))
    }
  }
  
  
  formula_p <- parse(text = formula)
  
  y <- eval(formula_p, x) + rnorm(N, 0, 1)
  
  data <- data.frame(cbind(x, y))
  
  # Estimated Variable Means
  x_bar <- colMeans(x)
  # True Variable Means and Standard Deviations
  true_mean <- rep(x = 0, times = n_vars)
  names(true_mean) <- unique(unlist(str_extract_all(formula,"x.[0-9]+")))
  true_s <- rep(x = 1, times = n_vars)
  names(true_s) <- unique(unlist(str_extract_all(formula,"x.[0-9]+")))
  
  
  list_names <- c('effects.rf','effects.f.true', 'effects.se.rf')
  
  result_list <- vector(mode = 'list', length = length(list_names))
  result_list <- lapply(X = 1:length(list_names), FUN = function(x){
    init <- matrix(data = 0, nrow = reps, ncol = n_coefs)
    colnames(init) <- extract_terms(formula)
    result_list[[x]] <- init
  })
  
  names(result_list) <- list_names
  

  
  for (rep in 1:reps) {
    rf <- ranger( formula = y ~ .,
                  data = data,
                  num.trees = num.trees,
                  keep.inbag = T, 
                  oob.error = F, # Save computational time
                  min.node.size = node_size)
    
    
    var_effects <- numeric(n_coefs)
    var_effects_se <- numeric(n_coefs)
    var_effects_true <- numeric(n_coefs)
    

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
                                        k_ratio = k, 
                                        save = FALSE,
                                        predict.all = T,
                                        inbag.counts = rf$inbag.counts)
      # Using f_hat and sample mean&variance of X
      predictions.rf <- rowMeans(rf.predict$predictions)
      
      # Using f and true mean&variance of X
      predictions.f.true <- eval(formula_p, new_data.true)
      
      # Using f_hat predictions based on sample mean&variance of X 
      var_effects[v] <- (predictions.rf[2] - predictions.rf[1]) / (2*k*sd(x[,v]))
      
      # Using f predictions based on true mean&variance of X 
      var_effects_true[v] <- (predictions.f.true[2] - predictions.f.true[1]) / (2*k*true_s[v])
      
  
      effect.se.direct <- rf.predict$se_effect 
      effect.var <- pmax((rf.predict$cov[1,1] + rf.predict$cov[2,2]
                          - 2*rf.predict$cov[1,2]) /  (2*k*sd(x[,v]))^2, 0)
      
      var_effects_se[v] <- effect.se.direct #sqrt(effect.var)
      
      if (exists('interaction_term')) {
        
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
          
          
          rf.predict <- RangerForestPredict(rf$forest,
                                            new_data,
                                            type='se',
                                            k_ratio = k,
                                            save = FALSE, 
                                            se.method = 'jack_cov',
                                            predict.all = T,
                                            inbag.counts = rf$inbag.counts)
          
          
          ### Computing Predictions
          
          # Using f_hat and sample mean&variance of X
          predictions.rf <- rowMeans(rf.predict$predictions)
          # Using f and true mean&variance of X
          predictions.f.true <- eval(formula_p, new_data.true)
          
          var_effects[n_coefs] <- ((predictions.rf[5] - predictions.rf[3]) - (predictions.rf[6] - predictions.rf[4])) / (-4*k^2*sd(x[,interaction_term[1]])*sd(x[,interaction_term[2]]))
          var_effects_true[n_coefs] <- ((predictions.f.true[5] - predictions.f.true[3]) - (predictions.f.true[6] - predictions.f.true[4])) / (-4*k^2*true_s[interaction_term[1]]*true_s[interaction_term[2]])
          
          effect.se.direct <- rf.predict$se_effect
          effect.var <- pmax((sum(diag(rf.predict$cov))
                              - 2*rf.predict$cov[2,1]
                              - 2*rf.predict$cov[3,1] + 2*rf.predict$cov[3,2]
                              + 2*rf.predict$cov[4,1] - 2*rf.predict$cov[4,2] - 2*rf.predict$cov[4,3]) /  ((-4*k^2*sd(x[,interaction_term[1]])*sd(x[,interaction_term[2]])))^2, 0)
          
          var_effects_se[n_coefs] <- effect.se.direct #sqrt(effect.var)
        }
      }
      
    }
    
    result_list[['effects.rf']][rep,] <- var_effects
    result_list[['effects.f.true']][rep,] <- var_effects_true
    result_list[['effects.se.rf']][rep,] <- var_effects_se
    
  }
  
  result_list[['N']] <- scenario[['N']]
  result_list[['k']] <- scenario[['k']]
  result_list[['num.trees']] <- scenario[['N_Trees']]
  result_list[['cor']] <- scenario[['Correlation']] 
  result_list[['formula']] <- scenario[['Formula']]
  result_list[['node_size']] <- scenario[['Node_Size']]
  result_list[['reps']] <- scenario[['Reps']]
  result_list[['longest_latex_formula']] <- scenario[['Longest_Latex_formula']]
  
  
  return(result_list)
}



