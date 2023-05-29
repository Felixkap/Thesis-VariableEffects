library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)
library(ggnewscale)
library(ggh4x)
library(latex2exp)
source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/RangerPredictFunction.R')


### Print Results for each Scenario
print_results <- function(result, lm_results=T){
  
  if (lm_results) {
    for (i in 1:length(result)) {
      cat(paste0('Setting ',i, ': N ='), result[[i]]$N, '; k =', result[[i]]$k, 
          'N_Trees =', result[[i]]$num.trees,'; Correlation =', result[[i]]$cor,
          '; Minimum Node Size =', result[[i]]$node_size,
          sprintf(paste0(';%', nchar('Setting: ')*2, 's'), '\nFormula ='), result[[i]]$formula,
          '\nMean(s) of simulated RF Variable Effect(s):\n ',
          colMeans(result[[i]]$effects.rf),
          '\nMean(s) of simulated LM Variable Effect(s):\n ',
          colMeans(result[[i]]$effects.lm),
          '\nTrue Variable Effect(s):\n ',
          colMeans(result[[i]]$effects.f.true),
          '\nStandard Error of simulated Variable Effects (RF):\n ',
          apply(result[[i]]$effects.rf, 2, sd),
          '.\nMean of Standard Errors Estimates of Variable Effects (RF):\n ',
          colMeans((result[[i]]$effects.se.rf)),
          '.\nNumber of Smaller Nulls:\n ',
          apply(result[[i]]$n.nulls,2 , sum),'\n\n')
    }
  } else {
    for (i in 1:length(result)) {
      cat(paste0('Setting ',i, ': N ='), result[[i]]$N, '; k =', result[[i]]$k, 
          'N_Trees =', result[[i]]$num.trees,'; Correlation =', result[[i]]$cor,
          '; Minimum Node Size =', result[[i]]$node_size,
          sprintf(paste0(';%', nchar('Setting: ')*2, 's'), '\nFormula ='), result[[i]]$formula,
          '\nMean(s) of simulated RF Variable Effect(s):\n ',
          colMeans(result[[i]]$effects.rf),
          # '\nMean(s) of simulated RF Variable Effect(s) with true Means/Variances:\n ',
          # colMeans(result[[i]]$effects.rf.true),
          '\nTrue Variable Effect(s):\n ',
          colMeans(result[[i]]$effects.f.true),
          '\nStandard Error of simulated Variable Effects (RF):\n ',
          apply(result[[i]]$effects.rf, 2, sd),
          '.\nMean of Standard Errors Estimates of Variable Effects (RF):\n ',
          colMeans((result[[i]]$effects.se.rf)),
          '.\nNumber of Smaller Nulls:\n ',
          apply(result[[i]]$n.nulls,2 , sum),'\n\n')
    }
  }
  
}



### Plot Distribution of Effect Estimates
plot_effects <- function(result, ax='all', compare_lm=F){
  
  if (compare_lm) {
    scenario_list <- data_per_scenario(var_name = 'x.', 
                                       effect_types = c('effects.rf', 
                                                        'effects.lm'), 
                                       result = result)
  } else {
    scenario_list <- data_per_scenario(var_name = 'x.', 
                                       effect_types = c('effects.rf'), 
                                       result = result)
  }
  
  
  effects_df <- get_true_effects(result)
  alldat <- do.call('rbind', scenario_list) 
  
  alldat$effect.type <- as.factor(alldat$effect.type)
  
  if (sum(sapply(alldat[['variable']], function(x){str_count(x)!=3})) == 0) {
    title = 'Estimating Variable Main Effects'
  } else {
    title = 'Estimating Variable Main and Interaction Effects'
  }
  
  if (compare_lm) {
    levels(alldat$effect.type) <- c(effects.lm='LM', 
                                    effects.rf='RF')
  } else {
    levels(alldat$effect.type) <- c(effects.rf='RF')
  }
  
  
  plot_result <- ggplot(alldat, aes(x = effect.type, y=value)) +
    ggh4x::facet_grid2(get_label_formula(alldat), scales="free", axes=ax, independent=ax)+
    geom_boxplot(aes(fill=effect.type)) + 
    geom_hline(data=effects_df, aes(yintercept=mean_val, color='TrueEffect')) + 
    labs(y = 'Effect Estimates', x = 'Model', fill = "Model", 
         title = title,
         caption = TeX(paste0('Remaining Settings: ', get_caption(result, alldat), '$'))) + 
    theme(axis.text.x = element_text(angle = 45),
          axis.title = element_text(size=18),
          strip.text.x = element_text(size = 9),
          strip.text.y = element_text(size = 12),
          legend.title = element_text(size=15), #change legend title font size
          legend.text = element_text(size=15),
          plot.title = element_text(hjust = 0.5, lineheight = 0.9, size = 26),
          plot.caption = element_text(hjust = 0, size = 15)) +
    scale_color_manual(name = "", values = c(TrueEffect = "orange"))
  
  return(plot_result)
}




### Get true Variable Effects 
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







### Plot distribution of Jackknife-after-Bootstrap Standard Errors
plot_se <- function(result, ax='all'){
  
  scenario_list <- data_per_scenario(var_name = 'x.', 
                                     effect_types = c('effects.se.rf'), 
                                     result = result)
  
  sd_df <- get_sd_effects(result)
  mean_sd_df <- get_mean_sd_effects(result)
  alldat <- do.call('rbind', scenario_list)
  

  plot_result <- ggplot(alldat, aes(x = value)) +
    ggh4x::facet_grid2(get_label_formula(alldat), scales="free", independent=ax)+
    geom_histogram(aes(y = after_stat(density), fill = variable)) +
    geom_density() + 
    geom_vline(data=sd_df, aes(xintercept=sd_val, color='Sd_Variable_Effects')) +
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
          plot.caption = element_text(hjust = 0, size = 15)) +
    scale_color_manual(name = "", 
                       values = c(Sd_Variable_Effects = "orange"), 
                       labels = 'Standard Deviation\nof simulated\nVariable Effects') + 
    new_scale_colour() + 
    geom_vline(data=mean_sd_df, aes(xintercept=mean_sd_val, color='Mean_Sd_Variable_Effects')) +
    scale_color_manual(name = "", 
                       values = c(Mean_Sd_Variable_Effects = "green"), 
                       labels = 'Mean of Standard\nError Estimates\nof Variable Effects')
  
  return(plot_result)
}



### Get Standard Deviation of Variable Effects 
get_sd_effects <- function(result){
  
  scenario_list <- data_per_scenario(var_name = 'x.', 
                                     effect_types = c('effects.rf'), 
                                     result = result)
  
  alldat <- do.call('rbind', scenario_list)
  
  sd_df <- alldat %>% 
    group_by(N, cor, k, num.trees, node_size, variable, n_vars) %>% 
    summarise(sd_val = sd(value))
  
  return(sd_df)
}


### Get Mean of Jackknife-after-Bootstrap SE estimates
get_mean_sd_effects <- function(result){
  
  scenario_list <- data_per_scenario(var_name = 'x.', 
                                     effect_types = c('effects.se.rf'), 
                                     result = result)
  
  alldat <- do.call('rbind', scenario_list)
  
  mean_sd_df <- alldat %>% 
    group_by(N, cor, k, num.trees, node_size, variable, n_vars) %>% 
    summarise(mean_sd_val = mean(value))
  
  return(mean_sd_df)
}




### Data per Scenario 
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
                   sigma_e = paste('sigma_e=', result[[scenario]]$sigma_e),
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
  
  return(unique(all_terms_extracted))
}


### Obtain all interaction terms from a formula
get_all_interactions <- function(f){
  
  all_terms <- unlist(strsplit(f, "(?<=.)(?=[+])|(?<=.)(?=[-])",perl = TRUE))
  cnt <- 1
  interaction_list <- list()
  for (term in all_terms) {
    if (length(unique(unlist(str_extract_all(term,"x.[0-9]+")))) > 1) {
      interaction_list[[cnt]] <- unique(unlist(str_extract_all(term,"x.[0-9]+")))
      cnt = cnt + 1
    }
  }
  all_interactions <- sapply(interaction_list, function(interaction_term){
    paste(interaction_term, collapse = ":")
  })
  
  return(all_interactions)
}


### Function to generate data set
generate_data <- function(N, n_vars, cor, formula, sigma_e){
  
  # Generate Data: Predictors
  if (n_vars == 1) {
    # Normal Distribution with Mean 0 and Variance 1
    x <- data.frame(x.1 = rnorm(N, 0, 1))
  } else{
    # Multivariate Normal Distribution with Mean Vector 0 and diagonal Variance Covariance Matrix 
    x <- data.frame(x = rmvnorm(N, mean = rep(0, n_vars), sigma = CSgenerate(n_vars, cor)))
  }
  
  # Estimated Means of all variables
  x_bar <- colMeans(x)
  
  # Generate Data: Outcome
  y <- eval(parse(text = formula), x) + rnorm(N, 0, sigma_e)
  
  # Generated Data
  data <- data.frame(cbind(x, y))
  
  return(list(x =x,
              x_bar = x_bar,
              y = y,
              data = data))
  
}




### Get correct formula for facet grid function (needed for plots). 
#   It is of the form variable~(all varying parameters)
get_label_formula <- function(alldat){ 
  
  tmp <- alldat[,!names(alldat) %in% c("variable", "effect.type", "value")]
  col_idx <- sapply(colnames(tmp), function(x){length(unique(tmp[[x]]))!=1})
  formula_labels <- as.formula(paste0('variable~',paste0( colnames(tmp)[col_idx], collapse = '+')))
  
  return(formula_labels)
  
}



### Get caption for plots
get_caption <- function(result, alldat){
  
  tmp <- alldat[,!names(alldat) %in% c("variable", "effect.type", "value")]
  col_idx <- sapply(colnames(tmp), function(x){length(unique(tmp[[x]]))!=1})
  longest_latex_formula <- result[[1]]$longest_latex_formula
  
  return(paste0(paste(tmp[1,][!col_idx], collapse = ";  "), ';  Formula= $', longest_latex_formula))
  
}



