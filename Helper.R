library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)
library(ggh4x)
library(latex2exp)
source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/RangerPredictFunction.R')


### Print Results for each Scenario
print_results <- function(result){
  
  for (i in 1:length(result)) {
    cat(paste0('Setting ',i, ': N ='), result[[i]]$N, '; k =', result[[i]]$k, 
        'N_Trees =', result[[i]]$num.trees,'; Correlation =', result[[i]]$cor,
        '; Minimum Node Size =', result[[i]]$node_size,
        sprintf(paste0(';%', nchar('Setting: ')*2, 's'), '\nFormula ='), result[[i]]$formula,
        '\nMean(s) of simulated RF Variable Effect(s):\n ',
        colMeans(result[[i]]$effects.rf),
        '\nMean(s) of simulated LM Variable Effect(s):\n ',
        colMeans(result[[i]]$effects.lm),
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



### Plot Distribution of Effect Size Estimates
### Plot Distribution of Effect Size Estimates
plot_effects <- function(result){
  
  scenario_list <- data_per_scenario(var_name = 'x.', 
                                     effect_types = c('effects.rf', 
                                                      'effects.lm'), 
                                     result = result)
  
  effects_df <- get_true_effects(result)
  alldat <- do.call('rbind', scenario_list) 
  
  alldat$effect.type <- as.factor(alldat$effect.type)
  
  if (sum(sapply(alldat[['variable']], function(x){str_count(x)!=3})) == 0) {
    title = 'Estimating Variable Main Effects'
  } else {
    title = 'Estimating Variable Main and Interaction Effects'
  }
  
  levels(alldat$effect.type) <- c(effects.lm='LM', 
                                  effects.rf='RF')
  
  
  plot_result <- ggplot(alldat, aes(x = effect.type, y=value)) +
    ggh4x::facet_grid2(get_label_formula(alldat), scales="free", independent="x")+
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







### Plot distribution of Variable Effect Standard Errors
plot_se <- function(result){
  
  scenario_list <- data_per_scenario(var_name = 'x.', 
                                     effect_types = c('effects.se.rf'), 
                                     result = result)
  
  sd_df <- get_sd_effects(result)
  alldat <- do.call('rbind', scenario_list)
  

  plot_result <- ggplot(alldat, aes(x = value)) +
    ggh4x::facet_grid2(get_label_formula(alldat), scales="free", independent="y")+
    geom_histogram(aes(y = after_stat(density), fill = variable), binwidth = 0.05) +
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
                       labels = 'Standard Deviation\nof simulated\nVariable Effects')
  
  return(plot_result)
}



### Get Standard Deviation of Variable Effects ###
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





### Plot PDPs
plot_average_pdps <- function(result){
  
  mean_df <- get_true_effects(result)

  
  scenarios <- length(result)
  plot_list <- vector(mode = 'list', length = scenarios)
  case <- 1
  for (s in 1:scenarios) {
    n_vars <- length(unique(unlist(str_extract_all(result[[s]]$formula,"x.[0-9]+"))))
    plot_list[[s]] <- vector(mode = 'list', length = n_vars)
    for (j in 1:n_vars) {
      true_effect <- as.numeric(mean_df[case, 'mean_val'])
      plot_list[[s]][[j]] <- ggplot(result[[s]]$plt[[j]], aes(x=paste0('x.', j), y='.value')) +
        geom_point() +
        geom_line() + 
        geom_abline(slope = true_effect, colour = 'red') + 
        labs(x = paste0('X', j), 
             y = paste0('PDP Values\nScenario ', s)) + 
        theme_bw()
        case = case + 1
    }
  }
  
  mygg <- function(x){
    ggarrange(plotlist =  x, nrow = 1)
  }
  
  return(ggarrange(plotlist = lapply(X = plot_list, FUN = mygg), nrow = length(plot_list)))
}
  



plot_marginal <- function(result){
  
  mean_df <- get_true_effects(result)
  idx_plot <- result[[1]][['marginal_idx']]
  
  
  scenarios <- length(result)
  
  k_vec <- numeric(length(scenarios))
  for (s in 1:scenarios) {
    k_vec[s] <- result[[s]]$k
  }
  idx <- which(unique(k_vec)[1] == k_vec)
  
  plot_list <- vector(mode = 'list', length = length(idx))
  case <- 1
  for (s in idx) {
    
    formula <- result[[s]][['formula']]
    n_vars <- length(unique(unlist(str_extract_all(result[[s]]$formula,"x.[0-9]+"))))
    reps <- length(result[[s]][['plt_list']])
    var_list <- vector(mode = 'list', length = n_vars)
    var_list <- lapply(X = 1:n_vars, FUN = function(X){
      var_list[[X]] = vector(mode = 'list', length = reps)
    })
    
    for (j in 1:n_vars) {
      for (rep in 1:reps) {
        dat <- result[[s]][['plt_list']][[rep]][[paste0('x', j)]]
        names(dat) <- c('x', 'values')
        var_list[[j]][[rep]] <- cbind(dat, rep = rep, variable = paste0('x', j))
      }
    }
    
    for (j in 1:n_vars) {
      
      true_effect <- as.numeric(mean_df[case, 'mean_val'])
      dat <- do.call('rbind', var_list[[j]])
      plot_list[[s]][[j]] <- ggplot(dat) +
        geom_line(aes(x=x, y=values,group=rep), alpha = 0.3) +
        geom_smooth(formula = y ~ x, aes(x=x, y=values), method="loess", se=F) +
        geom_abline(slope = true_effect, colour = 'red') + 
        labs(x = paste0('X', j), 
             y = idx_plot) + 
        ggtitle(paste0('N=', result[[s]]$N, ', Cor=', result[[s]]$cor, 
                       ', \nMin_Node_size=', result[[s]]$node_size,
                       ', \nf(x)=', formula)) +
        theme(plot.title = element_text(hjust = 0.5))
        theme_bw()
      case = case + 1
    }
  }
  mygg <- function(x){
    ggarrange(plotlist =  x, nrow = 1)
  }
  
  return(ggarrange(plotlist = lapply(X = plot_list, FUN = mygg), nrow = length(plot_list)))
}







### Sum up elemets listwise
listwise_summation <- function(l1, l2){
  Map('+', l1, l2)
}


average_pdp <- function(pdp_list){
  summed_lists <- Reduce(listwise_summation, pdp_list)
  out <- lapply(summed_lists, function(x){x/length(pdp_list)})
  return(out)
}


plt_idx <- function(pdp, ale){
  if (pdp) {
    out <- 'PDP Values'
  } else if (ale) {
    out <- 'ALE Values'
  }
  return(out)
}


get_label_formula <- function(alldat){
  tmp <- alldat[,!names(alldat) %in% c("variable", "effect.type", "value")]
  col_idx <- sapply(colnames(tmp), function(x){length(unique(tmp[[x]]))!=1})
  formula_labels <- as.formula(paste0('variable~',paste0( colnames(tmp)[col_idx], collapse = '+')))
  #formula_labels <- as.formula(paste0(paste0(colnames(tmp)[col_idx], collapse = '+'), '~variable'))
  return(formula_labels)
}




get_caption <- function(result, alldat){
  tmp <- alldat[,!names(alldat) %in% c("variable", "effect.type", "value")]
  col_idx <- sapply(colnames(tmp), function(x){length(unique(tmp[[x]]))!=1})
  longest_latex_formula <- result[[1]]$longest_latex_formula
  return(paste0(paste(tmp[1,][!col_idx], collapse = ";  "), ';  Formula= $', longest_latex_formula))
}



