library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)
library(ggh4x)


### Print Results for each Scenario
print_results <- function(result){
  
  for (i in 1:length(result)) {
    cat(paste0('Setting ',i, ': N ='), result[[i]]$N, '; k =', result[[i]]$k, 
        'N_Trees =', result[[i]]$num.trees,'; Correlation =', result[[i]]$cor,
        '; Minimum Node Size =', result[[i]]$node_size,
        sprintf(paste0(';\n%', nchar('Setting: ')*2, 's'), 'Formula ='), result[[i]]$formula, '; N_Trees =', result[[i]]$num.trees,
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
plot_effects <- function(result){
  
  scenario_list <- data_per_scenario(var_name = 'x.', 
                                     effect_types = c('effects.rf', 
                                                      'effects.lm'), 
                                     result = result)
            
  effects_df <- get_true_effects(result)
  alldat <- do.call('rbind', scenario_list) 
  
  alldat$effect.type <- as.factor(alldat$effect.type)
  
  levels(alldat$effect.type) <- c(effects.lm='LM', 
                                 effects.rf='RF')

  
  plot_result <- ggplot(alldat, aes(x = value, y=effect.type)) +
    ggh4x::facet_grid2(variable~N+cor+k+num.trees+node_size, scales="free_x", independent = "x")+
    geom_boxplot(aes(fill=effect.type)) + 
      geom_vline(data=effects_df, aes(xintercept=mean_val, color='TrueEffect')) + 
    labs(y = 'Model', x = 'Effect Estimates', fill = "Model") + 
    theme(axis.text.x = element_text(angle = 45)) +
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
    group_by(N, cor, k, num.trees, node_size, variable) %>% 
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
    ggh4x::facet_grid2(variable~N+cor+k+num.trees+node_size, scales="free", independent = "all")+
    geom_histogram(aes(y = after_stat(density), fill = variable), binwidth = 0.05) +
    geom_density() + 
    geom_vline(data=sd_df, aes(xintercept=sd_val, color='Sd_Variable_Effects')) +
    labs(x = 'SE Estimates of Variable Effects', fill = "Variable") + 
    theme(axis.text.x = element_text(angle = 45)) +
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
    group_by(N, cor, k, num.trees, node_size, variable) %>% 
    summarise(sd_val = sd(value))
  
  return(sd_df)
}




### Data per Scenario ###
data_per_scenario <- function(var_name = 'x.', effect_types, result){
  
  n_scenarios <- length(result)
  
  scenario_list <- vector(mode = 'list', length = n_scenarios*length(effect_types))
  cnt <- 1 ; max_vars <- 0
  for (scenario in 1:n_scenarios) {
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
                   num.trees = paste('N_Trees=', result[[scenario]]$num.trees),
                   node_size = paste('Node_Size=', result[[scenario]]$node_size),
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
    } else {
      all_terms_extracted[cnt] <- paste(unlist(str_extract_all(term,"x.[0-9]+")), collapse = "")
    }
    cnt <- cnt + 1
  }
  return(all_terms_extracted)
}



### Plot PDPs
plot_pdps <- function(result){
  
  
  for (res in 1:length(result)) {
    
    n_vars <- length(unique(unlist(str_extract_all(result$formula,"x.[0-9]+"))))
    col_plots <- floor(n_vars/2)
    row_plots <- n_vars - col_plots
    
    par(mfrow=row_plots,col_plots)
    plot(result$`1`$pdp$x1, 
         type = 'l', xlab = paste('Variable:', 'X1'))
    
    plot(result$`1`$pdp$x2, 
         type = 'l', xlab = paste('Variable:', 'X2'))
    
    plot(result$`1`$pdp$x3, 
         type = 'l', xlab = paste('Variable:', 'X3'))
    
    plot(result$`1`$pdp$x4, 
         type = 'l', xlab = paste('Variable:', 'X4'))
    
    
  }
  
}


### Plot PDPs
plot_pdps <- function(result){
  
  scenarios <- length(result)
  plot_list <- vector(mode = 'list', length = scenarios)
  for (s in 1:scenarios) {
    n_vars <- length(unique(unlist(str_extract_all(result[[s]]$formula,"x.[0-9]+"))))
    plot_list[[s]] <- vector(mode = 'list', length = n_vars)
    for (j in 1:n_vars) {
      plot_list[[s]][[j]] <- ggplot(result[[s]]$pdp[[j]], aes_string(x=paste0('x.', j), y='.value')) +
        geom_line(color = 'red') +
        geom_point() +
        labs(x = paste0('X', j), 
             y = paste0('PDP Values\nScenario ', s)) + 
        theme_bw()
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



### Plot PDPs - Option 2
scenario_pdps <- function(scenario, s){
  
  n_vars <- length(unique(unlist(str_extract_all(scenario$formula,"x.[0-9]+"))))
  col_plots <- floor(n_vars/2)
  row_plots <- n_vars - col_plots
  par(mfrow=c(row_plots,col_plots))
  for (v in 1:n_vars) {
    plot(scenario$pdp[[paste0('x', v)]], 
         type = 'l', xlab = paste('Variable:', paste0('X', v)))
  }
  mtext(paste("PDPs for Scenario", s), side = 3, line = -3, cex = 2, outer = TRUE)
}


plot_pdps2 <- function(result){
  
  scenarios <- length(result)
  for (s in 1:scenarios) {
    scenario_pdps(result[[s]], s)
  }
  
}



