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
library(ALEPlot)


source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/RangerPredictFunction.R')
source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/Helper.R')

##### Simulation --- Distribution of Variable Effect and Standard Error
# Fit Random Forest & Linear Regression Model N times
# Compute Variable Effects (as we defined it) N times
# Estimate Standard Error of Variable Effects N times directly
# ---> Based on estimated Variance-Covariance Matrix

# Compare Standard Error of simulated Variable Effects with mean of directly estimated Standard Errors of Variable Effects






sim_once <- function(k, N, num.trees, cor, formula, node_size, pdp, ale, k_idx){
  
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
  
  
  rf <- ranger( formula = y ~ .,
                data = data,
                num.trees = num.trees,
                keep.inbag = T, 
                oob.error = F, # Save computational time
                min.node.size = node_size)
  
  
  
  
  
  if (pdp & k_idx) {
    predictor <- Predictor$new(rf, data=x, y=y)
  }

  
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
                  'effects.se.rf', 'effects.se.lm', 
                  'smaller_nulls', 'plt_list')
  
  result_list <- vector(mode = 'list', length = length(list_names))
  result_list <- lapply(X = 1:length(list_names), FUN = function(X){
    if (X==length(list_names)) {
      result_list[[X]] = vector(mode = 'list', length = n_vars)
    } else {
      result_list[[X]] = numeric(length = n_vars)
    }})
  names(result_list) <- list_names
  names(result_list[['plt_list']]) <- paste0('x', 1:n_vars)
  
  
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
    
    lm.predict <- predict(linreg,
                          newdata = new_data,
                          se.fit = T)
    
    ### Computing Predictions
    
    # Using f_hat and sample mean&variance of X
    predictions.rf <- rowMeans(rf.predict$predictions)
    
    predictions.lm <- lm.predict$fit
    
    # Using f and true mean&variance of X
    predictions.f.true <- eval(formula_p, new_data.true)
    
    ### Computing Variable Effects
    
    # Using f_hat predictions based on sample mean&variance of X 
    result_list[['effects.rf']][v] <- (predictions.rf[2] - predictions.rf[1]) / (2*k*sd(x[,v]))
    result_list[['effects.lm']][v] <- (predictions.lm[2] - predictions.lm[1]) / (2*k*sd(x[,v]))
    
    
    # Using f predictions based on true mean&variance of X 
    result_list[['effects.f.true']][v] <- (predictions.f.true[2] - predictions.f.true[1]) / (2*k*true_s[v])
    
    
    # result_list[['smaller_nulls']][v] <- sum((rf.predict$cov[1,1] + rf.predict$cov[2,2]
    #                                           - 2*rf.predict$cov[1,2]) /  (2*k*sd(x[,v]))^2 < 0)
    # 
    # Standard Error of Variable Effect (Coviariance Estimation)
    # VAR[(f(B) - f(A)) / 2*sd(x)]
    # =( 1 / 4*sd(x)^2 ) * ( VAR[f(B)] + VAR[f(B)] - 2*COV[f(A), f(B)] )
    
    effect.se.direct <- rf.predict$se_effect
    # effect.var <- pmax((rf.predict$cov[1,1] + rf.predict$cov[2,2]
    #                     - 2*rf.predict$cov[1,2]) /  (2*k*sd(x[,v]))^2, 0)
    # 
    result_list[['effects.se.rf']][v] <- effect.se.direct #sqrt(effect.var)
    result_list[['effects.se.lm']][v] <- summary(linreg)$coef[,"Std. Error"][v+1]
    
    
    ############################################################################
    ############################################################################
    ############################################################################
    ######################        PDP PLOT       ###############################
    ############################################################################
    ############################################################################
    ############################################################################
    yhat <- function(X.model, newdata){
      
      rf.predict <- RangerForestPredict(X.model$forest,
                                        newdata,
                                        type='se',
                                        se.method = 'jack_cov',
                                        k_ratio = k, 
                                        save = FALSE,
                                        predict.all = T,
                                        inbag.counts = rf$inbag.counts)
      predictions.rf <- rowMeans(rf.predict$predictions)
      return(predictions.rf)
    } 
    
    if (pdp & k_idx) {
      pdp_v <- FeatureEffect$new(predictor = predictor, 
                                 feature = paste0('x.', v),
                                 method = "pdp", 
                                 grid.points = seq(-3,3, length.out=20))
      
      result_list[['plt_list']][[paste0('x', v)]] <- pdp_v$results[,-ncol(pdp_v$results)]
    } else if (ale & k_idx) {
      ale_vals <- ALEPlot(x, rf, pred.fun=yhat, J=v, K=20)
      result_list[['plt_list']][[paste0('x', v)]] <- data.frame(ale_vals$x.values, ale_vals$f.values)
    }
    

    
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
        
        lm.predict <- predict(linreg,
                              newdata = new_data,
                              se.fit = T)
        
        ### Computing Predictions
        
        # Using f_hat and sample mean&variance of X
        predictions.rf <- rowMeans(rf.predict$predictions)
        predictions.lm <- lm.predict$fit
        
        # Using f and true mean&variance of X
        predictions.f.true <- eval(formula_p, new_data.true)
        
        
        result_list[['effects.rf']][n_coefs] <- ((predictions.rf[5] - predictions.rf[3]) - (predictions.rf[6] - predictions.rf[4])) / (-4*k^2*sd(x[,interaction_term[1]])*sd(x[,interaction_term[2]]))
        result_list[['effects.lm']][n_coefs] <- ((predictions.lm[5] - predictions.lm[3]) - (predictions.lm[6] - predictions.lm[4])) / (-4*k^2*sd(x[,interaction_term[1]])*sd(x[,interaction_term[2]]))
        
        # Using f predictions based on true mean&variance of X 
        
        result_list[['effects.f.true']][n_coefs] <- ((predictions.f.true[5] - predictions.f.true[3]) - (predictions.f.true[6] - predictions.f.true[4])) / (-4*k^2*true_s[interaction_term[1]]*true_s[interaction_term[2]])
        # 
        # result_list[['smaller_nulls']][n_coefs] <- sum((sum(diag(rf.predict$cov))
        #                                                 - 2*rf.predict$cov[2,1]
        #                                                 - 2*rf.predict$cov[3,1] + 2*rf.predict$cov[3,2]
        #                                                 + 2*rf.predict$cov[4,1] - 2*rf.predict$cov[4,2] - 2*rf.predict$cov[4,3]) /  ((-4*k^2*sd(x[,interaction_term[1]])*sd(x[,interaction_term[2]])))^2 < 0)
        
        effect.se.direct <- rf.predict$se_effect
        # effect.var <- pmax((sum(diag(rf.predict$cov))
        #                     - 2*rf.predict$cov[2,1]
        #                     - 2*rf.predict$cov[3,1] + 2*rf.predict$cov[3,2]
        #                     + 2*rf.predict$cov[4,1] - 2*rf.predict$cov[4,2] - 2*rf.predict$cov[4,3]) /  ((-4*k^2*sd(x[,interaction_term[1]])*sd(x[,interaction_term[2]])))^2, 0)
        # 
        
        result_list[['effects.se.rf']][n_coefs] <- effect.se.direct #sqrt(effect.var)
        result_list[['effects.se.lm']][n_coefs] <- (summary(linreg)$coef[,"Std. Error"][v+1])
      }
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
  longest_latex_formula <- scenario[["Longest_Latex_formula"]]
  repeats <- scenario[["Repeats"]]
  cor <- scenario[["Correlation"]]
  node_size <- scenario[["Node_Size"]]
  pdp <- scenario[["pdp"]]
  ale <-  scenario[["ale"]]
  k_idx <- scenario[["k_idx"]]
  
  all_terms <- extract_terms(formula)
  
  res <- replicate(n = repeats,
                   expr = sim_once(k=k, N = N, num.trees = num.trees,
                                   cor = cor, formula = formula, 
                                   node_size = node_size, pdp = pdp,
                                   ale = ale, k_idx = k_idx))
  
  

  
  ### Distribution of Test Statistics
  res_scenario <- list(formula=formula, longest_latex_formula=longest_latex_formula, 
                       N=N, k=k, cor=cor, num.trees=num.trees, node_size=node_size,
                       effects.rf = do.call('rbind', res[1,]),
                       effects.lm = do.call('rbind', res[2,]),
                       #effects.f.sample = do.call('rbind', res[3,]),
                       #effects.rf.true = do.call('rbind', res[3,]),
                       effects.f.true = do.call('rbind', res[3,]),
                       effects.se.rf = do.call('rbind', res[4,]),
                       effects.se.lm = do.call('rbind', res[5,]),
                       n.nulls = do.call('rbind', res[6,]),
                       plt = if (pdp | ale) {
                         average_pdp(res[7,])
                       },
                       plt_list = if (pdp | ale) {
                         res[7,]
                       },
                       marginal_idx = if (pdp | ale) {
                         plt_idx(pdp, ale)
                       }
  )
  
  colnames(res_scenario$effects.rf) <- all_terms
  colnames(res_scenario$effects.lm) <- all_terms
  colnames(res_scenario$effects.f.true) <- all_terms
  colnames(res_scenario$effects.se.rf) <- all_terms
  colnames(res_scenario$effects.se.lm) <- all_terms
  
  return(res_scenario)
}




# # # ###### Simulation Setup
n <- c(40, 400, 4000) ; num.trees <- 2000 ; repeats <- 4; cor <- c(0, 0.8)
k <- c(0.2, 1); node_size <- c(1); pdp <- F; ale <- F
formulas <- c("2*x.1+4*x.2-3*x.3+2.2*x.4-x.3*x.4")
longest_latex_formula <- "2x_1+4x_2-3x_3+2.2x_4-x_3x_4"


#parallel::clusterExport(cl = clust, varlist = 'formulas')
scenarios <- data.frame(expand.grid(n, num.trees, formulas, repeats,
                                    cor, k, node_size, pdp, ale))
colnames(scenarios) = c("N", "N_Trees", "Formula", "Repeats",
                        "Correlation", "k", "Node_Size", "pdp", "ale")
scenarios$k_idx <- (scenarios$k == unique(scenarios$k)[1])
scenarios[,"Formula"] <- as.character(scenarios[,"Formula"]) ### Formula became Factor
scenarios["Longest_Latex_formula"] <- longest_latex_formula
scenarios <- split(scenarios, seq(nrow(scenarios)))
# #Run Simulation
system.time(result <- lapply(X = scenarios, FUN = sim_multi))

if (!pdp | !ale) {
  print_results(result)
}
effect_plots <- plot_effects(result)

se_plot <- plot_se(result)
effect_plots
se_plot



