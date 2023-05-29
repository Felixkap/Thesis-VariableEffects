library(tidyverse)
library(ggplot2)
library(ggpubr)
library(stringr)
library(caret)
library(pdp)
library(latex2exp)
library(MixMatrix)
library(mvtnorm)
library(Matrix)
source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/RangerPredictFunction.R')




### get_all_terms: Function to obtain all individual terms of a formula
get_all_terms <- function(f){
  
  ## @param f: Formula 
  
  if (!(is.character(f) & length(f)==1)) {
    stop('Argument f must be string of length 1')
  }
  
  all_terms <- unlist(strsplit(f, "(?<=.)(?=[+])|(?<=.)(?=[-])",perl = TRUE))
  
  return(all_terms)
}



### get_all_interactions: Function to obtain all interactions between variables
get_all_interactions <- function(f){
  
  ## @param f: Formula
  
  all_terms <- get_all_terms(f)
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




### generate_data:  Function to generate data set
generate_data <- function(N, n_vars, cor, formula, sigma_e, cat=0){
  
  ## @param N: Integer. Number of observations
  ## @param n_vars: Integer. Number of independent variables
  ## @param cor: Numeric Value. Correlation between independent variables
  ## @param formula: String. Functional Relationship between predictors and outcome
  ## @param sigma_e: Numeric Value. Noise in the data (standard deviation of errors)
  

  ### Warnings
  if (N%%1!=0 | n_vars%%1!=0) {
    stop('Argument N and argument n_vars must be integers!')
  }
  
  if (!is.numeric(cor)) {
    stop('Argument cor must be sepcified as numeric value with -1 < cor < +1.')
  } else if (!(cor>=-1 & cor<=1)) {
    stop('Argument cor must be sepcified as numeric value with -1 < cor < +1.')
  }
  
  if (!is.character(formula)) {
    stop('Argument formula must be a string.')
  }
  
  if (!is.numeric(sigma_e)) {
    stop('Argument sigma_e must be numeric value with sigma_e >= 0.')
  } else if (sigma_e < 0) {
    stop('Argument sigma_e must be numeric value with sigma_e >= 0.')
  }
  
  if (cat > n_vars) {
    stop(paste('Formula accomodates max.', n_vars, 'categorical variable(s).'))
  }
  
  ### Begin Function
  # Generate Data: Predictors
  if (n_vars == 1) {
    # Normal Distribution with Mean 0 and Variance 1
    if (cat) {
      x <- data.frame(x.1 = sample(c(0,1), N, replace = TRUE))
      x.fac <- data.frame(x.1 = lapply(x, as.factor))
    } else {
      x <- data.frame(x.1 = rnorm(N, 0, 1))
    }
  } else{
    # Multivariate Normal Distribution with Mean Vector 0 and diagonal Variance Covariance Matrix 
    if (cat==0) {
      x <- data.frame(x = rmvnorm(N, mean = rep(0, n_vars-cat), sigma = CSgenerate(n_vars-cat, cor)))
    } else {
      x1 <- data.frame(x = rmvnorm(N, mean = rep(0, n_vars-cat), sigma = CSgenerate(n_vars-cat, cor)))
      x2 <- data.frame(sapply((n_vars-cat+1):n_vars, FUN = function(v){sample(c(0,1), N, replace = TRUE)}))
      colnames(x2) <- paste0('x.',(n_vars-cat+1):n_vars)
      x2.fac <- data.frame(lapply(x2, as.factor))
      colnames(x2.fac) <- paste0('x.',(n_vars-cat+1):n_vars)
      x <- cbind(x1,x2)
      x.fac <- cbind(x1,x2.fac)
    }
  }
  

  # Estimated Means of all variables
  if (cat) {
    point_of_interest <- calculate_mean_mode(x.fac)
  } else {
    point_of_interest <- calculate_mean_mode(x)
  }
  
  
  # Generate Data: Outcome
  y <- eval(parse(text = formula), x) + rnorm(N, 0, sigma_e)
  
  # Generated Data
  if (cat) {
    data <- data.frame(cbind(x.fac, y))
    return(list(x =x.fac,
                point_of_interest = point_of_interest,
                y = y,
                data = data))
  } else {
    data <- data.frame(cbind(x, y))
    return(list(x =x,
                point_of_interest = point_of_interest,
                y = y,
                data = data))
  }
}





### plot_error_rate: Plot Loss against Number of trees
plot_error_rate <- function(fitted.rf, data, outcome='y', data_type='train',...){
  
  ## @param fitted.rf: Fitted Random Forest object of type ranger
  ## @param data: DataFrame. Contains independent and dependent variables
  ## @param data_type: String. Determines whether train, validation or test set is evaluated
  ## @param first_tree: How many trees of the fitted RF should be at least used to evaluate MSE?
  ## @param thinning: Increment of sequence first_tree:rf$num.trees

  ### Warnings
  if (!inherits(fitted.rf$forest, "ranger.forest")) {
    stop("Error: Invalid class of input object. Must provide object of type forest (RF fitted with ranger function)")
  } 
  
  if (!is.data.frame(data)) {
    stop('You must provide an object of type dataframe to argument data')
  }
  
  if (!(data_type %in% c('train', 'test', 'validation', 'oob'))) {
    stop('Argument data_type must be specified as "train", "test", "oob" or "validation"')
  }
  
  if (!is.character(outcome) | !(outcome %in% colnames(data))) {
    stop('Argument outcome must be of type character and must match one of the columns of the data set')
  }

  
  ### Handling Ellipses
  l <- list(...)
  # If thinning is not further specified, set it to 2 if data_type is not evaluated with oob
  # If thinning is not further specified, set it to 25 if data_type is evaluated with oob
  if (is.null(l$thinning) & data_type!='oob') {
    thinning = 2
  } else if (is.null(l$thinning) & data_type=='oob') {
    thinning = 25
  } else {
    thinning = l$thinning
  }
  
  
  # if first_tree is not further specified, set it to 30
  first_tree <- ifelse(is.null(l$first_tree), 1,l$first_tree)
  
  
  ### Begin Function
  if (data_type != 'oob') {
    move_avg <- function(vec, first_tree, thinning){
      sapply(seq(first_tree,length(vec), by=thinning), FUN = function(x){mean(vec[1:x])})
    }
    
    all_tree_predictions <- predict(fitted.rf, data, predict.all = T)$predictions
    moving_predictions <- t(apply(all_tree_predictions, 1, move_avg, thinning=thinning, first_tree=first_tree))
    rmse_vec <- sqrt(colMeans((moving_predictions - data[[outcome]])^2))
    
    return(plot(seq(first_tree,fitted.rf$num.trees, by=thinning), 
                rmse_vec, 
                type='l', 
                ylab = 'Root Mean Squared Error', 
                xlab = 'N_Trees', 
                col='red',
                main = paste('RMSE on', data_type, 'data vs. \nNumber of Trees')))
  } else {
    nt <- seq(first_tree, fitted.rf$num.trees, thinning)
    
    oob_mse <- vector("numeric", length(nt))
    
    
    for(i in 1:length(nt)){
      rf.tmp <- ranger( formula = fitted.rf$call$formula,
                        data = data,
                        num.trees = nt[i],
                        keep.inbag = T, 
                        oob.error = T, 
                        mtry = fitted.rf$mtry,
                        write.forest = FALSE,
                        importance = 'permutation',
                        min.node.size = fitted.rf$min.node.size)
      oob_mse[i] <- sqrt(rf.tmp$prediction.error)
    }
    
   
    return(plot(x = nt, 
                y = oob_mse, 
                col = "red", 
                type = "l",
                ylab = 'Root Mean Squared Error', 
                xlab = 'N_Trees',
                main = paste('OOB estimate for RMSE vs. \nNumber of Trees')))
  }
  
} 



### Function to obtain main and low-order interaction effects
summary.rf <- function(x, fitted.rf, k, 
                       point_of_interest,
                       interaction = 'None', 
                       interaction_terms = NULL,
                       moving_var_effects = FALSE,
                       moving_se_effects = FALSE,
                       ...){
  
  ## @param x: Dataframe. Contains independent variables
  ## @param fitted.rf: Fitted Random Forest object of type ranger
  ## @param k: Numeric Value. Range over which the variable effect is computed
  ## @param interaction: String. Must be specified as None, all or specific. 
  ##                     Determines whether only main effects or also low-order interaction effects (either all or only specific ones).
  ## @param interaction_terms: String or NULL.
  ##                           Only used if argument 'interaction' is set to 'specific'.
  ##                           Determines between which variables low-order interaction effects are estimated.
  ##                           Must be specified as variable_name:variable_name or NULL.
  ## @param moving_var_effects: Boolean. Determine if Moving VAR Variable Effect should be plotted
  ## @param moving_se_effects: Boolean. Determine if Moving SE Variable Effect should be plotted
  ## @param ...: Ellipses: Additional Arguments first_tree and thinning
                           # @param first_tree: How many trees of the fitted RF should be at least used to evaluate VAR or SE of effect?
                           # @param thinning: Increment of sequence first_tree:rf$num.trees


  ### Warnings
  if (!is.data.frame(x)) {
    stop('You must provide an object of type dataframe to argument x')
  }
  
  if (!inherits(fitted.rf$forest, "ranger.forest")) {
    stop("Error: Invalid class of input object. Must provide object of type forest (RF fitted with ranger function)")
  } 
  
  if (!(interaction %in% c('None', 'all', 'specific'))) {
    stop('Argument "interaction" must be specified as "None", "all" or "specific"')
  }
  
  if (interaction == 'None' |interaction == 'all') {
    if (!is.null(interaction_terms)) {
      stop('If argument interaction is specified as "None", you must set argument interaction_terms to NULL')
    }
  }
  
  if (is.character(point_of_interest)) {
    if (point_of_interest == 'mean_mode') {
      point_of_interest <- calculate_mean_mode(x)
    } else {
      stop('Argument "point_of_interest" must be either specified as "mean_mode" or as numerical vector')
    }
  } else if (length(point_of_interest) != ncol(x)) {
    stop('if argument "point_of_interest" is specified manually you must specifiy a value of interest for each variable')
  } else {
    names(point_of_interest) <- names(x)
  }
  
  
  # Check which variables are factor variables
  num.cat.idx <- sapply(x, is.factor)
  
  if (any(apply(x[,num.cat.idx, drop=F], 2, function(c){length(unique(c))!=2}))) {
    stop('For this workflow you can not provide a data set that contains factor or character variables with more than 2 levels.')
  }
  
  
  
  ### Handling Ellipses
  l <- list(...)
  # If thinning is not further specified, set it to 2
  thinning <- ifelse(is.null(l$thinning), 2,l$thinning)
  # if first_tree is not further specified, set it to 30
  first_tree <- ifelse(is.null(l$first_tree), 30,l$first_tree)


  ### Begin Function
  
  # Number of variables
  n_vars <- ncol(x)
  # Variable Names
  x_names <- names(x)
  # Convert character variables to factor variables
  x <- char_to_fac(x)

  
  
  
  # Matrix with effect estimates and corresponding SE
  effect.mat <- matrix(NA, nrow = n_vars, ncol = 2)
  rownames(effect.mat) <- x_names 
  colnames(effect.mat) <- c('Effect Estimates', 'SE')
  
  # List with data points used for computing variable effect and corresponding variance-covariance matrix
  cov_list <- vector(mode = 'list', length = n_vars)
  names(cov_list) <- x_names
  
  # List of plots for moving VAR Estimate of effects
  if (moving_var_effects) {
    plot_list_var <- vector(mode = 'list', length = n_vars)
    names(plot_list_var) <- x_names
  }
  
  # List of plots for moving SE Estimate of effects
  if (moving_se_effects) {
    plot_list_se <- vector(mode = 'list', length = n_vars)
    names(plot_list_se) <- x_names
  }

  #### Main Effects
  for (v in 1:n_vars) {
    
    # New data points needed for main effect
    if (num.cat.idx[v]) {
      x_a <- replace(point_of_interest, v, unique(sort(x[,v]))[1])
      x_b <- replace(point_of_interest, v, unique(sort(x[,v]))[2])
    } else {
      x_a <- replace(point_of_interest, v, point_of_interest[v]-k*sd(x[,v]))
      x_b <- replace(point_of_interest, v, point_of_interest[v]+k*sd(x[,v]))
    }
    
    newdata <- data.frame(rbind(x_a, x_b))
    
    # Predictions of Random Forest
    rf.predict.obj <- RangerForestPredict(object = fitted.rf$forest, 
                                          data = newdata,  
                                          type = 'se', 
                                          se.method = 'se_direct', 
                                          inbag.counts = fitted.rf$inbag.counts)
    
    
    # Estimate Main Effect
    rf.predictions <- rowMeans(rf.predict.obj$predictions)
    rf.effect <- (rf.predictions[2] - rf.predictions[1]) / sum(abs(as.numeric(newdata[2,])-as.numeric(newdata[1,])))
    effect.mat[v,1] <-ifelse(num.cat.idx[v], rf.effect/2, rf.effect)  
    
    # Estimate SE of Main Effect
    rf.se.effect <- rf.predict.obj$se_effect
    effect.mat[v,2] <- rf.se.effect
    

    # Get Variance Covariance Matrix of RF Predictions that were used to compute effect
    effect_points <- data.frame(rbind(x_a, x_b))
    
    cov_list[[v]] <- list(effect_points = effect_points,
                          cov = RangerForestPredict(object = fitted.rf$forest, 
                                                    data = newdata,
                                                    type = 'se', 
                                                    se.method = 'jack_cov', 
                                                    inbag.counts = fitted.rf$inbag.counts)$cov)
    
    
    if (moving_var_effects) {
      plot_list_var[[v]] <- plot_moving_var_effect(rf, v, newdata, k, 
                                                   first_tree, thinning)
    }
    
    if (moving_se_effects) {
      plot_list_se[[v]] <- plot_moving_se_effect(rf, v, newdata, k, 
                                                 first_tree, thinning)
    }
    
  }
  
 
  ### All possible interaction effects or specific interaction effects
  if (interaction == 'all') {
    
    if (!is.null(interaction_terms)) {
      stop('If argument interaction is specified as "all", you must set argument interaction_terms to NULL')
    }
    
    v_cnt <- n_vars + 1
    
    # Get all possible interaction terms
    all_interactions <- apply(t(combn(x_names, 2)), 1, function(r){
      paste0(r[1], ':', r[2])
    })
  
    interaction_term <- all_interactions
    # For each interaction term form 4 new points and estimate interaction effect
    for (interaction_term in all_interactions) {
      
      
      v_idx <-  which(x_names %in% str_split(interaction_term, ':')[[1]])

      
      # New data points needed for main effect
      if (any(num.cat.idx[v_idx])) {
        
        if (num.cat.idx[v_idx][1]) {
          x_a <- replace(point_of_interest, v_idx[1], unique(sort(x[,v_idx[1]]))[1])   
          x_b <- replace(point_of_interest, v_idx[1], unique(sort(x[,v_idx[1]]))[2])
        } else {
          x_a <- replace(point_of_interest, v_idx[1], point_of_interest[v_idx[1]] - k*sd(x[,v_idx[1]]))
          x_b <- replace(point_of_interest, v_idx[1], point_of_interest[v_idx[1]] + k*sd(x[,v_idx[1]]))
        }
        
        if (num.cat.idx[v_idx][2]) {
          x_a1 <- replace(x_a, v_idx[2], unique(sort(x[,v_idx[2]]))[1]) 
          x_b1 <- replace(x_b, v_idx[2], unique(sort(x[,v_idx[2]]))[1]) 
          x_a2 <- replace(x_a, v_idx[2], unique(sort(x[,v_idx[2]]))[2]) 
          x_b2 <- replace(x_b, v_idx[2], unique(sort(x[,v_idx[2]]))[2]) 
        } else {
          x_a1 <- replace(x_a, v_idx[2], point_of_interest[v_idx[2]] - k*sd(x[,v_idx[2]]))
          x_b1 <- replace(x_b, v_idx[2], point_of_interest[v_idx[2]] - k*sd(x[,v_idx[2]]))
          x_a2 <- replace(x_a, v_idx[2], point_of_interest[v_idx[2]] + k*sd(x[,v_idx[2]]))
          x_b2 <- replace(x_b, v_idx[2], point_of_interest[v_idx[2]] + k*sd(x[,v_idx[2]]))
        }

      } else {
        x_a <- replace(point_of_interest, v_idx[1], point_of_interest[v_idx[1]] - k*sd(x[,v_idx[1]]))
        x_b <- replace(point_of_interest, v_idx[1], point_of_interest[v_idx[1]] + k*sd(x[,v_idx[1]]))
        
        x_a1 <- replace(x_a, v_idx[2], point_of_interest[v_idx[2]] - k*sd(x[,v_idx[2]]))
        x_b1 <- replace(x_b, v_idx[2], point_of_interest[v_idx[2]] - k*sd(x[,v_idx[2]]))
        x_a2 <- replace(x_a, v_idx[2], point_of_interest[v_idx[2]] + k*sd(x[,v_idx[2]]))
        x_b2 <- replace(x_b, v_idx[2], point_of_interest[v_idx[2]] + k*sd(x[,v_idx[2]]))
      }
      
      
      newdata <- data.frame(rbind(x_a1, x_b1, 
                                  x_a2, x_b2))
      
      # Predictions of Random Forest
      rf.predict.obj <- RangerForestPredict(object = fitted.rf$forest, 
                                            data = newdata, 
                                            type = 'se', 
                                            se.method = 'se_direct', 
                                            inbag.counts = fitted.rf$inbag.counts)
      
      # Estimate low-order interaction Effect
      rf.predictions <- rowMeans(rf.predict.obj$predictions)
      rf.effect <- ((rf.predictions[3] - rf.predictions[1]) - (rf.predictions[4] - rf.predictions[2])) / prod(abs((as.numeric(newdata[3,]) - as.numeric(newdata[2,])))[which(abs(as.numeric(newdata[3,]) - as.numeric(newdata[2,])) != 0)])
    
      
      # Estimate SE of low-order interaction Effect
      rf.se.effect <- rf.predict.obj$se_effect

      effect.mat <- rbind(effect.mat, c(rf.effect, rf.se.effect))
      rownames(effect.mat)[v_cnt] <- interaction_term
      
      # Get Variance Covariance Matrix of RF Predictions that were used to compute effect
      effect_points = data.frame(rbind(x_a1, x_b1, x_a2, x_b2))
      
      
      cov_list[[interaction_term]] <- list(effect_points = effect_points,
                                           cov = RangerForestPredict(object = fitted.rf$forest, 
                                                                     data = newdata, 
                                                                     type = 'se', 
                                                                     se.method = 'jack_cov', 
                                                                     inbag.counts = fitted.rf$inbag.counts)$cov)
      
      if (moving_var_effects) {
        plot_list_var[[interaction_term]] <- plot_moving_var_effect(rf, v_cnt, newdata, k, 
                                                                    first_tree, thinning)
      }
      
      if (moving_se_effects) {
        plot_list_se[[interaction_term]] <- plot_moving_se_effect(rf, v_cnt, newdata, k, 
                                                                  first_tree, thinning)
      }
      
      v_cnt <- v_cnt + 1 
      
    }
    
  } else if (interaction == 'specific') {
    
    if (!is.character(interaction_terms)) {
        stop('If argument "interaction" is specified as "specific", argument "interaction_term" must be a string and specified as variable_name:variable_name')
    }
    
    
    for (interaction_term in interaction_terms) {
      if (sum(x_names %in% str_split(interaction_term, ':')[[1]])<2) {
        stop('Argument "interaction_term" must be specified as variable_name:variable_name')
      }
    }
    
    

    v_cnt <- n_vars + 1
    
    # For each interaction term form 4 new points and estimate interaction effect
    for (interaction_term in interaction_terms) {
      
      
      v_idx <-  which(x_names %in% str_split(interaction_term, ':')[[1]])
      
      
      # New data points needed for main effect
      if (any(num.cat.idx[v_idx])) {
        
        if (num.cat.idx[v_idx][1]) {
          x_a <- replace(point_of_interest, v_idx[1], unique(sort(x[,v_idx[1]]))[1])   
          x_b <- replace(point_of_interest, v_idx[1], unique(sort(x[,v_idx[1]]))[2])
        } else {
          x_a <- replace(point_of_interest, v_idx[1], point_of_interest[v_idx[1]] - k*sd(x[,v_idx[1]]))
          x_b <- replace(point_of_interest, v_idx[1], point_of_interest[v_idx[1]] + k*sd(x[,v_idx[1]]))
        }
        
        if (num.cat.idx[v_idx][2]) {
          x_a1 <- replace(x_a, v_idx[2], unique(sort(x[,v_idx[2]]))[1]) 
          x_b1 <- replace(x_b, v_idx[2], unique(sort(x[,v_idx[2]]))[1]) 
          x_a2 <- replace(x_a, v_idx[2], unique(sort(x[,v_idx[2]]))[2]) 
          x_b2 <- replace(x_b, v_idx[2], unique(sort(x[,v_idx[2]]))[2]) 
        } else {
          x_a1 <- replace(x_a, v_idx[2], point_of_interest[v_idx[2]] - k*sd(x[,v_idx[2]]))
          x_b1 <- replace(x_b, v_idx[2], point_of_interest[v_idx[2]] - k*sd(x[,v_idx[2]]))
          x_a2 <- replace(x_a, v_idx[2], point_of_interest[v_idx[2]] + k*sd(x[,v_idx[2]]))
          x_b2 <- replace(x_b, v_idx[2], point_of_interest[v_idx[2]] + k*sd(x[,v_idx[2]]))
        }
        
      } else {
        x_a <- replace(point_of_interest, v_idx[1], point_of_interest[v_idx[1]] - k*sd(x[,v_idx[1]]))
        x_b <- replace(point_of_interest, v_idx[1], point_of_interest[v_idx[1]] + k*sd(x[,v_idx[1]]))
        
        x_a1 <- replace(x_a, v_idx[2], point_of_interest[v_idx[2]] - k*sd(x[,v_idx[2]]))
        x_b1 <- replace(x_b, v_idx[2], point_of_interest[v_idx[2]] - k*sd(x[,v_idx[2]]))
        x_a2 <- replace(x_a, v_idx[2], point_of_interest[v_idx[2]] + k*sd(x[,v_idx[2]]))
        x_b2 <- replace(x_b, v_idx[2], point_of_interest[v_idx[2]] + k*sd(x[,v_idx[2]]))
      }
      
      
      newdata <- data.frame(rbind(x_a1, x_b1, 
                                  x_a2, x_b2))
      
      # Predictions of Random Forest
      rf.predict.obj <- RangerForestPredict(object = fitted.rf$forest, 
                                            data = newdata, 
                                            type = 'se', 
                                            se.method = 'se_direct', 
                                            inbag.counts = fitted.rf$inbag.counts)
      
      # Estimate low-order interaction Effect
      rf.predictions <- rowMeans(rf.predict.obj$predictions)
      rf.effect <- ((rf.predictions[3] - rf.predictions[1]) - (rf.predictions[4] - rf.predictions[2])) / prod(abs((as.numeric(newdata[3,]) - as.numeric(newdata[2,])))[which(abs(as.numeric(newdata[3,]) - as.numeric(newdata[2,])) != 0)])
      
      # Estimate SE of low-order interaction Effect
      rf.se.effect <- rf.predict.obj$se_effect
      
      effect.mat <- rbind(effect.mat, c(rf.effect, rf.se.effect))
      rownames(effect.mat)[v_cnt] <- interaction_term
      
      # Get Variance Covariance Matrix of RF Predictions that were used to compute effect
      effect_points = data.frame(rbind(x_a1, x_b1, x_a2, x_b2))
      
      cov_list[[interaction_term]] <- list(effect_points = effect_points,
                                           cov = RangerForestPredict(object = fitted.rf$forest, 
                                                                     data = newdata, 
                                                                     type = 'se', 
                                                                     se.method = 'jack_cov', 
                                                                     inbag.counts = fitted.rf$inbag.counts)$cov)
      
      if (moving_var_effects) {
        plot_list_var[[interaction_term]] <- plot_moving_var_effect(rf, v_cnt, newdata, k, 
                                                                    first_tree, thinning)
      }
      
      if (moving_se_effects) {
        plot_list_se[[interaction_term]] <- plot_moving_se_effect(rf, v_cnt, newdata, k, 
                                                                  first_tree, thinning)
      }
      
      v_cnt <- v_cnt + 1 
    }
  } 
  
  if (moving_var_effects & moving_se_effects) {
    res <- list(summary = effect.mat, 
                cov_list = cov_list,
                moving_var_effect = plot_list_var,
                moving_se_effect = plot_list_se)
  } else if (moving_var_effects) {
    res <- list(summary = effect.mat, 
                cov_list = cov_list,
                moving_var_effect = plot_list_var)
  } else if (moving_se_effects) {
    res <- list(summary = effect.mat,
                cov_list = cov_list, 
                moving_se_effect = plot_list_se)
  } else {
    res <- list(summary = effect.mat,
                cov_list = cov_list,)
  }
  
  return(res)
}




plot_moving_var_effect <- function(fitted.rf, v, newdata, k, 
                                  first_tree, thinning){
  
  ## @param fitted.rf: Fitted Random Forest of type ranger
  ## @param v: Index of variable for which moving var variable effect is to be computed
  ## @param newdata: New data point at which the main or low-order interaction effects are evaluated
  ## @param k: Range over which the variable effect is computed
  ## @param first_tree: How many trees of the fitted RF should be at least used to evaluate VAR of effect?
  ## @param thinning: Increment of sequence first_tree:rf$num.trees
  
  
  #### Warnings
  if (!inherits(fitted.rf$forest, "ranger.forest")) {
    print("Error: Invalid class of input object. Must provide RF fitted with ranger function")
  } 
  
  
  if (!is.numeric(k)) {
    stop('Argumentk must be of type numeric')
  }
  
  if (!is.numeric(thinning)) {
    stop('Argument thinning must be of type numeric')
  } else if (thinning %% 1 != 0 | thinning < 1) {
    stop('Argument thinning must be integer >= 1')
  }


  
  #### Estimate moving SE of Variable Effect 

  moving_var_effect <- sapply(X = seq(first_tree,fitted.rf$num.trees, by = thinning), FUN = function(b){
    
    RangerForestPredict(fitted.rf$forest, 
                        newdata, 
                        num.trees = b, 
                        type = 'se', 
                        se.method = 'se_direct',
                        inbag.counts = fitted.rf$inbag.counts)$var_effect
  })
  
  plot_data <- data.frame(n_trees = seq(first_tree,fitted.rf$num.trees, by = thinning),
                          moving_var = moving_var_effect,
                          variable = v)
  
  gg <- ggplot(plot_data, aes(x = n_trees, y = moving_var)) +
    geom_line()+
    labs(x= 'Number of Trees', 
         y=TeX(paste0('VAR Estimate of $\\beta_{', v, '}')))+
    theme_bw()
  
  return(gg)
}







plot_moving_se_effect <- function(fitted.rf, v, newdata, k, 
                                   first_tree, thinning){
  
  ## @param fitted.rf: Fitted Random Forest of type ranger
  ## @param v: Index of variable for which moving se variable effect is to be computed
  ## @param newdata: New data point at which the main or low-order interaction effects are evaluated
  ## @param k: Range over which the variable effect is computed
  ## @param first_tree: How many trees of the fitted RF should be at least used to evaluate SE of effect?
  ## @param thinning: Increment of sequence first_tree:rf$num.trees
  
  
  #### Warnings
  if (!inherits(fitted.rf$forest, "ranger.forest")) {
    print("Error: Invalid class of input object. Must provide RF fitted with ranger function")
  } 
  
  
  if (!is.numeric(k)) {
    stop('Argumentk must be of type numeric')
  }
  
  if (!is.numeric(thinning)) {
    stop('Argument thinning must be of type numeric')
  } else if (thinning %% 1 != 0 | thinning < 1) {
    stop('Argument thinning must be integer >= 1')
  }
  
  
  
  #### Estimate moving SE of Variable Effect 
  
  moving_se_effect <- sapply(X = seq(first_tree,fitted.rf$num.trees, by = thinning), FUN = function(b){
    
    RangerForestPredict(fitted.rf$forest, 
                        newdata,  
                        num.trees = b, 
                        type = 'se', 
                        se.method = 'se_direct',
                        inbag.counts = fitted.rf$inbag.counts)$se_effect
  })
  
  plot_data <- data.frame(n_trees = seq(first_tree,fitted.rf$num.trees, by = thinning),
                          moving_se = moving_se_effect,
                          variable = v)
  
  gg <- ggplot(plot_data, aes(x = n_trees, y = moving_se)) +
    geom_line()+
    labs(x= 'Number of Trees', 
         y=TeX(paste0('SE Estimate of $\\beta_{', v, '}')))+
    theme_bw()
  
  return(gg)
}




compute_pdp <- function(data, model, var_name, n_grid = 100){
  
  # Extract variable values
  var_values <- data %>% 
    select(all_of(var_name)) %>% 
    distinct() %>% 
    pull()
  
  if (!is.numeric(var_values)) var_values <- as.numeric(var_values)
  
  # Generate grid of values
  grid_values <- seq(min(var_values), max(var_values), length.out = n_grid)
  
  # Create dataframe for storing results
  pdp_df <- data.frame(var_name = rep(var_name, n_grid),
                       grid_value = grid_values,
                       pdp = NA)
  
  # Loop through grid values and compute PDP
  for (i in 1:n_grid){
    data_pdp <- data
    data_pdp[[var_name]] <- grid_values[i]
    pred <- predict(model, data = data_pdp)$predictions
    pdp_df$pdp[i] <- mean(pred)
  }
  
  return(pdp_df)
}


plot_pdp <- function(fitted.rf, x, vars, n_grid=100){
  
  ## @param fitted.rf: Fitted Random Forest of type ranger
  ## @param vars: Variables for which the PDP should be plotted
  
  #### Warnings
  if (!inherits(fitted.rf$forest, "ranger.forest")) {
    print("Error: Invalid class of input object. Must provide RF fitted with ranger function")
  } 
  
  plt_length <- length(vars)
  plt_list <- vector(mode = 'list', length = plt_length)
  names(plt_list) <- vars
  
  for (v in vars) {
   
    par.v <- compute_pdp(data = x, model = fitted.rf, var_name = v, n_grid = n_grid)
    
    
    plt_list[[v]] <- ggplot(par.v, aes(x = grid_value, y = pdp)) +
      geom_line() +
      labs(x= v, 
           y='Partial Dependence Values', 
           title = paste0('Partial Dependence Plot\n for ', v))+
      theme_bw()
  
  }
 
  return(plt_list)
}



plot_2pdp <- function(fitted.rf, vars){
  
  ## @param fitted.rf: Fitted Random Forest of type ranger
  ## @param vars: Two Variables for which the PDP should be plotted
  
  #### Warnings
  if (!inherits(fitted.rf$forest, "ranger.forest")) {
    print("Error: Invalid class of input object. Must provide RF fitted with ranger function")
  } 
  
  if (length(vars) != 2) {
    stop('Length of vector vars must be 2')
  }
  
  
  par.v1.v2 <- partial(fitted.rf, pred.var = vars, chull = TRUE)
  plt.v1.v2 <- autoplot(par.v1.v2, contour = TRUE, 
                                   legend.title = "Partial\ndependence") +
    labs(title = paste0('Partial Dependence Plot for ', vars[1], ' and ', vars[2]))+
    theme_bw()
  
  return(suppressWarnings(print(plt.v1.v2)))
}



plot_varimp <- function(fitted.rf){
  
  ## @param fitted.rf: Fitted Random Forest of type ranger
  
  #### Warnings
  if (!inherits(fitted.rf$forest, "ranger.forest")) {
    print("Error: Invalid class of input object. Must provide RF fitted with ranger function")
  }
  
  var_importance <- ifelse(rf$variable.importance > 0, rf$variable.importance, 0.001)
  
  imp_table <- data.frame(variable=names(var_importance), importance=sqrt(var_importance), row.names=NULL)
  
  varimp <- ggplot(imp_table, aes(x=reorder(variable,importance), y=importance,fill=importance))+ 
    geom_bar(stat="identity", position="dodge")+ coord_flip()+
    ylab("Increase in RMSE after Variable Permutation")+
    xlab("")+
    ggtitle("Permutation Variable Importance")+
    scale_fill_gradient(low="red", high="blue")
  
  return(varimp)
}




char_to_fac <- function(dat){
  
  ## @param dat: Dataframe
  
  idx <- sapply(dat,is.character)
  dat[,idx] <- lapply(dat[,idx, drop=F], function(c){as.factor(c)})
  
  return(dat)
}




calculate_mean_mode <- function(data) {
  
  ## @param data: Dataframe
  
  # Calculate mean of numerical variables
  mean_numerical <- data.frame(t(colMeans(data[, sapply(data, is.numeric)], na.rm = TRUE)))

  # Calculate mode of factor variables
  var_modes <- sapply(data[, sapply(data, is.factor), drop=F], function(x) {
    tab <- table(x)
    mode_val <- names(tab)[which.max(tab)]
    return(mode_val)
  })

  save_fac_levels <- lapply(data[, sapply(data, is.factor), drop=F], levels)

  mode_factor <- data.frame(t(var_modes))
  
  if (length(mode_factor)!=0) {
    for (fac in 1:ncol(mode_factor)) {
      mode_factor[,fac] <-  factor(mode_factor[,fac], levels = save_fac_levels[[fac]])
    }
  }


  # Combine mean and mode into a vector
  result <- cbind(mean_numerical, mode_factor)[colnames(data)]

  return(result)
}

