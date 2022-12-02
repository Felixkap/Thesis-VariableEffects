library(Rcpp)


rangerCpp <- function(treetype, input_x, input_y, variable_names, mtry, num_trees, verbose, seed, num_threads, write_forest, importance_mode_r, min_node_size, split_select_weights, use_split_select_weights, always_split_variable_names, use_always_split_variable_names, prediction_mode, loaded_forest, snp_data, sample_with_replacement, probability, unordered_variable_names, use_unordered_variable_names, save_memory, splitrule_r, case_weights, use_case_weights, class_weights, predict_all, keep_inbag, sample_fraction, alpha, minprop, holdout, prediction_type_r, num_random_splits, sparse_x, use_sparse_data, order_snps, oob_error, max_depth, inbag, use_inbag, regularization_factor, use_regularization_factor, regularization_usedepth) {
  .Call('_ranger_rangerCpp', treetype, input_x, input_y, variable_names, mtry, num_trees, verbose, seed, num_threads, write_forest, importance_mode_r, min_node_size, split_select_weights, use_split_select_weights, always_split_variable_names, use_always_split_variable_names, prediction_mode, loaded_forest, snp_data, sample_with_replacement, probability, unordered_variable_names, use_unordered_variable_names, save_memory, splitrule_r, case_weights, use_case_weights, class_weights, predict_all, keep_inbag, sample_fraction, alpha, minprop, holdout, prediction_type_r, num_random_splits, sparse_x, use_sparse_data, order_snps, oob_error, max_depth, inbag, use_inbag, regularization_factor, use_regularization_factor, regularization_usedepth)
}



RangerForestPredict <- function(object, data, predict.all = FALSE,
                                  num.trees = object$num.trees, 
                                  type = "response", se.method = "jack",
                                  seed = NULL, num.threads = NULL,
                                  verbose = TRUE, inbag.counts = NULL, ...) {
  
  
  ## GenABEL GWA data
  if (inherits(data, "gwaa.data")) {
    snp.names <- snp.names(data)
    snp.data <- data@gtdata@gtps@.Data
    data <- data@phdata[, -1, drop = FALSE]
    gwa.mode <- TRUE
  } else {
    snp.data <- as.matrix(0)
    gwa.mode <- FALSE
  }
  
  
  ## Check forest argument
  if (!inherits(object, "ranger.forest")) {
    stop("Error: Invalid class of input object.")
  } else {
    forest <- object
  }
  if (is.null(forest$num.trees) ||
      is.null(forest$child.nodeIDs) || is.null(forest$split.varIDs) ||
      is.null(forest$split.values) || is.null(forest$independent.variable.names) ||
      is.null(forest$treetype)) {
    stop("Error: Invalid forest object.")
  }
  if (forest$treetype == "Survival" && (is.null(forest$chf) || is.null(forest$unique.death.times))) {
    stop("Error: Invalid forest object.")
  }
  
  ## Check for old ranger version
  if (length(forest$child.nodeIDs) != forest$num.trees || length(forest$child.nodeIDs[[1]]) != 2) {
    stop("Error: Invalid forest object. Is the forest grown in ranger version <0.3.9? Try to predict with the same version the forest was grown.")
  }
  if (!is.null(forest$dependent.varID)) {
    warning("Forest grown in ranger version <0.11.5, converting ...")
    forest <- convert.pre.xy(forest)
  }
  
  ## Prediction type
  if (type == "response" || type == "se" || type == "cov") {
    prediction.type <- 1
  } else if (type == "terminalNodes") {
    prediction.type <- 2
  } else if (type == "quantiles") {
    stop("Error: Apply predict() to the ranger object instead of the $forest object to predict quantiles.")
  } else {
    stop("Error: Invalid value for 'type'. Use 'response', 'se', 'terminalNodes', or 'quantiles'.")
  }
  
  ## Type "se" only for certain tree types
  if (type == "se" && se.method == "jack" && forest$treetype != "Regression") {
    stop("Error: Jackknife standard error prediction currently only available for regression.")
  }
  if (type == "se" && se.method == "infjack") {
    if (forest$treetype == "Survival") {
      stop("Error: Infinitesimal jackknife standard error prediction not yet available for survival.")
    } else if (forest$treetype == "Classification") {
      stop("Error: Not a probability forest. Set probability=TRUE to use the infinitesimal jackknife standard error prediction for classification.")
    }
  }
  
  ## Type "se" requires keep.inbag=TRUE
  if (type == "se" && is.null(inbag.counts)) {
    stop("Error: No saved inbag counts in ranger object. Please set keep.inbag=TRUE when calling ranger.")
  }
  
  ## Set predict.all if type is "se"
  if (type == "se") {
    predict.all <- TRUE
  }
  
  if (type == "se" || type == "cov") {
    predict.all <- TRUE
  }
  
  x <- data
  
  if (sum(!(forest$independent.variable.names %in% colnames(x))) > 0) {
    stop("Error: One or more independent variables not found in data.")
  }
  
  ## Subset to same column as in training if necessary
  if (length(colnames(x)) != length(forest$independent.variable.names) || any(colnames(x) != forest$independent.variable.names)) {
    x <- x[, forest$independent.variable.names, drop = FALSE]
  }
  
  ## Recode characters
  if (!is.matrix(x) && !inherits(x, "Matrix")) {
    char.columns <- sapply(x, is.character)
    if (length(char.columns) > 0) {
      x[char.columns] <- lapply(x[char.columns], factor)
    }
  }
  
  ## Recode factors if forest grown 'order' mode
  if (!is.null(forest$covariate.levels) && !all(sapply(forest$covariate.levels, is.null))) {
    x <- mapply(function(xx, yy) {
      if(is.null(yy)) {
        xx
      } else {
        new.levels <- setdiff(levels(xx), yy)
        factor(xx, levels = c(yy, new.levels), exclude = NULL)
      }
    }, x, forest$covariate.levels, SIMPLIFY = !is.data.frame(x))
  }
  if (is.list(x) && !is.data.frame(x)) {
    x <- as.data.frame(x)
  }
  
  ## Convert to data matrix
  if (!is.matrix(x) & !inherits(x, "Matrix")) {
    x <- data.matrix(x)
  }
  
  ## Check missing values
  if (any(is.na(x))) {
    offending_columns <- colnames(x)[colSums(is.na(x)) > 0]
    stop("Missing data in columns: ",
         paste0(offending_columns, collapse = ", "), ".", call. = FALSE)
  }
  
  ## Num threads
  ## Default 0 -> detect from system in C++.
  if (is.null(num.threads)) {
    num.threads = 0
  } else if (!is.numeric(num.threads) || num.threads < 0) {
    stop("Error: Invalid value for num.threads")
  }
  
  ## Seed
  if (is.null(seed)) {
    seed <- runif(1 , 0, .Machine$integer.max)
  }
  
  if (forest$treetype == "Classification") {
    treetype <- 1
  } else if (forest$treetype == "Regression") {
    treetype <- 3
  } else if (forest$treetype == "Survival") {
    treetype <- 5
  } else if (forest$treetype == "Probability estimation") {
    treetype <- 9
  } else {
    stop("Error: Unknown tree type.")
  }
  
  ## Defaults for variables not needed
  mtry <- 0
  importance <- 0
  min.node.size <- 0
  split.select.weights <- list(c(0, 0))
  use.split.select.weights <- FALSE
  always.split.variables <- c("0", "0")
  use.always.split.variables <- FALSE
  prediction.mode <- TRUE
  write.forest <- FALSE
  replace <- TRUE
  probability <- FALSE
  unordered.factor.variables <- c("0", "0")
  use.unordered.factor.variables <- FALSE
  save.memory <- FALSE
  splitrule <- 1
  alpha <- 0
  minprop <- 0
  case.weights <- c(0, 0)
  use.case.weights <- FALSE
  class.weights <- c(0, 0)
  keep.inbag <- FALSE
  sample.fraction <- 1
  holdout <- FALSE
  num.random.splits <- 1
  order.snps <- FALSE
  oob.error <- FALSE
  max.depth <- 0
  inbag <- list(c(0,0))
  use.inbag <- FALSE
  y <- matrix(c(0, 0))
  regularization.factor <- c(0, 0)
  use.regularization.factor <- FALSE
  regularization.usedepth <- FALSE
  ## Use sparse matrix
  if (inherits(x, "dgCMatrix")) {
    sparse.x <- x
    x <- matrix(c(0, 0))
    use.sparse.data <- TRUE
  } else {
    sparse.x <- Matrix(matrix(c(0, 0)))
    use.sparse.data <- FALSE
    x <- data.matrix(x)
  }
  
  ## Call Ranger
  result <- rangerCpp(treetype, x, y, forest$independent.variable.names, mtry,
                      num.trees, verbose, seed, num.threads, write.forest, importance,
                      min.node.size, split.select.weights, use.split.select.weights,
                      always.split.variables, use.always.split.variables,
                      prediction.mode, forest, snp.data, replace, probability,
                      unordered.factor.variables, use.unordered.factor.variables, save.memory, splitrule,
                      case.weights, use.case.weights, class.weights, 
                      predict.all, keep.inbag, sample.fraction, alpha, minprop, holdout, 
                      prediction.type, num.random.splits, sparse.x, use.sparse.data,
                      order.snps, oob.error, max.depth, inbag, use.inbag, 
                      regularization.factor, use.regularization.factor, regularization.usedepth)
  

  
  if (length(result) == 0) {
    stop("User interrupt or internal error.")
  }
  
  ## Prepare results
  result$num.samples <- nrow(x)
  result$treetype <- forest$treetype
  
  if (predict.all) {
    if (forest$treetype %in% c("Classification", "Regression")) {
      if (is.list(result$predictions)) {
        result$predictions <- do.call(rbind, result$predictions)
      } else {
        result$predictions <- array(result$predictions, dim = c(1, length(result$predictions)))
      }
    } else {
      if (is.list(result$predictions) & length(result$predictions) >= 1 & is.numeric(result$predictions[[1]])) {
        # Fix for single test observation
        result$predictions <- list(result$predictions)
      }
      result$predictions <- aperm(array(unlist(result$predictions), 
                                        dim = rev(c(length(result$predictions), 
                                                    length(result$predictions[[1]]), 
                                                    length(result$predictions[[1]][[1]])))))
    }
  } else {
    if (is.list(result$predictions)) {
      result$predictions <- do.call(rbind, result$predictions)
    } 
  }
  
  if (type == "response") {
    if (forest$treetype == "Classification" && !is.null(forest$levels)) {
      if (!predict.all) {
        result$predictions <- integer.to.factor(result$predictions, forest$levels)
      }
    } else if (forest$treetype == "Regression") {
      ## Empty
    } else if (forest$treetype == "Survival") {
      result$unique.death.times <- forest$unique.death.times
      result$chf <- result$predictions
      result$predictions <- NULL
      result$survival <- exp(-result$chf)
    } else if (forest$treetype == "Probability estimation") {
      if (predict.all) {
        ## Set colnames and sort by levels
        if (!is.null(forest$levels)) {
          colnames(result$predictions) <- forest$levels[forest$class.values]
          result$predictions <- result$predictions[, forest$levels[sort(forest$class.values)], , drop = FALSE]
        }
      } else {
        if (is.vector(result$predictions)) {
          result$predictions <- matrix(result$predictions, nrow = 1)
        }
        
        ## Set colnames and sort by levels
        if (!is.null(forest$levels)) {
          colnames(result$predictions) <- forest$levels[forest$class.values]
          result$predictions <- result$predictions[, forest$levels[sort(forest$class.values)], drop = FALSE]
        }
      }
    }
  } else if (type == "terminalNodes") {
    if (is.vector(result$predictions)) {
      result$predictions <- matrix(result$predictions, nrow = 1)
    }
  }
  
  ## Compute Jackknife
  if (type == "se") {
    ## Aggregated predictions
    if (length(dim(result$predictions)) > 2) {
      yhat <- apply(result$predictions, c(1, 2), mean)
    } else {
      yhat <- rowMeans(result$predictions)
    }
    
    ## Get inbag counts, keep only observations that are OOB at least once
    inbag.counts <- simplify2array(inbag.counts) 
    if (is.vector(inbag.counts)) {
      inbag.counts <- t(as.matrix(inbag.counts))
    }
    inbag.counts <- inbag.counts[rowSums(inbag.counts == 0) > 0, , drop = FALSE] 
    n <- nrow(inbag.counts)
    oob <- inbag.counts == 0 # Is observation in bootstrap sample b yes or no?

    if (num.trees != object$num.trees) {
      oob <- oob[, 1:num.trees]
    }
    
    if (all(!oob)) {
      stop("Error: No OOB observations found, consider increasing num.trees or reducing sample.fraction.")
    }
    
    if (se.method == "jack" || se.method == "jack_cov" || se.method == "jack_effect") {
      ## Compute Jackknife
      oob.count <- rowSums(oob) # In how many bootstrap samples appears observation?
      jack.n <- sweep(tcrossprod(result$predictions, oob), 
                      2, oob.count, "/", check.margin = FALSE)
      if (is.vector(jack.n)) {
        jack.n <- t(as.matrix(jack.n))
      }
      if (any(oob.count == 0)) {
        n <- sum(oob.count > 0)
        jack.n <- jack.n[, oob.count > 0]
      } 
      
      jack <- (n - 1) / n * rowSums((jack.n - yhat)^2)
      bias <- (exp(1) - 1) * n / result$num.trees^2 * rowSums((result$predictions - yhat)^2)
      jab <- pmax(jack - bias, 0)
      result$se <- sqrt(jab)
      
      if (se.method == "jack_cov") {
        if (nrow(data)!=2) {
          stop('Error: Covariance can only be computed between two data points')
        } else {
          ## Compute Jackknife
          
          jack.n1 <- sweep(tcrossprod(result$predictions[1,], oob), # For Covariance --- Observation 1 
                           2, oob.count, "/", check.margin = FALSE) 
          jack.n2 <- sweep(tcrossprod(result$predictions[2,], oob), # For Covariance --- Observation 2
                           2, oob.count, "/", check.margin = FALSE)
          
          if (is.vector(jack.n1)) {
            jack.n1 <- t(as.matrix(jack.n1))
          } else if (is.vector(jack.n2)) {
            jack.n2 <- t(as.matrix(jack.n2))
          }
          
          
          #### Covariance
          jack_cov <- (n - 2) / n * (sum((jack.n1 - yhat[1]) * (jack.n2 - yhat[2])))
          result$cov <- matrix(c(jab[1], rep(jack_cov,2), jab[2]), 2, 2)
          
        }
      } else if (se.method == "jack_effect"){
        if (nrow(data)!=2) {
          stop('Error: Variable Effect must only be computed based on two data points')
        } else {
          ## Compute Jackknife
          
          jack.n1 <- sweep(tcrossprod(result$predictions[1,], oob), # For Covariance --- Observation 1 
                           2, oob.count, "/", check.margin = FALSE) 
          jack.n2 <- sweep(tcrossprod(result$predictions[2,], oob), # For Covariance --- Observation 2
                           2, oob.count, "/", check.margin = FALSE) 
          
          if (is.vector(jack.n1)) {
            jack.n1 <- t(as.matrix(jack.n1))
          } else if (is.vector(jack.n2)) {
            jack.n2 <- t(as.matrix(jack.n2))
          }
          
          #### Direct Estimate for Standard Error of Variable Effect
          jack_var <- (n - 1) / n * sum(((jack.n2 - jack.n1) - (yhat[2] - yhat[1]))^2) 
          bias_var <- (exp(1) - 1) * n / result$num.trees^2 * sum(((result$predictions[2,] - result$predictions[1,]) - (yhat[2] - yhat[1]))^2) 
          jab_var <- (pmax(jack_var - bias_var, 0)) / 4
          result$se_effect <- sqrt(jab_var)
      }
     }
    } else if (se.method == "infjack") {
      if (forest$treetype == "Regression") {
        infjack <- rInfJack(pred = result$predictions, inbag = inbag.counts, used.trees = 1:num.trees)
        result$se <- sqrt(infjack$var.hat)
        } else if (forest$treetype == "Probability estimation") {
          infjack <- apply(result$predictions, 2, function(x) {
            rInfJack(x, inbag.counts)$var.hat
            })
        result$se <- sqrt(infjack)
      } 
    } else {
      stop("Error: Unknown standard error method (se.method).")
    }
    

    if (forest$treetype == "Probability estimation") {
      ## Set colnames and sort by levels
      colnames(result$predictions) <- forest$levels[forest$class.values]
      result$predictions <- result$predictions[, forest$levels, drop = FALSE]
      
      if (!is.matrix(result$se)) {
        result$se <- matrix(result$se, ncol = length(forest$levels))
      }
      colnames(result$se) <- forest$levels[forest$class.values]
      result$se <- result$se[, forest$levels, drop = FALSE]
    }
   }   
  class(result) <- "ranger.prediction"
  return(result)
}