library(Rcpp)
library(SAEforest)

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
  if (type == "response" || type == "se") {
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
      yhat <- rowMeans(result$predictions, na.rm = T)
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
    
    if (se.method == "jack" || se.method == "jack_cov") {
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

        ## Compute Jackknife
        
        jack.n <- sweep(tcrossprod(result$predictions, oob), # For Covariance --- Observation 1 and 2
                         2, oob.count, "/", check.margin = FALSE) 
      
        
        if (is.vector(jack.n)) {
          jack.n <- t(as.matrix(jack.n))
        }
        
        
        ######################################################################
        ############### Estimate Variance - Covariance Matrix ################
        ######################################################################
        n_obs <- nrow(jack.n)
        obs_vec <- sapply(1:n_obs, function(i){jack.n[i,] - yhat[i]})
        
        
        
        result$cov <- matrix(NA, n_obs, n_obs)
        for (i in 1:n_obs) {
          for (j in 1:n_obs) {
            if (i==j) {
              result$cov[i,j] = jab[i]
            } else {
              jack_cov <- (n - 1) / n * sum(obs_vec[,i] * obs_vec[,j])
              bias <- (exp(1) - 1) * n / result$num.trees^2 * sum((result$predictions[i,] - yhat[i])*(result$predictions[j,] - yhat[j]))
              jack_cov <- jack_cov - bias
              result$cov[i,j] = jack_cov
            }
          }
        }
        
        
        
        ######################################################################
        ######## Direct Estimate for Standard Error of Variable Effect #######
        ######################################################################
        jack_var <- (n - 1) / n * sum(((jack.n[2,] - jack.n[1,]) - (yhat[2] - yhat[1]))^2) 
        bias_var <- (exp(1) - 1) * n / result$num.trees^2 * sum(((result$predictions[2,] - result$predictions[1,]) - (yhat[2] - yhat[1]))^2) 
        jab_var <- (pmax(jack_var - bias_var, 0)) / 4
        result$se_effect <- sqrt(jab_var)
    
          
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










########################## INFJACK Function ####################################
rInfJack = function(pred, inbag, calibrate = TRUE, used.trees = NULL) {
  # original: https://github.com/swager/randomForestCI/blob/master/R/infinitesimalJackknife.R
  
  if (is.null(used.trees)) {
    used.trees = 1:ncol(inbag)
  }
  pred = pred[, used.trees, drop=FALSE]
  
  # Check if sampling without replacement
  no.replacement = (max(inbag) == 1)
  
  # Extract tree-wise predictions and variable counts from random forest
  B = length(used.trees)
  n = nrow(inbag)
  s = sum(inbag) / ncol(inbag)
  
  y.hat = rowMeans(pred)
  pred.centered = pred - rowMeans(pred)
  
  N = Matrix::Matrix(inbag[, used.trees], sparse = TRUE)
  N.avg = Matrix::rowMeans(N)
  
  # Compute raw infinitesimal jackknife
  if (as.numeric(B)^2 > as.numeric(n) * as.numeric(nrow(pred))) {
    
    C = Matrix::tcrossprod(N, pred.centered) -
      Matrix::Matrix(N.avg, nrow(N), 1) %*%
      Matrix::Matrix(rowSums(pred.centered), 1, nrow(pred.centered))
    raw.IJ = Matrix::colSums(C^2) / B^2
    
  } else {
    # Faster implementation when n is large. Uses the fact that
    # colSums((A - B)^2) = T1 - 2 * T2 + T3,
    # where T1 = diag(A'A), T2 = diag(B'A), and T3 = diag(B'B)
    
    NTN = Matrix::crossprod(N, N)
    NTNPT_T = Matrix::tcrossprod(pred.centered, NTN)
    T1 = Matrix::rowSums(pred.centered * NTNPT_T)
    
    RS = rowSums(pred.centered)
    NbarTN = Matrix::crossprod(N.avg, N)
    T2 = RS * Matrix::tcrossprod(NbarTN, pred.centered)
    
    T3 = sum(N.avg^2) * RS^2
    raw.IJ = as.numeric(T1 - 2 * T2 + T3) / B^2
    
  }
  
  # Apply Monte Carlo bias correction
  N.var = mean(Matrix::rowMeans(N^2) - Matrix::rowMeans(N)^2)
  boot.var = rowSums(pred.centered^2) / B
  bias.correction = n * N.var * boot.var / B
  vars = raw.IJ - bias.correction
  
  # Finite sample correction
  if (no.replacement) {
    variance.inflation = 1 / (1 - mean(inbag))^2
    vars = variance.inflation * vars
  }
  
  results = data.frame(y.hat=y.hat, var.hat=vars)
  
  if (isTRUE(calibrate) && nrow(results) <= 20) {
    calibrate = FALSE
    warning("Sample size <=20, no calibration performed.")
  }
  
  # If appropriate, calibrate variance estimates; this step in particular
  # ensures that all variance estimates wil be positive.
  
  if (calibrate) {
    # Compute variance estimates using half the trees
    calibration.ratio = 2
    n.sample = ceiling(B / calibration.ratio)
    results.ss = rInfJack(pred, inbag, calibrate = FALSE, used.trees = sample(used.trees, n.sample))
    
    # Use this second set of variance estimates to estimate scale of Monte Carlo noise
    sigma2.ss = mean((results.ss$var.hat - results$var.hat)^2)
    delta = n.sample / B
    sigma2 = (delta^2 + (1 - delta)^2) / (2 * (1 - delta)^2) * sigma2.ss
    
    # Use Monte Carlo noise scale estimate for empirical Bayes calibration
    results = tryCatch(
      expr = {
        vars.calibrated = calibrateEB(vars, sigma2)
        results$var.hat = vars.calibrated
        results
      }, 
      error = function(e) {
        warning(sprintf("Calibration failed with error:\n%sFalling back to non-calibrated variance estimates.", e))
        results = rInfJack(pred, inbag, calibrate = FALSE, used.trees = used.trees)
        return(results)
      }
    )
  }
  
  return(results)
}

# Fit an empirical Bayes prior in the hierarchical model
#     mu ~ G, X ~ N(mu, sigma^2)
#
# @param X a vector of observations
# @param sigma noise estimate
# @param p tuning parameter -- number of parameters used to fit G
# @param nbin tuning parameter -- number of bins used for discrete approximation
# @param unif.fraction tuning parameter -- fraction of G modeled as "slab"
#
# @return posterior density estimate g
#
# @section References:
# For more details about "g-estimation", see: B Efron. Two modeling strategies for
# empirical Bayes estimation. Stat. Sci., 29: 285-301, 2014.
# @author Stefan Wager
gfit = function(X, sigma, p = 2, nbin = 1000, unif.fraction = 0.1) {
  
  xvals = seq(min(min(X) - 2 * sd(X), 0), max(max(X) + 2 * sd(X), sd(X)), length.out = nbin)
  binw = xvals[2] - xvals[1]
  
  zero.idx = max(which(xvals <= 0))
  noise.kernel = dnorm(xvals / sigma) * binw / sigma
  
  if (zero.idx > 1) {
    noise.rotate = noise.kernel[c(zero.idx:length(xvals), 1:(zero.idx - 1))]
  } else {
    noise.rotate = noise.kernel
  }
  
  XX = sapply(1:p, function(j) xvals^j * as.numeric(xvals >= 0))
  neg.loglik = function(eta) {
    g.eta.raw = exp(XX %*% eta) * as.numeric(xvals >= 0)
    if ((sum(g.eta.raw) == Inf) | (sum(g.eta.raw) <= 100 * .Machine$double.eps)) {
      return (1000 * (length(X) + sum(eta^2)))
    }
    g.eta.main = g.eta.raw / sum(g.eta.raw)
    g.eta = (1 - unif.fraction) * g.eta.main +
      unif.fraction * as.numeric(xvals >= 0) / sum(xvals >= 0)
    f.eta = convolve(g.eta, noise.rotate)
    sum(approx(xvals, -log(pmax(f.eta, 0.0000001)), X)$y)
  }
  
  eta.hat = nlm(neg.loglik, rep(-1, p))$estimate
  g.eta.raw = exp(XX %*% eta.hat) * as.numeric(xvals >= 0)
  g.eta.main = g.eta.raw / sum(g.eta.raw)
  g.eta = (1 - unif.fraction) * g.eta.main +
    unif.fraction * as.numeric(xvals >= 0) / sum(xvals >= 0)
  
  return(data.frame(x=xvals, g=g.eta))
}

# Bayes posterior estimation with Gaussian noise
#
# @param x0 an obsevation
# @param g.est a prior density, as returned by gfit
# @param sigma noise estimate
#
# @return posterior estimate E[mu | x0]
# @author Stefan Wager
gbayes = function(x0, g.est, sigma) {
  Kx = dnorm((g.est$x - x0) / sigma)
  post = Kx * g.est$g
  post = post / sum(post)
  sum(post * g.est$x)
}

# Empirical Bayes calibration of noisy variance estimates
#
# @param vars list of variance estimates
# @param sigma2 estimate of the Monte Carlo noise in vars
#
# @return calibrated variance estimates
# @author Stefan Wager
calibrateEB = function(vars, sigma2) {
  
  if(sigma2 <= 0 | min(vars) == max(vars)) {
    return(pmax(vars, 0))
  }
  
  sigma = sqrt(sigma2)
  eb.prior = gfit(vars, sigma)
  
  if (length(vars >= 200)) {
    # If there are many test points, use interpolation to speed up computations
    calib.x = unique(quantile(vars, q = seq(0, 1, by = 0.02)))
    calib.y = sapply(calib.x, function(xx) gbayes(xx, eb.prior, sigma))
    calib.all = approx(x=calib.x, y=calib.y, xout=vars)$y
  } else {
    calib.all = sapply(vars, function(xx) gbayes(xx, eb.prior, sigma))
  }
  return(calib.all)
}

