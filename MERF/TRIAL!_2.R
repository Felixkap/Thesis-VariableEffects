library(mvtnorm)
library(MASS)
library(MixMatrix) # for compound symmetry matrix (CSgenerate)
library(fabricatr)
library(lmerTest)
library(tictoc)
library(effects)
library(pdp)
library(data.table)
library(ggplot2)
library(nlme)
library(mgsub)
library(ggeffects)

#Simulate data with certain parameters 
simulate.data.gaussian <- function( #only random intercept, no random slope, only possible to generate standard normal data
  nr.group, #number of groups
  group.size, #number of observations per group
  nrX.between, #number of predictors between group
  nrX.within, #number of predictors within group
  cor = 0, #rho (correlation between predictors)
  sigma.resid = 1, #size residual
  sigma.b = 1, # size of random effect, i.e. random intercept of group
  fixed.formula # fixed effects formula
){
  X.groups <- rmvnorm(nr.group, mean = rep(0, nrX.between), sigma = CSgenerate(nrX.between, cor))
  X.b <- X.groups[rep(seq_len(nrow(X.groups)), each = group.size), ] #for each group repeat the group variable
  
  X.w <- rmvnorm(nr.group * group.size, mean = rep(0, nrX.within), sigma = CSgenerate(nrX.within, cor))
  # merge the two sets of variables
  X = data.frame(X=cbind(X.b, X.w))
  
  #generate the response
  X$y <- (eval(parse( text = fixed.formula), X) + rnorm(nr.group * group.size, mean = 0, sd = sigma.resid) 
          + rep(rnorm(nr.group, mean = 0, sd = sigma.b), each=group.size))
  # last two terms add noise and random intercept (per group) respectively
  X$grp <- rep(1:nr.group, each=group.size)
  
  #now the dataset is complete
  data = X
  return( data )
}



save.simulated.data <- function(input #(N, #number of simulations (for each scenario)
                                #nr.group, #number of groups
                                #group.size, #number of observations per group
                                #nrX.between, #number of predictors between group
                                #nrX.within, #number of predictors within group
                                #cor = 0, #rho (correlation between predictors)
                                #sigma.resid = 1, #size residual
                                #sigma.b = 1, # size of random effect, i.e. random intercept of group
                                #fixed.formula, #fixed effects formula
                                #scenario number
                                #seed every run another seed, makes it reproducible
                                #batch number
){
  #initialize list with dataframes for random effects estimates, i.e. residual and r intercept
  if (length(input) != 12) print('Wrong number of input parameters')
  
  N <- as.numeric(input[1])
  nr.group<-as.numeric(input[2])
  group.size<-as.numeric(input[3]) 
  nrX.between <- as.numeric(input[4])
  nrX.within <- as.numeric(input[5])
  cor <- as.numeric(input[6])
  sigma.resid <- as.numeric(input[7])
  sigma.b<- as.numeric(input[8])
  fixed.formula <- input[9]
  seed <- as.numeric(input[11])
  scenario <- as.numeric(input[10])
  batch <- as.numeric(input[12])
  
  
  dat_list = vector(mode = 'list', length=N)
  for(i in 1:N){ # N is number of simulation, smaller than 50 is necessary to avoid duplicates
    #SET SEED
    set.seed(seed = seed)
    data <- simulate.data.gaussian(nr.group,group.size,nrX.between, nrX.within,cor,sigma.resid,sigma.b,fixed.formula)
    dat_list[[i]] <- data
  }
  return(dat_list)
}
  
  
  
  


# main function, this function has as input the simulation scenarios, then generates the data and runs the models and writes the results down
run.models <- function( input #(N, #number of simulations (for each scenario)
                        #nr.group, #number of groups
                        #group.size, #number of observations per group
                        #nrX.between, #number of predictors between group
                        #nrX.within, #number of predictors within group
                        #cor = 0, #rho (correlation between predictors)
                        #sigma.resid = 1, #size residual
                        #sigma.b = 1, # size of random effect, i.e. random intercept of group
                        #fixed.formula, #fixed effects formula
                        #scenario number
                        #seed every run another seed, makes it reproducible
                        #batch number
){
  #initialize list with dataframes for random effects estimates, i.e. residual and r intercept
  if (length(input) != 12) print('Wrong number of input parameters')
  
  N <- as.numeric(input[1])
  nr.group<-as.numeric(input[2])
  group.size<-as.numeric(input[3]) 
  nrX.between <- as.numeric(input[4])
  nrX.within <- as.numeric(input[5])
  cor <- as.numeric(input[6])
  sigma.resid <- as.numeric(input[7])
  sigma.b<- as.numeric(input[8]) ;
  fixed.formula <- input[9]
  seed <- as.numeric(input[11])
  scenario <- as.numeric(input[10])
  batch <- as.numeric(input[12])
  
  models <- c('LM', 'LME')
  RES <- c("resid", "int")
  RE.estimates <- data.frame('LM_s_resid' = rep(0,N), 'LME_s_resid' = rep(0,N),'
                             LME_s_int' = rep(0,N), 'LME_s_int_hat'=rep(0,N),
                             'seed' = rep(0,N))
  FE.estimates <- pf.list2(N, models[1:2], 1 + 2*(nrX.between+nrX.within)) #for intercept and every first and second order term
  ncol.fixed <- 2+4*(nrX.between+nrX.within)
  FE.estimates <- data.frame(matrix(NA, nrow=N, ncol=ncol.fixed))
  
  run.times <- as.data.frame(matrix(NA, ncol=2, nrow=N))
  colnames(run.times) <- c('LM', 'LME')
  
  #dataframe for mse accuracies
  ACC <- data.frame("LM"=rep(0, N), "LME"=rep(0, N), 'seed' = rep(0, N))
  #to store convergence
  conv <- matrix(1, ncol=2, nrow=N) #1 if it converged, 0 otherwise (added at convergence check)
  if(nrX.between == 10){fixed.formula <- mgsub(c("X.3", "X.4", "X.5"), c("X.11","X.12", "X.13"), string = fixed.formula)}
  
  
  for(i in 1:N){ # N is number of simulation, smaller than 50 is necessary to avoid duplicates
    #SET SEED
    set.seed(seed = seed)
    data <- simulate.data.gaussian(nr.group,group.size,nrX.between, nrX.within,cor,sigma.resid,sigma.b,fixed.formula)
    
    data.wg <- data[, - ncol(data)] #so get rid of the last column, i.e. the group number
    
    if ( i == 1 ) { #only do this first time, create model formula (for lm and lme)
      formula <- as.formula(
        paste('y ~',paste('poly(',colnames(data.wg[-length(colnames(data.wg))]),',2, raw = T)',collapse = ' + ')
        )) #formula that's fitted for LM and LME models, it only include first and second order terms, no tf or interactions
      formula2 <- as.formula(
        paste('y ~',paste('poly(',colnames(data.wg[-length(colnames(data.wg))]),',2)',collapse = ' + ')
        )) #formula that's fitted for LME pdp plot (marginal effect plot. somehow needed, unclear why)
      
      effect.estimates1_1 <- make.FE.dataframe(data.wg, c('LM', 'LME'), N)
      effect.estimates0_1 <- make.FE.dataframe(data.wg, c('LM', 'LME'), N)
    }
    
    ####### LM #############
    t <- system.time(LM <- lm(formula, data=data.wg))
    run.times[i, 'LM'] <- unname(t[3])
    slm <- summary(LM)
    
    ####### LME ############
    
    t <- system.time(LME <- lmer(update(formula,    ~ . + (1|grp)), data = data))#add Random Intercept for mixed-effects model
    run.times[i, 'LME'] <- unname(t[3])
    slme <- summary(LME)
    
    #------------------------------------#
    ######## Fixed effects ###############
    #------------------------------------#
    
    # For LM and LME store all fixed effects coefficients, and at the end compute mean, sd, mse etc.
    FE.estimates[i, 1:(1 + 2*(nrX.between+nrX.within))] <- unname(coef(slm)[,1])
    FE.estimates[i, (2 + 2*(nrX.between+nrX.within)):(2 + 4*(nrX.between+nrX.within))] <- unname(coef(slme)[,1])
    # For all 6:combined dependence plots, only of last fitted models (so after this for-loop)
    
    s <- 1; t<- nrX.between + nrX.within
    for(model in list(LM, LME)){
      effect.estimates0_1[i,s:t] <- var.effect(model, data, value=1, sym=F)
      effect.estimates1_1[i,s:t] <- var.effect(model, data, value=1, sym=T)
      s <- s+nrX.between+nrX.within; t <- t+nrX.between+nrX.within
      
      #if fixed formula has step function, then compute for that X the effect differently, i.e. between 0.25 and 0.75
      if(scenario > 64){
        x <- 4
        if(nrX.between==10) { x <- 12}
        adj.data_min <- data; adj.data_plus <- data
        adj.data_plus[,x] <-  rep(0.75, length(data[,x]))
        adj.data_min[,x] <-  rep(0.25, length(data[,x]))
        pred.hat.y_min <- predict(model, newdata=adj.data_min)
        pred.hat.y_plus <- predict(model, newdata=adj.data_plus)
        effect.estimates0_1[i,x] <- (mean(pred.hat.y_plus) - mean(pred.hat.y_min)) / (0.5)
      }
    }
    
    
    #------------------------------------#
    ######## Random effects ##############
    #------------------------------------#
    
    # Store random intercept estimate (only applicable for MERF and LME)
    # Store residual sigma_resid^2 estimate (non-applicable for RF)
    RE.estimates[i,'LM_s_resid'] <- slm$sigma
    RE.estimates[i,'LME_s_resid'] <- slme$sigma
    
    #Store random intercept estimates
    RE.estimates[i, 'LME_s_int'] <- unname(sqrt(unlist(slme$varcor)))
    
    #also store the sd group estimates, which is quasi estimate
    RE.estimates[i, 'LME_s_int_hat'] <- sd(unname(unlist(ranef(LME))))
    RE.estimates[i, 'seed'] <- seed
    #------------------------------------#
    ###### Prediction accuracy ###########
    #------------------------------------#
    
    # Compute MSE of fitted values, MSE lm
    ACC[i, 'LM'] <- mean(slm$residuals^2)
    ACC[i, 'LME'] <- mean(slme$residuals^2)
    ACC[i, 'seed'] <- seed
    
    seed <- seed + 1
  }
  #j<-1
  #MERF.lDB <- data.wg #this is (somehow) needed for the pdp plots of the RF part of the MERF models
  
  #for(var in plot.var){
  #  plots[j] <- effect.plot(list(LM, MERF1$fit.rf, MERF2$fit.rf, RF1, RF2, LME), var, data.wg)
  #  j<- j+1
  #}
  
  # plot.var <- colnames(data)[1:(nrX.between+nrX.within)]
  # plots <- rep(0, length(plot.var)) #to store the ggplot pdp plots
  # plots <- lapply(plot.var, effect.plot, models = list(LM, MERF1$fit.rf, MERF2$fit.rf, RF1, RF2, LME2))
  
  
  
  input <- c('start seed' = input[11], 'scenario' = scenario, 'batch'=batch ,'N' = N, 'nr.group'= nr.group,
             'group.size' = group.size, 'nrX.between' = nrX.between, 'nrX.within' = nrX.within,
             'cor' = cor, 'sigma.resid' = sigma.resid, 'sigma.b' = sigma.b, 'formula' = input[9])
  output <- list(ACC, FE.estimates, run.times, input, effect.estimates0_1, effect.estimates1_1)
  names(output) <- c("ACC", "FE.estimates", "run.times", "input", 'effects1', 'effects2')
  
  name <- paste('scenario','_' ,scenario,'_',batch, sep='')
  ############--------------------###########
  output
}

# this function computes the fixed effects, in a partial-dependence-plot-like way
var.effect <- function(model, data, value, sym = T){ # difference between (-point) and (+point) or point and 0
  effect <- c()
  if(sym){
    for(j in 1:(ncol(data)-2)){
      adj.data_min <- data
      adj.data_plus <- data
      adj.data_plus[,j] <-  rep(abs(value), length(data[,j]))
      adj.data_min[,j] <-  rep(-abs(value), length(data[,j]))
      pred.hat.y_min <- predict(model, newdata=adj.data_min)
      pred.hat.y_plus <- predict(model, newdata=adj.data_plus)
      effect[j] <- (mean(pred.hat.y_plus) - mean(pred.hat.y_min)) / (2 * abs(value))
    }
  }
  if(!sym){
    for(j in 1:(ncol(data)-2)){
      adj.data_val <- data
      adj.data_0 <- data
      adj.data_val[,j] <-  rep(value, length(data[,j]))
      adj.data_0[,j] <-  rep(0, length(data[,j]))
      pred.hat.y_val <- predict(model, newdata=adj.data_val)
      pred.hat.y_0 <- predict(model, newdata=adj.data_0)
      effect[j] <- (mean(pred.hat.y_val) - mean(pred.hat.y_0)) / value
    }
  }
  effect
}


#some short function for generating data frames
MSE.lm <- function(LM) sum(LM$residuals^2)/LM$df.residual

pf.list <- function(N, names, parameters) {
  df <- as.data.frame(matrix(NA, nrow = N, ncol = length(parameters)))
  colnames(df) <- parameters
  setNames(replicate(length(names),df,simplify=FALSE), names)
}

pf.list2 <- function(N, names, ncol) {
  df <- as.data.frame(matrix(NA, nrow = N, ncol = ncol ))
  setNames(replicate(length(names),df,simplify=FALSE), names)
}


make.FE.dataframe <- function(data.wg, models, N){
  nrX <- ncol(data.wg) - 1
  df <- as.data.frame(matrix(NA,nrow=N, ncol= nrX*length(models)))
  s <- 1
  t <- nrX
  for(model in models){
    colnames(df)[s:t] <- paste(model,"_",colnames(data.wg[-length(colnames(data.wg))]), sep='')
    s <- s+nrX ; t <- t + nrX
  }
  df
}












