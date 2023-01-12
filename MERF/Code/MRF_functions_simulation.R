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
  print('X.groups')
  print(X.groups)
  X.b <- X.groups[rep(seq_len(nrow(X.groups)), each = group.size), ] #for each group repeat the group variable
  print('X.b')
  print(X.b)
  
  X.w <- rmvnorm(nr.group * group.size, mean = rep(0, nrX.within), sigma = CSgenerate(nrX.within, cor))
  print('X.w')
  print(X.w)
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
  
  N <- as.numeric(input[1]); nr.group<-as.numeric(input[2]);  group.size<-as.numeric(input[3]); 
  nrX.between <- as.numeric(input[4]); nrX.within <- as.numeric(input[5])
  cor <- as.numeric(input[6]); sigma.resid <- as.numeric(input[7]); sigma.b<- as.numeric(input[8]) ;
  fixed.formula <- input[9]; seed <- as.numeric(input[11]); scenario <- as.numeric(input[10])
  batch <- as.numeric(input[12])
  
  models <- c('LM', 'LME', 'MERF1', 'MERF2')
  RES <- c("resid", "int")
  RE.estimates <- data.frame('LM_s_resid' = rep(0,N), 'LME_s_resid' = rep(0,N), 'MERF1_s_resid' = rep(0,N),
                             'MERF2_s_resid' = rep(0,N), 'LME_s_int' = rep(0,N), 'MERF1_s_int' = rep(0,N),
                             'MERF2_s_int' = rep(0,N), 'LME_s_int_hat'=rep(0,N), 'MERF1_s_int_hat'=rep(0,N), 
                             'MERF2_s_int_hat'=rep(0,N), 'seed' = rep(0,N))
  FE.estimates <- pf.list2(N, models[1:2], 1 + 2*(nrX.between+nrX.within)) #for intercept and every first and second order term
  ncol.fixed <- 2+4*(nrX.between+nrX.within)
  FE.estimates <- data.frame(matrix(, nrow=N, ncol=ncol.fixed))
  
  run.times <- as.data.frame(matrix(, ncol=6, nrow=N))
  colnames(run.times) <- c( 'RF1', 'RF2', 'LM', 'LME', 'MERF1', 'MERF2')
  
  #dataframe for mse accuracies
  ACC <- data.frame("LM"=rep(0, N), "LME"=rep(0, N), "RF1_oob" = rep(0, N), "RF1" = rep(0, N),
                    "RF2_oob" = rep(0, N), "RF2" = rep(0, N), "MERF1_oob" =rep(0, N), 'MERF1'=rep(0, N),
                    "MERF2_oob" =rep(0, N), 'MERF2'=rep(0, N), 'seed' = rep(0, N))
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
      
      effect.estimates1_1 <- make.FE.dataframe(data.wg, c('RF1', 'RF2', 'LM', 'LME', 'MERF1', 'MERF2'), N)
      effect.estimates0_1 <- make.FE.dataframe(data.wg, c('RF1', 'RF2', 'LM', 'LME', 'MERF1', 'MERF2'), N)
    }
    ####### RF #############
    
    m <- max(floor((nrX.between + nrX.within)/3), 1) #default mtry
    m2 <- ceiling(0.8*(nrX.between + nrX.within))
    
    t <- system.time(RF1 <- randomForest(y ~ ., mtry=m, nodesize=5,  data= data.wg, ntree = 500, importance = TRUE))
    run.times[i, 'RF1'] <- unname(t[3])
    t <- system.time(RF2 <- randomForest(y ~ ., mtry=m2, nodesize=1, data= data.wg, ntree = 500, importance = TRUE))
    run.times[i, 'RF2'] <- unname(t[3])
    
    ####### LM #############
    t <- system.time(LM <- lm(formula, data=data.wg))
    run.times[i, 'LM'] <- unname(t[3])
    slm <- summary(LM)
    
    ####### LME ############
    
    t <- system.time(LME <- lmer(update(formula,    ~ . + (1|grp)), data = data))#add Random Intercept for mixed-effects model
    run.times[i, 'LME'] <- unname(t[3])
    slme <- summary(LME)
    
    #LME2 <- lmer(update(formula2,    ~ . + (1|grp)), data = data)
    
    
    ####### MERF models ###########
    
    t <- system.time(
      MERF1 <- MERF.fit(data, mtry = m, nodesize = 5, crit=1e-3,  nr.group, group.size, nrX.within, nrX.between, max.niter = 15))
    run.times[i, 'MERF1'] <- unname(t[3])
    t <-system.time(
      MERF2 <- MERF.fit(data, mtry = m2, nodesize = 2, crit=1e-3,  nr.group, group.size, nrX.within, nrX.between, max.niter = 15))
    run.times[i, 'MERF2'] <- unname(t[3])
    
    #Check if it converged
    if(all(is.na(MERF1$convergence.iter))) conv[i,1] <- 0
    if(all(is.na(MERF2$convergence.iter))) conv[i,2] <- 0
    
    
    #Next, combine output, model performace wrt fixed effects, 
    #random effect (residual and intercept) amd prediction accuracy (MSE)
    
    #------------------------------------#
    ######## Fixed effects ###############
    #------------------------------------#
    
    # For LM and LME store all fixed effects coefficients, and at the end compute mean, sd, mse etc.
    FE.estimates[i, 1:(1 + 2*(nrX.between+nrX.within))] <- unname(coef(slm)[,1])
    FE.estimates[i, (2 + 2*(nrX.between+nrX.within)):(2 + 4*(nrX.between+nrX.within))] <- unname(coef(slme)[,1])
    # For all 6:combined dependence plots, only of last fitted models (so after this for-loop)
    
    s <- 1; t<- nrX.between + nrX.within
    for(model in list(RF1, RF2, LM, LME, MERF1$fit.rf, MERF2$fit.rf)){
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
    RE.estimates[i, 'MERF1_s_resid'] <- sqrt(MERF1$`sigma.sq[r]`)
    RE.estimates[i, 'MERF2_s_resid'] <- sqrt(MERF2$`sigma.sq[r]`)
    
    #Store random intercept estimates
    RE.estimates[i, 'LME_s_int'] <- unname(sqrt(unlist(slme$varcor)))
    RE.estimates[i, 'MERF1_s_int'] <- sqrt(unlist(last(MERF1$D, n=1)))
    RE.estimates[i, 'MERF2_s_int'] <- sqrt(unlist(last(MERF2$D, n=1)))
    
    #also store the sd group estimates, which is quasi estimate
    RE.estimates[i, 'LME_s_int_hat'] <- sd(unname(unlist(ranef(LME))))
    RE.estimates[i, 'MERF1_s_int_hat'] <- sd(MERF1$`bi[[r]]`)
    RE.estimates[i, 'MERF2_s_int_hat'] <- sd(MERF2$`bi[[r]]`)
    RE.estimates[i, 'seed'] <- seed
    #------------------------------------#
    ###### Prediction accuracy ###########
    #------------------------------------#
    
    # Compute MSE of fitted values, MSE lm
    ACC[i, 'LM'] <- mean(slm$residuals^2)
    ACC[i, 'LME'] <- mean(slme$residuals^2)
    ACC[i, 'RF1_oob'] <- last(RF1$mse, n=1)
    ACC[i, 'RF1'] <- mean((predict(RF1, data.wg)-data.wg[,'y'])^2)
    ACC[i, 'RF2_oob'] <- last(RF2$mse, n=1)
    ACC[i, 'RF2'] <- mean((predict(RF2, data.wg)-data.wg[,'y'])^2)
    
    ACC[i, 'MERF1_oob'] <- mean((MERF1$fit.rf$predicted + rep(MERF1$`bi[[r]]`, each=group.size)-data[,'y'])^2)
    ACC[i, 'MERF1'] <-   mean((predict(MERF1$fit.rf, data.wg) + rep(MERF1$`bi[[r]]`, each=group.size)-data[,'y'])^2)
    
    ACC[i, 'MERF2_oob'] <- mean((MERF2$fit.rf$predicted + rep(MERF2$`bi[[r]]`, each=group.size)-data[,'y'])^2)
    ACC[i, 'MERF2'] <-   mean((predict(MERF2$fit.rf, data.wg) + rep(MERF2$`bi[[r]]`, each=group.size)-data[,'y'])^2)
    
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
  output <- list(ACC, RE.estimates, FE.estimates, run.times, input, conv, effect.estimates0_1, effect.estimates1_1)
  names(output) <- c("ACC", "RE.estimates", "FE.estimates", "run.times", "input", "conv",
                     'effects1', 'effects2')
  
  name <- paste('scenario','_' ,scenario,'_',batch, sep='')
  ###### write output to text-files ##############
  write.table( input,
               paste("C:/Users/rtn/Documents/simulationResults/ACC/",name,sep=''), 
               append = F,
               sep = "\t",
               row.names = F,
               col.names = T)
  write( "\n Mean squared error residuals \n", file = paste("C:/Users/rtn/Documents/simulationResults/ACC/",name,sep=''), append = T)
  write.table( ACC,
               paste("C:/Users/rtn/Documents/simulationResults/ACC/",name,sep=''), 
               append = T,
               sep = "\t",
               row.names = F,
               col.names = T)
  
  write.table( input,
               paste("C:/Users/rtn/Documents/simulationResults/FES/",name,sep=''), 
               append = F,
               sep = "\t",
               row.names = F,
               col.names = T)
  write( "\n Fixed effects estimates \n", file = paste("C:/Users/rtn/Documents/simulationResults/FES/",name,sep=''), append = T)
  write.table( FE.estimates,
               paste("C:/Users/rtn/Documents/simulationResults/FES/",name,sep=''), 
               append = T,
               sep = "\t",
               row.names = F,
               col.names = T)
  
  write.table( input,
               paste("C:/Users/rtn/Documents/simulationResults/RES/",name,sep=''), 
               append = F,
               sep = "\t",
               row.names = F,
               col.names = T)
  write( "\n Random effects estiamtes \n", file = paste("C:/Users/rtn/Documents/simulationResults/RES/",name,sep=''), append = T)
  write.table( RE.estimates,
               paste("C:/Users/rtn/Documents/simulationResults/RES/",name,sep=''), 
               append = T,
               sep = "\t",
               row.names = F,
               col.names = T)
  
  write.table( input,
               paste("C:/Users/rtn/Documents/simulationResults/convergence/",name,sep=''), 
               append = F,
               sep = "\t",
               row.names = F,
               col.names = T)
  write( "\n convergence of MERF models \n", file = paste("C:/Users/rtn/Documents/simulationResults/convergence/",name,sep=''), append = T)
  write.table( conv,
               paste("C:/Users/rtn/Documents/simulationResults/convergence/",name,sep=''), 
               append = T,
               sep = "\t",
               row.names = F,
               col.names = T)
  
  write.table( input,
               paste("C:/Users/rtn/Documents/simulationResults/runtimes/",name,sep=''), 
               append = F,
               sep = "\t",
               row.names = F,
               col.names = T)
  write( "\n Run times \n", file = paste("C:/Users/rtn/Documents/simulationResults/runtimes/",name,sep=''), append = T)
  write.table( run.times,
               paste("C:/Users/rtn/Documents/simulationResults/runtimes/",name,sep=''), 
               append = T,
               sep = "\t",
               row.names = F,
               col.names = T)
  
  write.table( input,
               paste("C:/Users/rtn/Documents/simulationResults/effectEstimates1/",name,sep=''), 
               append = F,
               sep = "\t",
               row.names = F,
               col.names = T)
  write( "\n Fixed effect estimates between -1 and 1 \n", file = paste("C:/Users/rtn/Documents/simulationResults/effectEstimates1/",name,sep=''), append = T)
  write.table( effect.estimates1_1,
               paste("C:/Users/rtn/Documents/simulationResults/effectEstimates1/",name,sep=''), 
               append = T,
               sep = "\t",
               row.names = F,
               col.names = T)
  
  write.table( input,
               paste("C:/Users/rtn/Documents/simulationResults/effectEstimates2/",name,sep=''), 
               append = F,
               sep = "\t",
               row.names = F,
               col.names = T)
  write( "\n Fixed effect estimates between 0 and 1 (or 0.25 and 0.75 for step function [X.4/X.11]) \n", file = paste("C:/Users/rtn/Documents/simulationResults/effectEstimates2/",name,sep=''), append = T)
  write.table( effect.estimates0_1,
               paste("C:/Users/rtn/Documents/simulationResults/effectEstimates2/",name,sep=''), 
               append = T,
               sep = "\t",
               row.names = F,
               col.names = T)
  ############--------------------###########
  output
}

# this function computes the fixed effects, in a partial-dependence-plot-like way
var.effect <- function(model, data, value, sym = T){ # difference between (-point) and (+point) or point and 0
  effect <- c()
  if(sym){
    for(j in 1:(ncol(data)-2)){
      adj.data_min <- data; adj.data_plus <- data
      adj.data_plus[,j] <-  rep(abs(value), length(data[,j]))
      adj.data_min[,j] <-  rep(-abs(value), length(data[,j]))
      pred.hat.y_min <- predict(model, newdata=adj.data_min)
      pred.hat.y_plus <- predict(model, newdata=adj.data_plus)
      effect[j] <- (mean(pred.hat.y_plus) - mean(pred.hat.y_min)) / (2 * abs(value))
    }
  }
  if(!sym){
    for(j in 1:(ncol(data)-2)){
      adj.data_val <- data; adj.data_0 <- data
      adj.data_val[,j] <-  rep(value, length(data[,j]))
      adj.data_0[,j] <-  rep(0, length(data[,j]))
      pred.hat.y_val <- predict(model, newdata=adj.data_val)
      pred.hat.y_0 <- predict(model, newdata=adj.data_0)
      effect[j] <- (mean(pred.hat.y_val) - mean(pred.hat.y_0)) / value
    }
  }
  effect
}

#function used to fit the MERF
MERF.fit <- function(data, mtry = 4, nodesize = 5, crit=1e-4, nr.group, group.size, nrX.within, nrX.between,  max.niter=15){
  xnam = names(data)[1:(nrX.between+nrX.within)] # the predictor variables needed for MERF fit
  ni = c(rep(group.size, nr.group))
  
  ### Yi ###
  Yi = split(data$y, data$grp)
  #str(Yi)
  
  ### Zi ###
  ## if q = 1 (i.e. random intercept)
  Zi.ri = list(); length(Zi.ri)= length(ni)
  for(i in 1:length(Zi.ri))    
    Zi.ri[[i]] <- matrix(rep(1,ni[i]),nrow=ni[i],ncol=1)             	         
  Zi = Zi.ri
  
  #m <- max(floor((nrX.between + nrX.within)/3), 1)
  
  ### MERF application ### ----
  MERF.fit <- MERF(
    xnam=xnam
    ,MERF.lDB = data			
    ,ni = ni
    ,Zi = Zi
    ,Yi = Yi
    ,grp = data$grp
    
    ,ntree= 500
    ,mtry= mtry		
    ,nodesize= nodesize
    
    ,sigmasqzero = NULL		
    ,Dzero = NULL			
    ,bizero = NULL
    
    ,F.niter = 0	
    ,max.niter = max.niter
    ,smallest.Jump.allowed = crit
    
    ,verbose= TRUE
  )
  return(MERF.fit) # has output: -- fitted RF -- sigma_squared (resid) -- sigma_INT (Random effect)
}


#some short function for generating data frames
MSE.lm <- function(LM) sum(LM$residuals^2)/LM$df.residual

pf.list <- function(N, names, parameters) {
  df <- as.data.frame(matrix(, nrow = N, ncol = length(parameters)))
  colnames(df) <- parameters
  setNames(replicate(length(names),df,simplify=FALSE), names)
}

pf.list2 <- function(N, names, ncol) {
  df <- as.data.frame(matrix(, nrow = N, ncol = ncol ))
  setNames(replicate(length(names),df,simplify=FALSE), names)
}


make.FE.dataframe <- function(data.wg, models, N){
  nrX <- ncol(data.wg) - 1
  df <- as.data.frame(matrix(,nrow=N, ncol= nrX*length(models)))
  s <- 1
  t <- nrX
  for(model in models){
    colnames(df)[s:t] <- paste(model,"_",colnames(data.wg[-length(colnames(data.wg))]), sep='')
    s <- s+nrX ; t <- t + nrX
  }
  df
}

#The actual R code for MERF, as provided by one of the authors

#############################################################################################################
#													#####
##				Ahlem HAJJEM, Fran?ois BELLAVANCE, and Denis LAROCQUE		 	 ####
###				Department of Management Sciences HEC Montr?al			  	  ###
####				HEC Montr?al, 3000, chemin de la C?te-Sainte-Catherine,	   	 	   ##
#####				Montr?al, QC, Canada H3T 2A7

###	Revised: 27st september 2017	###
#############################################################################################################

###########################################
###	Description of MERF arguments	###
###########################################

#xnam			#A charcter vector of p columns, corresponding to the names of the p fixed effect covariates (as they appear in the learning dataset).

#MERF.lDB 		#The learning dataset: a dataframe of N rows (i.e. N level-one observations) and (p+2) columns, 
#where the column "cluster.id" corresponds to the variable uniquely identifying the n clusters (or level-two obervations), 
#the column "Y" corresponds to the name of the continuous response variable,
#and the other p columns correspond to the p fixed effects covariates.	

#ni			#A vector of n columns, where ni[i] corresponds to the size of cluster i, 
#for i = 1, ..., n (ATTENTION: should keep the same order as in MERF.lDB).

#Zi			#A list of n matrices Zi[[i]] of dimension (ni[i] X q), where  
#Zi[[i]] corresponds to the q random effects covariates values of the ni[i] observations nested within cluster i, for i= 1, ..., n. 
#Note that q=1 when there is only a random intercept, i.e. no random effect covariates.

#Yi			#A list of n vectors Yi[i] of ni[i] rows, where 
#Yi[i] corresponds to the response values of the ni[i] observations nested within cluster i, for i= 1, ..., n.

#ntree			#ntree argument of randomForest function, i.e., number of trees to grow (e.g. 300).

#mtry			#mtry argument of randomForest function, i.e., number of variables randomly sampled as candidates at each split (e.g. 3).

#nodesize		#nodesize argument of randomForest function, i.e., minimum size of terminal nodes (e.g.5).

#sigmasqzero = NULL	#Starting value of s2, where the covariance matrix of the errors Ri = s2Ini , for i = 1, ..., n. 
#Default: if( is.null(sigmasqzero) ) sigmasqzero <- 1.

#Dzero = NULL		#Starting value of the covariance matrix of the random effects. 
#Default:
#if( is.null(Dzero) ){
#	Dzero <- diag(0.01,nrow=q,ncol=q)
#}

#bizero = NULL		#Starting values of the random effects: a list of n matrices of q ? 1 unknown vector of random effects. 
#Default:
#if( is.null(bizero) ){
#bizero <- list(); length(bizero)<- n
#	for(i in 1:n) bizero[[i]] <- matrix(0,nrow=q,ncol=1)}
#}

#F.niter		#The number of iterations forced to avoid early stopping (e.g. 100).

#max.niter		#Maximum number of iterations, in addition to the F.niter forced iterations (e.g. 300).

#smallest.Jump.allowed 	#A given small value (e.g. 1e-4). 

#verbose		#Logical. Should R report extra information on progress?

########################################################################################################################
###############################                   ### MERF function ###  	         ##############################
##################		fits a mixed effects random forest of regression trees model		###############
########################################################################################################################


library(randomForest)

MERF  <- function(
  
  xnam
  
  ,MERF.lDB 	
  
  ,ni
  ,Zi
  ,Yi
  ,grp
  
  ,ntree
  ,mtry
  ,nodesize
  
  ,sigmasqzero = NULL
  ,Dzero = NULL
  ,bizero = NULL
  
  ,F.niter
  ,max.niter
  ,smallest.Jump.allowed
  
  ,verbose = TRUE
  
){
  
  ####################
  ####	STEP 0	####
  #################### 
  
  #Memory Allocation and initialization:
  
  #Parameters values 
  n <- length(ni)
  N <- sum(ni)
  q <- dim(Zi[[1]])[2] 	# q=1 in random intercept case
  
  
  #Initial values of sigmasqzero, Dzero, and bizero
  if( is.null(sigmasqzero) ) sigmasqzero <- 1
  else sigmasqzero <- sigmasqzero
  
  if( is.null(Dzero) ){
    Dzero <- diag(1 ,nrow=q,ncol=q)
  }
  else Dzero <- Dzero
  
  if( is.null(bizero) ){
    bizero <- list(); length(bizero)<- n
    for(i in 1:n) bizero[[i]] <- matrix(0,nrow=q,ncol=1)
  }	
  else bizero <- bizero
  
  #iter number
  r <- 1
  if (verbose) 	
    message("MERF iter no: ", r)
  
  #transformed outcome, star.Yi[[r]][[i]], initialized with the original values
  star.Yi <- list() 
  for(i in 1:n){
    star.Yi[[i]] <- Yi[[i]] - Zi[[i]] %*% bizero[[i]]
  }
  
  #one STD random forest
  ######################
  
  MERF.lDB$star.Yi <- unlist(star.Yi) 
  rm(star.Yi) ;  gc(verbose=FALSE)
  
  fit.rf.formula <-  as.formula(paste("star.Yi ~ ", paste(xnam, collapse= "+")))
  
  fit.rf <- randomForest(
    formula=fit.rf.formula 
    ,data=MERF.lDB	
    ,ntree=ntree
    ,mtry = mtry 	
    ,replace=TRUE	
    ,nodesize = nodesize
    ,proximity=FALSE	
  )
  
  #fixed part
  #as vector			
  MERF.lDB$f.pred  <- predict(fit.rf, type="response")			#!!! use the out-of-bag predictions
  #in matrix format
  fixed.pred <- list()	
  fixed.pred <- split(MERF.lDB, grp) 	
  for(i in 1:n)fixed.pred[[i]] <- as.matrix(subset(fixed.pred[[i]] ,select=f.pred), ncol=1)
  
  #random	part
  ############
  #random effects parameters in list format
  bi <- list(list()) ; length(bi) <- r
  for(i in 1:n)bi[[r]][[i]] <- bizero[[i]] 	
  #print("bizero");print(bizero)
  rm(bizero) ; gc(verbose=FALSE)
  
  #level-1 variance component
  #residuals
  epsili <- list()
  for(i in 1:n)
    epsili[[i]] <- Yi[[i]] - fixed.pred[[i]] - Zi[[i]] %*% bi[[r]][[i]]
  
  sigma.sq <- vector(mode="numeric") ;length(sigma.sq) <- r
  sigma.sq[r] <- sigmasqzero
  #print("sigmasqzero") ;print(sigmasqzero)			
  rm(sigmasqzero) ; gc(verbose=FALSE)
  #message("sigmasq of current micro iter", sigma.sq[r] )
  
  #level-2 variance component
  D <- list() ;length(D) <- r
  D[[r]] <- Dzero		#!!!Dzero <- diag(x=0.01, nrow=q, ncol = q)
  #print("Dzero") ;print(Dzero)	
  rm(Dzero) ; gc(verbose=FALSE)
  #message("D of current micro iter: ", D[[r]] )
  
  #level-1 and level-2 variance components (or typical or total variance)
  Vi <- list() 
  inv.Vi <- list(list()) ; length(inv.Vi) <- r
  
  for(i in 1:n){
    Vi[[i]] <- Zi[[i]] %*% D[[r]] %*% t(Zi[[i]]) + sigma.sq[r]*diag(x = 1, nrow=ni[i], ncol = ni[i])
    if(q==1)
      inv.Vi[[r]][[i]] <- 
        (1/sigma.sq[r]) * (diag(rep(1,ni[i]))
                           -((as.numeric(D[[r]])/sigma.sq[r])/(1+ni[i]*(as.numeric(D[[r]])/sigma.sq[r])))*matrix(rep(1,(ni[i])^2)
                                                                                                                 , ncol=ni[i], nrow=ni[i]) )
    else inv.Vi[[r]][[i]] <- solve(Vi[[i]])
  }
  
  Vi <- list(NULL) 
  #inv.Vi[[r-1]] <- list(NULL) #not to run at step 0
  
  #the generalized log-likelihood (GLL) 
  GLL <- vector(mode="numeric") ; length(GLL) <- r
  term <- vector(mode="numeric",length=n)
  for(i in 1:n)
    term[i]<-t(epsili[[i]]) %*% solve(sigma.sq[r]*diag(x=1,nrow=ni[i],ncol=ni[i])) %*% epsili[[i]]
  + t(bi[[r]][[i]]) %*% solve(D[[r]]) %*% bi[[r]][[i]]
  + log(abs(D[[r]]))
  + log(abs(sigma.sq[r]*diag(x=1, nrow=ni[i], ncol = ni[i])))
  GLL[r] <- sum(term)
  rm(term)
  gc(verbose=FALSE)
  
  #convergence criterion
  Jump <- rep(NA,r) 		#at this first iteration Jump = NA
  convergence.iter<- rep(NA,r) 	#at this first convergence.iter = NA
  
  ####################
  ####	STEP 1	####        
  ####################
  
  #update iteration number r
  r <- r+1
  if (verbose) 	
    message("MERF iter no: ", r)
  
  #update the length of the different lists
  
  length(sigma.sq) <- r
  length(D) <- r
  length(inv.Vi) <- r
  length(bi) <- r
  length(GLL) <- r
  
  length(Jump) <- r
  length(convergence.iter) <- r
  
  #update the transformed outcome, star.Yi
  star.Yi <- list() 
  for(i in 1:n){
    star.Yi[[i]] <- Yi[[i]] - Zi[[i]] %*% bi[[r-1]][[i]]	
  }	
  
  #one STD random forest
  ######################
  
  MERF.lDB$star.Yi <- unlist(star.Yi) 
  rm(star.Yi) ;  gc(verbose=FALSE)
  
  
  fit.rf <- randomForest(
    formula=fit.rf.formula 
    ,data=MERF.lDB	
    ,ntree=ntree
    ,mtry = mtry 	
    ,replace=TRUE
    ,nodesize = nodesize
    ,proximity=FALSE	
  )
  
  #fixed part
  #as vector
  MERF.lDB$f.pred  <- predict(fit.rf, type="response")
  
  #in matrix format
  fixed.pred <- list()	
  fixed.pred <- split(MERF.lDB, grp) 	
  for(i in 1:n)fixed.pred[[i]] <- as.matrix(subset(fixed.pred[[i]] ,select=f.pred), ncol=1)
  
  #random	part
  ############
  for(i in 1:n)
    bi[[r]][[i]] <- D[[r-1]]%*%t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]] %*% (Yi[[i]] - fixed.pred[[i]])
  
  bi[r-1] <- list(NULL)		 
  
  ####################
  ####	STEP 2	####         
  ####################
  
  #level-1 variance component
  #residuals
  epsili <- list()
  for(i in 1:n)
    epsili[[i]] <- Yi[[i]] - fixed.pred[[i]] - Zi[[i]] %*% bi[[r]][[i]]
  
  
  term <- vector(mode="numeric",length=n)
  for(i in 1:n)
    term[i] <- crossprod(epsili[[i]]) + 
    sigma.sq[r-1] * (ni[i] - sigma.sq[r-1]* sum(diag(inv.Vi[[r-1]][[i]])))
  sigma.sq[r] <- (1/N)*(sum(term))
  rm(term) ;gc(verbose=FALSE)
  #message("sigmasq of current micro iter", sigma.sq[r] )
  
  #level-2 variance component
  term <- list()
  term[[1]] <- tcrossprod(bi[[r]][[1]]) + 
    (	D[[r-1]] - 
        D[[r-1]] %*% t(Zi[[1]])%*% inv.Vi[[r-1]][[1]] %*% Zi[[1]] %*% D[[r-1]]
    )
  for(i in 2:n) 
    term[[i]] <- term[[i-1]]+ tcrossprod(bi[[r]][[i]]) + 
    (	D[[r-1]] - 
        D[[r-1]] %*% t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]]%*% Zi[[i]]%*% D[[r-1]]
    )
  term <- term[[n]]
  D[[r]] <- (1/n)*term
  rm(term) ;gc(verbose=FALSE)	 
  #message("D of current micro iter: ", D[[r]] )
  
  #level-1 and level-2 variance components (or typical or total variance)
  inv.Vi[[r]] <-list()
  for(i in 1:n){
    Vi[[i]] <- Zi[[i]] %*% D[[r]] %*% t(Zi[[i]])+sigma.sq[r]*diag(x = 1, nrow=ni[i], ncol = ni[i])
    if(q==1)
      inv.Vi[[r]][[i]] <- 
        (1/sigma.sq[r]) * (diag(rep(1,ni[i]))
                           -((as.numeric(D[[r]])/sigma.sq[r])/(1+ni[i]*(as.numeric(D[[r]])/sigma.sq[r])))*matrix(rep(1,(ni[i])^2)
                                                                                                                 , ncol=ni[i], nrow=ni[i]) )
    else inv.Vi[[r]][[i]] <- solve(Vi[[i]])
  }
  Vi <- list(NULL) 
  inv.Vi[[r-1]] <- list(NULL) 	#not to run at step 0
  
  
  #the generalized log-likelihood (GLL) 
  term <- vector(mode="numeric",length=n)
  for(i in 1:n)
    term[i]<-t(epsili[[i]]) %*% solve(sigma.sq[r]*diag(x=1,nrow=ni[i],ncol=ni[i])) %*% epsili[[i]]
  + t(bi[[r]][[i]]) %*% solve(D[[r]]) %*% bi[[r]][[i]]
  + log(abs(D[[r]]))
  + log(abs(sigma.sq[r]*diag(x=1, nrow=ni[i], ncol = ni[i])))
  GLL[r] <- sum(term)
  rm(term)
  gc(verbose=FALSE)
  
  #update the value of the Jump in GLL
  Jump[r] <- abs( (GLL[r]- GLL[r-1])/GLL[r] )
  
  if(Jump[r] < smallest.Jump.allowed | Jump[r] == smallest.Jump.allowed) {
    convergence.iter[r] <- r
    if (verbose) message("Converg. at iter no: ", r)
  } 
  
  
  ####################
  ####	STEP 3	####         
  ####################
  
  
  ###################################################
  #Iterating "F.niter" times to avoid early stopping#
  ###################################################
  
  for(I in 3:F.niter){#repeat step 1 and 2
    
    ####################
    ####	STEP 1	####        
    ####################
    
    #update iteration number r
    r <- r+1
    if (verbose) 	
      message("MERF iter no: ", r)
    
    #update the length of the different lists
    
    length(sigma.sq) <- r
    length(D) <- r
    length(inv.Vi) <- r
    length(bi) <- r
    length(GLL) <- r
    
    length(Jump) <- r
    length(convergence.iter) <- r
    
    #update the transformed outcome, star.Yi
    star.Yi <- list() 
    for(i in 1:n){
      star.Yi[[i]] <- Yi[[i]] - Zi[[i]] %*% bi[[r-1]][[i]]	
    }
    
    #one STD random forest
    ######################
    
    MERF.lDB$star.Yi <- unlist(star.Yi) 
    rm(star.Yi) ;  gc(verbose=FALSE)
    
    fit.rf <- randomForest(
      formula=fit.rf.formula 
      ,data=MERF.lDB	
      ,ntree=ntree
      ,mtry = mtry 	
      ,replace=TRUE
      ,nodesize = nodesize
      ,proximity=FALSE	
    )
    
    #fixed part
    #as vector
    MERF.lDB$f.pred  <- predict(fit.rf, type="response")
    
    #in matrix format
    fixed.pred <- list()	
    fixed.pred <- split(MERF.lDB, grp) 	
    for(i in 1:n)fixed.pred[[i]] <- as.matrix(subset(fixed.pred[[i]] ,select=f.pred), ncol=1)
    
    #random	part
    ############
    for(i in 1:n)
      bi[[r]][[i]] <- D[[r-1]]%*%t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]] %*% (Yi[[i]] - fixed.pred[[i]])
    bi[r-1] <- list(NULL)		 
    
    
    ####################
    ####	STEP 2	####        
    ####################
    
    #level-1 variance component
    #residuals
    epsili <- list()
    for(i in 1:n)
      epsili[[i]] <- Yi[[i]] - fixed.pred[[i]] - Zi[[i]] %*% bi[[r]][[i]]
    
    
    term <- vector(mode="numeric",length=n)
    for(i in 1:n)
      term[i] <- crossprod(epsili[[i]]) + 
      sigma.sq[r-1] * (ni[i] - sigma.sq[r-1]* sum(diag(inv.Vi[[r-1]][[i]])))
    sigma.sq[r] <- (1/N)*(sum(term))
    rm(term) ;gc(verbose=FALSE)
    #message("sigmasq of current micro iter", sigma.sq[r] )
    
    #level-2 variance component
    term <- list()
    term[[1]] <- tcrossprod(bi[[r]][[1]]) + 
      (	D[[r-1]] - 
          D[[r-1]] %*% t(Zi[[1]])%*% inv.Vi[[r-1]][[1]] %*% Zi[[1]] %*% D[[r-1]]
      )
    for(i in 2:n) 
      term[[i]] <- term[[i-1]]+ tcrossprod(bi[[r]][[i]]) + 
      (	D[[r-1]] - 
          D[[r-1]] %*% t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]]%*% Zi[[i]]%*% D[[r-1]]
      )
    term <- term[[n]]
    D[[r]] <- (1/n)*term
    rm(term) ;gc(verbose=FALSE)	 
    #message("D of current micro iter: ", D[[r]] )
    
    #level-1 and level-2 variance components (or typical or total variance)
    inv.Vi[[r]] <-list()
    for(i in 1:n){
      Vi[[i]] <- Zi[[i]] %*% D[[r]] %*% t(Zi[[i]])+sigma.sq[r]*diag(x = 1, nrow=ni[i], ncol = ni[i])
      if(q==1)
        inv.Vi[[r]][[i]] <- 
          (1/sigma.sq[r]) * (diag(rep(1,ni[i]))
                             -((as.numeric(D[[r]])/sigma.sq[r])/(1+ni[i]*(as.numeric(D[[r]])/sigma.sq[r])))*matrix(rep(1,(ni[i])^2)
                                                                                                                   , ncol=ni[i], nrow=ni[i]) )
      else inv.Vi[[r]][[i]] <- solve(Vi[[i]])
    }
    Vi <- list(NULL) 
    inv.Vi[[r-1]] <- list(NULL) 	#not to run at step 0
    
    
    #the generalized log-likelihood (GLL) 
    term <- vector(mode="numeric",length=n)
    for(i in 1:n)
      term[i]<-t(epsili[[i]]) %*% solve(sigma.sq[r]*diag(x=1,nrow=ni[i],ncol=ni[i])) %*% epsili[[i]]
    + t(bi[[r]][[i]]) %*% solve(D[[r]]) %*% bi[[r]][[i]]
    + log(abs(D[[r]]))
    + log(abs(sigma.sq[r]*diag(x=1, nrow=ni[i], ncol = ni[i])))
    GLL[r] <- sum(term)
    rm(term)
    gc(verbose=FALSE)
    
    #update the value of the Jump in GLL
    
    Jump[r] <- abs( (GLL[r]- GLL[r-1])/GLL[r] )
    
    if(Jump[r] < smallest.Jump.allowed | Jump[r] == smallest.Jump.allowed) {
      convergence.iter[r] <- r
      if (verbose) message("Converg. at iter no: ", r)
    } 
    
  }
  #end for (I in 1: F.niter)
  ###########################
  
  ######################################
  #Iterating "max.niter" times at most #
  ######################################
  
  while( r < (F.niter + max.niter) ){
    
    if(Jump[r] > smallest.Jump.allowed){ #repeat step 1 and 2
      
      ####################
      ####	STEP 1	####        
      ####################
      
      #update iteration number r
      r <- r+1
      if (verbose) 	
        message("MERF iter no: ", r)
      
      #update the length of the different lists
      
      length(sigma.sq) <- r
      length(D) <- r
      length(inv.Vi) <- r
      length(bi) <- r
      length(GLL) <- r
      
      length(Jump) <- r
      length(convergence.iter) <- r
      
      #update the transformed outcome, star.Yi
      star.Yi <- list() 
      for(i in 1:n){
        star.Yi[[i]] <- Yi[[i]] - Zi[[i]] %*% bi[[r-1]][[i]]	
      }
      
      #one STD random forest
      ######################
      
      MERF.lDB$star.Yi <- unlist(star.Yi) 
      rm(star.Yi) ;  gc(verbose=FALSE)
      
      fit.rf <- randomForest(
        formula=fit.rf.formula 
        ,data=MERF.lDB	
        ,ntree=ntree
        ,mtry = mtry 	
        ,replace=TRUE
        ,nodesize = nodesize
        ,proximity=FALSE	
      )
      
      #fixed part
      #as vector
      MERF.lDB$f.pred  <- predict(fit.rf, type="response")
      
      #in matrix format
      fixed.pred <- list()	
      fixed.pred <- split(MERF.lDB, grp) 	
      for(i in 1:n)fixed.pred[[i]] <- as.matrix(subset(fixed.pred[[i]] ,select=f.pred), ncol=1)
      
      #random	part
      ############
      for(i in 1:n)
        bi[[r]][[i]] <- D[[r-1]]%*%t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]] %*% (Yi[[i]] - fixed.pred[[i]])
      bi[r-1] <- list(NULL)		
      
      
      ####################
      ####	STEP 2	####        
      ####################
      
      #level-1 variance component
      #residuals
      epsili <- list()
      for(i in 1:n)
        epsili[[i]] <- Yi[[i]] - fixed.pred[[i]] - Zi[[i]] %*% bi[[r]][[i]]
      
      
      term <- vector(mode="numeric",length=n)
      for(i in 1:n)
        term[i] <- crossprod(epsili[[i]]) + 
        sigma.sq[r-1] * (ni[i] - sigma.sq[r-1]* sum(diag(inv.Vi[[r-1]][[i]])))
      sigma.sq[r] <- (1/N)*(sum(term))
      rm(term) ;gc(verbose=FALSE)
      #message("sigmasq of current micro iter", sigma.sq[r] )
      
      #level-2 variance component
      term <- list()
      term[[1]] <- tcrossprod(bi[[r]][[1]]) + 
        (	D[[r-1]] - 
            D[[r-1]] %*% t(Zi[[1]])%*% inv.Vi[[r-1]][[1]] %*% Zi[[1]] %*% D[[r-1]]
        )
      for(i in 2:n) 
        term[[i]] <- term[[i-1]]+ tcrossprod(bi[[r]][[i]]) + 
        (	D[[r-1]] - 
            D[[r-1]] %*% t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]]%*% Zi[[i]]%*% D[[r-1]]
        )
      term <- term[[n]]
      D[[r]] <- (1/n)*term
      rm(term) ;gc(verbose=FALSE)	 
      #message("D of current micro iter: ", D[[r]] )
      
      #level-1 and level-2 variance components (or typical or total variance)
      inv.Vi[[r]] <-list()
      for(i in 1:n){
        Vi[[i]] <- Zi[[i]] %*% D[[r]] %*% t(Zi[[i]])+sigma.sq[r]*diag(x = 1, nrow=ni[i], ncol = ni[i])
        if(q==1)
          inv.Vi[[r]][[i]] <- 
            (1/sigma.sq[r]) * (diag(rep(1,ni[i]))
                               -((as.numeric(D[[r]])/sigma.sq[r])/(1+ni[i]*(as.numeric(D[[r]])/sigma.sq[r])))*matrix(rep(1,(ni[i])^2)
                                                                                                                     , ncol=ni[i], nrow=ni[i]) )
        else inv.Vi[[r]][[i]] <- solve(Vi[[i]])
      }
      Vi <- list(NULL) 
      inv.Vi[[r-1]] <- list(NULL) 	#not to run at step 0
      
      
      #the generalized log-likelihood (GLL) 
      term <- vector(mode="numeric",length=n)
      for(i in 1:n)
        term[i]<-t(epsili[[i]]) %*% solve(sigma.sq[r]*diag(x=1,nrow=ni[i],ncol=ni[i])) %*% epsili[[i]]
      + t(bi[[r]][[i]]) %*% solve(D[[r]]) %*% bi[[r]][[i]]
      + log(abs(D[[r]]))
      + log(abs(sigma.sq[r]*diag(x=1, nrow=ni[i], ncol = ni[i])))
      GLL[r] <- sum(term)
      rm(term)
      gc(verbose=FALSE)
      
      #update the value of the Jump in GLL
      
      Jump[r] <- abs( (GLL[r]- GLL[r-1])/GLL[r] )
      
      if(Jump[r] < smallest.Jump.allowed | Jump[r] == smallest.Jump.allowed) {
        convergence.iter[r] <- r
        if (verbose) message("Converg. at iter no: ", r)
      } 
      
      
    }
    #end if(Jump[r] > smallest.Jump.allowed) and STOP repeating step 1 and 2
    
    else break 
    #end while( r < (2 + F.niter + max.niter) )
    
    
    ############################
  }###	END OF STEP 3	####        
  ############################
  
  
  #output to be returned (MERF model is the one at the last iteration)
  ###################### 
  
  
  output <- list(
    Jump[r]
    ,GLL
    ,convergence.iter
    ,fit.rf 
    ,bi[[r]]
    ,sigma.sq[r]
    ,D	#,D[[r]]
    ,MERF.lDB
  )
  
  names(output) <- c(
    "Jump[r]"
    ,"GLL"
    ,"convergence.iter"
    ,"fit.rf"
    ,"bi[[r]]"
    ,"sigma.sq[r]"
    ,"D"	#,"D[[r]]"
    ,"fit.data"
  )
  
  #clean memory
  #############
  rm(
    xnam
    
    ,MERF.lDB	
    
    ,ni,n,N
    ,Zi,q
    ,Yi	
    
    ,ntree
    ,mtry
    ,nodesize
    ,fit.rf.formula 
    ,fit.rf
    ,fixed.pred ,epsili
    
    ,sigma.sq,D,bi,Vi,inv.Vi
    
    ,F.niter ,max.niter
    ,smallest.Jump.allowed ,GLL ,Jump ,convergence.iter
    ,r,i,I
    
    ,verbose 		
  )
  gc(verbose=FALSE)
  
  
  #return
  #######
  output
  
}#end of MERF function
############################################################################################################
############################################################################################################
############################################################################################################
### end of the script	###
#sink()






MERF.short  <- function(
  
  xnam
  
  ,MERF.lDB 	
  
  ,ni
  ,Zi
  ,Yi
  ,grp
  
  ,ntree
  ,mtry
  ,nodesize
  
  ,sigmasqzero = NULL
  ,Dzero = NULL
  ,bizero = NULL
  
  ,F.niter
  ,max.niter
  ,smallest.Jump.allowed
  
  ,verbose = TRUE
  
){
  
  ####################
  ####	STEP 0	####
  #################### 
  
  #Memory Allocation and initialization:
  
  #Parameters values 
  n <- length(ni)
  N <- sum(ni)
  q <- dim(Zi[[1]])[2] 	# q=1 in random intercept case
  
  
  #Initial values of sigmasqzero, Dzero, and bizero
  if( is.null(sigmasqzero) ) sigmasqzero <- 1
  else sigmasqzero <- sigmasqzero
  
  if( is.null(Dzero) ){
    Dzero <- diag(1 ,nrow=q,ncol=q)
  }
  else Dzero <- Dzero
  
  if( is.null(bizero) ){
    bizero <- list(); length(bizero)<- n
    for(i in 1:n) bizero[[i]] <- matrix(0,nrow=q,ncol=1)
  }	
  else bizero <- bizero
  
  #iter number
  r <- 1
  if (verbose) 	
    message("MERF iter no: ", r)
  
  #transformed outcome, star.Yi[[r]][[i]], initialized with the original values
  star.Yi <- list() 
  for(i in 1:n){
    star.Yi[[i]] <- Yi[[i]] - Zi[[i]] %*% bizero[[i]]
  }
  
  #one STD random forest
  ######################
  
  MERF.lDB$star.Yi <- unlist(star.Yi) 
  rm(star.Yi) ;  gc(verbose=FALSE)
  
  fit.rf.formula <-  as.formula(paste("star.Yi ~ ", paste(xnam, collapse= "+")))
  
  fit.rf <- randomForest(
    formula=fit.rf.formula 
    ,data=MERF.lDB	
    ,ntree=ntree
    ,mtry = mtry 	
    ,replace=TRUE	
    ,nodesize = nodesize
    ,proximity=FALSE	
  )
  
  #fixed part
  #as vector			
  MERF.lDB$f.pred  <- predict(fit.rf, type="response")			#!!! use the out-of-bag predictions
  #in matrix format
  fixed.pred <- list()	
  fixed.pred <- split(MERF.lDB, grp) 	
  for(i in 1:n)fixed.pred[[i]] <- as.matrix(subset(fixed.pred[[i]] ,select=f.pred), ncol=1)
  
  #random	part
  ############
  #random effects parameters in list format
  bi <- list(list()) ; length(bi) <- r
  for(i in 1:n)bi[[r]][[i]] <- bizero[[i]] 	
  #print("bizero");print(bizero)
  rm(bizero) ; gc(verbose=FALSE)
  
  #level-1 variance component
  #residuals
  epsili <- list()
  for(i in 1:n)
    epsili[[i]] <- Yi[[i]] - fixed.pred[[i]] - Zi[[i]] %*% bi[[r]][[i]]
  
  sigma.sq <- vector(mode="numeric") ;length(sigma.sq) <- r
  sigma.sq[r] <- sigmasqzero
  #print("sigmasqzero") ;print(sigmasqzero)			
  rm(sigmasqzero) ; gc(verbose=FALSE)
  #message("sigmasq of current micro iter", sigma.sq[r] )
  
  #level-2 variance component
  D <- list() ;length(D) <- r
  D[[r]] <- Dzero		#!!!Dzero <- diag(x=0.01, nrow=q, ncol = q)
  #print("Dzero") ;print(Dzero)	
  rm(Dzero) ; gc(verbose=FALSE)
  #message("D of current micro iter: ", D[[r]] )
  
  #level-1 and level-2 variance components (or typical or total variance)
  Vi <- list() 
  inv.Vi <- list(list()) ; length(inv.Vi) <- r
  
  for(i in 1:n){
    Vi[[i]] <- Zi[[i]] %*% D[[r]] %*% t(Zi[[i]]) + sigma.sq[r]*diag(x = 1, nrow=ni[i], ncol = ni[i])
    if(q==1)
      inv.Vi[[r]][[i]] <- 
        (1/sigma.sq[r]) * (diag(rep(1,ni[i]))
                           -((as.numeric(D[[r]])/sigma.sq[r])/(1+ni[i]*(as.numeric(D[[r]])/sigma.sq[r])))*matrix(rep(1,(ni[i])^2)
                                                                                                                 , ncol=ni[i], nrow=ni[i]) )
    else inv.Vi[[r]][[i]] <- solve(Vi[[i]])
  }
  
  Vi <- list(NULL) 
  #inv.Vi[[r-1]] <- list(NULL) #not to run at step 0
  
  #the generalized log-likelihood (GLL) 
  GLL <- vector(mode="numeric") ; length(GLL) <- r
  term <- vector(mode="numeric",length=n)
  for(i in 1:n)
    term[i]<-t(epsili[[i]]) %*% solve(sigma.sq[r]*diag(x=1,nrow=ni[i],ncol=ni[i])) %*% epsili[[i]]
  + t(bi[[r]][[i]]) %*% solve(D[[r]]) %*% bi[[r]][[i]]
  + log(abs(D[[r]]))
  + log(abs(sigma.sq[r]*diag(x=1, nrow=ni[i], ncol = ni[i])))
  GLL[r] <- sum(term)
  rm(term)
  gc(verbose=FALSE)
  
  #convergence criterion
  Jump <- rep(NA,r) 		#at this first iteration Jump = NA
  convergence.iter<- rep(NA,r) 	#at this first convergence.iter = NA
  
  ####################
  ####	STEP 1	####        
  ####################
  
  #update iteration number r
  r <- r+1
  if (verbose) 	
    message("MERF iter no: ", r)
  
  #update the length of the different lists
  
  length(sigma.sq) <- r
  length(D) <- r
  length(inv.Vi) <- r
  length(bi) <- r
  length(GLL) <- r
  
  length(Jump) <- r
  length(convergence.iter) <- r
  
  #update the transformed outcome, star.Yi
  star.Yi <- list() 
  for(i in 1:n){
    star.Yi[[i]] <- Yi[[i]] - Zi[[i]] %*% bi[[r-1]][[i]]	
  }	
  
  #one STD random forest
  ######################
  
  MERF.lDB$star.Yi <- unlist(star.Yi) 
  rm(star.Yi) ;  gc(verbose=FALSE)
  
  
  fit.rf <- randomForest(
    formula=fit.rf.formula 
    ,data=MERF.lDB	
    ,ntree=ntree
    ,mtry = mtry 	
    ,replace=TRUE
    ,nodesize = nodesize
    ,proximity=FALSE	
  )
  
  #fixed part
  #as vector
  MERF.lDB$f.pred  <- predict(fit.rf, type="response")
  
  #in matrix format
  fixed.pred <- list()	
  fixed.pred <- split(MERF.lDB, grp) 	
  for(i in 1:n)fixed.pred[[i]] <- as.matrix(subset(fixed.pred[[i]] ,select=f.pred), ncol=1)
  
  #random	part
  ############
  for(i in 1:n)
    bi[[r]][[i]] <- D[[r-1]]%*%t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]] %*% (Yi[[i]] - fixed.pred[[i]])
  
  bi[r-1] <- list(NULL)		 
  
  ####################
  ####	STEP 2	####         
  ####################
  
  #level-1 variance component
  #residuals
  epsili <- list()
  for(i in 1:n)
    epsili[[i]] <- Yi[[i]] - fixed.pred[[i]] - Zi[[i]] %*% bi[[r]][[i]]
  
  
  term <- vector(mode="numeric",length=n)
  for(i in 1:n)
    term[i] <- crossprod(epsili[[i]]) + 
    sigma.sq[r-1] * (ni[i] - sigma.sq[r-1]* sum(diag(inv.Vi[[r-1]][[i]])))
  sigma.sq[r] <- (1/N)*(sum(term))
  rm(term) ;gc(verbose=FALSE)
  #message("sigmasq of current micro iter", sigma.sq[r] )
  
  #level-2 variance component
  term <- list()
  term[[1]] <- tcrossprod(bi[[r]][[1]]) + 
    (	D[[r-1]] - 
        D[[r-1]] %*% t(Zi[[1]])%*% inv.Vi[[r-1]][[1]] %*% Zi[[1]] %*% D[[r-1]]
    )
  for(i in 2:n) 
    term[[i]] <- term[[i-1]]+ tcrossprod(bi[[r]][[i]]) + 
    (	D[[r-1]] - 
        D[[r-1]] %*% t(Zi[[i]]) %*% inv.Vi[[r-1]][[i]]%*% Zi[[i]]%*% D[[r-1]]
    )
  term <- term[[n]]
  D[[r]] <- (1/n)*term
  rm(term) ;gc(verbose=FALSE)	 
  #message("D of current micro iter: ", D[[r]] )
  
  #level-1 and level-2 variance components (or typical or total variance)
  inv.Vi[[r]] <-list()
  for(i in 1:n){
    Vi[[i]] <- Zi[[i]] %*% D[[r]] %*% t(Zi[[i]])+sigma.sq[r]*diag(x = 1, nrow=ni[i], ncol = ni[i])
    if(q==1)
      inv.Vi[[r]][[i]] <- 
        (1/sigma.sq[r]) * (diag(rep(1,ni[i]))
                           -((as.numeric(D[[r]])/sigma.sq[r])/(1+ni[i]*(as.numeric(D[[r]])/sigma.sq[r])))*matrix(rep(1,(ni[i])^2)
                                                                                                                 , ncol=ni[i], nrow=ni[i]) )
    else inv.Vi[[r]][[i]] <- solve(Vi[[i]])
  }
  Vi <- list(NULL) 
  inv.Vi[[r-1]] <- list(NULL) 	#not to run at step 0
  
  
  #the generalized log-likelihood (GLL) 
  term <- vector(mode="numeric",length=n)
  for(i in 1:n)
    term[i]<-t(epsili[[i]]) %*% solve(sigma.sq[r]*diag(x=1,nrow=ni[i],ncol=ni[i])) %*% epsili[[i]]
  + t(bi[[r]][[i]]) %*% solve(D[[r]]) %*% bi[[r]][[i]]
  + log(abs(D[[r]]))
  + log(abs(sigma.sq[r]*diag(x=1, nrow=ni[i], ncol = ni[i])))
  GLL[r] <- sum(term)
  rm(term)
  gc(verbose=FALSE)
  
  #update the value of the Jump in GLL
  Jump[r] <- abs( (GLL[r]- GLL[r-1])/GLL[r] )
  
  if(Jump[r] < smallest.Jump.allowed | Jump[r] == smallest.Jump.allowed) {
    convergence.iter[r] <- r
    if (verbose) message("Converg. at iter no: ", r)
  } 
  
  
  #output to be returned (MERF model is the one at the last iteration)
  ###################### 
  
  
  output <- list(
    Jump[r]
    ,GLL
    ,convergence.iter
    ,fit.rf 
    ,bi[[r]]
    ,sigma.sq[r]
    ,D	#,D[[r]]
    ,MERF.lDB
  )
  
  names(output) <- c(
    "Jump[r]"
    ,"GLL"
    ,"convergence.iter"
    ,"fit.rf"
    ,"bi[[r]]"
    ,"sigma.sq[r]"
    ,"D"	#,"D[[r]]"
    ,"fit.data"
  )
  
  #clean memory
  #############
  rm(
    xnam
    
    ,MERF.lDB	
    
    ,ni,n,N
    ,Zi,q
    ,Yi	
    
    ,ntree
    ,mtry
    ,nodesize
    ,fit.rf.formula 
    ,fit.rf
    ,fixed.pred ,epsili
    
    ,sigma.sq,D,bi,Vi,inv.Vi
    
    ,F.niter ,max.niter
    ,smallest.Jump.allowed ,GLL ,Jump ,convergence.iter
    ,r,i,I
    
    ,verbose 		
  )
  gc(verbose=FALSE)
  
  
  #return
  #######
  output
  
}#end of MERF function
############################################################################################################
############################################################################################################
############################################################################################################
### end of the script	###
#sink()












MERF.fit.short <- function(data, mtry = 4, nodesize = 5, crit=1e-4, nr.group, group.size, nrX.within, nrX.between,  max.niter=15, d0){
  xnam = names(data)[1:(nrX.between+nrX.within)] # the predictor variables needed for MERF fit
  ni = c(rep(group.size, nr.group))
  
  ### Yi ###
  Yi = split(data$y, data$grp)
  #str(Yi)
  
  ### Zi ###
  ## if q = 1 (i.e. random intercept)
  Zi.ri = list(); length(Zi.ri)= length(ni)
  for(i in 1:length(Zi.ri))    
    Zi.ri[[i]] <- matrix(rep(1,ni[i]),nrow=ni[i],ncol=1)             	         
  Zi = Zi.ri
  
  #m <- max(floor((nrX.between + nrX.within)/3), 1)
  
  ### MERF application ### ----
  MERF.fit <- MERF.short(
    xnam=xnam
    ,MERF.lDB = data			
    ,ni = ni
    ,Zi = Zi
    ,Yi = Yi
    ,grp = data$grp
    
    ,ntree= 500
    ,mtry= mtry		
    ,nodesize= nodesize
    
    ,sigmasqzero = NULL		
    ,Dzero =  Dzero <- diag(d0 ,nrow=1,ncol=1)		
    ,bizero = NULL
    
    ,F.niter = 0	
    ,max.niter = max.niter
    ,smallest.Jump.allowed = crit
    
    ,verbose= TRUE
  )
  return(MERF.fit) # has output: -- fited RF -- sigma_squared (resid) -- sigma_INT (Random effect)
}

