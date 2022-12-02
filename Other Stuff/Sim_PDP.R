rm(list=ls())
library(mvtnorm)
library(matrixcalc)
library(ggplot2)
library(ggpubr)
library(cowplot)
library(mlr)
library(iml)
library(MixMatrix)
################################################################################
################################################################################
########################     Simulate Once      ################################
################################################################################
################################################################################
################################################################################

simulate_once <- function(obs){
  
  #### Simulate data once from normal distribution
  sigma <- CSgenerate(3, 0)
  data_once <- as.data.frame(rmvnorm(n = obs, 
                                     mean = rep(0, times = 3), 
                                     sigma = sigma)) # variance=1 -> cov(X,Y)=corr(X,Y)
  colnames(data_once) <- c("X1", "X2", "X3")
  
  ### Dataset independent features ###
  data_once_ind <- data_once
  std <- abs(mean(sin(3*data_once_ind$X1) + data_once_ind$X2 + data_once_ind$X3))*0.1
  data_once_ind$Y <- sin(3*data_once_ind$X1) + data_once_ind$X2 + data_once_ind$X3 + rnorm(n = obs, mean = 0, sd = std)
 
  ### Dataset correlated features ###
  sigma <- CSgenerate(3, 0.9)
  if(is.positive.semi.definite(sigma) == FALSE) 
    warning("Covariance Matrix is not positive semidefinite!")
  data_once_cor <- as.data.frame(rmvnorm(n = obs, 
                                         mean = rep(0, times = 3), 
                                         sigma = sigma))
  colnames(data_once_cor) <- c("X1", "X2", "X3")
  std <- abs(mean(sin(3*data_once_cor$X1) + data_once_cor$X2 + data_once_cor$X3))*0.1
  data_once_cor$Y <- sin(3*data_once_cor$X1) + data_once_cor$X2 + data_once_cor$X3 + rnorm(n = obs, mean = 0, sd = std)
  list(data_once_ind, data_once_cor)
}


################################################################################
################################################################################
######################     Repeat Simulation      ##############################
################################################################################
################################################################################
################################################################################
simulate_repeated <- function(no_sim, obs){
  replicate_list <- replicate(n = no_sim, simulate_once(obs=obs))
  # Rows: Independent, Correlated
  # Columns: Predictors
  

  # Get all Datasets with independent features
  data_replicated_ind <- vector(mode = 'list', length = no_sim)
  data_replicated_cor <- vector(mode = 'list', length = no_sim)
  l1 <- seq(from = 1, by=2, length = no_sim)
  cnt <- 1
  for (j in l1){
    data_replicated_ind[[cnt]] <- replicate_list[[j]]
    data_replicated_cor[[cnt]] <- replicate_list[[j+1]]
    cnt <- cnt + 1
  }  
  
  list(data_replicated_ind, data_replicated_cor)
}


################################################################################
################################################################################
######################     Get Simulated Data     ##############################
################################################################################
################################################################################
################################################################################
obs <- 500
no_sim <- 4
lrn = makeLearner("regr.randomForest")
grid_size=30
simulation <- simulate_repeated(no_sim, obs)
independent <- simulation[[1]]
correlated <- simulation[[2]]




################################################################################
################################################################################
########################     Make Predictions     ##############################
################################################################################
################################################################################
################################################################################
prediction <- function(sim_datasets){
  
  ### Derive prediction for each simulated sample ###
  multi_sim_X1 <- data.frame()
  multi_sim_X2 <- data.frame()
  multi_sim_X3 <- data.frame()
  
  sim_datasets <- independent
  j <- 1
  
  f <- function(data){
    t(sapply(data, range))
  }
  
  tmp <- lapply(sim_datasets, f)
  range_mat <- matrix(NA, nrow = 3, ncol = 2)
  range_mat[,1] <- do.call('pmin', tmp)[1:3,1]
  range_mat[,2] <- do.call('pmax', tmp)[1:3,2]
  
  for (j in 1:no_sim){
    task <- makeRegrTask(data=sim_datasets[[j]], target="Y")
    trained <- train(learner = lrn, task = task)
    pred <- Predictor$new(trained, data = sim_datasets[[j]])
    mod1 <- FeatureEffect$new(pred, feature = "X1", method = "pdp", grid.size = grid_size)
    multi_sim_X1[((j-1)*grid_size+1):(j*grid_size),1] <- mod1$results[,1]
    multi_sim_X1[((j-1)*grid_size+1):(j*grid_size),2] <- mod1$results[,2]
    multi_sim_X1$Simulation[((j-1)*grid_size+1):(j*grid_size)] <- j 
    mod2 <- FeatureEffect$new(pred, feature = "X2", method = "pdp", grid.size = grid_size)
    multi_sim_X2[((j-1)*grid_size+1):(j*grid_size),1] <- mod2$results[,1]
    multi_sim_X2[((j-1)*grid_size+1):(j*grid_size),2] <- mod2$results[,2]
    multi_sim_X2$Simulation[((j-1)*grid_size+1):(j*grid_size)] <- j 
    
    mod3 <- FeatureEffect$new(pred, feature = "X3", method = "pdp", grid.size = grid_size)
    multi_sim_X3[((j-1)*grid_size+1):(j*grid_size),1] <- mod3$results[,1]
    multi_sim_X3[((j-1)*grid_size+1):(j*grid_size),2] <- mod3$results[,2]
    multi_sim_X3$Simulation[((j-1)*grid_size+1):(j*grid_size)] <- j }
  
  colnames(multi_sim_X1) <- c("X1", "Y_hat", "Simulation")
  colnames(multi_sim_X2) <- c("X2", "Y_hat", "Simulation")
  colnames(multi_sim_X3) <- c("X3", "Y_hat", "Simulation")
  multi_sim_X1$Simulation <- factor(multi_sim_X1$Simulation)
  multi_sim_X2$Simulation <- factor(multi_sim_X2$Simulation)
  multi_sim_X3$Simulation <- factor(multi_sim_X3$Simulation)
  
  list(multi_sim_X1, multi_sim_X2, multi_sim_X3)
}



ind <- prediction(independent)

#### Function to add an "overall grid" to produce mean PDP and error bars
add_avg_grid <- function(input, element){
  pdp <- as.data.frame(input[[element]])
  pdp$Grid <- rep(0, nrow(pdp))
  for (j in 1:grid_size){
    for (i in 1:nrow(pdp)){
      abst <- (max(pdp[,1])-min(pdp[,1]))/grid_size
      if (min(pdp[,1])+(j-1)*abst <= pdp[,1][i] & pdp[,1][i] < min(pdp[,1])+j*abst) {
        pdp$Grid[i] <- min(pdp[,1])+(j-0.5)*abst
        pdp$Grid[pdp[,1]==max(pdp[,1])] <- max(pdp[,1])-0.5*abst
      }
    }
  }
  pdp
}

pdp_ind_X1 <- add_avg_grid(input=ind, element=1)
pdp_ind_X2 <- add_avg_grid(input=ind, element=2)
pdp_ind_X3 <- add_avg_grid(input=ind, element=3)


#pdp_dep_X1 <- add_avg_grid(input=dep, element=1)
#pdp_dep_X2 <- add_avg_grid(input=dep, element=2)
#pdp_dep_X3 <- add_avg_grid(input=dep, element=3)

#### Function to derive mean and error bars per grid interval
add_mean_sd <- function(data){
  average <- data %>% 
    group_by(Grid) %>% 
    summarise(mean = mean(Y_hat), 
              sd = sd(Y_hat))
  average
}

avg_ind_X1 <- add_mean_sd(pdp_ind_X1)
avg_ind_X2 <- add_mean_sd(pdp_ind_X2)
avg_ind_X3 <- add_mean_sd(pdp_ind_X3)


#avg_dep_X1 <- add_mean_sd(pdp_dep_X1)
#avg_dep_X2 <- add_mean_sd(pdp_dep_X2)
#avg_dep_X3 <- add_mean_sd(pdp_dep_X3)


#### Function to produce PDP
pdp_plot <- function(pdp_data, feature){
  ggplot(data.frame(x = c(-3, 3)), aes(x = x))+
    geom_line(data = pdp_data, aes(x = pdp_data[,1], y = Y_hat, group = Simulation), alpha = 0.3)+
    theme_bw()+
    theme(plot.title = element_text(size=16, hjust = 0))+
    xlab(feature)+
    labs(title = "")+
    ylab("")+
    xlim(c(-3,3))+
    ylim(c(-3,3))
}

# Plots Independent Case
p11 <- pdp_plot(pdp_ind_X1, feature = "X1") +
  stat_function(fun = function(x) sin(3*x), geom = "line", col="red", size=1)
p12 <- pdp_plot(pdp_ind_X2, feature = "X2") + geom_abline(slope = 1, col = "red")
p13 <- pdp_plot(pdp_ind_X3, feature = "X3") + geom_abline(slope = 1, col = "red")

ga_ind <- ggarrange(p11+labs(title = "PDPs for independent Features"), 
                    p12, 
                    p13,
                    nrow = 1, ncol = 3)



list(ga_ind)
