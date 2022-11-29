rm(list=ls())
library(mvtnorm)
library(matrixcalc)
library(tidyverse)
library(stringr)
library(mlr)
library(iml)
library(ggplot2)
library(ggpubr)
## Setting 1: Linear Dependence of all Variables ## 
#### PDP - 1 Feature

simulate_once <- function(obs){
    
    #### Simulate data once from normal distribution:
    data_once <- as.data.frame(rmvnorm(n = obs, 
                                       mean = rep(0, times = 3), 
                                       sigma = diag(1, nrow = 3))) # variance=1 -> cov(X,Y)=corr(X,Y)
    colnames(data_once) <- c("X1", "X2", "X3")
    
    ### Dataset independent features ###
    data_once_ind <- data_once
    std <- abs(mean(data_once_ind$X1 + data_once_ind$X2 + data_once_ind$X3))*0.1
    data_once_ind$Y <- sin(3*data_once_ind$X1) + data_once_ind$X2 + data_once_ind$X3 + rnorm(n = obs, mean = 0, sd = std)
    ### Dataset dependent features ###
    data_once_dep <- data_once
    data_once_dep$X2 <- data_once_dep$X1
    std <- abs(mean(data_once_dep$X1 + data_once_dep$X2 + data_once_dep$X3))*0.1
    data_once_dep$Y <- data_once_dep$X1 + data_once_dep$X2 + data_once_dep$X3 + rnorm(n = obs, mean = 0, sd = std)
    ### Dataset correlated features ###
    sigma <- diag(1, nrow = 3)
    sigma[1,2] <- sigma[2,1] <- 0.9
    if(is.positive.semi.definite(sigma) == FALSE) 
      warning("Covariance Matrix is not positive semidefinite!")
    data_once_cor <- as.data.frame(rmvnorm(n = obs, 
                                           mean = rep(0, times = 3), 
                                           sigma = sigma))
    colnames(data_once_cor) <- c("X1", "X2", "X3")
    std <- abs(mean(data_once_cor$X1 + data_once_cor$X2 + data_once_cor$X3))*0.1
    data_once_cor$Y <- data_once_cor$X1 + data_once_cor$X2 + data_once_cor$X3 + rnorm(n = obs, mean = 0, sd = std)
    list(data_once_ind, data_once_dep, data_once_cor)
}
  
simulate_repeated <- function(no_sim, obs){
    replicate_list <- replicate(n = no_sim, simulate_once(obs=obs))
    
    # Dataset independent features - repeated simulation
    data_replicated_ind <- data.frame(c(rep(0, obs)))
    l1 <- seq(from = 1, by=3, length = no_sim)
    for (j in l1){
      for (i in 1:4){
        data_replicated_ind[,ncol(data_replicated_ind)+1] <- c(replicate_list[[j]][i])}}
    data_replicated_ind <- data_replicated_ind[, -1]
    colnames(data_replicated_ind) <- c(paste0(rep(c("X1", "X2", "X3", "Y"),times=no_sim), "_", rep(c(1:no_sim), each = 4)))
    
    # Dataset dependent features - repeated simulation
    data_replicated_dep <- data.frame(c(rep(0, obs)))
    l2 <- seq(from = 2, by=3, length = no_sim)
    for (k in l2){
      for (l in 1:4){
        data_replicated_dep[,ncol(data_replicated_dep)+1] <- c(replicate_list[[k]][l])}}
    data_replicated_dep <- data_replicated_dep[, -1]
    colnames(data_replicated_dep) <- c(paste0(rep(c("X1", "X2", "X3", "Y"),times=no_sim), "_", rep(c(1:no_sim), each = 4)))
    
    # Dataset correlated features - repeated simulation
    data_replicated_cor <- data.frame(c(rep(0, obs)))
    l3 <- seq(from = 3, by=3, length = no_sim)
    for (k in l3){
      for (l in 1:4){
        data_replicated_cor[,ncol(data_replicated_cor)+1] <- c(replicate_list[[k]][l])}}
    data_replicated_cor <- data_replicated_cor[, -1]
    colnames(data_replicated_cor) <- c(paste0(rep(c("X1", "X2", "X3", "Y"),times=no_sim), "_", rep(c(1:no_sim), each = 4)))
    
    list(data_replicated_ind, data_replicated_dep, data_replicated_cor)
}
  
no_sim=20
obs=500
lrn = makeLearner("regr.randomForest")
grid_size=30
independent <- simulate_repeated(no_sim, obs)[[1]]
dependent <- simulate_repeated(no_sim, obs)[[2]]
correlated <- simulate_repeated(no_sim, obs)[[3]] 
dataset <- independent
#### Function to produce predictions based on repeated simulations
prediction <- function(dataset){
    
    # Convert simulated datasets to long-format
    data <- dataset %>% 
      gather("key", "Value") %>% 
      mutate(Variable = str_sub(key, 1, 2),
             Simulation = str_sub(key, -2, -1)) %>% 
      select(Value:Simulation)
    data$Variable = sub("_$", "", data$Variable)
    data$Simulation = sub("_", "", data$Simulation)
    
    # Prepare Data for Prediction
    set <- data.frame()
    set[1:obs, 1:4] <- data[data$Simulation == 1,] %>% 
      mutate(Index = rep(1:obs, 4)) %>% 
      spread(key = Variable, value = Value) %>% 
      select(X1:Y)
    for (i in 1:no_sim){
      set[((i-1)*obs+1):(i*obs), 1:4] <- data[data$Simulation == i,] %>% 
        mutate(Index = rep(1:obs, 4)) %>% 
        spread(key = Variable, value = Value) %>% 
        select(X1:Y)}
    
    ### Derive prediction for each simulated sample ###
    multi_sim_X1 <- data.frame()
    multi_sim_X2 <- data.frame()
    multi_sim_X3 <- data.frame()
    
    for (j in 1:no_sim){
      task <- makeRegrTask(data=set[((j-1)*obs+1):(j*obs), ], target="Y")
      trained <- train(learner = lrn, task = task)
      pred <- Predictor$new(trained, data = set[((j-1)*obs+1):(j*obs), ])
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
pdp_plot <- function(pdp_data, avg_data, feature){
  ggplot(data.frame(x = c(-3, 3)), aes(x = x))+
    geom_line(data = pdp_data, aes(x = pdp_data[,1], y = Y_hat, group = Simulation), alpha = 0.3)+
    geom_line(data = avg_data, aes(x = Grid, y = mean), col="black", size = 1.2)+
    geom_errorbar(data = avg_data, aes(x=Grid, ymin = mean-sd, ymax=mean+sd), col = "black")+
    theme_bw()+
    theme(plot.title = element_text(size=16, hjust = 0))+
    xlab(feature)+
    labs(title = "")+
    ylab("")+
    xlim(c(-3,3))+
    ylim(c(-3,3))
}

# Plots Independent Case
p11 <- pdp_plot(pdp_ind_X1, avg_ind_X1, feature = "X1") +
  stat_function(fun = function(x) sin(3*x), geom = "line", col="red", size=1)
p12 <- pdp_plot(pdp_ind_X2, avg_ind_X2, feature = "X2") + geom_abline(slope = 1, col = "red")
p13 <- pdp_plot(pdp_ind_X3, avg_ind_X3, feature = "X3") + geom_abline(slope = 1, col = "red")

ga_ind <- ggarrange(p11+labs(title = "PDPs for independent Features"), 
                    p12, 
                    p13,
                    nrow = 1, ncol = 3)



list(ga_ind)