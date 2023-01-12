################################################################################
################################################################################
########################        MY SIMULATION       ############################
################################################################################
################################################################################

rm(list=ls())
setwd('C://Users/feix_//iCloudDrive//Studium Master//CQM - Thesis Internship//Thesis-VariableEffects')

source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/SE_Variable_Effect_RF.R')
library(parallel)

cores <- detectCores()
clust <- makeCluster(cores)
# Export all objects in R-Script 'Covariance Predictions R' 
parallel::clusterEvalQ(clust,
                       expr = {source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/SE_Variable_Effect_RF.R')})



###### Simulation Setup
formulas <- c("2*x+2*-x^2+e")
n <- c(100, 200, 500, 1000)
num.trees <- 200
repeats <- 1e3
scenarios <- data.frame(expand.grid(n, num.trees, formulas, repeats))
colnames(scenarios) = c("N", "N_Trees", "Formula", "Repeats")
scenarios[,"Formula"] <- as.character(scenarios[,"Formula"]) ### Formula became Factor
scenarios <- split(scenarios, seq(nrow(scenarios)))





# Simulate
system.time(result <- parLapply(cl = clust, 
                                X = scenarios, 
                                fun = sim_multi))

test <- do.call(rbind, result)
colnames(test) <- c("Formula", "N", "N_Trees",
                     "Mean Effect", "SD Effect", "SE Estimate Effect", 
                     "SE2 Estimate Effect", "N_Nulls")



### Write Table
write.table(test, file = 'NoBiasCorrection - 1e3 - cbindxe.csv', sep = ",", row.names = F)