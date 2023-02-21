rm(list=ls())
library(batch)
library(parallel)

detectCores()

### Create clusters for parallel simulation
cl.n <- min(16, detectCores()/2+4)
cl.n <- 6
cl <- makeCluster(cl.n)
source("C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/MERF/Code/parallel_simulation.R")

#read functions and load necessary libraries from other file
parallel::clusterEvalQ(cl,
                       expr = {source("C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/MERF/Code/parallel_simulation.R")})



# an input list for the simulation (creation of this list shown below)
input_list  <- split(sim.seta, seq(nrow(sim.seta)))
output_list <- parLapply(cl, input_list, run.models) #run.models in the main function of the simulation


###### creation of input list: order of factors is important, since this is also the order the factors are assigned in run.models

N <- 10 #number of simulation per batch
nr.group <- c(25,100)
group.size <- c(3,15)
nrX.between <- c(2, 10)
nrX.within <- c(3, 15)
cor <- c(0,0.8)
sigma.b <- c(0.7, 3)
sigma.resid <- c(1, 2)
fixed.formula <- c("X.1+X.2-2*X.2^2-0.5*X.3+1.5*X.4-2*X.5+X.5^2",
                   "X.1-2*X.1^3+0.5*exp(X.2)+2*log(abs(X.3*X.5))-4*(abs(X.4 )> 0.5)"
)

co <- expand.grid(N=N, nr.group = nr.group, group.size = group.size, nrX.between = nrX.between, 
                  nrX.within = nrX.within, cor = cor, sigma.resid = sigma.resid, sigma.b = sigma.b,
                  fixed.formula = fixed.formula,
                  stringsAsFactors = FALSE)
sim.seta <- subset(co, (nrX.between == 2 & nrX.within == 3) | (nrX.between == 10 & nrX.within == 15))


#add batches and scenarios
sim.seta$scenario <- 1:128

sim.seta <- sim.seta[rep(row.names(sim.seta), each = 10),] #repeat every scenario 10 times, so we get 100 simulation in total
sim.seta$seed <- rep(seq(from = 100, by = 50, length.out = 128), 10) #create seeds
sim.seta$batch <- rep(c(1:10), 128) #add batch number






