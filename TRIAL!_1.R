rm(list=ls())
library(batch)
library(parallel)
install.packages('pdp')

detectCores()

### Create clusters for parallel simulation
cl.n <- min(16, detectCores()/2+4)
cl <- length(makeCluster(cl.n))
source("C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thomas_Thesis/Code/TRIAL!_2.R")



### Define Data Setup

N <- 10 # Number of simulation per batch
n_batches <- 5 # For parallelization
nr.group <- c(5) # Number of groups
group.size <- c(15) # Group Size
nrX.between <- c(2) # Variables to vary between groups
nrX.within <- c(3) # Variables to vary within groups
cor <- c(0) # Correlation between predictors
sigma.b <- c(0.7) # random group effect
sigma.resid <- c(1) # random noise
fixed.formula <- c("X.1+X.2-2*X.2^2-0.5*X.3+1.5*X.4-2*X.5+X.5^2",
                   "2*X.1^3-0.5*exp(X.2)+X.3-2*log(abs(X.5))-4*(abs(X.4 )> 0.5)"
)

### Collect all Data Setups
sim.seta <- expand.grid(N=N, 
                  nr.group = nr.group, 
                  group.size = group.size, 
                  nrX.between = nrX.between, 
                  nrX.within = nrX.within, 
                  cor = cor, 
                  sigma.resid = sigma.resid, 
                  sigma.b = sigma.b,
                  fixed.formula = fixed.formula,
                  stringsAsFactors = FALSE)


n_scenarios <- nrow(sim.seta)
sim.seta$scenario <- 1:n_scenarios # add scenarios

sim.seta <- sim.seta[rep(row.names(sim.seta), each = n_batches),] # Repeat every scenario n_batches times (for parallelization)
sim.seta$seed <- rep(seq(from = 100, by = 50, length.out = n_batches), n_scenarios) # Add seeds
sim.seta$batch <- rep(c(1:n_batches), n_scenarios) # Add batch number


input_list  <- split(sim.seta, seq(nrow(sim.seta)))
saved_data <- save.simulated.data(input_list[[1]])




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
  
  
  for(i in 1:N){ # For each Simulation...
    set.seed(seed = seed) 
    data <- simulate.data.gaussian(nr.group,
                                   group.size,
                                   nrX.between, 
                                   nrX.within,
                                   cor,
                                   sigma.resid,
                                   sigma.b,
                                   fixed.formula)
    
    data.wg <- data[, - ncol(data)] # Get rid of group column
    
    if ( i == 1) { #only do this first time, create model formula (for lm and lme)
      formula <- as.formula(
        paste('y ~',paste('poly(',colnames(data.wg[-length(colnames(data.wg))]),',2)',collapse = ' + ')
        )) #formula that's fitted for LM and LME models, it only include first and second order terms, no tf or interactions
    }
    
    ####### LM #############
    #t <- system.time(LM <- lm(formula, data=data))
    #run.times[i, 'LM'] <- unname(t[3])
    LM <- lm(formula, data=data)
    slm <- summary(LM)
    
    ####### LME ############
    
    #t <- system.time(LME <- lmer(update(formula,    ~ . + (1|grp)), data = data))#add Random Intercept for mixed-effects model
    #run.times[i, 'LME'] <- unname(t[3])
    LME <- lmer(update(formula,    ~ . + (1|grp)), data = data) #add Random Intercept for mixed-effects model
    slme <- summary(LME)
    
  }
}

data <- saved_data[[1]]
form <-  paste('y ~',paste('poly(',colnames(data[-length(colnames(data))]),',2)',collapse = ' + '))

LM <- lm(form, data)
###X1
range(data['X.1'])
x1_range <- seq(from=-3, to=3, length.out=100)
predictions_new <- numeric(length(x1_range))
for (val in seq_along(x1_range)) {
  data_new <- data
  data_new['X.1'] <- x1_range[val]
  print(mean(predict(LM, data_new)))
  predictions_new[val] <- mean(predict(LM, data_new))
}



