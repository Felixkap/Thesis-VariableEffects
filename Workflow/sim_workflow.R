rm(list=ls())


################################################################################
################################################################################
#####################      Import Packages     #################################
################################################################################
################################################################################
library(ranger)
library(mlr)

# Import Helper File: Change to directory where Helper_workflow.R is saved
source('C:/Users/kpl/Thesis Project/Thesis-VariableEffects/Workflow/Helper_Workflow.R')



################################################################################
################################################################################
#####################       Generate Data      #################################
################################################################################
################################################################################

## @param seed: Seed Number to obtain same data set
## @param N_train: Number of training observations
## @param N_test: Number of test observations
## @param cor: Correlation between features
## @param sigma_e: Noise in the data (standard deviation of errors)
## @param formula: Functional Relationship between predictors and outcome
## @param cat: Number of binary categorical variables to include as covariates
## @param k: Range over which the variable effect is to be computed

seed <- 12391
set.seed(seed)
N_train <- 500
N_test <- 200
cor <- 0
sigma_e <- 1
formula <- c("4*x.1+0.2*x.2-5*x.3+2*x.4-2*x.3*x.4")
cat <- 2
k <- 1 

# Obtain number of variables: n_vars
n_vars <- length(unique(unlist(str_extract_all(formula,"x.[0-9]+")))) 


train_data <- generate_data(N_train, n_vars, cor, formula, sigma_e, cat=cat)

# Vector containing mean of all numerical and mode of all categorical variables
print(train_data$point_of_interest) 


################################################################################
################################################################################
#####################   RF - Hyperparameters  #################################
################################################################################
################################################################################

## @param num.trees: Number of trees used to fit the Random Forest
## @param mtry: Number of variables to possibly split at in each node 
#               ---> Smaller value -> more randomness <=> Uncorrelated trees ->
## @param node_size: Minimal node size to split at. Default 1 for classification, 5 for regression



## Set Hyperparameters manually (Otherwise Hyperparameter Tuning)
num.trees <- 2500
mtry <- 1
node_size <- 5


################################################################################
################################################################################
#################              Fit Models             ##########################
################################################################################
################################################################################


# Random Forest
rf <- ranger( formula = y ~ .,
              data = train_data$data,
              num.trees = num.trees,
              keep.inbag = T, 
              oob.error = T, 
              importance = 'permutation',
              min.node.size = node_size)


# Linear Model
linreg <- lm( formula = y ~ .^2,
              data = train_data$data)



################################################################################
################################################################################
#################       Performance on Train Data     ##########################
################################################################################
################################################################################

# Random Forest
all_tree_predictions <- predict(rf, train_data$x, predict.all = T)$predictions
predictions.rf <- rowMeans(all_tree_predictions)

# Linear Model
predictions.lm <- predict(linreg,
                          newdata = train_data$x,
                          se.fit = T)$fit

cat("RMSE on train data of Linear Model:", sqrt(mean((train_data$y - predictions.lm)^2)), 
    ". RMSE on train data of Random Forest:", sqrt(mean((train_data$y - predictions.rf)^2)),
    ".\nR^2 of Linear Model:", 1 - mean((train_data$y - predictions.lm)^2) / var(train_data$y), 
    ". R^2 of Random Forest:", 1 - mean((train_data$y - predictions.rf)^2) / var(train_data$y),
    ".\nOOB estimate of RF for generalization error:", sqrt(rf$prediction.error), ".")



# Plot RMSE vs Number of Trees (train data)
plot_error_rate(rf, train_data$data, data_type='train')

################################################################################
################################################################################
#################       Performance on Test Data      ##########################
################################################################################
################################################################################

test_data <- generate_data(N_test, n_vars, cor, formula, sigma_e, cat=cat)


all_tree_predictions <- predict(rf, test_data$x, predict.all = T)$predictions
predictions.rf <- rowMeans(all_tree_predictions)

predictions.lm <- predict(linreg,
                          newdata = test_data$x,
                          se.fit = T)$fit

cat("RMSE on test data of Linear Model:", sqrt(mean((test_data$y - predictions.lm)^2)), 
    ". RMSE on test data of Random Forest:", sqrt(mean((test_data$y - predictions.rf)^2)))


# Plot RMSE vs Number of Trees (train data)
plot_error_rate(rf, test_data$data, data_type='test')




################################################################################
################################################################################
#######                        RF Summary                      #################
################################################################################
#######    Contains: - Effect Estimates                        #################
#######              - SE of Effects                           #################
#######              - Data Points used to compute effects     #################
#######                and Variance-Covariance of              #################
#######                corresponding RF Prediction             #################
#######              - Moving SE/VAR Effect vs Trees           #################
################################################################################
################################################################################

### Definition of Main Effect for variable x_j (Example for point_of_interest = 'mean_mode')
# (i) Form two points 
#     A:= (mean(x_1), mean(x_2), ... ,mean(x_j) - k * sd(x_j), ... mean(x_P) )
#     B:= (mean(x_1), mean(x_2), ... ,mean(x_j) + k * sd(x_j), ... mean(x_P) )
# (ii) Compute Main Effect: beta_j = (f(B) - f(A)) / (2 * k * sd(x_j))


### Definition for low-order interaction effect between variables x_j and x_h
# (i) Form four points 
#     A0:= (mean(x_1), mean(x_2), ... ,mean(x_j) - k * sd(x_j), mean(x_h) - k * sd(x_h), ... mean(x_P) )
#     A1:= (mean(x_1), mean(x_2), ... ,mean(x_j) - k * sd(x_j), mean(x_h) + k * sd(x_h), ... mean(x_P) )
#     B0:= (mean(x_1), mean(x_2), ... ,mean(x_j) + k * sd(x_j), mean(x_h) - k * sd(x_h), ... mean(x_P) )
#     B1:= (mean(x_1), mean(x_2), ... ,mean(x_j) + k * sd(x_j), mean(x_h) + k * sd(x_h), ... mean(x_P) )
# (ii) Compute low-order Interaction Effect: gamma_{j,h} = ((f(B1) - f(B0)) - (f(A1) - f(A0))) / (4 * k^2 * sd(x_j) * sd(x_h))



system.time(result <- summary.rf(x = train_data$x, 
               fitted.rf = rf, 
               point_of_interest = 'mean_mode',
               k = k, 
               interaction = 'specific', 
               interaction_terms = 'x.1:x.2', 
               moving_var_effects = T, 
               moving_se_effects = F, 
               first_tree=100, 
               thinning=100))

result$summary
result$cov_list
result$moving_var_effect



### Use different data point where we evaluate effects
#   --> change value of numeric variable
#   --> change level of binary categorical variable
point_of_interest <- train_data$point_of_interest
print(point_of_interest)

point_of_interest[,1] <- 1
point_of_interest <- other_level(train_data$point_of_interest, 4)
print(point_of_interest)


system.time(result <- summary.rf(x = train_data$x, 
                                 fitted.rf = rf, 
                                 point_of_interest = point_of_interest,
                                 k = k, 
                                 interaction = 'specific', 
                                 interaction_terms = 'x.1:x.3', 
                                 moving_var_effects = T, 
                                 moving_se_effects = F, 
                                 first_tree=100, 
                                 thinning=100))

result$summary



################################################################################
################################################################################
##################    PDP's and Variable Importance   ##########################
################################################################################
################################################################################

# PDP of each single variables
plot_pdp(rf, train_data$x, names(train_data$x))

# PDP of two specific variable
plot_2pdp(rf, c('x.3', 'x.1'))

# Visualizing Variable Importance
plot_varimp(rf)






################################################################################
################################################################################
##########                      ILLUSTRATION                         ###########
################################################################################
##########           Obtaining Variance-Covariance Matrix            ###########
##########                  of any RF Predictions                    ###########
################################################################################
################################################################################
################################################################################

# Step 1: Define Data points of interest
newdata <- data.frame(x.1 = c(2, 1.5, 0.333, -0.8), 
                      x.2 = c(1, 1.7, -0.6, 3), 
                      x.3 = c(0.5, -2.9, 0.1, 0.9), 
                      x.4 = factor(c(0, 1, 1, 1), levels = c(0,1)))

# Step 2: Apply ranger object with new data to RangerForestPredict() 
# - se.method = 'jack_cov': Jackknife-after-Bootstrap to estimate standard errors of 
#                          and covariance between RF predictions
# - se.method = 'se_direct': Estimate SE of variable effect (based on data points provided by newdata) directly

obj <- RangerForestPredict(object = rf$forest, 
                           data = newdata,  
                           type = 'se', 
                           se.method = 'jack_cov', 
                           inbag.counts = rf$inbag.counts)

# Predictions: Average over tree predictions
rowMeans(obj$predictions)

# Variance Covariance Matrix of random forest predictions
obj$cov



