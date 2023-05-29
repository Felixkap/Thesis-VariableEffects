################################################################################
################################################################################
#####################      Import Packages     #################################
################################################################################
################################################################################
rm(list=ls())
set.seed(1234)
library(tidyverse)
library(dplyr)
library(ranger)

# Import Helper File: Change to directory where Helper_workflow.R is saved
source('C:/Users/kpl/Thesis Project/Thesis-VariableEffects/Workflow/Helper_Workflow.R')



################################################################################
################################################################################
#################      Import and Preprocess Data      #########################
################################################################################
################################################################################
data <- read.csv('C:/Users/kpl/Thesis Project/Thesis-VariableEffects/Workflow/Apply_Workflow_To_Data/fietsdata.csv', )
data <- char_to_fac(data)
head(data)
### Variable Selection
# - Use techniques such as correlation analysis, feature importance from a preliminary model 
#   or domain knowledge to guide your feature selection process


x <- data %>% select(-fietskans)



### Outliers: -Base model used in RF is a large decision tree. 
#             -Decision trees are robust to outliers, because they isolate them in small regions of the feature space. 
#             -Since prediction for each leaf is the average (for regression) or the majority class (for classification), 
#              being isolated in separate leaves, outliers won't influence the rest of the predictions. 
#             -Bottom line: We don't care about outliers in RF. 
#             -Only remove them if they are 'wrong' observations (recording errors etc)

### NA-Values: - Ranger Function can not handle NA-Values. Hence use Imputation method
#              - Function impute_missing imputes missing values of numerical variables with mean and 
#                of categorical variable with most common observation




################################################################################
################################################################################
#####################   RF - Hyperparameters  #################################
################################################################################
################################################################################

## Set Hyperparameters manually


## @param num.trees: Number of trees used to fit the Random Forest
## @param mtry: Number of variables to possibly split at in each node 
#               ---> Smaller value -> more randomness <=> Uncorrelated trees ->
## @param node_size: Minimal node size to split at. Default 1 for classification, 5 for regression
## @k: Range over which the effect of numerical variables is computed 
## @sample.fraction: Fraction of observations to sample. Default is 1 for sampling with replacement and 
#                    0.632 for sampling without replacement. 
#                    For classification, this can be a vector of class-specific values.

num.trees<- 3000
mtry <- 2
node_size <- 5
k <- 0.5
sample.fraction <- 1


################################################################################
################################################################################
#################              Fit Models             ##########################
################################################################################
################################################################################

# Random Forest
rf <- ranger( formula = fietskans ~ .,
              data = data,
              num.trees = num.trees,
              keep.inbag = T, 
              oob.error = T, 
              mtry = mtry,
              sample.fraction = sample.fraction,
              importance = 'permutation',
              min.node.size = node_size)

# Linear Model
linreg <- lm( formula = fietskans ~ .^2,
              data = data)
summary(linreg)


################################################################################
################################################################################
#################       Performance on Train Data     ##########################
################################################################################
################################################################################

# Random Forest
all_tree_predictions <- predict(rf, x, predict.all = T)$predictions
predictions.rf <- rowMeans(all_tree_predictions)

# Linear Model
predictions.lm <- predict(linreg,
                          newdata = data,
                          se.fit = T)$fit
n <- nrow(data)
cat("RMSE on train data of Linear Model:", sqrt(mean((data$fietskans - predictions.lm)^2)), 
    ". RMSE on train data of Random Forest:", sqrt(mean((data$fietskans - predictions.rf)^2)),
    ".\nR^2 of Linear Model:", 1 - mean((data$fietskans - predictions.lm)^2) / ((n-1) / (n) * var(data$fietskans)), 
    ". R^2 of Random Forest:", 1 - mean((data$fietskans - predictions.rf)^2) / ((n-1) / (n) * var(data$fietskans)),
    ".\nOOB estimate of RF for generalization error:", sqrt(rf$prediction.error), ".")


# Add optional Arguments: thinning, first_tree
plot_error_rate(rf, data, outcome='fietskans', data_type='train')
plot_error_rate(rf, data, outcome='fietskans', data_type='oob')



################################################################################
################################################################################
#################              RF Summary             ##########################
################################################################################
################# Contains: - Effect Estimates        ##########################
#################           - SE of Effects           ##########################
################################################################################
################################################################################


point_of_interest <- calculate_mean_mode(x)
system.time(result <- summary.rf(x = x, 
                                 fitted.rf = rf, 
                                 point_of_interest = 'mean_mode', # Default <=> calculate_mean_mode(x)
                                 k = k, 
                                 interaction = 'all', 
                                 interaction_terms = NULL, 
                                 moving_var_effects = T, 
                                 moving_se_effects = F, 
                                 first_tree=100, 
                                 thinning=100))

result$summary
#r$cov_list
result$moving_var_effect



################################################################################
################################################################################
##################    PDP's and Variable Importance   ##########################
################################################################################
################################################################################

# PDP of each single variables
plot_pdp(rf, x, names(x))

# PDP of two specific variable
#plot_2pdp(rf, c('geslacht', 'leeftijd'))

# Visualizing Variable Importance
plot_varimp(rf)





# Corrected for all non-linear effects --> wee just fix all other variables at some value
# smaller k, less bias --> condinitional distribution is narrower (can I use plots?)
# Fietsdata: description of data


# Future work