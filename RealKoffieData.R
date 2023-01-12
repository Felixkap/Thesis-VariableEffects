rm(list=ls())
library(ranger)
library(creditmodel)
library(mlr)
library(partykit)
library(iml)
library(Matrix)
library(ggplot2)
library(ALEPlot)
library(caret)

source('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/RangerPredictFunction.R')



#### Predict Function
yhat <- function(X.model, newdata) {
  
  ### Predict RF Responses (Training Data)
  rf.predict <- RangerForestPredict(X.model$forest, 
                                    newdata, 
                                    predict.all = T,
                                    inbag.counts = X.model$inbag.counts)
  
  
  # For each observation we get all predictions from all trees
  rf.all.predictions <- rf.predict$predictions
  
  # For each observation average predctions over all trees
  rf.predictions <- rowMeans(rf.all.predictions)
  
  return(rf.predictions)
}



##### Import Data
data <- read.csv('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/example project and data/06 example data.csv', )
# Only numerical data
data_main = data[,-which(names(data) %in% c("machine", "test", "setting"))]

data_onehot = data[,-which(names(data) %in% 'test')]
data_onehot[,c('setting', 'machine')] <- sapply(data_onehot[,c('setting', 'machine')], as.character)
dmy <- dummyVars("~.", data = data_onehot)
data_onehot <- data.frame(predict(dmy, newdata = data_onehot))


X_main = data_main[,-which(names(data_main) == "y")]
X_onehot = data_onehot[,-which(names(data_onehot) == "y")]


### Fit Random Forest Model from 'ranger'
rf <- ranger( formula = y ~ ., 
              data = data_main, 
              num.trees = 200, 
              keep.inbag = T, 
              min.node.size = 5) # min.node.size to smooth forest fit

rf_onehot <- ranger( formula = y ~ ., 
                     data = data_onehot, 
                     num.trees = 200, 
                     keep.inbag = T, 
                     min.node.size = 5) # min.node.size to smooth forest fit


# For each observation average predctions over all trees
rf.predictions <- yhat(rf, X_main)
rf.predictions.onehot <- yhat(rf_onehot, X_onehot)

#### Standard Error of Random Forest Predictions
rf.se <- RangerForestPredict(rf$forest, 
                             data_main, 
                             type = 'se', 
                             se.method = 'jack', 
                             inbag.counts = rf$inbag.counts)$se

rf.se.onehot <- RangerForestPredict(rf_onehot$forest, 
                             data_onehot, 
                             type = 'se', 
                             se.method = 'jack', 
                             inbag.counts = rf$inbag.counts)$se


data_main <- cbind(data_main, rf.predictions, rf.se)
data_onehot <- cbind(data_onehot, rf.predictions.onehot, rf.se.onehot)


gg_se <- ggplot(data = data, aes(x=y,y=rf.predictions)) + 
  geom_point(alpha=0.3) +
  geom_errorbar(aes(ymin=rf.predictions - 0.5*rf.se,
                    ymax=rf.predictions + 0.5*rf.se), width=.1)+
  geom_abline(slope = 1, linetype = 2 , color = 'red')+
  labs(x = 'True Y', y = 'Predicted Y', 
       title = 'Random Forest Prediction with Standard Errors 
       vs. True Outcome Variable')+
  theme_bw()


gg_se



gg_se.onehot <- ggplot(data = data_onehot, aes(x=y,y=rf.predictions.onehot)) + 
  geom_point(alpha=0.3) +
  geom_errorbar(aes(ymin=rf.predictions.onehot - 0.5*rf.se.onehot,
                    ymax=rf.predictions.onehot + 0.5*rf.se.onehot), width=.1)+
  geom_abline(slope = 1, linetype = 2 , color = 'red')+
  labs(x = 'True Y', y = 'Predicted Y', 
       title = 'Random Forest Prediction with Standard Errors 
       vs. True Outcome Variable')+
  theme_bw()


gg_se.onehot



### PDP and ALE Plots for all x variables
par("mar")
par(mar=c(1,1,1,1))
par(mfrow=c(4,2))
PD.1=PDPlot(X, rf, pred.fun=yhat, J="x1", K=100)
ALE.2=ALEPlot(X, rf, pred.fun=yhat, J="x1", K=100)
PD.1=PDPlot(X, rf, pred.fun=yhat, J="x2", K=100)
ALE.2=ALEPlot(X, rf, pred.fun=yhat, J="x2", K=100)
PD.1=PDPlot(X, rf, pred.fun=yhat, J="x3", K=100)
ALE.2=ALEPlot(X, rf, pred.fun=yhat, J="x3", K=100)
PD.1=PDPlot(X, rf, pred.fun=yhat, J="x4", K=100)
ALE.2=ALEPlot(X, rf, pred.fun=yhat, J="x4", K=100)




