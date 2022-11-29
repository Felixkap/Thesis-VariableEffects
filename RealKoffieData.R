library(randomForest)
library(creditmodel)
library(mlr)
library(partykit)
library(iml)
data <- read.csv('C:/Users/feix_/iCloudDrive/Studium Master/CQM - Thesis Internship/Thesis-VariableEffects/example project and data/06 example data.csv')
rf <- randomForest(y ~., data = data)
#data.split <- train_test_split(data, prop = 0.75)
#nrow(data.split$test)


task = makeRegrTask(data = data, target = "y")
mod.bike = mlr::train(mlr::makeLearner(cl = 'regr.svm'), task)
pred.bike = Predictor$new(mod.bike, data = data[, names(data) != "y"])
tree = TreeSurrogate$new(pred.bike) 
plot(tree)
pred.tree  = predict(tree, data)
pred.svm = getPredictionResponse(predict(mod.bike, task))

#n_machines <- unique(data$machine)
#n_settings <- unique(data$setting)
#combs <- expand.grid(n_settings, n_machines)

#cnt <- 0
#for (comb in 1:nrow(combs)) {
#  machine_setting <- combs[comb,]
#  idx <- data$setting==as.numeric(machine_setting[1]) & data$machine==as.numeric(machine_setting[2])
#  cnt = cnt + sum(idx)
#}
#cnt

#install.packages("ALEPlot")
#library(ALEPlot)
#library(ggplot2)
#library(randomForest)

#set.seed(4543)
#data(mtcars)
#nrow(mtcars)
#rf.fit <- randomForest(mpg ~ ., data=mtcars)
#yhat <- function(X.model, newdata) as.numeric(predict(X.model, newdata))
#X=mtcars[,-which(names(mtcars) == "mpg")]
#PD.1=PDPlot(X, rf.fit, pred.fun=yhat, J="vs", K=100)
#ALE.2=ALEPlot(X, rf.fit, pred.fun=yhat, J="vs", K=40)


#library(pdp)
#set.seed(4543)
#N=500
#x1 <- sample(0:3,N,replace=TRUE)
#x2 <- runif(N, min=0, max=1)
#x3 <- runif(N, min=0, max=1)
#y = x1 + log(x2)^2 + rnorm(N, 0, 0.1)
#DAT = data.frame(y, x1, x2, x3)
#rf.DAT = randomForest(y ~ ., DAT)

## Define the predictive function (easy in this case, since \code{lm} has 
## a built-in predict function that suffices)
#yhat <- function(X.model, newdata) as.numeric(predict(X.model, newdata))
## Calculate and plot the PD main effects and second-order interaction effects of x1, x2, x3
#plotPartial(pdp::partial(rf.DAT, pred.var="x1"))
#ALE.1=ALEPlot(DAT[,2:4], rf.DAT, pred.fun=yhat, J=1, K=50)
#plotPartial(pdp::partial(rf.DAT, pred.var="x2"))
#ALE.2=ALEPlot(DAT[,2:4], rf.DAT, pred.fun=yhat, J=2, K=400)
#length(ALE.2$x.values)


data(cervical)
cervical.task = makeClassifTask(data = cervical, target = "Biopsy")
mod.cervical = mlr::train(mlr::makeLearner(cl = 'classif.randomForest', predict.type = "prob"), cervical.task)
pred.cervical = Predictor$new(mod.cervical, data = cervical[names(cervical) != "Biopsy"], type = "prob")
tree.cervical = TreeSurrogate$new(pred.cervical, maxdepth = 2) 
plot(tree.cervical) + 
  theme(strip.text.x = element_text(size = 8))
pred.tree.cervical  = predict(tree.cervical, cervical)["Cancer"]
pred.cervical = getPredictionProbabilities(predict(mod.cervical, cervical.task))


