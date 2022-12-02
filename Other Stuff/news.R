data <- read.csv("OnlineNewsPopularity.csv")

predictors <- c(2:19, 29:44, 45:50, 57:58)

## s <- sample(1:nrow(data), 20000)
## datas <- data[s, ]

require(randomForest)
news.rf <- randomForest(data[, predictors], log(data$shares), nodesize = 20)
mean((log(data$shares[-s]) - predict(news.rf, data[-s, ]))^2)

require(gbm)
news.gbm <- gbm(log(shares) ~ . , data = cbind(data[, predictors], shares = data$shares), distribution = "gaussian", n.trees = 1000, interaction.depth = 10, shrinkage = 0.005, cv.folds = 3)
best.iter <- gbm.perf(mpg.gbm, method="cv")
print(best.iter)

require(BayesTree)
mpg.bart <- pdbart(data[, predictors], log(data$shares), levquants = seq(0.01, 0.99, 0.05))

variable <- "title_sentiment_polarity"
variable.lab <- "title sentiment polarity"
ylim <- c(7.25, 7.6)
## variable <- "title_subjectivity"
## variable <- "n_tokens_title"
variable <- "num_keywords"
variable.lab <- "number of keywords"
ylim <- c(7.1, 7.5)

plot(data[, variable], log(data$shares), ann = FALSE)
lines(loess.smooth(data[, variable], log(data$shares), span = 2/3), col = "red")
title(main = "Scatter plot", xlab = variable.lab, ylab = "log(shares)")

pdf(paste0("news_", variable, ".pdf"), width = 8, height = 5)
par(mfrow = c(1, 2))
plot(loess.smooth(data[, variable], log(data$shares), span = 2/3), type = "l", main = "LOESS (conditional expectation)", xlab = variable.lab, ylab = "log(shares)", cex.lab = 1.5, ylim = ylim)
partialPlot(news.rf, data, x.var = eval(variable), ann = FALSE, ylim = ylim)
title(main = "Random Forest (partial dependence)", xlab = variable.lab, ylab = "log(shares)", cex.lab = 1.5)
dev.off()

plot(news.gbm, n.trees = best.iter, variable, main = "GBM", ann = FALSE)
title(main = "GBM", xlab = variable, ylab = "shares")


## partialPlot(news.rf, data, x.var = "num_hrefs", ann = FALSE)
## title(main = "Random Forest", xlab = "num_hrefs", ylab = "log(shares)")
