library(ggplot2)
data <- mpg
colnames(data) <- c("mpg", "cylinders", "displacement", "horsepower", "weight", "acceleration", "year", "origin", "name")
data$horsepower <- as.numeric(data$horsepower)
data$origin <- as.factor(data$origin)
data$origin.us <- as.numeric(data$origin == 1)
## data$origin.eu <- as.numeric(data$origin == 2)
data$origin.jp <- as.numeric(data$origin == 3)
data$origin <- NULL
data$name <- NULL

require(randomForest)
mpg.rf <- randomForest(data[, 2:9], data$mpg, nodesize = 5)

require(gbm)
mpg.gbm <- gbm(mpg ~ . , data = data, distribution = "gaussian", n.trees = 5000, interaction.depth = 5, shrinkage = 0.01)
best.iter <- gbm.perf(mpg.gbm, method="OOB")
print(best.iter)

require(BayesTree)
mpg.bart <- pdbart(data[, 2:9], data$mpg, levquants = seq(0.01, 0.99, 0.05))
save(mpg.bart, file = "bhd_bart.rda")

variable <- "acceleration"

pdf(paste0("mpg_", variable, ".pdf"))

par(mfrow = c(2, 2), mar=c(4, 4, 3, 1) + 0.1)

if (class(data[, variable]) %in% c("integer", "factor")) {
    boxplot(data$mpg ~ data[, variable], ann = FALSE)
    title(main = "Box plot", xlab = variable, ylab = "mpg", cex.lab = 1.5, cex.main = 1.5)
} else {
    plot(data[, variable], data$mpg, ann = FALSE, ylim = c(10, 40), xlim = c(5, 25))
    lines(loess.smooth(data[, variable], data$mpg, span = 2/3), col = "red")
    title(main = "Scatter plot", xlab = variable, ylab = "mpg", cex.lab = 1.5, cex.main = 1.5)
}

partialPlot(mpg.rf, data, x.var = eval(variable), ann = FALSE, ylim = c(20, 28), xlim = c(5, 25))
title(main = "Random Forest", xlab = variable, ylab = "mpg", cex.lab = 1.5, cex.main = 1.5)

plot(mpg.gbm, n.trees = best.iter, variable, main = "GBM", ann = FALSE, ylim = c(20, 28), xlim = c(5, 25))
title(main = "GBM", xlab = variable, ylab = "mpg", cex.lab = 1.5, cex.main = 1.5)

plot(mpg.bart, xind = which(colnames(data) == variable) - 1, main = "BART", ann = FALSE, ylim = c(20, 28), xlim = c(5, 25))
title(main = "BART", xlab = variable, ylab = "mpg", cex.lab = 1.5, cex.main = 1.5)

dev.off()


pdf("mpg_acceleration_ice.pdf")
require(ICEbox)
rf.ice = ice(object = mpg.rf, X = data[, 2:9], y = data$mpg, predictor = "acceleration", frac_to_build = .5)
plot(rf.ice, ann = FALSE)
title(main = "ICE of Random Forest", xlab = "acceleration", ylab = "mpg", cex.lab = 1.5, cex.main = 1.5)
dev.off()

pdf("mpg_origin_ice.pdf")
par(mfrow = c(1, 2))
require(ICEbox)
rf.ice = ice(object = mpg.rf, X = data[data$origin.jp == 0, 2:9], y = data$mpg[data$origin.jp == 0], predictor = "origin.us")
plot(rf.ice, ann = FALSE, frac_to_plot = 0.5)
title(main = "ICE of Random Forest", xlab = "origin is US", ylab = "mpg", cex.lab = 1.5, cex.main = 1.5)
rf.ice = ice(object = mpg.rf, X = data[data$origin.us == 0, 2:9], y = data$mpg[data$origin.us == 0], predictor = "origin.jp")
plot(rf.ice, ann = FALSE, frac_to_plot = 1)
title(main = "ICE of Random Forest", xlab = "origin is JP", ylab = "mpg", cex.lab = 1.5, cex.main = 1.5)
dev.off()



pdf("mpg_origin_scatter.pdf")
levels(data$origin) <- c("US", "EU", "JP")
boxplot(data$mpg ~ data$origin, ann = FALSE)
title(main = "Box plot", xlab = "origin", ylab = "mpg", cex.lab = 1.5, cex.main = 1.5)
dev.off()
