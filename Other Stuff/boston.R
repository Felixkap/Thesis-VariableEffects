require(ICEbox)
require(randomForest)
require(MASS)

data(Boston) #Boston Housing data
X = Boston
y = X$medv
X$medv = NULL


## build a RF:
bhd.rf = randomForest(Boston[, 1:13], Boston$medv, nodesize = 5)
## rf.ice = ice(object = bhd.rf, X = Boston[, 1:13], y = Boston$medv, predictor = "nox", frac_to_build = .2)
## plot(rf.ice)

## build a gbm
require(gbm)
bhd.gbm = gbm(medv ~ ., data = Boston, distribution = "gaussian", n.trees = 5000, interaction.depth = 5, shrinkage = 0.01, cv.folds = 5)
best.iter <- gbm.perf(bhd.gbm, method="cv")
print(best.iter)

## build a bart
require(BayesTree)
bhd.bart <- pdbart(Boston[, 1:13], Boston$medv, levquants = seq(0.01, 0.99, 0.05))
save(bhd.bart, file = "bhd_bart.rda")

pdf("bhd_dis.pdf")
par(mfrow = c(2, 2), mar=c(4, 4, 3, 1) + 0.1)
## plot 1: conditional expectation
plot(Boston$dis, Boston$medv, ann = FALSE, xlim = c(1, 13), ylim = c(5, 50))
lines(loess.smooth(Boston$dis, Boston$medv, span = 1/3), col = "red")
title(main = "Scatter plot", xlab = "nitrix oxides concentration (DIS)", ylab = "median housing price (MEDV)", cex.lab = 1.3, cex.main = 1.5)
## plot 2: random forest
partialPlot(bhd.rf, Boston[, 1:13], x.var = "dis", xlim = c(1, 13), ylim = c(5, 50), ann = FALSE)
title(main = "Random Forest", xlab = "nitrix oxides concentration (DIS)", ylab = "median housing price (MEDV)", cex.lab = 1.3, cex.main = 1.5)
## plot 3: GBM
plot(bhd.gbm, n.trees = best.iter, "dis", main = "GBM", xlim = c(1, 13), ylim = c(5, 50), ann = FALSE)
title(main = "GBM", xlab = "nitrix oxides concentration (DIS)", ylab = "median housing price (MEDV)", cex.lab = 1.3, cex.main = 1.5)
## plot 4: BART
plot(bhd.bart, xind = 8, main = "BART", xlim = c(1, 13), ylim = c(5, 50), ann = FALSE)
title(main = "BART", xlab = "nitrix oxides concentration (DIS)", ylab = "median housing price (MEDV)", cex.lab = 1.3, cex.main = 1.5)
dev.off()


pdf("bhd_nox_ice.pdf")
require(ICEbox)
rf.ice = ice(object = bhd.rf, X = Boston[, 1:13], y = Boston$medv, predictor = "nox", frac_to_build = .2)
plot(rf.ice, ann = FALSE)
title(main = "ICE of Random Forest", xlab = "nitrix oxides concentration (NOX)", ylab = "median housing price (MEDV)", cex.lab = 1.5, cex.main = 1.5)
dev.off()


## scatterplot
pdf("bhd_dis_scatter.pdf")
plot(Boston$dis, Boston$medv, ann = FALSE, xlim = c(1, 13))
lines(loess.smooth(Boston$dis, Boston$medv, span = 1/4), col = "red")
title(main = "Scatter plot", xlab = "distance to city center (DIS)", ylab = "median housing price (MEDV)", cex.lab = 1.5, cex.main = 1.5)
dev.off()

plot(ksmooth(Boston$dis, Boston$medv, bandwidth = 1))

## ## weighted regression
## require(CBPS)
## require(mgcv)
## fit <- CBPS(nox ~ . - medv, data = Boston)

## require(ipw)
## fit <- ipwpoint(nox, family = "gaussian", data = Boston, numerator = ~ 1, denominator = ~ . - medv)
## ## plot(gam(medv ~ s(nox), data = Boston, weights = fit$weights))
## svydata <- Boston
## svydata$weights <- fit$ipw.weights
## require(survey)
## plot(svysmooth(medv ~ nox, design = svydesign(~ 1, weights = ~ weights, data = svydata), bandwidth = 0.05))

## plot(Boston$nox, Boston$medv, xlim = c(1, 0.86), ann = FALSE)
## lines(loess.smooth(Boston$nox, Boston$medv, weights = fit$ipw.weights, span = 1/5))



## mediation plot (total effect)
var <- "dis"
med <- c("zn", "indus", "nox", "rad", "chas")
conf <- setdiff(names(X), c(var, med))
yy <- rep(0, nrow(X))
for (i in 1:nrow(X)) {
    xs <- X[i, var]
    M.pop <- X[order(abs(X[, var] - xs))[1:20], med]
    Z.pop <- X[, conf]
    X.mimic <- cbind(X[i, var, drop = FALSE], Z.pop, M.pop[sample.int(nrow(M.pop), nrow(Z.pop), replace = TRUE), ])
    yy[i] <- mean(predict(bhd.rf, X.mimic))
}

pdf("bhd_dis_tot.pdf")
plot(X[, var], yy, xlab = "dis", ylab = "medv", main = "Total effect of dis")
lines(loess.smooth(X[, var], yy, span = 1/10), col = "red")
dev.off()

## ## mediation analysis (controlled direct effect)
## var <- "dis"
## med <- c("crim", "indus", "nox", "rad", "lstat")
## ## med <- c("rad", "nox")
## conf <- setdiff(names(X), c(var, med))
## yy <- rep(0, nrow(X))
## for (i in 1:nrow(X)) {
##     xs <- X[i, var]
##     M.pop <- colMeans(X[, med])
##     Z.pop <- X[, conf]
##     X.mimic <- cbind(X[i, var, drop = FALSE], Z.pop, t(M.pop))
##     yy[i] <- mean(predict(bhd.rf, X.mimic))
## }

## plot(X[, var], yy)

## ## require(neuralnet)
## ## data <- Boston
## ## maxs <- apply(data, 2, max)
## ## mins <- apply(data, 2, min)
## ## scaled <- as.data.frame(scale(data, center = mins, scale = maxs - mins))
## ## data <- scaled
## ## n <- names(data)
## ## f <- as.formula(paste("medv ~", paste(n[!n %in% "medv"], collapse = " + ")))
## ## nn <- neuralnet(f,data=data,hidden=c(5,3),linear.output=T)

## ## predict.nn <- function(object, ...) {compute(object, ...)$net.result * (maxs["medv"] - mins["medv"]) + mins["medv"]}

## ## nn.ice = ice(object = nn, X = data[, 1:13], y = y, predictor = "dis", frac_to_build = 1)
## ## plot(nn.ice)
