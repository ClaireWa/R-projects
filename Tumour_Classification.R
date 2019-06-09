library(readr)
library(caret)
library(glmnet)
library(MASS)
library(pROC)
library(QuantPsyc)
library(data.table)

# Part A

# read in file
wdbc2 <- read.csv("wdbc2.csv")

# set seed so data partition remains constant each time code is ran
set.seed(1)
# create data partition s.t. training set contains 70% of the observations
train.idx <- createDataPartition(wdbc2$diagnosis, p=0.7)$Resample1

# Fit both a ridge regression model and a lasso model on the training set to diagnose
# the type of tumour from the 30 biomarkers
# Ridge regression and lasso are implemented in the glmnet package
# use the following function to facilitate the tranformation of a dataframe to a 
# matrix as expected by the glmnet package.

prepare.glmnet <- function(data, formula=~ .) {
  
  ## create the design matrix to deal correctly with factor variables,
  ## without losing rows containing NAs
  old.opts <- options(na.action='na.pass')
  x <- model.matrix(formula, data)
  options(old.opts)
  
  ## remove the intercept column, as glmnet will add one by default
  x <- x[, -match("(Intercept)", colnames(x))]
  
  return(x)
}

# By default, the function uses all existing columns in the dataframe 
# However, we do not want the outcome variable to be in the matrix of predictors: 
# so first remove it, then convert the rest of the dataframe to a matrix.
y.wdbc <- wdbc2$diagnosis[train.idx]      # store the outcome separately
x.wdbc <- prepare.glmnet(wdbc2[train.idx,], ~ . - diagnosis-id) # exclude the outcome

# check for missing values
table(!is.na(wdbc2))
# check if columns in x.wdbc are numeric
table(apply(x.wdbc, 2, is.numeric))

# Fit the ridge regression and lasso model, performing cross-validation 
# for each model. by default, function cv.glmnet() performs cross=validation with 
# lasso penalty. To change it to ridge regression, set the alpha option to 0.

# fit ridge
# set seed again so fitted values remain constant 
set.seed(1)
fit.cv.ridge <- cv.glmnet(x.wdbc, y.wdbc, alpha = 0, type.measure = "auc", family = "binomial")

# fit lasso
# set seed again so fitted values remain constant 
set.seed(1)
fit.cv.lasso <- cv.glmnet(x.wdbc, y.wdbc, family = "binomial", type.measure = "auc")

# Get the lambda that maximises the AUC for each model
ridge.opt <- fit.cv.ridge$lambda.min
lasso.opt <- fit.cv.lasso$lambda.min

####################################################################################

# Part B

# Extract the AUC corresponding to the optimal lambda
AUC.opt.r <- fit.cv.ridge$cvm[fit.cv.ridge$lambda == ridge.opt]
AUC.opt.l <- fit.cv.lasso$cvm[fit.cv.lasso$lambda == lasso.opt]

# Extract the AUC corresponding to the lambda such that the AUC is within
# 1 standard error of the maximum
ridge.1se <- fit.cv.ridge$lambda.1se
lasso.1se <- fit.cv.lasso$lambda.1se
AUC.r.1se <- fit.cv.ridge$cvm[fit.cv.ridge$lambda == ridge.1se]
AUC.l.1se <- fit.cv.lasso$cvm[fit.cv.lasso$lambda == lasso.1se]

####################################################################################

# Part C

# define the coefficients for each of the lambdas w.r.t both models
ridge.coeffs <- fit.cv.ridge$glmnet.fit$beta[,which(fit.cv.ridge$lambda == ridge.opt)]
lasso.coeffs <- fit.cv.lasso$glmnet.fit$beta[,which(fit.cv.lasso$lambda == lasso.opt)]
ridge.coeffs2 <- fit.cv.ridge$glmnet.fit$beta[,which(fit.cv.ridge$lambda == ridge.1se)]
lasso.coeffs2 <- fit.cv.lasso$glmnet.fit$beta[,which(fit.cv.lasso$lambda == lasso.1se)]

# Get the model sizes
l1 <- length(ridge.coeffs[abs(ridge.coeffs)>0])
l2 <- length(lasso.coeffs[abs(lasso.coeffs)>0])
l3 <- length(ridge.coeffs2[abs(ridge.coeffs2)>0])
l4 <- length(lasso.coeffs2[abs(lasso.coeffs2)>0])

# create a data table for each type of lambda
dt1 <- data.table("Model" = c("Ridge Regression", "Lasso Regression"), "Lambda Type" = c("Optimal"),
                 "Lambda Value" = signif(c(ridge.opt, lasso.opt), 3), "Model Size" = c(l1, l2), 
                 "AUC" = signif(c(AUC.opt.r, AUC.opt.l), 3))

dt2 <- data.table("Model" = c("Ridge Regression", "Lasso Regression"), "Lambda Type" = c("AUC within 1 se of max"),
                 "Lambda Value" = signif(c(ridge.1se, lasso.1se),3),
                 "Model Size" = c(l3, l4), 
                 "AUC" = signif(c(AUC.r.1se, AUC.l.1se),3))

# combine the data tables into one
lambdas.dt <- rbind(dt1, dt2)

####################################################################################

# Part D

# Perform backward elimination on the same training set derived in Part A
# define the training set
wdbc.training <- wdbc2[train.idx,]
wdbc.training$id <- NULL

# check for missing values
table(!is.na(wdbc.training))
# check if columns which are used as dependent variables are not numeric
table(apply(wdbc.training[2:length(wdbc.training)], 2, is.numeric))

# Use all possible variables in the model
full.model <- glm(diagnosis ~ ., data = wdbc.training, family = "binomial")
 
# perform backward elimination
B <- stepAIC(full.model, direction = "back")
summary(B)

# standardize the regression coefficients using lm.beta() function
# Order the coefficients in decreasing order of absolute value
B.std.coefs <- lm.beta(B)
B.std.coefs <- B.std.coefs[order(abs(B.std.coefs), decreasing = TRUE)]
B.std.coefs

#############################################################################################

# Part E
# Perform stepwise selection starting from the null model
null.model <- glm(diagnosis ~ 1, data = wdbc.training, family = "binomial")

S <- stepAIC(null.model, scope = list(upper=full.model), direction = "both")
summary(S)

S.std.coefs <- lm.beta(S)
S.std.coefs <- S.std.coefs[order(abs(S.std.coefs), decreasing = TRUE)]
S.std.coefs
# check which variable(s) entered the model and were later discarded

###############################################################################################

# Part G
# Use models S and B to predict obs in the training set and from that compute the AUC
B.pred <- predict(B, data = wdbc2[train.idx,], type = "response")
S.pred <- predict(S, data = wdbc2[train.idx,], type = "response")

# compute the training AUC for each model
roc.B <- roc(wdbc2$diagnosis[train.idx], B.pred)
roc.S <- roc(wdbc2$diagnosis[train.idx], S.pred)
roc.B$auc
roc.S$auc

##############################################################################################

#Part H
# Use the four models to predict the outcome for the observations in the test set
# use lambda at 1 se for the penalised models
x.pred <- prepare.glmnet(wdbc2[-train.idx,], ~ . - diagnosis - id)
y.pred <- wdbc2$diagnosis[-train.idx]

# Predict the outcome for the test set
test.pred.B <- predict(B, wdbc2[-train.idx,])
test.pred.S <- predict(S, wdbc2[-train.idx,])
test.pred.lasso <- predict(fit.cv.lasso, newx = x.pred, type = "response", s = lasso.1se)
test.pred.ridge <- predict(fit.cv.ridge, newx = x.pred, type = "response", s = ridge.1se)

# Plot the ROC curves for each model
# on the same plot using different colours
# and report the test AUCs
roc.test.B <- roc(wdbc2$diagnosis[-train.idx], test.pred.B,
                  plot = TRUE, col = "black")
roc.test.S <- roc(wdbc2$diagnosis[-train.idx], test.pred.S,
                  plot = TRUE, add = TRUE, col = "red")
roc.test.lasso <- roc(y.pred, as.numeric(test.pred.lasso),
                      plot = TRUE, add = TRUE, col = "green")
roc.test.ridge <- roc(y.pred, as.numeric(test.pred.ridge),
                      plot = TRUE, add = TRUE, col = "blue")
legend = c("B Model", "S Model", "Lasso Model", "Ridge Model")
legend("right", legend, pch = 4, col = 1:4)
roc.test.B$auc
roc.test.S$auc
roc.test.lasso$auc
roc.test.ridge$auc
