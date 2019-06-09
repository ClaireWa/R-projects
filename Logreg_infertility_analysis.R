library(readr)
library(data.table)
library(caret)

# Part A

# Import the dataset and convert to data table
infert.dt <- data.table(infert)

# Fit a logistic regression model (M1) to predict secondary infertility using age and 
# parity as predictors
M1 <- glm(case ~ age + parity, family = "binomial", data = infert)

# Use deviance to judge the goodness of fit and report a a p-value to 3 significant figures
# degrees of freedom = 2 (difference in parameters of the two models)
summary(M1)
signif(pchisq(M1$null.deviance - M1$deviance, df=2, lower.tail = FALSE), 3)

#####################################################################################

# Part B

# Fit a second model (M2) by adding the number of spontaneous abortions
# to the set of predictors used in model M1
M2 <- glm(case ~ age + parity + spontaneous, family = "binomial", data = infert)
summary(M2)

# Report the odds ratio and the 95% confidence interval for the spontaneous abortions
# variable
or.spontaneous <- exp(coef(M2)[4])
ci.spontaneous <- exp(confint(M2))[4,]

# Note: report confint() function answers instead of the following way below as generally 
# gives better results especially for smaller sample sizes (ref: Lab 3, pg 2.)
round(exp(coef(M2)[4] + 1.96 * (coef(summary(M2))[4,2]) * c(-1,1)), 3)

# Perform a likelihood ratio test to compare model M2 to M1
# and report the p-value for the test
# df = 1 (difference in parameters of the two models)
pval <- pchisq(M1$deviance - M2$deviance, df = 1, lower.tail = FALSE)
pval
signif(pval, 3)

# Given that the p-value is < 0.05, the model that includes 
# spontaneous variable is significantly better.

####################################################################################

# Part C

# Implement a function that computes the binomial log-likelihood

loglik.binom <- function(y.obs, y.pred) {
  
  # INPUT:
  #     y.obs: A vector of observed outcomes with values either 0 or 1 to represent
  #            controls and cases, respectively
  #     y.pred: A vector of fitted probabilities learnt from a logistic regression model
  #
  # OUTPUT: 
  #     The binomial log-likelihood for the outcome variable, y.
  
  case <- y.obs == 1
  case.p <- y.pred[case]
  L1 <- sum(log(case.p))
  
  ctrl <- y.obs == 0
  ctrl.p <- y.pred[ctrl]
  L2 <- sum(log(1 - ctrl.p))
  
  L <- L1 + L2
  
  return(L)
}

# Use the function to compute the deviance and null deviance for model M2
# To calculate the deviance for M2, use the following formula:
# D = 2*(log(Lstar) - log(L(Beta)))
# Here, log(Lstar) = 0
# The deviance of a model is minus twice the log-likelihood
# The null deviance is minus twice the log-likelihood of a model which only uses 
# the intercept term and no other predictor

dev.M2 <- 2*(- loglik.binom(infert$case, M2$fitted.values))

null.dev.M2 <- 2*(- loglik.binom(infert$case, rep(mean(infert$case), length(infert$case))))

summary(M2)

# These results match the deviance and null deviance in the summary output of M2
# when rounded to two decimal places

#####################################################################################

# Part D

# Function glm.cv() from Lab 3
glm.cv <- function(formula, data, folds) {
  regr.cv <- NULL
  for (fold in 1:length(folds)) {
    regr.cv[[fold]] <- glm(formula, data=data[-folds[[fold]], ],
                           family="binomial")
  }
  return(regr.cv)
}

# Function predict.cv() from Lab 3
predict.cv <- function(regr.cv, data, outcome, folds) {
  pred.cv <- NULL
  for (fold in 1:length(folds)) {
    test.idx <- folds[[fold]]
    pred.cv[[fold]] <- data.frame(obs=outcome[test.idx],
                                  pred=predict(regr.cv[[fold]], newdata=data,
                                               type="response")[test.idx])
  }
  return(pred.cv)
}

# perform 10-folds cross-validation for model M2
# set the random seed to 1 before creating the folds
set.seed(1)
folds <- createFolds(infert$case, k=10)

# create the set of fitted models using cv.M2
cv.M2 <- glm.cv(case ~ age + parity + spontaneous, infert, folds)

# Use predict.cv() function to create a list of dataframes containing the 
# observed and predicted outcomes for the fitted models
pred.M2 <- predict.cv(cv.M2, infert, infert$case, folds)

# Use logkik.binom() to compute the log-likelihood of the predicted probabilites 
# for each test fold.
# Create an empty vector to store the log-likelihoods in
loglik.per.fold <- NULL

for(i in 1:length(pred.M2)) {
  y <- pred.M2[[i]]
  loglik.per.fold[i] <- loglik.binom(y$obs, y$pred)
}

# Calculate the sum of the test log-likelihoods over all folds
loglik.sum <- sum(loglik.per.fold)
loglik.sum

# The sum of the test log-likelihoods over all the folds is -144.5131


