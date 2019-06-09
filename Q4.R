library(corrplot)
library(glmnet)
library(readr)
library(data.table)

# Part A
# read in csv file
nki <- read.csv("nki.csv", header = TRUE)

corr <- cor(nki[, sapply(nki, is.numeric)], use = "pairwise.complete")

# Identify unique pairs of variables with correlation coefficient > 0.8 in absolute value
highly.correlated <- NULL
corr.8 <- NULL
for(i in 1:(ncol(corr)-1)) {
  for(j in (i+1):ncol(corr)) {
    if(abs(corr[i,j]) > 0.8) {
      highly.correlated <- rbind(highly.correlated, c(i,j))
      corr.8 <- rbind(corr.8, corr[i,j])
    }
  }
}
# Make a subset of the column names that refer to numeric values
number.cols <- colnames(nki[, sapply(nki, is.numeric)])

# Find the variables that are highly correlated
corr.var <- cbind(number.cols[highly.correlated[,1]], number.cols[highly.correlated[,2]])

# Compute matrix of correlations between gene expression variables and display it so that a block structure
# is highlighted
corrplot(corr[3:72,3:72], order = "hclust", diag=FALSE, tl.cex = 0.3, tl.col = "black",
         title = "Correlation matrix", mar = c(0,1,2,0))

# create table of pairs that have corr coeff > 0.8 
cc.80 <- data.table(corr.var, corr.8)

################################################################################################
# Part B

# Run PCA over the columns containing gene expressions (7 to 76) so that it is possible to identify
# variable clusters

# Use only the columns containing gene expressions
nki.genes <- nki[,7:length(nki[1,])]

# Note that PCA makes sense only over numerical column, and missing values are not allowed. 
# Check both
stopifnot(!(is.na(nki.genes)))
stopifnot(!(is.numeric(nki.genes))||is.integer(nki.genes))

# Run the PCA so that it is possible to identify variable clusters
# use transpose of nki.genes
pca <- prcomp(t(nki.genes), scale. = TRUE)

plot(pca$x[,1:2], main = "Projection of predictors on first two Principal Components", cex = 0.7,
     xlab = "First Principal Component", ylab ="Second Principal Component")

# Percentage of variance explained by first two components
var.perc <- pca$sdev^2/sum(pca$sdev^2)
sum(var.perc[1:2])

# Devise a simple rule that identifies the four genes that are most different from the rest 
# There are four clear outliers, those with (PC1,PC2) < (x, -9)
outliers <- which(pca$x[,2] <= -9)
outliers

##############################################################################################
# Part C

# Run PCA again as in B, this time in order to derive a patient-wise summary of all gene expressions
# Use regular matrix instead of the transpose to complete this
pca2 <- prcomp(nki.genes, scale. = TRUE)

# Scatter plot of the projection of patients on the first 3 PCs
plot(pca2$x[,1:3], main = "Projection of predictors on the first 3 Principal Components")

# Test if these principal components (independently) are associated with the outcome
# in unadjusted log reg models 
glm1 <- glm(nki$Event ~ pca2$x[,1], family = "binomial" )
glm2 <- glm(nki$Event ~ pca2$x[,2], family = "binomial" )
glm3 <- glm(nki$Event ~ pca2$x[,3], family = "binomial" )

summary(glm1)
summary(glm2)
summary(glm3)

# in models adjusted for age and estrogen receptor
glm.adj1 <- glm(nki$Event ~ pca2$x[,1] + nki$Age + nki$EstrogenReceptor + nki$Grade , family = "binomial" )
glm.adj2 <- glm(nki$Event ~ pca2$x[,2] + nki$Age + nki$EstrogenReceptor + nki$Grade, family = "binomial" )
glm.adj3 <- glm(nki$Event ~ pca2$x[,3] + nki$Age + nki$EstrogenReceptor + nki$Grade, family = "binomial" )

summary(glm.adj1)
summary(glm.adj2)
summary(glm.adj3)

#############################################################################################
# Part D
# Fit a lasso model to predict the binary outcome using all available predictors
# Transform to a matrix as expected by the glmnet package
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

# Remove 'Event' variable in the matrix of predictors
y.nki <- nki$Event
x.nki <- prepare.glmnet(nki, ~ . - Event)

# set seed so fitted values remain constant 
set.seed(1)
# fit the lasso model
fit.lasso <- cv.glmnet(x.nki, y.nki, family = "binomial", type.measure = "auc")

# Get the optimal lambda
opt.par <- fit.lasso$lambda.min

# Repeat procedure but only penalise gene expression variables
set.seed(1)
# Penalise only the gene expression variables
# for first 6 set penalty factor to 0 and for the rest set penalty factor to 1
fit.lasso.adj <- cv.glmnet(x.nki, y.nki, family = "binomial", type.measure = "auc", penalty.factor = c(rep(0,6), rep(1,70)))

# Get the optimal lambda
opt.par2 <- fit.lasso.adj$lambda.min

# Extract the AUC corresponding to the optimal lambda for each model
fit.lasso$cvm[fit.lasso$lambda == opt.par]
fit.lasso.adj$cvm[fit.lasso.adj$lambda == opt.par2]
