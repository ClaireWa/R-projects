library(readr)
library(data.table)

# Part A
# Read in file and store as data table
gdm.dt <- read.delim("GDM.raw.txt", header = TRUE)
gdm.dt <- data.table(gdm.dt)

# Store rsID and coded allele in separate data table
# get list of column names
x <- colnames(gdm.dt)
# split SNP names in x into rsID and reference allele
y <- strsplit(x[4:length(x)],"_")
z <- unlist(y)
rsID <- z[ c(TRUE,FALSE) ]
ref.allele <- z[ c(FALSE,TRUE) ]

# Create a data table snp.allele.dt which stores rsID and the coded allele
snp.allele.dt <- as.data.table(cbind(rsID, ref.allele))

# Check for missing values in gdm.dt
table(is.na(gdm.dt))

# Impute missing values according to SNP-wise average allele count
# Loop over columns with missing values
# sex shows as false, loop will update these, so ignore this column and other columns
# that are not SNPs so to not change their values
for (col in colnames(gdm.dt[,4:length(gdm.dt)])) {
  gdm.dt[[col]][is.na(gdm.dt[[col]])] <- mean(gdm.dt[[col]], na.rm = T)
}

# stop if any missing values remain
stopifnot(!is.na(gdm.dt))

##############################################################################################

# Part B

# Write a function that fits a logistic regression model for each SNP in x, where x
# is a data table of SNPs, y is a binary oytcine vector and order is either
# True or False
univ.glm.test <- function(x, y, order = FALSE) {
  # check vectors have same number of elements
  if (!(nrow(x) == length(y)))
    stop("The vector lengths in x do not match the vector length of y")
  # check y is numeric and binary vector
  if (!(is.numeric(y)||is.integer(y))) 
    stop("Only numeric vectors are accepted.")
  vSet = unique(y)
  if (length(vSet) > 2)
    stop("y must be a binary outcome vector")
  
  
  crude <- NULL
  
  for(snp in colnames(x)) {
    fit <- glm(y ~ x[[snp]], family = "binomial")
    output <- as.data.table(coef(summary(fit)))[-1,-3]
    crude <- rbind(crude, output)
  }
  colnames(crude) <- c("beta", "std.error", "p.value")
  # Add the column of SNP names and odds ratios
  crude <- cbind(crude, odds.ratio = exp(crude$beta))
  crude <- cbind(snp=colnames(x), crude)
  
  # If order is set to TRUE, output is ordered by increasing p-value
  if(order==TRUE){
    crude <- crude[order(crude$p.value), ]
  }
  
  # Get rid of the confusing row names
  rownames(crude) <- c()
  return(crude)
}

##############################################################################################

# PART C

# Use the function from part B to run an association study for all SNPs in gdm.dt against
# having gestational diabetes ("pheno")

# The set of snps is from the 4th to 179th column in gdm.dt, so exclude the first 3 columns
fits.dt <- univ.glm.test(gdm.dt[,4:length(gdm.dt)], gdm.dt$pheno)

# SNP which is most strongly associated with "pheno" will have smallest p-value
signif.ID <- which(min(fits.dt$p.value)==fits.dt$p.value)

# Summary statistics for the most strongly associate SNP
fits.dt[signif.ID,]

# 95 and 99 percent CI on the odds ratio
#95%
exp(cbind(fits.dt$beta[signif.ID]- 1.96*fits.dt$std.error[signif.ID], fits.dt$beta[signif.ID] + 1.96*fits.dt$std.error[signif.ID]))
#99%
exp(cbind(fits.dt$beta[signif.ID]- 2.58*fits.dt$std.error[signif.ID], fits.dt$beta[signif.ID] + 2.58*fits.dt$std.error[signif.ID]))

# SNP which most significant protective effect will have largest p-value
signif.ID2 <- which(max(fits.dt$p.value)==fits.dt$p.value)

# Summary statistics for the most strongly associate SNP
fits.dt[signif.ID2,]

# 95 and 99 percent CI on the odds ratio
#95%
exp(cbind(fits.dt$beta[signif.ID2]- 1.96*fits.dt$std.error[signif.ID2], fits.dt$beta[signif.ID2] + 1.96*fits.dt$std.error[signif.ID2]))
#99%
exp(cbind(fits.dt$beta[signif.ID2]- 2.58*fits.dt$std.error[signif.ID2], fits.dt$beta[signif.ID2] + 2.58*fits.dt$std.error[signif.ID2]))

#############################################################################################
# Part D

gdm.annot <- read.delim("GDM.annot.txt", header = TRUE)
gdm.annot.dt <- data.table(gdm.annot)

fits.dt$rsID <- snp.allele.dt$rsID

fits.dt <- merge(fits.dt, gdm.annot.dt, by.x = "rsID", by.y = "snp")

# Report SNPs that have a p.value less than 10^-4
hit <- fits.dt[fits.dt$p.value<1e-4,]

# Report genes that are within a 1MB window of each hit snp position
unique(fits.dt[fits.dt$pos > (hit$pos[1]-1024^2) & fits.dt$pos < (hit$pos[1] + 1024^2),]$gene)
unique(fits.dt[fits.dt$pos > (hit$pos[2]-1024^2) & fits.dt$pos < (hit$pos[2] + 1024^2),]$gene)

##############################################################################################
# Part E

#subset of snps with p value < 10^-4
snps.grs1 <- subset(fits.dt, p.value < 1e-4)
gdm.grs1 <- gdm.dt[, .SD, .SDcols = snps.grs1[p.value < 1e-4]$snp]

# Ensure ordering is respected
stopifnot(colnames(gdm.grs1) == snps.grs1$snp)

weighted.score1 <- as.matrix(gdm.grs1) %*% snps.grs1$beta

#subset of snps with p value < 10^-3
snps.grs2 <- subset(fits.dt, p.value < 1e-3)
gdm.grs2 <- gdm.dt[, .SD, .SDcols = snps.grs2[p.value < 1e-3]$snp]

# Ensure ordering is respected
stopifnot(colnames(gdm.grs2) == snps.grs2$snp)

weighted.score2 <- as.matrix(gdm.grs2) %*% snps.grs2$beta


#subset of snps on the FTO gene
snps.grs3 <- subset(fits.dt, gene == "FTO")
gdm.grs3 <- gdm.dt[, .SD, .SDcols = snps.grs3[gene == "FTO"]$snp]

# Ensure ordering is respected
stopifnot(colnames(gdm.grs3) == snps.grs3$snp)

weighted.score3 <- as.matrix(gdm.grs3) %*% snps.grs3$beta

# Add the three scores to the gdm.dt dataframe
gdm.dt$score1 <- weighted.score1
gdm.dt$score2 <- weighted.score2
gdm.dt$score3 <- weighted.score3


# Fits each score in a separate logistic regression model

fit.score1 <- glm(pheno ~ score1, data = gdm.dt, family = "binomial")
Odds1 <- exp(coef(summary(fit.score1))[2,1])
confint(fit.score1)
p.val1 <- coef(summary(fit.score1))[2,4]

fit.score2 <- glm(pheno ~ score2, data = gdm.dt, family = "binomial")
Odds2 <- exp(coef(summary(fit.score2))[2,1])
confint(fit.score2)
p.val2 <- coef(summary(fit.score2))[2,4]

fit.score3 <- glm(pheno ~ score3, data = gdm.dt, family = "binomial")
Odds3 <- exp(coef(summary(fit.score3))[2,1])
confint(fit.score3)
p.val3 <- coef(summary(fit.score3))[2,4]

#############################################################################################
# PART F

gdm.test <- read.delim("GDM.test.txt", header = TRUE)

# GRS including snps with p value < 10^-4
gdm.test.grs1 <- gdm.test[, colnames(gdm.test) %in% snps.grs1$rsID]

# Ensure ordering is respected
stopifnot(colnames(gdm.test.grs1) == snps.grs1$rsID)

weighted.score.test1 <- as.matrix(gdm.test.grs1) %*% snps.grs1$beta

# GRS including  snps with p value < 10^-3
gdm.test.grs2 <- gdm.test[, colnames(gdm.test) %in% snps.grs2$rsID]

# Order correctly
gdm.test.grs2 <- gdm.test.grs2[, order(colnames(gdm.test.grs2))]

# Ensure ordering is respected
stopifnot(colnames(gdm.test.grs2) == snps.grs2$rsID)

weighted.score.test2 <- as.matrix(gdm.test.grs2) %*% snps.grs2$beta


# GRS including snps on the FTO gene
gdm.test.grs3 <- gdm.test[, colnames(gdm.test) %in% snps.grs3$rsID]

gdm.test.grs3 <- gdm.test.grs3[, order(colnames(gdm.test.grs3))]

# Ensure ordering is respected
stopifnot(colnames(gdm.test.grs3) == snps.grs3$rsID)

weighted.score.test3 <- as.matrix(gdm.test.grs3) %*% snps.grs3$beta

# Add the three scores to the gdm.test dataframe
gdm.test$score1 <- weighted.score.test1
gdm.test$score2 <- weighted.score.test2
gdm.test$score3 <- weighted.score.test3

gdm.test <- data.table(gdm.test)

#############################################################################################
# PART G
# Use the log regression models fitted in E to predict the outcome of patients in gdm.test
predict.gdm1 <- predict(fit.score1, gdm.test, type = "response")
predict.gdm2 <- predict(fit.score2, gdm.test, type = "response")
predict.gdm3 <- predict(fit.score3, gdm.test, type = "response")

# Implement a function to compute the binomial loglikelihood

loglik.binom <- function(y.obs, y.pred) {
  
  # INPUT
  #    y.obs: The vector of observed outcomes that we have from our data (0 or 1)
  #   y.pred: The vector of fitted probabilities that our model provides
  #
  
  # locate the controls in y.obs and obtain the corresponding probabilities in y.pred
  ctrl.idx <- y.obs == 0
  ctrl.prob <- y.pred[ctrl.idx]
  L1 <- sum(log(1 - ctrl.prob))
  
  case.idx <- y.obs == 1
  case.prob <- y.pred[case.idx]
  L2 <- sum(log(case.prob))
  
  L <- L1 + L2
  
  # OUTPUT: The binomial log-likelihood for the outcome variable
  
  return(L)
}
# Compute the test log-likelihoods for the predicted probs from the three risk score models
loglik.binom(gdm.test$pheno, predict.gdm1)
loglik.binom(gdm.test$pheno, predict.gdm2)
loglik.binom(gdm.test$pheno, predict.gdm3)

##############################################################################################
# Part H
# Read in file (contains summary stats from a different study on same set of SNPs)
gdm.study <- read.delim("GDM.study2.txt", header = TRUE, stringsAsFactors = FALSE)

# Perform a meta-analysis with the results obtained in part C
fits.dt <- univ.glm.test(gdm.dt[,4:ncol(gdm.dt)], gdm.dt$pheno)
# remove score rows
fits.dt <- fits.dt[-c(177:179)]

# split snp names into ID and allele to match the format in gdm.study
y <- strsplit(fits.dt$snp,"_")
z <- unlist(y)
rsID <- z[ c(TRUE,FALSE) ]
ref.allele <- z[ c(FALSE,TRUE) ]

fits.dt$snp <- rsID
fits.dt$ref.allele <- ref.allele

gwas1 <- fits.dt
gwas2 <- data.table(gdm.study)

# harmonize the two datasets
# Check which SNP IDs are in both sets
snps.in.both <- intersect(gwas1$snp, gwas2$snp)  
# They contain all the same snps

# To ensure ordering is the same, order by snp ID
gwas1 <- gwas1[order(gwas1$snp), ]
gwas2 <- gwas2[order(gwas2$snp), ]
stopifnot(all.equal(gwas1$snp, gwas2$snp))

# order snp.allele.dt so that it matches gwas1
snp.allele <- snp.allele.dt[order(snp.allele.dt$rsID), ]

# check which snps have alleles matching
matches <- (snp.allele$ref.allele == gwas2$effect.allele)

# Check which SNPs have the other allele matching
flipped <- (snp.allele$ref.allele == gwas2$other.allele)

# Create a table to see which alleles do and do not match
table(matches, flipped) 

# Change the sign for of the coefficient for the SNPs that have flipped alleles
# before entering the meta-analysis
gwas2$beta[flipped] <- -gwas2$beta[flipped]

# The alleles which do not match cannot be meta-analysed
# get subset of matching alleles for meta-analysis
all.match <- matches != flipped

gwas.match1 <- gwas1[all.match,]
gwas.match2 <- gwas2[all.match,]

beta1 <- gwas.match1$beta
beta2 <- gwas.match2$beta

# using inverse variance weighting to perform fixed effect meta-analysis
weight.gwas1 <- 1/gwas.match1$std.error^2
weight.gwas2 <- 1/gwas.match2$se^2

# By looking at the weights assigned to the two studies, it appears that in general the 
# second study is less powered.
head(weight.gwas1)
head(weight.gwas2)

# The compuation of the meta-analysis effect size is a weighted sum of the effect sizes 
# from each study, weighted according to the weight just derived.
beta.ma <- (weight.gwas1*beta1 + weight.gwas2*beta2)/(weight.gwas1 + weight.gwas2)
se.ma <- sqrt(1 / (weight.gwas1 + weight.gwas2))

# Note that by combining the two studies through a meta-analysis we have increased the power.

# calculate meta-analysis p-values
pval.ma <- 2*(pnorm(abs(beta.ma / se.ma), lower.tail = FALSE))
# get the set of SNP IDs that have matching alleles
rsID.ma <- as.character(snp.allele$rsID[all.match])

meta.summary <- as.data.frame(cbind(beta.ma, se.ma, pval.ma))
meta.summary$rsID <- rsID.ma

# Find the SNPs with meta-analysis p-value < 1e-4
signif.snps <- subset(meta.summary, pval.ma < 1e-4)
# Sort by increasing p-value
signif.snps <- signif.snps[order(signif.snps$pval.ma),]

# Summary of meta-analysis
signif.snps
