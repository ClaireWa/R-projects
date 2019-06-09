library(readr)
library(data.table)

#Part A

# import the two lipids files
lipids <- read.delim("lipids.txt", stringsAsFactors=FALSE)
lipid.classes <- read.delim("lipid-classes.txt", header=FALSE, stringsAsFactors=FALSE)

# create the results data table and include a column (character type)
# where the classes will be added
results.dt <- data.table(lipids, "Class" = "0")

# annotate each lipid species with its corresponding full class name
for (i in (1:nrow(lipids))) {
  for (j in (1:length(lipid.classes[,1]))) {
    if  (grepl(paste("\\b", lipid.classes[j,1], "\\b", sep=""), lipids[i,1], ignore.case = TRUE) == "TRUE") {
      results.dt[i,4] <- lipid.classes[j,2]
    }
  }
}

# check how variables are coded
str(results.dt)
# convert class column to factor
results.dt <- results.dt[, Class:= as.factor(results.dt$Class)]

# Check for Nas
table(is.na(results.dt))

# summary of lipids per class
summary(results.dt$Class)

#Ceramides            Cholesterol esters 
#13                             9 
#Diacylglycerols      Lysophosphatidylcholines 
#16                            13 
#Lysophosphatidylethanolamines          Phosphatidylcholines 
#3                            42 
#Phosphatidylethanolamines           Phosphatidylserines 
#25                             8 
#Triacylglycerols 
#147 

#####################################################################################

# Part B

# Compute the Wald test statistic for each lipid species
# Wald stat = (regressesion coefficient (beta))/se(regression coeff (beta))
# odds ratio = exp(X*beta), log(odds) = X*beta
# Include a column (numeric type) where the Wald test statistics will be added
results.dt <- data.table(results.dt, "Wald Test Statistic" = c(log(results.dt$oddsratio)/results.dt$se))

# Assuming 288 patients in dataset, derive a p-value for each lipid species using
# t distribution and append to results.dt
N <- 288
# Four predictors in each fitted model
p <- 4
# include p-value for each species in data table
results.dt<- data.table(results.dt, "P-Value" = 2*pt(-abs(results.dt$`Wald Test Statistic`), df=N-p-1))

# Repeat calculation of p-values using the normal distribution
p.val.norm <- 2*pnorm(-abs(results.dt$`Wald Test Statistic`))

# compare values
pvals.dt <- data.table("T Dist" = results.dt$`P-Value`, "Norm Dist" = p.val.norm)
pvals.dt$Difference <- pvals.dt$`T Dist` - pvals.dt$`Norm Dist`
diff.range <- paste(signif(min(pvals.dt$Difference), 3), signif(max(pvals.dt$Difference), 3), sep = ",")
plot(pvals.dt$`T Dist`, pvals.dt$`Norm Dist`, main = "T vs. Normal Distribution", xlab="T Distribution", ylab = "Normal Distribution")

####################################################################################

# Part C

holm.bonferroni <- function(results.dt, alpha) {

  # A function which implements the Holm-Bonferroni method to control the family
  # -wise error rate, it returns the subset of the input data table that are significant
  # according to this method. All hypotheses with index at most k are rejected.
  
  # Check that p-values are in numeric form
  if (!(is.numeric(results.dt$`P-Value`) || is.integer(results.dt$`P-Value`)))
    stop("The function works only for numerical p-value vectors")
  # check that provided alpha is in numeric form
  if ((!(is.numeric(alpha) || is.integer(alpha)))||(!(length(alpha)==1)))
    stop("The function works only for single numerical alpha value")

  # m is the total number of p-values which is equivalent to the total number of rows
  m <- nrow(results.dt)
  
  # Sort the results according to p-values
  results.sort <- results.dt[order(results.dt$`P-Value`),]
  pval.sort <- results.sort$`P-Value`
  
  # find the largest index k for which pval_k < alpha/(m+1-k)
  # where pval_k is the k-th p-value after sorting
  k <- 0
  for(i in 1:length(pval.sort)) {
    if(pval.sort[i] < alpha/(m+1-i)) {
      k <- i
    }
    else
      break
  }
  
  # Return the subset of results that are significant
  results.significant <- results.sort[1:k,]
  print(paste("Largest index k =", k))
  return(results.significant)
}

###################################################################################

# PART D

benjamini.hochberg <- function(results.dt, q) {
  
  # A function which implements the Benjamini-Hochberg method to control the false
  # discovery rate, it returns the subset of the input datatable that are significant
  # according to this method. All hypotheses with index at most k are rejected.
  
  # Check that p-values are in numeric form
  if (!(is.numeric(results.dt$`P-Value`) || is.integer(results.dt$`P-Value`)))
    stop("The function works only for numerical p-value vectors")
  # check that provided q is in numeric form
  if ((!(is.numeric(q) || is.integer(q)))||(!(length(q)==1)))
    stop("The function works only for single numerical q value")
  
  # m is the total number of p-values
  m <- length(col(results.dt)[,1])
  
  # Sort the results according to p-values
  results.sort <- results.dt[order(results.dt$`P-Value`),]
  pval.sort <- results.sort$`P-Value`
  
  # find the largest k for which pval_k < (k/m)*q
  # where q is the false discovery rate desired
  k <- 0
  for(i in 1:length(pval.sort)) {
    if(pval.sort[i] <= (i/m)*q) {
      k <- i
    }
  }
  
  # Return the subset of results that are significant
  results.significant <- results.sort[1:k,]
  print(paste("Largest index k =", k))
  return(results.significant)
}

#################################################################################

# Part E

# Produce a volcano plot with a different colour for lipid species that are 
# significant after controlling for the family-wise error rate given a nominal 
# significance threshold alpha = 0.05

# Get the subset of significant points
sgnf.Holm <- holm.bonferroni(results.dt, 0.05)

# Do the same procedure for the BH procedure
sgnf.Benj <- benjamini.hochberg(results.dt, 0.01)

# Note: In logisitic regression, a volcano plot should use the log-odds ratio on 
# the x-axis. Using odds ratio would produce an asymmetric plot, as odds ratios are
# bounded below by zero but are unbounded above.

# volcano plot
# choose overall symbol
symbol <- rep(0,nrow(results.dt))
# choose overall colour
c <- rep(1,nrow(results.dt))
# change colour for significant values according to Holm-Bonf method
c[(results.dt$`P-Value` < max(sgnf.Holm$`P-Value`))] <- 2

# produce volcano plot using different colours
plot(log(results.dt$oddsratio), -log10(results.dt$`P-Value`), main="Volcano Plot",
     xlab="Coefficient (Log of odds ratio)", ylab="-log10(p-value)", pch=symbol, cex=0.7,
     col=c)

# create another plot over this plot to highlight signif points according to Benj-Hoch method
par(new=TRUE)
# use different symbols for significant points according to this method
symbol[(results.dt$`P-Value` < max(sgnf.Benj$`P-Value`))] <- 3
# use different colout to make it clearer
c[(results.dt$`P-Value` < max(sgnf.Benj$`P-Value`))] <- 3

plot(log(results.dt$oddsratio), -log10(results.dt$`P-Value`), main="Volcano Plot",
     xlab="Coefficient (Log of odds ratio)", ylab="-log10(p-value)", pch=symbol, cex=0.7,
     col=c)

# add legend
legend(x=0.293, y=1.2, legend=c('Non-significant', 'Holm-Bonferroni Signigicant', 'Benjamini-Hochberg Signigicant'),
       col=c('black', 'red', 'green'), pch=c(0,0,3),  cex = 0.55)

##################################################################################

# PART F

# Get the significant subsets for each procedure
sgnf.holm <- holm.bonferroni(results.dt, 0.05)
sgnf.benj <- benjamini.hochberg(results.dt, 0.05)

# inner merge the subsets to get the lipid species that are considered 
# significant by eihter method
sub.signif <- merge(sgnf.holm, sgnf.benj)

# Now we can see the significant subset of lipid scpecies
sub.signif$lipid.species

