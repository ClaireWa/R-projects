library(readr)
library(data.table)

# Part A
# load in the three files provided on learn as data tables
cohort <- data.table(read_csv("cohort.csv"))
lab1 <- data.table(read_csv("lab1.csv"))
linker <- data.table(read_csv("linker.csv"))

# merge the three files into one data table
cohort.dt <- merge(lab1, linker, by.x="LABID", by.y="LABID", all= TRUE)
cohort.dt <- merge(cohort, cohort.dt, by.x=c("id", "yob"), by.y=c("id", "yob"), all= TRUE)

# re-order original cohort so that the order matches that in cohort.dt
# so that they can be comparted using assertion checks
cohort <- cohort[order(cohort$id)] 

# check if the columns in original cohort file are in cohort.dt
stopifnot(all(colnames(cohort) %in% colnames(cohort.dt)))

# check if the columns in original lab1 file are in cohort.dt
stopifnot(all(colnames(lab1) %in% colnames(cohort.dt)))

# check if all ids from original cohort file match cohort.dt
stopifnot(cohort$id == cohort.dt$id)

# check if all data from original cohort file is in cohort.dt
stopifnot(all(cohort %in% cohort.dt))

# check if yob fields match 
stopifnot(cohort.dt$yob.x == cohort.dt$yob.y)

# check number of rows match
stopifnot(nrow(cohort) == nrow(cohort.dt))
stopifnot(nrow(lab1) == nrow(cohort.dt))


# Drop LABID column from cohort.dt
cohort.dt <- subset(cohort.dt, select = -c(LABID))

# Convert albumin field to factor and order as 1=normo, 2=micro, 3=macro
cohort.dt$albumin <- as.factor(cohort.dt$albumin)
cohort.dt$albumin<- relevel(cohort.dt$albumin, "micro") #re-orders micro = 1
cohort.dt$albumin<- relevel(cohort.dt$albumin, "normo") # re-orders normo = 1 and micro pushed to level 2

#########################################################################

# Part B
# convert yob to integer to calculate ages
cohort.dt$yob <- as.integer(cohort.dt$yob)

# Count the number of missing ages
table(is.na(cohort.dt$age))[2]

# Based on age/year offset being consistent for NA ages, calculate missing ages
cohort.dt$age[which(is.na(cohort.dt$age))] <- 2019 - cohort.dt$yob[which(is.na(cohort.dt$age))[1]]

##########################################################################

# Part C
# convert diabetes status to categorical variable
cohort.dt$diabetes <- as.factor(cohort.dt$diabetes)

# create shell summary data table
summary.df <- data.frame("Variable name" = character(), "Median (IQR) or N (%)" = character(),
                         "Total N" = numeric(), "Missing N (%)" = character(), stringsAsFactors = FALSE)

# function to impute any missing values to the mean value of vector
impute.to.mean <- function(x) {
  # check the type of objects we have been given
  if (!(is.numeric(x) || is.integer(x)))
    return(x)
  # find which values are missing
  na.idx <- is.na(x)
  # replace NAs with the mean computed over the observed values
  x[na.idx] <- mean(x, na.rm=TRUE)
  # return the vector with imputed values
  return(x)
}

# function to return median and IQR for cts variables
summary.cts <- function(x, impute.to.mean = FALSE) {
  # check impute to mean 
  if (!(impute.to.mean == FALSE))
    x <- impute.to.mean(x)
  # compute median
  med <- round(median(x, na.rm = impute.to.mean), 2)
  # compute IQR
  IQR <- stats::quantile(x, probs = c(1, 3) / 4, na.rm = impute.to.mean)
  
  ans <-  paste0(med, " (", round(IQR[1L], 2), ", ", round(IQR[2L], 2), ")")
  return(ans)
}

add.summary.to.table.cts <- function(x, var, na.rm = FALSE){
  # count number of values (missing and not missing)
  n.missing <- sum(is.na(x))
  n.existing <- sum(!is.na(x))
  #total n
  n <- n.missing+n.existing
  # if vector contains missing values impute to mean
  if (!(na.rm == FALSE))
    x <- impute.to.mean(x)

  #
  summary.df <- rbind.data.frame(summary.df,
                      data.frame("Variable name" = deparse(substitute(var)), 
                                 "Median (IQR) or N (%)" = summary.cts(x),
                                 "Total N" = n, "Missing N (%)" = 
                                 paste0(n.missing, " (", round((n.missing/n)*100, 2), ")"), stringsAsFactors = FALSE))
  #return(summary.df)
}
# add continuous variable fields to summary.df
summary.df <- add.summary.to.table.cts(cohort.dt$yob, yob)
summary.df <- add.summary.to.table.cts(cohort.dt$age, age)
summary.df <- add.summary.to.table.cts(cohort.dt$bp, bp, na.rm = TRUE)
summary.df <- add.summary.to.table.cts(cohort.dt$urea, urea, na.rm = TRUE)
summary.df <- add.summary.to.table.cts(cohort.dt$creatinine, creatinine, na.rm = TRUE)
summary.df <- add.summary.to.table.cts(cohort.dt$glucose, glucose, na.rm = TRUE)


# convert to data table
summary.dt <- data.table(summary.df)
colnames(summary.dt) <- c("Variable name", "Median (IQR) or N (%)", "Total N", "Missing N (%)" )

##########################################################################

# Part D
# function that takes categorical variablesand returns tot number of non-missing
# rows of factor and overall percentage of field based on N values of input vector

summary.cat <- function(x, use.str) {
  # check the type of objects we have been given
  if (!(is.character(x) || is.factor(x)))
    stop("The vector must contain categorical type variables")

  # count non-missing rows
  n.existing <- summary(x)[use.str]
  
  # N values of input
  N <- length(x)
  
  # compute percentage 
  percentage <- round((n.existing/N)*100, 2)
  
  ans <-  paste0(n.existing, " (", percentage, ")")
  return(ans)
}


add.summary.to.table.cat <- function(x, var){
  # count number of values (missing and not missing)
  n.missing <- sum(is.na(x))
  n.existing <- sum(!is.na(x))
  #total n
  n <- n.missing+n.existing

  #
  summary.df <- rbind.data.frame(summary.df,
                                 data.frame("Variable name" = var, 
                                            "Median (IQR) or N (%)" = summary.cat(x, var),
                                            "Total N" = n, "Missing N (%)" = 
                                              paste0(n.missing, " (", round((n.missing/n)*100, 2), ")"), stringsAsFactors = FALSE))
}

# add categorical variable fields to summary.df
summary.df <- add.summary.to.table.cat(cohort.dt$albumin, var = "normo")
summary.df <- add.summary.to.table.cat(cohort.dt$albumin, var = "micro")
summary.df <- add.summary.to.table.cat(cohort.dt$albumin, var = "macro")
summary.df <- add.summary.to.table.cat(cohort.dt$diabetes, var = "0")
summary.df <- add.summary.to.table.cat(cohort.dt$diabetes, var = "1")


# convert to data table
summary.dt <- data.table(summary.df)
colnames(summary.dt) <- c("Variable name", "Median (IQR) or N (%)", "Total N", "Missing N (%)" )
