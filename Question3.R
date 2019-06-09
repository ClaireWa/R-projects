library(readr)
library(data.table)

# Part A

egfr.mdrd4 <- function(scr, age, sex, ethnic) {
  
  # INPUT:
  #     scr: A vector containing serum creatinine data as positive real value(s)
  #     age: A vector containing age data as positive real value(s)
  #     sex: A factor variable with levels "Male" and "Female"
  #     ethnic: A factor variable with levels "Black" and "Other"
  #
  # OUTPUT: 
  #     A vector which returns the estimated eGFR obtained by using the MDRD4 equation
  
  sex.coeff <- (sex=="Female")*0.742 + (sex=="Male")
  ethnic.coeff <- (ethnic=="Black")*1.212 + (ethnic=="Other")
  x <- 175*((scr/88.42)^-1.154)*(age^-0.203)*sex.coeff*ethnic.coeff
  
  return(x)
}

egfr.ckdepi <- function(scr, age, sex, ethnic) {
  
  # INPUT:
  #     scr: A vector containing serum creatinine data as positive real value(s)
  #     age: A vector containing age data as positive real value(s)
  #     sex: A factor variable with levels "Male" and "Female"
  #     ethnic: A factor variable with levels "Black" and "Other"
  #
  # OUTPUT: 
  #     A vector which returns the estimated eGFR obtained by using the CKD-EPI equation
  
  k <- (sex=="Female")*0.7 + (sex=="Male")*0.9
  alpha <- (sex=="Female")*-0.329 + (sex=="Male")*-0.411 
  min <- (pmin((scr/88.42)/k, 1))^alpha
  max <- (pmax((scr/88.42)/k, 1))^-1.209
  sex.coeff <- (sex=="Female")*1.018 + (sex=="Male")
  ethnic.coeff <- (ethnic=="Black")*1.159 + (ethnic=="Other")
  x <- 141*min*max*(0.993^age)*sex.coeff*ethnic.coeff
  
  return(x)
}

#####################################################################################

# Part B

# import scr2.csv dataset 
scr2.dt <- data.table(read_csv("Biomedical Data Science/assessment1/scr2.csv"), stringsAsFactors = TRUE)

# check for missing data
table(is.na(scr2.dt))

# Remove the missing data
scr2.dt <- na.omit(scr2.dt)

# Compute the eGFR according to the two equations
mdrd4 <- egfr.mdrd4(scr2.dt$scr, scr2.dt$age, scr2.dt$sex, scr2.dt$ethnic)
ckdepi <- egfr.ckdepi(scr2.dt$scr, scr2.dt$age, scr2.dt$sex, scr2.dt$ethnic)

# Calculate the mean and standard deviation of mdrd4 and ckdepi rounded to two decimal places
m.mdrd4 <- round(mean(mdrd4), 2)
m.ckdepi <- round(mean(ckdepi), 2)

sd.mdrd4 <- round(sd(mdrd4), 2)
sd.ckdepi <- round(sd(ckdepi), 2)

# Calculate their Pearson correlation coefficient
Pearson.CorrCoef <- cor(mdrd4, ckdepi)

# Report same quantities according to strata
# split mdrd4 into strata as per instructions and create a data table of mdrd4 estimates and the range the fall in
mdrd4.split <- cut(mdrd4,  c(0, 60, 90, Inf))
mdrd4.dt <- data.table(mdrd4, mdrd4.split, stringsAsFactors = TRUE)

# split data table into ranges
mr1.dt <- mdrd4.dt[mdrd4.dt$mdrd4.split == "(0,60]"]
mr2.dt <- mdrd4.dt[mdrd4.dt$mdrd4.split == "(60,90]"]
mr3.dt <- mdrd4.dt[mdrd4.dt$mdrd4.split == "(90,Inf]"]

# Calculate the means and standard deviations
m.mdrd4.r1 <- round(mean(mr1.dt$mdrd4), 2)
m.mdrd4.r2 <- round(mean(mr2.dt$mdrd4), 2)
m.mdrd4.r3 <- round(mean(mr3.dt$mdrd4), 2)

sd.mdrd4.r1 <- round(sd(mr1.dt$mdrd4), 2)
sd.mdrd4.r2 <- round(sd(mr2.dt$mdrd4), 2)
sd.mdrd4.r3 <- round(sd(mr3.dt$mdrd4), 2)


# Report same quantities according to strata
# split ckdepi into strata as per instructions and create a data table of mdrd4 estimates and the range the fall in
ckdepi.split <- cut(ckdepi,  c(0, 60, 90, Inf))
ckdepi.dt <- data.table(ckdepi, ckdepi.split, stringsAsFactors = TRUE)

# split data table into ranges
cr1.dt <- ckdepi.dt[ckdepi.dt$ckdepi.split == "(0,60]"]
cr2.dt <- ckdepi.dt[ckdepi.dt$ckdepi.split == "(60,90]"]
cr3.dt <- ckdepi.dt[ckdepi.dt$ckdepi.split == "(90,Inf]"]

# Calculate the means and standard deviations
m.ckdepi.r1 <- round(mean(cr1.dt$ckdepi), 2)
m.ckdepi.r2 <- round(mean(cr2.dt$ckdepi), 2)
m.ckdepi.r3 <- round(mean(cr3.dt$ckdepi), 2)

sd.ckdepi.r1 <- round(sd(cr1.dt$ckdepi), 2)
sd.ckdepi.r2 <- round(sd(cr2.dt$ckdepi), 2)
sd.ckdepi.r3 <- round(sd(cr3.dt$ckdepi), 2)

# cannot calculate their Pearson correlation coefficients for each strata as dimensions do not match

####################################################################################

# part C

# Produce a scatter plot of the two eGFR vectors
plot(mdrd4, ckdepi, main="Comparison of eGFR Measurements",
     xlab="MDRD4", ylab="CKDEPI")

# Add vertical and horizontal lines corresponding to median, first and third quartiles
qt.mdrd4 <- quantile(mdrd4, probs = seq(0, 1, 0.25))
qt.ckdepi <- quantile(ckdepi, probs = seq(0, 1, 0.25))

abline(h = median(ckdepi), col = "blue")
abline(h = qt.ckdepi[2], col = "red")
abline(h = qt.ckdepi[4], col = "red")
abline(v = median(mdrd4), col = "blue")
abline(v = qt.mdrd4[2], col = "red")
abline(v = qt.mdrd4[4], col = "red")

