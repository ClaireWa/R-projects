library(readr)
library(data.table)
library(tibble)
library(car)
library(lsr)

num.dat <- data.table(read.csv("numerics-data.csv"))

# clean the data
# remove respondents who answered check questions incorrectly
# sebis.Check question should be answered "always" which corresponds to 5
# get subset of data table where sebis.Check answered corretly
data.dt <- num.dat[num.dat$sebisLikert.sebisCheck. == 5]
# 4 obs removed

# hackersLikert.folkHCheck. should be answered "disagree" which corresponds to 2
# get subset of data table where hackersLikert.folkHCheck. answered corretly
data.dt <- num.dat[num.dat$hackersLikert.folkHCheck. == 2]
# 3 obs removed
# 7 obs removed in total

# check for NAs in data
table(is.na(data.dt))

# create westin column
data.dt <- add_column(data.dt, westin = c("0"), .after = 2)

# sort each participant into privacy view group in westin column
for (i in (1:nrow(data.dt))) {
  if (((data.dt$segmentation.westinC.[i] == 1)||(data.dt$segmentation.westinC.[i] == 2))
       && ((data.dt$segmentation.westinM.[i] == 3)||(data.dt$segmentation.westinM.[i] == 4))
       &&((data.dt$segmentation.westinE.[i] == 3)||(data.dt$segmentation.westinE.[i] == 4))){
    data.dt$westin[i] <- "Fundamentalist"
  } else if (((data.dt$segmentation.westinC.[i] == 3)||(data.dt$segmentation.westinC.[i] == 4))
      && ((data.dt$segmentation.westinM.[i] == 1)||(data.dt$segmentation.westinM.[i] == 2))
      &&((data.dt$segmentation.westinE.[i] == 1)||(data.dt$segmentation.westinE.[i] == 2))){
    data.dt$westin[i] <- "Unconcerned"
  } else {
    data.dt$westin[i] <- "Pragmatist"
  }
}

# check if westin is in factor format
is.factor(data.dt$westin)

# change westin to factor and define ordering
data.dt$westin <- as.factor(data.dt$westin)
data.dt$westin <- ordered(data.dt$westin, levels = c("Fundamentalist", "Unconcerned", "Pragmatist"))

# Calculate the SeBIS subscales
# Password Generation
# Adjust for inverted questions
data.dt <- add_column(data.dt, sebisLikert.sebisF12.inverted = 6-data.dt$sebisLikert.sebisF12., .after = "sebisLikert.sebisF12.")
data.dt <- add_column(data.dt, sebisLikert.sebisF14.inverted = 6-data.dt$sebisLikert.sebisF14., .after = "sebisLikert.sebisF14.")

# calculate the password average
password <- (data.dt$sebisLikert.sebisF12.inverted + data.dt$sebisLikert.sebisF13.
                     + data.dt$sebisLikert.sebisF14.inverted + data.dt$sebisLikert.sebisF15.)/4
data.dt <- add_column(data.dt,password, .after = 3)

# Device Securement has no inverted scales
device <- (data.dt$sebisLikert.sebisF04.+ data.dt$sebisLikert.sebisF06.
             + data.dt$sebisLikert.sebisF03. + data.dt$sebisLikert.sebisF05.)/4
data.dt <- add_column(data.dt, device, .after = 4)

# Proactive Awareness
# F8, F11, F16 abd F7 are reverse scored questions
data.dt <- add_column(data.dt, sebisLikert.sebisF08.inverted = 6-data.dt$sebisLikert.sebisF08., .after = "sebisLikert.sebisF08.")
data.dt <- add_column(data.dt, sebisLikert.sebisF11.inverted = 6-data.dt$sebisLikert.sebisF11., .after = "sebisLikert.sebisF11.")
data.dt <- add_column(data.dt, sebisLikert.sebisF16.inverted = 6-data.dt$sebisLikert.sebisF16., .after = "sebisLikert.sebisF16.")
data.dt <- add_column(data.dt, sebisLikert.sebisF07.inverted = 6-data.dt$sebisLikert.sebisF07., .after = "sebisLikert.sebisF07.")

awareness <- (data.dt$sebisLikert.sebisF08.inverted + data.dt$sebisLikert.sebisF11.inverted +
                data.dt$sebisLikert.sebisF16.inverted + data.dt$sebisLikert.sebisF10. + data.dt$sebisLikert.sebisF07.inverted)/5
data.dt <- add_column(data.dt, awareness, .after = 5)

# Updating
# no inverted scales
# round numbers to match number of decimal places in other sebis scales
updating <- round(((data.dt$sebisLikert.sebisF01.+ data.dt$sebisLikert.sebisF02.
           + data.dt$sebisLikert.sebisF09.)/3), 2)
data.dt <- add_column(data.dt, updating, .after = 6)

w <- data.dt$westin
# password Descriptive Stats
pass <- data.dt$password
# look at the means, standard deviations and count for each group
tapply(pass, w, mean) 
tapply(pass, w, sd)
tapply(pass, w, length) 

m1 <- pass ~ w
# boxplot
boxplot(m1, main = "Password by Westin Category", xlab = "Privacy View", ylab = "Password")

# device Descriptive Stats
dev <- data.dt$device
# look at the means, standard deviations and count for each group
tapply(dev, w, mean) 
tapply(dev, w, sd)
tapply(dev, w, length) 

m2 <- dev ~ w
# boxplot
boxplot(m2, main = "Device by Westin Category", xlab = "Privacy View", ylab = "Device")

# awareness Descriptive Stats
aware <- data.dt$awareness
# look at the means, standard deviations and count for each group
tapply(aware, w, mean) 
tapply(aware, w, sd)
tapply(aware, w, length) 

m3 <- aware~w
# boxplot
boxplot(m3, main = "Awareness by Westin Category", xlab = "Privacy View", ylab = "Awareness")

# updating Descriptive Stats
up <- data.dt$updating
# look at the means, standard deviations and count for each group
tapply(up, w, mean)
tapply(up, w, sd)
tapply(up, w, length)

m4 <- up ~ w
# boxplot
boxplot(m4, main = "Updating by Westin Category", xlab = "Privacy View", ylab = "Updating")

# ANOVA tests
# aov assumes equal variances and provides more info
a.pass<- aov(m1)
a.dev <- aov(m2)
a.aware<- aov(m3)
a.up <- aov(m4)
summary(a.pass)
summary(a.dev)
summary(a.aware)
summary(a.up)

# oneway.test assumes homeogeneity of variance violated
# Welch's correction applied 
a.dev2 <- oneway.test(m2)
a.dev2

# none report statistially significant results

# effect sizes
lsr::etaSquared(a.pass)
lsr::etaSquared(a.dev)
lsr::etaSquared(a.aware)
lsr::etaSquared(a.up)

# Test Assumptions
# Homogenity of Variance
# Levene's test is less sensitive to departures from normal distribution than Bartlett's test
leveneTest(m1)
leveneTest(m2)
leveneTest(m3)
leveneTest(m4)
# non-significant result, no evidence to suggest that
# variances cannot be assumed to be equal
# Test is sensitive to unequal samples however so check model plots

# Model Checking Plots
plot(a.pass,1, main = "Password")
plot(a.dev,1, main = "Device")
plot(a.aware,1, main = "Awareness")
plot(a.up,1, main = "Updating")

# Normality check
# Run Shapiro-Wilk test on ANOVA residuals
# again sensitive to unequal group samples so rely on plots
shapiro.test(residuals(object = a.pass))
shapiro.test(residuals(object = a.dev))
shapiro.test(residuals(object = a.aware))
shapiro.test(residuals(object = a.up))
# Q-Q plots
plot(a.pass,2, main = "Password")
plot(a.dev,2, main = "Device")
plot(a.aware,2, main = "Awareness")
plot(a.up,2, main = "Updating")

#Levene's Test and Shapiro-Wilk testing are considered here as alternatives to each plotting 
#methods but as they are sensitive to unequal sample sizes and small sample sizes it is decided 
#not to use them for checking the homogeneity of variance and normality assumptions 
#respectively.

