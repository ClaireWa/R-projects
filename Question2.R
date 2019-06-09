library(readr)
library(data.table)

# Part A
# load in the files provided on learn as data tables
l1.dt <- data.table(read_csv("Biomedical Data Science/assessment1/ltegfr1.csv"))
l2.dt <- data.table(read_csv("Biomedical Data Science/assessment1/ltegfr2.csv"))

# check intersection between identifiers
length(intersect(l1.dt$id, l2.dt$ID))

# check patient records that are in l2.dt and not in l1.dt
l2.dt[!ID %in% l1.dt$id]$ID

# check patient records that are in l1.dt and not in l2.dt
l1.dt[!id %in% l2.dt$ID]$id

# merge the files into one data table, keeping all observations from both datasets, 
# ordered by identifier and follow-up time
l.dt <- merge(l1.dt, l2.dt, by.x = c("id", "fu.years"), by.y = c("ID", "fu.years"), 
              all = TRUE, sort = TRUE)

# assertion to check ordering is correct
# check order of patient ids first
stopifnot(l.dt$id[order(l.dt$id)] == l.dt$id)

# check order of follow-up times for each patient id second
for(i in 1:length(unique(l.dt$id))) {
  times <- l.dt$fu.years[l.dt$id == unique(l.dt$id)[i]]
  order <- order(times)
  stopifnot(times[order] == times)
}

####################################################################################
# Part B
# compute average eGFR and length of follow-up for each patient

ids.num <- unique(l.dt$id)
egfr.avg <- NULL
fup.avg <- NULL

# calculation
for(i in 1:length(ids.num)) {
  egfr.avg[i] <- mean(l.dt$egfr[l.dt$id==ids.num[i]], na.rm=TRUE)
  fup.avg[i] <- mean(l.dt$fu.years[l.dt$id==ids.num[i]])
}

# create table of averages
pat.avgs <- cbind(ids.num, egfr.avg, fup.avg)
pat.avgs.dt <- data.table(pat.avgs)

# tabulate number of patients with average eGFR in multiple ranges
ranges <- cut(pat.avgs.dt$egfr.avg,  c(0, 15, 30, 60, 90, Inf))
pat.avgs.split.dt <- data.table(pat.avgs.dt, ranges, stringsAsFactors = TRUE)
ranges <- data.table("Ranges" = levels(pat.avgs.split.dt$ranges), "Number of Patients in Range" = summary(pat.avgs.split.dt$ranges)[-6], keep.rownames = TRUE)

# count number of patients with missing average eGFR
pat.no.egfr.avg <- table(is.na(pat.avgs))["TRUE"]
pat.no.egfr.avg
# patient IDs who have no eGFR average
unique(l.dt$id[l.dt$egfr.avg == "NaN"])

# merge with l.dt
l.dt <- merge(l.dt, pat.avgs.split.dt, by.x = c("id"), by.y = c("ids.num"), 
              all = TRUE, sort = TRUE)

# assertion checks
# check order of patient ids first
stopifnot(l.dt$id[order(l.dt$id)] == l.dt$id)

# check order of follow-up times for each patient id second
for(i in 1:length(unique(l.dt$id))) {
  times <- l.dt$fu.years[l.dt$id == unique(l.dt$id)[i]]
  order <- order(times)
  stopifnot(times[order] == times)
}

# fit linear regression models
# determine patients with at least 3 eGFR readings

pat.for.reg.models <- NULL

for(i in 1:length(unique(l.dt$id))) {
  if(sum(!is.na(l.dt$egfr[l.dt$id == unique(l.dt$id)[i]])) >= 3) {
    pat.for.reg.models <- append(pat.for.reg.models, i)
  }
}

# subset l.dt to only contain data on patients with atleast 3 eGFR readings
l.reg.models.dt <- l.dt[l.dt$id %in% pat.for.reg.models, ]

# Remove observations with missing values 
l.reg.models.dt <- na.omit(l.reg.models.dt)

n <- length(pat.for.reg.models)

# fit lm for each patient with atleast 3 eGFR readings against follow-up time
lms <- lapply(1:n, function(x) lm(l.reg.models.dt$egfr[l.reg.models.dt$id == unique(l.reg.models.dt$id)[x]] ~ 
                                     l.reg.models.dt$fu.years[l.reg.models.dt$id == unique(l.reg.models.dt$id)[x]], na.action=na.exclude))

# store coefficients (intercept and slope) for each lin reg model
coeff <- sapply(lms, coef)

# create a data frame containing coefficients for each lin reg model and 
# corresponding Patient ID 
lin.reg.models.dt <- data.frame("Patient ID" = pat.for.reg.models, "Intercept" = coeff[1,], "Slope" = coeff[2,])

# count number of patients with slope in multiple ranges
# create ranges
r.1 <- sum(lin.reg.models.dt$Slope < -3)
r.2 <- sum(lin.reg.models.dt$Slope >= -3 & lin.reg.models.dt$Slope < 0)
r.3 <- sum(lin.reg.models.dt$Slope >= 0 & lin.reg.models.dt$Slope <= 3)
r.4 <- sum(lin.reg.models.dt$Slope > 3)

# combine ranges and number of patients per range in table
pat.per.slope.dt <- data.table("Range < -3" = r.1, "Range [-3, 0)" = r.2, 
                     "Range [0, 3]" = r.3, "Range > 3" = r.4)


##################################################################################

# Part C

# use patient averages from part b to create a data table of patients in range(0,15]
pat.egfr.avg.r1 <- pat.avgs.dt[pat.avgs.dt$egfr.avg >0 & pat.avgs.dt$egfr.avg<= 15]

# Store the rest of the patient details for patients in range(0,15]
ID <- NULL
Sex <- NULL
Age <- NULL 
Fup.last <- NULL 
eGFR.Measurements <- NULL
j <- 0
for (i in pat.egfr.avg.r1$ids.num) {
  j <- j + 1
  ID[j] <- unique(l.dt$id)[i]
  Sex[j] <- mean(l.dt$sex[l.dt$id == ID[j]], na.rm=TRUE)
  Age[j] <- mean(l.dt$baseline.age[l.dt$id == ID[j]], na.rm=TRUE)
  Fup.last[j] <- max(l.dt$fu.years[l.dt$id == ID[j]])  
  eGFR.Measurements[j] <- sum(l.dt$id == ID[j])
}

temp.pat.dt <- data.table(ID, Sex, Age, Fup.last, eGFR.Measurements)
pat.egfr.avg.r1 <- merge(temp.pat.dt, pat.egfr.avg.r1, by.x = "ID", by.y = "ids.num")
pat.egfr.avg.r1 <- subset(pat.egfr.avg.r1, select= -c(fup.avg))

setnames(pat.egfr.avg.r1, old = c("Age", "Fup.last", "eGFR.Measurements", "egfr.avg"), 
         new = c("Age at baseline", "Time of last eGFR reading", "Number of eGFR measurements taken", 
                 "Average eGFR"))

##################################################################################

# Part D

# Patient 3
# collect patient 3 data by subsetting l.dt
pat3.dt <- l.dt[l.dt$id ==3]

# Plot eGFR wrt time
plot(pat3.dt$fu.years, pat3.dt$egfr, main = "Patient 3", xlab = "Follow-up Time", ylab = "eGFR Measurements")

# Fit a linear regression model and add the regression line to the plot
lm3 <- lm(pat3.dt$egfr~ pat3.dt$fu.years)
abline(lm3, col = "blue")

# Compute the 95% confidence intervals
confint(lm3)

# Remove extreme values and plot the regression line again in different colour
pat3.new.dt <- pat3.dt[pat3.dt$egfr > min(pat3.dt$egfr) & pat3.dt$egfr < max(pat3.dt$egfr),]
lm3.new <- lm(pat3.new.dt$egfr ~ pat3.new.dt$fu.years)
abline(lm3.new, col = "red")
abline()

# Patient 37
# collect patient 3 data by subsetting l.dt
pat37.dt <- l.dt[l.dt$id ==37]

# Plot eGFR wrt time
plot(pat37.dt$fu.years, pat37.dt$egfr, main = "Patient 37", xlab = "Follow-up Time", ylab = "eGFR")

# Fit a linear regression model and add the regression line to the plot
lm37 <- lm(pat37.dt$egfr~ pat37.dt$fu.years)
abline(lm37, col = "blue")

# Compute the 95% confidence intervals
confint(lm37)

# Remove extreme values and plot the regression line again in different colour
pat37.new.dt <- pat37.dt[pat37.dt$egfr > min(pat37.dt$egfr) & pat37.dt$egfr < max(pat37.dt$egfr),]
lm37.new <- lm(pat37.new.dt$egfr ~ pat37.new.dt$fu.years)
abline(lm37.new, col = "red")

# Patient 162
# collect patient 3 data by subsetting l.dt
pat162.dt <- l.dt[l.dt$id ==162]

# Plot eGFR wrt time
plot(pat162.dt$fu.years, pat162.dt$egfr, main = "Patient 162", xlab = "Follow-up Time", ylab = "eGFR")

# Fit a linear regression model and add the regression line to the plot
lm162 <- lm(pat162.dt$egfr~ pat162.dt$fu.years)
abline(lm162, col = "blue")

# Compute the 95% confidence intervals
confint(lm162)

# Remove extreme values and plot the regression line again in different colour
pat162.new.dt <- pat162.dt[pat162.dt$egfr > min(pat162.dt$egfr) & pat162.dt$egfr < max(pat162.dt$egfr),]
lm162.new <- lm(pat162.new.dt$egfr ~ pat162.new.dt$fu.years)
abline(lm162.new, col = "red")

# Patient 223
# collect patient 3 data by subsetting l.dt
pat223.dt <- l.dt[l.dt$id ==223]
pat223.dt <- na.omit(pat223.dt)

# Plot eGFR wrt time
plot(pat223.dt$fu.years, pat223.dt$egfr, main = "Patient 223", xlab = "Follow-up Time", ylab = "eGFR")

# Fit a linear regression model and add the regression line to the plot
lm223 <- lm(pat223.dt$egfr~ pat223.dt$fu.years)
abline(lm223, col = "blue")

# Compute the 95% confidence intervals
confint(lm223)

# Remove extreme values and plot the regression line again in different colour
pat223.new.dt <- pat223.dt[pat223.dt$egfr > min(pat223.dt$egfr) & pat223.dt$egfr < max(pat223.dt$egfr),]
lm223.new <- lm(pat223.new.dt$egfr ~ pat223.new.dt$fu.years)
abline(lm223.new, col = "red")


