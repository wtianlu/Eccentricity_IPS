# *******************************************************************************************
# PROJECT: Eccentricity and VF modulations
# AUTHOR: Lulu Wang
# INSTITUTION: KU Leuven
# CONTENT:  Analyse raw accuracy data and TVA weights
# *******************************************************************************************
# Initialisation
# install.packages("rstatix")
cat("\f") # clear console (or Ctrl+L)

# Libraries
library(dplyr)
library(lme4)
library(reshape2) # wide to long for QCM
library(stats) # for aov
library(rstatix)

# *******************************************************************************************
# Repeated measures ANOVA - raw performance ----
# *******************************************************************************************

# Load data
fn_data = "/Users/Lulu/Documents/Experiments/WP2c_eccentricity/Data_proc_beh/WP2_beh_accuracyR.txt"
rawdata = data.frame(read.table(fn_data,header = TRUE, sep = ','))

# Organise data
factorcols = c("pid","ecc","vf","cond")
rawdata[factorcols] <- lapply(rawdata[factorcols],factor)
summary(rawdata)
rawdata %>%  group_by(ecc, vf, cond) %>%  get_summary_stats(perf, type = "mean_sd")
rawdata %>%  group_by(ecc, vf, cond) %>%  identify_outliers(perf)

# Calculate ANOVA
res.aov <- anova_test(data = rawdata, dv = perf, wid = pid,  within = c(ecc, vf, cond))
get_anova_table(res.aov)

# comparison between ecc at level cond
ecc.effect <- rawdata %>%  group_by(cond,vf) %>%  anova_test(dv = perf, wid = pid, within = ecc)
ecc.effect
get_anova_table(ecc.effect) %>%  filter(cond == 1)
get_anova_table(ecc.effect) %>%  filter(cond == 2)


# *******************************************************************************************
# Multiple Linear Regression  dPSC ----
# *******************************************************************************************

# Load data
fn_data = "/Volumes/LACIE SHARE/WP3/DataDerived/Group/WP3_dpsctab.txt"
rawdata = data.frame(read.table(fn_data,header = TRUE, sep = ','))

# Organise data
factorcols = c("ses","run","istestr") 
rawdata[factorcols] <- lapply(rawdata[factorcols],factor)
summary(rawdata)

# Loop over participants
for (part in c(1,2,3,4,5,6)) {
  print(sprintf("Participant %d",part))
  subdata <- rawdata  %>% filter(pid==part) %>% filter(istestr==0)
  summary(subdata)
  
  # Define model
  fit.dpsc <- lm(dpsc ~ ses + run + ses*run, data=subdata)
  print(summary(fit.dpsc))
}

# *******************************************************************************************
#  Linear Regression  TVA weights ----
# *******************************************************************************************
# Load data
fn_data = "/Volumes/LACIE SHARE/WP3/DataDerived/Group/fitted_awCu_processed.txt"
rawdata = data.frame(read.table(fn_data,header = TRUE, sep = ','))

# Loop over participants
for (part in c(1,2,3,4,5,6)) {
  print(sprintf("Participant %d",part))
  subdata <- rawdata$w_index[grepl(sprintf("sub-%d",part),rawdata$ID)]
  sessions <- c(1,2,3,4,5)
  # Define model
  fit.weights <- lm(subdata~sessions)
  print(summary(fit.weights))
}


# *******************************************************************************************
#  T-test  Functional connectivity ----
# *******************************************************************************************
# Load data
fn_data = "/Volumes/LACIE SHARE/WP3/DataDerived/Group/FCscores.csv"
rawdata = data.frame(read.table(fn_data,header = TRUE, sep = ','))
# Calculate paired t-test
result <-t.test(rawdata$IPSdiff, rawdata$AMYdiff, paired = TRUE, alternative = "two.sided")

# *******************************************************************************************
#  rmANOVA questionnaire data ----
# *******************************************************************************************
# Load data
fn_data = "/Volumes/LACIE SHARE/WP3/DataDerived/Group/AllQuestResponses.csv"
rawdata = data.frame(read.table(fn_data,header = TRUE, sep = ','))

for (qcmcat in c("Challenge","Interest","success","Anxiety")){
  subdata <- rawdata[grepl(qcmcat, names(rawdata))]
  subdata$id = c(1,2,3,4,5,6)
  subdata <- melt(subdata,id.vars=c("id"))
  subdata$variable <- factor(subdata$variable)
  subdata$id <- factor(subdata$id)
  print(subdata)

  my.aov <- aov(value ~ variable + Error(id),data = subdata)
  print(summary(my.aov))
  
}












