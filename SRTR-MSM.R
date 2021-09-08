# TABLE OF CONTENTS
# 1. Required packages
# 2. Data loading, subsetting and cleaning
## 2.1 Load data
## 2.2 Subsetting
## 2.3 Data cleanup
# 3. Crude survival analyses
## 3.1 Kaplan-Meier
## 3.2 Competing risk
## 3.3 Multi-state models
# 4. Additional plots






# 1. REQUIRED PACKAGES
library(haven)
library(lubridate)
library(dplyr)
library(Epi)
library(ggpubr)
library(patchwork)
library(tidyr)
library(Hmisc)


########################################################################################
# 2. DATA LOADING, SUBSETTING AND CLEANING
# 2.1 LOAD DATA
df <- read_sas("~/SRTR/2020/pubsaf2006/tx_ki.sas7bdat") # Only 'tx_ki.sas7bdat' file will be used

# 2.2 SUBSETTING
# Subset relevant variables (this includes variables that would be useful for regression, but will not be used here)
df <- select(df, "CAN_AGE_AT_LISTING","CAN_ABO","CAN_ANGINA","CAN_CEREB_VASC","CAN_DIAB","CAN_DRUG_TREAT_COPD","CAN_DRUG_TREAT_HYPERTEN",
              "CAN_ETHNICITY_SRTR","CAN_GENDER","CAN_HGT_CM","CAN_INIT_STAT","CAN_LAST_ALLOC_PRA","CAN_LAST_SRTR_PEAK_PRA",
              "CAN_LAST_STAT","CAN_LISTING_DT","CAN_MALIG",    
              "CAN_PREV_TX","CAN_RACE_SRTR","CAN_REM_CD","DON_AGE","DON_ANTI_HCV",      
              "DON_CAD_DON_COD","DON_CONT_CIGARETTE","DON_CREAT","DON_ETHNICITY_SRTR","DON_EXPAND_DON_KI","DON_GENDER","DON_HGT_CM",
             "DON_HIGH_CREAT","DON_HIST_DIAB",          
              "DON_HTN","DON_NON_HR_BEAT","DON_RACE_SRTR","DON_TY","DON_WGT_KG",             
              "DONOR_ID","ORG_TY","PERS_ID","PERS_OPTN_DEATH_DT","PERS_RELIST","PERS_RETX","PERS_RETX_TRR_ID",       
              "PERS_SSA_DEATH_DT","PX_ID","REC_AGE_AT_TX","REC_BMI","REC_COLD_ISCH_TM","REC_CTR_ID","REC_DGN",                
              "REC_DIAL_DT","REC_FAIL_DT","REC_FIRST_WEEK_DIAL","REC_GRAFT_STAT","REC_HGT_CM","REC_HISTO_TX_ID",             
              "REC_MM_EQUIV_TX","REC_PRETX_DIAL","REC_RESUM_MAINT_DIAL","REC_RESUM_MAINT_DIAL_DT",
              "REC_TX_DT","REC_TX_ORG_TY","TFL_DEATH_DT","TFL_GRAFT_DT","TFL_LAFUDATE","TFL_LASTATUS","TRR_ID","TX_ID")

# Code factors as such
df[,c("CAN_ANGINA","CAN_ABO","CAN_CEREB_VASC","CAN_DIAB","CAN_DRUG_TREAT_COPD","CAN_DRUG_TREAT_HYPERTEN","CAN_ETHNICITY_SRTR","CAN_GENDER","CAN_INIT_STAT","CAN_LAST_STAT","CAN_MALIG",
       "CAN_PREV_TX","CAN_RACE_SRTR","CAN_REM_CD","DON_ANTI_HCV","DON_CAD_DON_COD","DON_CONT_CIGARETTE","DON_ETHNICITY_SRTR","DON_EXPAND_DON_KI",
       "DON_GENDER","DON_HIGH_CREAT","DON_HIST_DIAB","DON_HTN","DON_NON_HR_BEAT","DON_RACE_SRTR","DON_TY","ORG_TY","REC_DGN",
       "REC_FIRST_WEEK_DIAL","REC_GRAFT_STAT","REC_PRETX_DIAL","REC_RESUM_MAINT_DIAL","REC_TX_ORG_TY",
       "TFL_LASTATUS")] <- lapply(df[,c("CAN_ANGINA","CAN_ABO","CAN_CEREB_VASC","CAN_DIAB","CAN_DRUG_TREAT_COPD","CAN_DRUG_TREAT_HYPERTEN","CAN_ETHNICITY_SRTR","CAN_GENDER","CAN_INIT_STAT","CAN_LAST_STAT","CAN_MALIG",
                                         "CAN_PREV_TX","CAN_RACE_SRTR","CAN_REM_CD","DON_ANTI_HCV","DON_CAD_DON_COD","DON_CONT_CIGARETTE","DON_ETHNICITY_SRTR","DON_EXPAND_DON_KI",
                                         "DON_GENDER","DON_HIGH_CREAT","DON_HIST_DIAB","DON_HTN","DON_NON_HR_BEAT","DON_RACE_SRTR","DON_TY","ORG_TY","REC_DGN",
                                         "REC_FIRST_WEEK_DIAL","REC_GRAFT_STAT","REC_PRETX_DIAL","REC_RESUM_MAINT_DIAL","REC_TX_ORG_TY",
                                         "TFL_LASTATUS")], factor) 



# Subsetting
df <- subset(df, REC_AGE_AT_TX>=18) # Only adults
df <- subset(df, REC_TX_ORG_TY=="KI") # No multi-organ
start_date <- as.Date("2000-01-01") # Start date of analysis
stop_date <- as.Date("2017-12-31") # Include patients transplanted to this date (follow-up is until June 2020 regardless)
df <- df[df$REC_TX_DT>=start_date,] # Subset transplants from this date forward
df <- df[df$REC_TX_DT<=stop_date,] # Subset until stop date
df <- subset(df, CAN_PREV_TX == 0) # Only first kidneys

# There are still some duplicate patients after subsetting first kidneys
temp <- df$PERS_ID[duplicated(df$PERS_ID)]
temp2 <- df[df$PERS_ID %in% temp,]
# Some patients have 2 entries on list. Rarely, there are 2 donors for 1 recipient on 1 day with differing outcomes 
# for the different donor rows. 
# Remove the latter cases (since I don't know which is correct) from the data frame. 
df <- df[!duplicated(df[c("PERS_ID","REC_TX_DT")]),]

# More often it's because of 2nd (or 3rd) transplants that are labeled as first transplant.
# In case of mislabeled 2nd kidneys, subset to keep only the real first kidney 
temp <- df[order(df$PERS_ID,df$REC_TX_DT),]
df <- temp[!duplicated(temp[c("PERS_ID")]),]



########################################################################################
# 2.3 DATA CLEANUP

# Group multiracial, pacific islander and native american in 'MNP' or 'other' category
df$CAN_RACE_SRTR <- recode_factor(df$CAN_RACE_SRTR, 'WHITE' = "WHITE", 'BLACK' = "BLACK", 'ASIAN' = "ASIAN", 
                         'MULTI' = "MNP", 'NATIVE' = "MNP", 'PACIFIC' = "MNP")

# Create variable transplant year
df$tx_year <- lubridate::year(df$REC_TX_DT)

# Create different age group variables
df$age_group <- cut(df$REC_AGE_AT_TX, breaks = c(0, 50, 55, 60, 65, 70, 75, 80, 100))
df$age_group_2 <- cut(df$REC_AGE_AT_TX, breaks = c(0, 50, 60, 70, 100))
df$age_group_3 <- cut(df$REC_AGE_AT_TX, breaks = c(0, 55, 65, 75, 100)) # This is the one mainly used
df$age_group_4 <- cut(df$REC_AGE_AT_TX, breaks = c(0, 30, 40, 50, 60, 70, 100))
df$age_group_5 <- cut(df$REC_AGE_AT_TX, breaks = c(0, 25,30,35,40,45,50,55,60,65,70,75,80,100))


# Create dummy coded variable for RETRANSPLANT (1/0), using 'PERS_RETX' date
# (status 'R' for retransplant in variable #'TFL_LASTATUS' is unrealistically low, so will not use this - likely many deaths after 
# retransplant are the final status here, obscuring whether this person was ever retransplanted)
df$retx <- as.factor(as.character(ifelse(!is.na(df$PERS_RETX), 1, 0)))
# Create time of retransplant (retx_yrs)
df$retx_yrs <- (as.numeric(difftime(df$PERS_RETX, df$REC_TX_DT, units = "days"))/365.25)


# DEATHS AND GRAFT LOSSES
length(which(df$TFL_LASTATUS == "D")) # Died according to TFL_LASTATUS
length(which(!is.na(df$PERS_SSA_DEATH_DT))) # SSDMF deaths
length(which(!is.na(df$PERS_OPTN_DEATH_DT))) # OPTN deaths, higher number

# Discrepancies TFL death date and SSA_DEATH_DT are uncommon, but some cases where SSMDF date is years earlier or (more often) later
plot(subset(df, TFL_LASTATUS == "D")$TFL_LAFUDATE, subset(df, TFL_LASTATUS == "D")$PERS_SSA_DEATH_DT) 

# Concordance TFL and OPTN is almost perfect (in cases where there is a TFL date)
plot(subset(df, TFL_LASTATUS == "D")$TFL_LAFUDATE, subset(df, TFL_LASTATUS == "D")$PERS_OPTN_DEATH_DT)

# Use SSA date if available (should be most accurate), otherwise OPTN date
df$death_comp_dt <- ifelse(!is.na(df$PERS_SSA_DEATH_DT), df$PERS_SSA_DEATH_DT, df$PERS_OPTN_DEATH_DT)
df$death_comp_dt <- as_date(df$death_comp_dt)
length(which(!is.na(df$death_comp_dt))) # Has added about 20.000 deaths to the SS count (and 6.000 to OPTN count)

# Create variable 'years until death'
df$death_comp_yrs <- (as.numeric(difftime(df$death_comp_dt, df$REC_TX_DT, units = "days"))/365.25)

# Create variable death_comp_status (1 for dead, 0 for censored), corresponding event time is death_comp_yrs
df$death_comp_status <- as.factor(ifelse(!is.na(df$death_comp_dt), 1, 0))

# Calculate years follow up corresponding with TFL_LAFUDATE
df$tfl_lafu_yrs <- (as.numeric(difftime(df$TFL_LAFUDATE, df$REC_TX_DT, units = "days"))/365.25)

# As noted above for TFL_LASTATUS, tfl_lafu_yrs is very inaccurate for retransplant
length(which(df$retx_yrs > df$tfl_lafu_yrs)) # many patients 'retransplanted after last FU' - should be impossible
# But we need a good 'last follow up' time variable for the multi-state models. 
# So we will update the tfl_lafu_yrs with correct retx time and SSA death information
# New variable will be called 'tfl_lafu_comp_yrs' (comp for composite)
df$tfl_lafu_comp_yrs <- if_else(df$retx == 1, df$retx_yrs, df$tfl_lafu_yrs, missing = df$tfl_lafu_yrs)
df$tfl_lafu_comp_yrs <- if_else(df$death_comp_status == 1, df$death_comp_yrs, df$tfl_lafu_comp_yrs, missing = df$tfl_lafu_comp_yrs)

# And update TFL_LASTATUS with this new info (in new variable 'tfl_lafu_comp' because some info about retransplant is lost here)
# We will mostly work with 'tfl_lafu_comp' from now on
df$tfl_lafu_comp_status <- as.factor(if_else(df$retx == 1, "R", as.character(df$TFL_LASTATUS), 
                                             missing = as.character(df$TFL_LASTATUS)))
df$tfl_lafu_comp_status <- as.factor(if_else(df$death_comp_status == 1, "D", as.character(df$tfl_lafu_comp_status), 
                                             missing = as.character(df$tfl_lafu_comp_status)))


# Have look at REC_FAIL_DT (graft failure date)
summary(as.numeric(difftime(df$REC_FAIL_DT, df$REC_TX_DT, units = "days"))/365.25) # It's almost always just the transplant date
# Which makes no sense, this will not be a usable parameter

# Using TFL_GRAFT_DT (should contain same info as REC_FAIL_DT) instead
summary(as.numeric(difftime(df$TFL_GRAFT_DT, df$REC_TX_DT, units = "days"))/365.25) # This is believable
df$tfl_graft_yrs <- as.numeric(difftime(df$TFL_GRAFT_DT, df$REC_TX_DT, units = "days"))/365.25

summary(df$tfl_graft_yrs - df$tfl_lafu_yrs) # 'tfl last FU' can be after graft has failed, plot below:
plot(df$tfl_graft_yrs, df$tfl_lafu_yrs)
# Importantly, graft never fails after last FU time, as expected

# Create binary variable 'tfl_graft_status' (1 graft failed, 0 failed) corresponding with 'tfl_graft_yrs'
df$tfl_graft_status <- as_factor(if_else(!is.na(df$tfl_graft_yrs), 1, 0))



# Calculate years from transplant to relisting
df$relist_yrs <- (as.numeric(difftime(df$PERS_RELIST, df$REC_TX_DT, units = "days"))/365.25)
length(which(df$relist_yrs < df$tfl_graft_yrs)) # Many patients were relisted before graft loss, which is possible
summary(df$relist_yrs - df$tfl_graft_yrs)
# For the purpose of later multi-state model, we need the relist time to start right after (1 day after) graft loss
df$relist_yrs <- if_else(df$relist_yrs <= df$tfl_graft_yrs, df$tfl_graft_yrs+0.0027, df$relist_yrs, missing = df$relist_yrs)
length(which(!is.na(df$relist_yrs)))
length(which(df$relist_yrs >= df$death_comp_yrs)) 
# Set 'relist_yrs' to NA if relisting happened on or after time of death
is.na(df$relist_yrs) <- which(df$relist_yrs >= df$death_comp_yrs)


# Create binary event variable 'relist' (1=yes, 0=no) corresponding to the 'relist_yrs' time variable
df$relist <- as.factor(if_else(!is.na(df$relist_yrs), 1, 0, missing = 0))

length(which((subset(df, relist == 0)$retx == 1))) # Only very few patients were retransplanted without having been relisted
table(subset(df, relist == 0 & retx == 1)$DON_TY) # Partly living donor, as expected, but some deceased (should be impossible)

is.na(df$relist) <- which(df$retx_yrs < df$relist_yrs) # If retx prior to relisting, set relist to NA
is.na(df$relist_yrs) <- which(df$retx_yrs < df$relist_yrs) # and time of relist to NA



# Discrepancies time of death and graft failure time
plot(df$tfl_graft_yrs, df$death_comp_yrs) 
df$temp2 <- (df$death_comp_yrs - df$tfl_graft_yrs)
df$temp2 <- (df$death_comp_yrs - df$tfl_graft_yrs)
length(subset(df$temp2, df$temp2 < (-0.25))) # Handful of cases of 'graft survival' more than 3 months after death
length(subset(df$temp2, df$temp2 < (0))) # Bit more cases graft survival any time after death
summary(subset(df, df$temp2 < (0))$TFL_LASTATUS) #but TFL last status is correct ('D') in all of these
# So TFL_LASTATUS is correct, but tfl_graft_yrs (derived from TFL_GRAFT_DT) is not (in these rare cases)

# Replace these tfl_graft_yrs with time of death
df$tfl_graft_yrs <- if_else(((df$death_comp_yrs - df$tfl_graft_yrs) < 0), df$death_comp_yrs, df$tfl_graft_yrs, 
                            missing = df$tfl_graft_yrs)


length(which(df$death_comp_yrs < df$tfl_graft_yrs)) # No cases of death prior to graft loss, as expected
length(which(df$death_comp_yrs == df$tfl_graft_yrs)) # >1000 cases of death same day as graft loss, technically possible but
# unclear what exactly happened in these cases. 
# Perhaps cases where graft loss was known to have occurred prior to death but date was unclear, and just set to same as death date
# Will resolve the tied cases of graft loss and death with least possible manipulation, by setting time of death to 1 day later.
# (because multi-state model would not be able to handle these ties).
# So these people are taken to have had graft failure prior to death. 
# Will change the '_yrs' variables only, as we'll be working with those (and not with the dates). 
length(which((df$death_comp_yrs - df$tfl_graft_yrs)==0)) # Confirms >1000 cases of death same day as graft loss
df$death_comp_yrs <- if_else((df$death_comp_yrs - df$tfl_graft_yrs)==0, (df$death_comp_yrs + 0.0027), df$death_comp_yrs, 
                             missing = df$death_comp_yrs)

# Do same for 'tfl_lafu_yrs' (although this will not really be used)
length(which(subset(df, TFL_LASTATUS == "D")$tfl_lafu_yrs == subset(df, TFL_LASTATUS == "D")$tfl_graft_yrs))
df$tfl_lafu_yrs <- if_else((df$tfl_lafu_yrs - df$tfl_graft_yrs)==0, (df$tfl_lafu_yrs + 0.0027), df$tfl_lafu_yrs, 
                           missing = df$tfl_lafu_yrs)

# Create variable 'ki_status' for which dead=2, graft failure=1, otherwise=0 (loss to FU or alive with functioning graft) - 
# whichever happened first - the corresponding event time variable being the new variable 'ki_yrs'. The kidney lasted until this point. 
df$ki_status <- as.factor(if_else(df$tfl_graft_status == 1, 1, ifelse(df$death_comp_status == 1, 2, 0), missing = 0))
df$temp <- if_else(df$tfl_graft_status == 1, df$tfl_graft_yrs, 
                      if_else(df$death_comp_status == 1, df$death_comp_yrs, df$tfl_lafu_comp_yrs))
df$ki_yrs <- ifelse(is.na(df$temp), df$tfl_lafu_comp_yrs, df$temp)


# Create ki_status_ACGF (1 = graft failure or death with functioning graft) to be used for KM 'all cause graft failure'
df$ki_status_ACGF <- as.factor(if_else(df$ki_status == 0, 0, 1))

length(which((df$tfl_graft_yrs-df$retx_yrs)>0)) # Rarely, time of retransplant was prior to graft loss, which should be impossible
# Assuming that retransplant is not an error, we will set 'graft loss time' to 'time of retransplant' in these cases
df$tfl_graft_yrs <- if_else((df$tfl_graft_yrs-df$retx_yrs)>0, df$retx_yrs, df$tfl_graft_yrs, missing = df$tfl_graft_yrs)



length(which((df$tfl_graft_yrs - df$retx_yrs)==0))
# Some cases where the time of tfl_graft_yrs is exactly the same as retx_yrs
# Maybe exact time of graft loss is not known and it was set as the same day, and maybe some patients where retransplanted pre-emptively
# Without having been relisted. The retx is unlikely to be an error, so set time of graft failure as 1 day prior to retransplant
# Because later multi-state model cannot deal with ties (events on same day)
# Change tfl_graft_yrs to 1 day earlier in these cases to resolve the ties
df$tfl_graft_yrs <- if_else((df$tfl_graft_status == 1 & (df$retx_yrs - df$tfl_graft_yrs) == 0), (df$tfl_graft_yrs - 0.0027), 
                     df$tfl_graft_yrs, missing = df$tfl_graft_yrs)


# Handful of patients had retransplant same day of death (indicating they died on day of surgery)
with(df, length(which((death_comp_status == 1 & retx == 1 & (retx_yrs - death_comp_yrs) == 0)))) # Number of pts
# Again, to resolve ties, set death time to 1 day after retransplant here
df$death_comp_yrs <- if_else(df$death_comp_status == 1 & df$retx == 1 & (df$retx_yrs - df$death_comp_yrs) == 0, 
                   df$death_comp_yrs + 0.0027, df$death_comp_yrs, missing = df$death_comp_yrs)

# For later multi-state: Use tfl_lafu_comp_yrs as last follow-up time variable for tmerge
length(df[df$tfl_lafu_comp_yrs==0,])
# Small number of cases of tfl_lafu_comp_yrs (last FU) on day 0, MSM will not accept this
# Add 0.001 to tfl_lafu_comp_yrs of '0' and call this 'last_fu'
df$last_fu <- ifelse(df$tfl_lafu_comp_yrs == 0, df$tfl_lafu_comp_yrs + 0.001, df$tfl_lafu_comp_yrs)

# Final data checks:
# The following should all be impossible (so vector should be length 0):
length(which(df$tfl_graft_yrs > df$tfl_lafu_yrs)) # Graft loss after last follow-up
length(which(df$retx_yrs > df$death_comp_yrs)) # retransplanted after death
length(which(df$retx_yrs > df$tfl_lafu_yrs)) # retransplanted after last FU (see above, this is why tfl_lafu_comp_yrs was created)
length(which(df$retx_yrs > df$tfl_lafu_comp_yrs)) # Same using tfl_lafu_comp_yrs
length(which(df$death_comp_yrs < df$tfl_graft_yrs)) # Death before graft loss
length(which(df$tfl_lafu_comp_yrs > df$death_comp_yrs)) # Last FU after death
length(which(df$relist_yrs >= df$death_comp_yrs)) # Relisted on day of death or after death
length(which(df$relist_yrs < df$tfl_graft_yrs)) # Relisted prior to graft loss
length(which(df$retx_yrs < df$relist_yrs)) # Retransplanted prior to relisting

# And the following could be possible but had to be addressed to avoid ties:
length(which((df$tfl_graft_yrs - df$retx_yrs)==0)) # graft loss same day as retransplant
length(which(df$death_comp_yrs == df$tfl_graft_yrs)) # death same day as graft loss
length(which(df$retx_yrs <= df$tfl_graft_yrs)) # retransplanted before or day of graft failure
length(which(df$relist_yrs == df$tfl_graft_yrs)) # Relisted on the day of graft loss

summary(df)



#################################################################################
# 3. CRUDE SURVIVAL ANALYSES
## 3.1 KAPLAN-MEIER
# Estimate survival rates with Kaplan-Meier
summary(survfit(Surv(ki_yrs, ki_status==1) ~ 1, data = df), times = c(1,3,5,10,15,20)) # Death-censored graft survival
summary(survfit(Surv(ki_yrs, ki_status==2) ~ 1, data = df), times = c(1,3,5,10,15,20)) # death w functioning graft
summary(survfit(Surv(ki_yrs, ki_status_ACGF==1) ~ 1, data = df), times = c(1,3,5,10,15,20)) # all-cause graft failure

# Same by age group
summary(survfit(Surv(ki_yrs, ki_status==1) ~ age_group_3, data = df), times = c(1,3,5,10,15,20)) # DC graft survival
summary(survfit(Surv(ki_yrs, ki_status==2) ~ age_group_3, data = df), times = c(1,3,5,10,15,20)) # death w functioning graft
summary(survfit(Surv(ki_yrs, ki_status_ACGF==1) ~ age_group_3, data = df), times = c(1,3,5,10,15,20)) # all-cause graft failure


# Kaplan-Meier (KM) for DC graft survival, by age group
km.g <- survfit(Surv(ki_yrs, ki_status==1) ~ age_group_3, data = df)

# Same, subset living donor
km.g.liv <- survfit(Surv(ki_yrs, ki_status==1) ~ age_group_3, data = subset(df, DON_TY == "L"))

# Same, subset deceased donor
km.g.dec <- survfit(Surv(ki_yrs, ki_status==1) ~ age_group_3, data = subset(df, DON_TY == "C"))
summary(km.g.dec, times = c(1,2))
summary(km.g.liv, times = c(1,2))


# Death with functioning graft (DWFG)
km.d <- survfit(Surv(ki_yrs, ki_status==2) ~ age_group_3, data = df)


# KM for all cause graft failure, by age group
km.acgf <- survfit(Surv(ki_yrs, ki_status_ACGF==1) ~ age_group_3, data = df)
summary(km.acgf, times = c(0,5,10,15,20))
# Plot
oldpar <- par(mar = c(5.1, 4.1, 4.1, 10))
plot(km.acgf, 
     col = c("#d7191c","#fdae61","#abd9e9","#2c7bb6"),
     xlim = c(0, 15),
     xlab = "Years", 
     ylab = "1 - Graft failure or death with functioning graft",
     main = "All-cause graft failure")
minor.tick(nx=5, ny=2, tick.ratio = 0.5)
legend("bottomleft", c("18 - 54","55 - 64","65 - 74","\u2265 75"), 
       col = c("#d7191c","#fdae61","#abd9e9","#2c7bb6"), 
       lty=c(1,1,1,1), lwd = 2, bty = "n", inset = c(1,0), xpd = NA)
par(oldpar)



# Overall patient survival (OS) Kaplan-Meier
# We must use tfl_lafu_comp_yrs rather than death_comp_yrs here, because the time of last follow-up in everyone 
# who didn't die is NA in death_comp_yrs (ie you don't know how long they are known to have survived before
# loss to FU or retransplant). This is NOT an issue for the multi-state models (where you start by defining 
# the time of last FU, so you can use death_comp_status to indicate when a patient died)
summary(survfit(Surv(tfl_lafu_comp_yrs, tfl_lafu_comp_status == 'D') ~ 1, data = df), times = c(1,3,5,10,15))
summary(survfit(Surv(ki_yrs, ki_status==2) ~ 1, data = df), times = c(1,3,5,10,15)) # Compare with DC graft survival

# OS by age group
summary(survfit(Surv(tfl_lafu_comp_yrs, tfl_lafu_comp_status == 'D') ~ age_group_3, data = df), times = c(1,3,5,10,15))
summary(survfit(Surv(ki_yrs, ki_status==2) ~ age_group_3, data = df), times = c(1,3,5,10,15)) # Compare with DC graft survival


# Plot for OS
# Whole population
km.os.all <- survfit(Surv(tfl_lafu_comp_yrs, tfl_lafu_comp_status == 'D') ~ 1, data = df)

# By age group
km.os <- survfit(Surv(tfl_lafu_comp_yrs, tfl_lafu_comp_status == 'D') ~ age_group_3, data = df)

# Extract list of the survival probability and time to use in ggplot
# DCGF
f.kmg.plotdat <- data.frame(survival = numeric(), lower = numeric(), upper = numeric(),  time = numeric(), 
                            age_group_3 = numeric())
datalist <- list()
for (i in 1:4) {
  tmp <- data.frame(survival = km.g[i]$surv, lower = km.g[i]$lower, upper = km.g[i]$upper,
                    time = km.g[i]$time, age_group_3 = rep(i, length(km.g[i]$time)))
  datalist[[i]] <- tmp
}
f.kmg.plotdat <- bind_rows(datalist)

# ACGF
f.kmacgf.plotdat <- data.frame(survival = numeric(), lower = numeric(), upper = numeric(),  
                               time = numeric(), age_group_3 = numeric())
datalist <- list()
for (i in 1:4) {
  tmp <- data.frame(survival = km.acgf[i]$surv, lower = km.acgf[i]$lower, upper = km.acgf[i]$upper,
                    time = km.acgf[i]$time, age_group_3 = rep(i, length(km.acgf[i]$time)))
  datalist[[i]] <- tmp
}
f.kmacgf.plotdat <- bind_rows(datalist)


# OS
f.kmos.plotdat <- data.frame(survival = numeric(), lower = numeric(), upper = numeric(),  
                             time = numeric(), age_group_3 = numeric())
datalist <- list()
for (i in 1:4) {
  tmp <- data.frame(survival = km.os[i]$surv, lower = km.os[i]$lower, upper = km.os[i]$upper,
                    time = km.os[i]$time, age_group_3 = rep(i, length(km.os[i]$time)))
  datalist[[i]] <- tmp
}
f.kmos.plotdat <- bind_rows(datalist)

# DWFG
f.kmd.plotdat <- data.frame(survival = numeric(), lower = numeric(), upper = numeric(),  
                            time = numeric(), age_group_3 = numeric())
datalist <- list()
for (i in 1:4) {
  tmp <- data.frame(survival = km.d[i]$surv, lower = km.d[i]$lower, upper = km.d[i]$upper,
                    time = km.d[i]$time, age_group_3 = rep(i, length(km.d[i]$time)))
  datalist[[i]] <- tmp
}
f.kmd.plotdat <- bind_rows(datalist)



# Plot KM curves manually (I prefer this over ggsurvplot)
p.kmg <- ggplot() + 
  geom_step(aes(x=time, y=survival), data = subset(f.kmg.plotdat, age_group_3 == 4), lwd = 0.7, col = "#2c7bb6", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmg.plotdat, age_group_3 == 3), lwd = 0.7, col = "#abd9e9", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmg.plotdat, age_group_3 == 2), lwd = 0.7, col = "#fdae61", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmg.plotdat, age_group_3 == 1), lwd = 0.7, col = "#d7191c", lty = 1) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0.5,1.0), breaks = seq(0.5,1,0.05),
                     labels = c(0.5,"",0.6,"",0.7,"",0.8,"",0.9,"","1.0")) +
  labs(x = "", y = "Survival", title = "Graft survival, death-censored") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

p.kmd <- ggplot() + 
  geom_step(aes(x=time, y=survival), data = subset(f.kmd.plotdat, age_group_3 == 1), lwd = 0.7, col = "#d7191c", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmd.plotdat, age_group_3 == 2), lwd = 0.7, col = "#fdae61", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmd.plotdat, age_group_3 == 3), lwd = 0.7, col = "#abd9e9", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmd.plotdat, age_group_3 == 4), lwd = 0.7, col = "#2c7bb6", lty = 1) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "Years after transplant", y = "Survival", title = "Death with functioning graft") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

p.kmos <- ggplot() + 
  geom_step(aes(x=time, y=survival), data = subset(f.kmos.plotdat, age_group_3 == 1), lwd = 0.7, col = "#d7191c", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmos.plotdat, age_group_3 == 2), lwd = 0.7, col = "#fdae61", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmos.plotdat, age_group_3 == 3), lwd = 0.7, col = "#abd9e9", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmos.plotdat, age_group_3 == 4), lwd = 0.7, col = "#2c7bb6", lty = 1) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "Years after transplant", y = "", title = "Overall patient survival") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

p.kmacgf <- ggplot() + 
  geom_step(aes(x=time, y=survival), data = subset(f.kmacgf.plotdat, age_group_3 == 1), lwd = 0.7, col = "#d7191c", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmacgf.plotdat, age_group_3 == 2), lwd = 0.7, col = "#fdae61", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmacgf.plotdat, age_group_3 == 3), lwd = 0.7, col = "#abd9e9", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmacgf.plotdat, age_group_3 == 4), lwd = 0.7, col = "#2c7bb6", lty = 1) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "", y = "", title = "All-cause (uncensored) graft survival") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

pw.km <- (p.kmg+p.kmacgf)/(p.kmd+p.kmos)
pw.km + plot_annotation(tag_levels = "A")


# Same but 5 yr plots with confidence intervals
p.kmg <- ggplot() + 
  geom_step(aes(x=time, y=survival), data = subset(f.kmg.plotdat, age_group_3 == 4), lwd = 0.7, col = "#2c7bb6", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmg.plotdat, age_group_3 == 3), lwd = 0.7, col = "#abd9e9", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmg.plotdat, age_group_3 == 2), lwd = 0.7, col = "#fdae61", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmg.plotdat, age_group_3 == 1), lwd = 0.7, col = "#d7191c", lty = 1) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmg.plotdat, age_group_3 == 4), 
              alpha = 0.1, col = "#2c7bb6", fill = "#2c7bb6", lty = 1, lwd = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmg.plotdat, age_group_3 == 3), 
              alpha = 0.1, col = "#abd9e9", fill = "#abd9e9", lty = 1, lwd = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmg.plotdat, age_group_3 == 2), 
              alpha = 0.1, col = "#fdae61", fill = "#fdae61", lty = 1, lwd = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmg.plotdat, age_group_3 == 1), 
              alpha = 0.1, col = "#d7191c", fill = "#d7191c", lty = 1, lwd = 0.2) +
  scale_x_continuous(limits = c(0,5), breaks = seq(0,5,1)) +
  scale_y_continuous(limits = c(0.8,1.0), breaks = seq(0.8,1,0.02),
                     labels = c(0.8,rep("",4), 0.9, rep("",4), "1.0")) +
  labs(x = "", y = "Survival", title = "Graft survival, death-censored") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

p.kmd <- ggplot() + 
  geom_step(aes(x=time, y=survival), data = subset(f.kmd.plotdat, age_group_3 == 1), lwd = 0.7, col = "#d7191c", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmd.plotdat, age_group_3 == 2), lwd = 0.7, col = "#fdae61", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmd.plotdat, age_group_3 == 3), lwd = 0.7, col = "#abd9e9", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmd.plotdat, age_group_3 == 4), lwd = 0.7, col = "#2c7bb6", lty = 1) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmd.plotdat, age_group_3 == 4), 
              alpha = 0.1, col = "#2c7bb6", fill = "#2c7bb6", lty = 1, lwd = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmd.plotdat, age_group_3 == 3), 
              alpha = 0.1, col = "#abd9e9", fill = "#abd9e9", lty = 1, lwd = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmd.plotdat, age_group_3 == 2), 
              alpha = 0.1, col = "#fdae61", fill = "#fdae61", lty = 1, lwd = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmd.plotdat, age_group_3 == 1), 
              alpha = 0.1, col = "#d7191c", fill = "#d7191c", lty = 1, lwd = 0.2) +
  scale_x_continuous(limits = c(0,5), breaks = seq(0,5,1)) +
  scale_y_continuous(limits = c(0.5,1.0), breaks = seq(0.5,1,0.05),
                     labels = c(0.5,"",0.6,"",0.7,"",0.8,"",0.9,"","1.0")) +
  labs(x = "Years after transplant", y = "Survival", title = "Death with functioning graft") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

p.kmos <- ggplot() + 
  geom_step(aes(x=time, y=survival), data = subset(f.kmos.plotdat, age_group_3 == 1), lwd = 0.7, col = "#d7191c", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmos.plotdat, age_group_3 == 2), lwd = 0.7, col = "#fdae61", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmos.plotdat, age_group_3 == 3), lwd = 0.7, col = "#abd9e9", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmos.plotdat, age_group_3 == 4), lwd = 0.7, col = "#2c7bb6", lty = 1) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmos.plotdat, age_group_3 == 4), 
              alpha = 0.1, col = "#2c7bb6", fill = "#2c7bb6", lty = 1, lwd = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmos.plotdat, age_group_3 == 3), 
              alpha = 0.1, col = "#abd9e9", fill = "#abd9e9", lty = 1, lwd = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmos.plotdat, age_group_3 == 2), 
              alpha = 0.1, col = "#fdae61", fill = "#fdae61", lty = 1, lwd = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmos.plotdat, age_group_3 == 1), 
              alpha = 0.1, col = "#d7191c", fill = "#d7191c", lty = 1, lwd = 0.2) +
  scale_x_continuous(limits = c(0,5), breaks = seq(0,5,1)) +
  scale_y_continuous(limits = c(0.5,1.0), breaks = seq(0.5,1,0.05),
                     labels = c(0.5,"",0.6,"",0.7,"",0.8,"",0.9,"","1.0")) +
  labs(x = "Years after transplant", y = "", title = "Overall patient survival") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

p.kmacgf <- ggplot() + 
  geom_step(aes(x=time, y=survival), data = subset(f.kmacgf.plotdat, age_group_3 == 1), lwd = 0.7, col = "#d7191c", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmacgf.plotdat, age_group_3 == 2), lwd = 0.7, col = "#fdae61", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmacgf.plotdat, age_group_3 == 3), lwd = 0.7, col = "#abd9e9", lty = 1) +
  geom_step(aes(x=time, y=survival), data = subset(f.kmacgf.plotdat, age_group_3 == 4), lwd = 0.7, col = "#2c7bb6", lty = 1) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmacgf.plotdat, age_group_3 == 4), 
              alpha = 0.1, col = "#2c7bb6", fill = "#2c7bb6", lty = 1, lwd = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmacgf.plotdat, age_group_3 == 3), 
              alpha = 0.1, col = "#abd9e9", fill = "#abd9e9", lty = 1, lwd = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmacgf.plotdat, age_group_3 == 2), 
              alpha = 0.1, col = "#fdae61", fill = "#fdae61", lty = 1, lwd = 0.2) +
  geom_ribbon(aes(ymin = lower, ymax=upper, x=time), data = subset(f.kmacgf.plotdat, age_group_3 == 1), 
              alpha = 0.1, col = "#d7191c", fill = "#d7191c", lty = 1, lwd = 0.2) +
  scale_x_continuous(limits = c(0,5), breaks = seq(0,5,1)) +
  scale_y_continuous(limits = c(0.5,1.0), breaks = seq(0.5,1,0.05),
                     labels = c(0.5,"",0.6,"",0.7,"",0.8,"",0.9,"","1.0")) +
  labs(x = "", y = "", title = "All-cause (uncensored) graft survival") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

pw.km <- (p.kmg+p.kmacgf)/(p.kmd+p.kmos)
pw.km + plot_annotation(tag_levels = "A")

# ggsave("KM CI.tiff", device = "tiff", dpi = 300, width = 20, height = 16, units = "cm", 
       path = "C:/Users/Thomas Vanhove/Box Sync/R/SRTR/Plots")







##################################################################################################
## 3.2 COMPETING RISK 
# Using Survfit, so Aalen-Johansen
cr.fit <- survfit(Surv(ki_yrs, ki_status) ~ 1, data=df)
print(cr.fit)
cr.fit$transitions # Check that no pt leaves state 1 or 2 (they are both treated as terminal states)

plot(cr.fit, col = c(1,2), lwd = 2, xlab="Years", ylab="Probability in State", ylim = c(0,1), xlim = c(0,15))
legend("topleft", c("Graft failure", "Death with functioning graft"), col = c(1,2), lwd = 2, bty = "n")


# By age group
cr.fit2 <- survfit(Surv(ki_yrs, ki_status) ~ age_group_3, data=df)
print(cr.fit2)

oldpar <- par(mar = c(5.1, 4.1, 4.1, 12))
plot(cr.fit2, col = c("#d7191c","#fdae61","#abd9e9","#2c7bb6"), lty=c(2,2,2,2,1,1,1,1), lwd = 2, 
     xlab="Years", ylab="Probability in State", xlim = c(0,15), ylim = c(0,1))
minor.tick(nx=5, ny=2, tick.ratio = 0.5)
legend("topleft", c("Death : 18-55","Death : 56-65",
                    "Death : 66-75","Death : >75",
                    "Graft failure : 18-55", "Graft failure : 56-65","Graft failure : 66-75","Graft failure : >75"), 
       col = c("#d7191c","#fdae61","#abd9e9","#2c7bb6"), 
       lty=c(1,1,1,1,2,2,2,2), lwd = 2, bty = "n", inset = c(1,0), xpd = NA)
par(oldpar)



# Plot Cumulative incidence function (CIF) from competing risk analysis
# Graft failure, CR estimate
plotCIF(cr.fit2, event = 1, 
        col = c("#d7191c","#fdae61","#abd9e9","#2c7bb6"), 
        lty=c(1,1,1,1), 
        lwd = 2, 
        xlab="Years", 
        ylab="Cumulative incidence",
        main = "Graft failed, alive",
        xlim = c(0,15))
legend("topleft", c("18 - 54","55 - 64","65 - 74","\u2265 75"), 
       col = c("#d7191c","#fdae61","#abd9e9","#2c7bb6"), 
       lty=c(1,1,1,1), lwd = 2, bty = "n")


# Plot just CR estimate of DWFG
plotCIF(cr.fit2, event = 2, 
        col = c("#d7191c","#fdae61","#abd9e9","#2c7bb6"), 
        lty=c(1,1,1,1), 
        lwd = 2, 
        xlab="Years", 
        ylab="Cumulative incidence",
        main = "Death with functioning graft",
        xlim = c(0,15))
legend("topleft", c("18 - 54","55 - 64","65 - 74","\u2265 75"), 
       col = c("#d7191c","#fdae61","#abd9e9","#2c7bb6"), 
       lty=c(1,1,1,1), lwd = 2, bty = "n")


# Graft failure AND DWFG in stacked CIF plot
stackedCIF(cr.fit, 
        lwd = 2, 
        fill = c("#f0f0f0","#bdbdbd"),
        xlab="Years", 
        ylab="Stacked cumulative incidence",
        main = "",
        xlim = c(0,15))
legend("topleft", c("Death with functioning graft","Graft failed, alive"), 
       col = c("#bdbdbd","#f0f0f0"),
       lwd = 2, 
       bty = "n")


# Same, stratified by age groups
oldpar <- par(mai=c(0.5, 0.5, 0.5, 0.5), oma=c(0.5,2,0.5,0.5), mgp=c(2,0.5,0) , mfrow = c(2,2))
stackedCIF(cr.fit2, 
           group = 1,
           lwd = 2, 
           fill = c("grey90","grey60"),
           xlab="", 
           ylab="Cumulative incidence",
           main = "Age 18 - 54",
           xlim = c(0,15), ylim = c(0,1), las = 1)
minor.tick(nx=5, ny=2, tick.ratio = 0.5)
legend("topleft", c("Alive with functioning graft","Death with functioning graft","Graft failed"),
       col = c("black","grey60","grey90"), 
       bty = "n",
       pch = c(0,15,15), pt.cex = 2)

stackedCIF(cr.fit2, 
           group = 2,
           lwd = 2, 
           fill = c("grey90","grey60"),
           xlab="", 
           ylab="",
           main = "Age 55 - 64",
           xlim = c(0,15), ylim = c(0,1), las = 1)
minor.tick(nx=5, ny=2, tick.ratio = 0.5)

stackedCIF(cr.fit2, 
           group = 3,
           lwd = 2, 
           fill = c("grey90","grey60"),
           xlab="Years after transplant", 
           ylab="Cumulative incidence",
           main = "Age 65 - 74",
           xlim = c(0,15), ylim = c(0,1), las = 1)
minor.tick(nx=5, ny=2, tick.ratio = 0.5)

stackedCIF(cr.fit2, 
           group = 4,
           lwd = 2, 
           fill = c("grey90","grey60"),
           xlab="Years after transplant", 
           ylab="",
           main = "Age \u2265 75",
           xlim = c(0,15), ylim = c(0,1), las = 1)
minor.tick(nx=5, ny=2, tick.ratio = 0.5)
par(oldpar)


# Same plot with color
oldpar <- par(mai=c(0.5, 0.5, 0.5, 0.5), oma=c(0.5,2,0.5,0.5), mgp=c(2,0.5,0) , mfrow = c(2,2))
stackedCIF(cr.fit2, 
           group = 1,
           lwd = 2, 
           fill = c("#56B4E9","#E69F00"),
           xlab="", 
           ylab="Cumulative incidence",
           main = "Age 18 - 54",
           xlim = c(0,15), ylim = c(0,1), las = 1)
minor.tick(nx=5, ny=2, tick.ratio = 0.5)
legend("topleft", c("Alive with functioning graft","Death with functioning graft","Graft failed"), 
       col = c("black","#E69F00","#56B4E9"), 
       bty = "n",
       pch = c(0,15,15), pt.cex = 2)

stackedCIF(cr.fit2, 
           group = 2,
           lwd = 2, 
           fill = c("#56B4E9","#E69F00"),
           xlab="", 
           ylab="",
           main = "Age 55 - 64",
           xlim = c(0,15), ylim = c(0,1), las = 1)
minor.tick(nx=5, ny=2, tick.ratio = 0.5)

stackedCIF(cr.fit2, 
           group = 3,
           lwd = 2, 
           fill = c("#56B4E9","#E69F00"),
           xlab="Years after transplant", 
           ylab="Cumulative incidence",
           main = "Age 65 - 74",
           xlim = c(0,15), ylim = c(0,1), las = 1)
minor.tick(nx=5, ny=2, tick.ratio = 0.5)

stackedCIF(cr.fit2, 
           group = 4,
           lwd = 2, 
           fill = c("#56B4E9","#E69F00"),
           xlab="Years after transplant", 
           ylab="",
           main = "Age \u2265 75",
           xlim = c(0,15), ylim = c(0,1), las = 1)
minor.tick(nx=5, ny=2, tick.ratio = 0.5)
par(oldpar)












##################################################################################################
# 3.3 MULTI-STATE MODEL (MSM)
# Create appropriate 'start-stop' dataframe using 'tmerge'
df.ss <- tmerge(data1 = df[,c("PERS_ID", "DONOR_ID", "REC_TX_DT","REC_AGE_AT_TX","DON_AGE","REC_MM_EQUIV_TX","REC_COLD_ISCH_TM",
                              "DON_TY","CAN_GENDER", "age_group","age_group_2","age_group_3","age_group_4",
                              "CAN_DIAB", "CAN_RACE_SRTR","CAN_ETHNICITY_SRTR","REC_BMI",
                              "CAN_CEREB_VASC","CAN_ANGINA")], 
                 data2 = df, id = PERS_ID, tstop = last_fu)
df.ss <- tmerge(data1 = df.ss, data2 = df, id = PERS_ID, graftloss = event(tfl_graft_yrs, tfl_graft_status),
                 death_comp = event(death_comp_yrs, death_comp_status))
df.ss <- tmerge(df.ss, df.ss, id = PERS_ID, enum=cumtdc(tstart))
temp <- if_else(df.ss$graftloss == 1, 1, if_else(df.ss$death_comp == 1, 2, 0))
df.ss$event <- factor(temp, 0:2, labels=c("censor", "gl", "death")) # 'event' = 3 possible outcomes
df.ss$tstate <- with(df.ss, tstop - tstart)

dim(df.ss)
summary(df.ss)
attr(df.ss, "tcount")
with(df.ss, table(graftloss, death_comp)) # No ties (no graft loss and death on same row)

msm.1 <- survfit(Surv(tstart, tstop, event) ~ 1, data=df.ss, id=PERS_ID)
plot(msm.1, ylim = c(0,1), xlim = c(0,20))
print(msm.1)


# Create event 'death after graft loss' 
d2 <- with(df.ss, ifelse(enum==2 & event=='death', 4, as.numeric(event)))
e2 <- factor(d2, labels=c("censor", "gl", "death w functioning graft",
                            "death after gl"))


msm.2 <- survfit(Surv(tstart, tstop, e2) ~ 1, data=df.ss, id=PERS_ID)

plot(msm.2, lty=c(3,1,5), lwd = 2, 
     xlab="Years", ylab="Probability in State", xlim = c(0,15), ylim = c(0,0.5))
minor.tick(nx=5, ny=2, tick.ratio = 0.5)
legend("topleft", c("Failed graft, alive", "Death with functioning graft", "Death after graft failure"), 
       lty=c(3,1,5), lwd = 2, bty = "n", xpd = NA)


# Use noplot=NULL option to also show the baseline state of 'alive with functioning graft'
# Although this is redundant because all probabilities must sum to 1
plot(msm.2, lty=c(6,3,1,5), lwd = 2, 
     xlab="Years", ylab="Probability in State", xlim = c(0,15), noplot = NULL)
minor.tick(nx=5, ny=1, tick.ratio = 0.5)
legend("topright", c("Alive with functioning graft","Failed graft, alive", "Death with functioning graft", 
                     "Death after graft failure"), 
       lty=c(6,3,1,5), lwd = 2, bty = "n", xpd = NA)


# Create individual plots msm.2 for the age groups and arrange together
msm.2a <- survfit(Surv(tstart, tstop, e2) ~ 1, data=df.ss, subset = df.ss$age_group_3 == "(0,55]", id=PERS_ID)
msm.2b <- survfit(Surv(tstart, tstop, e2) ~ 1, data=df.ss, subset = df.ss$age_group_3 == "(55,65]", id=PERS_ID)
msm.2c <- survfit(Surv(tstart, tstop, e2) ~ 1, data=df.ss, subset = df.ss$age_group_3 == "(65,75]", id=PERS_ID)
msm.2d <- survfit(Surv(tstart, tstop, e2) ~ 1, data=df.ss, subset = df.ss$age_group_3 == "(75,100]", id=PERS_ID)

msm.2a.plotdat <- as_tibble(msm.2a$pstate)
msm.2a.plotdat <- cbind(msm.2a.plotdat, msm.2a$time)
colnames(msm.2a.plotdat) <- c("censor", "gl", "dwfg", "dagf","time")

msm.2b.plotdat <- as_tibble(msm.2b$pstate)
msm.2b.plotdat <- cbind(msm.2b.plotdat, msm.2b$time)
colnames(msm.2b.plotdat) <- c("censor", "gl", "dwfg", "dagf","time")

msm.2c.plotdat <- as_tibble(msm.2c$pstate)
msm.2c.plotdat <- cbind(msm.2c.plotdat, msm.2c$time)
colnames(msm.2c.plotdat) <- c("censor", "gl", "dwfg", "dagf","time")

msm.2d.plotdat <- as_tibble(msm.2d$pstate)
msm.2d.plotdat <- cbind(msm.2d.plotdat, msm.2d$time)
colnames(msm.2d.plotdat) <- c("censor", "gl", "dwfg", "dagf","time")



p.msm.2.1 <- ggplot() + 
  geom_step(aes(x=time, y=censor), data = msm.2a.plotdat, col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=censor), data = msm.2b.plotdat, col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=censor), data = msm.2c.plotdat, col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=censor), data = msm.2d.plotdat, col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "", y = "Probability in state", title = "Alive with functioning graft") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 

p.msm.2.2 <- ggplot() + 
  geom_step(aes(x=time, y=dwfg), data = msm.2a.plotdat, col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=dwfg), data = msm.2b.plotdat, col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=dwfg), data = msm.2c.plotdat, col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=dwfg), data = msm.2d.plotdat, col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "", y = "", title = "Death with functioning graft") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 


p.msm.2.3 <- ggplot() + 
  geom_step(aes(x=time, y=gl), data = msm.2a.plotdat, col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=gl), data = msm.2b.plotdat, col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=gl), data = msm.2c.plotdat, col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=gl), data = msm.2d.plotdat, col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,0.5), breaks = seq(0,1,0.1)) +
  labs(x = "Years", y = "Probability in state", title = "Graft failed, alive") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 


p.msm.2.4 <- ggplot() + 
  geom_step(aes(x=time, y=dagf), data = msm.2a.plotdat, col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=dagf), data = msm.2b.plotdat, col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=dagf), data = msm.2c.plotdat, col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=dagf), data = msm.2d.plotdat, col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,0.5), breaks = seq(0,1,0.1)) +
  labs(x = "Years", y = "", title = "Death after graft failure") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))  

pw.msm.2 <- (p.msm.2.1+p.msm.2.2)/(p.msm.2.3+p.msm.2.4)
pw.msm.2 + plot_annotation(tag_levels = "A")





# Create second start-stop file with a 'retransplant' event
# Death after retx will be coded as death with functioning graft (dwfg) in event 'e3' 
# In event 'e4', all deaths (dwfg, dagf, death after retx) will be lumped together in 'os' (overall survival)
df.ss2 <- tmerge(data1 = df[,c("PERS_ID", "DONOR_ID", "REC_TX_DT","REC_AGE_AT_TX","DON_AGE","REC_MM_EQUIV_TX","REC_COLD_ISCH_TM",
                               "DON_TY","CAN_GENDER", "age_group","age_group_2","age_group_3","age_group_4",
                               "CAN_DIAB", "CAN_RACE_SRTR","CAN_ETHNICITY_SRTR","REC_BMI",
                               "CAN_CEREB_VASC","CAN_ANGINA","retx","death_comp_status","relist","ki_status")], 
                  data2 = df, id = PERS_ID, tstop = last_fu)
df.ss2 <- tmerge(data1 = df.ss2, data2 = df, id = PERS_ID, graftloss = event(tfl_graft_yrs, tfl_graft_status), 
                  retx = event(retx_yrs, retx), death_comp = event(death_comp_yrs, death_comp_status))
df.ss2 <- tmerge(df.ss2, df.ss2, id = PERS_ID, enum=cumtdc(tstart))

temp <- if_else(df.ss2$graftloss == 1, 1, if_else(df.ss2$retx == 1, 2, if_else(df.ss2$death_comp == 1, 3, 0)), missing = 0)
table(temp)
df.ss2$event <- factor(temp, 0:3, labels=c("censor", "gl", "retx", "death")) # 'event' = 4 possible outcomes

dim(df.ss2)
summary(df.ss2)
attr(df.ss2, "tcount")
with(df.ss2, table(graftloss, death_comp, retx)) # No ties (no graft loss and death on same row)

msm.retx <- survfit(Surv(tstart, tstop, event) ~ 1, data=df.ss2, id=PERS_ID)
plot(msm.retx, xlim = c(0,15))
print(msm.retx)

# Create event 'death after graft loss' 
d3 <- with(df.ss2, ifelse(enum==2 & event=='death', 5, as.numeric(event)))
df.ss2$e3 <- factor(d3, labels=c("censor", "gl", "retx", "death w functioning graft",
                          "death after gl"))
# Create e3alt, which has the events in a different order (useful later for stacked CIF plots)
df.ss2$e3alt <- factor(as.character(df.ss2$e3, levels=c("censor","death after gl","death w functioning graft","retx","gl")))
df.ss2$e4 <- recode_factor(df.ss2$e3, 'censor'="censor", 'gl'="gl", 'retx'="retx", 'death w functioning graft'="death",
                            'death after gl'="death") # event 'e4' with only 1 state for death


# Plot stacked CIF for endpoint e3
msm.retx.2 <- survfit(Surv(tstart, tstop, e3) ~ 1, data=df.ss2, id=PERS_ID)

stackedCIF(msm.retx.2, 
           fill = c("#f7f7f7", "#cccccc", "#969696", "#525252"),
           xlab="Years", ylab="Stacked cumulative incidence", xlim = c(0,15))
minor.tick(nx=5, ny=2, tick.ratio = 0.5)






# Create individual plots endpoint e3 for the age groups and arrange together
# Run separate models
msm.retx.2a <- survfit(Surv(tstart, tstop, e3) ~ 1, data=df.ss2, subset = df.ss2$age_group_3 == "(0,55]", id=PERS_ID)
msm.retx.2b <- survfit(Surv(tstart, tstop, e3) ~ 1, data=df.ss2, subset = df.ss2$age_group_3 == "(55,65]", id=PERS_ID)
msm.retx.2c <- survfit(Surv(tstart, tstop, e3) ~ 1, data=df.ss2, subset = df.ss2$age_group_3 == "(65,75]", id=PERS_ID)
msm.retx.2d <- survfit(Surv(tstart, tstop, e3) ~ 1, data=df.ss2, subset = df.ss2$age_group_3 == "(75,100]", id=PERS_ID)


# Create dataframe of data to plot
fit.list <- list(msm.retx.2a, msm.retx.2b, msm.retx.2c, msm.retx.2d)

p.msm.retx.plotdat <- data.frame(alive = numeric(), gfa = numeric(), retx = numeric(), dwfg = numeric(), 
                                 dagf = numeric(), time = numeric())
datalist <- list()
for (i in 1:4) {
  tmp <- data.frame(alive = fit.list[[i]]$pstate[,1], gfa = fit.list[[i]]$pstate[,2], retx = fit.list[[i]]$pstate[,3], 
                    dwfg = fit.list[[i]]$pstate[,4], dagf = fit.list[[i]]$pstate[,5], time = fit.list[[i]]$time,
                    age_group_3 = rep(i, length(fit.list[[i]]$time)))
  datalist[[i]] <- tmp
}
p.msm.retx.plotdat <- bind_rows(datalist)


p.msm.retx.1 <- ggplot() + 
  geom_step(aes(x=time, y=alive), data = subset(p.msm.retx.plotdat, age_group_3 == 1), col = "#d7191c", lwd = 0.5) +
  geom_step(aes(x=time, y=alive), data = subset(p.msm.retx.plotdat, age_group_3 == 2), col = "#fdae61", lwd = 0.5) +
  geom_step(aes(x=time, y=alive), data = subset(p.msm.retx.plotdat, age_group_3 == 3), col = "#abd9e9", lwd = 0.5) +
  geom_step(aes(x=time, y=alive), data = subset(p.msm.retx.plotdat, age_group_3 == 4), col = "#2c7bb6", lwd = 0.5) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1), labels = c(0,"",0.2,"",0.4,"",0.6,"",0.8,"","1.0")) +
  labs(x = "", y = "Probability in state", title = "Alive with functioning graft") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 

p.msm.retx.2 <- ggplot() + 
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.retx.plotdat, age_group_3 == 1), col = "#d7191c", lwd = 0.5) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.retx.plotdat, age_group_3 == 2), col = "#fdae61", lwd = 0.5) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.retx.plotdat, age_group_3 == 3), col = "#abd9e9", lwd = 0.5) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.retx.plotdat, age_group_3 == 4), col = "#2c7bb6", lwd = 0.5) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,0.2), breaks = seq(0,0.2,0.02), 
                     labels = c(0,rep("", 4), 0.1,rep("", 4),0.2)) +
  labs(x = "", y = "Probability in state", title = "Graft failed, alive") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 

p.msm.retx.3 <- ggplot() + 
  geom_step(aes(x=time, y=retx), data = subset(p.msm.retx.plotdat, age_group_3 == 1), col = "#d7191c", lwd = 0.5) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.retx.plotdat, age_group_3 == 2), col = "#fdae61", lwd = 0.5) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.retx.plotdat, age_group_3 == 3), col = "#abd9e9", lwd = 0.5) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.retx.plotdat, age_group_3 == 4), col = "#2c7bb6", lwd = 0.5) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,0.2), breaks = seq(0,0.2,0.02), 
                     labels = c(0,rep("", 4), 0.1,rep("", 4),0.2)) +
  labs(x = "Years after transplant", y = "Probability in state", title = "Retransplanted") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 

p.msm.retx.4 <- ggplot() + 
  geom_step(aes(x=time, y=dwfg), data = subset(p.msm.retx.plotdat, age_group_3 == 1), col = "#d7191c", lwd = 0.5) +
  geom_step(aes(x=time, y=dwfg), data = subset(p.msm.retx.plotdat, age_group_3 == 2), col = "#fdae61", lwd = 0.5) +
  geom_step(aes(x=time, y=dwfg), data = subset(p.msm.retx.plotdat, age_group_3 == 3), col = "#abd9e9", lwd = 0.5) +
  geom_step(aes(x=time, y=dwfg), data = subset(p.msm.retx.plotdat, age_group_3 == 4), col = "#2c7bb6", lwd = 0.5) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1), labels = c(0,"",0.2,"",0.4,"",0.6,"",0.8,"","1.0")) +
  labs(x = "", y = "", title = "Death with functioning graft") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")  

p.msm.retx.5 <- ggplot() + 
  geom_step(aes(x=time, y=dagf), data = subset(p.msm.retx.plotdat, age_group_3 == 1), col = "#d7191c", lwd = 0.5) +
  geom_step(aes(x=time, y=dagf), data = subset(p.msm.retx.plotdat, age_group_3 == 2), col = "#fdae61", lwd = 0.5) +
  geom_step(aes(x=time, y=dagf), data = subset(p.msm.retx.plotdat, age_group_3 == 3), col = "#abd9e9", lwd = 0.5) +
  geom_step(aes(x=time, y=dagf), data = subset(p.msm.retx.plotdat, age_group_3 == 4), col = "#2c7bb6", lwd = 0.5) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,0.2), breaks = seq(0,0.2,0.02), 
                     labels = c(0,rep("", 4), 0.1,rep("", 4),0.2)) +
  labs(x = "Years after transplant", y = "", title = "Death after graft failure") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")  

p.blank <- ggplot() + geom_blank() + theme_minimal() + labs(title = "Model diagram") +
  theme(plot.title = element_text(hjust = 0.5)) 


pw.msm.retx.e3 <- (p.msm.retx.1+p.msm.retx.4)/(p.msm.retx.2+p.msm.retx.5)/(p.msm.retx.3+p.blank)
pw.msm.retx.e3 + plot_annotation(tag_levels = "A")
# ggsave("MSM 5.tiff", device = "tiff", dpi = 300, width = 16, height = 20, units = "cm", 
       path = "C:/Users/Thomas Vanhove/Box Sync/R/SRTR/Plots")





# Same with event 'e4', in which there is only 1 state for death
# Create individual plots msm.retx for the age groups and arrange together
msm.retx.3a <- survfit(Surv(tstart, tstop, e4) ~ 1, data=df.ss2, subset = df.ss2$age_group_3 == "(0,55]", id=PERS_ID)
msm.retx.3b <- survfit(Surv(tstart, tstop, e4) ~ 1, data=df.ss2, subset = df.ss2$age_group_3 == "(55,65]", id=PERS_ID)
msm.retx.3c <- survfit(Surv(tstart, tstop, e4) ~ 1, data=df.ss2, subset = df.ss2$age_group_3 == "(65,75]", id=PERS_ID)
msm.retx.3d <- survfit(Surv(tstart, tstop, e4) ~ 1, data=df.ss2, subset = df.ss2$age_group_3 == "(75,100]", id=PERS_ID)

fit.list <- list(msm.retx.3a, msm.retx.3b, msm.retx.3c, msm.retx.3d)

p.msm.retx2.plotdat <- data.frame(alive = numeric(), gfa = numeric(), retx = numeric(), death = numeric(), 
                                 time = numeric())
datalist <- list()
for (i in 1:4) {
  tmp <- data.frame(alive = fit.list[[i]]$pstate[,1], gfa = fit.list[[i]]$pstate[,2], retx = fit.list[[i]]$pstate[,3], 
                    death = fit.list[[i]]$pstate[,4], time = fit.list[[i]]$time,
                    age_group_3 = rep(i, length(fit.list[[i]]$time)))
  datalist[[i]] <- tmp
}
p.msm.retx2.plotdat <- bind_rows(datalist)

p1.e4 <- ggplot() + 
  geom_step(aes(x=time, y=alive), data = subset(p.msm.retx2.plotdat, age_group_3 == 1), col = "#d7191c", lwd = 0.5) +
  geom_step(aes(x=time, y=alive), data = subset(p.msm.retx2.plotdat, age_group_3 == 2), col = "#fdae61", lwd = 0.5) +
  geom_step(aes(x=time, y=alive), data = subset(p.msm.retx2.plotdat, age_group_3 == 3), col = "#abd9e9", lwd = 0.5) +
  geom_step(aes(x=time, y=alive), data = subset(p.msm.retx2.plotdat, age_group_3 == 4), col = "#2c7bb6", lwd = 0.5) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "", y = "Probability in state", title = "Alive with functioning graft") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 

p2.e4 <- ggplot() + 
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.retx2.plotdat, age_group_3 == 1), col = "#d7191c", lwd = 0.5) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.retx2.plotdat, age_group_3 == 2), col = "#fdae61", lwd = 0.5) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.retx2.plotdat, age_group_3 == 3), col = "#abd9e9", lwd = 0.5) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.retx2.plotdat, age_group_3 == 4), col = "#2c7bb6", lwd = 0.5) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,0.2), breaks = seq(0,0.2,0.02),
                     labels = c(0, rep("",4), 0.1, rep("",4), 0.2)) +
  labs(x = "Years", y = "Probability in state", title = "Graft failed, alive") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 


p3.e4 <- ggplot() + 
  geom_step(aes(x=time, y=retx), data = subset(p.msm.retx2.plotdat, age_group_3 == 1), col = "#d7191c", lwd = 0.5) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.retx2.plotdat, age_group_3 == 2), col = "#fdae61", lwd = 0.5) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.retx2.plotdat, age_group_3 == 3), col = "#abd9e9", lwd = 0.5) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.retx2.plotdat, age_group_3 == 4), col = "#2c7bb6", lwd = 0.5) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,0.2), breaks = seq(0,0.2,0.02),
                     labels = c(0, rep("",4), 0.1, rep("",4), 0.2)) +
  labs(x = "Years", y = "", title = "Retransplanted") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 


p4.e4 <- ggplot() + 
  geom_step(aes(x=time, y=death), data = subset(p.msm.retx2.plotdat, age_group_3 == 1), col = "#d7191c", lwd = 0.5) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.retx2.plotdat, age_group_3 == 2), col = "#fdae61", lwd = 0.5) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.retx2.plotdat, age_group_3 == 3), col = "#abd9e9", lwd = 0.5) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.retx2.plotdat, age_group_3 == 4), col = "#2c7bb6", lwd = 0.5) +
  scale_x_continuous(limits = c(0,15), breaks = seq(0,15,1), 
                     labels = c(0, rep("",4), 5, rep("",4),10, rep("",4),15)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "", y = "", title = "Death") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom")  


pw.msm.retx.e4 <- (p1.e4+p4.e4)/(p2.e4+p3.e4)
pw.msm.retx.e4 + plot_annotation(tag_levels = "A")



# Summary of this model (e4) for all age groups (mean time in state)
msm.retx.3.com <- survfit(Surv(tstart, tstop, e4) ~ age_group_3, data=df.ss2, id=PERS_ID)
print(msm.retx.3.com)





# Stacked CIF endpoint e3 by age group - color
msm.retx.2.com <- survfit(Surv(tstart, tstop, e3) ~ age_group_3, data=df.ss2, id=PERS_ID)
oldpar <- par(mai=c(0.5, 0.5, 0.5, 0.5), oma=c(0.5,2,0.5,0.5), mgp=c(2,0.5,0) , mfrow = c(2,2))
stackedCIF(msm.retx.2.com, group = 1,
           fill = c("#CC79A7","#009E73","#E69F00","#000000"),
           xlab="", ylab="Cumulative Incidence", main = "Age 18-54", xlim = c(0,15))
minor.tick(nx=5, ny=2, tick.ratio = 0.5)
stackedCIF(msm.retx.2.com, group = 2,
           fill = c("#CC79A7","#009E73","#E69F00","#000000"),
           xlab="", ylab="", main = "Age 55-64", xlim = c(0,15))
minor.tick(nx=5, ny=2, tick.ratio = 0.5)
stackedCIF(msm.retx.2.com, group = 3,
           fill = c("#CC79A7","#009E73","#E69F00","#000000"),
           xlab="Years after transplant", ylab="Cumulative Incidence", main = "Age 65-74", xlim = c(0,15))
minor.tick(nx=5, ny=2, tick.ratio = 0.5)
stackedCIF(msm.retx.2.com, group = 4,
           fill = c("#CC79A7","#009E73","#E69F00","#000000"),
           xlab="Years after transplant", ylab="", main = "Age \u2265 75", xlim = c(0,15))
minor.tick(nx=5, ny=2, tick.ratio = 0.5)
par(oldpar)


# Legend (will be pasted in with graphics editor) - color
png(file="C:/Users/Thomas Vanhove/Box Sync/R/SRTR/Plots/legend MSM 5 CIF.png",
    width=10, height=10, units = "cm", res = 1200)
plot.new()
legend("topleft", c("Alive with functioning graft", "Death after graft failure",
                    "Death with functioning graft","Retransplanted","Graft failed, alive"), 
       fill = c("white","#000000","#E69F00","#009E73","#CC79A7"), bty = "y", xpd = NA)
dev.off()






###################################################################################
# MSM of states AFTER graft loss: 'alive, graft failed', 'death', 'retransplanted'
df$age_gl <- df$REC_AGE_AT_TX + df$tfl_graft_yrs # create variable age at graft loss (rather than at transplant)
df$age_group_3_gl <- cut(df$age_gl, breaks = c(0, 55, 65, 75, 100)) # And cut into same age_group_3 bins
df$age_group_4_gl <- cut(df$age_gl, breaks = c(0, 30, 40, 50, 60, 70, 100)) # Same with age_group_4

temp <- subset(df, tfl_graft_status == 1) # subset everyone with death censored graft loss
temp$last_fu <- (temp$last_fu - temp$tfl_graft_yrs) # new time zero is time of graft loss, so subtract from other outcome times
temp$retx_yrs <- (temp$retx_yrs - temp$tfl_graft_yrs)
temp$death_comp_yrs <- (temp$death_comp_yrs - temp$tfl_graft_yrs)
temp$tfl_lafu_comp_yrs <- (temp$tfl_lafu_comp_yrs - temp$tfl_graft_yrs)
temp$relist_yrs <- (temp$relist_yrs - temp$tfl_graft_yrs)
temp$tfl_graft_yrs <- 0
temp$last_fu <- ifelse(temp$last_fu == 0, (temp$last_fu + 0.001), temp$last_fu) # Last follow up can't be 0
temp <- subset(temp, !is.na(last_fu)) # And it can't be NA
df.agl <- temp # Will use this dataframe (not in start-stop format) for competing risk analysis later
df.agl$pt_retx <- df.agl$retx # Variable for whether patient was retransplanted, to be used later in subsetting them 
# (as 'retx' will become the event variable in the multi-state dataframe)

df.ss3 <- tmerge(data1 = df.agl[,c("PERS_ID", "DONOR_ID", "REC_TX_DT","REC_AGE_AT_TX","DON_AGE","REC_MM_EQUIV_TX","REC_COLD_ISCH_TM",
                                 "DON_TY","CAN_GENDER", "age_group","age_group_2","age_group_3","age_group_4",
                                 "age_group_3_gl", "age_group_4_gl",
                                 "CAN_DIAB", "CAN_RACE_SRTR","CAN_ETHNICITY_SRTR","REC_BMI",
                                 "CAN_CEREB_VASC","CAN_ANGINA", "relist", 'pt_retx')], 
                  data2 = df.agl, id = PERS_ID, tstop = last_fu)
df.ss3 <- tmerge(data1 = df.ss3, data2 = df.agl, id = PERS_ID,  
                  retx = event(retx_yrs, retx), death_comp = event(death_comp_yrs, death_comp_status))
df.ss3 <- tmerge(df.ss3, df.ss3, id = PERS_ID, enum=cumtdc(tstart))
temp2 <- if_else(df.ss3$retx == 1, 1, if_else(df.ss3$death_comp == 1, 2, 0), missing = 0)
table(temp2)
df.ss3$event <- factor(temp2, 0:2, labels=c("censor", "retx", "death"))

dim(df.ss3)
summary(df.ss3)
attr(df.ss3, "tcount")
with(df.ss3, table(retx, death_comp)) # No ties (no retx and death on same row)

msm.agl <- survfit(Surv(tstart, tstop, event) ~ 1, data=df.ss3, id=PERS_ID) #agl for 'after graft loss'

plot(msm.agl, lty=c(3,1,5), lwd = 2, 
     xlab="Years", ylab="Probability in State", xlim = c(0,15), noplot = NULL, 
     main = "Survival after graft failure")
minor.tick(nx=5, ny=2, tick.ratio = 0.5)
legend("right", c("Failed graft, alive", "Retransplanted", "Death"), 
       lty=c(3,1,5), lwd = 2, bty = "n", xpd = NA)


# First, Kaplan-Meier estimates of death after graft loss
d.agl <- summary(survfit(Surv(tfl_lafu_comp_yrs, death_comp_status==1) ~ age_group_3_gl, data = df.agl), times = c(2))
1-d.agl$surv
1-d.agl$lower
1-d.agl$upper


# Create individual plots msm.agl for the age groups and arrange together
msm.agl.a <- survfit(Surv(tstart, tstop, event) ~ 1, data=df.ss3, subset = df.ss3$age_group_3_gl == "(0,55]", id=PERS_ID)
msm.agl.b <- survfit(Surv(tstart, tstop, event) ~ 1, data=df.ss3, subset = df.ss3$age_group_3_gl == "(55,65]", id=PERS_ID)
msm.agl.c <- survfit(Surv(tstart, tstop, event) ~ 1, data=df.ss3, subset = df.ss3$age_group_3_gl == "(65,75]", id=PERS_ID)
msm.agl.d <- survfit(Surv(tstart, tstop, event) ~ 1, data=df.ss3, subset = df.ss3$age_group_3_gl == "(75,100]", id=PERS_ID)

msm.agl.a # Mean time in states (s0 being 'graft failed, alive') in the 4 age groups
msm.agl.b
msm.agl.c
msm.agl.d

summary(msm.agl.a[,3], times = 2) # 2y mortality after graft loss in youngest age group
summary(msm.agl.b[,3], times = 2)
summary(msm.agl.c[,3], times = 2)
summary(msm.agl.d[,3], times = 2)

# Create dataframe with data to plot
fit.list <- list(msm.agl.a, msm.agl.b, msm.agl.c, msm.agl.d)

p.msm.agl.plotdat <- data.frame(gfa = numeric(), retx = numeric(), death = numeric(), 
                                 time = numeric())
datalist <- list()
for (i in 1:4) {
  tmp <- data.frame(gfa = fit.list[[i]]$pstate[,1], retx = fit.list[[i]]$pstate[,2], 
                    death = fit.list[[i]]$pstate[,3], time = fit.list[[i]]$time,
                    age_group_3 = rep(i, length(fit.list[[i]]$time)))
  datalist[[i]] <- tmp
}
p.msm.agl.plotdat <- bind_rows(datalist)

# Plot
p1.agl <- ggplot() + 
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.agl.plotdat, age_group_3 == 1), col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.agl.plotdat, age_group_3 == 2), col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.agl.plotdat, age_group_3 == 3), col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.agl.plotdat, age_group_3 == 4), col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), 
                     labels = seq(0,10,1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "", y = "Probability in state", title = "Graft failed, alive") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 

p2.agl <- ggplot() + 
  geom_step(aes(x=time, y=retx), data = subset(p.msm.agl.plotdat, age_group_3 == 1), col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.agl.plotdat, age_group_3 == 2), col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.agl.plotdat, age_group_3 == 3), col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.agl.plotdat, age_group_3 == 4), col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), 
                     labels = seq(0,10,1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "Years", y = "Probability in state", title = "Retransplanted") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 


p3.agl <- ggplot() + 
  geom_step(aes(x=time, y=death), data = subset(p.msm.agl.plotdat, age_group_3 == 1), col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.agl.plotdat, age_group_3 == 2), col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.agl.plotdat, age_group_3 == 3), col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.agl.plotdat, age_group_3 == 4), col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), 
                     labels = seq(0,10,1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "Years", y = "Probability in state", title = "Death") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 



pw.agl <- (p1.agl + p2.agl)/(p3.agl + plot_spacer())
pw.agl + plot_annotation(tag_levels = "A")

# Same without the baseline state of 'alive with failed graft'
pw.agl.2 <- (p3.agl + p2.agl)
pw.agl.2 + plot_annotation(tag_levels = "A")



# Single plot of states 'death' and 'retransplanted' after graft failure
msm.agl.com <- survfit(Surv(tstart, tstop, event) ~ age_group_3, data=df.ss3, id=PERS_ID)
print(msm.agl.com)
oldpar <- par(mar = c(5.1, 4.1, 4.1, 18), mfrow=c(1,1))
plot(msm.agl.com, col = rep(c("#d7191c","#fdae61","#abd9e9","#2c7bb6"), times = 2), lty=c(3,3,3,3,1,1,1,1), lwd = 2, 
     xlab="Years", ylab="Probability in State", xlim = c(0,10)) 
minor.tick(nx=2, ny=2, tick.ratio = 0.5)
legend("topleft", c("Death : >75","Death : 66-75", "Death : 56-65","Death : 18-55",
                    "Retransplanted : 18-55", "Retransplanted : 56-65",
                    "Retransplanted : 66-75","Retransplanted : >75"), 
       col = c("#2c7bb6","#abd9e9","#fdae61","#d7191c","#d7191c","#fdae61","#abd9e9","#2c7bb6"), 
       lty=c(1,1,1,1,3,3,3,3), lwd = 2, cex = 0.75,
       bty = "n", inset = c(1,0), xpd = NA)
par(oldpar)






# Stacked CIF outcome after graft failure by age group - black and white
oldpar <- par(mai=c(0.5, 0.5, 0.5, 0.5), oma=c(0.5,2,0.5,0.5), mgp=c(2,0.5,0) , mfrow = c(2,2))
stackedCIF(msm.agl.com, group = 1,
           fill = c("gray80", "gray40"),
           xlab="", ylab="Stacked Cumulative Incidence", main = "Age 18-54", xlim = c(0,10))
minor.tick(nx=2, ny=2, tick.ratio = 0.5)
stackedCIF(msm.agl.com, group = 2,
           fill = c("gray80", "gray40"),
           xlab="", ylab="", main = "Age 55-64", xlim = c(0,10))
minor.tick(nx=2, ny=2, tick.ratio = 0.5)
stackedCIF(msm.agl.com, group = 3,
           fill = c("gray80", "gray40"),
           xlab="Years", ylab="Stacked Cumulative Incidence", main = "Age 65-74", xlim = c(0,10))
minor.tick(nx=2, ny=2, tick.ratio = 0.5)
stackedCIF(msm.agl.com, group = 4,
           fill = c("gray80", "gray40"),
           xlab="Years", ylab="", main = "Age \u2265 75", xlim = c(0,10))
minor.tick(nx=2, ny=2, tick.ratio = 0.5)
par(oldpar)

# Legend (will be pasted in with graphics editor) - black and white
plot.new()
legend("topleft", c("Graft failed, alive", "Death after graft failure", "Retransplanted"), 
       fill = c("white", "gray40", "gray80"), bty = "n", xpd = NA)





# Create same stackedCIF plots, but with age groups as 'age at graft loss'
msm.agl.agegl <- survfit(Surv(tstart, tstop, event) ~ age_group_3_gl, data=df.ss3, id=PERS_ID) 
oldpar <- par(mai=c(0.5, 0.5, 0.5, 0.5), oma=c(0.5,2,0.5,0.5), mgp=c(2,0.5,0) , mfrow = c(2,2))
stackedCIF(msm.agl.agegl, group = 1,
           fill = c("gray80", "gray40"),
           xlab="", ylab="Stacked Cumulative Incidence", main = "Age 18-54", xlim = c(0,10))
minor.tick(nx=2, ny=2, tick.ratio = 0.5)
stackedCIF(msm.agl.agegl, group = 2,
           fill = c("gray80", "gray40"),
           xlab="", ylab="", main = "Age 55-64", xlim = c(0,10))
minor.tick(nx=2, ny=2, tick.ratio = 0.5)
stackedCIF(msm.agl.agegl, group = 3,
           fill = c("gray80", "gray40"),
           xlab="Years", ylab="Stacked Cumulative Incidence", main = "Age 65-74", xlim = c(0,10))
minor.tick(nx=2, ny=2, tick.ratio = 0.5)
stackedCIF(msm.agl.agegl, group = 4,
           fill = c("gray80", "gray40"),
           xlab="Years", ylab="", main = "Age \u2265 75", xlim = c(0,10))
minor.tick(nx=2, ny=2, tick.ratio = 0.5)
par(oldpar)





# ADDING STATE FOR 'WAITLISTED AFTER GRAFT LOSS'
df.ss4 <- tmerge(data1 = df.agl[,c("PERS_ID", "DONOR_ID", "REC_TX_DT","REC_AGE_AT_TX","DON_AGE","REC_MM_EQUIV_TX","REC_COLD_ISCH_TM",
                                 "DON_TY","CAN_GENDER", "age_group","age_group_2","age_group_3","age_group_4",
                                 "age_group_3_gl", "age_group_4_gl",
                                 "CAN_DIAB", "CAN_RACE_SRTR","CAN_ETHNICITY_SRTR","REC_BMI",
                                 "CAN_CEREB_VASC","CAN_ANGINA")], 
                 data2 = df.agl, id = PERS_ID, tstop = last_fu)
df.ss4 <- tmerge(data1 = df.ss4, data2 = df.agl, id = PERS_ID,  
                 retx = event(retx_yrs, retx), death_comp = event(death_comp_yrs, death_comp_status),
                 relist = event(relist_yrs, relist))
df.ss4 <- tmerge(df.ss4, df.ss4, id = PERS_ID, enum=cumtdc(tstart))
temp2 <- if_else(df.ss4$relist == 1, 1, if_else(df.ss4$retx == 1, 2, if_else(df.ss4$death_comp == 1, 3, 0, missing = 0), 
                                                missing = 0), missing = 0)
table(temp2)
df.ss4$event <- factor(temp2, 0:3, labels=c("censor", "relist", "retx", "death"))

dim(df.ss4)
summary(df.ss4)
attr(df.ss4, "tcount")

msm.relist <- survfit(Surv(tstart, tstop, event) ~ 1, data=df.ss4, id=PERS_ID) 
plot(msm.relist)
print(msm.relist)


# Plot 3 states (relisted, retransplanted, death) by 4 age groups
# Below is with age groups AT GRAFT LOSS (could do same with age at transplant by using age_group_3)
msm.relist.a <- survfit(Surv(tstart, tstop, event) ~ 1, data=df.ss4, subset = df.ss4$age_group_3_gl == "(0,55]", id=PERS_ID)
msm.relist.b <- survfit(Surv(tstart, tstop, event) ~ 1, data=df.ss4, subset = df.ss4$age_group_3_gl == "(55,65]", id=PERS_ID)
msm.relist.c <- survfit(Surv(tstart, tstop, event) ~ 1, data=df.ss4, subset = df.ss4$age_group_3_gl == "(65,75]", id=PERS_ID)
msm.relist.d <- survfit(Surv(tstart, tstop, event) ~ 1, data=df.ss4, subset = df.ss4$age_group_3_gl == "(75,100]", id=PERS_ID)

fit.list <- list(msm.relist.a, msm.relist.b, msm.relist.c, msm.relist.d)

p.msm.relist.plotdat <- data.frame(gfa = numeric(), relist = numeric(), retx = numeric(), death = numeric(), 
                                time = numeric())
datalist <- list()
for (i in 1:4) {
  tmp <- data.frame(gfa = fit.list[[i]]$pstate[,1], relist = fit.list[[i]]$pstate[,2], 
                    retx = fit.list[[i]]$pstate[,3], 
                    death = fit.list[[i]]$pstate[,4], time = fit.list[[i]]$time,
                    age_group_3_gl = rep(i, length(fit.list[[i]]$time)))
  datalist[[i]] <- tmp
}
p.msm.relist.plotdat <- bind_rows(datalist)


p1.agl <- ggplot() + 
  geom_step(aes(x=time, y=relist), data = subset(p.msm.relist.plotdat, age_group_3_gl == 1), col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=relist), data = subset(p.msm.relist.plotdat, age_group_3_gl == 2), col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=relist), data = subset(p.msm.relist.plotdat, age_group_3_gl == 3), col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=relist), data = subset(p.msm.relist.plotdat, age_group_3_gl == 4), col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), 
                     labels = seq(0,10,1)) +
  scale_y_continuous(limits = c(0,0.5), breaks = seq(0,1,0.1)) +
  labs(x = "", y = "Probability in state", title = "Relisted") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 

p2.agl <- ggplot() + 
  geom_step(aes(x=time, y=retx), data = subset(p.msm.relist.plotdat, age_group_3_gl == 1), col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.relist.plotdat, age_group_3_gl == 2), col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.relist.plotdat, age_group_3_gl == 3), col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.relist.plotdat, age_group_3_gl == 4), col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), 
                     labels = seq(0,10,1)) +
  scale_y_continuous(limits = c(0,0.5), breaks = seq(0,1,0.1)) +
  labs(x = "Years after graft failure", y = "", title = "Retransplanted") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 


p3.agl <- ggplot() + 
  geom_step(aes(x=time, y=death), data = subset(p.msm.relist.plotdat, age_group_3_gl == 1), col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.relist.plotdat, age_group_3_gl == 2), col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.relist.plotdat, age_group_3_gl == 3), col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.relist.plotdat, age_group_3_gl == 4), col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), 
                     labels = seq(0,10,1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "Years after graft failure", y = "Probability in state", title = "Death") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 
p.blank <- ggplot() + geom_blank() + theme_minimal() + labs(title = "Model diagram") +
  theme(plot.title = element_text(hjust = 0.5)) 

pw.agl <- (p1.agl + p2.agl)/(p3.agl + p.blank)
pw.agl + plot_annotation(tag_levels = "A")

# ggsave("MSM AGL.tiff", device = "tiff", dpi = 300, width = 16, height = 16, units = "cm", 
       path = "C:/Users/Thomas Vanhove/Box Sync/R/SRTR/Plots")



# Same in stackedCIF plots, age groups as 'age at graft loss'
msm.relist.agegl <- survfit(Surv(tstart, tstop, event) ~ age_group_3_gl, data=df.ss4, id=PERS_ID) 
oldpar <- par(mai=c(0.5, 0.5, 0.5, 0.5), oma=c(0.5,2,0.5,0.5), mgp=c(2,0.5,0) , mfrow = c(2,2))
stackedCIF(msm.relist.agegl, group = 1,
           fill = c("gray40", "gray80", "gray60"),
           xlab="", ylab="Stacked Cumulative Incidence", main = "Age 18-54", xlim = c(0,10))
minor.tick(nx=2, ny=2, tick.ratio = 0.5)
stackedCIF(msm.relist.agegl, group = 2,
           fill = c("gray40", "gray80", "gray60"),
           xlab="", ylab="", main = "Age 55-64", xlim = c(0,10))
minor.tick(nx=2, ny=2, tick.ratio = 0.5)
stackedCIF(msm.relist.agegl, group = 3,
           fill = c("gray40", "gray80", "gray60"),
           xlab="Years after graft failure", ylab="Stacked Cumulative Incidence", main = "Age 65-74", xlim = c(0,10))
minor.tick(nx=2, ny=2, tick.ratio = 0.5)
stackedCIF(msm.relist.agegl, group = 4,
           fill = c("gray40", "gray80", "gray60"),
           xlab="Years after graft failure", ylab="", main = "Age \u2265 75", xlim = c(0,10))
minor.tick(nx=2, ny=2, tick.ratio = 0.5)
par(oldpar)





# MSM looking only at relisted patients (outcomes retx and death)
msm.relist.a <- survfit(Surv(tstart, tstop, event) ~ 1, data=subset(df.ss3, age_group_3_gl == "(0,55]" & relist == 1), id=PERS_ID)
msm.relist.b <- survfit(Surv(tstart, tstop, event) ~ 1, data=subset(df.ss3, age_group_3_gl == "(55,65]" & relist == 1), id=PERS_ID)
msm.relist.c <- survfit(Surv(tstart, tstop, event) ~ 1, data=subset(df.ss3, age_group_3_gl == "(65,75]" & relist == 1), id=PERS_ID)
msm.relist.d <- survfit(Surv(tstart, tstop, event) ~ 1, data=subset(df.ss3, age_group_3_gl == "(75,100]" & relist == 1), id=PERS_ID)

msm.relist.a # Mean time in states (s0 being 'graft failed, alive') in the 4 age groups
msm.relist.b
msm.relist.c
msm.relist.d

# Mean time in state retransplanted
survfit(Surv(tstart, tstop, event) ~ age_group_3_gl, data=subset(df.ss3, pt_retx == 1), id=PERS_ID) 

# Mortality after graft loss with CI
summary(msm.relist.a[,3], times = 2)
summary(msm.relist.b[,3], times = 2)
summary(msm.relist.c[,3], times = 2)
summary(msm.relist.d[,3], times = 2)

# create data for plots
fit.list <- list(msm.relist.a, msm.relist.b, msm.relist.c, msm.relist.d)

p.msm.relist.plotdat <- data.frame(gfa = numeric(), retx = numeric(), death = numeric(), 
                                time = numeric())
datalist <- list()
for (i in 1:4) {
  tmp <- data.frame(gfa = fit.list[[i]]$pstate[,1], retx = fit.list[[i]]$pstate[,2], 
                    death = fit.list[[i]]$pstate[,3], time = fit.list[[i]]$time,
                    age_group_3_gl = rep(i, length(fit.list[[i]]$time)))
  datalist[[i]] <- tmp
}
p.msm.relist.plotdat <- bind_rows(datalist)


p1.relist <- ggplot() + 
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.relist.plotdat, age_group_3_gl == 1), col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.relist.plotdat, age_group_3_gl == 2), col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.relist.plotdat, age_group_3_gl == 3), col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.relist.plotdat, age_group_3_gl == 4), col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), 
                     labels = seq(0,10,1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "", y = "Probability in state", title = "Graft failed, alive") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 

p2.relist <- ggplot() + 
  geom_step(aes(x=time, y=retx), data = subset(p.msm.relist.plotdat, age_group_3_gl == 1), col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.relist.plotdat, age_group_3_gl == 2), col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.relist.plotdat, age_group_3_gl == 3), col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.relist.plotdat, age_group_3_gl == 4), col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), 
                     labels = seq(0,10,1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "Years after graft failure", y = "", title = "Retransplanted") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 


p3.relist <- ggplot() + 
  geom_step(aes(x=time, y=death), data = subset(p.msm.relist.plotdat, age_group_3_gl == 1), col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.relist.plotdat, age_group_3_gl == 2), col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.relist.plotdat, age_group_3_gl == 3), col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.relist.plotdat, age_group_3_gl == 4), col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), 
                     labels = seq(0,10,1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "Years after graft failure", y = "Probability in state", title = "Death") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 



pw.relist <- (p1.relist + p2.relist)/(p3.relist + plot_spacer())
pw.relist + plot_annotation(tag_levels = "A")


# Same in patients who were never relisted
# MSM looking only at relisted patients (outcomes retx and death)
msm.nolist.a <- survfit(Surv(tstart, tstop, event) ~ 1, data=subset(df.ss3, age_group_3_gl == "(0,55]" & relist == 0), id=PERS_ID)
msm.nolist.b <- survfit(Surv(tstart, tstop, event) ~ 1, data=subset(df.ss3, age_group_3_gl == "(55,65]" & relist == 0), id=PERS_ID)
msm.nolist.c <- survfit(Surv(tstart, tstop, event) ~ 1, data=subset(df.ss3, age_group_3_gl == "(65,75]" & relist == 0), id=PERS_ID)
msm.nolist.d <- survfit(Surv(tstart, tstop, event) ~ 1, data=subset(df.ss3, age_group_3_gl == "(75,100]" & relist == 0), id=PERS_ID)

msm.nolist.a # Mean time in states (s0 being 'graft failed, alive') in the 4 age groups
msm.nolist.b
msm.nolist.c
msm.nolist.d

summary(msm.nolist.a[,3], times = 2)
summary(msm.nolist.b[,3], times = 2)
summary(msm.nolist.c[,3], times = 2)
summary(msm.nolist.d[,3], times = 2)

fit.list <- list(msm.nolist.a, msm.nolist.b, msm.nolist.c, msm.nolist.d)

p.msm.nolist.plotdat <- data.frame(gfa = numeric(), retx = numeric(), death = numeric(), 
                                   time = numeric())
datalist <- list()
for (i in 1:4) {
  tmp <- data.frame(gfa = fit.list[[i]]$pstate[,1], retx = fit.list[[i]]$pstate[,2], 
                    death = fit.list[[i]]$pstate[,3], time = fit.list[[i]]$time,
                    age_group_3_gl = rep(i, length(fit.list[[i]]$time)))
  datalist[[i]] <- tmp
}
p.msm.nolist.plotdat <- bind_rows(datalist)


p1.nolist <- ggplot() + 
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 1), col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 2), col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 3), col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=gfa), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 4), col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), 
                     labels = seq(0,10,1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "", y = "Probability in state", title = "Graft failed, alive") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 

p2.nolist <- ggplot() + 
  geom_step(aes(x=time, y=retx), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 1), col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 2), col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 3), col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 4), col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), 
                     labels = seq(0,10,1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "Years after graft failure", y = "", title = "Retransplanted") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 


p3.nolist <- ggplot() + 
  geom_step(aes(x=time, y=death), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 1), col = "#d7191c", lwd = 0.7) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 2), col = "#fdae61", lwd = 0.7) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 3), col = "#abd9e9", lwd = 0.7) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 4), col = "#2c7bb6", lwd = 0.7) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), 
                     labels = seq(0,10,1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "Years after graft failure", y = "Probability in state", title = "Death") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 



pw.nolist <- (p1.nolist + p2.nolist)/(p3.nolist + plot_spacer())
pw.nolist + plot_annotation(tag_levels = "A")



# Plot of death after graft loss in relisted and non-relisted patients
p.list.d <- ggplot() + 
  geom_step(aes(x=time, y=death), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 1), col = "#d7191c", lwd = 1, lty = 2) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 2), col = "#fdae61", lwd = 1, lty = 2) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 3), col = "#abd9e9", lwd = 1, lty = 2) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.nolist.plotdat, age_group_3_gl == 4), col = "#2c7bb6", lwd = 1, lty = 2) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.relist.plotdat, age_group_3_gl == 1), col = "#d7191c", lwd = 1, lty = 1) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.relist.plotdat, age_group_3_gl == 2), col = "#fdae61", lwd = 1, lty = 1) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.relist.plotdat, age_group_3_gl == 3), col = "#abd9e9", lwd = 1, lty = 1) +
  geom_step(aes(x=time, y=death), data = subset(p.msm.relist.plotdat, age_group_3_gl == 4), col = "#2c7bb6", lwd = 1, lty = 1) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), 
                     labels = seq(0,10,1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "", y = "Probability in state", title = "Death") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 

# Plot of retx after graft loss in relisted patients
p.list.retx <- ggplot() + 
  geom_step(aes(x=time, y=retx), data = subset(p.msm.relist.plotdat, age_group_3_gl == 1), col = "#d7191c", lwd = 1, lty = 1) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.relist.plotdat, age_group_3_gl == 2), col = "#fdae61", lwd = 1, lty = 1) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.relist.plotdat, age_group_3_gl == 3), col = "#abd9e9", lwd = 1, lty = 1) +
  geom_step(aes(x=time, y=retx), data = subset(p.msm.relist.plotdat, age_group_3_gl == 4), col = "#2c7bb6", lwd = 1, lty = 1) +
  scale_x_continuous(limits = c(0,10), breaks = seq(0,10,1), 
                     labels = seq(0,10,1)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.1)) +
  labs(x = "Years after graft failure", y = "Probability in state", title = "Retransplanted") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5)) 

p.blank <- ggplot() + geom_blank() + theme_minimal() + labs(title = "Model diagram") +
  theme(plot.title = element_text(hjust = 0.5)) 

# combine
(p.list.d + plot_spacer())/(p.list.retx + p.blank) + plot_annotation(tag_levels = "A")


# ggsave("relist or no.tiff", device = "tiff", dpi = 300, width = 20, height = 16, units = "cm", 
       path = "C:/Users/Thomas Vanhove/Box Sync/R/SRTR/Plots")


# Legend
png(file="C:/Users/Thomas Vanhove/Box Sync/R/SRTR/Plots/legend relisted.png",
    width=10, height=10, units = "cm", res = 1200)
plot.new()
legend("topleft", c("\u2265 75","65 - 74 ","55 - 64","18 - 54"),
       col = c("#2c7bb6","#abd9e9","#fdae61","#d7191c"), 
       lty=1, lwd = 2, cex = 1,
       bty = "y", inset = c(0.5,0), xpd = NA,
       title = "Relisted", title.adj = 0.5)
dev.off()

png(file="C:/Users/Thomas Vanhove/Box Sync/R/SRTR/Plots/legend never relisted.png",
    width=10, height=10, units = "cm", res = 1200)
plot.new()
legend("topleft", c("\u2265 75","65 - 74 ","55 - 64","18 - 54"),
       col = c("#2c7bb6","#abd9e9","#fdae61","#d7191c"), 
       lty=2, lwd = 2, cex = 1,
       bty = "y", inset = c(0.5,0), xpd = NA,
       title = "Never relisted", title.adj = 0.5)
dev.off()


# Summary of numbers of older patients relisted and retransplanted
summary(subset(df, ki_status == 1)$age_group_3_gl)
summary(subset(df, ki_status == 1 & relist == 1)$age_group_3_gl)
summary(subset(df, ki_status == 1 & relist == 0)$age_group_3_gl)
summary(subset(df, ki_status == 1 & retx == 1)$age_group_3_gl)
summary(subset(df, ki_status == 1 & relist == 1 & retx == 1)$age_group_3_gl)
summary(subset(df, ki_status == 1 & relist == 0 & retx == 1)$age_group_3_gl)

# Patient survival after graft loss in different subgroups
summary(survfit(Surv(tfl_lafu_comp_yrs, tfl_lafu_comp_status == 'D') ~ age_group_3, 
                data = subset(df.agl, ki_status == 1 & relist == 0)), times = c(1,2,3,5,10))
























###################################################################################
# ADDENDUM: ADDITIONAL PLOTS

# Recipient age by 5y period
ggplot(transform(df, contrib.bins=cut(REC_TX_DT, breaks = "5 years"), include.lowest=T)) + 
  geom_violin(aes(x = contrib.bins, y = REC_AGE_AT_TX)) +
  geom_boxplot(aes(x = contrib.bins, y = REC_AGE_AT_TX), width = 0.15) +
  scale_y_continuous(breaks = seq(20, 90, 5)) +
  labs(x = "Tx date", y = "Age")


# Recipient age by 5y period
ggplot() + 
  geom_violin(aes(x = contrib.bins, y = REC_AGE_AT_TX), 
              data =  transform(df, contrib.bins=cut(REC_TX_DT, breaks = "5 years"), include.lowest=T), scale = "count") +
  geom_boxplot(aes(x = contrib.bins, y = REC_AGE_AT_TX), width = 0.1, 
               data =  transform(df, contrib.bins=cut(REC_TX_DT, breaks = "5 years"), include.lowest=T), outlier.shape = NA) +
  geom_hline(yintercept = 65, linetype = 2) +
  scale_y_continuous(breaks = seq(20, 90, 5)) +
  scale_x_discrete(labels = c("1990-1994", "1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2019")) +
  labs(x = "Year of transplant", y = "Recipient age") + 
  theme_bw() + 
  theme(axis.title.x = element_text(margin = margin(10,0,0,0)), axis.title.y = element_text(margin = margin(0,10,0,0)))  

summary(df)

# Recipient age at graft loss by 5y period (!required subsetting a different part of the original dataframe)
ggplot() + 
  geom_violin(aes(x = contrib.bins, y = age_gl), 
              data =  transform(df, contrib.bins=cut(TFL_GRAFT_DT, breaks = "5 years"), include.lowest=T), scale = "count") +
  geom_boxplot(aes(x = contrib.bins, y = age_gl), width = 0.1, 
               data =  transform(df, contrib.bins=cut(TFL_GRAFT_DT, breaks = "5 years"), include.lowest=T), outlier.shape = NA) +
  geom_hline(yintercept = 65, linetype = 2) +
  scale_y_continuous(breaks = seq(20, 90, 5)) +
  scale_x_discrete(labels = c("1990-1994", "1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2019")) +
  labs(x = "Year of graft failure", y = "Recipient age at graft failure") + 
  theme_bw() + 
  theme(axis.title.x = element_text(margin = margin(10,0,0,0)), axis.title.y = element_text(margin = margin(0,10,0,0)))


# Recipient age at death by 5y period (!required subsetting a different part of the original dataframe)
df$age_death <- df$REC_AGE_AT_TX + df$death_comp_yrs # create variable 'age at death'
ggplot() + 
  geom_violin(aes(x = contrib.bins, y = age_death), 
              data =  transform(subset(df, !is.na(death_comp_dt)), contrib.bins=cut(death_comp_dt, breaks = "5 years"), 
                                        include.lowest=T), scale = "count") +
  geom_boxplot(aes(x = contrib.bins, y = age_death), width = 0.1, 
               data =  transform(subset(df, !is.na(death_comp_dt)), contrib.bins=cut(death_comp_dt, breaks = "5 years"), 
                                         include.lowest=T), outlier.shape = NA) +
  geom_hline(yintercept = 65, linetype = 2) +
  scale_y_continuous(breaks = seq(20, 90, 5)) +
  scale_x_discrete(labels = c("1990-1994", "1995-1999", "2000-2004", "2005-2009", "2010-2014", "2015-2019")) +
  labs(x = "Year of death", y = "Recipient age at death") + 
  theme_bw() + 
  theme(axis.title.x = element_text(margin = margin(10,0,0,0)), axis.title.y = element_text(margin = margin(0,10,0,0)))


# Donor age by 5y period
ggplot(transform(df, contrib.bins=cut(REC_TX_DT, breaks = "5 years"), include.lowest=T)) + 
  geom_violin(aes(x = contrib.bins, y = DON_AGE)) +
  labs(x = "Tx date", y = "Donor age")


# Time of last follow up
ggplot(data = df, aes(x = tfl_lafu_comp_yrs)) +
  stat_ecdf() + 
  scale_y_continuous(breaks = seq(0,1,0.1)) +
  scale_x_continuous(breaks = seq(0,20,1), labels = c(0,rep("",4), 5,rep("",4),10,rep("",4),15,rep("",4),20)) +
  labs(x = "Years", y = "Cumulative probability", title = "Time of last recorded follow up") +
  theme_bw() + 
  theme(axis.title.x = element_text(margin = margin(10,0,0,0), size = 12), 
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 12),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

# Time of last follow up by age group
ggplot() +
  stat_ecdf(aes(x = tfl_lafu_comp_yrs), data = subset(df, age_group_3 == "(0,55]"), color = "#d7191c") + 
  stat_ecdf(aes(x = tfl_lafu_comp_yrs), data = subset(df, age_group_3 == "(55,65]"), color = "#fdae61") +
  stat_ecdf(aes(x = tfl_lafu_comp_yrs), data = subset(df, age_group_3 == "(65,75]"), color = "#abd9e9") +
  stat_ecdf(aes(x = tfl_lafu_comp_yrs), data = subset(df, age_group_3 == "(75,100]"), color = "#2c7bb6") +
  scale_y_continuous(breaks = seq(0,1,0.1)) +
  scale_x_continuous(breaks = seq(0,20,1), labels = c(0,rep("",4), 5,rep("",4),10,rep("",4),15,rep("",4),20)) +
  labs(x = "Years", y = "Cumulative probability", title = "Time of last recorded follow up") +
  theme_bw() + 
  theme(axis.title.x = element_text(margin = margin(10,0,0,0), size = 12), 
        axis.title.y = element_text(margin = margin(0,10,0,0), size = 12),
        axis.text = element_text(size = 12),
        plot.title = element_text(hjust = 0.5, size = 15),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"),
        panel.grid.minor = element_blank())




# MSM plot where there is a state for 'lost to follow-up'
# First, create event variable for lost to FU, corresponding time is still tfl_lafu_comp_yrs
df$lostfu <- factor(if_else(df$tfl_lafu_comp_status == "L" | df$tfl_lafu_comp_status == "N", 1, 0, missing = 1))

df.ss.lost <- tmerge(data1 = df[,c("PERS_ID", "DONOR_ID", "REC_TX_DT","REC_AGE_AT_TX","DON_AGE","REC_MM_EQUIV_TX","REC_COLD_ISCH_TM",
                              "DON_TY","CAN_GENDER", "age_group","age_group_2","age_group_3","age_group_4",
                              "CAN_DIAB", "CAN_RACE_SRTR","CAN_ETHNICITY_SRTR","REC_BMI",
                              "CAN_CEREB_VASC","CAN_ANGINA")], 
                data2 = df, id = PERS_ID, tstop = last_fu)
df.ss.lost <- tmerge(data1 = df.ss.lost, data2 = df, id = PERS_ID, lostfu = event(tfl_lafu_comp_yrs, lostfu),
                     death_comp = event(death_comp_yrs, death_comp_status),
                     retx = event(retx_yrs, retx),
                     graftloss = event(tfl_graft_yrs, tfl_graft_status))
df.ss.lost <- tmerge(df.ss.lost, df.ss.lost, id = PERS_ID, enum=cumtdc(tstart))
temp <- if_else(df.ss.lost$lostfu == 1, 1, 
                if_else(df.ss.lost$death_comp == 1, 2, 
                        if_else(df.ss.lost$retx == 1, 3,
                                if_else(df.ss.lost$graftloss == 1, 4, 0))))
df.ss.lost$event <- factor(temp, 0:4, labels=c("censor", "lostfu", "death", "retx", "gl")) # 'event' = 5 possible outcomes
df.ss.lost$tstate <- with(df.ss.lost, tstop - tstart)


dim(df.ss.lost)
summary(df.ss.lost)
attr(df.ss.lost, "tcount")
with(df.ss.lost, table(graftloss, death_comp)) # No ties (no graft loss and death on same row)

msm.lost <- survfit(Surv(tstart, tstop, event) ~ age_group_4, data=df.ss.lost, id=PERS_ID)
msm.lost <- survfit(Surv(tstart, tstop, event) ~ 1, data=df.ss.lost, subset = df.ss.lost$age_group_3 == "(0,55]", id=PERS_ID)
print(msm.lost)
plot(msm.lost)


oldpar <- par(mai=c(0.5, 0.5, 0.5, 0.5), oma=c(0.5,2,0.5,0.5), mgp=c(2,0.5,0) , mfrow = c(3,2))
stackedCIF(msm.lost, 
           group = 1,
           lwd = 1, 
           fill = c("#0072B2","#E69F00","#009E73","#CC79A7"),
           xlab="", 
           ylab="Cumulative incidence",
           main = "Age 18 - 29",
           xlim = c(0,15), ylim = c(0,1), las = 1)
minor.tick(nx=5, ny=2, tick.ratio = 0.5)

stackedCIF(msm.lost, 
           group = 2,
           lwd = 1, 
           fill = c("#0072B2","#E69F00","#009E73","#CC79A7"),
           xlab="", 
           ylab="Cumulative incidence",
           main = "Age 30 - 39",
           xlim = c(0,15), ylim = c(0,1), las = 1)
minor.tick(nx=5, ny=2, tick.ratio = 0.5)

stackedCIF(msm.lost, 
           group = 3,
           lwd = 1, 
           fill = c("#0072B2","#E69F00","#009E73","#CC79A7"),
           xlab="", 
           ylab="Cumulative incidence",
           main = "Age 40 - 49",
           xlim = c(0,15), ylim = c(0,1), las = 1)
minor.tick(nx=5, ny=2, tick.ratio = 0.5)

stackedCIF(msm.lost, 
           group = 4,
           lwd = 1, 
           fill = c("#0072B2","#E69F00","#009E73","#CC79A7"),
           xlab="", 
           ylab="",
           main = "Age 50 - 59",
           xlim = c(0,15), ylim = c(0,1), las = 1)
minor.tick(nx=5, ny=2, tick.ratio = 0.5)

stackedCIF(msm.lost, 
           group = 5,
           lwd = 1, 
           fill = c("#0072B2","#E69F00","#009E73","#CC79A7"),
           xlab="Years after transplant", 
           ylab="Cumulative incidence",
           main = "Age 60 - 69",
           xlim = c(0,15), ylim = c(0,1), las = 1)
minor.tick(nx=5, ny=2, tick.ratio = 0.5)

stackedCIF(msm.lost, 
           group = 6,
           lwd = 1, 
           fill = c("#0072B2","#E69F00","#009E73","#CC79A7"),
           xlab="Years after transplant", 
           ylab="",
           main = "Age \u2265 70",
           xlim = c(0,15), ylim = c(0,1), las = 1)
minor.tick(nx=5, ny=2, tick.ratio = 0.5)
par(oldpar)


# Legend 
png(file="C:/Users/Thomas Vanhove/Box Sync/R/SRTR/Plots/legend MSM lostfu.png",
    width=10, height=10, units = "cm", res = 1200)
plot.new()
legend("topleft", c("Alive with functioning graft","Graft failed, alive","Retransplanted","Death","Lost to follow-up"), 
       fill = c("white","#CC79A7","#009E73","#E69F00","#0072B2"), bty = "y", xpd = NA)
dev.off()






# Percent stacked barchart of proportion living vs deceased donor by recipient age
ggplot(df) +
  geom_bar(aes(x=REC_AGE_AT_TX, fill = DON_TY), position = "fill") + 
  scale_x_continuous(breaks = seq(20,90,10), limits = c(18,91)) + 
  scale_y_continuous(breaks = seq(0,1,0.1), labels = scales::percent_format(accuracy = 1)) + 
  labs(x = "Age at transplant", y = "") + 
  theme_bw() +
  scale_fill_manual(name = "Donor type", labels = c("Deceased", "Living"), values = c("grey90", "grey60")) +
  theme(panel.grid.major = element_blank(),
           panel.grid.minor = element_blank())
# ggsave("don type by age.tiff", device = "tiff", dpi = 300, width = 20, height = 16, units = "cm", 
       path = "C:/Users/Thomas Vanhove/Box Sync/R/SRTR/Plots")


# Percent stacked barchart of proportion race by recipient age
ggplot(subset(df, !is.na(CAN_RACE_SRTR))) +
  geom_bar(aes(x=REC_AGE_AT_TX, fill = CAN_RACE_SRTR), position = "fill") + 
  scale_x_continuous(breaks = seq(20,90,10), limits = c(18,90)) + 
  scale_y_continuous(breaks = seq(0,1,0.1), labels = scales::percent_format(accuracy = 1)) + 
  labs(x = "Age at transplant", y = "") + 
  theme_bw() +
  scale_fill_manual(name = "Recipient race", labels = c("White", "Black", "Asian", "Other"),
                    values = c("grey90","black","grey80","grey60")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# ggsave("rec race by age.tiff", device = "tiff", dpi = 300, width = 20, height = 16, units = "cm", 
       path = "C:/Users/Thomas Vanhove/Box Sync/R/SRTR/Plots")


# Percent stacked barchart of sex by recipient age
ggplot(df) +
  geom_bar(aes(x=REC_AGE_AT_TX, fill = CAN_GENDER), position = "fill") + 
  scale_x_continuous(breaks = seq(20,90,10), limits = c(18,90)) + 
  scale_y_continuous(breaks = seq(0,1,0.1), labels = scales::percent_format(accuracy = 1)) + 
  labs(x = "Age at transplant", y = "") + 
  theme_bw() +
  scale_fill_manual(name = "Sex", labels = c("Female", "Male"),
                    values = c("grey90","grey60")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
# ggsave("sex by age.tiff", device = "tiff", dpi = 300, width = 20, height = 16, units = "cm", 
       path = "C:/Users/Thomas Vanhove/Box Sync/R/SRTR/Plots")

# Percent stacked barchart of proportion relisted by age graft loss
ggplot(subset(df.agl, !is.na(relist))) + 
  geom_bar(aes(x=cut(age_gl, breaks = c(0,30,40,50,60,70,80,100), 
                     labels = c("18-29","30-39","40-49","50-59","60-69","70-79","\u2265 80")), 
               fill = relist), position = "fill", color = "black") +
  scale_y_continuous(breaks = seq(0,1,0.1), labels = scales::percent_format(accuracy = 1)) + 
  labs(x = "Age at graft failure", y = "", fill = "") +
  theme_bw() +
  scale_fill_manual(values = c("grey60", "grey80"), labels = c("Never relisted", "Relisted")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())



