# ------------------------------------------------------------------------------------ #
# Program:    PNAS_Social_Mobility_public_v2.R                                         #  
# Author:     Damien Morris                                                            #
# Date:       29 09 2020                                                               #  
# Re-analaysis of twin and social mobility data in Engzell and Tropf 2019              #
# ------------------------------------------------------------------------------------ #

# ---------------------------------------CONTENTS------------------------------------- #

#--------------PART 1: INITIAL SETUP------------------#                

#1.  CLEAR WORKSPACE AND LOAD PACKAGES

#--------------PART 2: PREPARE DATA-------------------#               

#2.  IMPORT GDIM SOCIAL MOBILITY DATA AND PREPARE KEY VARIABLES
#3.  IMPORT BRANIGAN ET AL (2013) TWIN ACE DATA AND PREPARE KEY VARIABLES
#4.  MERGE SOCIAL MOBILITY DATA AND ACE DATA
#5.  CALCULATE INVERSE VARIANCE WEIGHTS FOR PRECISION WEIGHTING ACE ESTIMATES
#6.  SUBSET THE DATA (IN ACCORD WITH E&T METHODOLOGY)
#7.  CREATE ALT_NORWAY DATA AND PERFORM PART 2 AGAIN

#--------------PART 3: REPRODUCE E&T RESULTS----------#                

#8.  REPRODUCE ENGZELL & TROPF TABLE 1 RESULTS (CONTROLLING FOR GENDER IN SMOOTHED RESULTS)
#9.  ALTERNATIVE METHOD OF OBTAINING TABLE 1 RESULTS USING COEFTEST
#10. REPRODUCE ENGZELL & TROPF FIGURES 1B AND 1C

#--------------PART 4: REVISED RESULTS----------------#                

#11. ALT_NORWAY RESULTS. REVISED TABLE 1 RESULTS
#12. ALT_NORWAY RESULTS. REVISED CHARTS
#13. PRECISION-WEIGHTED RESULTS. REVISED CHARTS
#14. PRECISION-WEIGHTED RESULTS. REVISED TABLE 1 DATA
#15. PLOT WITHIN-COUNTRY TRENDS FOR EACH GENDER
#16. CREATE MULTI-PANEL FIGURE FOR LETTER
#17. PREPARE PRECISION-WEIGHTED CHARTS FOR SUPPLEMENTARY MATERIALS
#18. PRODUCE WITHIN-COUNTRY C2 AND E2 CHARTS FOR SUPPLEMENTARY MATERIALS
#19. CREATE MULTI-PANEL FIGURES FOR SUPPLEMENTARY MATERIALS

# ------------------------------------------------------------------------------------ #

# --------------------------------PART 1: INITIAL SETUP------------------------------- #

#1. CLEAR WORKSPACE AND LOAD PACKAGES

#Clear the console (or Control+L)
cat("\014")

#Empty the environment
rm(list=ls())

#Open libraries
library(readstata13) #for reading in stata files
library(miceadds) #linear models with clustered standard errors
library(sandwich) #for vcovCL analysis to corroborate clustered SEs
library(ggplot2)  #for plots
library(dplyr) #for pipe operator and mutate function etc
library(lmtest) #for coeftest
library(scales) #for checking colours for manual recolouring of ggplot charts
library(cowplot) #for plot_grid in multipanel figure

#Set working directory
setwd("Enter_your_filepath")

# --------------------------------PART 2: PREPARE DATA-------------------------------- #

#INSTRUCTIONS:

#Run Part 2 of script twice.
# - 1st time: run with alternative Norwegian parent-offspring correlations. Then create Alt_Norway data object and comment out Norwegian substitution code
# - 2nd time: run straight through to create ETmerged_max_exclusions_nozero data for reproducing E&T's main results.

#2. IMPORT GDIM SOCIAL MOBILITY DATA AND PREPARE KEY VARIABLES------------------------------------------

#Import E&T's GDIM social mobility file.
#This can be downloaded from: https://osf.io/qhcjz/
StataGDIM <- read.dta13("Enter_your_filepath")

#Substitute in endogenous Norwegian rIG values from Heath et al (1985) to create Alt_Norway data. Then comment this out and run Part 2 of the script a 2nd time.
StataGDIM[429, "cor"] <- 0.625 #male 1940s Norway  (GDIM value was 0.23)
StataGDIM[441, "cor"] <- 0.530 #male 1950s Norway  (GDIM value was 0.29)
StataGDIM[428, "cor"] <- 0.680 #female 1940s Norway (GDIM value was 0.45)
StataGDIM[440, "cor"] <- 0.575 #female 1950s Norway (GDIM value was 0.33)

#Reorder the GDIM data 
StataGDIM <- StataGDIM[order(StataGDIM$child, StataGDIM$parent, StataGDIM$wbcode, StataGDIM$year),]

#Create lag and lead variable by group for rIG 
StataGDIM <- StataGDIM %>% group_by(child, parent, wbcode) %>% mutate(rIG_lag = lag(cor))
StataGDIM <- StataGDIM %>% group_by(child, parent, wbcode) %>% mutate(rIG_lead = lead(cor))

#Create lag and lead variable by group for sdc 
StataGDIM <- StataGDIM %>% group_by(child, parent, wbcode) %>% mutate(sdc_lag = lag(sdc))
StataGDIM <- StataGDIM %>% group_by(child, parent, wbcode) %>% mutate(sdc_lead = lead(sdc))

#Create smooth variable for rIG and sdc
StataGDIM$rIGsmooth <-rowMeans(StataGDIM[,c("cor", "rIG_lag", "rIG_lead")], na.rm = TRUE)
StataGDIM$sdc_smooth <-rowMeans(StataGDIM[,c("sdc", "sdc_lag", "sdc_lead")], na.rm = TRUE)

#Rename child variable to gender variable in StataGDIM 
StataGDIM <- StataGDIM %>% rename(gender = child)

#In renamed "gender" variable replace "son" with "Male" etc.
StataGDIM$gender <- gsub('son', 'Male', StataGDIM$gender)
StataGDIM$gender <- gsub('daughter', 'Female', StataGDIM$gender)
StataGDIM$gender <- gsub('all', 'Mixed', StataGDIM$gender)

#3.  IMPORT BRANIGAN ET AL (2013) TWIN ACE DATA AND PREPARE KEY VARIABLES

#Import E&T's version of Branigan et al (2013) twin ACE data.
#This can be downloaded from: https://osf.io/fwqa5/
StataACE <- read.dta13("Enter_your_filepath")

#Create dob_midrange variable
StataACE <- StataACE %>% mutate(dob_midrange = (yr_max + yr_min)/2)

#Create cohort (decade of birth) variable in ACE data
StataACE <- StataACE %>% mutate(cohort = as.integer(dob_midrange/10)*10)

#Strip out all cohorts with no max year before 1940
StataACE <-subset(StataACE, yr_max >= 1940)

#Recode all remaining sub_1940 cohorts as 1940
StataACE$cohort <- gsub(1920, 1940, StataACE$cohort)
StataACE$cohort <- gsub(1930, 1940, StataACE$cohort)

#4. MERGE SOCIAL MOBILITY DATA AND ACE DATA

ETmerged <- merge(StataACE,StataGDIM,by=c("countryname","gender", "cohort"))

#Create cohort ID by concatenating two character vectors
ETmerged$cohortID <- paste(ETmerged$wbcode, ETmerged$cohort, sep = "", collapse = NULL)

#Create short cohort ID by concatenating two character vectors.
#cohortIDshort is created to reduce clutter in multipanel figures
ETmerged$cohort_short <- substring(ETmerged$cohort, 3)
ETmerged$cohortIDshort <- paste(ETmerged$wbcode, ETmerged$cohort_short, sep = "", collapse = NULL) 

#Create de-standardised and smoothed ACE values
ETmerged$h2primesmooth <- ETmerged$h2 * (ETmerged$sdc_smooth)^2 
ETmerged$c2primesmooth <- ETmerged$c2 * ETmerged$sdc_smooth^2
ETmerged$e2primesmooth <- ETmerged$e2 * ETmerged$sdc_smooth^2

#5. CALCULATE INVERSE VARIANCE WEIGHTS FOR PRECISION WEIGHTING ACE ESTIMATES

#Sampling variance of twin correlations (method from Borenstein et al 2009 Intro to Meta-analysis)
ETmerged<-ETmerged %>% mutate (rMZvar = (1-r_mz^2)^2/(n_mz-1))
ETmerged<-ETmerged %>% mutate (rDZvar = (1-r_dz^2)^2/(n_dz-1))

#Sampling variance of ACE variance components (method from Branigan et al, 2013)
ETmerged<-ETmerged %>% mutate (h2var = (4*(rMZvar+rDZvar)))
ETmerged<-ETmerged %>% mutate (c2var = (4*rDZvar+rMZvar))
ETmerged<-ETmerged %>% mutate (e2var = rMZvar)

#Inverse variance weightings for variance of h2, c2 and e2
ETmerged<-ETmerged %>% mutate (h2weight = 1/h2var)
ETmerged<-ETmerged %>% mutate (c2weight = 1/c2var)
ETmerged<-ETmerged %>% mutate (e2weight = 1/e2var)

#Remove excess variables which are no longer needed.
ETmerged <- select(ETmerged,-c(year, yr_min, yr_max, wbcode, survey, ginip, ginic, rIG_lag, rIG_lead, sdc_lag, sdc_lead, rMZvar, rDZvar, h2var,c2var, e2var, cohort_short)) 

#6. SUBSET THE DATA (IN ACCORD WITH E&T METHODOLOGY)

#Subset data for correlations with highest parent only. The main analysis in Engzell Tropf uses this particular rIG subset.
ETmerged_max <-subset(ETmerged, parent == "max")

#Remove excluded studies (i.e. Minnesota, SRI, Vietnam Vets)
ETmerged_max_exclusions <- subset(ETmerged_max, sample != "Minnesota")
ETmerged_max_exclusions <- subset(ETmerged_max_exclusions, sample != "SRI")
ETmerged_max_exclusions <- subset(ETmerged_max_exclusions, sample != "Vietnam Veterans")

#Set sub-zero variance components to zero (using mutate command)
ETmerged_max_exclusions_nozero <- ETmerged_max_exclusions %>% mutate(h2 = replace(h2, h2 <0, 0))
ETmerged_max_exclusions_nozero <- ETmerged_max_exclusions_nozero %>% mutate(c2 = replace(c2, c2 <0, 0))
ETmerged_max_exclusions_nozero <- ETmerged_max_exclusions_nozero %>% mutate(h2primesmooth = replace(h2primesmooth, h2primesmooth <0, 0))
ETmerged_max_exclusions_nozero <- ETmerged_max_exclusions_nozero %>% mutate(c2primesmooth = replace(c2primesmooth, c2primesmooth <0, 0))

#7.  CREATE ALT_NORWAY DATA AND RUN PART 2 AGAIN

Alt_Norway <-  ETmerged_max_exclusions_nozero

#After creating Alt_Norway data comment this command out and run part 2 a second time

# ---------------------------PART 3: REPRODUCE E&T RESULTS---------------------------- #

#8. REPRODUCE ENGZELL & TROPF TABLE 1 RESULTS (CONTROLLING FOR GENDER IN SMOOTHED RESULTS)

#Reproduce the left side of Table 1 in E&T (smoothed, standardised and gender-controlled results).
summary(lm.cluster(data = ETmerged_max_exclusions_nozero, formula = scale(h2) ~ scale(rIGsmooth) + gender, cluster= ETmerged_max_exclusions_nozero$cohortID)) #h2
summary(lm.cluster(data = ETmerged_max_exclusions_nozero, formula = scale(c2) ~ scale(rIGsmooth) + gender, cluster= ETmerged_max_exclusions_nozero$cohortID)) #c2
summary(lm.cluster(data = ETmerged_max_exclusions_nozero, formula = scale(e2) ~ scale(rIGsmooth) + gender, cluster= ETmerged_max_exclusions_nozero$cohortID)) #e2

#Manually derive p-values from t-values (for a two-sided t-test):
2*pt(abs(-3.61417454), df=14,lower.tail=FALSE) 
2*pt(abs(4.7448572), df=14,lower.tail=FALSE)   
2*pt(abs(-0.06274453), df=14,lower.tail=FALSE) 

#Reproduction results:
#h2= -0.512 se= 0.142 t= -3.61 p= 0.000   Man-p = 0.003  #MATCH (with manual p)
#c2=  0.574 se= 0.121 t=  4.74 p= 0.000   Man-p = 0.000  #MATCH
#e2= -0.011 se= 0.179 t= -0.06 p= 0.950   Man-p = 0.951  #MATCH (with manual p)

#Reproduce the right side of Table 1 in E&T (smoothed, destandardised and gender-controlled results).
summary(lm.cluster(data = ETmerged_max_exclusions_nozero, formula = scale(h2primesmooth) ~ scale(rIGsmooth) + gender, cluster= ETmerged_max_exclusions_nozero$cohortID)) #de_h2
summary(lm.cluster(data = ETmerged_max_exclusions_nozero, formula = scale(c2primesmooth) ~ scale(rIGsmooth) + gender, cluster= ETmerged_max_exclusions_nozero$cohortID)) #de_c2
summary(lm.cluster(data = ETmerged_max_exclusions_nozero, formula = scale(e2primesmooth) ~ scale(rIGsmooth) + gender, cluster= ETmerged_max_exclusions_nozero$cohortID)) #de_e2

#Manually derive p-values from t-values (for a two-sided t-test):
2*pt(abs(-0.1957366), df=14,lower.tail=FALSE) 
2*pt(abs(4.9388966), df=14,lower.tail=FALSE)  
2*pt(abs(1.652004), df=14,lower.tail=FALSE)   

#Destandardised Results:
#de_h2= -0.056 se= 0.284 t= -0.20 p=0.845  Man-p = 0.848  #MATCH (with manual p))
#de_c2=  0.726 se= 0.147 t=  4.94 p=0.000  Man-p = 0.000  #MATCH
#de_e2=  0.342 se= 0.207 t=  1.65 p=0.099  Man-p = 0.121  #MATCH (with manual p))

#9. ALTERNATIVE METHOD OF OBTAINING MAIN RESULTS USING COEFTEST

#Reproduce the left side of Table 1 in E&T (smoothed, standardised and gender-controlled results).

#Create named models
h2model <-lm(data=ETmerged_max_exclusions_nozero, formula = scale(h2) ~ scale(rIGsmooth) + gender)
c2model <-lm(data=ETmerged_max_exclusions_nozero, formula = scale(c2) ~ scale(rIGsmooth) + gender)
e2model <-lm(data=ETmerged_max_exclusions_nozero, formula = scale(e2) ~ scale(rIGsmooth) + gender)

#Apply clustered SEs (#HC1 is the decfault small sample adjustment used in Stata)
coeftest(h2model,vcovCL, cluster=ETmerged_max_exclusions_nozero$cohortID, type="HC1")
coeftest(c2model,vcovCL, cluster=ETmerged_max_exclusions_nozero$cohortID, type="HC1")
coeftest(e2model,vcovCL, cluster=ETmerged_max_exclusions_nozero$cohortID, type="HC1")

#Manually derive p-values from t-values (for a two-sided t-test):
2*pt(abs(-3.6142), df=14,lower.tail=FALSE)
2*pt(abs(4.7449), df=14,lower.tail=FALSE) 
2*pt(abs(-0.0627), df=14,lower.tail=FALSE)

#Results
#h2= -0.512 SE= 0.142 t= -3.61 p= 0.002 Man-p= 0.003 #MATCH (with man-p)
#c2=  0.574 SE= 0.121 t=  4.74 p= 0.000 Man-p= 0.000 #MATCH
#e2= -0.011 SE= 0.179 t=-0.06  p= 0.951 Man-p= 0.951 #MATCH   

#Reproduce the right side of Table 1 in E&T (smoothed, destandardised and gender-controlled results).
#Create named models
h2DSmodel <-lm(data=ETmerged_max_exclusions_nozero, formula = scale(h2primesmooth) ~ scale(rIGsmooth) + gender)
c2DSmodel <-lm(data=ETmerged_max_exclusions_nozero, formula = scale(c2primesmooth) ~ scale(rIGsmooth) + gender)
e2DSmodel <-lm(data=ETmerged_max_exclusions_nozero, formula = scale(e2primesmooth) ~ scale(rIGsmooth) + gender)

#Apply clustered SEs (#HC1 is the decfault small sample adjustment used in Stata)
coeftest(h2DSmodel,vcovCL, cluster=ETmerged_max_exclusions_nozero$cohortID, type="HC1")
coeftest(c2DSmodel,vcovCL, cluster=ETmerged_max_exclusions_nozero$cohortID, type="HC1")
coeftest(e2DSmodel,vcovCL, cluster=ETmerged_max_exclusions_nozero$cohortID, type="HC1")

#Manually derive p-values from t-values (for a two-sided t-test):
2*pt(abs(-0.1957), df=14,lower.tail=FALSE) 
2*pt(abs(4.9389), df=14,lower.tail=FALSE)  
2*pt(abs(1.6520), df=14,lower.tail=FALSE)  

#Destandardised Results:
#de_h2= -0.056 se= 0.284 t= -0.20  p=0.847  Man-p = 0.848  #MATCH (with manual p))
#de_c2=  0.726 se= 0.147 t=  4.94  p=0.000  Man-p = 0.000  #MATCH
#de_e2=  0.342 se= 0.207 t=  1.65  p=0.113  Man-p = 0.121  #MATCH (with manual p))

#10.  REPRODUCE ENGZELL & TROPF FIGURES 1B AND 1C

#H2 on rIG
#Base chart with  resized labels and points
h2_on_rIG <- ggplot(ETmerged_max_exclusions_nozero, aes(rIGsmooth, h2))+geom_point(aes(color=countryname, shape=gender), size=1)+geom_text(hjust = -0.2, vjust = 0.4, aes(label=cohortIDshort, color=countryname), size=2.7)+theme_bw()+ stat_smooth(method=lm, alpha=0.3) +labs(x = "parent-offspring correlation", y = "heritability for educational attainment") + xlim(0.21, 0.69) +ylim(0,0.8) + guides(color=FALSE) 

#Reformat gender legend
h2_on_rIG <- h2_on_rIG + theme(legend.justification=c(1,1), legend.position=c(0.98,0.98), 
                               panel.border = element_rect(colour = "black", fill=NA),
                               legend.box.background = element_rect(colour = "black", size = 0.75),
                               legend.title = element_blank(),
                               legend.spacing.y = unit(0, "mm"))
#Resize text
h2_on_rIG <- h2_on_rIG + theme(text = element_text(size = 8), legend.text = element_text(size = 8), axis.text.y = element_text(size = 7.5), axis.text.x = element_text(size = 7.5))  

h2_on_rIG
#Save output to PDF 
#ggsave("h2_on_rIG.pdf", width=12, height=10, units ="cm")

#C2 on rIG
#Base chart with  resized labels and points
c2_on_rIG <- ggplot(ETmerged_max_exclusions_nozero, aes(rIGsmooth, c2))+geom_point(aes(color=countryname, shape=gender), size=1)+geom_text(hjust = -0.2, vjust = 0.4, aes(label=cohortIDshort, color=countryname), size=2.7)+theme_bw()+ stat_smooth(method=lm, alpha=0.3) +labs(x = "parent-offspring correlation", y = "shared environmental influence") + xlim(0.21, 0.69) +ylim(0,0.8) + guides(color=FALSE)

c2_on_rIG <- c2_on_rIG + theme(legend.justification=c(0,1), legend.position=c(0.02,0.98), 
                               panel.border = element_rect(colour = "black", fill=NA),
                               legend.box.background = element_rect(colour = "black", size = 0.75),
                               legend.title = element_blank(),
                               legend.spacing.y = unit(0, "mm"))
#Resize text
c2_on_rIG <- c2_on_rIG + theme(text = element_text(size = 8), legend.text = element_text(size = 8), axis.text.y = element_text(size = 7.5), axis.text.x = element_text(size = 7.5))  

c2_on_rIG
#Save output to PDF
#ggsave("c2_on_rIG.pdf", width=12, height=10, units ="cm") 

#E2 on rIG
#Base chart with  resized labels and points
e2_on_rIG <- ggplot(ETmerged_max_exclusions_nozero, aes(rIGsmooth, e2))+geom_point(aes(color=countryname, shape=gender), size=1)+geom_text(hjust = -0.2, vjust = 0.4, aes(label=cohortIDshort, color=countryname), size=2.7)+theme_bw()+ stat_smooth(method=lm, alpha=0.3) +labs(x = "parent-offspring correlation", y = "non-shared environmental influence") + xlim(0.21, 0.69) +ylim(0,0.8) + guides(color=FALSE) 

#Reformat gender legend
e2_on_rIG <- e2_on_rIG + theme(legend.justification=c(1,1), legend.position=c(0.98,0.98), 
                               panel.border = element_rect(colour = "black", fill=NA),
                               legend.box.background = element_rect(colour = "black", size = 0.75),
                               legend.title = element_blank(),
                               legend.spacing.y = unit(0, "mm"))
#Resize text
e2_on_rIG <- e2_on_rIG + theme(text = element_text(size = 8), legend.text = element_text(size = 8), axis.text.y = element_text(size = 7.5), axis.text.x = element_text(size = 7.5))  

e2_on_rIG
#Save output to PDF
#ggsave("e2_on_rIG.pdf", width=12, height=10, units ="cm")

# ------------------------------PART 4: REVISED RESULTS------------------------------- #

#11. ALT_NORWAY RESULTS. REVISED TABLE 1 RESULTS

#New results for left side of Table 1 in E&T (smoothed, standardised and gender-controlled).
summary(lm.cluster(data = Alt_Norway, formula = scale(h2) ~ scale(rIGsmooth) + gender, cluster= Alt_Norway$cohortID)) #Nor_h2
summary(lm.cluster(data = Alt_Norway, formula = scale(c2) ~ scale(rIGsmooth) + gender, cluster= Alt_Norway$cohortID)) #Nor_c2
summary(lm.cluster(data = Alt_Norway, formula = scale(e2) ~ scale(rIGsmooth) + gender, cluster= Alt_Norway$cohortID)) #Nor_e2

#Manually derive p-values from t-values for Alt_Norway (for a two-sided t-test):
2*pt(abs(-1.0320512), df=14,lower.tail=FALSE)
2*pt(abs(2.48367345), df=14,lower.tail=FALSE)
2*pt(abs(-2.828329), df=14,lower.tail=FALSE) 

#Results with Alt_Norway values
#Nor_h2= -0.208 se= 0.202 t= -1.03 p= 0.302 Man-p = 0.320  #ATTENUATED. NO LONGER SIGNIFICANT.
#Nor_c2=  0.423 se= 0.170 t=  2.48 p= 0.013 Man-p = 0.026  #ATTENUATED. STILL SIGNIFICANT.
#Nor_e2= -0.433 se= 0.153 t= -2.83 p= 0.005 Man-p = 0.013  #STRENGTHENED. NOW SIGNIFICANT.

#New results for right side of Table 1 in E&T (smoothed, destandardised and gender-controlled).
summary(lm.cluster(data = Alt_Norway, formula = scale(h2primesmooth) ~ scale(rIGsmooth) + gender, cluster= Alt_Norway$cohortID)) #Nor_de_h2
summary(lm.cluster(data = Alt_Norway, formula = scale(c2primesmooth) ~ scale(rIGsmooth) + gender, cluster= Alt_Norway$cohortID)) #Nor_de_c2
summary(lm.cluster(data = Alt_Norway, formula = scale(e2primesmooth) ~ scale(rIGsmooth) + gender, cluster= Alt_Norway$cohortID)) #Nor_de_e2

#Manually derive p-values from t-values for Alt_Norway (for a two-sided t-test):
2*pt(abs(1.2811099), df=14,lower.tail=FALSE)
2*pt(abs(3.4124371), df=14,lower.tail=FALSE)
2*pt(abs(0.6433726), df=14,lower.tail=FALSE)

#Destandardised results with Alt_Norway values
#Nor_de_h2=  0.263 se= 0.205 t= 1.28 p=0.200 Man-p= 0.221  #SIGN REVERSED. BUT STILL NON-SIGNIFICANT
#Nor_de_c2=  0.660 se= 0.193 t= 3.41 p=0.001 Man-p= 0.004  #ATTENUATED. BUT STILL SIGNIFICANT
#Nor_de_e2=  0.151 se= 0.235 t= 0.64 p=0.520 Man-p= 0.530  #ATTENUATED. BUT STILL NON-SIGNIFICANT

#ALTERNATIVE METHOD OF OBTAINING TABLE 1 RESULTS USING COEFTEST

#New results for left side of Table 1 in E&T (smoothed, standardised and gender-controlled).
#Create named models
h2model_AltNor <-lm(data=Alt_Norway, formula = scale(h2) ~ scale(rIGsmooth) + gender)
c2model_AltNor <-lm(data=Alt_Norway, formula = scale(c2) ~ scale(rIGsmooth) + gender)
e2model_AltNor <-lm(data=Alt_Norway, formula = scale(e2) ~ scale(rIGsmooth) + gender)

#Apply clustered SEs (#HC1 is the decfault small sample adjustment used in Stata)
coeftest(h2model_AltNor,vcovCL, cluster=Alt_Norway$cohortID, type="HC1")
coeftest(c2model_AltNor,vcovCL, cluster=Alt_Norway$cohortID, type="HC1")
coeftest(e2model_AltNor,vcovCL, cluster=Alt_Norway$cohortID, type="HC1")

#Manually derive p-values from t-values for Alt_Norway (for a two-sided t-test):
2*pt(abs(-1.0321), df=14,lower.tail=FALSE) 
2*pt(abs(2.4837), df=14,lower.tail=FALSE)  
2*pt(abs(-2.8283), df=14,lower.tail=FALSE) 

#Alt_Norway results
#Nor_h2= -0.208 se= 0.202 t= -1.03 p= 0.313 Man-p = 0.320  #ATTENUATED. NO LONGER SIGNIFICANT.
#Nor_c2=  0.423 se= 0.170 t=  2.48 p= 0.021 Man-p = 0.026  #ATTENUATED. STILL SIGNIFICANT.
#Nor_e2= -0.433 se= 0.153 t= -2.83 p= 0.010 Man-p = 0.013  #STRENGTHENED. NOW SIGNIFICANT.

#New results for right side of Table 1 in E&T (smoothed, destandardised and gender-controlled).
#Create named models
h2DSmodel_AltNor <-lm(data=Alt_Norway, formula = scale(h2primesmooth) ~ scale(rIGsmooth) + gender)
c2DSmodel_AltNor <-lm(data=Alt_Norway, formula = scale(c2primesmooth) ~ scale(rIGsmooth) + gender)
e2DSmodel_AltNor <-lm(data=Alt_Norway, formula = scale(e2primesmooth) ~ scale(rIGsmooth) + gender)

#Apply clustered SEs (#HC1 is the decfault small sample adjustment used in Stata)
coeftest(h2DSmodel_AltNor,vcovCL, cluster=Alt_Norway$cohortID, type="HC1")
coeftest(c2DSmodel_AltNor,vcovCL, cluster=Alt_Norway$cohortID, type="HC1")
coeftest(e2DSmodel_AltNor,vcovCL, cluster=Alt_Norway$cohortID, type="HC1")

#Manually derive p-values from t-values (for a two-sided t-test):
2*pt(abs(1.2811), df=14,lower.tail=FALSE) 
2*pt(abs(3.4124), df=14,lower.tail=FALSE) 
2*pt(abs(0.6434), df=14,lower.tail=FALSE)  

#Destandardised results with Alt_Norway values
#Nor_de_h2=  0.263 se= 0.205 t= 1.28 p=0.214 Man-p= 0.221  #SIGN REVERSED. BUT STILL NON-SIGNIFICANT
#Nor_de_c2=  0.660 se= 0.193 t= 3.41 p=0.002 Man-p= 0.004  #ATTENUATED. BUT STILL SIGNIFICANT
#Nor_de_e2=  0.151 se= 0.235 t= 0.64 p=0.527 Man-p= 0.530  #ATTENUATED. BUT STILL NON-SIGNIFICANT

#12: ALT_NORWAY RESULTS: REVISED CHARTS

#Create reference ET models for inclusion in Alt_Norway versions of charts

#Create named linear models of ET model (this is different to the scaled mode with gender controls used on coeftest earlier)
ET_h2model<-lm(ETmerged_max_exclusions_nozero$h2~ETmerged_max_exclusions_nozero$rIGsmooth)
ET_c2model<-lm(ETmerged_max_exclusions_nozero$c2~ETmerged_max_exclusions_nozero$rIGsmooth)
ET_e2model<-lm(ETmerged_max_exclusions_nozero$e2~ETmerged_max_exclusions_nozero$rIGsmooth)

#Create named equations from these models
#Equations drawn from data
ET_h2equation <- function(x){coef(ET_h2model)[2]*x+coef(ET_h2model)[1]} #original h2 slope 
ET_c2equation <- function(x){coef(ET_c2model)[2]*x+coef(ET_c2model)[1]} #original C2 slope 
ET_e2equation <- function(x){coef(ET_e2model)[2]*x+coef(ET_e2model)[1]} #original e2 slope 

#Check ggplot colour scales (for 10 countries) to use for manual recolouring below.
show_col(hue_pal()(10))

#F8766D  Australia  (in within-country plot)
#D89000  Denmark
#A3A500  Finland    
#39B600  Germany
#00BF7D  Italy
#00BFC4  Norway     (in within-country plot)
#00B0F6  Spain
#9590FF  Sweden
#E76BF3  United Kingdom  (in within-country plot)
#FF62BC  United States   (in within-country plot)

#H2 on rIG (Alt_Norway)
#Base chart with resized labels and points
h2_on_rIG_Alt_Nor <- ggplot(Alt_Norway, aes(rIGsmooth, h2))+geom_point(aes(color=countryname, shape=gender), size=1)+geom_text(hjust = -0.2, vjust = 0.4, aes(label=cohortIDshort, color=countryname), size=2.7)+theme_bw()+ stat_smooth(method=lm, alpha=0.3) +labs(x = "parent-offspring correlation", y = "heritability for educational attainment") + xlim(0.21, 0.69) +ylim(0,0.8) + guides(color=FALSE) + stat_function(fun=ET_h2equation,geom="line", size=0.5, aes(colour = "original E&T slope"), linetype="solid")

#Change colours manually:
h2_on_rIG_Alt_Nor <-h2_on_rIG_Alt_Nor + scale_colour_manual(values = c("red", "#F8766D", "#D89000", "#A3A500", "#39B600", "#00BF7D", "#00BFC4", "#00B0F6", "#9590FF", "#E76BF3", "#FF62BC"), limits = c("original E&T slope", "Australia", "Denmark", "Finland", "Germany", "Italy", "Norway", "Spain", "Sweden", "United Kingdom", "United States"))

#Reformat gender legend
h2_on_rIG_Alt_Nor <- h2_on_rIG_Alt_Nor + theme(legend.justification=c(0,0), legend.position=c(0.02,0.02), 
                               panel.border = element_rect(colour = "black", fill=NA),
                               legend.box.background = element_rect(colour = "black", size = 0.75),
                               legend.title = element_blank(),
                               legend.spacing.y = unit(0, "mm"))
#Resize text
h2_on_rIG_Alt_Nor <- h2_on_rIG_Alt_Nor + theme(text = element_text(size = 8), legend.text = element_text(size = 8), axis.text.y = element_text(size = 7.5), axis.text.x = element_text(size = 7.5))  

h2_on_rIG_Alt_Nor
#Save output to PDF 
#ggsave("h2_on_rIG_Alt_Nor2_BL.pdf", width=12, height=10, units ="cm") #sets pdf output size (in cms)

#C2 on rIG (Alt_Norway)
#Base chart with  resized labels and points
c2_on_rIG_Alt_Nor <- ggplot(Alt_Norway, aes(rIGsmooth, c2))+geom_point(aes(color=countryname, shape=gender), size=1)+geom_text(hjust = -0.2, vjust = 0.4, aes(label=cohortIDshort, color=countryname), size=2.7)+theme_bw()+ stat_smooth(method=lm, alpha=0.3) +labs(x = "parent-offspring correlation", y = "shared environmental influence") + xlim(0.21, 0.69) +ylim(0,0.8) + guides(color=FALSE) + stat_function(fun=ET_c2equation,geom="line", size=0.5, aes(colour = "original E&T slope"), linetype="solid")

#Change colours manually:
c2_on_rIG_Alt_Nor <- c2_on_rIG_Alt_Nor + scale_colour_manual(values = c("red", "#F8766D", "#D89000", "#A3A500", "#39B600", "#00BF7D", "#00BFC4", "#00B0F6", "#9590FF", "#E76BF3", "#FF62BC"), limits = c("original E&T slope", "Australia", "Denmark", "Finland", "Germany", "Italy", "Norway", "Spain", "Sweden", "United Kingdom", "United States"))
#Reformat gender legend

c2_on_rIG_Alt_Nor <- c2_on_rIG_Alt_Nor + theme(legend.justification=c(0,1), legend.position=c(0.02,0.98), 
                               panel.border = element_rect(colour = "black", fill=NA),
                               legend.box.background = element_rect(colour = "black", size = 0.75),
                               legend.title = element_blank(),
                               legend.spacing.y = unit(0, "mm"))
#Resize text
c2_on_rIG_Alt_Nor <- c2_on_rIG_Alt_Nor + theme(text = element_text(size = 8), legend.text = element_text(size = 8), axis.text.y = element_text(size = 7.5), axis.text.x = element_text(size = 7.5))  

c2_on_rIG_Alt_Nor
#Save output to PDF
#ggsave("c2_on_rIG_Alt_Nor.pdf", width=12, height=10, units ="cm")

#E2 on rIG
#Base chart with  resized labels and points
e2_on_rIG_Alt_Nor <- ggplot(Alt_Norway, aes(rIGsmooth, e2))+geom_point(aes(color=countryname, shape=gender), size=1)+geom_text(hjust = -0.2, vjust = 0.4, aes(label=cohortIDshort, color=countryname), size=2.7)+theme_bw()+ stat_smooth(method=lm, alpha=0.3) +labs(x = "parent-offspring correlation", y = "non-shared environmental influence") + xlim(0.21, 0.69) +ylim(0,0.8) + guides(color=FALSE) + stat_function(fun=ET_e2equation,geom="line", size=0.5, aes(colour = "original E&T slope"), linetype="solid")

#Change colours manually:
e2_on_rIG_Alt_Nor <- e2_on_rIG_Alt_Nor + scale_colour_manual(values = c("red", "#F8766D", "#D89000", "#A3A500", "#39B600", "#00BF7D", "#00BFC4", "#00B0F6", "#9590FF", "#E76BF3", "#FF62BC"), limits = c("original E&T slope", "Australia", "Denmark", "Finland", "Germany", "Italy", "Norway", "Spain", "Sweden", "United Kingdom", "United States"))

#Reformat gender legend
e2_on_rIG_Alt_Nor <- e2_on_rIG_Alt_Nor + theme(legend.justification=c(1,1), legend.position=c(0.98,0.98), 
                               panel.border = element_rect(colour = "black", fill=NA),
                               legend.box.background = element_rect(colour = "black", size = 0.75),
                               legend.title = element_blank(),
                               legend.spacing.y = unit(0, "mm"))
#Resize text
e2_on_rIG_Alt_Nor <- e2_on_rIG_Alt_Nor + theme(text = element_text(size = 8), legend.text = element_text(size = 8), axis.text.y = element_text(size = 7.5), axis.text.x = element_text(size = 7.5))  

e2_on_rIG_Alt_Nor
#Save output to PDF 
#ggsave("e2_on_rIG_Alt_Nor.pdf", width=12, height=10, units ="cm") #sets pdf output size (in cms)

#13. PRECISION-WEIGHTED RESULTS. REVISED PLOTS

#Base chart with resized labels and points
h2_on_rIG_WGHT <- ggplot(ETmerged_max_exclusions_nozero, aes(rIGsmooth, h2))+geom_point(aes(color=countryname, shape=gender, size=h2weight)) + geom_text(hjust = -0.2, vjust = 0.4, aes(label=cohortIDshort, color=countryname, size=h2weight)) + theme_bw()+ stat_smooth(method=lm, formula=y~x, alpha=0.3, aes(weight=h2weight)) +labs(x = "parent-offspring correlation", y = "heritability for educational attainment") + xlim(0.21, 0.69) +ylim(0,0.8) + guides(color=FALSE, size=FALSE) + stat_function(fun=ET_h2equation,geom="line", size=0.5, aes(colour = "original E&T slope"), linetype="solid")

#Change colours manually:
h2_on_rIG_WGHT <- h2_on_rIG_WGHT + scale_colour_manual(values = c("red", "#F8766D", "#D89000", "#A3A500", "#39B600", "#00BF7D", "#00BFC4", "#00B0F6", "#9590FF", "#E76BF3", "#FF62BC"), limits = c("original E&T slope", "Australia", "Denmark", "Finland", "Germany", "Italy", "Norway", "Spain", "Sweden", "United Kingdom", "United States"))

#Reformat gender legend
h2_on_rIG_WGHT <- h2_on_rIG_WGHT + theme(legend.justification=c(1,1), legend.position=c(0.98,0.98), 
                               panel.border = element_rect(colour = "black", fill=NA),
                               legend.box.background = element_rect(colour = "black", size = 0.75),
                               legend.title = element_blank(),
                               legend.spacing.y = unit(0, "mm"))
#Text resizing
h2_on_rIG_WGHT <- h2_on_rIG_WGHT + theme(text = element_text(size = 8), legend.text = element_text(size = 8), axis.text.y = element_text(size = 7.5), axis.text.x = element_text(size = 7.5))  

h2_on_rIG_WGHT

#Save output to PDF 
#ggsave("h2_on_rIG_WGHT.pdf", width=12, height=10, units ="cm") #sets pdf output size (in cms)

#14: PRECISION WEIGHTED RESULTS. REVISED TABLE 1 DATA

#Apply lm.cluster model with precison weights, which controls for gender and uses scaled variables. 
summary(lm.cluster(data = ETmerged_max_exclusions_nozero, formula = scale(h2) ~ scale(rIGsmooth) + gender, cluster= ETmerged_max_exclusions_nozero$cohortID, weights = ETmerged_max_exclusions_nozero$h2weight)) #h2
summary(lm.cluster(data = ETmerged_max_exclusions_nozero, formula = scale(c2) ~ scale(rIGsmooth) + gender, cluster= ETmerged_max_exclusions_nozero$cohortID, ETmerged_max_exclusions_nozero$c2weight)) #c2
summary(lm.cluster(data = ETmerged_max_exclusions_nozero, formula = scale(e2) ~ scale(rIGsmooth) + gender, cluster= ETmerged_max_exclusions_nozero$cohortID, ETmerged_max_exclusions_nozero$e2weight)) #e2

#Manually check t-values (for two-sided t-test):
2*pt(abs(-1.1773333), df=14,lower.tail=FALSE) 
2*pt(abs(5.268264), df=14,lower.tail=FALSE)   
2*pt(abs(-0.346513), df=14,lower.tail=FALSE)  

#RESULTS
#           Estimate       SE        t value      Pr(>|t|)       Man p
#h2weights  -0.3007034   0.2554106   -1.1773333   0.2390625      0.2586906     #NO LONGER SIGNIFICANT
#c2weights   0.4881143   0.09265184   5.268264    1.377202e-07   0.0001188292  #STILL SIGNIFICANT
#e2weights  -0.1157921   0.3341638   -0.346513    7.289572e-01   0.7341101     #STILL NON-SIGNIFICANT

#Now Corroborate using lm() and with vcovCLfor clustered SEs
#First create named models
h2modelWEIGHTED <-lm(data=ETmerged_max_exclusions_nozero, formula = scale(h2) ~ scale(rIGsmooth) + gender, weights = h2weight)
c2modelWEIGHTED <-lm(data=ETmerged_max_exclusions_nozero, formula = scale(c2) ~ scale(rIGsmooth) + gender, weights = c2weight)
e2modelWEIGHTED <-lm(data=ETmerged_max_exclusions_nozero, formula = scale(e2) ~ scale(rIGsmooth) + gender, weights = e2weight)

#Second apply clustered SEs (#HC1 is the decfault small sample adjustment used in Stata)
coeftest(h2modelWEIGHTED,vcovCL, cluster=ETmerged_max_exclusions_nozero$cohortID, type="HC1")
coeftest(c2modelWEIGHTED,vcovCL, cluster=ETmerged_max_exclusions_nozero$cohortID, type="HC1")
coeftest(e2modelWEIGHTED,vcovCL, cluster=ETmerged_max_exclusions_nozero$cohortID, type="HC1")

#Manually check t-values (for two-sided t-test):
2*pt(abs(-1.1773), df=14,lower.tail=FALSE) 
2*pt(abs(5.2683), df=14,lower.tail=FALSE)  
2*pt(abs(-0.3465), df=14,lower.tail=FALSE) 

#RESULTS
#      Estimate       SE      t value    Pr(>|t|)       Man p
#h2=   -0.30070    0.25541   -1.1773     0.2516         0.2587035    #MATCH ON ESTIMATE, SE and T. DIFFERENT P
#c2=    0.488114   0.092652   5.2683     2.758e-05 ***  0.0001188214 #MATCH ON ESTIMATE, SE and T. DIFFERENT P
#e2=   -0.11579    0.33416   -0.3465     0.7322495      0.7341196    #MATCH ON ESTIMATE, SE and T. DIFFERENT P

#15. PLOT WITHIN-COUNTRY TRENDS FOR EACH GENDER

#Create country_gender variable to use as factor
ETmerged_max_exclusions_nozero$country_gender <- paste0(ETmerged_max_exclusions_nozero$countryname, ETmerged_max_exclusions_nozero$gender)

#Base chart with resized labels and points
h2WithinCountrychart <-ggplot(ETmerged_max_exclusions_nozero, aes(rIGsmooth, h2, color=factor(country_gender)))+geom_point(aes(shape=gender), size=1)+geom_text(hjust = -0.2, vjust = 0.4, aes(label=cohortIDshort), size=2.7)+theme_bw()+ stat_smooth(method=lm, alpha=0.0, size=0.5) + labs(x = "parent-offspring correlation", y= "heritability for educational attainment") +xlim(0.21, 0.69) + ylim(0,0.8) + scale_shape_manual(values=c(16,17,15)) + stat_function(fun=ET_h2equation,geom="line", size=0.5, aes(colour = "original E&T slope"), linetype="solid") 

#Change colours manually:
h2WithinCountrychart <-h2WithinCountrychart + scale_colour_manual(values = c("red", "#F8766D", "#F8766D", "#00BFC4","#00BFC4", "#E76BF3", "#FF62BC", "#FF62BC", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey"), limits = c("original E&T slope", "AustraliaMale", "AustraliaFemale", "NorwayMale", "NorwayFemale", "United KingdomMixed", "United StatesMale", "United StatesFemale", "DenmarkMale", "FinlandMale", "FinlandFemale", "GermanyMale", "GermanyFemale", "ItalyMale", "ItalyFemale", "SpainMale", "SpainFemale", "SwedenMixed"))+ guides(color=FALSE)

#Edit and reformat legend
h2WithinCountrychart <- h2WithinCountrychart + theme(legend.justification=c(1,1), legend.position=c(0.98,0.98), 
                                 panel.border = element_rect(colour = "black", fill=NA),
                                 legend.box.background = element_rect(colour = "black", size = 0.75),
                                 legend.title = element_blank(),
                                 legend.spacing.y = unit(0, "mm"))

#Resize text
h2WithinCountrychart <- h2WithinCountrychart + theme(text = element_text(size = 8), legend.text = element_text(size = 8), axis.text.y = element_text(size = 7.5), axis.text.x = element_text(size = 7.5))  

h2WithinCountrychart
#ggsave("WithinCountry.pdf", width=12, height=10, units ="cm") #sets pdf output size (in cms)

#16. CREATE MULTI-PANEL FIGURE FOR LETTER

#Create new named objects for each panel in the figure
Fig1A <- h2_on_rIG_Alt_Nor + guides(shape = FALSE) + ggtitle("A") +ylab("heritability for educational attainment") + xlab(" ") + xlim(0.21, 0.71)
Fig1B <- h2_on_rIG_WGHT + guides(shape = FALSE) + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) + ggtitle("B") + xlim(0.21, 0.71)
Fig1C <- h2WithinCountrychart + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) + ggtitle("C") + xlab(" ") + xlim(0.21, 0.71) + theme(legend.key.width=unit(0.4,"cm"),legend.key.height=unit(0.4,"cm"))

Grob1A <-ggplotGrob(Fig1A)
Grob1B <-ggplotGrob(Fig1B)
Grob1C <-ggplotGrob(Fig1C)

plot_grid(Grob1A, Grob1B, Grob1C, align = "h", ncol = 3, rel_widths = c(.37, .315, .315))
#ggsave("Figure_1.pdf", width=17.4, height=7.5, units ="cm") #sets pdf output size (in cms)

#17. PREPARE PRECISION-WEIGHTED CHARTS FOR SUPPLEMENTARY MATERIALS

#Precision-weighted c2 results
#Base chart with resized labels and points
c2_on_rIG_WGHT <- ggplot(ETmerged_max_exclusions_nozero, aes(rIGsmooth, c2))+geom_point(aes(color=countryname, shape=gender, size=c2weight)) + geom_text(hjust = -0.2, vjust = 0.4, aes(label=cohortIDshort, color=countryname, size=c2weight)) + theme_bw()+ stat_smooth(method=lm, formula=y~x, alpha=0.3, aes(weight=c2weight)) +labs(x = "parent-offspring correlation", y = "shared environmental influence") + xlim(0.21, 0.69) +ylim(0,0.8) + guides(color=FALSE, size=FALSE) + stat_function(fun=ET_c2equation,geom="line", size=0.5, aes(colour = "original E&T slope"), linetype="solid")

#Change colours manually:
c2_on_rIG_WGHT <- c2_on_rIG_WGHT + scale_colour_manual(values = c("red", "#F8766D", "#D89000", "#A3A500", "#39B600", "#00BF7D", "#00BFC4", "#00B0F6", "#9590FF", "#E76BF3", "#FF62BC"), limits = c("original E&T slope", "Australia", "Denmark", "Finland", "Germany", "Italy", "Norway", "Spain", "Sweden", "United Kingdom", "United States"))

#Reformat gender legend
c2_on_rIG_WGHT <- c2_on_rIG_WGHT + theme(legend.justification=c(0,1), legend.position=c(0.02,0.98), 
                                         panel.border = element_rect(colour = "black", fill=NA),
                                         legend.box.background = element_rect(colour = "black", size = 0.75),
                                         legend.title = element_blank(),
                                         legend.spacing.y = unit(0, "mm"))
#Text resizing
c2_on_rIG_WGHT <- c2_on_rIG_WGHT + theme(text = element_text(size = 8), legend.text = element_text(size = 8), axis.text.y = element_text(size = 7.5), axis.text.x = element_text(size = 7.5))  

c2_on_rIG_WGHT

#Precision-weighted e2 results
#Base chart with resized labels and points
e2_on_rIG_WGHT <- ggplot(ETmerged_max_exclusions_nozero, aes(rIGsmooth, e2))+geom_point(aes(color=countryname, shape=gender, size=e2weight)) + geom_text(hjust = -0.2, vjust = 0.4, aes(label=cohortIDshort, color=countryname, size=e2weight)) + theme_bw()+ stat_smooth(method=lm, formula=y~x, alpha=0.3, aes(weight=e2weight)) +labs(x = "parent-offspring correlation", y = "un-shared environmental influence") + xlim(0.21, 0.69) +ylim(0,0.8) + guides(color=FALSE, size=FALSE) + stat_function(fun=ET_e2equation,geom="line", size=0.5, aes(colour = "original E&T slope"), linetype="solid")

#Change colours manually:
e2_on_rIG_WGHT <- e2_on_rIG_WGHT + scale_colour_manual(values = c("red", "#F8766D", "#D89000", "#A3A500", "#39B600", "#00BF7D", "#00BFC4", "#00B0F6", "#9590FF", "#E76BF3", "#FF62BC"), limits = c("original E&T slope", "Australia", "Denmark", "Finland", "Germany", "Italy", "Norway", "Spain", "Sweden", "United Kingdom", "United States"))

#Reformat gender legend
e2_on_rIG_WGHT <- e2_on_rIG_WGHT + theme(legend.justification=c(1,1), legend.position=c(0.98,0.98), 
                                         panel.border = element_rect(colour = "black", fill=NA),
                                         legend.box.background = element_rect(colour = "black", size = 0.75),
                                         legend.title = element_blank(),
                                         legend.spacing.y = unit(0, "mm"))
#Text resizing
e2_on_rIG_WGHT <- e2_on_rIG_WGHT + theme(text = element_text(size = 8), legend.text = element_text(size = 8), axis.text.y = element_text(size = 7.5), axis.text.x = element_text(size = 7.5))  

e2_on_rIG_WGHT

#18. PRODUCE WITHIN-COUNTRY C2 AND E2 CHARTS FOR SUPPLEMENTARY MATERIALS

#C2 WITHIN-COUNTRY CHART
#Base chart with  resized labels and points
c2WithinCountrychart <-ggplot(ETmerged_max_exclusions_nozero, aes(rIGsmooth, c2, color=factor(country_gender)))+geom_point(aes(shape=gender), size=1)+geom_text(hjust = -0.2, vjust = 0.4, aes(label=cohortIDshort), size=2.7)+theme_bw()+ stat_smooth(method=lm, alpha=0.0, size=0.5) + labs(x = "parent-offspring correlation", y= "shared environmental influence") +xlim(0.21, 0.69) + ylim(0,0.8) + scale_shape_manual(values=c(16,17,15)) + stat_function(fun=ET_c2equation,geom="line", size=0.5, aes(colour = "original E&T slope"), linetype="solid") 

#Change colours manually:
c2WithinCountrychart <-c2WithinCountrychart + scale_colour_manual(values = c("red", "#F8766D", "#F8766D", "#00BFC4","#00BFC4", "#E76BF3", "#FF62BC", "#FF62BC", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey"), limits = c("original E&T slope", "AustraliaMale", "AustraliaFemale", "NorwayMale", "NorwayFemale", "United KingdomMixed", "United StatesMale", "United StatesFemale", "DenmarkMale", "FinlandMale", "FinlandFemale", "GermanyMale", "GermanyFemale", "ItalyMale", "ItalyFemale", "SpainMale", "SpainFemale", "SwedenMixed"))+ guides(color=FALSE)

#Edit and reformat legend
c2WithinCountrychart <- c2WithinCountrychart + theme(legend.justification=c(0,1), legend.position=c(0.02,0.98), 
                                                     panel.border = element_rect(colour = "black", fill=NA),
                                                     legend.box.background = element_rect(colour = "black", size = 0.75),
                                                     legend.title = element_blank(),
                                                     legend.spacing.y = unit(0, "mm"))

#Resize text
c2WithinCountrychart <- c2WithinCountrychart + theme(text = element_text(size = 8), legend.text = element_text(size = 8), axis.text.y = element_text(size = 7.5), axis.text.x = element_text(size = 7.5))  

c2WithinCountrychart

#E2 WITHIN-COUNTRY CHART

#Layer Within-Country chart
e2WithinCountrychart <-ggplot(ETmerged_max_exclusions_nozero, aes(rIGsmooth, e2, color=factor(country_gender)))+geom_point(aes(shape=gender), size=1)+geom_text(hjust = -0.2, vjust = 0.4, aes(label=cohortIDshort), size=2.7)+theme_bw()+ stat_smooth(method=lm, alpha=0.0, size=0.5) + labs(x = "parent-offspring correlation", y= "non-shared environmental influence") +xlim(0.21, 0.69) + ylim(0,0.8) + scale_shape_manual(values=c(16,17,15)) + stat_function(fun=ET_e2equation,geom="line", size=0.5, aes(colour = "original E&T slope"), linetype="solid") 

#Change colours manually:
e2WithinCountrychart <-e2WithinCountrychart + scale_colour_manual(values = c("red", "#F8766D", "#F8766D", "#00BFC4","#00BFC4", "#E76BF3", "#FF62BC", "#FF62BC", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey", "darkgrey"), limits = c("original E&T slope", "AustraliaMale", "AustraliaFemale", "NorwayMale", "NorwayFemale", "United KingdomMixed", "United StatesMale", "United StatesFemale", "DenmarkMale", "FinlandMale", "FinlandFemale", "GermanyMale", "GermanyFemale", "ItalyMale", "ItalyFemale", "SpainMale", "SpainFemale", "SwedenMixed"))+ guides(color=FALSE)

#Edit and reformat legend
e2WithinCountrychart <- e2WithinCountrychart + theme(legend.justification=c(1,1), legend.position=c(0.98,0.98), 
                                                     panel.border = element_rect(colour = "black", fill=NA),
                                                     legend.box.background = element_rect(colour = "black", size = 0.75),
                                                     legend.title = element_blank(),
                                                     legend.spacing.y = unit(0, "mm"))
#Resize text
e2WithinCountrychart <- e2WithinCountrychart + theme(text = element_text(size = 8), legend.text = element_text(size = 8), axis.text.y = element_text(size = 7.5), axis.text.x = element_text(size = 7.5))  

e2WithinCountrychart

#19. CREATE MULTI-PANEL FIGURES FOR SUPPLEMENTARY MATERIALS

#Fig S1: revised results for association between shared environmental influence and social mobility
FigS1A <- c2_on_rIG_Alt_Nor + guides(shape = FALSE) + ggtitle("A") +ylab("shared environmental influence") + xlab(" ") + xlim(0.21, 0.71)
FigS1B <- c2_on_rIG_WGHT + guides(shape = FALSE) + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) + ggtitle("B") + xlim(0.21, 0.71)
FigS1C <- c2WithinCountrychart + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) + ggtitle("C") + xlab(" ") + xlim(0.21, 0.71) + theme(legend.key.width=unit(0.4,"cm"),legend.key.height=unit(0.4,"cm"))

GrobS1A <-ggplotGrob(FigS1A)
GrobS1B <-ggplotGrob(FigS1B)
GrobS1C <-ggplotGrob(FigS1C)

plot_grid(GrobS1A, GrobS1B, GrobS1C, align = "h", ncol = 3, rel_widths = c(.37, .315, .315))
#ggsave("Figure_S1.pdf", width=17.4, height=7.5, units ="cm")

#Fig S2: revised results for association between non-shared environmental influence and social mobility
FigS2A <- e2_on_rIG_Alt_Nor + guides(shape = FALSE) + ggtitle("A") +ylab("non-shared environmental influence") + xlab(" ") + xlim(0.21, 0.71)
FigS2B <- e2_on_rIG_WGHT + guides(shape = FALSE) + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) + ggtitle("B") + xlim(0.21, 0.71)
FigS2C <- e2WithinCountrychart + theme(axis.title.y=element_blank(), axis.text.y=element_blank()) + ggtitle("C") + xlab(" ") + xlim(0.21, 0.71) + theme(legend.key.width=unit(0.4,"cm"),legend.key.height=unit(0.4,"cm"))

GrobS2A <-ggplotGrob(FigS2A)
GrobS2B <-ggplotGrob(FigS2B)
GrobS2C <-ggplotGrob(FigS2C)

plot_grid(GrobS2A, GrobS2B, GrobS2C, align = "h", ncol = 3, rel_widths = c(.37, .315, .315))
#ggsave("Figure_S2.pdf", width=17.4, height=7.5, units ="cm")

# ------------------------------------SCRIPT ENDS------------------------------------- #