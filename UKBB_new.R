
##    dir<-"~/Desktop/UKBB_project/"
##    setwd(dir)
#####   ukb_test0 <- ukb4908[which(ukb4908$eid %in% ukb4908_add$eid),]
#####   ukb_test <- merge(ukb_test0,ukb4908_add,by="eid")
#####   write.table(ukb_test,"ukb_test.txt")
##    ukb4908 <- read.table("ukb_test.txt")
##    ukb4908 <- read.delim("app1372_standard_data_2016Oct15.txt")
##    ukb4908_add <- read.table("ukbb_5000.txt")
##    ukb6142 <- read.csv("~/Desktop/UKBB_project/ukb6142.csv") # chromosome intensities
##    ukb6030 <- read.csv("~/Desktop/UKBB_project/ukb6030.csv") 
##    ukb4908 <- read.csv("~/Desktop/UKBB_project/ukb4908.csv") 
##    EXCL <- read.csv("w1372_20160406_EXCUDE_THESE_INDIVIDUALS.csv")

ukb4908_add <- read.delim("/home/chrin/glob/private/app1372_standard_data_2016Oct15.txt")   ## new FU file 20 Oct 2016
ukb4908 <- read.csv("/home/chrin/glob/private/ukb4908.csv")                                 ## all other phenotypes
EXCL <- read.csv("/home/chrin/glob/private/w1372_20160406_EXCUDE_THESE_INDIVIDUALS.csv")    ## exclude non-consenters

ukb_test <- merge(ukb4908, ukb4908_add, by="eid")
ukb_test <- ukb_test[which( ukb_test$eid %in% EXCL !=T),] 

  ## Assessment centre nation, needed as different FU censoring dates
Wales <- c(11003,11022,11023)
Scotland <-c(11005,11004)  
ukb_test$NATION <- 1
ukb_test$NATION[which(ukb_test$X54.0.0 %in% Scotland)] <- 2
ukb_test$NATION[which(ukb_test$X54.0.0 %in% Wales)] <- 3
ukb_test$NATION  <- factor(ukb_test$NATION, labels = c("Engld", "Scotld", "Wales"))

  ## variable denotion X54.0.0 = original phenotype files
  ## n_40000_0_0 or s_40001_0_0 or ts_40007_0_0 etc. is updated FU file

##########################################################################################
############### PHENOTYPE DEFINITIONS ####################################################
##########################################################################################

##########################################################################################
############### TYPE 2 DIABETES MELLITUS #################################################
##########################################################################################
## self-reported "DM" or in-hospital diagnosis ICD-10 (E11) or ICD-9 (250 = any DM)
## on anti-diabetics medication incl. insulin
## Diabetes diagnosed by doctor
## age at diagnosis >30yrs - because in self-reported (rather: nurse interview-recorded)
## diseases, only a few are "T2D", many more are "DM"; only including "T2D" would exclude
## many false negatives with "DM" codes - to minimise confouding by T1D, all under-30-yo
## at DM-diagnosis are excluded 

base_noncancerselfreport <- grep("^n_20002",names(ukb_test))  
base_treatmentcode <- grep("^n_20003",names(ukb_test))

  ## 1) self-reported DM
temp_1 <- c()   # self-reported baseline "T2D"
for (i in base_noncancerselfreport){
  temp_1 <- append(temp_1,(which(ukb_test[,i]=="1220"))) 
  ## 1220 = "Diab Mell", 1223 = "Type 2 Diab" but only n=883
}

  ## 2) anti-diabetic medication
antidiab <- c(1140883066,  1140884600, 1141157284,  1141157284,  1141152590,  1140910566,  1140857138,  1141153254, 
              1141171646,  1141177600,  1141189090,  1141168660,  1141173882,  1140874674,  1140874664,  1140874706) 
    ## 1140883066 insulin product  # 1140884600 metformin  # 1141157284 glipizide  # 1141157284 glipizide product
    ## 1141152590 glimepiride  # 1140910566 glyclizide # 1140857138 sulfaurea  # 1141153254 troglitazone # 1141171646 pioglitazone 
    ## 1141177600 rosiglitazone  # 1141189090 rosiglit/metf  # 1141168660 repaglinide  # 1141173882 nateglinide
    ## 1140874674 tolbutamide  # 1140874664 tolazamide # 1140874706 chlorpropamide
temp_2<-c() 
for (i in base_treatmentcode){
  temp_2 <- append(temp_2,(which(ukb_test[,i] %in% antidiab))) 
}

  ## 3) ICD hospital-record DM
diagmainicd10 <- grep("^s_41202",names(ukb_test)) # E11-14 = DM
diagmainicd9 <- grep("^s_41203",names(ukb_test)) # 250 = "diabetes mellitus" 

temp_3 <- c()   # in-hospital Dx of ICD10-T2D
for (i in diagmainicd10){
  temp_3 <- append(temp_3,(grep("^E11",ukb_test[,i])))
}
temp_4 <- c() # in-hospital Dx of ICD9-"DM"
for (i in min(diagmainicd9):max(diagmainicd9)){
  temp_4 <- append(temp_4,(grep("^2",ukb_test[,i])))
}

  ## 4) physician-diagnosed 
temp_5 <- c(which(ukb_test$n_2443_0_0 == 1),which(ukb_test$n_2443_1_0 == 1)) 

a <- (unique(c(temp_1,temp_2,temp_3,temp_4,temp_5))) 
b <- (which(ukb_test$X2976.0.0 >30 | ukb_test$X2976.1.0 >30)) # DM-diagnosis age
t2d_id <- unique(c(a,b))
ukb_test$T2D <- 0
ukb_test[t2d_id,]$T2D <- 1

##########################################################################################
############### All-cause mortality  #####################################################
##########################################################################################

date_death <- grep("^ts_40000",names(ukb_test)) 

dead<-c()
for (i in date_death){
  dead <- append(dead,which(!is.na(as.Date(ukb_test[,i],format="%d%b%Y")))) 
}
ukb_test$died <- 0
ukb_test$died[dead] <- 1  ## dummy for anyone who died during FU

for (i in date_death){
  ukb_test[,i] <- as.Date(ukb_test[,i], format="%d%b%Y")
}
    ## as death date in 3 variables for 3 instances - merge into 1
ukb_test$m = ukb_test$ts_40000_0_0
ukb_test$m[!is.na(ukb_test$ts_40000_1_0)] = ukb_test$ts_40000_1_0[!is.na(ukb_test$ts_40000_1_0)] 
ukb_test$m[!is.na(ukb_test$ts_40000_2_0)] = ukb_test$ts_40000_2_0[!is.na(ukb_test$ts_40000_2_0)] 

  ## FU times for death = death date or censored at last complete FU according to UKBB announcement - differ for E&W and Scotland
ukb_test$deathtime <- ukb_test$m - as.Date(ukb_test$ts_53_0_0, format="%d%b%Y")
  ## E&W 31/1/2016 Scotland 30/11/2015
ukb_test$deathtime[is.na(ukb_test$deathtime) & !ukb_test$NATION %in% c("Engld","Wales")] <- as.Date("2015-11-30",format="%Y-%m-%d") - as.Date(ukb_test$ts_53_0_0[is.na(ukb_test$deathtime)  & !ukb_test$NATION %in% c("Engld","Wales")],format="%d%b%Y")
ukb_test$deathtime[is.na(ukb_test$deathtime) & ukb_test$NATION %in% c("Engld","Wales")] <- as.Date("2016-01-31",format="%Y-%m-%d") - as.Date(ukb_test$ts_53_0_0[is.na(ukb_test$deathtime)  & ukb_test$NATION %in% c("Engld","Wales")],format="%d%b%Y")
ukb_test$deathtime <- ukb_test$deathtime/365.25

ukb_test$FUage <- ukb_test$X21022.0.0 +  ukb_test$deathtime

ukb_test<- ukb_test[ukb_test$deathtime>0,] ## 2 have negative death dates... 

##########################################################################################
############################ CVD mortality (1y/main CoD unless stated otherwise) #########
##########################################################################################

  ## CVD_death_v1: due to any circulatory disease - ICD10 "I" category
ukb_test$CVD_death_v1 <- 0
ukb_test[c(grep("^I",ukb_test$s_40001_0_0),grep("^I",ukb_test$s_40001_1_0),grep("^I",ukb_test$s_40001_2_0)),]$CVD_death_v1 <- 1

  ## CVD_death_v2: death due to circulatory disease excluding I26 (PE) and I89 (other
  ## noninfectious lymphatic..); based on PMID 27050205 - 2016 NEJM paper
exclude <- c(grep("^I26",ukb_test$s_40001_0_0), grep("^I89",ukb_test$s_40001_0_0),
             grep("^I26",ukb_test$s_40001_1_0), grep("^I89",ukb_test$s_40001_1_0),
             grep("^I26",ukb_test$s_40001_2_0), grep("^I89",ukb_test$s_40001_2_0))
ukb_test$CVD_death_v2 <- ukb_test$CVD_death_v1
ukb_test[exclude,]$CVD_death_v2 <- 0

  ## CVD_death_v3: death due to IHD = ICD10 I20 to I25
chddeath <- c(grep("^I20",ukb_test$s_40001_0_0),grep("^I20",ukb_test$s_40001_1_0),grep("^I20",ukb_test$s_40001_2_0),
              grep("^I21",ukb_test$s_40001_0_0),grep("^I21",ukb_test$s_40001_1_0),grep("^I21",ukb_test$s_40001_2_0),
              grep("^I22",ukb_test$s_40001_0_0),grep("^I22",ukb_test$s_40001_1_0),grep("^I22",ukb_test$s_40001_2_0),
              grep("^I23",ukb_test$s_40001_0_0),grep("^I23",ukb_test$s_40001_1_0),grep("^I23",ukb_test$s_40001_2_0),
              grep("^I24",ukb_test$s_40001_0_0),grep("^I24",ukb_test$s_40001_1_0),grep("^I24",ukb_test$s_40001_2_0),
              grep("^I24",ukb_test$s_40001_0_0),grep("^I24",ukb_test$s_40001_1_0),grep("^I24",ukb_test$s_40001_2_0) )
ukb_test$CVD_death_v3 <- 0
ukb_test[chddeath,]$CVD_death_v3 <- 1

  ## CVD_death_v4: death due to cerebrovascular disease = ICD10 I60-I69
cerdeath <- c(grep("^I6",ukb_test$s_40001_0_0),grep("^I6",ukb_test$s_40001_1_0),grep("^I6",ukb_test$s_40001_2_0))
ukb_test$CVD_death_v4 <- 0
ukb_test[cerdeath,]$CVD_death_v4 <- 1

###########################################################################################
############################ CVD diagnosis  ###############################################
###########################################################################################
## no FU/censoring time required because will run log_regression for lifetime risk
## or CVD rather than incidental CVD during FU only (too few events for powerful MR)

loopId_icd10 <- grep("^s_41202",names(ukb_test))
loopId_icd9 <- grep("^s_41203",names(ukb_test))

  ## CVD_v1: IHD = I20-I25; 410-414
loopId_icd10 <- grep("^s_41202",names(ukb_test))
loopId_icd9 <- grep("^s_41203",names(ukb_test))

id_icd10_cvd_v1 <- c()
for (i in loopId_icd10){
  id_icd10_cvd_v1 <- append(id_icd10_cvd_v1, c(grep("^I20",ukb_test[,i]),grep("^I21",ukb_test[,i]),grep("^I22",ukb_test[,i]),
                                               grep("^I23",ukb_test[,i]),grep("^I24",ukb_test[,i]),grep("^I25",ukb_test[,i])))
}
id_icd10_cvd_v1 <- unique(id_icd10_cvd_v1)

id_icd9_cvd_v1 <- c()
for (i in loopId_icd9){
  id_icd9_cvd_v1 <- append(id_icd9_cvd_v1, c(grep("^410",ukb_test[,i]),grep("^411",ukb_test[,i]),grep("^412",ukb_test[,i]),
                                             grep("^413",ukb_test[,i]),grep("^414",ukb_test[,i])))
}
id_icd9_cvd_v1 <- unique(id_icd9_cvd_v1)

ukb_test$CVD_v1 <- 0
ukb_test[unique(c(id_icd10_cvd_v1,id_icd9_cvd_v1)),]$CVD_v1 <- 1

  ##  CVD_v2: IHD as defined as per EI's grant application
  ## I20.0, I21, I22; 410, 411; surgical CABG/PTCA
loopId_surg <- grep("^s_41200",names(ukb_test))

id_surg <- c()
for (i in loopId_surg){
  id_surg <- append(id_surg, c(grep("^K40",ukb_test[,i]),grep("^K41",ukb_test[,i]),grep("^K42",ukb_test[,i]),
                               grep("^K43",ukb_test[,i]),grep("^K44",ukb_test[,i]),grep("^K49",ukb_test[,i]),
                               grep("^K75",ukb_test[,i])))
}
id_surg <- unique(id_surg)

id_icd10_cvd_v2 <- c()
for (i in loopId_icd10){
  id_icd10_cvd_v2 <- append(id_icd10_cvd_v2, c(grep("^I200",ukb_test[,i]),grep("^I21",ukb_test[,i]),grep("^I22",ukb_test[,i])))
}
id_icd10_cvd_v2 <- unique(id_icd10_cvd_v2)

id_icd9_cvd_v2 <- c()
for (i in loopId_icd9){
  id_icd9_cvd_v2 <- append(id_icd9_cvd_v2, c(grep("^410",ukb_test[,i]),grep("^411",ukb_test[,i])))
}
id_icd9_cvd_v2 <- unique(id_icd9_cvd_v2)

ukb_test$CVD_v2 <- 0
ukb_test[unique(c(id_icd10_cvd_v2,id_icd9_cvd_v2,id_surg)),]$CVD_v2 <- 1

  ## CVD_v3: any circulatory dis
id_icd10_cvd_v3 <- c()
for (i in loopId_icd10){
  id_icd10_cvd_v3 <- append(id_icd10_cvd_v3, c(grep("^I",ukb_test[,i])))
}
id_icd10_cvd_v3 <- unique(id_icd10_cvd_v3)

id_icd9_cvd_v3 <- c()
for (i in loopId_icd9){
  id_icd9_cvd_v3 <- append(id_icd9_cvd_v3, c(grep("^39",ukb_test[,i]),grep("^40",ukb_test[,i]),grep("^41",ukb_test[,i]),grep("^42",ukb_test[,i]),
                                             grep("^43",ukb_test[,i]),grep("^44",ukb_test[,i]),grep("^45",ukb_test[,i])))
}
id_icd9_cvd_v3 <- unique(id_icd9_cvd_v3)

ukb_test$CVD_v3 <- 0
ukb_test[unique(c(id_icd10_cvd_v3,id_icd9_cvd_v3)),]$CVD_v3 <- 1

 ## CVD_v4: any circulatory dis excluding Hypertension
id_icd10_cvd_v4 <- c()
for (i in loopId_icd10){
  id_icd10_cvd_v4 <- append(id_icd10_cvd_v4, c(grep("^I0",ukb_test[,i]),grep("^I2",ukb_test[,i]),grep("^I3",ukb_test[,i]),grep("^I4",ukb_test[,i]),
                                               grep("^I5",ukb_test[,i]),grep("^I6",ukb_test[,i]),grep("^I7",ukb_test[,i]),grep("^I8",ukb_test[,i]),
                                               grep("^I9",ukb_test[,i])))
}
id_icd10_cvd_v4 <- unique(id_icd10_cvd_v4)

id_icd9_cvd_v4 <- c()
for (i in loopId_icd9){
  id_icd9_cvd_v4 <- append(id_icd9_cvd_v4, c(grep("^39",ukb_test[,i]),grep("^41",ukb_test[,i]),grep("^42",ukb_test[,i]),
                                             grep("^43",ukb_test[,i]),grep("^44",ukb_test[,i]),grep("^45",ukb_test[,i])))
}
id_icd9_cvd_v4 <- unique(id_icd9_cvd_v4)

ukb_test$CVD_v4 <- 0
ukb_test[unique(c(id_icd10_cvd_v4,id_icd9_cvd_v4)),]$CVD_v4 <- 1

 ## CerebroVasc I60-I69; 430-438
id_icd10_cerebrovasc <- c()
for (i in loopId_icd10){
  id_icd10_cerebrovasc <- append(id_icd10_cerebrovasc, c(grep("^I6",ukb_test[,i])))
}
id_icd10_cerebrovasc <- unique(id_icd10_cerebrovasc)

id_icd9_cerebrovasc <- c()
for (i in loopId_icd9){
  id_icd9_cerebrovasc <- append(id_icd9_cerebrovasc, c(grep("^43",ukb_test[,i])))
}
id_icd9_cerebrovasc <- unique(id_icd9_cerebrovasc)

ukb_test$CerebroVasc <- 0
ukb_test[unique(c(id_icd10_cerebrovasc,id_icd9_cerebrovasc)),]$CerebroVasc <- 1

  ## HaemStroke I60-I62; 430-432
id_icd10_haemstroke <- c()
for (i in loopId_icd10){
  id_icd10_haemstroke <- append(id_icd10_haemstroke, c(grep("^I60",ukb_test[,i]),grep("^I61",ukb_test[,i]),grep("^I62",ukb_test[,i])))
}
id_icd10_haemstroke <- unique(id_icd10_haemstroke)

id_icd9_haemstroke <- c()
for (i in loopId_icd9){
  id_icd9_haemstroke <- append(id_icd9_haemstroke, c(grep("^430",ukb_test[,i]),grep("^431",ukb_test[,i]),grep("^432",ukb_test[,i])))
}
id_icd9_haemstroke <- unique(id_icd9_haemstroke)

ukb_test$HaemStroke <- 0
ukb_test[unique(c(id_icd10_haemstroke,id_icd9_haemstroke)),]$HaemStroke <- 1

  ## IschStroke I63-I65; 433-435

id_icd10_ischstroke <- c()
for (i in loopId_icd10){
  id_icd10_ischstroke <- append(id_icd10_ischstroke, c(grep("^I63",ukb_test[,i]),grep("^I64",ukb_test[,i]),grep("^I65",ukb_test[,i])))
}
id_icd10_ischstroke <- unique(id_icd10_ischstroke)

id_icd9_ischstroke <- c()
for (i in loopId_icd9){
  id_icd9_ischstroke <- append(id_icd9_ischstroke, c(grep("^433",ukb_test[,i]),grep("^434",ukb_test[,i]),grep("^435",ukb_test[,i])))
}
id_icd9_ischstroke <- unique(id_icd9_ischstroke)

ukb_test$IschStroke <- 0
ukb_test[unique(c(id_icd10_ischstroke,id_icd9_ischstroke)),]$IschStroke <- 1

###
write.table(ukb_test,"ukb_w_phenotypes.txt")
###