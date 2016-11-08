###################################################################################################################
#################### MR analysis Type 2 Diabetes ##################################################################
###################################################################################################################

library(survival)
library(metafor)
library(ggplot2)
library(MendelianRandomization)
library(grid)
library(gridExtra) 
library(fmsb)

## This Function provides for association testing between the outcome (e.g. BMI) and all risk SNPs
## returns matrix with Betas (Row 1) and SEs (Row 2) for all SNPs
## Adjusted for Age, Sex, UKBB-genotyping platform, UKBB assessment centre, first 10 genetic PCs

association_function <- function(data, outcome, method = "linear" ){
  ## assign variable index for outcome
  ind <- which(names(data)==outcome)
  ## log regression for binary outcome 
  if (method == "binary") {
    l3 <- c(); l4 <- c(); o <- 1
    for(i in grep("^rs",names(data))){
      ## for(i in GRS_obesity)){  ## replace to get associations based on 53 obesity-unrelated SNPs  
      temp <- summary(glm(data[,ind] ~ data[,i] +  data$FUage + data$sex + data$gen_platform + data$ctr_1 + data$ctr_2 + data$ctr_3 + 
                            data$ctr_4 + data$ctr_5 + data$ctr_6 + data$ctr_7 + data$ctr_8 + data$ctr_9 + data$ctr_10 + data$ctr_11 +
                            data$ctr_12 + data$ctr_13 + data$ctr_14 + data$ctr_15 + data$ctr_16 + data$ctr_17 + data$ctr_18 + data$ctr_19 +
                            data$ctr_20 + data$ctr_21 + data$X22009.0.1 + data$X22009.0.2 + data$X22009.0.3 + data$X22009.0.4 + data$X22009.0.5 +
                            data$X22009.0.6 + data$X22009.0.7 + data$X22009.0.8 + data$X22009.0.9 + data$X22009.0.10, family="binomial"))
      l3[o] <- temp$coef[2,1]
      l4[o] <- temp$coef[2,2]
      o<-o+1
    }
  } 
  ## linear regression for continuous measures
  else if (method == "linear") {
    l3 <- c(); l4 <- c(); o <- 1
    for(i in grep("^rs",names(data))){
      ## for(i in GRS_obesity)){  ## replace to get associations based on 53 obesity-unrelated SNPs
      temp <- summary(glm(data[,ind] ~ data[,i] +  data$FUage + data$sex + data$gen_platform + data$ctr_1 + data$ctr_2 + data$ctr_3 +
                            data$ctr_4 + data$ctr_5 + data$ctr_6 + data$ctr_7 + data$ctr_8 + data$ctr_9 + data$ctr_10 + data$ctr_11 +
                            data$ctr_12 + data$ctr_13 + data$ctr_14 + data$ctr_15 + data$ctr_16 + data$ctr_17 + data$ctr_18 + data$ctr_19 +
                            data$ctr_20 + data$ctr_21 + data$X22009.0.1 + data$X22009.0.2 + data$X22009.0.3 + data$X22009.0.4 + data$X22009.0.5 +
                            data$X22009.0.6 + data$X22009.0.7 + data$X22009.0.8 + data$X22009.0.9 + data$X22009.0.10))
      l3[o] <- temp$coef[2,1] 
      l4[o] <- temp$coef[2,2]
                            o<-o+1
    }
  }
  ## Cox regression
  else if (method == "Cox"){
    l3 <- c(); l4 <- c(); o <- 1
    for(i in grep("^rs",names(data))){
      ## for(i in GRS_obesity)){  ## replace to get associations based on 53 obesity-unrelated SNPs
      temp <- summary(coxph(Surv(as.numeric(data$deathtime),data[,ind]) ~ data[,i] +  data$FUage + data$sex + data$gen_platform + data$ctr_1 + data$ctr_2 + data$ctr_3 +
                      data$ctr_4 + data$ctr_5 + data$ctr_6 + data$ctr_7 + data$ctr_8 + data$ctr_9 + data$ctr_10 + data$ctr_11 +
                      data$ctr_12 + data$ctr_13 + data$ctr_14 + data$ctr_15 + data$ctr_16 + data$ctr_17 + data$ctr_18 + data$ctr_19 +
                      data$ctr_20 + data$ctr_21 + data$X22009.0.1 + data$X22009.0.2 + data$X22009.0.3 + data$X22009.0.4 + data$X22009.0.5 +
                      data$X22009.0.6 + data$X22009.0.7 + data$X22009.0.8 + data$X22009.0.9 + data$X22009.0.10))
      l3[o] <- temp$coef[1,1]  # in survival model, first coeff is SNP (not as in GLM the intercept)
      l4[o] <- temp$coef[1,2]
      o<-o+1
    }
  }
  return(rbind(l3,l4))
} 

############################################################
#### UKBB upload and file processing (pheno/geno merge) ####
############################################################

#### Upload & align alleles ####

ukb_test <- read.table("/home/chrin/glob/private/ukb_w_phenotypes.txt") ## UKBB phenotypes
t2d <- read.table("/home/chrin/glob/private/ch_short_t2d.txt") ## UKBB extracted genotypes for all 120 SNPs reported by DIAGv3 and in the UKBB (out of 121)
t2d2 <- t(t2d)
snps <- read.delim("/home/chrin/glob/private/T2Dsnps_sign_associated.txt") ## 69 SNPs with genome-wide significant T2D associated and replicated in DIAGv3 whose OR 95% LCI does not incude < 1.000

snps_coded <- t2d2[max(nrow(t2d2)), 1:120] # as last row is UK_BB coded allele
t2d3 <- data.frame(t2d2[1:152249,1:120]) 
snps_names <- names(t2d3)
ids <- rownames(t2d3)
t <- data.frame(lapply(t2d3,function(x) as.numeric(as.character(x))))
rownames(t) <- ids
t2 <- t[,which(names(t) %in% as.character(snps$SNP_diagv3))] # select T2D-SNPs
snps <- (snps[order(as.character(snps$SNP_diagv3)),]) # make sure rs order is similar across files
t3 <- t2[,order(names(t2))] 

for (i in which(snps$risk_diagv3 != snps$risk_ukbb)){
  ## re-align if risk and ukbb_coded allele opposite
  t3[,i] <- -(t3[,i] - 2)
}

# write.table(t3,"t3.txt")

#### merge pheno and genotypes, genotypes exclusion criteria #### 

temp <- t3
temp$eid <- rownames(temp)
temp$eid <- gsub("^X","",temp$eid)
T2Dmaster <- merge(ukb_test, temp, by="eid")

T2Dmaster <- subset(T2Dmaster,X21000.0.0==1001) ## "white-British" included
T2Dmaster<-(T2Dmaster[-which(T2Dmaster$X22051.0.0==0),]) ## Affymetrix QC: exclude fails
T2Dmaster <- subset(T2Dmaster, T2Dmaster$X31.0.0 == T2Dmaster$X22001.0.0) ## genetic unequal reported sex excluded
T2Dmaster<-subset(T2Dmaster, X22006.0.0 == 1) ## genetic "Caucasians" included

names(T2Dmaster)[names(T2Dmaster)=="X21001.0.0"] <- "BMI"
names(T2Dmaster)[names(T2Dmaster)=="X31.0.0"] <- "sex"  
T2Dmaster$sex <- ordered(T2Dmaster$sex, levels=c(0,1), labels=c("female","male"))

T2Dmaster$gen_platform <- ifelse(T2Dmaster$X22000.0.0<0,0,1) ## Genotyping platform BliLEVE vs. Axiom dummy variable

centre_id <- unique(T2Dmaster$X54.0.0) ## assessment centre series of dummies called ctr_1, ctr_2 etc.
name <- paste("ctr",1:(length(centre_id)-1),sep="_") 
a <- list()
for (i in 1:(length(centre_id)-1)){
  a[[i]] <- ifelse(T2Dmaster$X54.0.0 == centre_id[i],1,0)
  T2Dmaster[name[i]] <- a[[i]]
}

#####################################################################
################ Analysis/MR ########################################
#####################################################################

#### produce table with N and % for outcomes ####

attach(T2Dmaster)
no <- c(table(died)[2],
        table(CVD_death_v1)[2], # Circ disease
        table(CVD_death_v2)[2], # Circ disease narrow
        table(CVD_death_v3)[2], # IHD
        table(CVD_death_v4)[2], # Cerebrovasc.
        table(CVD_v1)[2], # IHD
        table(CVD_v2)[2], # IHD_strict
        table(CVD_v3)[2], # Circulatory disease
        table(CVD_v4)[2], # Circulatory disease excluding HT
        table(CerebroVasc)[2],table(HaemStroke)[2],table(IschStroke)[2])
no_perc = 100*(no/nrow(T2Dmaster))
names <- c("Mortality", "Mortality_CirculatoryDis","Mortality_CirculatoryDis_narrow","Mortality_IHD","Mortality_Cerebrovasc",
           "IHD","IHD_strict","CirculatoryDisease","CirculatoryDis_wo_HT","CerebrovascDis","HaemStroke","IschStroke")
outcome_table <- data.frame(cbind(names,no,no_perc))

# write.table(outcome_table, "outcome_table.txt")

#### prepare sensitivity analysis excluding all obesity-associated SNPs (see separate script) #### 

obesity_snps <- c("rs11873305","rs12970134","rs1552224","rs1801282","rs7903146","rs9936385","rs7593730","rs13233731",
                  "rs13389219","rs6795735","rs7756992","rs780094","rs10203174","rs2943640","rs459193","rs8182584")
ob_pos <- which(names(T2Dmaster) %in% obesity_snps)
rs_pos <- grep("^rs",names(T2Dmaster))
GRS_obesity <- rs_pos[which(!(rs_pos %in% ob_pos))] ## has all variable-position for the 53 obesity-nonassociated SNPs in it
                                                    ## see adjustment in function above to run it for these SNPs only

#### to get the pseudoR2 from log regression for % explained variance in T2D by T2D-SNPs 

Nagel <- NagelkerkeR2(glm(T2Dmaster$T2D ~ (T2Dmaster[,GRS_which]) +  T2Dmaster$FUage + T2Dmaster$sex +  T2Dmaster$gen_platform+
                            T2Dmaster$ctr_1+T2Dmaster$ctr_2+T2Dmaster$ctr_3+T2Dmaster$ctr_4+T2Dmaster$ctr_5+T2Dmaster$ctr_6+T2Dmaster$ctr_7+
                            T2Dmaster$ctr_8+T2Dmaster$ctr_9+T2Dmaster$ctr_10+T2Dmaster$ctr_11+T2Dmaster$ctr_12+T2Dmaster$ctr_13+T2Dmaster$ctr_14+
                            T2Dmaster$ctr_15+T2Dmaster$ctr_16+T2Dmaster$ctr_17+T2Dmaster$ctr_18+T2Dmaster$ctr_19+T2Dmaster$ctr_20+T2Dmaster$ctr_21+
                            T2Dmaster$X22009.0.1+T2Dmaster$X22009.0.2+T2Dmaster$X22009.0.3+T2Dmaster$X22009.0.4+T2Dmaster$X22009.0.5+
                            T2Dmaster$X22009.0.6+T2Dmaster$X22009.0.7+T2Dmaster$X22009.0.8+T2Dmaster$X22009.0.9+T2Dmaster$X22009.0.10,family="binomial")) 

# write.table(Nagel,"T2D_nagelkerke.txt")

#### use custom function for get genetic-outcome associations

snp_names <- names(T2Dmaster)[grep("^rs",names(T2Dmaster))]

t2d_b <- association_function(T2Dmaster, outcome = "T2D", method= "binary")[1,]
t2d_se <- association_function(T2Dmaster, outcome = "T2D", method= "binary")[2,]
l3mort <- association_function(T2Dmaster, outcome = "died", method= "Cox")[1,]
l4mort <- association_function(T2Dmaster, outcome = "died", method= "Cox")[2,]
l3circ <- association_function(T2Dmaster, outcome = "CVD_v3", method= "binary")[1,]
l4circ <- association_function(T2Dmaster, outcome = "CVD_v3", method= "binary")[2,]
l3ihd <- association_function(T2Dmaster, outcome = "CVD_v1", method= "binary")[1,]
l4ihd <- association_function(T2Dmaster, outcome = "CVD_v1", method= "binary")[2,]
l3cere <- association_function(T2Dmaster, outcome = "CerebroVasc", method= "binary")[1,]
l4cere <- association_function(T2Dmaster, outcome = "CerebroVasc", method= "binary")[2,]
l3haem <- association_function(T2Dmaster, outcome = "HaemStroke", method= "binary")[1,]
l4haem <- association_function(T2Dmaster, outcome = "HaemStroke", method= "binary")[2,]
l3ischstroke <- association_function(T2Dmaster, outcome = "IschStroke", method= "binary")[1,]
l4ischstroke <- association_function(T2Dmaster, outcome = "IschStroke", method= "binary")[2,]

mr1 <- mr_allmethods(mr_input(bx = l1, bxse = l2, by = l3mort, byse = l4mort, exposure = "T2D_full_69snps", outcome = "Death", snps = snp_names))$Values
mr3 <- mr_allmethods(mr_input(bx = l1, bxse = l2, by = l3circ, byse = l4circ, exposure = "T2D_full_69snps", outcome = "Circ_dis", snps = snp_names))$Values
mr4 <- mr_allmethods(mr_input(bx = l1, bxse = l2, by = l3ihd, byse = l4ihd, exposure = "T2D_full_69snps", outcome = "IHD", snps = snp_names))$Values
mr5 <- mr_allmethods(mr_input(bx = l1, bxse = l2, by = l3cere, byse = l4cere, exposure = "T2D_full_69snps", outcome = "Cerebrovasc_Dis", snps = snp_names))$Values
mr6 <- mr_allmethods(mr_input(bx = l1, bxse = l2, by = l3ischstroke, byse = l4ischstroke, exposure = "T2D_full_69snps", outcome = "Isch_stroke", snps = snp_names))$Values
mr7 <- mr_allmethods(mr_input(bx = l1, bxse = l2, by = l3haem, byse = l4haem, exposure = "T2D_full_69snps", outcome = "Haem_stroke", snps = snp_names))$Values
 
list_MR <- list(mr1, mr3, mr4, mr5, mr6, mr7)
names(list_MR) <- c("Mortality", "Circ_Disease", "IHD", "Cerebrovasc_Disease", "Ischaemic_Stroke", "Haemorrhagic_Stroke")
 
# write.table(as.data.frame(do.call(rbind, list_MR)), "T2D_MR_results_full69snps_BurgessMethods.txt")

#### Scatter and Funnel plot ####
  
## below just an example for one outcome - see separate script with custom-made function that produces it automatically 

l3 <- l3mort
l4 <- l4mort
ivw_slope <- mr1[7,2]     # from penalized robust IVW
egger_slope <- mr1[14,2]  # from Penalised robust MR Egger
egger_intercept <- mr1[15,2]
exp <- "T2D_69snps"
out <- "Mort"
df <- data.frame(x = l1, y = l3,
                 ymin = l3 - (1.96*l4),ymax = l3 + (1.96*l4),
                 xmin = l1 - (1.96*(l2)),xmax = l1 + (1.96*(l2)))
 
p1 <- ggplot(data = df,aes(x = x,y = y)) + geom_point() + 
  geom_errorbar(aes(ymin = ymin,ymax = ymax)) + 
  geom_errorbarh(aes(xmin = xmin,xmax = xmax)) +
  geom_hline(aes(yintercept=0)) +
  geom_vline(aes(xintercept=0)) +
  ggtitle(paste(exp,out,sep = "_")) + 
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12)) +
  labs(x="Genetic association with risk factor",y="Genetic association with outcome") +
  geom_abline(slope = egger_slope, intercept = egger_intercept, col="red") +
  geom_abline(slope = ivw_slope, intercept = 0, col="blue")  

seiv <- (l4 / l1)
dfunnel <- data.frame(y = l1/l4,x = l3/l1,xmin = (l3/l1) - (1.96*seiv),xmax = (l3/l1) + (1.96*seiv))
p2 <- ggplot(data = dfunnel,aes(x = x,y = abs(y))) + 
  geom_point(size=2) + 
  geom_errorbarh(aes(xmin = xmin,xmax = xmax)) +
  geom_vline(xintercept=0) +
  geom_vline(xintercept=ivw_slope,col="blue") +
  geom_vline(xintercept=egger_slope,col="red") +
  ggtitle(paste(exp,out,sep = "_")) + 
  labs(y="IV strength",x="IV estimate")  +
  theme(axis.text=element_text(size=12),axis.title=element_text(size=12)) 

pdf(w=12, h=5.5,paste(exp,out,"pdf",sep = "."))
grid.arrange(p1, p2, ncol = 2)
dev.off()

################################################################################################## 
#### repeat for GRS excluding all obesity association variants... (the function code above) ######
##################################################################################################

#### Multivariable-MR ####

names(T2Dmaster)[which(names(T2Dmaster)=="X48.0.0")] <- "Waist"
names(T2Dmaster)[which(names(T2Dmaster)=="X49.0.0")] <- "Hip"
names(T2Dmaster)[which(names(T2Dmaster)=="X4080.0.0")] <- "SBP"
names(T2Dmaster)[which(names(T2Dmaster)=="X4079.0.0")] <- "DBP"
names(T2Dmaster)[which(names(T2Dmaster)=="X23099.0.0")] <- "Bodyfat_perc"

l1 <- association_function(T2Dmaster, outcome = "T2D", method= "binary")[1,]
l2 <- association_function(T2Dmaster, outcome = "T2D", method= "binary")[2,]
l5 <- association_function(T2Dmaster, outcome = "BMI", method= "linear")[1,]
l6 <- association_function(T2Dmaster, outcome = "BMI", method= "linear")[2,]
l7 <- association_function(T2Dmaster, outcome = "Waist", method= "linear")[1,]
l8 <- association_function(T2Dmaster, outcome = "Waist", method= "linear")[2,]
l9 <- association_function(T2Dmaster, outcome = "Hip", method= "linear")[1,]
l10 <- association_function(T2Dmaster, outcome = "Hip", method= "linear")[2,]
l11 <- association_function(T2Dmaster, outcome = "SBP", method= "linear")[1,]
l12 <- association_function(T2Dmaster, outcome = "SBP", method= "linear")[2,]
l13 <- association_function(T2Dmaster, outcome = "DBP", method= "linear")[1,]
l14 <- association_function(T2Dmaster, outcome = "DBP", method= "linear")[2,]
l15 <- association_function(T2Dmaster, outcome = "Bodyfat_perc", method= "linear")[1,]
l16 <- association_function(T2Dmaster, outcome = "Bodyfat_perc", method= "linear")[2,]

y_mort <- l3mort; y_se_mort <- l4mort
y_circdis <- l3circ; y_se_circdis <- l4circ
y_ihd <- l3ihd; y_se_ihd <- l4ihd
y_cerebro <- l3cere; y_se_cerebro <- l4cere
y_isch <- l3ischstroke; y_se_isch <- l4ischstroke
y_haem <- l3haem; y_se_haem <- l4haem

x1 <- l1 ; x1_se <- l2
x2 <-  l5  ; x3 <- l7 ; x4 <- l9  ; x5 <- l11 ; x6 <- l13 ; x7 <- l15 

#### MV-MR as per Burgess, Bowden et al. original paper and with robust SE published in addendum to it as a Letter
## adjusting for all, obesity-associated, and bp-associated variants
b_all_mort <- summary(lm(y_mort ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_mort^-2))$coef[1,1]
se_all_mort <- summary(lm(y_mort ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_mort^-2))$coef[1,2]/ summary(lm(y_mort ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_mort^-2))$sigma
b_ob_mort <- summary(lm(y_mort ~ x1+x2+x3+x4+x7-1, weights = y_se_mort^-2))$coef[1,1]
se_ob_mort <- summary(lm(y_mort ~ x1+x2+x3+x4+x7-1, weights = y_se_mort^-2))$coef[1,2]/ summary(lm(y_mort ~ x1+x2+x3+x4+x7-1, weights = y_se_mort^-2))$sigma
b_bp_mort <- summary(lm(y_mort ~ x1+x5+x6-1, weights = y_se_mort^-2))$coef[1,1]
se_bp_mort <- summary(lm(y_mort ~ x1+x5+x6-1, weights = y_se_mort^-2))$coef[1,2]/ summary(lm(y_mort ~ x1+x5+x6-1, weights = y_se_mort^-2))$sigma

b_all_circdis <- summary(lm(y_circdis ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_circdis^-2))$coef[1,1]
se_all_circdis <- summary(lm(y_circdis ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_circdis^-2))$coef[1,2]/ summary(lm(y_circdis ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_circdis^-2))$sigma
b_ob_circdis <- summary(lm(y_circdis ~ x1+x2+x3+x4+x7-1, weights = y_se_circdis^-2))$coef[1,1]
se_ob_circdis <- summary(lm(y_circdis ~ x1+x2+x3+x4+x7-1, weights = y_se_circdis^-2))$coef[1,2]/ summary(lm(y_circdis ~ x1+x2+x3+x4+x7-1, weights = y_se_circdis^-2))$sigma
b_bp_circdis <- summary(lm(y_circdis ~ x1+x5+x6-1, weights = y_se_circdis^-2))$coef[1,1]
se_bp_circdis <- summary(lm(y_circdis ~ x1+x5+x6-1, weights = y_se_circdis^-2))$coef[1,2]/ summary(lm(y_circdis ~ x1+x5+x6-1, weights = y_se_circdis^-2))$sigma

b_all_ihd <- summary(lm(y_ihd ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_ihd^-2))$coef[1,1]
se_all_ihd <- summary(lm(y_ihd ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_ihd^-2))$coef[1,2]/ summary(lm(y_ihd ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_ihd^-2))$sigma
b_ob_ihd <- summary(lm(y_ihd ~ x1+x2+x3+x4+x7-1, weights = y_se_ihd^-2))$coef[1,1]
se_ob_ihd <- summary(lm(y_ihd ~ x1+x2+x3+x4+x7-1, weights = y_se_ihd^-2))$coef[1,2]/ summary(lm(y_ihd ~ x1+x2+x3+x4+x7-1, weights = y_se_ihd^-2))$sigma
b_bp_ihd <- summary(lm(y_ihd ~ x1+x5+x6-1, weights = y_se_ihd^-2))$coef[1,1]
se_bp_ihd <- summary(lm(y_ihd ~ x1+x5+x6-1, weights = y_se_ihd^-2))$coef[1,2]/ summary(lm(y_ihd ~ x1+x5+x6-1, weights = y_se_ihd^-2))$sigma

b_all_cerebro <- summary(lm(y_cerebro ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_cerebro^-2))$coef[1,1]
se_all_cerebro <- summary(lm(y_cerebro ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_cerebro^-2))$coef[1,2]/ summary(lm(y_cerebro ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_cerebro^-2))$sigma
b_ob_cerebro <- summary(lm(y_cerebro ~ x1+x2+x3+x4+x7-1, weights = y_se_cerebro^-2))$coef[1,1]
se_ob_cerebro <- summary(lm(y_cerebro ~ x1+x2+x3+x4+x7-1, weights = y_se_cerebro^-2))$coef[1,2]/ summary(lm(y_cerebro ~ x1+x2+x3+x4+x7-1, weights = y_se_cerebro^-2))$sigma
b_bp_cerebro <- summary(lm(y_cerebro ~ x1+x5+x6-1, weights = y_se_cerebro^-2))$coef[1,1]
se_bp_cerebro <- summary(lm(y_cerebro ~ x1+x5+x6-1, weights = y_se_cerebro^-2))$coef[1,2]/ summary(lm(y_cerebro ~ x1+x5+x6-1, weights = y_se_cerebro^-2))$sigma

b_all_isch <- summary(lm(y_isch ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_isch^-2))$coef[1,1]
se_all_isch <- summary(lm(y_isch ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_isch^-2))$coef[1,2]/ summary(lm(y_isch ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_isch^-2))$sigma
b_ob_isch <- summary(lm(y_isch ~ x1+x2+x3+x4+x7-1, weights = y_se_isch^-2))$coef[1,1]
se_ob_isch <- summary(lm(y_isch ~ x1+x2+x3+x4+x7-1, weights = y_se_isch^-2))$coef[1,2]/ summary(lm(y_isch ~ x1+x2+x3+x4+x7-1, weights = y_se_isch^-2))$sigma
b_bp_isch <- summary(lm(y_isch ~ x1+x5+x6-1, weights = y_se_isch^-2))$coef[1,1]
se_bp_isch <- summary(lm(y_isch ~ x1+x5+x6-1, weights = y_se_isch^-2))$coef[1,2]/ summary(lm(y_isch ~ x1+x5+x6-1, weights = y_se_isch^-2))$sigma

b_all_haem <- summary(lm(y_haem ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_haem^-2))$coef[1,1]
se_all_haem <- summary(lm(y_haem ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_haem^-2))$coef[1,2]/ summary(lm(y_haem ~ x1+x2+x3+x4+x5+x6+x7-1, weights = y_se_haem^-2))$sigma
b_ob_haem <- summary(lm(y_haem ~ x1+x2+x3+x4+x7-1, weights = y_se_haem^-2))$coef[1,1]
se_ob_haem <- summary(lm(y_haem ~ x1+x2+x3+x4+x7-1, weights = y_se_haem^-2))$coef[1,2]/ summary(lm(y_haem ~ x1+x2+x3+x4+x7-1, weights = y_se_haem^-2))$sigma
b_bp_haem <- summary(lm(y_haem ~ x1+x5+x6-1, weights = y_se_haem^-2))$coef[1,1]
se_bp_haem <- summary(lm(y_haem ~ x1+x5+x6-1, weights = y_se_haem^-2))$coef[1,2]/ summary(lm(y_haem ~ x1+x5+x6-1, weights = y_se_haem^-2))$sigma


col <- c("outcome_69_SNPs","beta_all","se_all","beta_obesity","se_obesity","beta_bp","se_bp")
mort <- c("Mort",b_all_mort, se_all_mort, b_ob_mort, se_ob_mort, b_bp_mort, se_bp_mort)
cdis <- c("Circ_Dis",b_all_circdis, se_all_circdis, b_ob_circdis, se_ob_circdis, b_bp_circdis, se_bp_circdis)
ihd <- c("IHD",b_all_ihd, se_all_ihd, b_ob_ihd, se_ob_ihd, b_bp_ihd, se_bp_ihd)
cer <- c("Cerebrov",b_all_cerebro, se_all_cerebro, b_ob_cerebro, se_ob_cerebro, b_bp_cerebro, se_bp_cerebro)
isch <- c("Mort",b_all_isch, se_all_isch, b_ob_isch, se_ob_isch, b_bp_isch, se_bp_isch)
haem <- c("Mort",b_all_haem, se_all_haem, b_ob_haem, se_ob_haem, b_bp_haem, se_bp_haem)

r <- rbind(mort,cdis,ihd,cer,isch,haem)
colnames(r) <- col

# write.table(r, "MV_MR_results_69_SNPs.txt", quote=F,row.names=F)

######################################################################################################################################################################
 


