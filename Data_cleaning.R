#################
#################
#Code: Nicolo Foppa Pedretti
#Edited: Elena Colicino
#Date: 20201012
#################
#################

#################
#library
library(nhanesA)
library(data.table)
library(rstan)
library(rvest)
library(monomvn)
library(survey)
library(corrplot)
library(BWQS)
library(textreadr)

#setwd("/Results/")

#################
# Download demographic data with no NA
demo_H = data.table(nhanes('DEMO_H'),YEAR='2013-2014')
demo_H = demo_H[WTMEC2YR>0,.(SEQN,RIAGENDR,RIDAGEYR,WTMEC2YR,INDFMPIR,RIDRETH1,YEAR)]
demo_H[,SUM_MEC:=sum(demo_H$WTMEC2YR)]
demo_H[, check := rowSums(.SD), .SDcol=c('RIAGENDR','RIDAGEYR','INDFMPIR','RIDRETH1')]
demo_H = demo_H[!is.na(check)]
demo_H[,check:=NULL]

#################
# Download PFCs data with no NA
pfc_H = data.table(nhanes("PFAS_H"),YEAR='2013-2014')
pfc_list <- c("LBXPFHS", "LBXMPAH", "LBXPFDE",  "LBXPFBS", "LBXPFHP", "LBXPFNA", "LBXPFUA", "LBXPFDO")
pfc_H = pfc_H[,c("SEQN","WTSB2YR",pfc_list,"YEAR"),with=F]
pfc_H[, check := rowSums(.SD), .SDcol=pfc_list]
pfc_H = pfc_H[!is.na(check)]
pfc_H[,check:=NULL]


pfc_H_v2 = data.table(nhanes("SSPFAS_H"),YEAR='2013-2014')
pfc_list_v2 <- c("SSNPFOA", "SSBPFOA", "SSNPFOS", "SSMPFOS")
pfc_H_v2 = pfc_H_v2[,c("SEQN",pfc_list_v2),with=F]
pfc_H_v2[, check := rowSums(.SD), .SDcol=pfc_list_v2]
pfc_H_v2 = pfc_H_v2[!is.na(check)]
pfc_H_v2[,check:=NULL]

pfc_H = Reduce(function(...) merge(...,by=c("SEQN")), list(pfc_H,pfc_H_v2))

#################
# Download bones mineral density data with no NA
fem_H = data.table(nhanes("DXXFEM_H"),YEAR='2013-2014')
spn_H = data.table(nhanes("DXXSPN_H"),YEAR='2013-2014')
bones_H = merge(fem_H,spn_H,by=c("SEQN","YEAR"))
bones_H = bones_H[,.(SEQN,DXXNKBMD,DXXOFBMD,DXXOSBMD,YEAR)]
bones_H[, check := rowSums(.SD), .SDcol=c('DXXNKBMD','DXXOSBMD','DXXOFBMD')]
bones_H = bones_H[!is.na(check)]
bones_H[,check:=NULL]

core_H = Reduce(function(...) merge(...,by=c("SEQN","YEAR")), list(demo_H,pfc_H,bones_H))

#################
# Build covariates dataset
# Download hormones female data with no NA (menopause)
HF = rbindlist(lapply(c('RHQ_D','RHQ_E','RHQ_F','RHQ_H'),
                      function(k) as.data.table(nhanes(paste0(k)))), fill = T)
HF = HF[,.(SEQN,RHQ060,RHQ305,RHQ310,RHQ540)]
HF[RHQ305==1 | RHQ310==1, hysterectomy:=1]
HF[RHQ305!=1 | RHQ310!=1, hysterectomy:=0]
HF = HF[,.(SEQN,age_last_period=RHQ060,hormones=RHQ540,hysterectomy)]

HF[, check := rowSums(.SD), .SDcol=c('age_last_period','hysterectomy')]
HF = HF[!is.na(check)]
HF[,check:=NULL]
HF = HF[age_last_period<100]

#################
# Download Activity data
ACT = rbindlist(lapply(c("PAQ_H"),
                       function(k) as.data.table(nhanes(paste0(k)))), fill = T)
ACT = ACT[,.(SEQN,PAQ650,PAQ605)]#,PAD320,PAD200)]
ACT[!is.na(PAQ650) | !is.na(PAQ605),act := ifelse(PAQ650==1 | PAQ605==1, 1, 0)]
ACT = ACT[,.(SEQN,act)]

#################
# Download smoking status data

SMQ = rbindlist(lapply(c('SMQ_H'),
                       function(k) as.data.table(nhanes(paste0(k)))), fill = T)
SMQ = SMQ[,.(SEQN,SMQ020,SMQ040)]
SMQ = SMQ[!is.na(SMQ020)]
SMQ[,smoking_status:=ifelse(SMQ020==2,0,ifelse(SMQ040==1 | SMQ040==2,2,1))]
SMQ = SMQ[!is.na(smoking_status),.(SEQN,smoking_status)]
SMQ[,smoking_status_ever_never:=ifelse(smoking_status==0,0,1)]

COV = Reduce(function(...) merge(...,by="SEQN", all = T), list(HF,ACT,SMQ))

X_H = merge(core_H,COV,by="SEQN")

X = rbindlist(list(X_H))
X = X[RIDAGEYR>=20]

X[,':='(hispanic = ifelse(RIDRETH1==1 | RIDRETH1==2,1,0),
        white = ifelse(RIDRETH1==3,1,0),
        black = ifelse(RIDRETH1==4 | RIDRETH1==5,1,0))]
X[,':='(smoker = ifelse(smoking_status==0,0,1))]
X[,mod_act:=act]

#################
# Definition of the weights
# clean weight outliers
OutVals = boxplot(X$WTMEC2YR)$out
plot(density(X$WTMEC2YR))
which(X$WTMEC2YR %in% OutVals)
X = X[-which(X$WTMEC2YR %in% OutVals)]
#################
#final dataset: X
#################

#################
# Define menopause
X[,menopause:=ifelse((RIDAGEYR-age_last_period)!=0 | hysterectomy==1,1,0)]

#################
# subset dataset
#men over 50 =115
X_man_50 = X[RIDAGEYR>=50 & RIAGENDR==1]
#postmenopausal women =117
X_women = X[RIAGENDR==2]
X_meno_hormones = X_women[menopause==1] 
X_meno_nohormones = X_women[menopause==1 & hormones==2] 

z = data.table(nhanes('DR1TOT_H'))
z = z[,.(SEQN,DR1TVD,DR1TCALC,DR1TVK)]
z = z[complete.cases(z)]

X_sensitivity = merge(X,z,by="SEQN")
X_sensitivity = X_sensitivity[is.na(hormones) | hormones==2]

X[, PFOS_aggregated := (0.7*SSNPFOS + 0.3*SSMPFOS)]
X[,RIAGENDR:=as.numeric(RIAGENDR-1)]
X_meno_nohormones[,RIAGENDR:=as.numeric(RIAGENDR-1)]
X_man_50[,RIAGENDR:=as.numeric(RIAGENDR-1)]
X_sensitivity[,RIAGENDR:=as.numeric(RIAGENDR-1)]
#################
# End
#################

# saveRDS(X_man_50,
#         "J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/PFC_bones_M.rds")
# saveRDS(X,
#         "J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/PFC_bones_20.rds")
# saveRDS(X_meno_nohormones,
#         "J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/PFC_bones_F.rds")
# saveRDS(X_sensitivity,
#         "J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/PFC_bones_S.rds")


### Sensitivity 2009-2010 Nhanes
#################
# Download demographic data with no NA
demo_F = data.table(nhanes('DEMO_F'),YEAR='2009-2010')
demo_F = demo_H[WTMEC2YR>0,.(SEQN,RIAGENDR,RIDAGEYR,WTMEC2YR,INDFMPIR,RIDRETH1,YEAR)]
demo_F[,SUM_MEC:=sum(demo_F$WTMEC2YR)]
demo_F[, check := rowSums(.SD), .SDcol=c('RIAGENDR','RIDAGEYR','INDFMPIR','RIDRETH1')]
demo_F = demo_F[!is.na(check)]
demo_F[,check:=NULL]

#################
# Download PFCs data with no NA
pfc_F = data.table(nhanes("PFC_F"),YEAR='2009-2010')
pfc_list <- c("LBXPFOA", "LBXPFOS", "LBXPFHS",  "LBXPFNA")
pfc_F = pfc_F[,c("SEQN","WTSC2YR",pfc_list,"YEAR"),with=F]
pfc_F[, check := rowSums(.SD), .SDcol=pfc_list]
pfc_F = pfc_F[!is.na(check)]
pfc_F[,check:=NULL]

#################
# Download bones mineral density data with no NA
fem_F = data.table(nhanes("DXXFEM_F"),YEAR='2009-2010')
spn_F = data.table(nhanes("DXXSPN_F"),YEAR='2009-2010')
bones_F = merge(fem_F,spn_F,by=c("SEQN","YEAR"))
bones_F = bones_F[,.(SEQN,DXXNKBMD,DXXOFBMD,DXXOSBMD,YEAR)]
bones_F[, check := rowSums(.SD), .SDcol=c('DXXNKBMD','DXXOSBMD','DXXOFBMD')]
bones_F = bones_F[!is.na(check)]
bones_F[,check:=NULL]


core_F = Reduce(function(...) merge(...,by=c("SEQN","YEAR")), list(demo_F,pfc_F,bones_F))

#################
# Build covariates dataset
# Download hormones female data with no NA (menopause)
HF = rbindlist(lapply(c('RHQ_D','RHQ_E','RHQ_F','RHQ_H'),
                      function(k) as.data.table(nhanes(paste0(k)))), fill = T)
HF = HF[,.(SEQN,RHQ060,RHQ305,RHQ310,RHQ540)]
HF[RHQ305==1 | RHQ310==1, hysterectomy:=1]
HF[RHQ305!=1 | RHQ310!=1, hysterectomy:=0]
HF = HF[,.(SEQN,age_last_period=RHQ060,hormones=RHQ540,hysterectomy)]

HF[, check := rowSums(.SD), .SDcol=c('age_last_period','hysterectomy')]
HF = HF[!is.na(check)]
HF[,check:=NULL]
HF = HF[age_last_period<100]

#################
# Download Activity data
ACT = rbindlist(lapply(c("PAQ_F"),
                       function(k) as.data.table(nhanes(paste0(k)))), fill = T)
ACT = ACT[,.(SEQN,PAQ650,PAQ605)]#,PAD320,PAD200)]
ACT[!is.na(PAQ650) | !is.na(PAQ605),act := ifelse(PAQ650==1 | PAQ605==1, 1, 0)]
ACT = ACT[,.(SEQN,act)]

#################
# Download smoking status data

SMQ = rbindlist(lapply(c('SMQ_F'),
                       function(k) as.data.table(nhanes(paste0(k)))), fill = T)
SMQ = SMQ[,.(SEQN,SMQ020,SMQ040)]
SMQ = SMQ[!is.na(SMQ020)]
SMQ[,smoking_status:=ifelse(SMQ020==2,0,ifelse(SMQ040==1 | SMQ040==2,2,1))]
SMQ = SMQ[!is.na(smoking_status),.(SEQN,smoking_status)]
SMQ[,smoking_status_ever_never:=ifelse(smoking_status==0,0,1)]

COV = Reduce(function(...) merge(...,by="SEQN", all = T), list(HF,ACT,SMQ))

X_F = merge(core_F,COV,by="SEQN")

X = rbindlist(list(X_F))
X = X[RIDAGEYR>=20]

X[,':='(hispanic = ifelse(RIDRETH1==1 | RIDRETH1==2,1,0),
        white = ifelse(RIDRETH1==3,1,0),
        black = ifelse(RIDRETH1==4 | RIDRETH1==5,1,0))]
X[,':='(smoker = ifelse(smoking_status==0,0,1))]
X[,mod_act:=act]

#################
# Definition of the weights
# clean weight outliers
OutVals = boxplot(X$WTMEC2YR)$out
plot(density(X$WTMEC2YR))
which(X$WTMEC2YR %in% OutVals)
#X = X[-which(X$WTMEC2YR %in% OutVals)]
#################
#final dataset: X
#################

#################
# Define menopause
X[,menopause:=ifelse((RIDAGEYR-age_last_period)!=0 | hysterectomy==1,1,0)]

#################
# subset dataset
#men over 50 =115
X_man_50 = X[RIDAGEYR>=50 & RIAGENDR==1]
#postmenopausal women =117
X_women = X[RIAGENDR==2]
X_meno_hormones = X_women[menopause==1] 
X_meno_nohormones = X_women[menopause==1 & hormones==2]

z = Reduce(function(...) merge(...,by="SEQN"), list(nhanes('DR1TOT_F'),nhanes('BMX_F'),nhanes('PbCd_F')))
z = data.table(z)
z = z[,.(SEQN,DR1TVD,DR1TCALC,DR1TVK,BMXBMI,LBXBPB,LBDBPBLC)]
z = z[complete.cases(z)]
z = z[LBDBPBLC == 0]

X_meno_nohormones_20092010 = merge(X_meno_nohormones,z,by="SEQN")
#saveRDS(X_meno_nohormones_20092010,"J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/PFC_bones_20092010.rds")



