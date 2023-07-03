###############################
# libraries
###############################
library(data.table)
library(BWQS)
library(rstan)
library(parallel)

###############################
#BWQS regression for the overall population, men and women
###############################
model_bwqs_stan_w = "data {
int<lower=0> N;          // number of individual
int<lower=0> C;          // number of chemicals
int<lower=0> K;          // number of covariates
matrix[N,C] X;		       // matrix of indipendent variable
matrix[N,K] KV;		       // matrix of covariates
vector[C] Dalp;          // vector of the Dirichlet coefficients
vector[N] sw;            // sampling weights
real y[N];               // outcome continuos variable
}
parameters {
real beta0;              // intercepts
real beta1;              // overall effect
vector[K] delta;         // covariates coefficients
simplex[C] W;            // weights
real<lower=0> sigma;     // standard deviation of the model
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + beta1*(X*W) + KV*delta;
}
model {
beta0 ~ normal(0, 100);
beta1 ~ normal(0, 100);
delta ~ normal(0, 100);
W ~ dirichlet(Dalp);
sigma ~ inv_gamma(0.01, 0.01);

for(n in 1:N){
  target +=  normal_lpdf(y[n]| Xb[n], sigma) * sw[n];
}
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = normal_lpdf(y[nn]| Xb[nn], sigma);

}
"

dataset = list.files(path = "J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/",
           pattern = "PFC_bones", full.names = T)[-4]
# results = list.files(path = "J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/",
#                      pattern = "_w.rds",full.names = T)[4:12]
input_list = data.frame(dt = as.character(rep(dataset[c(1,3,2)],each = 3)),
                        output = rep(c("DXXNKBMD","DXXOFBMD","DXXOSBMD"),3))
covariates = c("RIDAGEYR","mod_act","INDFMPIR","hispanic","black","smoker","RIAGENDR")
input_list$covariates[1:3] = I(list(covariates))
input_list$covariates[4:9] = I(list(covariates[-7]))
input_list$title = c('all_femurneck','all_femur','all_spine',
                   'men_femurneck','men_femur','men_spine',
                   'women_femurneck','women_femur','women_spine')

no_cores = 10 #detectCores()/2
cl <- makeCluster(no_cores)

clusterEvalQ(cl, {
  library(BWQS)
  library(data.table)
})

clusterExport(cl, varlist = c("input_list","model_bwqs_stan_w"))

tmp = parLapplyLB(cl,1:9,function(i){
  cat(paste0("initialize model: ",input_list$title[i],"\n"))
  dt = readRDS(as.character(input_list$dt[i]))
  output = as.character(input_list$output[i])
  mean_w = mean(dt$WTMEC2YR)
  dt[, weight_new := WTMEC2YR/mean_w]
  chem = data.table(new_name=c('n-PFOA','Sb-PFOA','n-PFOS','Sm-PFOS','PFHxS','Me-PFOSA-AcOH','PFDeA','PFBuS',
                               'PFHpA','PFNA','PFUA','PFDoA'),
                    old_name = c('NPFOA','BPFOA','NPFOS','MPFOS','PFHS','MPAH',
                                 'PFDE','PFBS','PFHP','PFNA','PFUA','PFDO'),
                    label = c('SSNPFOA','SSBPFOA','SSNPFOS','SSMPFOS',
                              'LBXPFHS','LBXMPAH','LBXPFDE',
                              'LBXPFBS','LBXPFHP','LBXPFNA','LBXPFUA','LBXPFDO'),
                    include = c(1,0,1,1,1,1,1,0,0,1,1,0))
  chemical_name = chem[include == 1, label]
  X = quantile_split(as.data.frame(dt[,chemical_name,with=F]),
                     mix_name = chemical_name,
                     q = 4)
  data_reg <- list(
    N = NROW(dt),
    C = length(chemical_name),
    K = length(input_list$covariates[[i]]),
    X = X,
    KV = dt[,input_list$covariates[[i]],with=F],
    sw = as.vector(dt$weight_new),
    Dalp = rep(1,length(chemical_name)),
    y = as.vector(unlist(dt[,output,with=F]))
  )
  
  iter = 20000
  thin = 2
  fit <- stan(model_code = model_bwqs_stan_w,
              data = data_reg,
              chains = 1,
              warmup = iter/2,
              iter = iter,
              cores = 1,
              thin = thin,
              refresh = 0,
              algorithm = 'NUTS',
              control=list(max_treedepth = 20,
                           adapt_delta = 0.999999999999999))
})

saveRDS(tmp,"J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/last_results_weights.rds")

for(i in 1:9){
  if(i %in% c(1,2,3)){
    fit = summary(tmp[[i]])$summary[c(1,2,10:17),c(1,3,4,8,9,10)]
    fit = trimws(format(round(fit,3),nsmall=3))
    CRI = paste0("(",t(fit[,3]),"; ",t(fit[,4]),")")
    fit = cbind(Estimate = fit[,1],CRI,Rhat = fit[,6])
    rownames(fit) = c("beta0","beta1",chem[include==1,new_name])
    cat(paste0(input_list$title[i],'\n'))
    write.csv(fit,paste0("J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/new_weight_",
                         input_list$title[i],".csv"))
  } else {
    fit = summary(tmp[[i]])$summary[c(1,2,9:16),c(1,3,4,8,9,10)]
    fit = trimws(format(round(fit,3),nsmall=3))
    CRI = paste0("(",t(fit[,3]),"; ",t(fit[,4]),")")
    fit = cbind(Estimate = fit[,1],CRI,Rhat = fit[,6])
    rownames(fit) = c("beta0","beta1",chem[include==1,new_name])
    cat(paste0(input_list$title[i],'\n'))
    write.csv(fit,paste0("J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/new_weight_",
                         input_list$title[i],".csv"))
  }
}



###############################
#WQS regression
###############################
library(gWQS)
library(Bolstad)
library(asbio)

X = readRDS("J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111//PFC_bones_20.rds")
mean_w = mean(X$WTMEC2YR)
X[, weight_new := WTMEC2YR/mean_w]
X = data.frame(X)
chem = data.table(new_name=c('n-PFOA','Sb-PFOA','n-PFOS','Sm-PFOS','PFHxS','Me-PFOSA-AcOH','PFDeA','PFBuS',
                             'PFHpA','PFNA','PFUA','PFDoA'),
                  old_name = c('NPFOA','BPFOA','NPFOS','MPFOS','PFHS','MPAH',
                               'PFDE','PFBS','PFHP','PFNA','PFUA','PFDO'),
                  label = c('SSNPFOA','SSBPFOA','SSNPFOS','SSMPFOS',
                            'LBXPFHS','LBXMPAH','LBXPFDE',
                            'LBXPFBS','LBXPFHP','LBXPFNA','LBXPFUA','LBXPFDO'),
                  include = c(1,0,1,1,1,1,1,0,0,1,1,0))
chemical_name = chem[include == 1, label]

results_spine = gwqs(DXXOSBMD ~ wqs + RIDAGEYR + mod_act + INDFMPIR + hispanic + black + smoker + RIAGENDR,
               mix_name = chemical_name, data = X, q = 4, validation = 0.7,
               b = 100, b1_pos = FALSE, b1_constr = FALSE, family = gaussian, 
               weights = "weight_new",
               seed = 2019)

summary(results_spine$fit)
results_spine$final_weights

results_fn = gwqs(DXXNKBMD ~ wqs + RIDAGEYR + mod_act + INDFMPIR + hispanic + black + smoker + RIAGENDR,
               mix_name = chemical_name, data = X, q = 4, validation = 0.7,
               b = 100, b1_pos = FALSE, b1_constr = FALSE, family = gaussian, 
               weights = "weight_new",
               seed = 2019)

summary(results_fn$fit)
results_fn$final_weights

results_femur = gwqs(DXXOFBMD ~ wqs + RIDAGEYR + mod_act + INDFMPIR + hispanic + black + smoker + RIAGENDR,
               mix_name = chemical_name, data = X, q = 4, validation = 0.7,
               b = 100, b1_pos = FALSE, b1_constr = FALSE, family = gaussian, 
               weights = "weight_new",
               seed = 2019)

BB = rbind(c(Estimate = summary(results_spine$fit)$coefficients[2,1],
        confint(results_spine$fit)[2,]),
      c(Estimate = summary(results_femur$fit)$coefficients[2,1],
        confint(results_femur$fit)[2,]),
      c(Estimate = summary(results_fn$fit)$coefficients[2,1],
        confint(results_fn$fit)[2,]))
w_s  = data.table(t(results_spine$final_weights))[2]
w_f  = data.table(t(results_femur$final_weights))[2]
w_fn = data.table(t(results_fn$final_weights))[2]

w = rbindlist(list(w_s,w_f,w_fn), use.names = T)
colnames(w) = chem$new_name[c(6,11,5,7,3,10,1,4)]

###############################
# Individual analysis for chemicals
###############################
model_single_stan_w = "data {
int<lower=0> N;          // number of individual
int<lower=0> K;          // number of covariates
matrix[N,K+1] X;		     // matrix of indipendent variable
vector[N] sw;            // sampling weights
real y[N];               // outcome continuos variable
}
parameters {
real beta0;              // intercepts
vector[K+1] beta;         // covariates coefficients
real<lower=0> sigma;     // standard deviation of the model
}
transformed parameters {
vector[N] Xb;
Xb = beta0 + X*beta;
}
model {
beta0 ~ normal(0, 100);
beta ~ normal(0, 100);
sigma ~ inv_gamma(0.01, 0.01);

for(n in 1:N){
  target +=  normal_lpdf(y[n]| Xb[n], sigma) * sw[n];
}
}
generated quantities {
vector[N] log_lik;
for (nn in 1:N)
log_lik[nn] = normal_lpdf(y[nn]| Xb[nn], sigma);

}
"

input_list = expand.grid(output = c("DXXNKBMD","DXXOFBMD","DXXOSBMD"),
                chemicals = chem[include == 1,label])
covariates = c("RIDAGEYR","mod_act","INDFMPIR","hispanic","black","smoker","RIAGENDR")
input_list$covariates = I(list(covariates))


no_cores = 24 #detectCores()/2
cl <- makeCluster(no_cores)

clusterEvalQ(cl, {
  library(BWQS)
  library(data.table)
})

clusterExport(cl, varlist = c("input_list","model_single_stan_w"))

tmp = parLapplyLB(cl,1:24,function(i){
  dt = readRDS("J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/PFC_bones_20.rds")
  output = as.character(input_list$output[i])
  mean_w = mean(dt$WTMEC2YR)
  dt[, weight_new := WTMEC2YR/mean_w]
  chem = data.table(new_name=c('n-PFOA','Sb-PFOA','n-PFOS','Sm-PFOS','PFHxS','Me-PFOSA-AcOH','PFDeA','PFBuS',
                               'PFHpA','PFNA','PFUA','PFDoA'),
                    old_name = c('NPFOA','BPFOA','NPFOS','MPFOS','PFHS','MPAH',
                                 'PFDE','PFBS','PFHP','PFNA','PFUA','PFDO'),
                    label = c('SSNPFOA','SSBPFOA','SSNPFOS','SSMPFOS',
                              'LBXPFHS','LBXMPAH','LBXPFDE',
                              'LBXPFBS','LBXPFHP','LBXPFNA','LBXPFUA','LBXPFDO'),
                    include = c(1,0,1,1,1,1,1,0,0,1,1,0))
  data_reg <- list(
    N = NROW(dt),
    K = length(input_list$covariates[[i]]),
    X = dt[,c(as.character(input_list$chemicals[i]),
              input_list$covariates[[i]]),with=F],
    sw = as.vector(dt$weight_new),
    y = as.vector(unlist(dt[,output,with=F]))
  )
  
  iter = 20000
  thin = 2
  fit <- stan(model_code = model_single_stan_w,
              data = data_reg,
              chains = 1,
              warmup = iter/2,
              iter = iter,
              cores = 1,
              thin = thin,
              refresh = 0,
              algorithm = 'NUTS',
              control=list(max_treedepth = 20,
                           adapt_delta = 0.999999999999999))
})

saveRDS(tmp,"J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/single_analysis_weights.rds")

table_results = lapply(1:24,function(i){
  data.table(t(c(Output = as.character(input_list$output[i]),
    PFC = as.character(input_list$chemicals[i]),
    round(summary(tmp[[i]])$summary[2,c(1,4,8)],3))))
})

table_results = rbindlist(table_results)

##########################################################################
# Sensitivity analysis: BWQS exclude 60% LOD and add calcium + vitamin (D & K) #
##########################################################################

dataset = list.files(path = "J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/",
                     pattern = "PFC_bones", full.names = T)[c(1,4)]
#dk = readRDS(dataset)
# dt = dt[,-c("DS1IVD","DS1ICALC","DS1IVK"),with=F]
# intake = data.table(nhanes("DR1TOT_H"))
# intake = intake[,.(SEQN,DR1TVK,DR1TVD,DR1TCALC)]
# dt = merge(dt,intake,by="SEQN")
# results = list.files(path = "J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/",
#                      pattern = "_w.rds",full.names = T)[4:12]
input_list = data.frame(dt = c(as.character(rep(dataset[1],6)),as.character(rep(dataset[2],3))),
                        output = rep(c("DXXNKBMD","DXXOFBMD","DXXOSBMD"),3))
covariates = c("RIDAGEYR","mod_act","INDFMPIR","hispanic","black","smoker","RIAGENDR",
               "DR1TVD","DR1TCALC","DR1TVK")
input_list$covariates[7:9] = I(list(covariates))
input_list$covariates[1:6] = I(list(covariates[-c(8,9,10)]))

input_list$title = c('all_femurneck','all_femur','all_spine')
chem = data.table(new_name=c('n-PFOA','Sb-PFOA','n-PFOS','Sm-PFOS','PFHxS','Me-PFOSA-AcOH','PFDeA','PFBuS',
                             'PFHpA','PFNA','PFUA','PFDoA'),
                  old_name = c('NPFOA','BPFOA','NPFOS','MPFOS','PFHS','MPAH',
                               'PFDE','PFBS','PFHP','PFNA','PFUA','PFDO'),
                  label = c('SSNPFOA','SSBPFOA','SSNPFOS','SSMPFOS',
                            'LBXPFHS','LBXMPAH','LBXPFDE',
                            'LBXPFBS','LBXPFHP','LBXPFNA','LBXPFUA','LBXPFDO'),
                  include = c(1,0,1,1,1,1,1,0,0,1,1,0))
chemical_name = chem[include == 1, label]
input_list$mix_metals[1:3] = I(list(chem$label[c(1,3,4,5,7,10)]))
input_list$mix_metals[4:6] = I(list(c("PFOS_aggregated",chemical_name[-c(2,3)]))) 
input_list$mix_metals[7:9] = I(list(chemical_name))

no_cores = 9 #detectCores()/2
cl <- makeCluster(no_cores)

clusterEvalQ(cl, {
  library(BWQS)
  library(data.table)
})

clusterExport(cl, varlist = c("input_list","model_bwqs_stan_w"))

tmp = parLapplyLB(cl,1:9,function(i){
  dt = readRDS(as.character(input_list$dt[i]))
  dt = dt[complete.cases(dt[,input_list$covariates[[i]],with=F])]
  output = as.character(input_list$output[i])
  mean_w = mean(dt$WTMEC2YR)
  dt[, weight_new := WTMEC2YR/mean_w]
  X = quantile_split(as.data.frame(dt[,input_list$mix_metals[[i]],with=F]),
                     mix_name = input_list$mix_metals[[i]],
                     q = 4)
  data_reg <- list(
    N = NROW(dt),
    C = length(input_list$mix_metals[[i]]),
    K = length(input_list$covariates[[i]]),
    X = X,
    KV = dt[,input_list$covariates[[i]],with=F],
    sw = as.vector(dt$weight_new),
    Dalp = rep(1,length(input_list$mix_metals[[i]])),
    y = as.vector(unlist(dt[,output,with=F]))
  )
  
  iter = 20000
  thin = 2
  fit <- stan(model_code = model_bwqs_stan_w,
              data = data_reg,
              chains = 1,
              warmup = iter/2,
              iter = iter,
              cores = 1,
              thin = thin,
              refresh = 0,
              algorithm = 'NUTS',
              control=list(max_treedepth = 20,
                           adapt_delta = 0.999999999999999))
})

saveRDS(tmp,"J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/last_results_weights_sensitivity.rds")

for(i in 1:9){
  if(i %in% c(1,2,3)){
    fit = summary(tmp[[i]])$summary[c(1,2,10:15),c(1,4,8)]
    fit = trimws(format(round(fit,3),nsmall=3))
    CRI = paste0("(",t(fit[,2]),"; ",t(fit[,3]),")")
    fit = cbind(Estimate = fit[,1],CRI)
    rownames(fit) = c("beta0","beta1",chem[label %in% input_list$mix_metals[[i]],new_name])
    cat(paste0(input_list$title[i],'\n'))
    write.csv(fit,paste0("J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/sensitivity_60_",
                         input_list$title[i],".csv"))
  } else if(i %in% c(4,5,6)){
    fit = summary(tmp[[i]])$summary[c(1,2,10:16),c(1,4,8)]
    fit = trimws(format(round(fit,3),nsmall=3))
    CRI = paste0("(",t(fit[,2]),"; ",t(fit[,3]),")")
    fit = cbind(Estimate = fit[,1],CRI)
    rownames(fit) = c("beta0","beta1",c("PFOS",chem[label %in% input_list$mix_metals[[i]],new_name]))
    cat(paste0(input_list$title[i],'\n'))
    write.csv(fit,paste0("J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/sensitivity_PFOS_",
                         input_list$title[i],".csv"))
  } else {
    fit = summary(tmp[[i]])$summary[c(1,2,13:20),c(1,4,8)]
    fit = trimws(format(round(fit,3),nsmall=3))
    CRI = paste0("(",t(fit[,2]),"; ",t(fit[,3]),")")
    fit = cbind(Estimate = fit[,1],CRI)
    rownames(fit) = c("beta0","beta1",c(chem[label %in% input_list$mix_metals[[i]],new_name]))
    cat(paste0(input_list$title[i],'\n'))
    write.csv(fit,paste0("J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/sensitivity_diet_",
                         input_list$title[i],".csv"))
  }
}

## Sensitivity cicle 2009-2010 
## Postmenopausal women only (remove the women with hormones replacement)

dataset = "J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/PFC_bones_20092010.rds"
input_list = data.frame(dt = as.character(rep(dataset,3)),
                        output = c("DXXNKBMD","DXXOFBMD","DXXOSBMD"))
covariates = c("RIDAGEYR","mod_act","INDFMPIR","hispanic","black","smoker",
               "DR1TVD","DR1TCALC","DR1TVK","BMXBMI","LBXBPB")
input_list$covariates[1:3] = I(list(covariates))

input_list$title = c('all_femurneck','all_femur','all_spine')
chem = data.table(new_name=c('PFOA','PFOS','PFHxS','PFNA'),
                  old_name = c('PFOA','PFOS','PFHS','PFNA'),
                  label = c('LBXPFOA','LBXPFOS',
                            'LBXPFHS','LBXPFNA'),
                  include = c(1,1,1,1))
chemical_name = chem[include == 1, label]
input_list$mix_metals[1:3] = I(list(chemical_name))

no_cores = 3 #detectCores()/2
cl <- makeCluster(no_cores)

clusterEvalQ(cl, {
  library(BWQS)
  library(data.table)
})

clusterExport(cl, varlist = c("input_list","model_bwqs_stan_w"))

tmp = parLapplyLB(cl,1:3,function(i){
  dt = readRDS(as.character(input_list$dt[i]))
  dt = dt[complete.cases(dt[,input_list$covariates[[i]],with=F])]
  output = as.character(input_list$output[i])
  mean_w = mean(dt$WTMEC2YR)
  dt[, weight_new := WTMEC2YR/mean_w]
  X = quantile_split(as.data.frame(dt[,input_list$mix_metals[[i]],with=F]),
                     mix_name = input_list$mix_metals[[i]],
                     q = 4)
  data_reg <- list(
    N = NROW(dt),
    C = length(input_list$mix_metals[[i]]),
    K = length(input_list$covariates[[i]]),
    X = X,
    KV = dt[,input_list$covariates[[i]],with=F],
    sw = as.vector(dt$weight_new),
    Dalp = rep(1,length(input_list$mix_metals[[i]])),
    y = as.vector(unlist(dt[,output,with=F]))
  )
  
  iter = 20000
  thin = 2
  fit <- stan(model_code = model_bwqs_stan_w,
              data = data_reg,
              chains = 1,
              warmup = iter/2,
              iter = iter,
              cores = 1,
              thin = thin,
              refresh = 0,
              algorithm = 'NUTS',
              control=list(max_treedepth = 20,
                           adapt_delta = 0.999999999999999))
})

saveRDS(tmp,"J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/last_results_weights_20092010.rds")

for(i in 1:3){
    fit = summary(tmp[[i]])$summary[c(1,2,14:17),c(1,4,8)]
    fit = trimws(format(round(fit,3),nsmall=3))
    CRI = paste0("(",t(fit[,2]),"; ",t(fit[,3]),")")
    fit = cbind(Estimate = fit[,1],CRI)
    rownames(fit) = c("beta0","beta1",chem[label %in% input_list$mix_metals[[i]],new_name])
    cat(paste0(input_list$title[i],'\n'))
    write.csv(fit,paste0("J:/PM/Colicino_Lab/Manuscripts/2019_BWQS_PFCs_Bone/Results_20210111/sensitivity_20092010_",
                         input_list$title[i],".csv"))
}


