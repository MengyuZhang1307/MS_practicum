
### This script is for model with lagged X,Y and without confounding

library(tidyverse)
library(tseries)
library(changepoint)
library(KFAS)
library(dlm)
library(matrixStats)


source("rcode/APTE/change_seed_nocon.R")
source("rcode/ts_gene_nocon.R")
source("rcode/APTE/fit_whole_getres.R")
source("rcode/APTE/split_whole.R")
source("rcode/APTE/fit_sep_getres.R")


################### parameter ###################

## m_t
m = c(500, 800, 400, 600, 700)

## probability
p = c(6/7, 2/7, 5/7, 3/7, 1/7)

## variance of random error
var_kt = 0.1
var_vt = 0.2
var_wt = 0.15
var_ut = 0.05
var_ut2 = 0.05

## coefficient beta_{0,t}, beta_{1,t}^0, mu_t^0
beta_0t = 2/2
beta_1t = c(1.5, 2.0, 1.5, 1.0, 0.5)
mu = c(2, 0.5, 1.5, 1, 0.5)
mu_2 = c(1.5, 0.2, 1, 0.5, 1)

## V_t^0
v = c(1.5, 1, 0.8, 1.2, 0.9)

## stationary phi_it, rho_it
phi_1 = c(0.2, 0.2, 0.1, 0.3, 0.1)/2
phi_2 = c(-0.3, -0.2, -0.25, -0.25, -0.3)/2
rho_1 = c(0.45, 0.5, 0.55, 0.5, 0.4)/2
rho_2 = c(0.6, 0.5, 0.45, 0.55, 0.5)/2
rho_3 = c(0.6, 0.5, 0.45, 0.55, 0.5)/2
rho_4 = c(0.2, 0.3, 0.25, 0.25, 0.15)/2

## U_t^0
u = c(1.8,1.1,1.2,0.5,1)
rho_5 = c(0.2, 0.3, 0.15, 0.35, 0.25)
var_et = 0.4
var_lt = 0.1

## gamma_(0,t)
gamma_0t = c(log(6), log(0.4), log(2.5), log(0.75), log(1/6))
gamma_1t = c(0.1,0.05,0.15,0.08,0.15)
gamma_2t = c(0.05,0.025,0.1,0.05,0.1)
gamma_3t = 0.2
gamma_4t = 0.5





################### simulation (main) ###################


all_res <- list()

j = 1
# while(j<=100) {
#   print(paste0("Simulation: ", j))
#   seed_set = sample.int(1e+6, 1)
#   all_res[[paste0("sim", j)]] = change_seed_nocon(seed_set, mod = "small")
#   j = j+1
# }

while(j<=100) {
  print(paste0("Simulation: ", j))
  tmp = try(change_seed_nocon(j, mod = "small"),TRUE)
  if(!inherits(tmp,"try-error")){
    all_res[[paste0("simulation", j)]] = tmp
    j = j+1
  }
}


APTE1 <- NULL
length1 <- NULL
for (i in 1:length(all_res)) {
  APTE1[i] = mean(all_res[[i]][[1]]$mean[-c(1:50)])
  length1[i] = length(all_res[[i]][[1]]$mean)
}
APTE2 <- NULL
length2 <- NULL
for (i in 1:length(all_res)) {
  APTE2[i] = mean(all_res[[i]][[2]]$mean)
  length2[i] = length(all_res[[i]][[2]]$mean)
}
APTE3 <- NULL
length3 <- NULL
for (i in 1:length(all_res)) {
  APTE3[i] = mean(all_res[[i]][[3]]$mean)
  length3[i] = length(all_res[[i]][[3]]$mean)
}
APTE4 <- NULL
length4 <- NULL
for (i in 1:length(all_res)) {
  APTE4[i] = mean(all_res[[i]][[4]]$mean)
  length4[i] = length(all_res[[i]][[4]]$mean)
}
APTE5 <- NULL
length5 <- NULL
for (i in 1:length(all_res)) {
  APTE5[i] = mean(all_res[[i]][[5]]$mean)
  length5[i] = length(all_res[[i]][[5]]$mean)
}
# lag
APTE_1 <- NULL
length_1 <-NULL
for (i in 1:length(all_res)) {
  APTE_1[i] = mean(all_res[[i]][[1]]$mean_lag[-c(1:50)])
  length_1[i] = length(all_res[[i]][[1]]$mean_lag)
}
APTE_2 <- NULL
length_2 <-NULL
for (i in 1:length(all_res)) {
  APTE_2[i] = mean(all_res[[i]][[2]]$mean_lag)
  length_2[i] = length(all_res[[i]][[2]]$mean_lag)
}
APTE_3 <- NULL
length_3 <-NULL
for (i in 1:length(all_res)) {
  APTE_3[i] = mean(all_res[[i]][[3]]$mean_lag)
  length_3[i] = length(all_res[[i]][[3]]$mean_lag)
}
APTE_4 <- NULL
length_4 <-NULL
for (i in 1:length(all_res)) {
  APTE_4[i] = mean(all_res[[i]][[4]]$mean_lag)
  length_4[i] = length(all_res[[i]][[4]]$mean_lag)
}
APTE_5 <- NULL
length_5 <-NULL
for (i in 1:length(all_res)) {
  APTE_5[i] = mean(all_res[[i]][[5]]$mean_lag)
  length_5[i] = length(all_res[[i]][[5]]$mean_lag)
}


APTE = cbind(APTE1,APTE2,APTE3,APTE4,APTE5)
APTE_lag = cbind(APTE_1,APTE_2,APTE_3,APTE_4,APTE_5)

cbind(Length = c(mean(length1),mean(length2), mean(length3), mean(length4), mean(length5)), 
      m_t = m,
      Estimated = colMeans(APTE), True = mu, Variance = colVars(APTE))
cbind(Length = c(mean(length_1),mean(length_2), mean(length_3), mean(length_4), mean(length_5)), 
      m_t = m,
      estimated = colMeans(APTE_lag), true = mu_2, variance = colVars(APTE_lag))

# ####################
# 
# # alpha for 12 model and different seeds
# table_alpha <- NULL
# for (i in 1:12) {
#   temp = as_table(means,length(seed_set),i)
#   colnames(temp) <- paste0("seed_",seed_set)
#   table_alpha[[names(means[[1]][i])]] = temp
# }
# 
# # APTE for each seed 
# APTE = map(table_alpha, ~colMeans(.x))
# 
# 
# 
# miss_alpha = table_alpha$mis_whole
# 
# 
# 
# change_num <-NULL
# for (i in seq(1,200, 1)) {
#   cptm_stationary <- cpt.mean(miss_alpha[,1], penalty='Asymptotic', pen.value = 0.05, method='PELT')
#   change_num[i] <- length(cpts(cptm_stationary)) # number of point
#   #plot(cptm_stationary)
# }
# 
# plot(seq(1,200, 1), y = change_num)
# lines(predict(loess(change_num~seq(1,200, 1))), col='red', lwd=2)
# 
# 
# cptm_stationary <- cpt.mean(miss_alpha[,1], penalty='Asymptotic', pen.value =0.05, method='PELT')
# cpts(cptm_stationary)
# # 500 1300 1700 2300 3000
# 

# ################### generate data (for debug) ###################
# source("rcode/ts_gene_nocon.R")
# data = ts_gene_nocon(seed = 3)
# 
# ################### fit whole time series model and check periods (for debug) ###################
# 
# # fit whole model
# # source("rcode/fit_whole_getres.R")
# whole_res = fit_whole_getres(data, TRUE, model = "large")
# 
# # split data
# split_data = split_whole(data, whole_res)
# 
# 
# # var_est = 0
# # data = test
# # model.nocon.noxy <- SSModel(data$df$Y ~ SSMregression(~ data$df$V + data$df$X + data$df$X1 + data$df$Y1 + data$df$Y2, type = "distinct", Q = diag(c(rep(var_est,3),0,0))), H = NA)
# # res.nocon.noxy = ssm_fit(model.nocon.noxy,6)
# 
# ################### fit seperate model (for debug) ###################
# 
# # fit whole model
# # source("rcode/fit_sep_getres.R")
# sep_res = fit_sep_getres(split_data,FALSE,model = "large")
# 
# 
# # var_est = 0
# # model1.nocon.xy <- SSModel(split_data$Y[[1]] ~ SSMregression(~ data$V[[1]] + data$X[[1]] + data$X1[[1]] + data$Y1[[1]] + data$Y2[[1]], type = "distinct", Q = diag(c(rep(var_est,3),0,0))), H = NA)
# #
# # model2.nocon.xy <- SSModel(split_data$Y[[2]] ~ SSMregression(~ data$V[[2]] + data$X[[2]] + data$X1[[2]] + data$Y1[[2]] + data$Y2[[2]], type = "distinct", Q = diag(c(rep(var_est,3),0,0))), H = NA)
# #
# # model3.nocon.xy <- SSModel(split_data$Y[[3]] ~ SSMregression(~ data$V[[3]] + data$X[[3]] + data$X1[[3]] + data$Y1[[3]] + data$Y2[[3]], type = "distinct", Q = diag(c(rep(var_est,3),0,0))), H = NA)
# #
# # model4.nocon.xy <- SSModel(split_data$Y[[4]] ~ SSMregression(~ data$V[[4]] + data$X[[4]] + data$X1[[4]] + data$Y1[[4]] + data$Y2[[4]], type = "distinct", Q = diag(c(rep(var_est,3),0,0))), H = NA)
# #
# # model5.nocon.xy <- SSModel(split_data$Y[[5]] ~ SSMregression(~ data$V[[5]] + data$X[[5]] + data$X1[[5]] + data$Y1[[5]] + data$Y2[[5]], type = "distinct", Q = diag(c(rep(var_est,3),0,0))), H = NA)
# #
# # res1.nocon.xy = ssm_fit(model1.nocon.xy,6)
# # res2.nocon.xy = ssm_fit(model2.nocon.xy,6)
# # res3.nocon.xy = ssm_fit(model3.nocon.xy,6)
# # res4.nocon.xy = ssm_fit(model4.nocon.xy,6)
# # res5.nocon.xy = ssm_fit(model5.nocon.xy,6)
# 
# # ################### get APTE ###################
# # # with 50 burn out time points
# #
# # for (t in 1:length(m)) {
# #   est_APTE[t] = mean(sep_res[[t]]$mean[-c(1:50)])
# # }
# # est_APTE
# #
# # cbind(estimated = est_APTE, true = mu)
