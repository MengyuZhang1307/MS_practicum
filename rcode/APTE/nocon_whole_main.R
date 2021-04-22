# this script only for whole series 

library(tidyverse)
library(tseries)
library(changepoint)
library(KFAS)
library(dlm)
library(matrixStats)





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

##################

source("rcode/APTE/fit_whole_getres.R")
source("rcode/ts_gene_nocon.R")
source("rcode/APTE/split_whole.R")
source("rcode/APTE/fit_sep_getres.R")

simulation <- function(Seed) {
  
  # generate data
  data = ts_gene_nocon(seed = Seed) 
  
  # fit whole times series
  whole_res = fit_whole_getres(data, variance = TRUE)
  
  return (whole_res)
  
}

seed_set = 1:10
all_res <- list()
for (i in seed_set) {
  print(paste0("set.seed", i))
  all_res[[paste0("set.seed", i)]] = simulation(i)
}




APTE_hat <- NULL
for (i in 1:length(seed_set)) {
   # change point time points
  cpts_stationary = rep(0,6)
  j = 4
  while (length(cpts_stationary) > 4 & j<100){
    cptm_stationary <- cpt.mean(all_res[[i]]$mean[-c(1:50)], penalty='Manual', pen.value = j, method='PELT')
    cpts_stationary <- cpts(cptm_stationary) # change point
    j = j+2
  }
  # plot(cptm_stationary)
  APTE_hat[[i]] = cbind(c(cpts_stationary+50,3000),cptm_stationary@param.est$mean)
}

