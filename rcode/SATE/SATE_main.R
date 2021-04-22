
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
source("rcode/SATE/change_seed_sate.R")

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
mu_2 = c(1.5, 0.2, 1, 0.5, 0.2)

## V_t^0
v = c(1.5, 1, 0.8, 1.2, 0.9)

## stationary phi_it, rho_it
phi_1 = c(0.2, 0.2, 0.1, 0.3, 0.1)/2
phi_2 = c(-0.3, -0.2, -0.25, -0.25, -0.3)/2
rho_1 = c(0.45, 0.5, 0.55, 0.5, 0.4)/2
rho_2 = c(0.6, 0.5, 0.45, 0.55, 0.5)/2
rho_3 = c(0.6, 0.5, 0.45, 0.55, 0.5)/2
rho_4 = c(0.2, 0.3, 0.25, 0.25, 0.15)/2



seed_set = 3
all_sate <- list()
for (i in seed_set) {
  print(paste0("set.seed", i))
  all_sate[[paste0("set.seed", i)]] = change_seed_sate(i, mod = "large")
}


################### get tau (for debug) ###################

# period 1
period = 1
model = sep_res[[period]]
B = 100


m_split = map(split_data$Y, ~length(.x)) %>% unlist # number of time points in each period
Z = split_data
Z$beta0 = map(m_split, ~rep(1,.x))
Z_tilda = Z
Z_tilda$X = map(split_data$X, ~(1-.x))


z_tilda = cbind(beta0 = Z_tilda$beta0[[period]],
                V = Z_tilda$V[[period]],
                X = Z_tilda$X[[period]],
                X1 = Z_tilda$X1[[period]],
                Y1 = Z_tilda$Y1[[period]],
                Y2 = Z_tilda$Y2[[period]]) %>% as.matrix()

alpha_hat = model$fitres$alphahat %>% as.matrix() # smoothing state

alpha_telda = rowSums(z_tilda * alpha_hat) # mean of counterfactuals at each time point

V_telda <- vector() # variance of counterfactuals at each time point
for (i in 1:m_split[period]) {
  V_telda[i] = z_tilda[i,] %*% model$fitres$V[,,i] %*% z_tilda[i,] + var_kt
}

y_telda = data.frame()
set.seed(1)
for (i in 1:m_split[period]) {
  temp = rnorm(B, alpha_telda[i], V_telda[i])
  y_telda = rbind(y_telda, temp)
}
colnames(y_telda) <- paste0("B_", 1:B)

tau = (y_telda-split_data$Y[[period]])*(z_tilda[,match("X",colnames(z_tilda))]-Z$X[[period]])

APTE_sate = colMeans(tau)
APTE_sate

