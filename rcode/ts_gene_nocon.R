# func.3: data generation function reutrn V,X,Y, and df for all variables (no confounding)
source("rcode/newyy.R")
ts_gene_nocon <- function (methods = "cons_var", seed) {
  

  
  set.seed(seed)
  # V: simultanous cause of Y
  for (t in 1:5){
    V_temp = rnorm(m[t], v[t], sqrt(var_vt/(1-rho_4[t]^2)))
    error_v = rnorm(m[t], 0, sqrt(var_vt)) 
    for (j in 2:m[t]) V_temp[j] = v[t]*(1-rho_4[t]) + rho_4[t]*V_temp[j-1] + error_v[j]
    assign(paste0("V_", t), V_temp)
    #plot(V_temp, type = "l", main = c("period",t))
    #print(adf.test(V_temp)$p.value)
  }
  V.1 = list(V_1 = V_1, V_2 = V_2, V_3 = V_3, V_4 = V_4, V_5 = V_5, whole = c(V_1,V_2,V_3,V_4,V_5))
  
  # exposure without confounding
  X.1 = NULL
  for (t in 1:length(m)) {
    X_temp = as.numeric(runif(m[t]) < p[t])
    assign(paste0("X_", t), X_temp)
    #plot(X_temp, type = "l", main = c("period",t))
    #print(adf.test(X_temp)$p.value)
    #print(c(p[t],mean(X_temp)))
  }
  X.1 = list(X_1 = X_1, X_2 = X_2, X_3 = X_3, X_4 = X_4, X_5 = X_5, whole = c(X_1,X_2,X_3,X_4,X_5))
  
  # beta: effect between V and Y
  if (methods == "cons_var") {
    for (t in 1:length(m)){
      beta_temp = rnorm(m[t], beta_1t[t], sqrt(var_wt))
      assign(paste0("beta_", t), beta_temp)
      #plot(beta_temp, type = "l", main = c("period",t))
      #print(adf.test(beta_temp)$p.value)
    }
  }
  if (methods == "cons_novar") {
    for (t in 1:length(m)){
      assign(paste0("beta_", t), rep(beta_1t[t], m[t]))
    }
  }
  beta_t.1 = list(beta_1 = beta_1, beta_2 = beta_2, beta_3 = beta_3, beta_4 = beta_4, beta_5 = beta_5, whole = c(beta_1,beta_2,beta_3,beta_4,beta_5))
  
  # mu_1: causal effect between Xt and Y
  if (methods == "cons_var") {
    for (t in 1:5){
      mu_temp = rnorm(m[t], mu[t], sqrt(var_ut))
      assign(paste0("mu_1", t), mu_temp)
      #plot(mu_temp, type = "l", main = c("period",t))
      #print(adf.test(mu_temp)$p.value)
    }
  } 
  if (methods == "cons_novar") {
    for (t in 1:length(m)){
      assign(paste0("mu_1", t), rep(mu[t], m[t]))
    }
  }
  mu_1t.1 = list(mu_11 = mu_11, mu_12 = mu_12, mu_13 = mu_13, mu_14 = mu_14, mu_15 = mu_15, whole =  c(mu_11,mu_12,mu_13,mu_14,mu_15))
  
  # mu_1: causal effect between Xt and Y
  if (methods == "cons_var") {
    for (t in 1:5){
      mu_temp = rnorm(m[t], mu_2[t], sqrt(var_ut2))
      assign(paste0("mu_2", t), mu_temp)
      #plot(mu_temp, type = "l", main = c("period",t))
      #print(adf.test(mu_temp)$p.value)
    }
  } 
  if (methods == "cons_novar") {
    for (t in 1:length(m)){
      assign(paste0("mu_2", t), rep(mu_2[t], m[t]))
    }
  }
  mu_2t.1 = list(mu_21 = mu_21,mu_22 = mu_22,mu_23 = mu_23,mu_24 = mu_24,mu_25 = mu_25, whole = c(mu_21,mu_22,mu_23,mu_24,mu_25))
  
  # Y 
  for (t in 1:length(m)){
    assign(paste0("data1", t), newyy(t,V.1[[t]], mu_1t.1[[t]], mu_2t.1[[t]], X.1[[t]], beta_t.1[[t]]))
  }
  Y.1 = list(Y_1 = data11$Y.new, Y_2 = data12$Y.new, Y_3 = data13$Y.new, Y_4 = data14$Y.new, Y_5 = data15$Y.new, whole = c(data11$Y.new, data12$Y.new, data13$Y.new,data14$Y.new, data15$Y.new))
  
  # all data with two lags of Y and one lag of X
  df = list(Y = Y.1$whole, 
             V = V.1$whole, 
             X = X.1$whole, 
             Y1 = lag(Y.1$whole, n = 1, default = Y.1$whole[1]),
             Y2 = lag(Y.1$whole, n = 2, default = Y.1$whole[1]),
             X1 = lag(X.1$whole, n = 1, default = X.1$whole[1])) 

  
  return(list(V.1 = V.1, mu_1t.1 = mu_1t.1, mu_2t.1 = mu_2t.1, X.1 = X.1, beta_t.1 = beta_t.1, Y.1 = Y.1, df = df))
}
# debug
# df = cbind(Y = test$Y.1$whole, V = test$V.1$whole, X = test$X.1$whole) %>% as.data.frame()
# 
# df = df %>% 
#   mutate(Y1 = lag(Y, n = 1, default = Y[1]),
#          Y2 = lag(Y, n = 2, default = Y[1]),
#          X1 = lag(X, n = 1, default = X[1]))

