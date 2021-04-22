get_tau <- function(split_data, m_split, Z, Z_tilda, model,period,B){
  z_tilda = cbind(beta0 = Z_tilda$beta0[[period]], 
                  V = Z_tilda$V[[period]], 
                  X = Z_tilda$X[[period]], 
                  X1 = Z_tilda$X1[[period]], 
                  Y1 = Z_tilda$Y1[[period]], 
                  Y2 = Z_tilda$Y2[[period]])
  
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
  return(APTE_sate)
}