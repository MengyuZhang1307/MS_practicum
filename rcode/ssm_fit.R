# func.2: fit state space model return kalma state and smoothed state and other model results

ssm_fit <- function(model, n) {
  fit <- fitSSM(model, inits = rep(0,n), method = "BFGS") # MLE of variance
  model <- fit$model
  out <- KFS(model)
  mean_filter = out$a[,3]
  mean_smoothed = out$alphahat[,3]
  mean_filter_lag = out$a[,4]
  mean_smoothed_lag = out$alphahat[,4]
  
  return(list(mean = mean_filter, mean_smoothed = mean_smoothed, mean_lag=mean_filter_lag, mean_smoothed_lag = mean_smoothed_lag, fitres = out))
}