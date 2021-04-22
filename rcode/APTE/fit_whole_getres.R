# define model, fit whole model, and get results
source("rcode/ssm_fit.R")

fit_whole_getres <- function(data, variance, model){
  var_est = ifelse (variance == TRUE, NA, 0)
  if (model == "large"){
  model.nocon <- SSModel(data$df$Y ~ SSMregression(~ data$df$V + data$df$X + data$df$X1 + data$df$Y1 + data$df$Y2, type = "distinct", Q = diag(c(rep(var_est,3),0,0))), H = NA)
  res.nocon = ssm_fit(model.nocon,6)
  }
  if (model == "small") {
  model.nocon <- SSModel(data$df$Y ~ SSMregression(~ data$df$V + data$df$X, type = "distinct", Q = diag(rep(var_est,2))), H = NA)
    res.nocon = ssm_fit(model.nocon,3)
  }
  return(res.nocon)
}
