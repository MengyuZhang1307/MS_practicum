
# define model, fit seperate model, and get results
source("rcode/ssm_fit.R")

fit_sep_getres <- function (data, variance, model) {
  var_est = ifelse (variance == TRUE, NA, 0)
  
  if(model == "large"){
  model1.nocon <- SSModel(data$Y[[1]] ~ SSMregression(~ data$V[[1]] + data$X[[1]] + data$X1[[1]] + data$Y1[[1]] + data$Y2[[1]], type = "distinct", Q = diag(c(rep(var_est,3),0,0))), H = NA)
  
  model2.nocon <- SSModel(data$Y[[2]] ~ SSMregression(~ data$V[[2]] + data$X[[2]] + data$X1[[2]] + data$Y1[[2]] + data$Y2[[2]], type = "distinct", Q = diag(c(rep(var_est,3),0,0))), H = NA)
  
  model3.nocon <- SSModel(data$Y[[3]] ~ SSMregression(~ data$V[[3]] + data$X[[3]] + data$X1[[3]] + data$Y1[[3]] + data$Y2[[3]], type = "distinct", Q = diag(c(rep(var_est,3),0,0))), H = NA)
  
  model4.nocon <- SSModel(data$Y[[4]] ~ SSMregression(~ data$V[[4]] + data$X[[4]] + data$X1[[4]] + data$Y1[[4]] + data$Y2[[4]], type = "distinct", Q = diag(c(rep(var_est,3),0,0))), H = NA)
  
  model5.nocon <- SSModel(data$Y[[5]] ~ SSMregression(~ data$V[[5]] + data$X[[5]] + data$X1[[5]] + data$Y1[[5]] + data$Y2[[5]], type = "distinct", Q = diag(c(rep(var_est,3),0,0))), H = NA)
  
  res1.nocon = ssm_fit(model1.nocon,6)
  res2.nocon = ssm_fit(model2.nocon,6)
  res3.nocon = ssm_fit(model3.nocon,6)
  res4.nocon = ssm_fit(model4.nocon,6)
  res5.nocon = ssm_fit(model5.nocon,6)
  }
  
  if (model == "small") {
    model1.nocon <- SSModel(data$Y[[1]] ~ SSMregression(~ data$V[[1]] + data$X[[1]], type = "distinct", Q = diag(rep(var_est,2))), H = NA)
    
    model2.nocon <- SSModel(data$Y[[2]] ~ SSMregression(~ data$V[[2]] + data$X[[2]], type = "distinct", Q = diag(rep(var_est,2))), H = NA)
    
    model3.nocon <- SSModel(data$Y[[3]] ~ SSMregression(~ data$V[[3]] + data$X[[3]], type = "distinct", Q = diag(rep(var_est,2))), H = NA)
    
    model4.nocon <- SSModel(data$Y[[4]] ~ SSMregression(~ data$V[[4]] + data$X[[4]], type = "distinct", Q = diag(rep(var_est,2))), H = NA)
    
    model5.nocon <- SSModel(data$Y[[5]] ~ SSMregression(~ data$V[[5]] + data$X[[5]], type = "distinct", Q = diag(rep(var_est,2))), H = NA)
    res1.nocon = ssm_fit(model1.nocon,3)
    res2.nocon = ssm_fit(model2.nocon,3)
    res3.nocon = ssm_fit(model3.nocon,3)
    res4.nocon = ssm_fit(model4.nocon,3)
    res5.nocon = ssm_fit(model5.nocon,3)
  }
  
  
  
  return(list(res1.nocon, res2.nocon, res3.nocon, res4.nocon, res5.nocon))
}