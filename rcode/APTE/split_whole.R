splitAt <- function(x, pos) unname(split(x, cumsum(seq_along(x) %in% pos)))
# split the data 
split_whole <- function (data, whole_model_res) {
  cptm_stationary <- cpt.mean(whole_model_res$mean[-c(1:50)], penalty='Manual', pen.value = 10, method='PELT')
  cpts_stationary <- cpts(cptm_stationary)
  # cptm_stationary1 <- cpt.mean(whole_model_res$mean_lag[-c(1:50)], penalty='Manual', pen.value = 10, method='PELT')
  # period check
  i = 5
  while(length(cptm_stationary@param.est[["mean"]])!=5){
    i = i+1
    cptm_stationary <- cpt.mean(whole_model_res$mean[-c(1:50)], penalty='Manual', pen.value = i, method='PELT')
    cpts_stationary <- cpts(cptm_stationary)
  }
  # i = 5
  # while(length(cptm_stationary1@param.est[["mean"]])!=5){
  #   i = i+1
  #   cptm_stationary1 <- cpt.mean(whole_model_res$mean_lag[-c(1:50)], penalty='Manual', pen.value = i, method='PELT')
  #   cpts_stationary1 <- cpts(cptm_stationary1)
  # }
  
  # plot(cptm_stationary)
  # plot(cptm_stationary1)
  # cptm_stationary@param.est[["mean"]]
  cpts_stationary <- cpts_stationary + 50
  split_df = map(data$df, ~splitAt(.,cpts_stationary))

  return(split_df)
}