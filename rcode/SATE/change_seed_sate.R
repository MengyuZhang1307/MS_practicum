source("rcode/ts_gene_nocon.R")
source("rcode/APTE/fit_whole_getres.R")
source("rcode/APTE/split_whole.R")
source("rcode/APTE/fit_sep_getres.R")
source("rcode/SATE/get_tau.R")

change_seed_sate <- function(Seed, mod){
  # Seed = 35
  # mod = "large"
  # generate data
  data = ts_gene_nocon(seed = Seed)
  
  # fit whole times series
  whole_res = fit_whole_getres(data, variance = TRUE, model = mod) # if we change this to FALSE, it would not be able to detect different periods.
  
  # split data
  split_data = split_whole(data, whole_res)
  
  # fit seperate model
  sep_res = fit_sep_getres(split_data, variance = FALSE, model = mod)
  
  # get SATE
  m_split = map(split_data$Y, ~length(.x)) %>% unlist # number of time points in each period
  Z = split_data
  Z$beta0 = map(m_split, ~rep(1,.x))
  Z_tilda = Z
  Z_tilda$X = map(split_data$X, ~(1-.x))

  SATE <- list()
  for (i  in 1:length(m_split)) {
    SATE[[paste0("period",i)]] = get_tau(split_data,m_split, Z, Z_tilda, sep_res[[i]], i, 100)
  }
  
  mean_sate = map(SATE, ~mean(.x)) %>% unlist
  var_sate = map(SATE, ~var(.x)) %>% unlist
  
  return(list(mean_sate,var_sate))
}