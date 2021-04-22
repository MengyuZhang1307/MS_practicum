
# different seed for simulation inference
source("rcode/ts_gene_nocon.R")
source("rcode/APTE/fit_whole_getres.R")
source("rcode/APTE/split_whole.R")
source("rcode/APTE/fit_sep_getres.R")

change_seed_nocon <- function(Seed, mod) {
  
  # generate data
  data = ts_gene_nocon(seed = Seed) 
  
  # fit whole times series
  whole_res = fit_whole_getres(data, variance = TRUE, model = mod) # if we change this to FALSE, it would not be able to detect different periods.
  
  # split data
  split_data = split_whole(data, whole_res)
  
  # fit seperate model
  sep_res = fit_sep_getres(split_data, variance = FALSE, model = mod)
  
  return (sep_res)
  
}
