# func.1: function to generate y for unconfounded situation

newyy <- function(period, V, mu, mu_2, x, beta) {
  
  Y.new = rnorm(m[period], mean = 0, sd =sqrt((1-phi_2[period])*var_kt/(1+phi_2[period])/(1-phi_1[period]-phi_2[period])/(1+phi_1[period]-phi_2[period]))) # initialization
  
  for (j in 3:m[period]) {
    Y.new[j] = beta_0t + beta[j] * V[j] + mu[j] * x[j] + mu_2[j] * x[j-1] + phi_1[period] * Y.new[j-1] + phi_2[period]*Y.new[j-2] + rnorm(1, 0, sqrt(var_kt))
  }
  
  return(list(X.new = x, Y.new = Y.new, prob = p))
}