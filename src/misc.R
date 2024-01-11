## SNR estimator for erasures
ov_est = function(x, d) {
  sqrt(norm(x, '2')^2 - length(x)) / sqrt(d)
}

## Numerically compute SE overlap
rad_ov = function(rho, scale) {
  integrand = function(x) {
    exp( log( tanh(rho*x)^2 ) + dnorm(x, mean = rho, log = TRUE) ) 
  }
  scale * integrate( integrand, lower=-Inf, upper=Inf )
}