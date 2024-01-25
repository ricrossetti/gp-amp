## Expected overlap functions
rad_ov = function(snr) {
  integrand = function(x) {
    exp( log( tanh(x)^2 ) + dnorm(x, 
                                  mean = snr, 
                                  sd = sqrt(snr), 
                                  log = TRUE) ) 
  }
  integrate(integrand, -Inf, Inf)$value
}

gaus_ov = function(snr) {
  snr / (1 + snr)
}

sparse_rad_ov = function(snr, p) {
  snr = snr / p
  f = function(x) { sqrt(p) / tanh(sqrt(snr)*x) }
  g = function(x) {
    a = log( 2*(1-p) / sqrt(p) ) + snr / 2  
    b = log( exp(sqrt(snr)*abs(x)) - exp(-sqrt(snr)*abs(x)) )
    res = exp(a - b)
    return( sign(x) * res )
  }
  integrand0 = function(x) {
    exp( -log( (f(x) + g(x))^2 ) + dnorm(x, log = TRUE) )
  }
  integrand1 = function(x) {
    exp( -log( (f(x) + g(x))^2 ) + 
           dnorm(x, mean = sqrt(snr), log = TRUE) )
  }
  res = (1-p)*integrate(integrand0, -Inf, Inf)$value +
    p*integrate(integrand1, -Inf, Inf)$value
  return(res)
}

# SE mappings
se_erasures = function(initsnr, 
                       ov_fun,
                       lambda, 
                       beta = 1, 
                       tol = -6,
                       iter = NULL,
                       ...) {
  snr = c(initsnr, lambda^2 * ov_fun(initsnr, ...))
  improve = Inf 
  if (is.null(iter)) {
    while (improve > tol) {
      snr = c(snr,
              lambda^2 * beta * ov_fun(snr[length(snr)], ...) + 
                (1-beta) * snr[length(snr)] )
      improve = log10( abs(snr[length(snr)] - snr[length(snr)-1]) ) -
        log10( snr[length(snr)-1] )
    }
  } else {
    for (i in 1:(iter-1)) {
      snr = c(snr,
              lambda^2 * beta * ov_fun(snr[length(snr)], ...) + 
                (1-beta) * snr[length(snr)] )
    }
  }
  return(sqrt(snr[-1] / lambda^2))
}

se_robin = function(initsnr, 
                    ov_fun,
                    lambda, 
                    beta = 1, 
                    tol = -6,
                    iter = NULL,
                    ...) {
  snr = c(initsnr, lambda^2 * ov_fun(initsnr, ...))
  improve = Inf 
  cycle_len = 1/beta
  if (cycle_len %% 1 != 0) { 
    return("Problem size is not a multiple of subset size") 
    }
  if (is.null(iter)) {
    while (improve > tol) {
      delta = ov_fun(snr[length(snr)], ...) -
        ov_fun(snr[max(1,length(snr)-cycle_len)], ...) 
      snr = c(snr, snr[length(snr)] + lambda^2 * beta * delta )
      improve = log10( abs(snr[length(snr)] - snr[length(snr)-1]) ) -
        log10( snr[length(snr)-1] )
    } 
  } else {
    for (i in 1:(iter-1)) {
      delta = ov_fun(snr[length(snr)], ...) -
        ov_fun(snr[max(1,length(snr)-cycle_len)], ...) 
      snr = c(snr, snr[length(snr)] + lambda^2 * beta * delta )
    }
  }
  return(sqrt(snr[-1] / lambda^2))
}

se_erasures_pi = function(initsnr,
                          lambda,
                          beta = 1,
                          tol = -6,
                          iter = NULL) {
  snr = initsnr
  improve = Inf
  if (is.null(iter)) {
    while(improve > tol) {
      log_num = 2* (log(lambda) + log(beta) + 
                      log( sum( sqrt(snr) * (1-beta)^rev((1:length(snr)) - 1) ) ) )
      log_den = log( beta * sum( snr * (1-beta)^rev((1:length(snr)) - 1) ) +
                       (1-(1-beta)^(length(snr))) +
                       (1-beta)^(length(snr)) * initsnr)
      snr = c(snr,
              exp( log_num - log_den ))
      improve = log10( abs(snr[length(snr)] - snr[length(snr)-1]) ) -
        log10( snr[length(snr)-1] )
    }
  } else {
    for (i in 1:(iter-1)) {
      log_num = 2* (log(lambda) + log(beta) + 
                      log( sum( sqrt(snr) * (1-beta)^rev((1:length(snr)) - 1) ) ) )
      log_den = log( beta * sum( snr * (1-beta)^rev((1:length(snr)) - 1) ) +
                       (1-(1-beta)^(length(snr))) +
                       (1-beta)^(length(snr)) * initsnr)
      snr = c(snr,
              exp( log_num - log_den ))
    }
  }
  return(sqrt(snr / lambda^2))
}

se_erasures_pi_robin = function(initsnr,
                                lambda,
                                beta = 1,
                                tol = -6,
                                iter = NULL) {
  snr = initsnr
  cycle_len = 1/beta
  if (cycle_len %% 1 != 0) { 
    return("Problem size is not a multiple of subset size") 
  }
  improve = Inf
  remainder = 1
  if (is.null(iter)) {
    while(improve > tol) {
      remainder = max(0, remainder - beta)
      log_num = 2* (log(lambda) +
                      log( sum( rev(sqrt(snr))[1:min(cycle_len,length(snr))]) +
                             remainder * initsnr ) )
      log_den = log( sum( rev(snr)[1:min(cycle_len,length(snr))]) + 
                       remainder * initsnr + (1 - remainder) * cycle_len ) + 
        log( cycle_len)
      snr = c(snr,
              exp( log_num - log_den ))
      improve = log10( abs(snr[length(snr)] - snr[length(snr)-1]) ) -
        log10( snr[length(snr)-1] )
    }
  } else {
    for (i in 1:(iter-1)) {
      remainder = max(0, remainder - beta)
      log_num = 2* (log(lambda) +
                      log( sum( rev(sqrt(snr))[1:min(cycle_len,length(snr))]) +
                             remainder * initsnr ) )
      log_den = log( sum( rev(snr)[1:min(cycle_len,length(snr))]) + 
                       remainder * initsnr + (1 - remainder) * cycle_len ) + 
        log( cycle_len)
      snr = c(snr,
              exp( log_num - log_den ))
    }
  }
  return(sqrt(snr / lambda^2))
}
