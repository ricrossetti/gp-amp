### AMP for pure Gaussian or GOE matrices
source('./helpers.R')

## AMP recursions
## Matrix-version implementations
# Fully Gaussian matrix
amp_gaus = function(x, denoise, theta, spike = NULL, tol = 1e-4) {
  # Add spike if present
  if (is.null(spike)) {
    n = length(x)
    data = rgaus_mat(n)
  } else {
    n = length(x)
    data = spike + rgaus_mat(n)
  }
  
  improve = Inf
  t = 0
  maxt = 100
  iter = matrix(NA, n, 0)
  
  while( (improve > tol) & (t < maxt) ) {
    t = t + 1
    f = sapply(x, denoise, theta = theta)
    x = as.vector(data %*% f)
    iter = cbind(iter, x)
    if (t > 1) {
      improve = abs( var(iter[,ncol(iter)]) - var(iter[,(ncol(iter) - 1)]) )
    }
  }
  
  return(iter)
}

# GOE matrix
amp_goe = function(x, denoise, theta, spike = NULL, tol = 1e-4) {
  # Add spike if present
  if (is.null(spike)) {
    n = length(x)
    data = rgoe(n)
  } else {
    n = length(x)
    data = spike + rgoe(n)
  }
  
  deriv = Deriv(denoise, formalArgs(denoise)[1])
  
  improve = Inf
  t = 0
  maxt = 100
  iter = matrix(NA, n, 0)
  onsager = 0
  
  while ( (improve > tol) & (t < maxt) ) {
    t = t+1
    f = sapply(x, denoise, theta = theta)
    x = as.vector(data %*% f) - onsager
    onsager = sum( sapply(x, deriv, theta = theta) ) / n * f
    iter = cbind(iter, x)
    if (t > 1) {
      improve = abs( var(iter[,ncol(iter)]) - var(iter[,(ncol(iter) - 1)]) )
    }
  }
  
  return(iter)
}

## Diagnostics of interest
# Deviations from normality of the iterates
deviations = function(iterates) {
  evalpoints = ppoints(nrow(iterates))
  quants = qnorm(evalpoints)
  norm_iters = apply(iterates, 2, function(x) sort(x) / sd(x) )
  
  dev = apply(norm_iters, 2, function(x,q) x - q, q = quants)
  return(dev)
}

# Normalized overlaps
overlap = function(x, truth) {
  abs(cor(x,truth))
}

## TESTING: COMMENT OUT WHEN DONE
# Parameter settings
n = 1000
signal = rspike(n)
snr = 2
spike = sqrt(snr) * signal[[2]]
init = 0.1 * signal[[1]] + rnorm(n)

# Pure Gaussian noise 
gaus_iter = amp_gaus(init, dtanh, 1)
goe_fact = amp_goe(init, dtanh, 1)

# Gaussian signal + noise matrix factorization
gaus_fact = amp_gaus(init, dtanh, 1, spike = spike)
goe_fact = amp_goe(init, dtanh, 1, spike = spike)

