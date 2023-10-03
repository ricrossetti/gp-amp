### AMP for pure Gaussian or GOE matrices
source('./helpers.R')

## AMP recursions
# Matrix formulation implementation
# Gaussian noise
gp_amp_gaus = function(x, denoise, signal = NULL, memory = FALSE, tol = 1e-4, ...) {
  n = length(x)
  if (is.null(signal)) {
    data = rnorm(n^2, sd = 1 / sqrt(n))
  } else {
    data = kronecker(signal, signal) / n + rnorm(n^2, sd = 1 / sqrt(n))
  }
  
  improve = Inf
  t = 0
  maxt = 100
  iter = matrix(x, n, 1)
  est = matrix(NA, n, 0)
  
  while( (improve > tol) & (t < maxt) ) {
    t = t + 1
    f = if (memory == FALSE) {
      mapply(denoise, x, ...)
    } else {
      mapply(denoise, x, memory = lapply( seq_len( nrow(iter) ), function(i) iter[i, -ncol(iter)] ), ...)
    }
    x = vkmult(f, data)
    est = cbind(est, f)
    iter = cbind(iter, x)
    if (t > 1) {
      improve = abs( var(iter[,ncol(est)]) - var(iter[,(ncol(est) - 1)]) )
    }
  }
  
  return(list(iter = iter[,-1], est = est))
}

# GOE noise
gp_amp_goe = function(x, denoise, deriv, signal = NULL, memory = FALSE, tol = 1e-4, ...) {
  n = length(x)
  if (is.null(signal)) {
    data = c(rgoe(n))
  } else {
    data = kronecker(signal, signal) / n + c(rgoe(n))
  }
  
  improve = Inf
  t = 0
  maxt = 100
  iter = matrix(x, n, 1)
  est = matrix(NA, n, 0)
  onsager = 0
  
  while( (improve > tol) & (t < maxt) ) {
    t = t + 1
    f = if (memory == FALSE) {
      mapply(denoise, x, ...)
    } else {
      mapply(denoise, x, memory = lapply( seq_len( nrow(iter) ), function(i) iter[i, -ncol(iter)] ), ...)
    }
    est = cbind(est, f)
    x = vkmult(f, data) - onsager
    if (memory == FALSE) {
      ovect = mapply(deriv, x, ...)
      onsager = sum(ovect) / n * f
    } else {
      ovect = mapply(deriv, x, memory = lapply( seq_len( nrow(iter) ), function(i) iter[i, -ncol(iter)] ), ...)
      if (is.vector(ovect)) {
        onsager = as.vector(est * sum(ovect) / n)
      } else {
        onsager = as.vector(est %*% rowSums(ovect) / n)
      }
    }
    iter = cbind(iter, x)
    if (t > 1) {
      improve = abs( var(iter[,ncol(est)]) - var(iter[,(ncol(est) - 1)]) ) 
    }
  }
  
  return(list(iter = iter[,-1], est = est))
}

## TESTING: COMMENT OUT WHEN DONE
# Parameter settings
n = 1000
snr = 3
signal = sqrt(snr) * rbinom(n, 1, prob = .8)
init = 0.1 * signal + rnorm(n)

# Test GOE + spike
test_nospike = gp_amp_goe(init, dtanh, theta = 1)
test_spike = gp_amp_goe(init, dtanh, signal = signal, theta = 1)

### MATRIX VERSIONS OF RECURSIONS
# # Fully Gaussian matrix
# amp_gaus = function(x, denoise, theta, signal = NULL, tol = 1e-4) {
#   # Add signal if present
#   if (is.null(signal)) {
#     n = length(x)
#     data = rgaus_mat(n)
#   } else {
#     n = length(x)
#     data = tcrossprod(signal) / n + rgaus_mat(n)
#   }
#   
#   improve = Inf
#   t = 0
#   maxt = 100
#   iter = matrix(NA, n, 0)
#   
#   while( (improve > tol) & (t < maxt) ) {
#     t = t + 1
#     f = sapply(x, denoise, theta = theta)
#     x = as.vector(data %*% f)
#     iter = cbind(iter, x)
#     if (t > 1) {
#       improve = abs( var(iter[,ncol(iter)]) - var(iter[,(ncol(iter) - 1)]) )
#     }
#   }
#   
#   return(iter)
# }
# 
# # GOE matrix
# amp_goe = function(x, denoise, theta, signal = NULL, tol = 1e-4) {
#   # Add signal if present
#   n = length(x)
#   if (is.null(signal)) {
#     data = rgoe(n)
#   } else {
#     data = tcrossprod(signal) / n + rgaus_mat(n)
#   }
#   
#   deriv = Deriv(denoise, formalArgs(denoise)[1])
#   
#   improve = Inf
#   t = 0
#   maxt = 100
#   iter = matrix(NA, n, 0)
#   onsager = 0
#   
#   while ( (improve > tol) & (t < maxt) ) {
#     t = t+1
#     f = sapply(x, denoise, theta = theta)
#     x = as.vector(data %*% f) - onsager
#     onsager = sum( sapply(x, deriv, theta = theta) ) / n * f
#     iter = cbind(iter, x)
#     if (t > 1) {
#       improve = abs( var(iter[,ncol(iter)]) - var(iter[,(ncol(iter) - 1)]) )
#     }
#   }
#   
#   return(iter)
# }
