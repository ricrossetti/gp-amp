### AMP for pure Gaussian or GOE matrices
# source('./helpers.R')

## AMP recursions, matrix formulation
amp_gaus = function(x, denoise, signal = NULL, maxt = 40, ...) {
  n = length(x)
  
  data = rnorm(n^2, sd = 1 /  sqrt(n))
  if (!is.null(signal)) {
    data = outer(signal, signal) + data
  }
  
  iter = matrix(x, n, 1)
  est = matrix(NA, n, 0)
  
  for (t in 1:maxt) {
    f = apply(iter, 1, denoise, ...)
    est = cbind(est, f / norm(f, '2') )
    x = data %*% est[,t]
    iter = cbind(iter, x)
  }
  
  return( list( iter = iter[,-1], est = est ) )
}

amp_goe = function(x, denoise, deriv, signal = NULL, maxt = 40, ...) {
  n = length(x)
  
  z = rnorm(n^2, sd = 1 /  sqrt(2*n))
  if (!is.null(signal)) {
    data = outer(signal, signal) + z + t(z)
  }
  
  iter = matrix(x, n, 1)
  est = matrix(0, n, 0)
  
  for (t in 1:maxt) {
    f = apply(iter, 1, denoise, ...)
    est = cbind(est, f / norm(f, '2') )
    # Prepare correction term
    if (t != 1) {
      b = apply(1, deriv, 
    }
  }
}


## AMP recursions, GP formulation
# Gaussian noise
gp_amp_gaus = function(x, denoise, signal = NULL, memory = FALSE, scale = TRUE, maxt = 40, ...) {
  n = length(x)
  if (is.null(signal)) {
    data = rnorm(n^2, sd = 1 / sqrt(n))
  } else {
    if (scale) { signal = signal / norm(signal, '2') }
    data = kronecker(signal, signal) + rnorm(n^2, sd = 1 / sqrt(n))
  }
  
  t = 0
  iter = matrix(x, n, 1)
  est = matrix(NA, n, 0)
  
  while( t < maxt ) {
    t = t + 1
    f = if (memory == FALSE) {
      mapply(denoise, x, ...) 
    } else {
      mapply(denoise, x, memory = lapply( seq_len( nrow(iter) ), function(i) iter[i, -ncol(iter)] ), ...)
    }
    f = f / norm(f, '2')
    x = vkmult(f, data)
    est = cbind(est, f)
    iter = cbind(iter, x)
  }
  
  return(list(iter = iter[,-1], est = est))
}

# GOE noise functions
gp_amp_goe = function(x, denoise, deriv, signal = NULL, memory = FALSE, scale = TRUE, maxt = 40, ...) {
  n = length(x)
  if (is.null(signal)) {
    data = c(rgoe(n))
  } else {
    data = kronecker(signal, signal) / n + c(rgoe(n))
  }

  t = 0
  iter = matrix(x, n, 1)
  est = matrix(NA, n, 0)
  onsager = 0
  
  while( t < maxt ) {
    t = t + 1
    f = if (memory == FALSE) {
      mapply(denoise, x, ...)
    } else {
      mapply(denoise, x, memory = lapply( seq_len( nrow(iter) ), function(i) iter[i, -ncol(iter)] ), ...)
    }
    if (scale) {
      normf = f / norm(f, '2')
      est = cbind(est, normf)
      x = vkmult(normf, data) - onsager
    } else {
      est = cbind(est, f)
      x = vkmult(f, data) - onsager
    }
    if (memory == FALSE) {
      ovect = mapply(deriv, x, ...)
      if (!scale) {
        ovect = ovect / norm(f, '2') * (1 - normf^2)
        onsager = sum(ovect) / n * normf
      }
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
  }
  
  return(list(iter = iter[,-1], est = est))
}

# Memory, no scaling
amp_goe_noscale_mem = function(x, denoise, deriv, signal = NULL, maxt = 40, ...) {
  n = length(x)
  if (is.null(signal)) {
    data = c(rgoe(n))
  } else {
    data = kronecker(signal, signal) / n + c(rgoe(n))
  }
  
  t = 0
  iter = matrix(x, n, 1)
  est = matrix(NA, n, 0)
  onsager = 0
  
  while( t < maxt ) {
    t = t + 1
    f = mapply(denoise, x, memory = lapply( seq_len( nrow(iter) ), function(i) iter[i, -ncol(iter)] ), ...)
    est = cbind(est, f)
    x = vkmult(f, data) - onsager
    ovect = mapply(deriv, x, memory = lapply( seq_len( nrow(iter) ), function(i) iter[i, -ncol(iter)] ), ...)
    if (is.vector(ovect)) {
      onsager = as.vector(est * sum(ovect) / n)
    } else {
      onsager = as.vector(est %*% rowSums(ovect) / n)
    }
    iter = cbind(iter, x)
  }
  
  return(list(iter = iter[,-1], est = est))
}

# No memory, no scaling
amp_goe_noscale_nomem = function(x, denoise, deriv, signal = NULL, maxt = 40, ...) {
  n = length(x)
  if (is.null(signal)) {
    data = c(rgoe(n))
  } else {
    data = kronecker(signal, signal) / n + c(rgoe(n))
  }
  
  t = 0
  iter = matrix(x, n, 1)
  est = matrix(NA, n, 0)
  onsager = 0
  
  while( t < maxt ) {
    t = t + 1
    f = mapply(denoise, x, ...)
    est = cbind(est, f)
    x = vkmult(f, data) - onsager
    ovect = mapply(deriv, x, ...)
    onsager = sum(ovect) / n * f
    iter = cbind(iter, x)
  }
  
  return(list(iter = iter[,-1], est = est))
}

# No memory, scaled
amp_goe_scale_nomem = function(x, denoise, deriv, signal = NULL, maxt = 40, ...) {
  n = length(x)
  if (is.null(signal)) {
    data = c(rgoe(n))
  } else { 
    signal = signal / norm(signal, '2')
    data = kronecker(signal, signal) + c(rgoe(n))
  }
  
  t = 0
  iter = matrix(x, n, 1)
  est = matrix(NA, n, 0)
  onsager = 0
  
  while( t < maxt ) {
    t = t + 1
    f = mapply(denoise, x, ...)
    normf = f / norm(f, '2')
    if (t > 1) {
      onsager = sum( ovect / norm(f, '2') * (1 - normf^2) ) / n * est[, ncol(est)]
    }
    est = cbind(est, normf)
    x = vkmult(f, data) - onsager
    ovect = mapply(deriv, x, ...)
    iter = cbind(iter, x)
  }
  
  return(list(iter = iter[,-1], est = est))
}

# Memory, scale
amp_goe_scale_mem = function(x, denoise, deriv, signal = NULL, maxt = 40, ...) {
  n = length(x)
  if (is.null(signal)) {
    data = c(rgoe(n))
  } else { 
    signal = signal / norm(signal, '2')
    data = kronecker(signal, signal) + c(rgoe(n))
  }
  
  t = 0
  iter = matrix(x, n, 1)
  est = matrix(NA, n, 0)
  onsager = 0
  
  while( t < maxt ) {
    t = t + 1
    f = mapply(denoise, x, memory = lapply( seq_len( nrow(iter) ), function(i) iter[i, -ncol(iter)] ), ...)
    normf = f / norm(f, '2')
    if (t > 1) {
      if (is.vector(ovect)) {
        ovect_scale = ovect / norm(f, '2') * (1 - normf^2)
        onsager = as.vector(sum(ovect_scale) * est)
      } else {
        normf2_replicate = matrix( rep(normf^2, nrow(ovect)), ncol = n, byrow = TRUE )
        ovect_scale = ovect / norm(f, '2') * (1 - normf2_replicate )
        onsager = as.vector(est %*% rowSums(ovect_scale))
      }
    }
    est = cbind(est, normf)
    x = vkmult(f, data) - onsager
    ovect = mapply(deriv, x, memory = lapply( seq_len( nrow(iter) ), function(i) iter[i, -ncol(iter)] ), ...)
    fold = normf
    iter = cbind(iter, x)
  }
  
  return(list(iter = iter[,-1], est = est))
}

## TESTING: COMMENT OUT WHEN DONE
# Parameter settings
# n = 1000
# snr = 2
# signal = sqrt(snr) * rbinom(n, 1, prob = .5)
# init = 0.1 * signal + rnorm(n)
# 
# # Test GOE + spike, denoisers with memory
# test_nospike = gp_amp_goe(init, dlmr, dlmr_grad, memory = TRUE)
# test_spike = gp_amp_goe(init, dlmr, dlmr_grad, signal = signal, memory = TRUE)
