### PERFORMANCE CHECKS ###
## Various state evolution quantities
# Empirical covariance of iterates
emp_cov = function(iter) {
  crossprod(iter) / n
}
# SE estimator from estimates
emp_cov = function(est) {
  crossprod(est) / n
}

# Exact SE (Monte Carlo approximation)
se_mc = function(denoise, x, snr, signal = NULL, size = 40, memory = FALSE, ...) {
  n = length(x)
  if (!is.null(signal)) {
    m = vector('numeric', size)
  }
  sig = matrix(NA, size, size)
  
  gaussians = matrix(NA, n, 0)
  est = matrix(NA, n, 0)
  
  for (s in 1:size) {
    f = if (memory == FALSE) {
      mapply(denoise, x, ...)
    } else {
      mapply(denoise, x, memory = lapply( seq_len( nrow(gaussian) ), function(i) gaussians[i, ] ), ...)
    }
    
    if (s == 1) {
      sig[s,s] = crossprod(f) / n
      g = rnorm(n) * sqrt(sig[s,s])
      if (!is.null(signal)) {
        m[s] = snr * crossprod(signal, f) / n
        x = m[s] * signal + g
      } else {
        x = g
      }
    } else {
      sig[s,s] = crossprod(f) / n
      sig[s, 1:(s-1)] = crossprod(f, est) / n; sig[1:(s-1), s] = sig[s, 1:(s-1)]
      pseudoinv = solve(sig[1:(s-1), 1:(s-1)]) %*% sig[s, 1:(s-1)]
      condvar = sig[s,s] - crossprod(sig[s, 1:(s-1)] %*% pseudoinv)
      g = rnorm(n) * c(sqrt(condvar)) + gaussians %*% pseudoinv
      if (!is.null(signal)) {
        m[s] = snr * crossprod(signal, f) / n
        x = m[s] * signal + g
      } else {
        x = g
      }
    }
    
    gaussians = cbind(gaussians, g)
    est = cbind(est, f)
    
  }
  if (!is.null(signal)) {
    return(list(mu = m, sig = sig))
  } else {
    return(sig)
  }
}

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

# Discrepancy between SE and iterate covariance
cov_discrepancy = function(iter, est) {
  res = abs(cov(est) - cov(iter))
  return(res)
}

# Demean iterations
iter_demean = function(signal, iter, est, snr) {
  mu = sqrt(snr) * apply(est, 2, crossprod, y = signal) / n
  res = iter - tcrossprod(signal, mu)
  return(res)
}

