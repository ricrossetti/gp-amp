### PERFORMANCE CHECKS ###
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
iter_demean = function(signal, iter, est) {
  mu = apply(est, 2, crossprod, y = signal)
  res = iter - tcrossprod(iter, mu)
  return(res)
}