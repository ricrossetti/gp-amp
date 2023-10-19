### AMP for pure Gaussian or GOE matrices

## AMP recursions, matrix formulation
amp_gaus = function(f0, denoise, signal = NULL, snr = 1, maxt = 40, ...) {
  n = length(f0)
  f0 = f0 / norm(f0, '2') * sqrt(n)
  
  data = matrix(rnorm(n^2, sd = 1 /  sqrt(n)), n, n)
  if (!is.null(signal)) {
    signal = signal / norm(signal, '2')
    data = snr * outer(signal, signal) + data
  }
  
  iter = matrix(NA, n, 0)
  est = matrix(f0, n, 1)
  
  for (t in 1:maxt) {
    x = data %*% est[,t]
    iter = cbind(iter, x)
    f = apply(iter, 1, denoise, ...)
    est = cbind(est, f / norm(f, '2') * sqrt(n))
  }
  
  return( list( iter = iter, est = est[,-ncol(est)] ) )
}

amp_goe = function(f0, denoise, deriv, signal = NULL, snr = 1, maxt = 40, ...) {
  n = length(f0)
  f0 = f0 / norm(f0, '2') * sqrt(n)
  
  data = matrix(rnorm(n^2, sd = 1 /  sqrt(2*n)), n, n)
  if (!is.null(signal)) {
    signal = signal / norm(signal, '2')
    data = snr * outer(signal, signal) + data + t(data)
  }
  
  iter = matrix(NA, n, 0)
  est = matrix(f0, n, 1)
  
  for (t in 1:maxt) {
    x = data %*% est[,t]
    if (t > 1) {
      x = x - as.matrix(est[,-t]) %*% b
    }
    iter = cbind(iter, x)
    f = apply(iter, 1, denoise, ...)
    b = matrix(apply(iter, 1, deriv, ...), ncol = t, byrow = TRUE)
    b = colSums(b) / (sqrt(n) * norm(f, '2'))
    est = cbind(est, f / norm(f, '2') * sqrt(n))
  }
  
  return( list( iter = iter, est = est[,-ncol(est)] ) )
}
