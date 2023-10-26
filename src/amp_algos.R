### AMP for pure Gaussian or GOE matrices

## AMP recursions, matrix formulation
amp_gaus = function(f0, denoise, signal = NULL, snr = 1, maxt = 40, ...) {
  n = length(f0)
  f0 = f0 / norm(f0, '2')
  
  data = matrix(rnorm(n^2), n, n)
  if (!is.null(signal)) {
    signorm = norm(signal, '2')
    data = snr / signorm * outer(signal, signal) + data
  }
  
  est = matrix(f0, n, 1)
  iter = data %*% f0
  
  for (t in 1:(maxt-1)) {
    f = apply(iter, 1, denoise, ...)
    normf = norm(f, '2')
    x = data %*% f / normf
    est = cbind(est, f / normf)
    iter = cbind(iter, x)
  }
  
  Sigma = crossprod(est)
  
  if (is.null(signal)) {
    return( list( iter = iter, est = est, Sigma = Sigma ) )
  } else {
    overlap = c(crossprod(est, signal)) / signorm
    xcenter = iter - tcrossprod(signal, overlap * snr)
    return( list( iter = iter, est = est, Sigma = Sigma, overlap = overlap, xcenter = xcenter ) )
  }
}

amp_goe = function(f0, denoise, deriv, signal = NULL, snr = 1, maxt = 40, ...) {
  n = length(f0)
  f0 = f0 / norm(f0, '2')
  data = matrix(rnorm(n^2, sd = 1 /  sqrt(2)), n, n)
  if (!is.null(signal)) {
    signorm = norm(signal, '2')
    data = snr / signorm * outer(signal, signal) + data + t(data)
  }
  
  est = matrix(f0, n, 1)
  iter = data %*% f0
  
  for (t in 1:(maxt-1)) {
    f = apply(iter, 1, denoise, ...)
    normf = norm(f, '2')
    b = matrix(apply(iter, 1, deriv, ...), ncol = t, byrow = TRUE)
    b = colSums(b) / normf
    x = data %*% f / normf - est %*% b
    est = cbind(est, f / normf)
    iter = cbind(iter, x)
    # x = data %*% est[,t]
    # if (t > 1) {
    #   x = x - as.matrix(est[,-t]) %*% b
    # }
    # iter = cbind(iter, x)
    # f = apply(iter, 1, denoise, ...)
    # b = matrix(apply(iter, 1, deriv, ...), ncol = t, byrow = TRUE)
    # b = colSums(b) / (sqrt(n) * norm(f, '2'))
    # est = cbind(est, f / norm(f, '2') * sqrt(n))
  }
  
  Sigma = crossprod(est)
  
  if (is.null(signal)) {
    return( list( iter = iter, est = est, Sigma = Sigma ) )
  } else {
    overlap = c(crossprod(est, signal)) / signorm
    xcenter = iter - tcrossprod(signal, overlap * snr)
    return( list( iter = iter, est = est, Sigma = Sigma, overlap = overlap, xcenter = xcenter ) )
  }
}
