### AMP for pure Gaussian or GOE matrices

# Gaussian AMP without memory
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
    f = denoise(iter[,t])
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

# Gaussian AMP with memory (CONSIDER USING GIP INSTEAD)
amp_gaus_memory = function(f0, denoise, signal = NULL, snr = 1, maxt = 40, ...) {
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

# GOE AMP without memory
amp_goe = function(f0, denoise, deriv, signal = NULL, snr = 1, maxt = 40, seed = NULL, ...) {
  n = length(f0)
  f0 = f0 / norm(f0, '2')
  if(!is.null(seed)) { set.seed = seed }
  data = matrix(rnorm(n^2, sd = 1 /  sqrt(2)), n, n)
  if (!is.null(signal)) {
    signorm = norm(signal, '2')
    data = snr / signorm * outer(signal, signal) + data + t(data)
  }
  
  est = matrix(f0, n, 1)
  iter = data %*% f0
  
  for (t in 1:(maxt-1)) {
    f = denoise(iter[,t])
    normf = norm(f, '2')
    b = deriv(iter[,t])
    b = sum(b) / normf
    x = data %*% f / normf - est[,t] * b
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

# GOE AMP with memory
amp_goe_memory = function(f0, denoise, deriv, signal = NULL, snr = 1, maxt = 40, ...) {
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

### AMP recursions for masked matrices
## Row-wise erasures
amp_goe_erasures = function(f0,
                            denoise,
                            deriv,
                            sub_size = ceiling(length(f0) / 20),
                            signal = NULL,
                            snr = 1,
                            maxt = 100,
                            seed = NULL,
                            ...) {
  n = length(f0)
  f0 = f0 / norm(f0, '2')
  if(!is.null(seed)) { set.seed = seed }
  data = matrix(rnorm(n^2, sd = 1 /  sqrt(2)), n, n)
  if (!is.null(signal)) {
    signorm = norm(signal, '2')
    data = snr / signorm * outer(signal, signal) + data + t(data)
  }
  
  est = matrix(f0, n, 1)
  S = sort(sample(1:n, sub_size))
  iter = matrix(0, n, 1)
  iter[S,1] = data[S,] %*% f0
  iter[-S,1] = rnorm(n-sub_size)
  if (!is.null(signal)) {
    overlap = as.numeric(crossprod(f0, signal) / signorm)
    m = matrix(0, n, 1)
    m[S,1] = snr * overlap * signal[S]
  }
  
  for (t in 1:(maxt-1)) {
    if (t == 1) {
      f = denoise(iter[,t])
      normf = norm(f, '2')
      b = deriv(iter[,t])
      b = sum(b) / normf
      S = sort(sample(1:n, sub_size))
      x = data[S,] %*% f / normf - b * est[S,t]
      est = cbind(est, f / normf)
      iter = cbind(iter, 0)
      iter[S,(t+1)] = x; iter[-S,(t+1)] = iter[-S,t]
      if (!is.null(signal)) {
        overlap = c(overlap, as.numeric(crossprod(est[,t+1], signal) / signorm))
        m = cbind(m, 0)
        m[S, t+1] = snr * overlap[t+1] * signal[S]; m[-S, t+1] = m[-S, t] 
      }
    } else {
      f = denoise(iter[,t])
      normf = norm(f, '2')
      b = deriv(iter[,t]) 
      b1 = sum(b[S]) / normf
      b2 = sum(b[-S]) / normf
      S = sort(sample(1:n, sub_size))
      x = data[S,] %*% f / normf - b1 * est[S,t] - b2 * est[S,(t-1)]
      est = cbind(est, f / normf)
      iter = cbind(iter, 0)
      iter[S,(t+1)] = x; iter[-S,(t+1)] = iter[-S,t]
      if (!is.null(signal)) {
        overlap = c(overlap, as.numeric(crossprod(est[,t+1], signal) / signorm))
        m = cbind(m, 0)
        m[S, t+1] = snr * overlap[t+1] * signal[S]; m[-S, t+1] = m[-S, t] 
      }
    }
  }
    if (is.null(signal)) {
      return( list( iter = iter,
                    est = est) )
    } else {
      return( list( iter = iter,
                    est = est,
                    signal = signal,
                    snr = snr,
                    overlap = overlap,
                    mean = m) )
    }
}



