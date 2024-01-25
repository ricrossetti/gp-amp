### AMP for GOE matrices
## Standard AMP
amp_goe = function(f0, 
                   denoise, 
                   deriv,
                   data = NULL,
                   signal = NULL, 
                   snr = 1, 
                   maxt = 40,
                   ...) {
  n = length(f0)
  f0 = f0 / norm(f0, '2')
  if(!is.null(signal)) {
    signorm = norm(signal, '2')
  }
  
  if (is.null(data)) {
    data = matrix(rnorm(n^2, sd = 1 /  sqrt(2)), n, n)
    if (!is.null(signal)) {
      data = snr / signorm * outer(signal, signal) + data + t(data)
    } else {
      data = data + t(data)
    }
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
    return( list( iter = iter, 
                  est = est, 
                  Sigma = Sigma, 
                  overlap = overlap, 
                  xcenter = xcenter ) )
  }
}

## AMP with memory
amp_goe_memory = function(f0, 
                          denoise, 
                          deriv, 
                          data = NULL,
                          signal = NULL, 
                          snr = 1, 
                          maxt = 40, 
                          ...) {
  n = length(f0)
  f0 = f0 / norm(f0, '2')
  if(!is.null(signal)) {
    signorm = norm(signal, '2')
  }
  
  if (is.null(data)) {
    data = matrix(rnorm(n^2, sd = 1 /  sqrt(2)), n, n)
    if (!is.null(signal)) {
      data = snr / signorm * outer(signal, signal) + data + t(data)
    } else {
      data = data + t(data)
    }
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
    return( list( iter = iter, 
                  est = est, 
                  Sigma = Sigma, 
                  overlap = overlap, 
                  xcenter = xcenter ) )
  }
}

## Row-wise erasures
# With parameter tuning
amp_goe_erasures = function(data,
                            signal = NULL,
                            init = rep(1, nrow(data)),
                            denoise = function(x, ...) x,
                            deriv = function(x, ...) 1,
                            beta = 1,
                            lambda = 1,
                            maxt = floor(30 / beta),
                            ...) {
  n = length(init)
  normf = norm(init, '2')
  if(!is.null(signal)) {
    signorm = norm(signal, '2')
    if (crossprod(signal, init) < 0) {
      init = -init
    }
  }
  
  sub_size = ceiling(n*beta)

  est = matrix(init / normf, n, 1)
  S = 1:sub_size
  iter = matrix(0, n, 1)
  iter[S,1] = data[S,] %*% init 
  iter[-S,1] = rnorm(n-sub_size, sd = normf)
  if (!is.null(signal)) {
    snr = lambda / signorm * abs(as.numeric(crossprod(est[,1], signal)))
    m = matrix(0, n, 1)
    m[S,1] = snr * signal[S]
  }
  masklist = list(S)
  b = vector('numeric', n)
  
  for (t in 1:(maxt-1)) {
    f = vector('numeric', n)
    f[S] = denoise(iter[S,t] / normf[t], rho = snr[t])
    f[-S] = est[-S,t] * normf[t]
    normf = c(normf, norm(f, '2'))
    b[S] = deriv(iter[S,t] / normf[t], rho = snr[t])
    onsager = sapply(masklist, function(l, b) { sum(b[l]) }, 
                     b = b)
    S = ((S[sub_size]-1 + (1:sub_size)) %% n) + 1; notS = (1:n)[-S]
    masklist = lapply(masklist, function(l, s) { intersect(s, l) },
                      s = notS)
    masklist[[t+1]] = S
    x = data[S,] %*% f - rowSums(est[S,] %*% as.matrix(onsager))
    est = cbind(est, f / normf[t+1])
    iter = cbind(iter, 0)
    iter[S,(t+1)] = x; iter[-S,(t+1)] = iter[-S,t]
    if (!is.null(signal)) {
      snr = c(
        snr, 
        lambda / (signorm * normf[t+1]) * abs(as.numeric(crossprod(f, signal)))
      )
      m = cbind(m, 0)
      m[S, t+1] = snr[t+1] * signal[S]; m[-S, t+1] = m[-S, t] 
    }
    if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
  }
  if (is.null(signal)) {
    return( list( iter = iter,
                  est = est,
                  se = normf) )
  } else {
    return( list( iter = iter,
                  mean = m,
                  est = est,
                  se = normf,
                  overlap = abs(snr) / lambda
                  ) )
  }
}


amp_goe_erasures_random = function(data,
                                   signal = NULL,
                                   init = rep(1, nrow(data)),
                                   denoise = function(x, ...) x,
                                   deriv = function(x, ...) 1,
                                   beta = 1,
                                   lambda = 1,
                                   maxt = floor(30 / beta),
                                   ...) {
  n = length(init)
  sub_size = ceiling(n*beta)
  normf = norm(init, '2')
  if(!is.null(signal)) {
    signorm = norm(signal, '2')
    if (crossprod(signal, init) < 0) {
      init = -init
    }
  } 
  
  est = matrix(init / normf, n, 1)
  S = sort(sample(1:n, sub_size))
  iter = matrix(0, n, 1)
  iter[S,1] = data[S,] %*% init
  iter[-S,1] = rnorm(n-sub_size, sd = normf)
  if (!is.null(signal)) {
    snr = lambda / signorm * abs(as.numeric(crossprod(est[,1], signal)))
    m = matrix(0, n, 1)
    m[S,1] = snr * signal[S]
  }
  masklist = list(S)
  b = vector('numeric', n)
  
  for (t in 1:(maxt-1)) {
    f = vector('numeric', n)
    f[S] = denoise(iter[S,t] / normf[t], rho = snr[t])
    f[-S] = est[-S,t] * normf[t]
    normf = c(normf, norm(f, '2'))
    b[S] = deriv(iter[S,t] / normf[t], rho = snr[t])
    onsager = sapply(masklist, function(l, b) { sum(b[l]) }, 
                     b = b)
    S = sort(sample(1:n, sub_size)); notS = (1:n)[-S]
    masklist = lapply(masklist, function(l, s) { intersect(s, l) },
                      s = notS)
    masklist[[t+1]] = S
    x = data[S,] %*% f - rowSums(est[S,] %*% as.matrix(onsager))
    est = cbind(est, f / normf[t+1])
    iter = cbind(iter, 0)
    iter[S,(t+1)] = x; iter[-S,(t+1)] = iter[-S,t]
    if (!is.null(signal)) {
      snr = c(
        snr, 
        lambda / (signorm * normf[t+1]) * abs(as.numeric(crossprod(f, signal)))
      )
      m = cbind(m, 0)
      m[S, t+1] = snr[t+1] * signal[S]; m[-S, t+1] = m[-S, t] 
    }
    if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
  }
  if (is.null(signal)) {
    return( list( iter = iter,
                  est = est,
                  se = normf) )
  } else {
    return( list( iter = iter,
                  mean = m,
                  est = est,
                  se = normf,
                  overlap = abs(snr) / lambda
    ) )
  }
}

## Row-wise erasures
## No parameter tuning
# amp_goe_erasures = function(f0,
#                             denoise,
#                             deriv,
#                             sub_size = length(f0),
#                             data = NULL,
#                             signal = NULL,
#                             snr = 1,
#                             maxt = 100,
#                             ...) {
#   n = length(f0)
#   f0 = f0 / norm(f0, '2')
#   if(!is.null(signal)) {
#     signorm = norm(signal, '2')
#   }
#   
#   if (is.null(data)) {
#     data = matrix(rnorm(n^2, sd = 1 /  sqrt(2)), n, n)
#     if (!is.null(signal)) {
#       data = snr / signorm * outer(signal, signal) + data + t(data)
#     } else {
#       data = data + t(data)
#     }
#   } 
#   
#   est = matrix(f0, n, 1)
#   S = sort(sample(1:n, sub_size))
#   iter = matrix(0, n, 1)
#   iter[S,1] = data[S,] %*% f0
#   iter[-S,1] = rnorm(n-sub_size)
#   if (!is.null(signal)) {
#     overlap = as.numeric(crossprod(f0, signal) / signorm)
#     m = matrix(0, n, 1)
#     m[S,1] = snr * overlap * signal[S]
#   }
#   masklist = list(S)
#   
#   for (t in 1:(maxt-1)) {
#     f = denoise(iter[,t])
#     normf = norm(f, '2')
#     b = deriv(iter[,t])
#     onsager = sapply(masklist, function(l, b, d) { sum(b[l]) / d }, 
#                      b = b,
#                      d = normf )
#     S = sort(sample(1:n, sub_size)); notS = (1:n)[-S]
#     masklist = lapply(masklist, function(l, s) { intersect(s, l) },
#                       s = notS)
#     masklist[[t+1]] = S
#     x = data[S,] %*% f / normf - rowSums(est[S,] %*% as.matrix(onsager))
#     est = cbind(est, f / normf)
#     iter = cbind(iter, 0)
#     iter[S,(t+1)] = x; iter[-S,(t+1)] = iter[-S,t]
#     if (!is.null(signal)) {
#       overlap = c(overlap, 
#                   as.numeric(crossprod(est[,t+1], signal) / signorm))
#       m = cbind(m, 0)
#       m[S, t+1] = snr * overlap[t+1] * signal[S]; m[-S, t+1] = m[-S, t] 
#     }
#   }
#   if (is.null(signal)) {
#     return( list( iter = iter,
#                   est = est) )
#   } else {
#     return( list( iter = iter,
#                   est = est,
#                   signal = signal,
#                   snr = snr,
#                   overlap = overlap,
#                   mean = m) )
#   }
# }