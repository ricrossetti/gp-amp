### PRIOR SPECIFIC FUNCTIONS
## Denoiser and derivative
denoise = function(x, a, b) {
  tanh(a / b * x)
}

deriv = function(x, a, b) {
  (a / b) / cosh(a / b * x)^2
}

denoise_linear = function(x, a, b) {
  a / (a^2 + b) * x
}

deriv_linear = function(a, b) {
  a / (a^2 + b)
}
 
## Overlap function
ov_fun = function(snr) {
  integrand = function(x) {
    exp( log( tanh(x)^2 ) + dnorm(x, 
                                  mean = snr,
                                  sd = sqrt(snr),
                                  log = TRUE) ) 
  }
  integrate(integrand, -Inf, Inf)$value
}

ov_fun_2 = function(lambda, a, b) {
  integrand = function(x) {
    exp( log( tanh((lambda*a/b)*x)^2 ) + dnorm(x, 
                                  mean = lambda*a,
                                  sd = sqrt(b),
                                  log = TRUE) ) 
  }
  integrate(integrand, -Inf, Inf)$value
}

### AMP FUNCTIONS
## Row-wise erasures
# With parameter tuning
amp_bayes = function(data,
                     signal,
                     init,
                     beta,
                     lambda,
                     robin = TRUE,
                     maxt = 30 / beta) {
  n = length(signal)
  
  sub_size = ceiling(n*beta)
  r = as.numeric(crossprod(signal, init)) / n
  q = as.numeric(crossprod(init)) / n
  
  theory_r = se_bayes(r^2 / q, beta, lambda, robin, maxt)
  
  est = matrix(init, n, 1)
  iter = matrix(data %*% est, n, 1)
  
  S=1:n
  masklist = list(S)
  b = vector('numeric', n)
  
  for (t in 1:maxt) {
    f = vector('numeric', n)
    f[S] = denoise(iter[S,t], lambda*r[t], q[t])
    f[-S] = est[-S,t]
    b[S] = deriv(iter[S,t], lambda*r[t], q[t])
    onsager = sapply(masklist, function(l, b) { sum(b[l]) }, 
                     b = b) / n
    if (robin) {
      S = ((S[sub_size]-1 + (1:sub_size)) %% n) + 1
    } else {
      S = sort(sample(1:n, sub_size))
    }
    notS = (1:n)[-S]
    masklist = lapply(masklist, function(l, notS) { intersect(notS, l) },
                      notS = notS)
    masklist[[t+1]] = S
    x = data[S,] %*% f - est[S,] %*% as.matrix(onsager)
    est = cbind(est, f)
    iter = cbind(iter, 0)
    iter[S,(t+1)] = x; iter[-S,(t+1)] = iter[-S,t]
    r[t+1] = crossprod(f, signal) / n    
    q[t+1] = crossprod(f) / n 
    if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
  }
  return( list( iter = iter,
                est = est,
                r = r,
                q = q,
                r_th = theory_r
  ) )
}

amp_linear_sphere = function(data,
                             signal,
                             init,
                             beta,
                             lambda,
                             robin = TRUE,
                             maxt = 30 / beta) {
  n = length(signal)
  init = init / norm(init, '2') * sqrt(n)
  
  sub_size = ceiling(n*beta)
  r = as.numeric(crossprod(signal, init)) / n
  
  theory_r = se_linear(r, beta, lambda, robin, maxt)
  
  est = matrix(init, n, 1)
  iter = matrix(data %*% est, n, 1)
  
  S=1:n
  masklist = list(S)
  b = vector('numeric', n)
  
  for (t in 1:maxt) {
    f = iter[,t] / norm(iter[,t], '2') * sqrt(n)
    b[S] = sqrt(n) / (norm(iter[,t], '2'))
    onsager = sapply(masklist, function(l, b) { sum(b[l]) }, 
                     b = b) / n
    if (robin) {
      S = ((S[sub_size]-1 + (1:sub_size)) %% n) + 1
    } else {
      S = sort(sample(1:n, sub_size))
    }
    notS = (1:n)[-S]
    masklist = lapply(masklist, function(l, notS) { intersect(notS, l) },
                      notS = notS)
    masklist[[t+1]] = S
    x = data[S,] %*% f - est[S,] %*% as.matrix(onsager)
    est = cbind(est, f)
    iter = cbind(iter, 0)
    iter[S,(t+1)] = x; iter[-S,(t+1)] = iter[-S,t]
    r[t+1] = crossprod(f, signal) / n
    if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
  }
  return( list( iter = iter,
                est = est,
                r = r,
                r_th = theory_r
  ) )
}

amp_linear_sphere2 = function(data,
                              signal,
                              init,
                              beta,
                              lambda,
                              robin = TRUE,
                              maxt = 30 / beta) {
  n = length(signal)
  init = init / norm(init, '2') * sqrt(n)
  
  sub_size = ceiling(n*beta)
  r = as.numeric(crossprod(signal, init)) / n
  
  theory_r = se_linear_sphere(r, beta, lambda, robin, maxt)
  
  est = matrix(init, n, 1)
  iter = matrix(data %*% est, n, 1)
  
  S=1:n
  masklist = list(S)
  
  for (t in 1:maxt) {
    f = iter[,t] / norm(iter[,t], '2') * sqrt(n)
    b = sqrt(n) / (norm(iter[,t], '2'))
    onsager = b * sapply(masklist, function(l) { length(l) }) / n
    if (robin) {
      S = ((S[sub_size]-1 + (1:sub_size)) %% n) + 1
    } else {
      S = sort(sample(1:n, sub_size))
    }
    notS = (1:n)[-S]
    masklist = lapply(masklist, function(l, notS) { intersect(notS, l) },
                      notS = notS)
    masklist[[t+1]] = S
    x = data[S,] %*% f - est[S,] %*% as.matrix(onsager)
    est = cbind(est, f)
    iter = cbind(iter, 0)
    iter[S,(t+1)] = x; iter[-S,(t+1)] = iter[-S,t]
    r[t+1] = crossprod(f, signal) / n
    if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
  }
  return( list( iter = iter,
                est = est,
                r = r,
                r_th = theory_r
  ) )
}

### STATE EVOLUTION FUNCTIONS
se_bayes = function(r,
                    beta,
                    lambda,
                    robin = TRUE,
                    maxt) {
  phi = function(r, l) ov_fun(l^2*r)
  if (robin) {
    eta = function(t, beta) {
      if (t<1/beta) {
        c( 1 - beta * min(t, 1/beta), rep(beta, t) )
      } else {
        c( rep(0, t + 1 - 1/beta), rep(beta, 1/beta) )
      }
    }
  } else {
    eta = function(t, beta) (1-beta)^rev((0:t)) * beta^((0:t)>0)
  }
  for( t in 1:maxt ) {
    phivec = sapply(r, phi, l = lambda)
    w = eta(t-1, beta)
    r = c(r, sum(w*phivec))
  }
  return(r)
}

se_linear_sphere = function(r,
                            beta,
                            lambda,
                            robin = TRUE,
                            maxt) {
  if (robin) {
    eta = function(t, beta) {
      if (t<1/beta) {
        c( 1 - beta * min(t, 1/beta), rep(beta, t) )
      } else {
        c( rep(0, t + 1 - 1/beta), rep(beta, 1/beta) )
      }
    }
  } else {
    eta = function(t, beta) (1-beta)^rev((0:t)) * beta^((0:t)>0)
  }
  for( t in 1:maxt ) {
    w = eta(t-1, beta)
    r = c(r, lambda*sum(w*r) / sqrt(sum(w * (lambda^2*r^2 + 1)) ))
  }
  return(r)
}

se_linear = function(r,
                     q,
                     beta,
                     lambda,
                     robin = TRUE,
                     maxt) {
  if (robin) {
    eta = function(t, beta) {
      if (t<1/beta) {
        c( 1 - beta * min(t, 1/beta), rep(beta, t) )
      } else {
        c( rep(0, t + 1 - 1/beta), rep(beta, 1/beta) )
      }
    }
  } else {
    eta = function(t, beta) (1-beta)^rev((0:t)) * beta^((0:t)>0)
  }
  for( t in 1:maxt ) {
    p = eta(t-1, beta)
    r = c(r, lambda * sum(r * p) )
    q = c(q, sum( p * (lambda^2 * r[-(t+1)]^2 + q) ))
  }
  return(r / sqrt(q))
}

### OLDSTUFF ####################################
# amp_robin_bayes = function(data,
#                            signal,
#                            initsnr,
#                            beta,
#                            lambda,
#                            maxt) {
#   n = length(signal)
#   init_obs = sqrt(initsnr) * signal + rnorm(n)
#   init = denoise(init_obs, initsnr)
#   if(!is.null(signal)) {
#     if (crossprod(signal, init) < 0) {
#       init = -init
#     }
#   }
#   
#   sub_size = ceiling(n*beta)
#   theory_ov = se_rand_bayes(initsnr, beta, lambda)
#   theory_snr = theory_ov^2*lambda^2
#   MAXT = min(maxt, length(theory_snr))
#   
#   normf = norm(init, '2')
#   est = matrix(init / normf, n, 1)
#   iter = matrix(data %*% est, n, 1)
#   snr = lambda / sqrt(n) * as.numeric(crossprod(signal, est))
#   m = matrix(snr * signal, n, 1)
#   
#   S=1:n
#   masklist = list(S)
#   b = vector('numeric', n)
#   
#   for (t in 1:(MAXT-1)) {
#     f = vector('numeric', n)
#     f[S] = denoise(iter[S,t], snr[t]^2)
#     f[-S] = est[-S,t] * normf[t]
#     normf = c(normf, norm(f, '2'))
#     b[S] = deriv(iter[S,t], snr[t]^2) / normf[t+1]
#     onsager = sapply(masklist, function(l, b) { sum(b[l]) }, 
#                      b = b)
#     S = ((S[sub_size]-1 + (1:sub_size)) %% n) + 1; notS = (1:n)[-S]
#     masklist = lapply(masklist, function(l, notS) { intersect(notS, l) },
#                       notS = notS)
#     masklist[[t+1]] = S
#     x = data[S,] %*% f / normf[t+1] - est[S,] %*% as.matrix(onsager)
#     est = cbind(est, f / normf[t+1])
#     iter = cbind(iter, 0)
#     iter[S,(t+1)] = x; iter[-S,(t+1)] = iter[-S,t]
#     snr = c(
#       snr, 
#       lambda / sqrt(n) * as.numeric(crossprod(signal, est[,t+1]))
#     )
#     m = cbind(m, 0)
#     m[S, t+1] = sqrt(snr[t+1]) * signal[S]; m[-S, t+1] = m[-S, t] 
#     if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
#   }
#   return( list( iter = iter,
#                 mean = m,
#                 est = est,
#                 se = normf,
#                 ov = snr / lambda,
#                 snr = snr^2,
#                 theory_ov = theory_ov,
#                 theory_snr = theory_snr
#   ) )
# }
# 
# amp_random_bayes = function(data,
#                             signal,
#                             initsnr,
#                             beta,
#                             lambda,
#                             maxt) {
#   n = length(signal)
#   init_obs = sqrt(initsnr) * signal + rnorm(n)
#   init = denoise(init_obs, initsnr)
#   if(!is.null(signal)) {
#     if (crossprod(signal, init) < 0) {
#       init = -init
#     }
#   }
#   
#   sub_size = ceiling(n*beta)
#   theory_ov = se_robin_bayes(initsnr, beta, lambda)
#   theory_snr = theory_ov^2*lambda^2
#   MAXT = min(maxt, length(theory_snr))
#   
#   normf = norm(init, '2')
#   est = matrix(init / normf, n, 1)
#   iter = matrix(data %*% est, n, 1)
#   snr = lambda / sqrt(n) * as.numeric(crossprod(signal, est))
#   m = matrix(snr * signal, n, 1)
#   
#   S=1:n
#   masklist = list(S)
#   b = vector('numeric', n)
#   
#   for (t in 1:(MAXT-1)) {
#     f = vector('numeric', n)
#     f[S] = denoise(iter[S,t], snr[t]^2)
#     f[-S] = est[-S,t] * normf[t]
#     normf = c(normf, norm(f, '2'))
#     b[S] = deriv(iter[S,t], snr[t]^2) / normf[t+1]
#     onsager = sapply(masklist, function(l, b) { sum(b[l]) }, 
#                      b = b)
#     S = sort(sample(1:n, sub_size)); notS = (1:n)[-S]
#     masklist = lapply(masklist, function(l, notS) { intersect(notS, l) },
#                       notS = notS)
#     masklist[[t+1]] = S
#     x = data[S,] %*% f / normf[t+1] - est[S,] %*% as.matrix(onsager)
#     est = cbind(est, f / normf[t+1])
#     iter = cbind(iter, 0)
#     iter[S,(t+1)] = x; iter[-S,(t+1)] = iter[-S,t]
#     snr = c(
#       snr, 
#       lambda / sqrt(n) * as.numeric(crossprod(signal, est[,t+1]))
#     )
#     m = cbind(m, 0)
#     m[S, t+1] = sqrt(snr[t+1]) * signal[S]; m[-S, t+1] = m[-S, t] 
#     if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
#   }
#   return( list( iter = iter,
#                 mean = m,
#                 est = est,
#                 se = normf,
#                 ov = snr / lambda,
#                 snr = snr^2,
#                 theory_ov = theory_ov,
#                 theory_snr = theory_snr
#   ) )
# }
# 

# se_rand_bayes = function(r,
#                          beta,
#                          lambda,
#                          maxt) {
#   phi = function(r, l) ov_fun(l^2*r)
#   eta = function(t, beta) (1-beta)^rev((0:t)) * beta^((0:t)>0)
#   improve = Inf
#   for( t in 1:maxt ) {
#     phivec = sapply(r, phi, l = lambda)
#     w = eta(t-1, beta)
#     r = c(r, sum(w*phivec))
#   }
#   return(r)
# }
