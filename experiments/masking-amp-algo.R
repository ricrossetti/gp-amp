### AMP algorithm with masking
library(tikzDevice)
sources = list.files(path = './src', pattern = '*.R', full.names = TRUE)
sapply(sources, source)
rm(sources)

### EXPERIMENT 3
### Empirical performance of AMP
amp_pi_mc = function(signal, 
                     init,
                     denoiser = function(x, ...) x,
                     deriv = function(x, ...) 1, 
                     lambda, 
                     beta,
                     sweep = FALSE,
                     iter = 30 / beta,
                     resample = 25, 
                     ...) {
  out = list()
  n = length(signal)
  spike = lambda / norm(signal, '2') * tcrossprod(signal)
  if (sweep) {amp=amp_goe_erasures} else {amp=amp_goe_erasures_random}
  for (i in 1:length(lambda)) {
    for (j in 1:length(beta)) {
      res = matrix(0, iter, resample)
      for (r in 1:resample) {
        data = matrix(rnorm(n^2, sd = 1/sqrt(2)), n, n)
        data = spike + data + t(data)
        ov = amp(data, 
                 signal, 
                 init, 
                 denoiser, 
                 deriv, 
                 beta = beta[j], 
                 lambda = lambda[i],
                 maxt = iter, 
                 ...)$overlap
        res[,r] = ov
        print(paste(r, "resamplings complete."))
        attr(res, 'lambda') = i
        attr(res, 'beta') = j
        attr(res, 'sweep') = sweep
        
      }
    }
    out[[length(out)+1]] = res
  }
  return(out)
} 

ov_theory = function(initsnr,
                     ov_fun,
                     lambda,
                     beta,
                     sweep = FALSE,
                     iter = 30,
                     ...) {
  res = list()
  if (sweep) {f=se_robin; g=se_erasures_pi_robin} 
  else {f=se_erasures; g=se_erasures_pi}
  for (i in 1:length(lambda)) {
    for (j in 1:length(beta)) {
      dyn = f(initsnr,
              ov_fun,
              lambda[i],
              beta[j],
              iter = iter / beta[j],
              ...)
      attr(dyn, 'lambda') = i
      attr(dyn, 'beta') = j
      attr(dyn, 'l') = ceiling(length(dyn) * beta[j])
      attr(dyn, 'sweep') = sweep
      attr(dyn, 'bayes') = TRUE
      dyn_pi = g(initsnr,
                 lambda[i],
                 beta[j],
                 iter = iter / beta[j])
      attr(dyn_pi, 'lambda') = i
      attr(dyn_pi, 'beta') = j
      attr(dyn_pi, 'l') = ceiling(length(dyn_pi) * beta[j])
      attr(dyn_pi, 'sweep') = sweep
      attr(dyn_pi, 'bayes') = FALSE
      res[[length(res)+1]] = dyn_pi
      res[[length(res)+1]] = dyn
    }
  }
  return( res )
}

# Rademacher signal
n = 5000
signal = sample(c(-1,1), n, TRUE)
initsnr = 1e-6
init = tanh_denoiser(sqrt(initsnr) * signal + rnorm(n), rho = initsnr)
init = rep(1,n)
la = 1.5
be = 1

a = ov_theory(initsnr, sparse_rad_ov, la, be, p=1)

a = amp_pi_mc(signal,
              rep(1,n),
              tanh_denoiser,
              tanh_denoiser_deriv,
              lambda = la,
              beta = be)

b = amp_pi_mc(signal,
              init,
              tanh_denoiser,
              tanh_denoiser_deriv,
              lambda = la,
              beta = be)

