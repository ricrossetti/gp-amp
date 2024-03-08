### AMP algorithm with masking
library(tikzDevice)
sources = list.files(path = './src', pattern = '*.R', full.names = TRUE)
sapply(sources, source)
rm(sources)

set.seed(1)

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
      res = matrix(0, iter[j], resample)
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
                 maxt = iter[j], 
                 ...)$overlap
        res[,r] = ov
        print(paste(r, "resamplings complete."))
      }
      attr(res, 'lambda') = i
      attr(res, 'beta') = j
      attr(res, 'sweep') = sweep
      out[[length(out)+1]] = res
    }
    
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

plot_amp = function(ov_list, amp_list, lambda, beta, legend = TRUE) {
  plot.new()
  xmax = max(sapply(ov_list, attr, which = 'l'))
  par(oma = c(0,0,0,0), mar = c(3,3,1,1))
  plot.window(xlim = c(0,xmax),
              ylim = c(0,1.05),
              xaxs = 'i', yaxs = 'i')
  box()
  lapply(ov_list, function(y) {
    lines( (0:(length(y)-1)) * beta[attr(y, 'beta')],
           y,
           col = attr(y, 'beta'),
           lwd = 1.5)
  })
  
  lapply(amp_list, function(x) {
    subset = (0:(xmax-1)) / beta[attr(x, 'beta')] + 1
    y = x[subset,]
    m = apply(y, 1, mean)
    stdev = apply(y, 1, sd)
    points((0:(nrow(y)-1)) + 5*1e-2 * (attr(x, 'beta')-1),
           m,
           col = attr(x, 'beta'),
           pch = attr(x, 'beta'))
    arrows((0:(nrow(y)-1)) + 5*1e-2 * (attr(x, 'beta')-1),
           y0 = m - stdev, y1 = m + stdev,
           col = attr(x, 'beta'),
           code = 3,
           length = 0.05,
           angle = 90)
  })
  axis(1, seq(0, xmax, by=5), cex.axis = .5, padj = -2.5, tcl = -.3)
  axis(2, seq(0,1, by = 0.1), cex.axis = .75, padj = 1.5, tcl = -.3)
  title(xlab = "Normalized iteration $t\\beta$", 
        ylab = "Normalized overlaps $o_t$", cex.lab = .75, line = 1.1)
  if (legend) {
    legend('topright',
           legend=paste0('$\\beta=', as.character(beta),"$"),
           col = 1:length(beta),
           lty = 1,
           lwd = 2)
  }
}

# Rademacher signal
n = 5000
signal = sample(c(-1,1), n, TRUE, prob = c(.5,.5))
maxt = 30
initsnr = 1e-2
init_oracle = tanh_denoiser(sqrt(initsnr) * signal + rnorm(n), rho = initsnr)
init = rep(1,n)
la = 1.5
be = c(1, .5)

# AMP ALGOS
amp_robin_bayes = amp_pi_mc(signal,
                            init_oracle,
                            tanh_denoiser,
                            tanh_denoiser_deriv,
                            la,
                            be,
                            sweep = TRUE,
                            iter = maxt / be)
amp_rand_bayes = amp_pi_mc(signal,
                           init_oracle,
                           tanh_denoiser,
                           tanh_denoiser_deriv,
                           la,
                           be,
                           sweep = FALSE,
                           iter = maxt / be)
amp_robin_linear = amp_pi_mc(signal,
                             init_oracle,
                             lambda = la,
                             beta = be,
                             sweep = TRUE,
                             iter = maxt / be)
amp_rand_linear = amp_pi_mc(signal,
                            init_oracle,
                            lambda = la,
                            beta = be,
                            sweep = FALSE,
                            iter = maxt / be)

# COMPARISON SE
se_comparison_robin = ov_theory(initsnr, sparse_rad_ov, la, be, 
                                iter = maxt, sweep = TRUE, p=1)
se_comparison_robin_linear = Filter(function(x) !attr(x, 'bayes'),
                                    se_comparison_robin)
se_comparison_robin_bayes = Filter(function(x) attr(x, 'bayes'),
                                    se_comparison_robin)

se_comparison_rand = ov_theory(initsnr, sparse_rad_ov, la, be,
                               iter = maxt, sweep = FALSE, p=1)
se_comparison_rand_linear = Filter(function(x) !attr(x, 'bayes'),
                                    se_comparison_rand)
se_comparison_rand_bayes = Filter(function(x) attr(x, 'bayes'),
                                   se_comparison_rand)



