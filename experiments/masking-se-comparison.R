### Comparison of Bayes-optimal SE and power method
library(tikzDevice)
sources = list.files(path = './src', pattern = '*.R', full.names = TRUE)
sapply(sources, source)
rm(sources)

### EXPERIMENT 2
## Bayes-optimal vs identity denoisers
ov_comparison = function(initsnr,
                         ov_fun,
                         lambda,
                         beta,
                         sweep = FALSE,
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
              ...)
      attr(dyn, 'lambda') = i
      attr(dyn, 'beta') = j
      attr(dyn, 'l') = ceiling(length(dyn) * beta[j])
      attr(dyn, 'sweep') = sweep
      attr(dyn, 'bayes') = TRUE
      dyn_pi = g(initsnr,
                 lambda[i],
                 beta[j])
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

pi_bayes_comparison_plot = function(ov_list, lambda, beta, legend = TRUE) {
  plot.new()
  xmax = max(sapply(ov_list, attr, which = 'l'))
  par(oma = c(0,0,0,0), mar = c(3,3,1,1))
  plot.window(xlim = c(1,xmax),
              ylim = c(0,1.05),
              xaxs = 'i', yaxs = 'i')
  box()
  lapply(ov_list, function(y) {
    l = if (attr(y, 'bayes')==TRUE) 1 else 2
    lines( (0:(length(y)-1)) * beta[attr(y, 'beta')],
           y,
           col = attr(y, 'beta'),
           lty = l,
           lwd = 1.5)
  })
  lapply(ov_list, function(y) {
    points( (length(y)-1) * beta[attr(y, 'beta')],
            y[length(y)],
            col = attr(y, 'beta'),
            pch = 13,
            cex = 1.5)
  })
  axis(1, seq(0, xmax, by=5), cex.axis = .5, padj = -2.5, tcl = -.3)
  axis(2, seq(0,1, by = 0.1), cex.axis = .75, padj = 1.5, tcl = -.3)
  # sapply( 1:length(lambda), function(x) {
  #   a = Filter( function(y) attr(y,'lambda') == x, ov_list)
  #   xlpos = max( sapply( a, function(y) attr(y, 'l') ) )
  #   ylpos = max(a[[1]])
  #   text( xlpos, ylpos, 
  #         labels = paste0('$\\lambda=',lambda[x],'$'),
  #         adj = c(1.2,-.75) )
  #   return(list(xlpos,ylpos))
  # } )
  if (legend) {
    legend('bottomright',
           legend=paste0('$\\beta=', as.character(beta),"$"),
           col = 1:length(beta),
           lty = 1,
           lwd = 2)
  }
  title(xlab = "Normalized iteration $t\\beta$", 
        ylab = "Normalized overlaps $o_t$", cex.lab = .75, line = 1.1)
}

# Rademacher spike
la = 1.5
be = c(1, .1, .01)
sweep_comp_rad = ov_comparison(1e-6, sparse_rad_ov, la, be, sweep = TRUE, p = 1)
tikz(file = './fig/se_roundrobin_comparison_rad.tex', width = 3, height = 3,
     documentDeclaration = '\\documentclass{standalone}\n', standAlone = TRUE)
pi_bayes_comparison_plot(sweep_comp_rad, la, be)
dev.off()

rand_comp_rad = ov_comparison(1e-6, sparse_rad_ov, la, be, sweep = FALSE, p = 1)
tikz(file = './fig/se_random_comparison_rad.tex', width = 3, height = 3,
     documentDeclaration = '\\documentclass{standalone}\n', standAlone = TRUE)
pi_bayes_comparison_plot(rand_comp_rad, la, be)
dev.off()

# Gaussian spike
la = 1.5
be = c(1, .1, .01)
sweep_comp_gaus = ov_comparison(1e-6, gaus_ov, la, be, sweep = TRUE)
tikz(file = './fig/se_roundrobin_comparison_gaus.tex', width = 3, height = 3,
     documentDeclaration = '\\documentclass{standalone}\n', standAlone = TRUE)
pi_bayes_comparison_plot(sweep_comp_gaus, la, be)
dev.off()

rand_comp_gaus = ov_comparison(1e-6, gaus_ov, la, be, sweep = FALSE)
tikz(file = './fig/se_random_comparison_gaus.tex', width = 3, height = 3,
     documentDeclaration = '\\documentclass{standalone}\n', standAlone = TRUE)
pi_bayes_comparison_plot(rand_comp_gaus, la, be)
dev.off()
