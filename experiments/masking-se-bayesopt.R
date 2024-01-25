### Bayes-optimal SE experiments with masking
library(tikzDevice)
sources = list.files(path = './src', pattern = '*.R', full.names = TRUE)
sapply(sources, source)
rm(sources)

### EXPERIMENT 1
### Overlap dynamics at varying batch sizes
# Function to compute trajectories
ov_dyn = function(initsnr,
                  ov_fun,
                  lambda,
                  beta,
                  sweep = FALSE,
                  ...) {
  res = list()
  if (sweep) {f=se_robin} else {f=se_erasures}
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
      res[[length(res)+1]] = dyn
    }
  }
  return( res )
}
# Function to plot results
ov_plot = function(ov_list, lambda, beta, legend = TRUE) {
  plot.new()
  xmax = max(sapply(ov_list, attr, which = 'l'))
  par(oma = c(0,0,0,0), mar = c(2,2,1,1))
  plot.window(xlim = c(1,xmax),
              ylim = c(0,1.05),
              xaxs = 'i', yaxs = 'i')
  box()
  lapply(ov_list, function(y) {
    l = if (beta[attr(y, 'beta')]==1) {
      1
    } else if (attr(y, 'sweep')==TRUE) {
      2
    } else 3
    lines( (0:(length(y)-1)) * beta[attr(y, 'beta')],
           y,
           col = attr(y, 'beta'),
           lty = l,
           lwd = 1.5)
  })
  lapply(ov_list, function(y) {
    pt = if (beta[attr(y, 'beta')]==1) {
      18
    } else if (attr(y, 'sweep')==TRUE) {
      13
    } else 7
    points( (length(y)-1) * beta[attr(y, 'beta')],
           y[length(y)],
           col = attr(y, 'beta'),
           pch = pt,
           cex = 1.5)
  })
  axis(1, seq(0, xmax, by=10), cex.axis = .5, padj = -2.5, tcl = -.3)
  axis(2, seq(0,1, by = 0.1), cex.axis = .75, padj = 1.5, tcl = -.3)
  sapply( 1:length(lambda), function(x) {
    a = Filter( function(y) attr(y,'lambda') == x, ov_list)
    xlpos = max( sapply( a, function(y) attr(y, 'l') ) )
    ylpos = max(a[[1]])
    text( xlpos, ylpos, 
          labels = paste0('$\\lambda=',lambda[x],'$'),
          adj = c(1.2,-.75) )
    return(list(xlpos,ylpos))
  } )
  if (legend) {
    legend('topright',
           legend=paste0('$\\beta=', as.character(beta),"$"),
           col = 1:length(beta),
           lty = 1,
           lwd = 2)
  }
}
# Initial perturbation size
eps = 1e-6
# Batch sizes
be = c(1, .5, .1, .01)
# SNR levels
la = c(1.1, 1.2, 1.5)

## Gaussian signal
gaus_dynamics = c(ov_dyn(eps, gaus_ov, la, be), 
                  ov_dyn(eps, gaus_ov, la, be, TRUE))
tikz(file = './fig/se_bayes_optimal_gaus.tex', width = 3, height = 3,
     documentDeclaration = '\\documentclass{standalone}\n', standAlone = TRUE)
ov_plot(gaus_dynamics, la, be)
dev.off()

## Rademacher signal
rad_dynamics = c(ov_dyn(eps, sparse_rad_ov, la, be, p=1), 
                  ov_dyn(eps, sparse_rad_ov, la, be, TRUE, p=1))
tikz(file = './fig/se_bayes_optimal_rad.tex', width = 3, height = 3,
     documentDeclaration = '\\documentclass{standalone}\n', standAlone = TRUE)
ov_plot(rad_dynamics, la, be, legend = FALSE)
dev.off()

## Sparse Rademacher signal
sparse_rad_dynamics = c(ov_dyn(eps, sparse_rad_ov, la, be, p=.1), 
                 ov_dyn(eps, sparse_rad_ov, la, be, TRUE, p=.1))
tikz(file = './fig/se_bayes_optimal_sparse_rad.tex', width = 3, height = 3,
     documentDeclaration = '\\documentclass{standalone}\n', standAlone = TRUE)
ov_plot(sparse_rad_dynamics, la, be, legend = FALSE)
dev.off()