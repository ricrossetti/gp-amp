### Numerics with Rademacher spike
source("./src/numerics-rad-new.R")
library(tikzDevice)

if (file.exists("./rad-numerics.RData")) {
  load("./rad-numerics.RData")
} else {
  
  set.seed(27)
  n = 5000
  la = sqrt(2)
  r0 = 0.1
  be = c(1,.1)
  maxt = 15
  resample = 100
  signal = sample(c(-1,1), n, TRUE)
  init = denoise( r0*signal + rnorm(n), r0, 1)
  
  ### Bayes-optimal experiments
  bayes_experiments_norm = function(signal,
                               init,
                               beta,
                               lambda,
                               robin,
                               maxt) {
    res = list()
    for (j in 1:length(beta)) {
      rmat = matrix(0, maxt / beta[j] + 1, resample)
      rmat2 = matrix(0, maxt / beta[j] + 1, resample)
      for (r in 1:resample) {
        data = matrix(rnorm(n^2, sd = 1/sqrt(2)),n)
        data = lambda / n * tcrossprod(signal) + (data + t(data)) / sqrt(n)
        amp = amp_bayes(data, 
                        signal, 
                        init, 
                        beta[j], 
                        la, 
                        robin, 
                        maxt / beta[j])
        rmat[,r] = (amp$r / sqrt(amp$q))^2
        print(paste(r, 'resamples complete.'))
      }
      th = amp$r_th
      attr(th, 'beta') = j; attr(th, 'theory') = TRUE; attr(th, 'robin') = robin
      attr(rmat, 'beta') = j; attr(rmat, 'theory') = FALSE; attr(rmat, 'robin') = robin
      res[[2*(j-1) + 1]] = th
      res[[2*(j-1) + 2]] = rmat
    }
    return(res)
  }
  
  # Run code
  amp_bayes_norm_robin = bayes_experiments_norm(signal, init, be, la, TRUE, maxt)
  amp_bayes_norm_rand = bayes_experiments_norm(signal, init, be, la, FALSE, maxt)[-c(1:2)]
  amp_bayes_norm_full = c(amp_bayes_norm_robin, amp_bayes_norm_rand)
  
  ### Bayes-optimal experiments, raw iteration count
  ### Raw iteration count power iteration experiments
  bayes_experiments = function(signal,
                               init,
                               beta,
                               lambda,
                               robin,
                               maxt) {
    res = list()
    for (j in 1:length(beta)) {
      rmat = matrix(0, maxt + 1, resample)
      rmat2 = matrix(0, maxt + 1, resample)
      for (r in 1:resample) {
        data = matrix(rnorm(n^2, sd = 1/sqrt(2)),n)
        data = lambda / n * tcrossprod(signal) + (data + t(data)) / sqrt(n)
        amp = amp_bayes(data, 
                        signal, 
                        init, 
                        beta[j], 
                        la, 
                        robin, 
                        maxt)
        rmat[,r] = (amp$r / sqrt(amp$q))^2
        print(paste(r, 'resamples complete.'))
      }
      th = amp$r_th
      attr(th, 'beta') = j; attr(th, 'theory') = TRUE; attr(th, 'robin') = robin
      attr(rmat, 'beta') = j; attr(rmat, 'theory') = FALSE; attr(rmat, 'robin') = robin
      res[[2*(j-1) + 1]] = th
      res[[2*(j-1) + 2]] = rmat
    }
    return(res)
  }
  
  ### Run code
  amp_bayes_robin = bayes_experiments(signal, init, be, la, TRUE, 4*maxt)
  amp_bayes_rand = bayes_experiments(signal, init, be, la, FALSE, 4*maxt)[-c(1:2)]
  amp_bayes_full = c(amp_bayes_robin, amp_bayes_rand)
  
  ### Iteration-normalized power iteration experiments
  pi_experiments_norm = function(signal,
                                 init,
                                 beta,
                                 lambda,
                                 robin,
                                 maxt) {
    res = list()
    for (j in 1:length(beta)) {
      rmat = matrix(0, maxt / beta[j] + 1, resample)
      rmat2 = matrix(0, maxt / beta[j] + 1, resample)
      for (r in 1:resample) {
        data = matrix(rnorm(n^2, sd = 1/sqrt(2)),n)
        data = lambda / n * tcrossprod(signal) + (data + t(data)) / sqrt(n)
        amp = amp_linear_sphere2(data, 
                         signal, 
                         init, 
                         beta[j], 
                         la, 
                         robin, 
                         maxt / beta[j])
        rmat[,r] = amp$r
        print(paste(r, 'resamples complete.'))
      }
      th = amp$r_th
      attr(th, 'beta') = j; attr(th, 'theory') = TRUE; attr(th, 'robin') = robin
      attr(rmat, 'beta') = j; attr(rmat, 'theory') = FALSE; attr(rmat, 'robin') = robin
      res[[2*(j-1) + 1]] = th
      res[[2*(j-1) + 2]] = rmat
    }
    return(res)
  }
  
  # Run code
  amp_pi_norm_robin = pi_experiments_norm(signal, init, be, la, TRUE, maxt)
  amp_pi_norm_rand = pi_experiments_norm(signal, init, be, la, FALSE, maxt)[-c(1:2)]
  amp_pi_norm_full = c(amp_pi_norm_robin, amp_pi_norm_rand)
  
  ### Raw iteration count power iteration experiments
  pi_experiments = function(signal,
                            init,
                            beta,
                            lambda,
                            robin,
                            maxt) {
    res = list()
    for (j in 1:length(beta)) {
      rmat = matrix(0, maxt + 1, resample)
      rmat2 = matrix(0, maxt + 1, resample)
      for (r in 1:resample) {
        data = matrix(rnorm(n^2, sd = 1/sqrt(2)),n)
        data = lambda / n * tcrossprod(signal) + (data + t(data)) / sqrt(n)
        amp = amp_linear_sphere2(data, 
                         signal, 
                         init, 
                         beta[j], 
                         la, 
                         robin, 
                         maxt)
        rmat[,r] = amp$r
        print(paste(r, 'resamples complete.'))
      }
      th = amp$r_th
      attr(th, 'beta') = j; attr(th, 'theory') = TRUE; attr(th, 'robin') = robin
      attr(rmat, 'beta') = j; attr(rmat, 'theory') = FALSE; attr(rmat, 'robin') = robin
      res[[2*(j-1) + 1]] = th
      res[[2*(j-1) + 2]] = rmat
    }
    return(res)
  }
  
  ### Run code
  amp_pi_robin = pi_experiments(signal, init, be, la, TRUE, 4*maxt)
  amp_pi_rand = pi_experiments(signal, init, be, la, FALSE, 4*maxt)[-c(1:2)]
  amp_pi_full = c(amp_pi_robin, amp_pi_rand)
  
  save.image("rad-numerics.RData")
}

### plotting functions
plot_bayes_norm = function(data, lambda, beta, legend = TRUE, maxt, maxline = NULL) {
  plot.new()
  xmax = maxt
  par(oma = c(0,0,0,0), mar = c(3,3,1,1))
  plot.window(xlim = c(0,xmax),
              ylim = c(0,1),
              xaxs = 'i', yaxs = 'i')
  box()
  cols = c(1,2,4)
  pts = c(1,0,5)
  ov = Filter(function(x) attr(x, 'theory'), data)
  for (i in 1:length(ov)) {
    attr(ov[[i]], 'c') = i
  }
  if (!is.null(maxline)) {
    abline(h = maxline, lty = 3, lwd = .3)
  }

  lapply(ov, function(y) {
    lines( (0:(length(y)-1)) * beta[attr(y, 'beta')],
           sqrt(y),
           col = cols[attr(y, 'c')],
           lwd = .5,
           lty = 1)
  })

  amp = Filter(function(x) !attr(x, 'theory'), data)
  for (i in 1:length(amp)) {
    attr(amp[[i]], 'c') = i
  }
  lapply(amp, function(x) {
    subset = (0:xmax) / beta[attr(x, 'beta')] + 1
    y = sqrt(x[subset,])
    m = apply(y, 1, mean)
    stdev = apply(y, 1, function(y) sd(y) / sqrt(length(y)) )
    points((0:(nrow(y)-1)),
           m,
           col = cols[attr(x, 'c')],
           pch = pts[attr(x, 'c')],
           cex = 1,
           lwd = .5)
    arrows((0:(nrow(y)-1)),
           y0 = m - stdev, y1 = m + stdev,
           col = cols[attr(x, 'c')],
           code = 3,
           length = 0.02,
           angle = 90,
           lwd = .5)
  })
  axis(1, seq(0, xmax, by=1), cex.axis = .85, padj = -1.5, tcl = -.2)
  axis(2, seq(0,1, by = 0.1), cex.axis = .85, padj = 1.5, tcl = -.2)
  title(xlab = "$\\# \\{$row updates$\\} / n$",
        ylab = "Overlap", cex.lab = 1, line = 1.2)
  if (legend) {
    lge = sapply(amp, function(x, beta) {
      if (attr(x, 'beta') == 1) {
        paste0('$\\gamma=', beta[attr(x, 'beta')],"$")
      } else if (attr(x, 'robin')) {
        paste0('$\\gamma=', beta[attr(x, 'beta')],"$, Round Robin")
      } else {
        paste0('$\\gamma=', beta[attr(x, 'beta')],"$, Random Update")
      }
    }, beta = beta)
    legend('bottomright',
           legend=lge,
           col = cols,
           lty = 1,
           pch = pts,
           lwd = 1,
           cex = .9,
           bty = 'n')
  }
}
plot_pi_norm = function(data, lambda, beta, legend = TRUE, maxt, maxline = NULL) {
  plot.new()
  xmax = maxt
  par(oma = c(0,0,0,0), mar = c(3,3,1,1))
  plot.window(xlim = c(0,xmax),
              ylim = c(0,1),
              xaxs = 'i', yaxs = 'i')
  box()
  cols = c(1,2,4)
  pts = c(1,0,5)
  ov = Filter(function(x) attr(x, 'theory'), data)
  for (i in 1:length(ov)) {
    attr(ov[[i]], 'c') = i
  }
  if (!is.null(maxline)) {
    abline(h = maxline, lty = 3, lwd = .3)
  }
  lapply(ov, function(y) {
    lines( (0:(length(y)-1)) * beta[attr(y, 'beta')],
           y,
           col = cols[attr(y, 'c')],
           lwd = .5,
           lty = 1)
  })

  amp = Filter(function(x) !attr(x, 'theory'), data)
  for (i in 1:length(amp)) {
    attr(amp[[i]], 'c') = i
  }
  lapply(amp, function(x) {
    subset = (0:xmax) / beta[attr(x, 'beta')] + 1
    y = x[subset,]
    m = apply(y, 1, mean)
    stdev = apply(y, 1, function(y) sd(y) / sqrt(length(y)) )
    points((0:(nrow(y)-1)),
           m,
           col = cols[attr(x, 'c')],
           pch = pts[attr(x, 'c')],
           lwd = .5)
    arrows((0:(nrow(y)-1)),
           y0 = m - stdev , y1 = m + stdev,
           col = cols[attr(x, 'c')],
           code = 3,
           length = 0.02,
           angle = 90,
           lwd = .5)
  })
  axis(1, seq(0, xmax, by=1), cex.axis = .85, padj = -1.5, tcl = -.2)
  axis(2, seq(0,1, by = 0.1), cex.axis = .85, padj = 1.5, tcl = -.2)
  title(xlab = "$\\# \\{$row updates$\\} / n$",
        ylab = "Overlap", cex.lab = 1, line = 1.2)
  if (legend) {
    lge = sapply(amp, function(x, beta) {
      if (attr(x, 'beta') == 1) {
        paste0('$\\gamma=', beta[attr(x, 'beta')],"$")
      } else if (attr(x, 'robin')) {
        paste0('$\\gamma=', beta[attr(x, 'beta')],"$, Round Robin")
      } else {
        paste0('$\\gamma=', beta[attr(x, 'beta')],"$, Random Update")
      }
    }, beta = beta)
    legend('bottomright',
           legend=lge,
           col = cols,
           lty = 1,
           pch = pts,
           lwd = 1,
           cex = .9,
           bty = 'n')
  }
}
plot_pi = function(data, lambda, beta, legend = TRUE, maxt, maxline = NULL) {
  plot.new()
  xmax = maxt
  par(oma = c(0,0,0,0), mar = c(3,3,1,1))
  plot.window(xlim = c(0,xmax),
              ylim = c(0,1),
              xaxs = 'i', yaxs = 'i')
  box()
  cols = c(1,2,4)
  pts = c(1,0,5)
  ov = Filter(function(x) attr(x, 'theory'), data)
  for (i in 1:length(ov)) {
    attr(ov[[i]], 'c') = i
  }
  if (!is.null(maxline)) {
    abline(h = maxline, lty = 3, lwd = .3)
  }
  lapply(ov, function(y) {
    lines( 0:(length(y)-1),
           y,
           col = cols[attr(y, 'c')],
           lwd = .5,
           lty = 1)
  })

  amp = Filter(function(x) !attr(x, 'theory'), data)
  for (i in 1:length(amp)) {
    attr(amp[[i]], 'c') = i
  }
  lapply(amp, function(x) {
    subset = (0:(xmax/4)) * 4 + 1
    y = x[subset,]
    m = apply(y, 1, mean)
    stdev = apply(y, 1, function(y) sd(y) / sqrt(length(y)) )
    points(subset-1,
           m,
           col = cols[attr(x, 'c')],
           pch = pts[attr(x, 'c')],
           lwd = .5)
    arrows(subset-1,
           y0 = m - stdev, y1 = m + stdev,
           col = cols[attr(x, 'c')],
           code = 3,
           length = 0.02,
           angle = 90,
           lwd = .5)
  })
  axis(1, seq(0, xmax, by=5), cex.axis = .85, padj = -1.5, tcl = -.2)
  axis(2, seq(0,1, by = 0.1), cex.axis = .85, padj = 1.5, tcl = -.2)
  title(xlab = "Iteration $t$",
        ylab = "Overlap", cex.lab = 1, line = 1.1)
  if (legend) {
    lge = sapply(amp, function(x, beta) {
      if (attr(x, 'beta') == 1) {
        paste0('$\\gamma=', beta[attr(x, 'beta')],"$")
      } else if (attr(x, 'robin')) {
        paste0('$\\gamma=', beta[attr(x, 'beta')],"$, Round Robin")
      } else {
        paste0('$\\gamma=', beta[attr(x, 'beta')],"$, Random Update")
      }
    }, beta = beta)
    legend('bottomright',
           legend=lge,
           col = cols,
           lty = 1,
           pch = pts,
           lwd = 1,
           cex = .9,
           bty = 'n')
  }
}

maxline = sqrt(max(amp_bayes_full[[1]]))

tikz(file = './isit-fig/bayes-optimal-comp.tex', width = 4, height = 2.5,documentDeclaration = '\\documentclass{standalone}\n', standAlone = TRUE)
plot_bayes_norm(amp_bayes_norm_full, la, be, TRUE, maxt, maxline)
dev.off()

tikz(file = './isit-fig/power-iteration-comp.tex', width = 4, height = 2.5,documentDeclaration = '\\documentclass{standalone}\n', standAlone = TRUE)
plot_pi_norm(amp_pi_norm_full, la, be, T, maxt, maxline)
dev.off()

tikz(file = './isit-fig/power-iteration-unnorm-iter-comp.tex', width = 4, height = 2.5, documentDeclaration = '\\documentclass{standalone}\n', standAlone = TRUE)
plot_pi(amp_pi_full, la, be, T, 4*maxt, maxline)
dev.off()


### export data as csv
# Bayes-optimal with normalized iterations
bayes_norm_iter = data.frame(iter = 0:maxt)
bayes_norm_iter$se_full = sqrt(amp_bayes_norm_full[[1]])
bayes_norm_iter$amp_full_mean = apply(sqrt(amp_bayes_norm_full[[2]]), 1, mean)
bayes_norm_iter$amp_full_sd = apply(sqrt(amp_bayes_norm_full[[2]]), 1, sd)
bayes_norm_iter$se_robin = sqrt(amp_bayes_norm_full[[3]][10*(0:15) + 1])
bayes_norm_iter$amp_robin_mean = apply(sqrt(amp_bayes_norm_full[[4]][10*(0:15) + 1,]), 1, mean)
bayes_norm_iter$amp_robin_sd = apply(sqrt(amp_bayes_norm_full[[4]][10*(0:15) + 1,]), 1, sd)
bayes_norm_iter$se_rand = sqrt(amp_bayes_norm_full[[5]][10*(0:15) + 1])
bayes_norm_iter$amp_rand_mean = apply(sqrt(amp_bayes_norm_full[[6]][10*(0:15) + 1,]), 1, mean)
bayes_norm_iter$amp_rand_sd = apply(sqrt(amp_bayes_norm_full[[6]][10*(0:15) + 1,]), 1, sd)
write.csv(bayes_norm_iter, "./data/bayes_norm_iter.csv", row.names = FALSE)

# Bayes-optimal with unnormalized iterations
bayes_iter = data.frame(iter = 0:(maxt*4))
bayes_iter$se_full = sqrt(amp_bayes_full[[1]])
bayes_iter$amp_full_mean = apply(sqrt(amp_bayes_full[[2]]), 1, mean)
bayes_iter$amp_full_sd = apply(sqrt(amp_bayes_full[[2]]), 1, sd)
bayes_iter$se_robin = sqrt(amp_bayes_full[[3]])
bayes_iter$amp_robin_mean = apply(sqrt(amp_bayes_full[[4]]), 1, mean)
bayes_iter$amp_robin_sd = apply(sqrt(amp_bayes_full[[4]]), 1, sd)
bayes_iter$se_rand = sqrt(amp_bayes_full[[5]])
bayes_iter$amp_rand_mean = apply(sqrt(amp_bayes_full[[6]]), 1, mean)
bayes_iter$amp_rand_sd = apply(sqrt(amp_bayes_full[[6]]), 1, sd)
write.csv(bayes_iter, "./data/bayes_iter.csv", row.names = FALSE)

# Power iteration with normalized iterations
pi_norm_iter = data.frame(iter = 0:maxt)
pi_norm_iter$se_full = amp_pi_norm_full[[1]]
pi_norm_iter$amp_full_mean = apply(amp_pi_norm_full[[2]], 1, mean)
pi_norm_iter$amp_full_sd = apply(amp_pi_norm_full[[2]], 1, sd)
pi_norm_iter$se_robin = amp_pi_norm_full[[3]][10*(0:15) + 1]
pi_norm_iter$amp_robin_mean = apply(amp_pi_norm_full[[4]][10*(0:15) + 1,], 1, mean)
pi_norm_iter$amp_robin_sd = apply(amp_pi_norm_full[[4]][10*(0:15) + 1,], 1, sd)
pi_norm_iter$se_rand = amp_pi_norm_full[[5]][10*(0:15) + 1]
pi_norm_iter$amp_rand_mean = apply(amp_pi_norm_full[[6]][10*(0:15) + 1,], 1, mean)
pi_norm_iter$amp_rand_sd = apply(amp_pi_norm_full[[6]][10*(0:15) + 1,], 1, sd)
write.csv(pi_norm_iter, "./data/pi_norm_iter.csv", row.names = FALSE)

# Power iteration with unnormalized iterations
pi_iter = data.frame(iter = 0:(maxt*4))
pi_iter$se_full = amp_pi_full[[1]]
pi_iter$amp_full_mean = apply(amp_pi_full[[2]], 1, mean)
pi_iter$amp_full_sd = apply(amp_pi_full[[2]], 1, sd)
pi_iter$se_robin = amp_pi_full[[3]]
pi_iter$amp_robin_mean = apply(amp_pi_full[[4]], 1, mean)
pi_iter$amp_robin_sd = apply(amp_pi_full[[4]], 1, sd)
pi_iter$se_rand = amp_pi_full[[5]]
pi_iter$amp_rand_mean = apply(amp_pi_full[[6]], 1, mean)
pi_iter$amp_rand_sd = apply(amp_pi_full[[6]], 1, sd)
write.csv(pi_iter, "./data/pi_iter.csv", row.names = FALSE)

###################################################################

## Scale-invariance for linear denoisers
scale_uscale_comp = function(lambda, beta, r = 0.01, robin = TRUE) {
  scaled = se_linear_sphere(r, beta, lambda, robin, maxt = 100)[100]
  unscaled = se_linear(r, 1, beta, lambda, robin, maxt = 100)[100]
  return(log(scaled) - log(unscaled))
}

# # Run experiments
# sparsity = c(1,1/2,1/5,1/10)
# lambdagrid = seq(0, 5, length.out = 1000)
# # Round robin
# mat_robin = matrix(0, length(lambdagrid), length(sparsity))
# for (i in 1:length(sparsity)) {
#   mat_robin[,i] = sapply(lambdagrid, function(l,b) scale_uscale_comp(l, b), b = sparsity[i] )
# }
# matplot(lambdagrid, mat_robin, type = 'l')
# # IID diagonal
# mat_rand = matrix(0, length(lambdagrid), length(sparsity))
# for (i in 1:length(sparsity)) {
#   mat_rand[,i] = sapply(lambdagrid, function(l,b) scale_uscale_comp(l, b, robin = F), b = sparsity[i] )
# }
# matplot(lambdagrid, mat_rand, type = 'l')

