### Plotting functions
library(latex2exp)
library(tikzDevice)
# Heatmap with annotated diagonals
heat_diag = function(m, diags,...) {

  l = ncol(m)
  if (l != length(diags)) {
    stop('Scale vector and matrix have incompatible dimensions')
  }
  
  # Create a heatmap without annotations
  fields::image.plot( 1:l, 1:l, m[,rev(seq_len(l))],
         col = colorRampPalette(terrain.colors(4, rev = TRUE))(20),
         xlab = NA,
         ylab = NA,
         axes = FALSE,
         ...)
  
  # Add annotations to the diagonal
  labs = signif(diags, 1)
  text(seq_len(l), rev(seq_len(l)), labels = labs, col = "black", cex = .8)
  
}

# Pairs plot for assessing normality and SE
pairs_gaus_diagnostics = function(iter, est, title = NULL, nfacets = 5, ...) {
  n = ncol(iter)
  l = nrow(iter)
  s = floor(seq(1, n, length.out = nfacets))
  iter.sub = iter[,s]
  est.sub = est[,s]
  data = rbind(iter.sub, est.sub)
  pairs(data, 
        labels = paste0('$x_{',s,'}$'), 
        lower.panel = function(x, y, digits = 2, prefix = "cor = ", ...) {
          subsample = sample(1:l, min(l,100))
          x = x[seq_len(l)][subsample]
          y = y[seq_len(l)][subsample]
          usr = par("usr"); on.exit(par("usr" = usr))
          par(usr = c(-4, 4, -4, 4))
          points(x/sd(x), y/sd(y), xaxp = c(-4,4,7), ...)
          r = cor(x, y)
          txt = format(r, digits = digits, nsmall = digits)
          txt = paste0(prefix, txt)
          coords = par("usr")
          text(coords[2]-0.05*(coords[2]-coords[1]), 
               coords[3]+0.05*(coords[4]-coords[3]), 
               txt,
               adj = c(1,0))
          if (par("mfg")[2]==1) axis(2, outer = TRUE)
        }, 
        upper.panel = function(x, y, digits = 2, prefix = "", cex.cor = 2, ...) {
          x = x[(l+1):length(x)]
          y = y[(l+1):length(y)]
          usr.curr <- par("usr"); on.exit(par("usr" = usr.curr))
          par(usr = c(0, 1, 0, 1))
          r <- abs(as.numeric(crossprod(x,y))/(norm(x, '2')*norm(y, '2')))
          txt <- format(r, digits = digits, nsmall = digits)
          txt <- paste0("$\\langle f_s, f_t\\rangle$ \n", txt)
          text(0.5, 0.5, txt, cex = cex.cor)
        }, 
        diag.panel = function(x, ...) {
          x = x[seq_len(l)]
          usr = par("usr"); on.exit(par("usr" = usr))
          h = hist(x / sd(x), plot = FALSE)
          breaks = h$breaks; nB = length(breaks)
          y = h$density
          par(usr = c(-4, 4, 0, 1.5 * max(y)))
          rect(breaks[-nB], 0, breaks[-1], y, xaxp = c(-4,4,7), ...)
          curve(dnorm, breaks[1], breaks[length(breaks)], add = TRUE)
          axis(1, outer = TRUE)
        }, 
        xaxt = 'n', yaxt = 'n', oma = c(2,2,2,2), ...)
  if (!is.null(title)) title(title, line = 3)
}

# Plot effective SNR trajectories
snr_trajectories_plot = function(snr, overlaps, ...) {
  n = length(snr)
  if (length(overlaps) != n) {stop("SNR vector and overlap lists don't match!")}
  maxt = max(sapply(overlaps, length)) - 1
  ubound = max(sapply(overlaps, max))
  cols = colorRampPalette(c('blue','orange'))(n)
  par(mar = c(2,2,2,10), xpd = TRUE)
  plot(1:length(overlaps[[1]]) - 1, overlaps[[1]], type = 'l', col = cols[1], 
       xaxt = 'n', yaxt = 'n', xlim = c(0, maxt), ylim =  c(0, 1.05 * ubound ),
       xaxs = 'i', yaxs = 'i', ...)
  for (i in 2:n) {
    lines(1:length(overlaps[[i]]) - 1, overlaps[[i]], col = cols[i])
  }
  legend('right', inset = c(-0.3,0), legend = paste0('$',snr,'$'), col = cols, 
         lty = 1, title = "SNR $\\lambda$", xjust = 1)
  axis(1, at = floor(seq(0, maxt, length.out = 10)))
  axis(2)
}

# Plot overlap trajectory across iteration
overlap_plot = function(est, signal, ...) {
  iter = 1:ncol(est)
  scales = apply(est, 2, norm, type = '2')
  overlap = as.vector(crossprod(est, signal)) / scales / norm(signal, '2')
  plot(iter, overlap, type = 'l', ...)
}

### UNSUSED FUNCTIONS
# panel.SE <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
#   usr <- par("usr"); on.exit(par("usr" = usr))
#   par(usr = c(0, 1, 0, 1))
#   r <- abs(as.numeric(crossprod(x,y))/(norm(x, '2')*norm(y, '2')))
#   txt <- format(r, digits = digits, nsmall = digits)
#   txt <- paste0(prefix, txt)
#   if(missing(cex.cor)) cex.cor <- 0.5/strwidth(txt)
#   text(0.5, 0.5, txt, cex = cex.cor)
# }
# 
# panel.scatter = function(x, y, digits = 2, prefix = "cor = ", ...) {
#   points(x/sd(x), y/sd(y))
#   r = cor(x, y)
#   txt = format(r, digits = digits, nsmall = digits)
#   txt = paste0(prefix, txt)
#   coords = par("usr")
#   text(coords[2]-0.05*(coords[2]-coords[1]),
#        coords[3]+0.05*(coords[4]-coords[3]),
#        txt,
#        adj = c(1,0))
# }
# 
# panel.hist = function(x, ...) {
#   usr = par("usr"); on.exit(par("usr" = usr))
#   h = hist(x / sd(x), plot = FALSE)
#   breaks = h$breaks; nB = length(breaks)
#   y = h$density
#   par(usr = c(-4, 4, 0, 1.5 * max(y)))
#   rect(breaks[-nB], 0, breaks[-1], y, xaxp = c(-4,4,7), ...)
#   curve(dnorm, breaks[1], breaks[length(breaks)], add = TRUE)
#   axis(1, outer = TRUE)
# }

