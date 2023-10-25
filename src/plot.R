### Plotting functions
library(latex2exp)
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
panel.SE <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par("usr" = usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(as.numeric(crossprod(x,y))/(norm(x, '2')*norm(y, '2')))
  txt <- format(r, digits = digits, nsmall = digits)
  txt <- paste0(prefix, txt)
  if(missing(cex.cor)) cex.cor <- 0.5/strwidth(txt)
  text(0.5, 0.5, txt, cex = cex.cor)
}

panel.scatter = function(x, y, digits = 2, prefix = "cor = ", ...) {
  points(x, y)
  r = cor(x, y)
  txt = format(r, digits = digits, nsmall = digits)
  txt = paste0(prefix, txt)
  coords = par("usr")
  text( coords[2]-0.05*(coords[2]-coords[1]), coords[3]+0.05*(coords[4]-coords[3]), txt, adj = c(1,0))
}

panel.hist = function(x, ...) {
  usr = par("usr"); on.exit(par("usr" = usr))
  rho = sd(x)
  h = hist(x, plot = FALSE)
  breaks = h$breaks; nB = length(breaks)
  y = h$density
  par(usr = c(breaks[1], breaks[length(breaks)], 0, 1.5 * max(y)))
  rect(breaks[-nB], 0, breaks[-1], y, ...)
  curve(dnorm(x, 0, rho), breaks[1], breaks[length(breaks)], add = TRUE)
}

pairs_gaus_diagnostics = function(data, subset.spacing = 4) {
  n = ncol(data)
  s = seq(1, n, subset.spacing)
  data.sub = data[,s]
  pairs(data.sub, labels = paste0('x',s), lower.panel = panel.scatter, upper.panel = NULL, diag.panel = panel.hist, xaxt = 'n', yaxt = 'n')
}
