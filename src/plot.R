### Plotting functions
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
