### HELPER FUNCTIONS ###

## Sample GOE
rgoe = function( n, scale = TRUE ) {
  out = matrix(0, n, n)
  vec = rnorm( (n+1)*n / 2, sd = n^(scale * -1/2))
  out[lower.tri(out, diag = TRUE)] = vec
  out = out + t(out)
  return(out)
}
