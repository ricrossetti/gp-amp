### HELPER FUNCTIONS ###
library(Deriv)

## Make noise matrices
# Sample Gaussian matrix
rgaus_mat = function(n, spherical = TRUE) {
  matrix(rnorm(n^2, sd = n^(spherical * -1/2)), n, n)
}

# Sample GOE
rgoe = function( n, scale = TRUE ) {
  out = matrix(0, n, n)
  vec = rnorm( (n+1)*n / 2, sd = n^(scale * -1/2))
  out[lower.tri(out, diag = TRUE)] = vec
  out = out + t(out)
  return(out)
}

## Some spike distributions (renormalized for second moment)
# Bernoulli spike
bspike = function( n, p, scale = TRUE) {
  v = rbinom(n, 1, p) / sqrt(p)
  if (scale) spike = tcrossprod(v) / n else spike = tcrossprod(v) / (sqrt(n))
  return(list(v, spike))
}

# Rademacher spike
rspike = function( n, scale = TRUE ) {
  v = sample(c(-1,1), n, TRUE, c(0.5,0.5))
  if (scale) spike = tcrossprod(v) / n else spike = tcrossprod(v) / (sqrt(n))
  return(list(v, spike))
}

# Sparse Rademacher spike
srspike = function( n, p, scale = TRUE) {
  v = sample(c(-1,0,1), n, TRUE, c(p/2, 1-p, p/2)) / sqrt(p)
  if (scale) spike = tcrossprod(v) / n else spike = tcrossprod(v) / (sqrt(n))
  return(list(v, spike))
}

# Gaussian spike
gspike = function(n, scale = TRUE) {
  v = rnorm(n)
  if (scale) spike = tcrossprod(v) / n else spike = tcrossprod(v) / (sqrt(n))
  return(list(v, spike))
}

# Bernoulli-Gaussian spike
bgspike = function(n, p, scale = TRUE) {
  v = rbinom(n, 1, p) * rnorm(n, sd = 1/sqrt(p))
  if (scale) spike = tcrossprod(v) / n else spike = tcrossprod(v) / (sqrt(n))
  return(list(v, spike))
}

## Some denoiser funcions
# Soft-thresholding
dst = function(x, theta) {
    if (x > theta) {
      x - theta 
    } else if (x <= -theta) {
      x + theta
    } else {
      0
    }
}

# Tanh (optimal for Rademacher)
dtanh = function(x, theta) {
  tanh( sqrt(theta) * x )
}

# Geometric-memory ReLU

