### HELPER FUNCTIONS ###
library(Deriv)

## Fast(-ish, this is still R) identity vec-Kronecker multiply
vkmult = function(kron, vec) {
  if ( length(vec) %% length(kron) != 0 ) {
    stop('Non-conformable dimensions!')
  }
  n = length(kron)
  d = length(vec) %/% n
  
  vsplit = split(vec, ceiling(seq_along(vec)/d))
  
  out = mapply(function(x,y) {
    sapply(x, function(a,b) a*b, b = y)
  }, vsplit, kron)
  
  out = as.vector(rowSums(out))
}

## Make noise matrices
# Sample Gaussian matrix
rgaus_mat = function(n, scal = TRUE) {
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

# ReLU
relu = function(x) {
  max(0, x)
}

relu_deriv = function(x) {
  as.numeric(x > 0)
}

# Long-memory ReLU (geometric decay)
dlmr = function(x, memory) {
  if (length(memory) == 0) {
    res = max(0, x) / 2
  } else {
    geom_memory = 2^(-rev(seq_along(memory) + 1)) * memory
    res = max(0, x) / 2 + sum( sapply(geom_memory, function(x) max(0,x) ) )
  }
  return(res)
}

dlmr_grad = function(x, memory) {
  if (length(memory) == 0) {
    res = as.numeric(x > 0) / 2
  } else {
    res = c((memory > 0) * 2^(-rev(seq_along(memory) + 1)), as.numeric(x > 0) / 2) 
  }
  return(res)
}
