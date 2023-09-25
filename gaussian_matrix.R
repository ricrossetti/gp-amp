### Pure Gaussian noise: is the unadjusted recursion accurately tracked by a Gaussian?
## Define some reasonable denoisers

# Soft-thresholding
ST = function(y,theta) {
  if (y > theta) 
    return(y - theta)
  else if (y <= -theta)
    return(theta + y)
  else
    return(0)
}

ST_deriv = function(y, theta) {
  if (y > theta) 
    return(1)
  else if (y <= -theta)
    return(-1)
  else
    return(0)
}

# Tanh nonlinearity
TH = function(y, theta) {
  tanh(sqrt(theta) * y)
}

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

## AMP recursions
# Fully Gaussian matrix
amp_gaus = function(rmat, x, denoise, theta, tol=1e-4) {
  if (ncol(rmat) != length(init)) {
    print('Non conformable d8imensions!')
    return(FALSE)
  }
  
  n = length(init)
  improve = tol + 1e2
  t = 0
  maxt = 100
  iter = matrix(NA, n, 0)
  
  while ((improve > tol) & (t < maxt)) {
    t = t+1
    f = sapply(x, denoise, theta = theta)
    x = as.vector(rmat %*% f)
    iter = cbind(iter, x)
    if (t > 1) {
      improve = sum(abs( iter[,ncol(iter)] - iter[,(ncol(iter) - 1)] ))
    }
  }
  
  return(iter)
}
