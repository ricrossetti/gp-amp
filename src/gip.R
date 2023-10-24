### Gaussian innovations processes
gip_goe_sampler = function(f0, denoise, deriv, maxt = 50, ...) {
  # Initialize with zeroth iteration
  n = length(f0)
  bigF = matrix(f0 / norm(f0, '2'), n, 1)
  
  bigPsi = matrix(0, n, 1)
  
  Sig_inv = 1
  bigX = matrix(rnorm(n), n)
  FtX = crossprod(bigF, bigX)
  
  for (t in 1:maxt) {
    # Compute conditional mean function
    f = apply(bigX, 1, denoise, ...)
    normf = norm(f, '2')
    f = f / normf
    b = matrix(apply(bigX, 1, deriv, ...), ncol = t, byrow = TRUE)
    b = colSums(b) / normf
    psi = bigF %*% b
    rho = crossprod(bigF, f)
    s = crossprod(bigX, f)
    crosscov = Sig_inv %*% rho
    
    cond_exp = bigX %*% crosscov + bigF %*% ( Sig_inv %*% s - psi ) + ( bigF %*% Sig_inv %*% FtX - bigPsi ) %*% crosscov
    
    # Sample innovation Gaussian
    z1 = rnorm(n)
    z2 = rnorm(ncol(bigF))
    z3 = rnorm(1)
    
    r = crossprod(rho, crosscov)
    eigdec = eigen(Sig_inv, symmetric = TRUE)
    invroot = eigdec$vectors %*% sqrt(eigdec$values) %*% t(eigdec$vectors)
    projhalf = bigF %*% invroot
    
    rankone = f - bigF %*% crosscov
    
    innovation = sqrt(r) * (z1 - projhalf %*% z2) + rankone * z3
    x = cond_exp + innovation
    
    # Update quantities
    bigF = cbind(bigF, f)
    bigX = cbind(bigX, x)
    schur = 1 - r
    Sig_inv = Sig_inv + tcrossprod(crosscov) / schur
    Sig_inv = cbind(Sig_inv, - crosscov / schur) 
    Sig_inv = rbind(Sig_inv, c(- crosscov / schur, 1 / schur))
    bigPsi = cbind(bigPsi, psi)
  }
  
  return(bigX)
} 