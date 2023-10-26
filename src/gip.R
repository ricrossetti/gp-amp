### Gaussian innovations processes
gip_iid_sampler = function(f0, denoise, signal = FALSE, snr = 1, maxt = 50, tol = 10e-6, ...) {
  # Initialize with zeroth iteration
  n = length(f0)
  bigF = matrix(f0 / norm(f0, '2'), n, 1)
  if (!isFALSE(signal)) {
    signal_scale = signal / norm(signal, '2')
    overlap = c(crossprod(signal_scale, f0)) / norm(f0, '2')
  } else {
    overlap = FALSE
  }
  Sig_inv = 1
  bigX = matrix(rnorm(n), n, 1)
  
  for (t in 1:maxt) {
    # Compute next innovation
    f = apply(bigX + if(!isFALSE(signal)) {
      snr * tcrossprod(signal, overlap)
    } else {
        0
    },
    1,
    denoise, 
    ...)
    normf = norm(f, '2')
    f = f / normf
    
    rho = crossprod(bigF, f)
    schur = 1 - c(crossprod(rho, Sig_inv %*% rho))
    condmean = bigX %*% Sig_inv %*% rho
    
    innovation = rnorm(n, sd = sqrt(schur))
    x = condmean + innovation
    
    # Update recursive quantities
    if (schur < tol) {
      print(paste0('Covariance collapsed to singular after ',t,' iterations.'))
      return(list(iter = bigX, est = bigF, ov = overlap))
    }
    bigF = cbind(bigF, f); bigX = cbind(bigX, x)
    if (!isFALSE(signal)) { 
      overlap = c(overlap, c(crossprod(signal_scale, f))) 
    }
    sig_add = matrix(0, t+1, t+1)
    sig_add[(1:t), (1:t)] = tcrossprod(Sig_inv %*% rho) / schur + Sig_inv
    sig_add[(1:t), (t+1)] = - Sig_inv %*% rho / schur
    sig_add[(t+1), (1:t)] = - t(Sig_inv %*% rho) / schur
    sig_add[(t+1), (t+1)] = 1 / schur
    Sig_inv = sig_add
  }
  return(list(iter = bigX, est = bigF, ov = overlap))
}

gip_goe_sampler = function(f0, denoise, deriv, maxt = 50, ...) {
  # Initialize with zeroth iteration
  n = length(f0)
  bigF = matrix(f0 / norm(f0, '2'), n, 1)
  
  bigPsi = matrix(0, n, 1)
  
  Sig_inv = 1
  bigX = matrix(rnorm(n), n, 1)
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

### Related functions
overlap_simulate = function(f0, denoiser, signal, snrvec, mcruns = 25, ...) {
  mean_overlap_list = list()
  sd_overlap_list = list()
  for (i in 1:length(snrvec)) {
    ovlist = list()
    for (m in 1:mcruns) {
      ov = gip_iid_sampler(f0, denoiser, signal, snr = snrvec[i], ...)$ov
      ovlist = c(ovlist, list(ov))
    }
    mint = min(sapply(ovlist, length))
    ovlist_cut = sapply(ovlist, function(x,mint) x[1:mint], mint = mint)
    ov_mean = rowMeans(ovlist_cut)
    ov_sd = apply(ovlist_cut, 1, sd)
    mean_overlap_list = c(mean_overlap_list, list(ov_mean))
    sd_overlap_list = c(sd_overlap_list, list(ov_sd))
  }
  return(list(avg_ov = mean_overlap_list, sd_ov = sd_overlap_list))
}