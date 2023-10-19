### Numerical experiments
sources = list.files(path = './src', pattern = '*.R', full.names = TRUE)
sapply(sources, source)

## Paramter settings
n = 1000
snr = 2

### CASE 1: No signal, no memory, tanh denoiser
denoise = dtanh
deriv = Deriv::Deriv(dtanh, 'x')
init = rnorm(n)

# Gaussian and GOE AMP
case1_gaus = gp_amp_gaus(init, denoise, maxt = 20, theta = 1)
case1_goe = gp_amp_goe(init, denoise, deriv, maxt = 20, theta = 1)
# Compute SE estimates and empirical covariances
case1_gaus_se = lapply(case1_gaus, cov )
case1_goe_se = lapply(case1_goe, cov )
par(mfrow = c(1,2))
heat_diag(cov2cor(case1_gaus_se$est), diag(case1_gaus_se$est), main = "Gaussian covariance")
heat_diag(cov2cor(case1_goe_se$est), diag(case1_goe_se$est), main = "GOE covariance")
par(mfrow = c(1,1))

# Monte Carlo estiates of SE 
case1_mc_se = se_mc(denoise, init, snr, size = 20, theta = 1)

# Do SE estimates and empirical covariance of iterates match?
case1_scale_diff_gaus = cov_scale_discrepancy(case1_gaus$iter, case1_gaus$est)
case1_cor_diff_gaus = cor_discrepancy(case1_gaus$iter, case1_gaus$est)
heat_diag(case1_cor_diff_gaus, case1_scale_diff_gaus, main = "Covariance Discrepancies", sub = "No signal, tanh, Gaussian noise")

case1_scale_diff_goe = cov_scale_discrepancy(case1_goe$iter, case1_goe$est)
case1_cor_diff_goe = cor_discrepancy(case1_goe$iter, case1_goe$est)
heat_diag(case1_cor_diff_goe, case1_scale_diff_goe, main = "Covariance Discrepancies", sub = "No signal, tanh, GOE noise")

# Are the two SE estimates comparable?
case1_scale_compare = cov_scale_discrepancy(case1_gaus$est, case1_goe$est)
case1_cor_compare = cor_discrepancy(case1_gaus$est, case1_goe$est)
heat_diag(case1_cor_compare, case1_scale_compare, main = "GOE vs Gaussian comparison", sub = "No signal, tanh")
dev.off()

### CASE 2: No signal, no memory, ReLU denoiser
denoise = relu
deriv = relu_deriv
init = rnorm(n)

# Gaussian and GOE AMP
case2_gaus = gp_amp_gaus(init, denoise, maxt = 20)
case2_goe = gp_amp_goe(init, denoise, deriv, maxt = 20)
# Compute SE estimates and empirical covariances
case2_gaus_se = lapply(case2_gaus, cov )
case2_goe_se = lapply(case2_goe, cov )
par(mfrow = c(1,2))
heat_diag(cov2cor(case2_gaus_se$est), diag(case2_gaus_se$est), main = "Gaussian covariance")
heat_diag(cov2cor(case2_goe_se$est), diag(case2_goe_se$est), main = "GOE covariance")
par(mfrow = c(1,1))
# Monte Carlo estiates of SE 
case2_mc_se = se_mc(denoise, init, snr, size = 20)

# Do SE estimates and empirical covariance of iterates match?
case2_scale_diff_gaus = cov_scale_discrepancy(case2_gaus$iter, case2_gaus$est)
case2_cor_diff_gaus = cor_discrepancy(case2_gaus$iter, case2_gaus$est)
heat_diag(case2_cor_diff_gaus, case2_scale_diff_gaus, main = "Covariance Discrepancies", sub = "No signal, ReLU, Gaussian noise")

case2_scale_diff_goe = cov_scale_discrepancy(case2_goe$iter, case2_goe$est)
case2_cor_diff_goe = cor_discrepancy(case2_goe$iter, case2_goe$est)
heat_diag(case2_cor_diff_goe, case2_scale_diff_goe, main = "Covariance Discrepancies", sub = "No signal, ReLU, GOE noise")

# Are the two SE estimates comparable?
case2_scale_compare = cov_scale_discrepancy(case2_gaus$est, case2_goe$est)
case2_cor_compare = cor_discrepancy(case2_gaus$est, case2_goe$est)
heat_diag(case2_cor_compare, case2_scale_compare, main = "GOE vs Gaussian comparison", sub = "No signal, ReLU")
dev.off()

### CASE 3: No Signal, denoiser with memory (geometric sum of ReLUs)
denoise = dlmr
deriv = dlmr_grad
init = rnorm(n)

# Gaussian and GOE AMP
case3_gaus = gp_amp_gaus(init, denoise, maxt = 20, memory = TRUE)
case3_goe = gp_amp_goe(init, denoise, deriv, maxt = 20, memory = TRUE)
# Compute SE estimates and empirical covariances
case3_gaus_se = lapply(case3_gaus, cov )
case3_goe_se = lapply(case3_goe, cov )
par(mfrow = c(1,2))
heat_diag(cov2cor(case3_gaus_se$est), diag(case3_gaus_se$est), main = "Gaussian covariance")
heat_diag(cov2cor(case3_goe_se$est), diag(case3_goe_se$est), main = "GOE covariance")
par(mfrow = c(1,1))
# Monte Carlo estiates of SE 
case3_mc_se = se_mc(denoise, init, snr, size = 20, memory = TRUE)

# Do SE estimates and empirical covariance of iterates match?
case3_scale_diff_gaus = cov_scale_discrepancy(case3_gaus$iter, case3_gaus$est)
case3_cor_diff_gaus = cor_discrepancy(case3_gaus$iter, case3_gaus$est)
heat_diag(case3_cor_diff_gaus, case3_scale_diff_gaus, main = "Covariance Discrepancies", sub = "No signal, ReLU + memory, Gaussian noise")

case3_scale_diff_goe = cov_scale_discrepancy(case3_goe$iter, case3_goe$est)
case3_cor_diff_goe = cor_discrepancy(case3_goe$iter, case3_goe$est)
heat_diag(case3_cor_diff_goe, case3_scale_diff_goe, main = "Covariance Discrepancies", sub = "No signal, ReLU + memory, GOE noise")

# Are the two SE estimates comparable?
case3_scale_compare = cov_scale_discrepancy(case3_gaus$est, case3_goe$est)
case3_cor_compare = cor_discrepancy(case3_gaus$est, case3_goe$est)
heat_diag(case3_cor_compare, case3_scale_compare, main = "GOE vs Gaussian comparison", sub = "No signal, ReLU + memory")
dev.off()
