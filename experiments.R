### Numerical experiments
sources = list.files(path = './src', pattern = '*.R', full.names = TRUE)
sapply(sources, source)

## Positively biased spike, ReLU denoising
# Parameter settings
n = 2000
snr = 2
signal = rbinom(n, 1, c(0.5,.5))
init = rep(1,n)

# Run AMP
relu_amp_iid = amp_gaus(init, relu_denoiser, signal, snr, maxt = 20)
relu_amp_goe = amp_goe(init, relu_denoiser, relu_denoiser_deriv, signal, snr, maxt = 20)

# Compute SE quantities
relu_ftf_iid = crossprod(relu_amp_iid$est) / n
relu_ftf_goe = crossprod(relu_amp_goe$est) / n
relu_snr_iid = crossprod(relu_amp_iid$est, signal) / n * snr
relu_snr_goe = crossprod(relu_amp_goe$est, signal) / n * snr

# Demeaned iterates
relu_x_center_iid = relu_amp_iid$iter - tcrossprod(signal, relu_snr_iid)
relu_x_center_goe = relu_amp_iid$iter - tcrossprod(signal, relu_snr_goe)
relu_xtx_iid = cov2cor(crossprod(relu_x_center_iid) / n)
relu_xtx_goe = crossprod(relu_x_center_goe) / n