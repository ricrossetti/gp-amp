### Numerical experiments
sources = list.files(path = './src', pattern = '*.R', full.names = TRUE)
sapply(sources, source)
rm(sources)

## Positively biased spike, ReLU denoising
# Parameter settings
n = 2000
snr = 1.5
signal = rbinom(n, 1, c(.5,.5))
scalefactor = norm(signal, '2')^2
init = rep(1,n)

# Run AMP
relu_amp_iid = amp_gaus(init, relu_denoiser, signal, snr, maxt = 20)
relu_amp_goe = amp_goe(init, relu_denoiser, relu_denoiser_deriv, signal, snr, maxt = 20)

# Compute SE quantities
relu_ftf_iid = crossprod(relu_amp_iid$est) / n
relu_ftf_goe = crossprod(relu_amp_goe$est) / n
relu_snr_iid = crossprod(relu_amp_iid$est, signal) * snr / scalefactor
relu_snr_goe = crossprod(relu_amp_goe$est, signal) * snr / scalefactor

# Demeaned iterates
relu_x_center_iid = relu_amp_iid$iter - tcrossprod(signal, relu_snr_iid)
relu_x_center_goe = relu_amp_iid$iter - tcrossprod(signal, relu_snr_goe)
relu_xtx_iid = cov2cor(crossprod(relu_x_center_iid) / n)
relu_xtx_goe = cov2cor(crossprod(relu_x_center_goe) / n)

## Positively biased spike, ReLU+memory denoising
# Parameter settings
n = 2000
snr = 1.5
signal = rbinom(n, 1, c(.5,.5))
scalefactor = norm(signal, '2')^2
init = rep(1,n)

# Run AMP
mrelu_amp_iid = amp_gaus(init, memory_relu, signal, snr, maxt = 20)
mrelu_amp_goe = amp_goe(init, memory_relu, memory_relu_deriv, signal, snr, maxt = 20)

# Compute SE quantities
mmrelu_ftf_iid = crossprod(mrelu_amp_iid$est) / n
mrelu_ftf_goe = crossprod(mrelu_amp_goe$est) / n
mrelu_snr_iid = crossprod(mrelu_amp_iid$est, signal) snr / scalefactor
mrelu_snr_goe = crossprod(mrelu_amp_goe$est, signal) snr / scalefactor

# Demeaned iterates
mrelu_x_center_iid = mrelu_amp_iid$iter - tcrossprod(signal, mrelu_snr_iid)
mrelu_x_center_goe = mrelu_amp_iid$iter - tcrossprod(signal, mrelu_snr_goe)
mrelu_xtx_iid = cov2cor(crossprod(mrelu_x_center_iid) / n)
mrelu_xtx_goe = cov2cor(crossprod(mrelu_x_center_goe) / n)
