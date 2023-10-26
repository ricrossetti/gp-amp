### Numerical experiments: ReLU without memory
sources = list.files(path = './src', pattern = '*.R', full.names = TRUE)
sapply(sources, source)
rm(sources)

# Parameter settings
n = 3000
snr = 1
signal = rbinom(n, 1, c(.1,.9))
init = rnorm(n, .1)

## AMP recursions
relu_amp_iid = amp_gaus(init, relu_denoiser, signal, snr, maxt = 20)
relu_amp_goe = amp_goe(init, relu_denoiser, relu_denoiser_deriv, signal, snr, maxt = 20)

tikz(file = './fig/relu_iid_pairs.tex', width = 6, height = 6)
pairs_gaus_diagnostics(relu_amp_iid$xcenter, relu_amp_iid$est)
dev.off()
tikz(file = './fig/relu_goe_pairs.tex', width = 6, height = 6)
pairs_gaus_diagnostics(relu_amp_goe$xcenter, relu_amp_goe$est)
dev.off()

## GIP recursions and SE overlaps
snrvec = round(seq(0.1, 2, length.out = 10), digits = 1)
overlap_list = overlap_simulate(init, relu_denoiser, signal, snrvec)

tikz(file = './fig/relu_overlaps.tex', width = 6, height = 5)
snr_trajectories_plot(snrvec, overlap_list$avg_ov)
dev.off()
