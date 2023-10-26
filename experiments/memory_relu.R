### Numerical experiments: ReLU + memory
sources = list.files(path = './src', pattern = '*.R', full.names = TRUE)
sapply(sources, source)
rm(sources)

# Parameter settings
n = 3000
snr = 1
signal = rbinom(n, 1, c(.1,.9))
init = rnorm(n, .1)

## AMP recursions
mrelu_amp_iid = amp_gaus(init, memory_relu, signal, snr, maxt = 20)
mrelu_amp_goe = amp_goe(init, memory_relu, memory_relu_deriv, signal, snr, maxt = 20)

tikz(file = './fig/mrelu_iid_pairs.tex', width = 6, height = 6)
pairs_gaus_diagnostics(mrelu_amp_iid$xcenter, mrelu_amp_iid$est)
dev.off()
tikz(file = './fig/mrelu_goe_pairs.tex', width = 6, height = 6)
pairs_gaus_diagnostics(mrelu_amp_goe$xcenter, mrelu_amp_goe$est)
dev.off()


## GIP recursions and SE overlaps
snrvec = round(seq(0.1, 2, length.out = 10), digits = 1)
overlap_list = overlap_simulate(init, memory_relu, signal, snrvec)

tikz(file = './fig/mrelu_overlaps.tex', width = 6, height = 5)
snr_trajectories_plot(snrvec, overlap_list$avg_ov)
dev.off()
