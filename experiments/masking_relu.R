### Numerical experiments: ReLU without memory
sources = list.files(path = './src', pattern = '*.R', full.names = TRUE)
sapply(sources, source)
rm(sources)

# Parameter settings
n = 3000
snr = 1
signal = rbinom(n, 1, c(.1,.9))
init = rnorm(n, .1)
er_prob = .75
maxt = 50
resamples = 25

# Initialize quantities
erasures_center_iter = array(NA, dim = c(n,maxt,resamples))
erasures_ov = matrix(NA, maxt, resamples)
full_center_iter = array(NA, dim = c(n,maxt,resamples))
full_ov = matrix(NA, maxt, resamples)

### NUMERICAL EXPERIMENTS
for (r in 1:resamples) {
  data = matrix(rnorm(n^2, sd = 1 /  sqrt(2)), n, n)
  data = snr * outer(signal, signal) / norm(signal, '2') + data + t(data)
  ## AMP recursions
  amp_erasures = amp_goe_erasures(init, 
                                  relu_denoiser, 
                                  relu_denoiser_deriv,
                                  sub_size = ceiling(length(init)*(1-er_prob)),
                                  data = data,
                                  signal = signal,
                                  snr = snr,
                                  maxt = maxt)
  erasures_center_iter[,,r] = amp_erasures$iter - amp_erasures$mean
  erasures_ov[,r] = amp_erasures$overlap
  amp_full = amp_goe(init, 
                     relu_denoiser, 
                     relu_denoiser_deriv, 
                     data = data,
                     signal = signal,
                     snr = snr,
                     maxt = maxt)
  full_center_iter[,,r] = amp_full$xcenter
  full_ov[,r] = amp_full$overlap
  print(paste0(r, " resamplings done"))
}


## Plotting
# Overlap with truth
par(mfrow=c(1,2))
matplot(1:maxt, full_ov, type = 'l', ylim = c(0,1), 
        main = "Full matrix overlap")
matplot(1:maxt, erasures_ov, type = 'l', ylim = c(0,1),
        main = "Masked overlap")
par(mfrow=c(1,1))
ov_gap_mean = apply((full_ov - erasures_ov), 1, mean) 
ov_gap_sd = apply((full_ov - erasures_ov), 1, sd)
plot(1:maxt, ov_gap_mean, 
     type = 'o',
     col = 'blue',
     lty = 3,
     axes = FALSE, 
     frame.plot = TRUE, 
     ylim = c(-.05,max(ov_gap_mean + 0.05)),
     ann = FALSE)
arrows(x0 = 1:maxt,  y0 = ov_gap_mean - ov_gap_sd, y1 = ov_gap_mean + ov_gap_sd,
       angle = 90,
       length = .05,
       code = 3,
       col = 'blue')
abline(h = 0)
axis(1, 1:maxt)
axis(2)
title(main = "Overlap gap from masking", 
      xlab = "Iteration $t$", 
      ylab = "$\\omega_t^{\\mathrm{full}} - \\omega_t$")

matplot(1:maxt, full_ov - erasures_ov, type = 'l', ylim = c(-.1,.2),
        main = "\\omega_t^{\\mathrm{full}} - \\omega_t")
abline(h = 0)

# Self-overlap evolution
full_sd = apply(full_center_iter, c(2,3), sd)
erasures_sd = apply(erasures_center_iter, c(2,3), sd)
# full_so = apply(full_center_iter, c(2,3), function(x) sqrt(sum(x^2)/length(x)) )
# erasures_so = apply(erasures_center_iter, c(2,3), function(x) sqrt(sum(x^2)/length(x)) )
par(mfrow=c(1,2))
matplot(1:maxt, full_sd,
        type = 'l',
        ylim = c(.9,1.1),
        axes = FALSE,
        frame.plot = TRUE,
        ann = FALSE)
axis(1, 1:maxt)
axis(2)
title(main = "Full matrix", 
      xlab = "Iteration $t$", 
      ylab = "$\\| \\bx_t \\|$")
matplot(1:maxt, erasures_sd, 
        type = 'l', 
        ylim = c(.9,1.1),
        axes = FALSE,
        frame.plot = TRUE,
        ann = FALSE)
axis(1, 1:maxt)
axis(2)
title(main = "Masked matrix", 
      xlab = "Iteration $t$")

# Centering diagnostics
full_mean = apply(full_center_iter, c(2,3), mean)
erasures_mean = apply(erasures_center_iter, c(2,3), mean)
par(mfrow=c(1,2))
matplot(1:maxt, full_mean,
        type = 'l',
        ylim = c(-.1,.1),
        axes = FALSE,
        frame.plot = TRUE,
        ann = FALSE)
axis(1, 1:maxt)
axis(2)
title(main = "Full matrix", 
      xlab = "Iteration $t$", 
      ylab = "$m_t$")
matplot(1:maxt, erasures_mean, 
        type = 'l', 
        ylim = c(-.1,.1),
        axes = FALSE,
        frame.plot = TRUE,
        ann = FALSE)
axis(1, 1:maxt)
axis(2)
title(main = "Masked matrix", 
      xlab = "Iteration $t$")

# ## AMP recursions
# amp_erasures = amp_goe_erasures(init, 
#                                 relu_denoiser, 
#                                 relu_denoiser_deriv,
#                                 sub_size = ceiling(length(init) * visible_frac),
#                                 signal = signal,
#                                 maxt = maxt,
#                                 seed = seed)
# # Center iterates
# amp_erasures_centered = amp_erasures$iter - amp_erasures$mean
# #self-overlaps and mean
# erasures_so_mean = apply(amp_erasures_centered, 2, mean)
# erasures_so_sd = apply(amp_erasures_centered, 2, sd)
# 
# amp_oneshot = amp_goe(init, 
#                       relu_denoiser, 
#                       relu_denoiser_deriv, 
#                       signal = signal,
#                       maxt = maxt,
#                       seed = seed)
# #self-overlaps and mean
# oneshot_so_mean = apply(amp_oneshot$xcenter, 2, mean)
# oneshot_so_sd = apply(amp_oneshot$xcenter, 2, sd)

# ## Plot stuff
# # estimate overlap with truth
# plot(1:maxt, amp_erasures$overlap, type = 'l', main = "Overlap with signal", ylim = c(0,1))
# lines(amp_oneshot$overlap, col = 'red',)
# 
# # iterate self-overlap (should be always about 1)
# plot(1:maxt, erasures_so_sd, type = 'l', main = "Centered iterate self-overlap", ylim = c(.9,1.1))
# lines(oneshot_so_sd, col = 'red')
# 
# # centered iterate mean (should be always about 0)
# plot(1:maxt, erasures_so_mean, type = 'l', main = "Centered iterate mean", ylim = c(-.1,.1))
# lines(oneshot_so_mean, col = 'red')
 