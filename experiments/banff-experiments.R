source("./src/amp-power-methods.R")

## Parameter settings
set.seed(27)
n = 5000
la = sqrt(2)
ga = 0.1
maxt = 100
resample = 100
signal = sample(c(-1,1), n, TRUE)
init = tanh( 0.1 * (0.1 * signal + rnorm(n)) )
r0 = crossprod(signal, init) / (sqrt(n) * norm(init, '2'))

## Function to evaluate recursions
mc_trials = function( resample, f, sig, init, J, lambda, maxt ) {
  res = matrix(0, maxt + 1, resample)
  for (t in 1:resample) {
    set.seed(t)
    data = matrix(rnorm(n^2, sd = 1/sqrt(2)), n)
    data = lambda / n * tcrossprod(signal) + (data + t(data)) / sqrt(n)
    amp = f( data, sig, init, J, maxt )
    res[,t] = amp
    if (t %% 5 == 0) print(paste(t, 'resamplings complete.'))
  }
  return(res)
}

## AMP SE predictions
datatable = data.frame( iter = 0:maxt,
                        fullmat = se_robin(r0, 1, la, maxt),
                        robin = se_robin(r0, ga, la, maxt),
                        rand = se_rand(r0, ga, la, maxt) )

## Build AMP data


