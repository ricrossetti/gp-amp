source("./src/amp-power-methods.R")

## Parameter settings
set.seed(27)
n = 5000
la = sqrt(2)
ga = 0.1
maxt = 60
resample = 50
signal = sample(c(-1,1), n, TRUE)
init = tanh( 0.1 * (0.1 * signal + rnorm(n)) )
r0 = crossprod(signal, init) / (sqrt(n) * norm(init, '2'))

## Function to evaluate recursions
mc_trials = function( f, sig, init, J, lambda, maxt, resample ) {
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

mc_trials_iidnoise = function( f, sig, init, J, lambda, maxt, resample ) {
  res = matrix(0, maxt + 1, resample)
  for (t in 1:resample) {
    set.seed(t)
    data = matrix(rnorm(n^2, sd = 1/sqrt(n)), n)
    data = lambda / n * tcrossprod(signal) + data
    amp = f( data, sig, init, J, maxt )
    res[,t] = amp
    if (t %% 5 == 0) print(paste(t, 'resamplings complete.'))
  }
  return(res)
}

mc_trials_sgnoise = function( f, sig, init, J, lambda, maxt, resample ) {
  res = matrix(0, maxt + 1, resample)
  for (t in 1:resample) {
    set.seed(t)
    noise = sample(c(-1,1), n*(n-1)/2, replace = TRUE)
    data = matrix(0, n, n); data[lower.tri(data)] = noise
    data = lambda / n * tcrossprod(signal) + (data + t(data) + diag(1, n)) / sqrt(n)
    amp = f( data, sig, init, J, maxt )
    res[,t] = amp
    if (t %% 5 == 0) print(paste(t, 'resamplings complete.'))
  }
  return(res)
}

## AMP SE predictions
se_table = data.frame( iter = 0:maxt,
                       se_fullmat = se_robin(r0, 1, la, maxt),
                       se_robin = se_robin(r0, ga, la, maxt),
                       se_rand = se_rand(r0, ga, la, maxt) )

## Build AMP data
method_list = list(pi_fullmat, 
                   pi_naive_rand, 
                   pi_naive_robin,
                   pi_memory_rand,
                   pi_memory_robin,
                   amp_fullmat,
                   fullamp_rand,
                   fullamp_robin,
                   opamp_rand,
                   opamp_robin)
res_list = lapply(method_list, mc_trials, 
                  sig = signal,
                  init = init,
                  J = 1/ga,
                  lambda = la,
                  maxt = maxt,
                  resample = resample)
saveRDS(res_list, file = 'banff-gaus.rds')

amp_table = data.frame( lapply(res_list, function(x) {
  apply(x, 1, median)
}) )
names(amp_table) = c('pi_fullmat', 
                     'pi_naive_rand', 
                     'pi_naive_robin',
                     'pi_memory_rand',
                     'pi_memory_robin',
                     'amp_fullmat',
                     'fullamp_rand',
                     'fullamp_robin',
                     'opamp_rand',
                     'opamp_robin')
gaus_data_table = cbind(se_table, amp_table)
write.csv(gaus_data_table, file = './data/banff-gaus.csv', row.names = FALSE)
gaus_data_table_sub3 = gaus_data_table[(0:floor(nrow(gaus_data_table)/3))*3 + 1,]
write.csv(gaus_data_table_sub3, file = './data/banff-gaus-sub3.csv', row.names = FALSE)

## AMP with iid matrices
method_list_iid = list(pi_fullmat,
                       pi_memory_rand,
                       pi_memory_robin)
res_list_iid = lapply(method_list_iid, mc_trials_iidnoise, 
                      sig = signal,
                      init = init,
                      J = 1/ga,
                      lambda = la,
                      maxt = maxt,
                      resample = resample)
saveRDS(res_list_iid, file = 'banff-gaus-iid.rds')
amp_table_iid = data.frame( lapply(res_list_iid, function(x) {
  apply(x, 1, median)
}) )
names(amp_table_iid) = c('amp_fullmat',
                         'opamp_rand',
                         'opamp_robin')
iid_data_table = cbind(se_table, amp_table_iid)
write.csv(iid_data_table, file = './data/banff-iid.csv', row.names = FALSE)
iid_data_table_sub3 = iid_data_table[(0:floor(nrow(iid_data_table)/3))*3 + 1,]
write.csv(iid_data_table_sub3, file = './data/banff-iid-sub3.csv', row.names = FALSE)

## AMP with iid matrices
method_list_sg = list(amp_fullmat,
                      opamp_rand,
                      opamp_robin)
res_list_sg = lapply(method_list_sg, mc_trials_sgnoise, 
                      sig = signal,
                      init = init,
                      J = 1/ga,
                      lambda = la,
                      maxt = maxt,
                      resample = resample)
saveRDS(res_list_sg, file = 'banff-gaus-sg.rds')
amp_table_sg = data.frame( lapply(res_list_sg, function(x) {
  apply(x, 1, median)
}) )
names(amp_table_sg) = c('amp_fullmat',
                         'opamp_rand',
                         'opamp_robin')
sg_data_table = cbind(se_table, amp_table_sg)
write.csv(sg_data_table, file = './data/banff-sg.csv', row.names = FALSE)
sg_data_table_sub3 = sg_data_table[(0:floor(nrow(sg_data_table)/3))*3 + 1,]
write.csv(sg_data_table_sub3, file = './data/banff-sg-sub3.csv', row.names = FALSE)