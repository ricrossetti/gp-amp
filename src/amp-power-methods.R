### Linear AMP and Power Iteration Algos
### All denoisers project on the sphere

## AMP SE
se_robin = function(r,
                    gamma,
                    lambda,
                    maxt) {
  
  eta = function(t, gamma) {
    if (t<1/gamma) {
      c( 1 - gamma * min(t, 1/gamma), rep(gamma, t) )
    } else {
      c( rep(0, t + 1 - 1/gamma), rep(gamma, 1/gamma) )
    }
  }
  for( t in 1:maxt ) {
    w = eta(t-1, gamma)
    r = c(r, lambda*sum(w*r) / sqrt(sum(w * (lambda^2*r^2 + 1)) ))
  }
  return(r)
}

se_rand = function(r,
                   gamma,
                   lambda,
                   maxt) {

  eta = function(t, gamma) (1-gamma)^rev((0:t)) * gamma^((0:t)>0)
    
  for( t in 1:maxt ) {
    w = eta(t-1, gamma)
    r = c(r, lambda*sum(w*r) / sqrt(sum(w * (lambda^2*r^2 + 1)) ))
  }
  return(r)
}

## TRUE CORRECTION AMP
opamp_robin = function(data,
                       signal,
                       init,
                       J,
                       maxt) {
  # Setup initial quantities
  n = length(signal)
  init = init / norm(init, '2') * sqrt(n)
  r = rep(as.numeric(crossprod(signal, init)) / n, maxt + 1)
  est = matrix(init, n, maxt + 1)
  sub_size = floor(n / J)
  
  # Setup first iteration
  x = as.vector( data %*% init )
  
  # Past weight function
  p = function(t, J) {
    if (t < J) {
      c( 1 - 1/J * min(t, J), rep(1/J, t) )
    } else {
      c( rep(0, t + 1 - J), rep(1/J, J) )
    }
  }
  
  # AMP recursion
  for (t in 1:maxt) {
    b = sqrt(n) / (norm(x, '2'))
    S = ((t-1) %% J) * sub_size + (1:sub_size)
    est[,t+1] = b * x
    x[S] = b * as.vector( data[S,] %*% x 
                          - as.matrix(est[S,(1:t)]) %*% p(t-1, J) )
    r[t+1] = crossprod(est[,t+1], signal) / n
    if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
  }
  return(r)
}

opamp_rand = function(data,
                      signal,
                      init,
                      J,
                      maxt) {
  # Setup initial quantities
  n = length(signal)
  init = init / norm(init, '2') * sqrt(n)
  r = rep(as.numeric(crossprod(signal, init)) / n, maxt + 1)
  est = matrix(init, n, maxt + 1)
  sub_size = floor(n / J)
  
  # Setup first iteration
  x = as.vector( data %*% init )
  
  # Past weight function
  p = function(t, J) (1-(1/J))^rev((0:t)) * (1/J)^((0:t)>0)
  
  # AMP recursion
  for (t in 1:maxt) {
    b = sqrt(n) / (norm(x, '2'))
    S = sort(sample( 1:n, sub_size ))
    est[,t+1] = b * x
    x[S] = b * as.vector( data[S,] %*% x 
                          - as.matrix(est[S,(1:t)]) %*% p(t-1, J) )
    r[t+1] = crossprod(est[,t+1], signal) / n
    if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
  }
  return(r)
}

## FULL MATRIX CORRECTION AMP
fullamp_robin = function(data,
                         signal,
                         init,
                         J,
                         maxt) {
  # Setup initial quantities
  n = length(signal)
  init = init / norm(init, '2') * sqrt(n)
  r = rep(as.numeric(crossprod(signal, init)) / n, maxt + 1)
  est = matrix(init, n, maxt + 1)
  sub_size = floor(n / J)
  
  # Setup first iteration
  x = as.vector( data %*% init )
  
  # AMP recursion
  for (t in 1:maxt) {
    b = sqrt(n) / (norm(x, '2'))
    S = ((t-1) %% J) * sub_size + (1:sub_size)
    est[,t+1] = b * x
    x[S] = b * as.vector( data[S,] %*% x - est[S,t])
    r[t+1] = crossprod(est[,t+1], signal) / n
    if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
  }
  return(r)
}

fullamp_rand = function(data,
                        signal,
                        init,
                        J,
                        maxt) {
  # Setup initial quantities
  n = length(signal)
  init = init / norm(init, '2') * sqrt(n)
  r = rep(as.numeric(crossprod(signal, init)) / n, maxt + 1)
  est = matrix(init, n, maxt + 1)
  sub_size = floor(n / J)
  
  # Setup first iteration
  x = as.vector( data %*% init )
  
  # AMP recursion
  for (t in 1:maxt) {
    b = sqrt(n) / (norm(x, '2'))
    S = sort(sample(1:n, sub_size))
    est[,t+1] = b * x
    x[S] = b * as.vector( data[S,] %*% x - est[S,t] )
    r[t+1] = crossprod(est[,t+1], signal) / n
    if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
  }
  return(r)
}

amp_fullmat = function(data,
                       signal,
                       init,
                       J,
                       maxt) {
  fullamp_robin(data,
                signal,
                init,
                1,
                maxt)
}

## POWER METHODS WITHOUT CORRECTIONS
pi_naive_robin = function(data,
                          signal,
                          init,
                          J,
                          maxt) {
  # Setup initial quantities
  n = length(signal)
  init = init / norm(init, '2') * sqrt(n)
  r = rep(as.numeric(crossprod(signal, init)) / n, maxt + 1)
  est = matrix(init, n, maxt + 1)
  sub_size = floor(n / J)
  
  # Setup first iteration
  x = as.vector( data %*% init )
  
  # Recursion
  for (t in 1:maxt) {
    b = sqrt(n) / (norm(x, '2'))
    S = ((t-1) %% J) * sub_size + (1:sub_size)
    est[,t+1] = b * x
    x[S] = b * as.vector( data[S,] %*% x )
    x[-S] = 0
    r[t+1] = crossprod(est[,t+1], signal) / n
    if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
  }
  return(r)
}

pi_naive_rand = function(data,
                         signal,
                         init,
                         J,
                         maxt) {
  # Setup initial quantities
  n = length(signal)
  init = init / norm(init, '2') * sqrt(n)
  r = rep(as.numeric(crossprod(signal, init)) / n, maxt + 1)
  est = matrix(init, n, maxt + 1)
  sub_size = floor(n / J)
  
  # Setup first iteration
  x = as.vector( data %*% init )
  
  # Recursion
  for (t in 1:maxt) {
    b = sqrt(n) / (norm(x, '2'))
    S = sort(sample(1:n, sub_size))
    est[,t+1] = b * x
    x[S] = b * as.vector( data[S,] %*% x )
    x[-S] = 0
    r[t+1] = crossprod(est[,t+1], signal) / n
    if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
  }
  return(r)
}

pi_fullmat = function(data,
                      signal,
                      init,
                      J,
                      maxt) {
  pi_naive_robin(data,
                 signal,
                 init,
                 1,
                 maxt)
}

## POWER METHODS WITH PAST MEMORY
pi_memory_robin = function(data,
                           signal,
                           init,
                           J,
                           maxt) {
  # Setup initial quantities
  n = length(signal)
  init = init / norm(init, '2') * sqrt(n)
  r = rep(as.numeric(crossprod(signal, init)) / n, maxt + 1)
  est = matrix(init, n, maxt + 1)
  sub_size = floor(n / J)
  
  # Setup first iteration
  x = as.vector( data %*% init )
  
  # Recursion
  for (t in 1:maxt) {
    b = sqrt(n) / (norm(x, '2'))
    S = ((t-1) %% J) * sub_size + (1:sub_size)
    est[,t+1] = b * x
    x[S] = b * as.vector( data[S,] %*% x )
    r[t+1] = crossprod(est[,t+1], signal) / n
    if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
  }
  return(r)
}

pi_memory_rand = function(data,
                          signal,
                          init,
                          J,
                          maxt) {
  # Setup initial quantities
  n = length(signal)
  init = init / norm(init, '2') * sqrt(n)
  r = rep(as.numeric(crossprod(signal, init)) / n, maxt + 1)
  est = matrix(init, n, maxt + 1)
  sub_size = floor(n / J)
  
  # Setup first iteration
  x = as.vector( data %*% init )
  
  # Recursion
  for (t in 1:maxt) {
    b = sqrt(n) / (norm(x, '2'))
    S = sort(sample(1:n, sub_size))
    est[,t+1] = b * x
    x[S] = b * as.vector( data[S,] %*% x )
    r[t+1] = crossprod(est[,t+1], signal) / n
    if (t %% 100 == 0) print(paste(t, 'iterations complete.'))
  }
  return(r)
}