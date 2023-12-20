# ReLU denoiser
relu_denoiser = function(x, scale = 1) {
  out = sapply( x, function(x, scale) {max(0, scale*x)}, scale = scale )
  return(out)
}

relu_denoiser_deriv = function(x, scale = 1) {
  out = sapply( x, function(x, scale) {max(0, sign(x)) * scale}, scale = scale)
  return(out)
}

# ReLU with memory
memory_relu = function(x, memory = Inf, rate = 2) {
  input_size = length(x)
  
  arg = rev(x[ max(1, input_size - memory + 1):input_size ])
  decay = rate^(-(1:length(arg)))
  
  arg_relu = sapply(arg, function(x) max(0,x))
  out = as.numeric(crossprod(decay, arg_relu))
  return(out)
}

memory_relu_deriv = function(x, memory = Inf, rate = 2) {
  
  input_size = length(x)
  
  arg = x[ max(1, input_size - memory + 1):input_size ]
  decay = rev(rate^(-(1:length(arg))))
  
  out = decay * sapply(arg, function(x) max(0,sign(x)))
  return(out)
}

# Tanh function
tanh_denoiser = function(x, scale = 1) {
  out = tanh(sqrt(scale) * x)
  return(out)
}

tanh_denoiser_deriv = function(x, scale = 1) {
  out = sqrt(scale) / cosh(sqrt(scale) * x)^2
  return(out)
}

# Soft-thresholding
st_denoiser = function(x, tau = 0.5) {
  out = sapply(x, function(x, tau) {
    if (abs(x) < tau) {
      return(0)
    } else {
      return(sign(x) * (abs(x) - tau))
    }
  },
  tau = tau)
  return(out)
}

st_denoiser_deriv = function(x, tau = 0.5) {
  out = sapply(x, function(x, tau) {
    if (abs(x) < tau) {
      return(0)
    } else {
      return(1)
    }
  },
  tau = tau)
  return(out)
}
