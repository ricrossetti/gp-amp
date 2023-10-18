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
