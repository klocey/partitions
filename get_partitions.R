
source('./part_func.R')

load_c()

cl_args = commandArgs(trailingOnly=TRUE)

len_cl_args = length(cl_args)

if (len_cl_args > 0) {
  method = cl_args[1]
  Q = as.numeric(cl_args[2])
  N = as.numeric(cl_args[3])
  sample_size = as.numeric(cl_args[4])
  zeros = as.logical(cl_args[5])
  use_c = as.logical(cl_args[6])
  use_dict = as.logical(cl_args[7])
}
if (len_cl_args == 0) {
  method = 'best'
  Q = 100
  N = 10
  sample_size = 10
  zeros = FALSE
  use_c = TRUE
  use_dict = FALSE
}

print(method)
print(Q)
print(N)
print(sample_size)
print(zeros)
print(use_c)
print(use_dict)


fix_arg = function(arg) {
  if (is.na(arg))
    arg = FALSE
  return(arg)
}

zeros = fix_arg(zeros)
use_dict = fix_arg(use_dict)
use_c = fix_arg(use_c)

eval(parse(text = 
     'rand_parts(Q, N, sample_size, method, hash(), zeros, use_c, use_dict)'))



