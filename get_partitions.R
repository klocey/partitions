
source('./part_func.R')

load_c()

cl_args = commandArgs(trailingOnly=TRUE)

len_cl_args = length(cl_args)

if (len_cl_args  > 0) {
  method = cl_args[1]
  Q = cl_args[2]
  N = cl_args[3]
  sample_size = cl_args[4]
  zeros = cl_args[5]
  use_dict = cl_args[6]
  use_c = cl_args[7]
}
if (len_cl_args == 0) {
  method = 'bottom_up'
  Q = 100
  N = 10
  sample_size = 10
  zeros = FALSE
  use_dict = FALSE
  use_c = TRUE
}


fix_arg = function(arg) {
  if (is.na(arg))
    arg = FALSE
  return(arg)
}

zeros = fix_arg(zeros)
use_dict = fix_arg(use_dict)
use_c = fix_arg(use_c)

eval(parse(text = paste(method,'(',paste('NULL', Q, N, sample_size, zeros,
                        use_dict, use_c, sep=','), ')', sep='')))



