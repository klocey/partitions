#!/usr/bin/env sage -python

from sage.all import *
import sys
import partitions as parts


if len(sys.argv) > 1:
    method = sys.argv[1]
    Q = int(sys.argv[2])
    N = int(sys.argv[3])
    sample_size = int(sys.argv[4])
    zeros = cl_args[5]
else:
    print('Error: must specify arguments at the command line')
    

# insert whatever the Sage command is to compute a random partition