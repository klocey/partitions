#!/usr/bin/env sage -python
import sys
import numpy as np
from scipy import stats
import random
from random import choice
import re
import math
import random, decimal

""" Functions for integer partitioning. The main function here is random_partition();
    other function support it. random_partition() is a variant of other random partitioning
    algorithms found in this repo. More functions to come. """

def conjugate(part): # Find the conjugate of an integer partition
    # Slightly recoded from the Sage source code: http://www.sagenb.org/src/combinat/partition.py
        
    if part == []:
        return []
    else:
        l = len(part)
        conj =  [l]*part[-1]
        for i in xrange(l-1,0,-1):
            conj.extend([i]*(part[i-1] - part[i]))
        return conj
            

def NrParts(N,S): # Find the number of partition for a given total N and number of parts S
    # Recoded and modified from GAP source code: http://www.gap-system.org/
    
    s=0
    if N == S or S == 1:
        s = 1
    elif N < S or S == 0:
        s = 0
    else:
        n = int(N)
        k = int(S)
        p = [1]*n
        
        for i in range(2,k+1):  
            for m  in range(i+1,n-i+1+1):
                p[m] = p[m] + p[m-i]
            
        s = p[n-k+1]

    return s;


   
def random_partition(n,s): # Generate a uniform random partition of n having s parts.
    
    numparts = NrParts(n,s)
    which = random.randrange(1,numparts+1)
    
    part = [s] # first element of part must equal s (because the conjugate must have s parts)
    n-= s # having found the first element, n is decreased by s
    
    while n:
        for k in range(1,n+1):
            count = NrParts(n+k,k) # number of partitions of N having K or less as the largest part
            if count >= which:
                count = NrParts(n+k-1,k-1)
                break
            
        part.append(k)
        n -= k
        if n == 0: break
        which -= count
        
    part = conjugate(part)
    return part


def rand_parts(N,S,sample_size): # Generate a sample of partitions of N having S parts 
    parts = []
    while len(parts) < sample_size:
        part = random_partition(N,S)
        parts.append(part)
    
    return parts
