#!/usr/bin/env sage -python
import sys
import numpy as np
from scipy import stats
import random
from random import choice
import re
import math
import random, decimal

""" Functions for integer partitioning. Most apply to partitions of N having S parts """

# Find the conjugate of an integer partition
# Taken from the Sage source code: http://www.sagenb.org/src/combinat/partition.py

def conjugate(part):
        
    if part == []:
        return []
    else:
        l = len(part)
        conj =  [l]*part[-1]
        for i in xrange(l-1,0,-1):
            conj.extend([i]*(part[i-1] - part[i]))
        return conj
            

# Find the number of partition for a given total N and number of parts S
# Taken from GAP source code: http://www.gap-system.org/

def NrParts(N,S):
    
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


# Generate a uniform random partition of n having s parts. 
# Derived by Ken J. Locey
    
def random_partition(n,s):
    
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

# Generate a sample of partitions of N having S parts 
def rand_parts(N,S,sample_size):
    parts = []
    while len(parts) < sample_size:
        part = random_partition(N,S)
        parts.append(part)
    
    return parts


# Find the last lexical (i.e. most even) partition of N having S parts
# derived by Ken J. Locey
def most_even_partition(n,s):
    most_even = [int(floor(float(n)/float(s)))]*s
    _remainder = int(n%s)
    
    j = 0
    while _remainder > 0:
        most_even[j] += 1
        _remainder -= 1
        j += 1
    return most_even


# Find the next lexical partition of N having S parts
# Derived by Ken J. Locey, uses Sage but this will change

def portion(alist, indices):

    return [alist[i:j] for i, j in zip([0]+indices, indices+[None])]

def next_restricted_part(p,n,s):
    
    if p == most_even_partition(n,s):return Partitions(n,length=s).first()
    
    for i in enumerate(reversed(p)):
        if i[1] - p[-1] > 1:
            if i[0] == (s-1):
                return Partitions(n,length=s,max_part=(i[1]-1)).first()
            else:
                parts = portion(p,[s-i[0]-1]) # split p into the part that won't change and the part that will
                h1 = parts[0]
                h2 = parts[1]
                next = list(Partitions(sum(h2),length=len(h2),max_part=(h2[0]-1)).first())
                return h1+next
                

