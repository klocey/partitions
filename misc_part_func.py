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

  
def rand_parts1(N,S,sample_size): # Generate a uniform random partition of n having s parts.
    
    D = {}
    def P(n,x):
        if (n,x) not in D:
            D[(n,x)] = NrParts(n+x,x)
        return D[(n,x)] # number of partitions of n having x or less as the largest part
    
    parts = []
    numparts = P(N-S,S)
    
    while len(parts) < sample_size:
        
        n = N-S
        part = [S] # first element of part must equal s (because the conjugate must have s parts)
        
        which = random.randrange(1,numparts+1)
    
        while n:
            for k in range(1,n+1):
                count = P(n,k) # number of partitions of N having K or less as the largest part
                if count >= which:
                    count = P(n,k-1)
                    break
            
            part.append(k)
            n -= k
            if n == 0: break
            which -= count
        
        part = conjugate(part)
        parts.append(part)
    return parts

      
def rand_parts2(N,S,sample_size): # Generate a uniform random partition of n having k parts.
    
    D = {}
    def P(n,x):
        if (n,x) not in D:
            D[(n,x)] = NrParts(n+x,x)
        return D[(n,x)] # number of partitions of n with s parts having x or less as the largest part
    
    parts = []
    numparts = P(N-S,S)
    
    while len(parts) < sample_size:
    
        which = random.randrange(1,numparts+1)
        n = int(N-S)
        part = [S]
        _max = int(S)
        _min = int(1)
        
        while n > 0:
            k = random.randrange(_min, _max + 1)
            upper = int(P(n,k))
            lower = int(P(n,k-1))
            if lower < which and which <= upper: 
                part.append(k)
                n -= k
                _max = k
                _min = 1
                num = int(upper - lower)
                which = random.randrange(1, num + 1)
                
            elif which > upper:
                _min = k+1    
            elif which <= lower:
                _max = k-1            
            
        part = conjugate(part)
        parts.append(part)
    
    return parts
    
    
def most_even_partition(n,s): # Find the last lexical (i.e. most even) partition of N having S parts
    
    most_even = [int(math.floor(float(n)/float(s)))]*s
    _remainder = int(n%s)
    
    j = 0
    while _remainder > 0:
        most_even[j] += 1
        _remainder -= 1
        j += 1
    return most_even



def min_max(n,s): # Find the smallest possible maximum part a partition of N having S parts

    _min = int(math.floor(float(n)/float(s)))
    if int(n%s) > 0:
        _min +=1

    return _min
    
    
def firstpart(N,S,k): # Find the first lexical partition of N having S parts with k as the largest part
    
    if k == None:
        part = [N-S+1]
        ones = [1]*(S-1)
        part.extend(ones)
    return part
    
    if k < min_max(n,s):
        return None
        
    else:
        part = [k]
        N -= k
        S -= 1
        while N > 0:
            k = min(k,N-S+1)
            part.append(k)
            N -= k
            S -= 1
        
    return part
       
    
# The 2 functions below find the next lexical partition of N having S parts

def portion(alist, indices):

    return [alist[i:j] for i, j in zip([0]+indices, indices+[None])]

def next_restricted_part(p,n,s):
    
    #if p == most_even_partition(n,s):return firstpart(n,s,k=None).first()
    
    for i in enumerate(reversed(p)):
        if i[1] - p[-1] > 1:
            if i[0] == (s-1):
                return firstpart(n,s,k=(i[1]-1))
            else:
                parts = portion(p,[s-i[0]-1]) # split p into the part that won't change and the part that will
                h1 = parts[0]
                h2 = parts[1]
                next = firstpart(sum(h2),len(h2),k=h2[0]-1)
                return h1+next
                
