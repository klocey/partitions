#!/usr/bin/python

import sys
import numpy as np
from scipy import stats
import random, decimal
from random import choice
import re
import math
import itertools

""" Functions for integer partitioning. Most apply to using integer partitioning to examine distributions
    of wealth and abundance using the feasible set. The feasible set is the set all forms of the distribution
    having the same constraint values (e.g. total abundance N, species richness S). To my knowledge, no
    mathematical environments provide the last 4 functions listed below.

Included functions:

- conjugate(): get the conjugate of an integer partition (recoded from Sage, see below)
- NrParts(): Find the number of partitions for a given total N and number of parts S (modified and recoded from GAP, see below)

Three functions to generate uniform random integer partitions of Q having N parts.
Each allows the option to have summands with zero values
- bottum_up(): Starts at small end of the feasible set 
- divide_and_conquer(): Start at random points in the feasible set
- best(): Calls one of the above two functions depending on the value of Q and Q/N

Four partitioning functions not provided in other softwares
- most_even_partition(): Get the last lexical (i.e. most even) partition of Q having N parts (no zeros)
- min_max(): Get the smallest possible maximum part a partition of Q having N parts (no zeros)
- firstpart(): Get the first lexical partition of Q having N parts with k as the largest part (no zeros)
- next_restricted_part(): Get the next lexical partition of Q having N parts """


def P(D,q,k):
    """ number of partitions of q with k or less parts (or having k or less as the largest part),
        i.e. P(q+k,k). D is a dictionary of P(q+k,k) values. """
        
    if (q,k) not in D:
        D[(q,k)] = NrParts(q+k,k)
    return [D, D[(q,k)]] # return the updated dictionary and P(q+k,k).

def conjugate(part):
    """ Find the conjugate of an integer partition. Recoded (on 24-Apr-2013)
        from the Sage source code: www.sagenb.org/src/combinat/partition.py """
    
    if part == []:
        return []
    else:
        l = len(part)
        conj =  [l]*part[-1]
        for i in xrange(l-1,0,-1):
            conj.extend([i]*(part[i-1] - part[i]))
        return conj
            
def NrParts(*arg):
    """ Find the number of partition for a given total Q and number of parts N. Recoded
        (on 24-Apr-2013) and modified from GAP source code: www.gap-system.org
        
        Note: p(Q) = p(Q+Q,Q) ...so NrParts(Q) returns the same value as NrParts(Q+Q,Q)
        ...wish I would've realized that before translating the block below 'if len(arg)==1:'
        from GAP source code. Well, now there's two waybotts to get p(Q). """
    
    parts = 0
    if len(arg) == 1:
        Q = arg[0]
        parts = 1                             # p(0) = 1
        p = [1]*(Q+1)
        for i in range(1,Q+1):
            parts = 0
            k = 1
            l = 1                         # k*(3*k-1)/2
            while 0 <= i-(l+k):
                parts = parts - (-1)**k * (p[i-l] + p[i-(l+k)])
                k = k + 1
                l = l + 3*k - 2
            
            if 0 <= i-l:
                parts = parts - (-1)**k * p[i-l]
            p[i] = parts
    
    elif len(arg) == 2:    
        Q = arg[0]
        N = arg[1]
        parts=0
        if Q == N or N == 1:
            parts = 1
        elif Q < N or N == 0:
            parts = 0
        else:
            q = int(Q)
            k = int(N)
            p = [1]*q
        
            for i in range(2,k+1):  
                for m  in range(i+1,q-i+1+1):
                    p[m] = p[m] + p[m-i]
            
            parts = p[q-k+1]
    
    return parts;


def bottom_up(D,Q,N,sample_size,zeros):
    """ Generate uniform random partitions of Q having N parts. D is a library that is passed
    back and forth between functions and frequently updated; it holds values for the number of
    partitions of Q having N or less parts (or N or less as the largest part), i.e. P(Q,Q+N).
    zeros is either 'yes' (i.e. summands can have zero values) or 'no' (i.e. summands have only
    positive values).  """
        
    parts = []
    if zeros == 'no':
        _list = P(D,Q-N,N)
        D = _list[0]
        numparts = _list[1]
        
    elif zeros == 'yes':
        _list = P(D,Q,N)
        D = _list[0]
        numparts = _list[1]
    
    while len(parts) < sample_size:
        which = random.randrange(1,numparts+1)
        if zeros == 'no':
            q = int(Q-N)
            part = [N]
        elif zeros == 'yes': 
            q = int(Q)
            part = []
            
        while q:
            for k in range(1,q+1):
                _list = P(D,q,k) # number of partitions of q having k or less as the largest part
                D = _list[0]
                count = _list[1]
                
                if count >= which:
                    _list = P(D,q,k-1)
                    D = _list[0]
                    count = _list[1]
                    break
            
            part.append(k)
            q -= k
            if q == 0: break
            which -= count
        
        part = conjugate(part)
        if zeros == 'yes':
            Zs = [0]*(N - len(part))
            part.extend(Zs)
        parts.append(part)
    
    return parts

      
def divide_and_conquer(D,Q,N,sample_size,zeros):
    """ Generate uniform random partitions of Q having N parts. D is a library that is passed,
        which holds values for the number of partitions of Q having N or less parts (or N or
        less as the largest part), i.e. P(Q,Q+N). zeros is either 'yes' (i.e. summands can have
        zero values) or 'no' (i.e. summands only have positive values). """
   
    parts = []
    if zeros == 'no':
        _list = P(D,Q-N,N)
        D = _list[0]
        numparts = _list[1]
        
    elif zeros == 'yes':
        _list = P(D,Q,N)
        D = _list[0]
        numparts = _list[1]
    
    while len(parts) < sample_size:    
        _max = int(N)
        _min = int(1)
        which = random.randrange(1,numparts+1)
        
        if zeros == 'no':
            q = int(Q-N)
            part = [N]
        elif zeros == 'yes': 
            q = int(Q)
            part = []
        
        while q > 0:
            k = random.randrange(_min, _max + 1)
            _list = P(D,q,k)
            D = _list[0]
            upper = _list[1]
            
            _list = P(D,q,k-1)
            D = _list[0]
            lower = _list[1]
            if lower < which and which <= upper: 
                part.append(k)
                q -= k
                _max = k
                _min = 1
                num = int(upper - lower)
                which = random.randrange(1, num + 1)
                
            elif which > upper:
                _min = k+1    
            elif which <= lower:
                _max = k-1            
            
        part = conjugate(part)
        if zeros == 'yes':
            Zs = [0]*(N - len(part))
            part.extend(Zs)
        parts.append(part)
    
    return parts


def rand_parts(Q,N,sample_size,which,zeros):
    """ Generate uniform random partitions of Q having N parts. """
    
    D = {} 
    if which == 'divide_and_conquer': parts = divide_and_conquer(D,Q,N,sample_size,zeros)
    if which == 'bottom_up':          parts = bottom_up(D,Q,N,sample_size,zeros)
    if which == 'best':
        if Q < 250 or N >= Q/1.5:
            parts = bottom_up(D,Q,N,sample_size,zeros)
        else:
            parts = divide_and_conquer(D,Q,N,sample_size,zeros)  
                   
    return parts


        
def most_even_partition(Q,N):
    """ Find the last lexical (i.e. most even) partition of Q having N parts """
    
    most_even = [int(math.floor(float(Q)/float(N)))]*N
    _remainder = int(Q%N)
    
    j = 0
    while _remainder > 0:
        most_even[j] += 1
        _remainder -= 1
        j += 1
    
    return most_even


def min_max(Q,N):
    """ Find the smallest possible maximum part a partition of Q having N parts """
    
    _min = int(math.floor(float(Q)/float(N)))
    if int(Q%N) > 0:
        _min +=1
    
    return _min

    
def firstpart(Q,N,k):
    """ Find the first lexical partition of Q having N parts with k as the largest part """
    
    part = []
    if k == None:
        part.append(Q-N+1)
        ones = [1]*(N-1)
        part.extend(ones)
        return part
    
    elif k < min_max(Q,N):
        return None
        
    else:
        part.append(k)
        Q -= k
        N -= 1
        while Q > 0:
            k = min(k,Q-N+1)
            part.append(k)
            Q -= k
            N -= 1
        
    return part


def next_restricted_part(p):
    """ Find the next lexical partition of Q having N parts """
    Q = sum(p)
    N = len(p)
    if p == most_even_partition(Q,N):
        return firstpart(Q,N,None)

    for i in enumerate(reversed(p)):
        if i[1] - p[-1] > 1:
            if i[0] == (N-1):
                p = firstpart(Q,N,int(i[1]-1))
                return p
            else:
                parts = np.split(p,[N-i[0]-1])
                h1 = list(parts[0])
                h2 = list(parts[1])
                next = list(firstpart(int(sum(h2)),int(len(h2)),int(h2[0])-1))
                return h1+next
