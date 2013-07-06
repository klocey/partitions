#!/usr/bin/python

import sys
import numpy as np
from scipy import stats
import random, decimal
from random import choice
import re
import math
import itertools


""" 
    Functions for integer partitioning. Most apply to using integer partitioning
    to examine distributions of wealth and abundance using the feasible set. The
    feasible set is the set all forms of the distribution having the same
    constraint values (e.g. total abundance N, species richness S).
    
    Included functions:
    - conjugate(): get the conjugate of an integer partition (recoded from Sage,
      see below)
    - NrParts(): Find the number of partitions for a given total N and number of
      parts S (modified and recoded from GAP, see below)
    
    Five functions to generate uniform random integer partitions of Q having N
    parts.  Each allows the option to have partitions with zero values.
    - bottom_up(): Starts at smallest possible value of the largest possible
      part (K) 
    - divide_and_conquer(): Start at random points in the feasible set
    - top_down(): Starts at largest possible value of the largest possible part
      (K) 
    - multiplicity(): uses the top_down approach but builds the partition using
      multiplicities (i.e. multiples of parts)
    - best(): Calls one of the above functions depending on the value of Q and
      Q / N (restricted to bottom_up and divide_and_conquer for now)

    Four partitioning functions not provided in other softwares
    - most_even_partition(): Get the last lexical (i.e. most even) partition of
      Q having N parts (no zeros)
    - min_max(): Get the smallest possible maximum part a partition of Q having
      N parts (no zeros)
    - firstpart(): Get the first lexical partition of Q having N parts with k as
      the largest part (no zeros)
    - next_restricted_part(): Get the next lexical partition of Q having N parts
    
"""


def P(D, q, k):
    """
    number of partitions of q with k or less parts (or having k or less as the
    largest part), i.e. P(q + k, k).
    Arguments:
        D : a dictionary for the number of partitions of Q having N or less
            parts (or N or less as the largest part), i.e. P(Q, Q + N).   
        q : total sum of the set
        k : number of parts      
        
    """
    if (q, k) not in D:
        D[(q, k)] = NrParts(q + k, k)
    return [D, D[(q, k)]] # return the updated dictionary and P(q + k, k).


def conjugate(part):
    """
    Find the conjugate of an integer partition. Recoded (on 24-Apr-2013) from
    the Sage source code: www.sagenb.org/src/combinat/partition.py
    
    """
    
    if part == []:
        return []
    else:
        l = len(part)
        conj =  [l] * part[-1]
        for i in xrange(l - 1, 0, -1):
            conj.extend([i] * (part[i - 1] - part[i]))
        return conj


def NrParts(*arg):
    """
    Find the number of partition for a given total Q and number of parts N.
    Recoded on 24-Apr-2013 and modified from GAP source code: www.gap-system.org
    Note: p(Q) = p(Q + Q, Q) thus NrParts(Q) returns the same value as NrParts
    (Q + Q, Q)
    
    Arguments:
    *arg : either Q or Q and N
    
    """

    if len(arg) == 1:  # using p(Q) = p(Q + Q, Q)
        Q = arg[0] * 2
        N = arg[0]
    elif len(arg) == 2:    
        Q = arg[0]
        N = arg[1]
    parts = 0
    if Q == N or N == 1:
        parts = 1
    elif Q < N or N == 0:
        parts = 0
    else:
        p = [1] * Q
        for i in range(2, N + 1):  
            for m  in range(i + 1, Q - i + 1 + 1):
                p[m] = p[m] + p[m - i]
        parts = p[Q - N + 1]
    return parts

def rand_parts(Q, N, sample_size, method, D={}, zeros=False):
    """
    Generate uniform random partitions of Q having N parts.
    
    Arguments:
        Q : Total sum across parts
        N : Number of parts to sum over
        sample_size : number of random partitions to generate
        method : method to use for generating the partition, options include:
            'bottom_up', 'top_down', 'divide_and_conquer', 'multiplicity', and
            'best'
        D : a dictionary for the number of partitions of Q having N or less
            parts (or N or less as the largest part), i.e. P(Q, Q + N). Defaults
            to a blank dictionary.
        zeros : boolean if True partitions can have zero values, if False
            partitions have only positive values, defaults to False
    
    Notes:
        method == 'best' attempts to use the values of Q and N to infer what the 
        fastest method to comput the partition.
    
    """
    parts = []
    if zeros:
        _list = P(D, Q, N)    
    else:
        _list = P(D, Q - N, N)
    D = _list[0]
    numparts = _list[1]        
    while len(parts) < sample_size:
        which = random.randrange(1, numparts + 1)
        if zeros:
            q = int(Q)
            part = []
        else:
            q = int(Q - N)
            part = [N]
        if method == 'bottom_up':
            part = bottom_up(part, q, D, which)
        if method == 'top_down':
            part = top_down(part, q, D, which)
        if method == 'divide_and_conquer':
            part = divide_and_conquer(part, q, N, D, which)
        if method == 'multiplicity':
            part = multiplicity(part, q, D, which)
        if method == 'best':
            if Q < 250 or N >= Q / 1.5:
                parts = bottom_up(part, q, D, which)
            else:
                parts = divide_and_conquer(part, q, N, D, which)
        if zeros:
            Zs = [0] * (N - len(part))
            part.extend(Zs)
        parts.append(part)
    return parts


def bottom_up(part, q, D, which):
    """
    Bottom up method of generating uniform random partitions of Q having N parts.
    
    Arguments:
        part : a list to hold the partition
        q : The total sum of the partition
        D : a dictionary for the number of partitions of Q having N or less
            parts (or N or less as the largest part), i.e. P(Q, Q + N).        
        which : 

    """    
    while q > 0:
        for k in range(1, q + 1):
            _list = P(D, q, k) # number of partitions of q having k or less as the largest part
            D = _list[0]
            count = _list[1]
            if count >= which:
                _list = P(D, q, k - 1)
                D = _list[0]
                count = _list[1]
                break
        part.append(k)
        q -= k
        if q == 0:
            break
        which -= count
    part = conjugate(part)    
    return(part)


def top_down(part, q, D, which):
    """
    Top down method of generating uniform random partitions of Q having N parts.
    
    Arguments:
        part : a list to hold the partition
        q : The total sum of the partition
        D : a dictionary for the number of partitions of Q having N or less
            parts (or N or less as the largest part), i.e. P(Q, Q + N).        
        which : 

    """    
    while q > 1:   ## Shouldn't this be q > 0 ? DJM
        if part: 
            x = min(part)
        else: 
            x = q
        for k in reversed(range(1, x + 1)):
            _list = P(D, q, k) # number of partitions of q having k or less as the largest part
            D = _list[0]
            count = _list[1]
            if count < which:
                k += 1
                break
        which -= count
        part.append(k)
        q -= k
        if q == 1:
            part.append(1)
            break
        if q <= 0:
            break
    part = conjugate(part)
    return(part)


def divide_and_conquer(part, q, N, D, which):
    """
    Divide and conquer method of generating uniform random partitions of Q
    having N parts.
        
    Arguments:
        part : a list to hold the partition
        q : The total sum of the partition
        N : Number of parts to sum over
        D : a dictionary for the number of partitions of Q having N or less
            parts (or N or less as the largest part), i.e. P(Q, Q + N).        
        which : 

    """
    _max = int(N)
    _min = int(1)
    while q > 0:
        k = random.randrange(_min, _max + 1)
        _list = P(D, q, k)
        D = _list[0]
        upper = _list[1]
        _list = P(D, q, k - 1)
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
            _min = k + 1    
        elif which <= lower:
            _max = k - 1    
    part = conjugate(part)
    return part


def get_multiplicity(q, k, D, which, count): 
    """ 
    Find the number of times a value k occurs in a partition that is being
    generated at random by the multiplicity() function. The resulting
    multiplicity is then passed back to the multiplicity() function along with
    an updated value of count and an updated dictionary D
    
    Arguments:
        q : 
        k : 
        D : a dictionary for the number of partitions of Q having N or less
            parts (or N or less as the largest part), i.e. P(Q, Q + N).                
        which :
        count : count < which
    
    """
    
    multi = [] # the multiplicity 
    f = 1
    while f:
        _list = P(D, (q - k * f), k - 1)
        D = _list[0]
        count += _list[1]
        if count >= which:
            count -= _list[1]
            multi = [k] * f
            break                
        f += 1
    return [multi, count, D]


def multiplicity(part, q, D, which):
    """
    multiplicity method of generating uniform random partitions of Q having N
    parts.
    
    Arguments:
        part : a list to hold the partition
        q : The total sum of the partition
        D : a dictionary for the number of partitions of Q having N or less
            parts (or N or less as the largest part), i.e. P(Q, Q + N).        
        which : 

    """
    while q > 0:
        multi = []
        if part:
            x = min(part)
        else: 
            x = int(q)
        for k in reversed(range(1, x + 1)): # start with largest k
            _list = P(D, q, k) # number of partitions of q having k or less as the largest part
            D = _list[0]
            count = _list[1]
            if count == which and which == 1:
                multi = [1] * q
                q=0
                break
            if count < which: # k has been found
                k += 1
                _list = get_multiplicity(q, k, D, which, count) # now, find how many times k occurs, i.e. the multiplicity of k 
                multi = _list[0]
                count = _list[1]
                D = _list[2]
                break
        q -= sum(multi)
        part.extend(multi)
        which -= count    
    part = conjugate(part)
    return part

        
def most_even_partition(Q, N):
    """ Find the last lexical (i.e. most even) partition of Q having N parts """
    
    most_even = [int(math.floor(float(Q) / float(N)))] * N
    _remainder = int(Q % N)
    j = 0
    while _remainder > 0:
        most_even[j] += 1
        _remainder -= 1
        j += 1
    return most_even


def min_max(Q, N):
    """ Find the smallest possible maximum part a partition of Q having N parts """
    
    _min = int(math.floor(float(Q) / float(N)))
    if int(Q % N) > 0:
        _min +=1
    return _min

    
def firstpart(Q, N, k):
    """ Find the first lexical partition of Q having N parts with k as the largest part """
    
    part = []
    if k == None:
        part.append(Q - N + 1)
        ones = [1] * (N - 1)
        part.extend(ones)
        return part
    
    elif k < min_max(Q, N):
        return None
        
    else:
        part.append(k)
        Q -= k
        N -= 1
        while Q > 0:
            k = min(k, Q - N + 1)
            part.append(k)
            Q -= k
            N -= 1
        
    return part


def next_restricted_part(p):
    """ Find the next lexical partition of Q having N parts """
    
    Q = sum(p)
    N = len(p)
    if p == most_even_partition(Q, N):
        return firstpart(Q, N, None)

    for i in enumerate(reversed(p)):
        if i[1] - p[-1] > 1:
            if i[0] == (N - 1):
                p = firstpart(Q, N, int(i[1] - 1))
                return p
            else:
                parts = np.split(p, [N - i[0] - 1])
                h1 = list(parts[0])
                h2 = list(parts[1])
                next = list(firstpart(int(sum(h2)), int(len(h2)), int(h2[0]) - 1))
                return h1 + next
