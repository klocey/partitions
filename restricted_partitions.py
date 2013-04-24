from sage.all import *
import random

""" Three functions that, together, generate uniform random partitions of N having S parts.
This method slows considerably for some values of S once N is greater than a few hundred.
Overall, this method is less efficient than that located in misc_part_func.py of the partitions repository."""


""" First, a function to find the smallest maximum addend for a partition of n with s parts """

def min_max(n,s):

    _min = int(floor(float(n)/float(s)))
    if int(n%s) > 0:
        _min +=1

    return _min


""" Next, A function that uses a cache and memoiziation to find the number of partitions of n with
    s parts having x as the largest part. This is fast, but I think there's a more elegant solution
    to be had. e.g., Often: P(N,S,max=K) = P(N-K,S-1) Thanks to ante (http://stackoverflow.com/users/494076/ante)
    for helping me with this: Finding the number of integer partitions given a total, a number of parts,
    and a maximum summand """

D = {}
def P(n,s,x):
    if n > s*x or x <= 0: return 0
    if n == s*x: return 1
    if (n,s,x) not in D:
        D[(n,s,x)] = sum(P(n-i*x, s-i, x-1) for i in xrange(s))
    return D[(n,s,x)]


""" Finally, a function to find uniform random partitions of n with s parts, with no rejection rate! Each randomly chosen number codes for a specific partition of n having s parts. """

def random_partition(n,s):
    S = s
    partition = []
    _min = min_max(n,S)
    _max = n-S+1

    total = number_of_partitions(n,S)
    which = random.randrange(1,total+1) # random number

    while n:
        for k in range(_min,_max+1):
            count = P(n,S,k)
            if count >= which:
                count = P(n,S,k-1)
                break

        partition.append(k)
        n -= k
        if n == 0: break
        S -= 1
        which -= count
        _min = min_max(n,S)
        _max = k

    return partition
