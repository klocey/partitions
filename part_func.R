## Functions for integer partitioning. Most apply to using integer partitioning to examine distributions
# of wealth and abundance using the feasible set. The feasible set is the set all forms of the distribution
# having the same constraint values (e.g. total abundance N, species richness S). To my knowledge, no
# mathematical environments have functions 3 through 10.
# 
# Included functions:
# 
# 1. conjugate(): get the conjugate of an integer partition (recoded from Sage, see below)
# 2. NrParts(): Find the number of partitions for a given total N and number of parts S (modified and recoded from GAP, see below)
# 
# 3. rand_parts1(): Generate uniform random integer partitions of n having s parts. Starts at small end of the feasible set. 
# 4. rand_parts2(): Generate uniform random integer partitions of n having s parts. Starts at random points in the feasible set.
# 
# 5. rand_parts_zero1(): Generate uniform random partitions of n having s parts, where some parts may = 0. Starts at small end of the feasible set.
# 6. rand_parts_zero2(): Generate uniform random partitions of n having s parts, where some parts may = 0. Starts at random points in the feasible set.
# 
# 7. most_even_partition(): Get the last lexical (i.e. most even) partition of N having S parts (no zeros)
# 8. min_max(): Get the smallest possible maximum part a partition of N having S parts (no zeros)
# 9. firstpart(): Get the first lexical partition of N having S parts with k as the largest part (no zeros)
# 10. next_restricted_part(): Get the next lexical partition of N having S parts """

dyn.load("partitions.dll")

library(hash)

lappend <- function (lst, ...){
  lst <- c(lst, list(...))
  return(lst)
}

last = function(x) { tail(x, n = 1) }

conjugate = function(part){ # Find the conjugate of an integer partition
  # Recoded (orginally on 24-Apr-2013) from the Sage source code:
  # http://www.sagenb.org/src/combinat/partition.py
  if (is.null(part))
    conj = NULL
  else {
    l = length(part)
    conj = rep(l, last(part))
    for (i in (l - 1):1)
      conj = c(conj, rep(i, part[i] - part[i+1]))
  }
  return(conj)
}

NrParts = function(N, S, use_c = FALSE){ 
  # Find the number of partition for a given total N and number of parts S
  # Recoded (originally on 24-Apr-2013) and modified from GAP source code:
  # http://www.gap-system.org/
  s = 0
  if (N == S | S == 1) {
    s = 1
  }
  else if (N < S | S == 0) {
      s = 0
  }
  else {
    n = N
    k = S
    p = rep(1, n)
    if (use_c) {
      p = .C("NrParts", n = as.integer(n), k = as.integer(k), p = as.double(p))$p
    }
    else {
      for (i in 2:k) {
        if ((i + 1) <= (n - i + 1)) {
          for (m in (i + 1):(n - i + 1)) {
            p[m + 1] = p[m + 1] + p[m - i + 1]
          }
        }
      }  
    }  
    s = p[n - k + 2]
  }  
  return(s)
}

rand_parts1 = function(N, S, sample_size, use_c=FALSE) { 
  # Generate a uniform random partition of n having s parts.
  D = hash()
  P = function(n, x, use_c) {
    key = paste(n, x, sep=',')
    if (!has.key(key, D))
      D[key] = NrParts(n + x, x, use_c)
    # number of partitions of n having x or less as the largest part
    return(D[[key]])
  }
  parts = list()
  numparts = P(N - S, S, use_c)
  while (length(parts) < sample_size) {
    n = N - S
    part = S # first element of part must equal s (because the conjugate must have s parts)
    which = sample(numparts, 1)
    while (n > 0) {
      for (k in 1:n) {
        count = P(n, k, use_c) # number of partitions of N having K or less as the largest part
        if (count >= which) {
          count = P(n, k - 1, use_c)
          break
        }
      }  
      part = c(part, k)
      n = n - k
      if (n == 0)
        break
      which = which - count
    }
    part = conjugate(part)
    parts= lappend(parts, part)
  }  
  return(parts)
}

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

def rand_parts_zero1(N,S,sample_size):
  """ Generate a uniform random partition of n having s parts, where some parts
may have zero values """

D = {}
def P(n,x):
  if (n,x) not in D:
  D[(n,x)] = NrParts(n+x,x)
return D[(n,x)] # number of partitions of n with s parts having x or less as the largest part

parts = []
numparts = P(N,S)

while len(parts) < sample_size:
  
  n = int(N)
part = []
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
_len = len(part)
if _len < S:
  zeros = [0]*(S-_len)
part.extend(zeros)
parts.append(part)

return parts

def rand_parts_zero2(N,S,sample_size): # Generate a uniform random partition of n having k parts.
  
  D = {}
def P(n,x):
  if (n,x) not in D:
  D[(n,x)] = NrParts(n+x,x)
return D[(n,x)] # number of partitions of n with s parts having x or less as the largest part

parts = []
numparts = P(N,S)

while len(parts) < sample_size:
  
  n = int(N)
part = []
which = random.randrange(1,numparts+1)
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
_len = len(part)
if _len < S:
  zeros = [0]*(S-_len)
part.extend(zeros)
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
  
  part = []
if k == None:
  part.append(N-S+1)
ones = [1]*(S-1)
part.extend(ones)
return part

if k < min_max(N,S):
  return None

else:
  part.append(k)
N -= k
S -= 1
while N > 0:
  k = min(k,N-S+1)
part.append(k)
N -= k
S -= 1

return part


# The 2 functions below find the next lexical partition of N having S parts

def next_restricted_part(p):
  n = sum(p)
s = len(p)
if p == most_even_partition(n,s):
  return firstpart(n,s,None)

for i in enumerate(reversed(p)):
  if i[1] - p[-1] > 1:
  if i[0] == (s-1):
  p = firstpart(n,s,int(i[1]-1))
return p
else:
  parts = np.split(p,[s-i[0]-1])
h1 = list(parts[0])
h2 = list(parts[1])
next = list(firstpart(int(sum(h2)),int(len(h2)),int(h2[0])-1))
return h1+next

