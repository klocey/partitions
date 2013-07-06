# Functions for integer partitioning. Most apply to using integer partitioning
# to examine distributions of wealth and abundance using the feasible set. The
# feasible set is the set all forms of the distribution having the same
# constraint values (e.g. total abundance N, species richness S).
# 
# Included functions:
# - conjugate(): get the conjugate of an integer partition (recoded from Sage,
#   see below)
# - NrParts(): Find the number of partitions for a given total N and number of
#   parts S (modified and recoded from GAP, see below)
# 
# Five functions to generate uniform random integer partitions of Q having N
# parts.  Each allows the option to have partitions with zero values.
# - bottom_up(): Starts at smallest possible value of the largest possible
#   part (K) 
# - divide_and_conquer(): Start at random points in the feasible set
# - top_down(): Starts at largest possible value of the largest possible part
#   (K) 
# - multiplicity(): uses the top_down approach but builds the partition using
#   multiplicities (i.e. multiples of parts)
# - best(): Calls one of the above functions depending on the value of Q and
#   Q / N (restricted to bottom_up and divide_and_conquer for now)
# 
# Four partitioning functions not provided in other softwares
# - most_even_partition(): Get the last lexical (i.e. most even) partition of
#   Q having N parts (no zeros)
# - min_max(): Get the smallest possible maximum part a partition of Q having
#   N parts (no zeros)
# - firstpart(): Get the first lexical partition of Q having N parts with k as
#   the largest part (no zeros)
# - next_restricted_part(): Get the next lexical partition of Q having N parts
##

dyn.load("partitions.dll")

library(hash)

rand_int = function(min=0, max=1) {
  int = ceiling(runif(1, min - 1, max))
  return(int)
}

last = function(x) { tail(x, n = 1) }

P = function(D, q, k, use_c) {
  # number of partitions of q with k or less parts (or having k or less as the
  # largest part), i.e. P(q+k,k).
  # Arguments:
  #   D : lookup table for numbers of partitions, P(q+k,k) values. 
  #   q : total sum of the set
  #   k : number of parts
  #   use_c : logical, if TRUE the number of partitions is computed in c
  key = paste(q, k, sep=',')
  if (!has.key(key, D))
    D[key] = NrParts(q + k, k, use_c)
  return(list(D, D[[key]]))
}


conjugate = function(part, use_c=TRUE){ 
  # Find the conjugate of an integer partition
  # Recoded (orginally on 24-Apr-2013) from the Sage source code:
  # http://www.sagenb.org/src/combinat/partition.py
  # Arguments:
  # part : a vector that represents an integer partition
  # use_c : default is TRUE, the conjugate is computed in c
  if (is.null(part))
    conj = NULL
  else {
    l = length(part)
    if (use_c) {
      conj = rep(0, max(part))
      conj[1:last(part)] = rep(l, last(part))      
      j = last(part) + 1
      conj = .C("conjugate", l = as.integer(l), j = as.integer(j),
                part = as.integer(part), conj = as.integer(conj))$conj
    }            
    else {
      conj = rep(l, last(part))
      for (i in (l - 1):1)
        conj = c(conj, rep(i, part[i] - part[i + 1]))
    }  
  }
  return(conj)
}


NrParts = function(Q, N=NULL, use_c = TRUE){ 
  # Find the number of partition for a given total Q and number of parts N.
  # Recoded (on 24-Apr-2013) and modified from GAP source code: www.gap-system.org
  # Arguments:
  #   Q : Total sum
  #   N : Number of items to sum across
  #   use_c : default is TRUE, the number of partitions is computed in c  
  if (is.null(N)) {  # using p(Q) = p(Q + Q, Q)
    N = Q
    Q = Q * 2
  }
  parts = 0  
  if (Q == N | N == 1) {
    parts = 1
  }
  else if (Q < N | N == 0) {
    parts = 0
  }
  else {
    p = rep(1, Q)
    if (use_c) {
      p = .C("NrParts", Q = as.integer(Q), N = as.integer(N), p = as.double(p))$p
    }
    else {
      for (i in 2:N) {
        for (m in (i + 1):(Q - i + 1)) {
          p[m + 1] = p[m + 1] + p[m - i + 1]
        }
      }  
    }
    parts = p[Q - N + 2]
  }  
  return(parts)
}

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



bottom_up = function(D, Q, N, sample_size, zeros) {
  # Bottom up algo for generating uniform random partitions of Q having N parts.
  # Arguments:
  #   D: a dictionary for the number of partitions of Q having N or less parts
  #     (or N or less as the largest part), i.e. P(Q,Q+N).
  #   Q : Total sum across parts
  #   N : Number of parts to sum over
  #   sample_size : number of random partitions to generate
  #   zeros: boolean if TRUE summands can have zero values, if FALSE summands
  #     have only positive values
  parts = NULL
  if (!zeros) {
        
  }
}

rand_parts1 = function(N, S, sample_size, use_dict=FALSE, use_c=TRUE) { 
  # Generate a uniform random partition of n having s parts.
  # Arguments:
  # N : Total number of individuals
  # S : Total number of species
  # sample_size : number of random partitions to return
  # use_dict : default is FALSE, if TRUE then a dictionary/hash table is used to store
  #   numbers of partitions, if running c-code for computed number of partions
  #   then this tends to slow things down
  # use_c : default is TRUE, the number of partitions and the conjugate are computed
  #   in c rather than in R this makes for a big speed gain
  if (use_dict) {
    D = hash()
    P = function(n, x, use_c) {
      # compute number of partitions of n having x or less as the
      # largest part
      key = paste(n, x, sep=',')
      if (!has.key(key, D))
        D[key] = NrParts(n + x, x, use_c)
      return(D[[key]])
    }
  }  
  else {
    P = function(n, x, use_c) {
      numparts = NrParts(n + x, x, use_c)
      return(numparts)
    }
  }  
  part_matrix = matrix(NA, ncol=sample_size, nrow=S)
  numparts = P(N - S, S, use_c)
  ipart = 1
  while (ipart <= sample_size) {
    n = N - S
    part = S # first element of part must equal s (because the conjugate must have s parts)
    which = rand_int(1, numparts)
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
    part = conjugate(part, use_c)
    part_matrix[ , ipart] = part
    ipart = ipart + 1
  }  
  return(part_matrix)
}

rand_parts2 = function(N, S, sample_size, use_dict=FALSE, use_c=TRUE) {
  # Generate a uniform random partition of n having k parts.
  # Arguments:
  # N : Total number of individuals
  # S : Total number of species
  # sample_size : number of random partitions to return
  # use_dict : default is FALSE, if TRUE then a dictionary/hash table is used to store
  #   numbers of partitions, if running c-code for computed number of partions
  #   then this tends to slow things down
  # use_c : default is TRUE, the number of partitions and the conjugate are computed
  #   in c rather than in R this makes for a big speed gain
  if (use_dict) {
    D = hash()
    P = function(n, x, use_c) {
      # compute number of partitions of n having x or less as the
      # largest part
      key = paste(n, x, sep=',')
      if (!has.key(key, D))
        D[key] = NrParts(n + x, x, use_c)
      return(D[[key]])
    }
  }  
  else {
    P = function(n, x, use_c) {
      numparts = NrParts(n + x, x, use_c)
      return(numparts)
    }
  }  
  part_matrix = matrix(NA, ncol=sample_size, nrow=S)
  numparts = P(N - S, S, use_c)
  ipart = 1
  while (ipart <= sample_size) {
    which = rand_int(1, numparts)
    n = N - S
    part = S
    max_k = S
    min_k = 1
    while (n > 0) {
      k = rand_int(min_k, max_k)
      upper = P(n, k, use_c)
      lower = P(n, k - 1, use_c)
      if (lower < which & which <= upper) {
        part = c(part, k)
        n = n - k
        max_k = k
        min_k = 1
        num = upper - lower
        which = rand_int(1, num)
      }
      else if (which > upper)
        min_k = k + 1
      else if (which <= lower)
        max_k = k - 1
    }  
    part = conjugate(part, use_c)
    part_matrix[ , ipart] = part
    ipart = ipart + 1
  }
  return(part_matrix)
}

rand_parts_zero1 = function(N, S, sample_size, use_dict=FALSE, use_c=TRUE) { 
  # Generate a uniform random partition of n having s parts.
  # Arguments:
  # N : Total number of individuals
  # S : Total number of species
  # sample_size : number of random partitions to return
  # use_dict : default is FALSE, if TRUE then a dictionary/hash table is used to store
  #   numbers of partitions, if running c-code for computed number of partions
  #   then this tends to slow things down
  # use_c : default is TRUE, the number of partitions and the conjugate are computed
  #   in c rather than in R this makes for a big speed gain
  if (use_dict) {
    D = hash()
    P = function(n, x, use_c) {
      # compute number of partitions of n having x or less as the
      # largest part
      key = paste(n, x, sep=',')
      if (!has.key(key, D))
        D[key] = NrParts(n + x, x, use_c)
      return(D[[key]])
    }
  }  
  else {
    P = function(n, x, use_c) {
      numparts = NrParts(n + x, x, use_c)
      return(numparts)
    }
  }  
  part_matrix = matrix(NA, ncol=sample_size, nrow=S)
  numparts = P(N, S, use_c)
  ipart = 1
  while (ipart <= sample_size) {
    n = N 
    part = NULL
    which = rand_int(1, numparts)
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
    part = conjugate(part, use_c)
    len_part = length(part)
    if (len_part < S) {
      zeros = rep(0, S - len_part)
      part = c(part, zeros)
    }  
    part_matrix[ , ipart] = part
    ipart = ipart + 1
  }  
  return(part_matrix)
}

rand_parts_zero2 = function(N, S, sample_size, use_dict=FALSE, use_c=TRUE) {
  # Generate a uniform random partition of n having k parts.
  # Arguments:
  # N : Total number of individuals
  # S : Total number of species
  # sample_size : number of random partitions to return
  # use_dict : default is FALSE, if TRUE then a dictionary/hash table is used to store
  #   numbers of partitions, if running c-code for computed number of partions
  #   then this tends to slow things down
  # use_c : default is TRUE, the number of partitions and the conjugate are computed
  #   in c rather than in R this makes for a big speed gain  
  if (use_dict) {
    D = hash()
    P = function(n, x, use_c) {
      # compute number of partitions of n having x or less as the
      # largest part
      key = paste(n, x, sep=',')
      if (!has.key(key, D))
        D[key] = NrParts(n + x, x, use_c)
      return(D[[key]])
    }
  }  
  else {
    P = function(n, x, use_c) {
      numparts = NrParts(n + x, x, use_c)
      return(numparts)
    }
  }  
  part_matrix = matrix(NA, ncol=sample_size, nrow=S)
  numparts = P(N, S, use_c)
  ipart = 1
  while (ipart <= sample_size) {
    which = rand_int(1, numparts)
    n = N
    part = NULL
    max_k = S
    min_k = 1
    while (n > 0) {
      k = rand_int(min_k, max_k)
      upper = P(n, k, use_c)
      lower = P(n, k - 1, use_c)
      if (lower < which & which <= upper) {
        part = c(part, k)
        n = n - k
        max_k = k
        min_k = 1
        num = upper - lower
        which = rand_int(1, num)
      }
      else if (which > upper)
        min_k = k + 1
      else if (which <= lower)
        max_k = k - 1
    }  
    part = conjugate(part, use_c)
    len_part = length(part)
    if (len_part < S) {
      zeros = rep(0, S - len_part)
      part = c(part, zeros)
    }  
    part_matrix[ , ipart] = part
    ipart = ipart + 1
  }
  return(part_matrix)
}

rand_parts = function(N, S, sample_size, method, include_zeros=FALSE,
                      use_dict=FALSE, use_c=TRUE) {
  # Generate a uniform random partition of n having k parts.
  # Arguments:
  # N : Total number of individuals
  # S : Total number of species
  # sample_size : number of random partitions to return
  # method : which algo to use 1 or 2
  # include_zeros : default is FALSE, if TRUE zeros are included in the partitions
  # use_dict : default is FALSE, if TRUE then a dictionary/hash table is used to store
  #   numbers of partitions, if running c-code for computed number of partions
  #   then this tends to slow things down
  # use_c : default is TRUE, the number of partitions and the conjugate are computed
  #   in c rather than in R this makes for a big speed gain
  # 
  if (use_dict) {
    D = hash()
    P = function(n, x, use_c) {
      # compute number of partitions of n having x or less as the
      # largest part
      key = paste(n, x, sep=',')
      if (!has.key(key, D))
        D[key] = NrParts(n + x, x, use_c)
      return(D[[key]])
    }
  }  
  else {
    P = function(n, x, use_c) {
      numparts = NrParts(n + x, x, use_c)
      return(numparts)
    }
  }  
  part_matrix = matrix(NA, ncol=sample_size, nrow=S)
  if (include_zeros)
    numparts = P(N, S, use_c)
  else
    numparts = P(N - S, S, use_c)
  ipart = 1
  while (ipart <= sample_size) {
    which = rand_int(1, numparts)
    if (include_zeros) {
      n = N
      part = NULL 
    }
    else {
      n = N - S 
      part = S # first element of part must equal s (because the conjugate must have s parts)
    }
    if (method == 1) {
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
    }
    if (method == 2) {
      max_k = S
      min_k = 1
      while (n > 0) {
        k = rand_int(min_k, max_k)
        upper = P(n, k, use_c)
        lower = P(n, k - 1, use_c)
        if (lower < which & which <= upper) {
          part = c(part, k)
          n = n - k
          max_k = k
          min_k = 1
          num = upper - lower
          which = rand_int(1, num)
        }
        else if (which > upper)
          min_k = k + 1
        else if (which <= lower)
          max_k = k - 1
      }  
      
    }  
    part = conjugate(part, use_c)
    if (include_zeros) {
      len_part = length(part)
      if (len_part < S) {
        zeros = rep(0, S - len_part)
        part = c(part, zeros)
      }
    }  
    part_matrix[ , ipart] = part
    ipart = ipart + 1  
  }
  return(part_matrix)
}