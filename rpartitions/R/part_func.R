#' Generate a random integer between two integers
#'
#' @param min minimum value 
#' @param max maximum value
#' @export
#' @examples
#' get_rand_int(min=0, max=10)
get_rand_int = function(min=0, max=1) {
  int = ceiling(runif(1, min - 1, max))
  return(int)
}


#' Returns the last element of a vector
#' 
#' @param x a vector
#' @export
#' @examples
#' last(1:10)
#' last(letters[1:10])
last = function(x) { tail(x, n = 1) }


#' Number of partitions of q with k or less parts (or having k or less as the
#' largest part), i.e. P(q+k,k).
#' 
#' @param D lookup table for numbers of partitions, P(q+k,k) values. 
#' @param q total sum of the set
#' @param use_c logical, if TRUE the number of partitions is computed in c
#' @export
#' @examples
#' P(list(), 100, 10, FALSE, FALSE)
P = function(D, q, k, use_c, use_dict) {
  if (use_dict) {
    key = paste(q, k, sep=',')
    if (!has.key(key, D)) {
      D[key] = NrParts(q + k, k, use_c)
    }
    out = list(D, D[[key]])
  }
  else {
    D = NULL
    out = list(D, NrParts(q + k, k, use_c))
  }
  return(out)
}


#' Find the conjugate of an integer partition
#' Recoded (orginally on 24-Apr-2013) from the Sage source code:
#' http://www.sagenb.org/src/combinat/partition.py
#' 
#' @param part a vector that represents an integer partition
#' @param use_c logical, defaults to TRUE, the conjugate is computed in c
#' @export
#' @useDynLib 'rpartitions'
#' @examples
#' conjugate(c(3,3,1,1), F)
conjugate = function(part, use_c=TRUE){ 
  if (is.null(part)) {
    conj = NULL
  }  
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


#' Find the number of partitions for a given total Q and number of parts N.
#' Recoded and modified from GAP source code: www.gap-system.org
#'
#' @param Q Total sum
#' @param N Number of items to sum across
#' @param use_c logical, defaults to TRUE, the number of partitions is computed in c  
#' @export 
#' @useDynLib 'rpartitions'
#' @examples
#' NrParts(100, 10)
NrParts = function(Q, N=NULL, use_c=TRUE){ 
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


#' Generate uniform random partitions of Q having N parts.
#'
#' @param Q Total sum across parts
#' @param N : Number of parts to sum over
#' @param sample_size : number of random partitions to generate
#' @param method : method to use for generating the partition, options include:
#'       'bottom_up', 'top_down', 'divide_and_conquer', 'multiplicity', and
#'       'best'. Defaults to 'best'
#' @param D : a dictionary for the number of partitions of Q having N or less
#'        parts (or N or less as the largest part), i.e. P(Q, Q + N). Defaults
#'        to a blank dictionary.
#' @param zeros : boolean if True partitions can have zero values, if False
#'        partitions have only positive values, defaults to False
#' @param use_c : boolean if TRUE then compiled c code is used, defaults to TRUE
#' @param use_dict : boolean, if TRUE then a hash table is used, defaults to FALSE
#' @return A matrix where each column is a random partition
#' @note method == 'best' attempts to use the values of Q and N to infer what the 
#'         fastest method to compute the partition.
#' @export
#' @examples
#' rand_parts(100, 10, 5)
rand_parts = function(Q, N, sample_size, method='best', D=hash(), zeros=FALSE,
                      use_c=TRUE, use_dict=FALSE) {
  parts= matrix(NA, ncol=sample_size, nrow=N)
  if (zeros) {
    Plist = P(D, Q, N, use_c, use_dict)
  }  
  else {
    Plist = P(D, Q - N, N, use_c, use_dict)
  }  
  D = Plist[[1]]
  numparts = Plist[[2]]
  ipart = 1
  while (ipart <= sample_size) {
    rand_int = get_rand_int(1, numparts)
    if (zeros) {
      q = Q
      part = NULL
    }
    else {
      q = Q - N
      part = N
    }  
    if (method == 'bottom_up') {
      part = bottom_up(part, q, D, rand_int, use_c, use_dict)
    }  
    if (method == 'top_down') {
      part = top_down(part, q, D, rand_int, use_c, use_dict)
    }  
    if (method == 'divide_and_conquer') {
      part = divide_and_conquer(part, q, N, D, rand_int, use_c, use_dict)
    }  
    if (method == 'multiplicity') {
      part = multiplicity(part, q, D, rand_int, use_c, use_dict)
    }  
    if (method == 'best') { 
      if (Q < 250 | N >= Q / 1.5)
        part = bottom_up(part, q, D, rand_int, use_c, use_dict)
      else
        part = divide_and_conquer(part, q, N, D, rand_int, use_c, use_dict)
    }  
    if (zeros) {
      Zs = rep(0, N - length(part))
      part = c(part, Zs)
    }  
    parts[ , ipart] = part
    ipart = ipart + 1
  } 
  return(parts)
}


#' Bottom up method of generating uniform random partitions of Q having N parts.  
#'
#' @param part a list to hold the partition
#' @param q the total sum of the partition
#' @param D a dictionary for the number of partitions of Q having N or less
#'   parts (or N or less as the largest part), i.e. P(Q, Q + N).        
#' @param rand_int the random integer
#' @param use_c boolean if TRUE then compiled c code is used
#' @param use_dict boolean, if TRUE then hash dictionary is used
#' @export
#' @examples
#' bottom_up(c(5, 4), 4, hash(), 1, T, T)
bottom_up = function(part, q, D, rand_int, use_c, use_dict) {
  while (q > 0) {
    for (k in 1:q) {
      Plist = P(D, q, k, use_c, use_dict)
      D = Plist[[1]]
      count = Plist[[2]]
      if (count >= rand_int) {
        Plist = P(D, q, k - 1, use_c, use_dict)
        D = Plist[[1]]
        count = Plist[[2]]
        break
      }
    }  
    part = c(part, k) 
    q = q - k
    if (q == 0) {
      break
    }  
    rand_int = rand_int - count
  }
  part = conjugate(part, use_c)
  return(part)
}


#' Top down method of generating uniform random partitions of Q having N parts.  
#' 
#' @param part a list to hold the partition
#' @param q the total sum of the partition
#' @param D a dictionary for the number of partitions of Q having N or less
#'   parts (or N or less as the largest part), i.e. P(Q, Q + N).        
#' @param rand_int the random integer
#' @param use_c boolean if TRUE then compiled c code is used
#' @param use_dict boolean, if TRUE then hash dictionary is used
#' @export
#' @examples
#' top_down(c(5, 4), 4, hash(), 1, T, T)
top_down = function(part, q, D, rand_int, use_c, use_dict) {
  while (q > 0) {
    if (!is.null(part)) {
      x = min(part)
    }  
    else {
      x = q
    }  
    for (k in x:1) {
      Plist = P(D, q, k, use_c, use_dict) # number of partitions of q having k or less as the largest part
      D = Plist[[1]]
      count = Plist[[2]]
      if (count < rand_int) {
        k = k + 1
        break
      }
    }
    rand_int = rand_int - count
    part = c(part, k)
    q = q - k
  }
  part = conjugate(part, use_c)
  return(part)
}


#' Divide and conquer method of generating uniform random partitions of Q
#' having N parts.
#'
#' @param part a list to hold the partition
#' @param q the total sum of the partition
#' @param N Number of parts to sum over
#' @param D a dictionary for the number of partitions of Q having N or less
#'        parts (or N or less as the largest part), i.e. P(Q, Q + N).        
#' @param rand_int the random integer
#' @param use_c boolean if TRUE then compiled c code is used
#' @param use_dict boolean, if TRUE then hash dictionary is used
#' @export
#' @examples
#' divide_and_conquer(c(5, 4), 4, hash(), 1, T, T)
divide_and_conquer = function(part, q, N, D, rand_int, use_c, use_dict) {
  max_int = N
  min_int = 1 
  while (q > 0) {
    k = get_rand_int(min_int, max_int)
    Plist = P(D, q, k, use_c, use_dict)
    D = Plist[[1]]
    upper = Plist[[2]]
    Plist = P(D, q, k - 1, use_c, use_dict)
    D = Plist[[1]]
    lower = Plist[[2]]
    if (lower < rand_int & rand_int <= upper) {
      part = c(part, k)
      q = q - k
      max_int = k
      min_int = 1
      num = upper - lower
      rand_int = get_rand_int(1, num)
    }
    else if (rand_int > upper)
      min_int = k + 1
    else if (rand_int <= lower)
      max_int = k - 1
  }
  part = conjugate(part, use_c)
  return(part)
}


#' Find the number of times a value k occurs in a partition that is being
#' generated at random by the multiplicity() function. The resulting
#' multiplicity is then passed back to the multiplicity() function along with
#' an updated value of count and an updated dictionary D
#'
#' @param q 
#' @param k 
#' @param D a dictionary for the number of partitions of Q having N or less
#'        parts (or N or less as the largest part), i.e. P(Q, Q + N).                
#' @param rand_int the random integer
#' @param count count < rand_int
#' @param use_c boolean if TRUE then compiled c code is used
#' @param use_dict boolean, if TRUE then hash dictionary is used
#' @export
#' @examples
#' get_multiplicity(10, 5, hash(), 3, 2, T, T)
get_multiplicity = function(q, k, D, rand_int, count, use_c, use_dict){
  multi = NULL # the multiplicity 
  f = 1
  while (f > 0) {
    Plist = P(D, (q - k * f), k - 1, use_c, use_dict)
    D = Plist[[1]]
    count = count + Plist[[2]]
    if (count >= rand_int) {
      count = count - Plist[[2]]
      multi = rep(k, f)
      break
    }
    f = f + 1
  }
  return(list(D, count, multi))
}  


#' multiplicity method of generating uniform random partitions of Q having N
#' parts.
#'
#' @param part a vector to hold the partition
#' @param q the total sum of the partition
#' @param D a dictionary for the number of partitions of Q having N or less
#'        parts (or N or less as the largest part), i.e. P(Q, Q + N).        
#' @param rand_int random integer
#' @param use_c boolean if TRUE then compiled c code is used
#' @export
#' @examples
#' multiplicity(c(5, 4), 4, hash(), 1, T, T) 
multiplicity =  function(part, q, D, rand_int, use_c, use_dict){
  while (q > 0) {
    multi = NULL
    if (!is.null(part)) {
      x = min(part)
    } 
    else {
      x = q
    }
    for (k in x:1) { # start with largest k
      Plist = P(D, q, k, use_c, use_dict) # number of partitions of q having k or less as the largest part
      D = Plist[[1]]
      count = Plist[[2]]
      if (count == rand_int & rand_int == 1) {
        multi = rep(1, q)
        q = 0
        break
      }   
      if (count < rand_int) { # k has been found
        k = k + 1
        Mlist = get_multiplicity(q, k, D, rand_int, count, use_c) # now, find how many times k occurs, i.e. the multiplicity of k 
        D = Mlist[[1]]
        count = Mlist[[2]]
        multi = Mlist[[3]]
        break
      }
    }
    q = q - sum(multi)
    part = c(part, multi)
    rand_int = rand_int - count
  }  
  part = conjugate(part)
  return(part)
}
