\name{divide_and_conquer}
\alias{divide_and_conquer}
\title{Divide and conquer method of generating uniform random partitions of Q
having N parts.}
\usage{
  divide_and_conquer(part, Q, N, D, rand_int, use_c,
    use_hash)
}
\arguments{
  \item{part}{a list to hold the partition}

  \item{Q}{the total sum of the partition}

  \item{N}{Number of parts to sum over}

  \item{D}{a dictionary for the number of partitions of Q
  having N or less parts (or N or less as the largest
  part), i.e. P(Q, Q + N).}

  \item{rand_int}{a number representing a member of the
  feasible set}

  \item{use_c}{boolean if TRUE then compiled c code is
  used}

  \item{use_hash}{boolean, if TRUE then a hash table is
  used}
}
\description{
  Divide and conquer method of generating uniform random
  partitions of Q having N parts.
}
\examples{
divide_and_conquer(c(5, 4), 5, 4, hash(), 2, TRUE, FALSE)
}

