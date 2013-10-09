\name{divide_and_conquer}
\alias{divide_and_conquer}
\title{Divide and conquer method of generating uniform random partitions of Q
having N parts.}
\usage{
  divide_and_conquer(part, q, N, D, rand_int, use_c,
    use_dict)
}
\arguments{
  \item{part}{a list to hold the partition}

  \item{q}{the total sum of the partition}

  \item{N}{Number of parts to sum over}

  \item{D}{a dictionary for the number of partitions of Q
  having N or less parts (or N or less as the largest
  part), i.e. P(Q, Q + N).}

  \item{rand_int}{the random integer}

  \item{use_c}{boolean if TRUE then compiled c code is
  used}

  \item{use_dict}{boolean, if TRUE then hash dictionary is
  used}
}
\description{
  Divide and conquer method of generating uniform random
  partitions of Q having N parts.
}
\examples{
divide_and_conquer(c(5, 4), 4, hash(), 1, T, T)
}

