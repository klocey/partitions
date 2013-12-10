\name{multiplicity}
\alias{multiplicity}
\title{multiplicity method of generating uniform random partitions of Q having N
parts.}
\usage{
  multiplicity(part, Q, D, rand_int, use_c, use_hash)
}
\arguments{
  \item{part}{a vector to hold the partition}

  \item{Q}{the total sum of the partition}

  \item{D}{a dictionary for the number of partitions of Q
  having N or less parts (or N or less as the largest
  part), i.e. P(Q, Q + N).}

  \item{rand_int}{random integer}

  \item{use_c}{boolean if TRUE then compiled c code is
  used}

  \item{use_hash}{boolean, if TRUE then a hash table is
  used}
}
\description{
  multiplicity method of generating uniform random
  partitions of Q having N parts.
}
\examples{
multiplicity(c(5, 4), 4, hash(), 1, TRUE, FALSE)
}

