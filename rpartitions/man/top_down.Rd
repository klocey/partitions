\name{top_down}
\alias{top_down}
\title{Top down method of generating uniform random partitions of Q having N parts.}
\usage{
  top_down(part, Q, D, rand_int, use_c, use_hash)
}
\arguments{
  \item{part}{a list to hold the partition}

  \item{Q}{the total sum of the partition}

  \item{D}{a dictionary for the number of partitions of Q
  having N or less parts (or N or less as the largest
  part), i.e. P(Q + N, N).}

  \item{rand_int}{a number representing a member of the
  feasible set}

  \item{use_c}{boolean if TRUE then compiled c code is
  used}

  \item{use_hash}{boolean, if TRUE then a hash table is
  used}
}
\description{
  Top down method of generating uniform random partitions
  of Q having N parts.
}
\examples{
top_down(c(5, 4), 4, hash(), 1, TRUE, FALSE)
}

