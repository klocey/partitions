\name{P}
\alias{P}
\title{Number of partitions of q with k or less parts (or having k or less as the
largest part), i.e. P(q+k,k).}
\usage{
  P(D, q, k, use_c, use_dict)
}
\arguments{
  \item{D}{lookup table for numbers of partitions, P(q+k,k)
  values.}

  \item{q}{total sum of the set}

  \item{use_c}{logical, if TRUE the number of partitions is
  computed in c}
}
\description{
  Number of partitions of q with k or less parts (or having
  k or less as the largest part), i.e. P(q+k,k).
}
\examples{
P(list(), 100, 10, FALSE, FALSE)
}

