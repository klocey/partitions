\name{get_multiplicity}
\alias{get_multiplicity}
\title{Find the number of times a value k occurs in a partition that is being
generated at random by the multiplicity() function. The resulting
multiplicity is then passed back to the multiplicity() function along with
an updated value of count and an updated dictionary D}
\usage{
  get_multiplicity(Q, k, D, rand_int, count, use_c,
    use_hash)
}
\arguments{
  \item{Q}{the total sum of the partition}

  \item{k}{the size of the largest (and also first) part}

  \item{D}{a dictionary for the number of partitions of Q
  having N or less parts (or N or less as the largest
  part), i.e. P(Q, Q + N).}

  \item{rand_int}{the random integer}

  \item{count}{count < rand_int}

  \item{use_c}{boolean if TRUE then compiled c code is
  used}

  \item{use_hash}{boolean, if TRUE then a hash table is
  used}
}
\description{
  Find the number of times a value k occurs in a partition
  that is being generated at random by the multiplicity()
  function. The resulting multiplicity is then passed back
  to the multiplicity() function along with an updated
  value of count and an updated dictionary D
}
\examples{
get_multiplicity(10, 5, hash(), 3, 2, TRUE, FALSE)
}

