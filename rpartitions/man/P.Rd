\name{P}
\alias{P}
\title{Number of partitions of Q with k or less parts.}
\usage{
  P(D, Q, k, use_c, use_hash)
}
\arguments{
  \item{D}{lookup table for numbers of partitions of Q
  having k or less parts (or k or less as the largest
  part), i.e. P(Q, Q + k)}

  \item{Q}{total (i.e., sum across all k or n parts)}

  \item{k}{the number of parts and also the size of the
  largest part (congruency)}

  \item{use_c}{boolean, if TRUE the number of partitions is
  computed in c}

  \item{use_hash}{boolean, if TRUE then a hash table is
  used instead of R's native list to store the information}
}
\value{
  a two element list, the first element is D the lookup
  table and the second element is the number of partitions
  for the specified Q and k value.
}
\description{
  This function was derived using the following theorem and
  proposition. The number of partitions of Q with k or less
  parts equals the number of partitions of Q with k or less
  as the largest part (see Bona 2006). This is a
  mathematical symmetry, i.e. congruency. Additionally, the
  number of partitions of Q with k or less parts equals the
  number of partitions of Q+k with k as the largest part
  when k>0, i.e. P(Q + k, k). We do not have a source for
  this proposition, but it can be shown when enumerating
  the entire feasible set or using the Sage computing
  enviornment
}
\examples{
P(list(), 100, 10, FALSE, FALSE)
}
\references{
  Bona, M. (2006). A Walk Through Combinatorics: An
  Introduction to Enumeration and Graph Theory. 2nd Ed.
  World Scientific Publishing Co. Singapore.
}

