\name{bottom_up}
\alias{bottom_up}
\title{Bottom up method of generating uniform random partitions of Q having N parts.}
\usage{
  bottom_up(part, q, D, rand_int, use_c, use_dict)
}
\arguments{
  \item{part}{a list to hold the partition}

  \item{q}{the total sum of the partition}

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
  Bottom up method of generating uniform random partitions
  of Q having N parts.
}
\examples{
bottom_up(c(5, 4), 4, hash(), 1, T, T)
}

