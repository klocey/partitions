\name{conjugate}
\alias{conjugate}
\title{Find the conjugate of an integer partition
Recoded (orginally on 24-Apr-2013) from the Sage source code:
http://www.sagenb.org/src/combinat/partition.py}
\usage{
  conjugate(partition, use_c = TRUE)
}
\arguments{
  \item{partition}{a vector that represents an integer
  partition}

  \item{use_c}{logical, defaults to TRUE, the conjugate is
  computed in c}
}
\description{
  Find the conjugate of an integer partition Recoded
  (orginally on 24-Apr-2013) from the Sage source code:
  http://www.sagenb.org/src/combinat/partition.py
}
\examples{
conjugate(c(3,3,1,1), FALSE)
}

