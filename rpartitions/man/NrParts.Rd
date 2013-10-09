\name{NrParts}
\alias{NrParts}
\title{Find the number of partitions for a given total Q and number of parts N.
Recoded and modified from GAP source code: www.gap-system.org}
\usage{
  NrParts(Q, N = NULL, use_c = TRUE)
}
\arguments{
  \item{Q}{Total sum}

  \item{N}{Number of items to sum across}

  \item{use_c}{logical, defaults to TRUE, the number of
  partitions is computed in c}
}
\description{
  Find the number of partitions for a given total Q and
  number of parts N. Recoded and modified from GAP source
  code: www.gap-system.org
}
\examples{
NrParts(100, 10)
}

