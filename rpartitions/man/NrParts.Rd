\name{NrParts}
\alias{NrParts}
\title{Find the number of partitions for a given total Q and number of parts N.}
\usage{
  NrParts(Q, N = NULL, use_c = TRUE)
}
\arguments{
  \item{Q}{Total sum}

  \item{N}{Number of items to sum across, if not specified
  than all possible values are considered}

  \item{use_c}{logical, defaults to TRUE, the number of
  partitions is computed in c}
}
\description{
  This function was recoded and modified from GAP source
  code: www.gap-system.org. Modifications for speed were
  based on the proposition that the number of partitions of
  Q having N parts is equal to the number of partitions of
  Q having N parts is equal to the number of partitions of
  Q - N, if N > Q/2 (for odd Q) or if N >= Q/2 (for even Q)
}
\examples{
NrParts(100)
NrParts(100, 10)
}

