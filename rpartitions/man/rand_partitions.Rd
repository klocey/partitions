\name{rand_partitions}
\alias{rand_partitions}
\title{Generate uniform random partitions of Q having N parts.}
\usage{
  rand_partitions(Q, N, sample_size, method = "best",
    D = hash(), zeros = FALSE, use_c = TRUE,
    use_hash = FALSE)
}
\arguments{
  \item{Q}{Total sum across parts}

  \item{N}{Number of parts to sum over}

  \item{sample_size}{number of random partitions to
  generate}

  \item{method}{: method to use for generating the
  partition, options include: 'bottom_up', 'top_down',
  'divide_and_conquer', 'multiplicity', and 'best'.
  Defaults to 'best'}

  \item{D}{a dictionary for the number of partitions of Q
  having N or less parts (or N or less as the largest
  part), i.e. P(Q, Q + N). Defaults to a blank dictionary.}

  \item{zeros}{boolean if True partitions can have zero
  values, if False partitions have only positive values,
  defaults to False}

  \item{use_c}{boolean if TRUE then compiled c code is
  used, defaults to TRUE}

  \item{use_hash}{boolean, if TRUE then a hash table is
  used, defaults to FALSE}
}
\value{
  A matrix where each column is a random partition
}
\description{
  Generate uniform random partitions of Q having N parts.
}
\note{
  method 'best' attempts to use the values of Q and N to
  infer what the fastest method to compute the partition.

  if zeros are allowed, then we must ask whether Q >= N. if
  not, then the total Q is partitioned among a greater
  number of parts than there are, say, individuals. In
  which case, some parts must be zero. A random partition
  would then be any random partition of Q with zeros
  appended at the end. But, if Q >= N, then Q is
  partitioned among less number of parts than there are
  individuals. In which case, a random partition would be
  any random partition of Q having N or less parts.
}
\examples{
rand_partitions(100, 10, 5)
}

