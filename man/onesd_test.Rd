% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{onesd_test}
\alias{onesd_test}
\title{Computes the one-sample log Bayes factor in favour of H1}
\usage{
onesd_test(
  n,
  s,
  popsd,
  u = 0.5,
  null_interval = NULL,
  alternative_interval = c(0, Inf),
  nonoverlapping_interval = FALSE,
  logarithm = TRUE
)
}
\arguments{
\item{n}{sample size}

\item{s}{sample standard deviation (with n - 1 as denominator)}

\item{popsd}{population standard deviation we test against}

\item{u}{parameter of the prior}

\item{null_interval}{interval for the null hypothesis (e.g., c(0.9, 1.1))}

\item{alternative_interval}{interval for the alternative hypothesis (e.g., c(1, Inf) and c(0, 1) give directed tests)}

\item{nonoverlapping_interval}{a logical specifying whether the intervals are non-overlapping (if TRUE, ignores alternative_interval and uses complement of null_interval)}

\item{logarithm}{a logical specifying whether the log should be taken}
}
\value{
The one-sample log Bayes factor in favour of H1 with \phi = popsd / s
}
\description{
Computes the one-sample log Bayes factor in favour of H1
}
\examples{
onesd_test(100, 1, 1, 0.50)

# one-sided test
onesd_test(100, 1, 1, 0.50, alternative_interval = c(1, Inf))

# Interval Bayes factor (overlapping)
onesd_test(100, 1, 1, 0.50, null_interval = c(0.9, 1.1))

# Interval Bayes factor (non-overlapping)
onesd_test(100, 1, 1, 0.50, null_interval = c(0.9, 1.1), nonoverlapping_interval = TRUE)

# one-sided interval Bayes factor
onesd_test(100, 1, 1, 0.50,
          alternative_interval = c(1.1, Inf),
          null_interval = c(0.9, 1.1))
}
