% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{twosd_test}
\alias{twosd_test}
\title{Computes the two-sample log Bayes factor in favour of H1}
\usage{
twosd_test(
  n1,
  n2,
  sd1,
  sd2,
  u = 0.5,
  null_interval = NULL,
  alternative_interval = c(0, Inf),
  nonoverlapping_interval = FALSE,
  logarithm = TRUE
)
}
\arguments{
\item{n1}{sample size of group 1}

\item{n2}{sample size of group 2}

\item{sd1}{sample standard deviation of group 1 (with n1 - 1 as denominator)}

\item{sd2}{sample standard deviation of group 2 (with n2 - 1 as denominator)}

\item{u}{parameter of the prior}

\item{null_interval}{interval for the null hypothesis (e.g., c(0.9, 1.1))}

\item{alternative_interval}{interval for the alternative hypothesis (e.g., c(1, Inf) and c(0, 1) give directed tests)}

\item{nonoverlapping_interval}{a logical specifying whether the intervals are non-overlapping (if TRUE, ignores alternative_interval and uses complement of null_interval)}

\item{logarithm}{a logical specifying whether the log should be taken}
}
\value{
The two-sample log Bayes factor in favour of H1 with \phi = sd2 / sd1
}
\description{
Computes the two-sample log Bayes factor in favour of H1
}
\examples{

# Interval null Bayes factor
twosd_test(100, 200, 1, 2, 4.5, null_interval = c(0.9, 1.1))

# One-sided tests
twosd_test(100, 200, 1, 2, 4.5, alternative_interval = c(1, Inf))
twosd_test(100, 200, 1, 2, 4.5, alternative_interval = c(0, 1))

# Interval Bayes factor (overlapping)
twosd_test(100, 100, 1, 1, 0.50, null_interval = c(0.9, 1.1))

# Interval Bayes factor (non-overlapping)
twosd_test(100, 100, 1, 1, 0.50, null_interval = c(0.9, 1.1), nonoverlapping_interval = TRUE)

# One-sided interval Bayes factor
twosd_test(100, 200, 1, 2, 4.5,
           alternative_interval = c(1.1, Inf),
           null_interval = c(0.9, 1.1))
}
