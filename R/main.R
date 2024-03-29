#' Computes the one-sample log Bayes factor in favour of H1
#'
#' @export
#' @param n sample size
#' @param s sample standard deviation (with n - 1 as denominator)
#' @param popsd population standard deviation we test against
#' @param u parameter of the prior
#' @param alternative_interval interval for the alternative hypothesis (e.g., c(1, Inf) and c(0, 1) give directed tests)
#' @param null_interval interval for the null hypothesis (e.g., c(0.9, 1.1))
#' @param nonoverlapping_interval a logical specifying whether the intervals are non-overlapping (if TRUE, ignores alternative_interval and uses complement of null_interval)
#' @param logarithm a logical specifying whether the log should be taken
#' @return The one-sample log Bayes factor in favour of H1 with \eqn{\phi = popsd / s}
#' @examples
#' onesd_test(100, 1, 1, 0.50)
#'
#' # one-sided test
#' onesd_test(100, 1, 1, 0.50, alternative_interval = c(1, Inf))
#'
#' # Interval Bayes factor (overlapping)
#' onesd_test(100, 1, 1, 0.50, null_interval = c(0.9, 1.1))
#'
#' # Interval Bayes factor (non-overlapping)
#' onesd_test(100, 1, 1, 0.50, null_interval = c(0.9, 1.1), nonoverlapping_interval = TRUE)
#'
#' # one-sided interval Bayes factor
#' onesd_test(100, 1, 1, 0.50,
#'           alternative_interval = c(1.1, Inf),
#'           null_interval = c(0.9, 1.1))
onesd_test <- function(
  n, s, popsd, u = 0.50, null_interval = NULL,
  alternative_interval = c(0, Inf), nonoverlapping_interval = FALSE, logarithm = TRUE
  ) {
  .check_interval_input(alternative_interval, null_interval)

  # Convert to sample sum of squares
  s2 <- s^2
  popvar <- popsd^2
  tau0 <- 1 / popvar

  # degrees of freedom
  v <- n - 1

  if (is.null(null_interval)) {
    logml0 <- v / 2 * log(tau0) - tau0 * v * s2 / 2
  } else {
    logml0 <- .compute_logml_restr_k1(v, s2, popsd, interval = null_interval, u = u)
  }

  if (nonoverlapping_interval) {
    if (is.null(null_interval)) {
      stop('Non-overlapping interval Bayes factors require a specification of the null interval')
    }

    logml1 <- .compute_logml_restr_k1(
      v, s2, popsd, interval = null_interval, u = u, nonoverlapping_interval = TRUE
    )

  } else {

    logml1 <- .compute_logml_restr_k1(v, s2, popsd, interval = alternative_interval, u = u)
  }


  ifelse(logarithm, logml1 - logml0, exp(logml1 - logml0))
}


#' Computes the two-sample log Bayes factor in favour of H1
#'
#' @export
#' @param n1 sample size of group 1
#' @param n2 sample size of group 2
#' @param sd1 sample standard deviation of group 1 (with n1 - 1 as denominator)
#' @param sd2 sample standard deviation of group 2 (with n2 - 1 as denominator)
#' @param u parameter of the prior
#' @param null_interval interval for the null hypothesis (e.g., c(0.9, 1.1))
#' @param alternative_interval interval for the alternative hypothesis (e.g., c(1, Inf) and c(0, 1) give directed tests)
#' @param nonoverlapping_interval a logical specifying whether the intervals are non-overlapping (if TRUE, ignores alternative_interval and uses complement of null_interval)
#' @param logarithm a logical specifying whether the log should be taken
#' @return The two-sample log Bayes factor in favour of H1 with \eqn{\phi = sd2 / sd1}
#' @examples
#'
#' # Interval null Bayes factor
#' twosd_test(100, 200, 1, 2, 4.5, null_interval = c(0.9, 1.1))
#'
#' # One-sided tests
#' twosd_test(100, 200, 1, 2, 4.5, alternative_interval = c(1, Inf))
#' twosd_test(100, 200, 1, 2, 4.5, alternative_interval = c(0, 1))
#'
#' # Interval Bayes factor (overlapping)
#' twosd_test(100, 100, 1, 1, 0.50, null_interval = c(0.9, 1.1))
#'
#' # Interval Bayes factor (non-overlapping)
#' twosd_test(100, 100, 1, 1, 0.50, null_interval = c(0.9, 1.1), nonoverlapping_interval = TRUE)
#'
#' # One-sided interval Bayes factor
#' twosd_test(100, 200, 1, 2, 4.5,
#'            alternative_interval = c(1.1, Inf),
#'            null_interval = c(0.9, 1.1))
twosd_test <- function(
  n1, n2, sd1, sd2, u = 0.50, null_interval = NULL,
  alternative_interval = c(0, Inf), nonoverlapping_interval = FALSE, logarithm = TRUE) {
  .check_interval_input(alternative_interval, null_interval)

  # convert to sample sum of squares
  s1 <- sd1^2
  s2 <- sd2^2

  # degrees of freedom
  v1 <- n1 - 1
  v2 <- n2 - 1
  v <- v1 + v2

  if (is.null(null_interval)) {
    logml0 <- (-v / 2) * log(v1 * s1 + v2 * s2) # proportional to logml0
  } else {
    logml0 <- .compute_logml_restr_k2(v1, v2, s1, s2, interval = null_interval, u = u) # proportional to logml1
  }

  if (nonoverlapping_interval) {
    if (is.null(null_interval)) {
      stop('Non-overlapping interval Bayes factors require a specification of the null interval')
    }

    logml1 <- .compute_logml_restr_k2(
      v1, v2, s1, s2, interval = null_interval, nonoverlapping_interval = TRUE, u = u
      )

  } else {

    logml1 <- .compute_logml_restr_k2(v1, v2, s1, s2, interval = alternative_interval, u = u)
  }


  val <- ifelse(logarithm, logml1 - logml0, exp(logml1 - logml0))[[1]]
  suppressWarnings(Rmpfr::asNumeric(val))
}


#' Computes the posterior density of \eqn{\phi} for the K = 2 case
#'
#' @export
#' @param x numerical value
#' @param n1 sample size of group 1
#' @param n2 sample size of group 2
#' @param sd1 sample standard deviation of group 1 (with n1 - 1 as denominator)
#' @param sd2 sample standard deviation of group 2 (with n2 - 1 as denominator)
#' @param u parameter of the prior
#' @param logarithm a logical specifying whether the log should be taken
#' @return The (log) density at x
#' @examples
#'
#' dphi2(seq(0, 2, .01), 100, 100, 1, 1)
dphi2 <- function(x, n1, n2, sd1, sd2, u = 0.50, logarithm = FALSE) {

  # convert to sample sum of squares
  s1 <- sd1^2
  s2 <- sd2^2

  # degrees of freedom
  v1 <- n1 - 1
  v2 <- n2 - 1
  v <- v1 + v2

  # Normalizing constant
  Z <- lbeta(v1/2 + u, v2/2 + u) +
    log(.Gauss2F1(2*u, v2/2 + u, v/2 + 2*u, 1 - (v1 * s1) / (v2 * s2))) +
    (-v1/2 + u) * log((v1 * s1) / (v2 * s2))

  # Unnormalized density
  val <- log(2) + (v1 + 2 * u - 1) * log(x) + (-2 * u) * log(1 + x^2) +
         (-v / 2) * log(x^2 * v1 * s1 / (v2 * s2) + 1)

  if (logarithm) {
    return(val - Z)
  } else {
    return(exp(val - Z))
  }
}


#' Computes the k-sample log Bayes factor for all hypotheses
#'
#' @export
#' @param hyp vector of hypotheses
#' @param ns vector of sample sizes
#' @param sds vector containing the sample standard deviations (with ns - 1 as denominators)
#' @param u numeric specifying the value of the parameters of the Dirichlet prior
#' @param logarithm logical specifying whether the log should be taken
#' @param compute_ml logical specifying whether the marginal likelihood should be computed
#' @param priors_only logical specifying whether we should only sample from the prior
#' @param ... arguments to rstan::sampling
#' @return A list of 'bfvar' objects, which include stanfit objects, log marginal likelihoods, and pairwise Bayes factors
#' @examples
#' ss <- c(3, 2, 1)
#' ns <- c(100, 100, 100)
#' hyp <- c('1=2=3', '1,2,3', '1>2>3', '1>2=3', '1>2,3')
#' ksd_test(hyp, ns, ss, 0.50)
ksd_test <- function(hyp, ns, sds, u = 0.50, logarithm = TRUE, compute_ml = TRUE, priors_only = FALSE, ...) {
  .check_user_input(hyp, ns, sds)

  # Translate hypothesis on variance / standard deviation into hypothesis on precision
  hyp2 <- gsub('>', '<', hyp)
  nr_hyp <- length(hyp)
  res <- list()

  for (i in seq(nr_hyp)) {
    res[[hyp[i]]] <- .create_bfvar_object(
      hyp2[i], ns, sds, u, priors_only = priors_only,
      compute_ml = ifelse(priors_only, FALSE, compute_ml), ...
    )
  }

  # Add Bayes factors
  if (!priors_only && compute_ml) {
    BF_matrix <- matrix(0, nrow = nr_hyp, ncol = nr_hyp, dimnames = list(hyp, hyp))

    for (i in seq(nr_hyp)) {
      for (j in seq(nr_hyp)) {
        if (i != j) {
          hi <- hyp[i]
          hj <- hyp[j]
          BF_matrix[i, j] <- res[[hi]]$logml - res[[hj]]$logml
        }
      }
    }

    if (logarithm) {
      res[['BF']] <- BF_matrix
    } else {
      res[['BF']] <- exp(BF_matrix)
    }
  }

  res
}
