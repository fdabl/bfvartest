#' Computes the one-sample log Bayes factor in favour of H1
#'
#' @param n sample size
#' @param s sample standard deviation (1/n * \sqrt{\sum_i (x_i - \mu)^2})
#' @param popsd population standard deviation we test against
#' @param alpha prior parameter
#' @param alternative_interval interval for the alternative hypothesis (e.g., c(1, Inf) and c(0, 1) give directed tests)
#' @param null_interval interval for the null hypothesis (e.g., c(0.9, 1.1))
#' @return the one-sample log Bayes factor in favour of H1 with \delta = popsd / s
#' @examples
#' one_sample(100, 1, 1, 0.50)
#' one_sample(100, 1, 1, 0.50, alternative_interval = c(1, Inf)) # one-sided test
#' one_sample(100, 1, 1, 0.50, alternative_interval = c(0.9, 1.1)) # interval Bayes factor
#' one_sample(100, 1, 1, 0.50, alternative_interval = c(1.1, Inf), null_interval = c(0.9, 1.1)) # one-sided interval Bayes factor
one_sample <- function(n, s, popsd, alpha = 0.50, alternative_interval = c(0, Inf), null_interval = NULL, logarithm = TRUE) {
  .check_interval_input(alternative_interval, null_interval)

  s2 <- s^2
  popvar <- popsd^2
  tau0 <- 1 / popvar

  if (is.null(null_interval)) {
    logml0 <- (n - 1)/2 * log(tau0) - 0.50 * tau0 * n * s2 # correct!
  } else {
    logml0 <- .compute_logml_restr_k1(n, s2, popsd, interval = null_interval, alpha = alpha)
  }

  logml1 <- .compute_logml_restr_k1(n, s2, popsd, interval = alternative_interval, alpha = alpha)
  ifelse(logarithm, logml1 - logml0, exp(logml1 - logml0))[[1]]
}


#' Computes the two-sample log Bayes factor in favour of H1
#'
#' @param n1 sample size of group 1
#' @param n2 sample size of group 2
#' @param s1 observed sum of squares of group 1
#' @param s2 observed sum of squares of group 2
#' @param alpha prior parameter
#' @param alternative_interval interval for the alternative hypothesis (e.g., c(1, Inf) and c(0, 1) give directed tests)
#' @param null_interval interval for the null hypothesis (e.g., c(0.9, 1.1))
#' @return the two-sample log Bayes factor in favour of H1
#' @examples
#' two_sample(100, 200, 1, 2, 4.5, alternative_interval = c(1, Inf)) # one-sided test
#' two_sample(100, 200, 1, 2, 4.5, alternative_interval = c(0, 1)) # one-sided test
#' two_sample(100, 200, 1, 2, 4.5, null_interval = c(0.9, 1.1)) # interval Bayes factor
#' two_sample(100, 200, 1, 2, 4.5, alternative_interval = c(1.1, Inf), null_interval = c(0.9, 1.1)) # one-sided interval Bayes factor
two_sample <- function(n1, n2, s1, s2, alpha = 0.50, alternative_interval = c(0, Inf), null_interval = NULL, logarithm = TRUE) {
  .check_interval_input(alternative_interval, null_interval)

  if (is.null(null_interval)) {
    logml0 <- ((2 - n1 - n2)/2) * log(n1*s1 + n2*s2)
  } else {
    logml0 <- .compute_logml_restr_k2(n1, n2, s1, s2, interval = null_interval, alpha = alpha)
  }

  logml1 <- .compute_logml_restr_k2(n1, n2, s1, s2, interval = alternative_interval, alpha = alpha)
  val <- ifelse(logarithm, logml1 - logml0, exp(logml1 - logml0))[[1]]
  suppressWarnings(Rmpfr::asNumeric(val))
}


#' Computes the k-sample log Bayes factor for all hypotheses
#'
#' @param hyp vector of hypotheses
#' @param ns vector of sample sizes
#' @params ss a vector containing sample sum of squares
#' @param alpha prior parameter
#' @param ... arguments to rstan::sampling
#' @return the log Bayes factors of all hypotheses against each other
#' @examples
#' ss <- c(1, 2, 3)
#' ns <- c(100, 100, 100)
#' hyp <- c('1=2=3', '1,2,3', '1>2>3')
#' k_sample(hyp, ns, ss, 0.50)
k_sample <- function(hyp, ns, ss, alpha = 0.50, logarithm = TRUE, compute_ml = TRUE, priors_only = FALSE, ...) {

  len <- c(length(hyp), length(ns), length(ss))

  if (len[1] < 2) {
    stop('Need at least two hypotheses to compare!')
  }
  if (len[2] != len[3]) {
    stop('Need information about sample size and observed sum of squares for all groups!')
  }

  res <- list()

  for (h in hyp) {
    res[[h]] <- .create_levene_object(
      h, ns, ss, alpha, priors_only = priors_only,
      compute_ml = ifelse(priors_only, !priors_only, compute_ml), ...
    )
  }

  # Add Bayes factors
  if (!priors_only && compute_ml) {
    nr_hyp <- length(hyp)
    BF_matrix <- diag(0, nr_hyp)
    colnames(BF_matrix) = rownames(BF_matrix) = hyp

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
