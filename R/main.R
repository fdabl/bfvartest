#' Computes the one-sample log Bayes factor in favour of H1
#'
#' @param n sample size
#' @param s2 observed sum of squares
#' @param popsd population standard deviation we test against
#' @param alpha prior parameter
#' @return the one-sample log Bayes factor in favour of H1
#' @examples
#' one_sample(100, 1, 1, 4.50)
one_sample <- function(n, s2, popsd, alpha = 0.50) {
  popvar <- popsd^2
  tau0 <- 1 / popvar

  value <- integrate(function(tau) {
    term <- ((n - 1)/2 + alpha) * log(tau / tau0) - log(tau) -2*alpha * log(1 + tau/tau0) - 0.5*n*s2*(tau - tau0)
    exp(term)
  }, 0, Inf)$value / beta(alpha, alpha)

  log(value)
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
#' two_sample(100, 200, 1, 2, 4.5, alternative_interval = c(1, Inf), null_interval = c(0.9, 1.1)) # one-sided interval Bayes factor
two_sample <- function(n1, n2, s1, s2, alpha = 0.50, alternative_interval = c(0, Inf), null_interval = NULL) {

  alt_condition <- (all(is.numeric(alternative_interval)) &&
                    alternative_interval[2] > alternative_interval[1] &&
                    all(alternative_interval >= 0))

  if (!alt_condition) {
    stop('Something is off with your specification of alternative_interval')
  }

  if (!is.null(null_interval)) {
    is_positive <- function(x) x >= 0

    null_condition <- (all(is.numeric(null_interval)) &&
                       all(null_interval >= 0) &&
                       null_interval[2] > null_interval[1] &&
                       all(is_positive(null_interval)))

    if (!null_condition) {
      stop('Something is off with your specification of null_interval')
    }
  }

  if (is.null(null_interval)) {
    logml0 <- ((2 - n1 - n2)/2) * log(n1*s1 + n2*s2)
  } else {
    logml0 <- .compute_logml_restr(n1, n2, s1, s2, interval = null_interval, alpha = alpha)
  }

  logml1 <- .compute_logml_restr(n1, n2, s1, s2, interval = alternative_interval, alpha = alpha)
  logml1 - logml0
}


.compute_logml_restr <- function(n1, n2, s1, s2, interval, alpha = 0.50) {
  n <- n1 + n2
  to_rho <- function(delta) ifelse(delta == Inf, 1, delta^2 / (1 + delta^2))
  lo <- to_rho(interval[1])
  hi <- to_rho(interval[2])
  Z <- pbeta(hi, alpha, alpha) - pbeta(lo, alpha, alpha)

  # Rmpfr provides arbitrary precision floating point arithmetic
  value <- Rmpfr::integrateR(function(rho) {
    rho <- Rmpfr::mpfr(rho, 100)
    llh <- ((n1 - 1)/2 + alpha - 1) * log(rho) + ((n2 - 1)/2 + alpha - 1) * log(1 - rho) + ((2 - n)/2) * log(rho*n1*s1 + (1 - rho)*n2*s2)

    exp(llh - log(Z) - lbeta(alpha, alpha))
  }, lo, hi)$value

  log(value)
}
