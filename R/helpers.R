# Checks whether the intervals for the null and the alternative are correctly specified
.check_interval_input <- function(alternative_interval, null_interval) {
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
}


# Computes the log marginal likelihood for 1 = 2 hypothesis
.compute_logml_restr_k1 <- function(n, s2, popsd, interval, alpha = 0.50) {
  popvar <- popsd^2
  tau0 <- 1 / popvar

  scaled_lbetaprime <- function(tau, tau0, alpha) {
    value <- (alpha - 1) * log(tau / tau0) - 2*alpha *
             log(1 + tau / tau0) - log(tau0) - lbeta(alpha, alpha)
    value
  }

  lo <- interval[1]
  hi <- interval[2]
  Z <- integrate(function(tau) exp(scaled_lbetaprime(tau, tau0, alpha)), 0, hi)$value -
       ifelse(
         lo == 0, 0, integrate(function(tau) exp(scaled_lbetaprime(tau, tau0, alpha)), 0, lo)$value
       )

  value <- integrate(function(tau) {
    llh <- (n - 1)/2 * log(tau) - 0.5 * tau * n * s2
    lprior <- scaled_lbetaprime(tau, tau0, alpha)

    exp(llh + lprior - log(Z))
  }, lo, hi)$value

  log(value)
}


# Computes the log marginal likelihood for K = 2 hypothesis
.compute_logml_restr_k2 <- function(n1, n2, s1, s2, interval, alpha = 0.50) {
  n <- n1 + n2
  to_rho <- function(delta) ifelse(delta == Inf, 1, delta^2 / (1 + delta^2))
  lo <- to_rho(interval[1])
  hi <- to_rho(interval[2])
  Z <- pbeta(hi, alpha, alpha) - pbeta(lo, alpha, alpha)

  # Rmpfr provides arbitrary precision floating point arithmetic
  value <- Rmpfr::integrateR(function(rho) {
    rho <- Rmpfr::mpfr(rho, 100)
    llh <- ((n1 - 1)/2 + alpha - 1) * log(rho) +
           ((n2 - 1)/2 + alpha - 1) * log(1 - rho) +
           ((2 - n)/2) * log(rho*n1*s1 + (1 - rho)*n2*s2)

    exp(llh - log(Z) - lbeta(alpha, alpha))
  }, lo, hi)$value

  log(value)
}


# computes the prior probability with which \rho follows the hypothesis
.compute_prior_restr <- function(hyp) {
  s <- strsplit(hyp, '<')[[1]]
  n <- length(s)

  ss <- strsplit(hyp, '')[[1]]
  ngroups <- as.numeric(ss[length(ss)])
  num <- factorial(ngroups)

  for (i in seq(n)) {
    num <- num / factorial(length(strsplit(s[i], ',')[[1]]))
  }

  1 / num
}


# '1,2<3=4,5<6' -> '1,2<3,4<5'
.reduce_hyp_equals <- function(string) {
  n <- length(string)
  res <- c()
  j <- 1
  for (i in seq(n)) {
    si <- strsplit(string[i], ',')[[1]]
    ni <- length(si)

    r <- j
    j <- j + 1
    if (ni > 1) {
      for (k in seq(2, ni)) {
        r <- paste0(r, ',', j)
        j <- j + 1
      }
    }

    res <- c(res, r)
  }

  res
}


# Parses user input such as '1,2<3=4,5<6'
# We reduce the dimensionality so that it results in '1,2<3,4<5'
.create_hyp_fn <- function(hyp) {
  s <- strsplit(hyp, '<')[[1]]
  sr <- .reduce_hyp_equals(s)
  dig <- gsub('[^0-9.-]', ',', sr)
  fn_string <- 'function(x) '

  # add order-constraints
  for (i in seq(2, length(dig))) {
    left <- paste0('max(x[c(', dig[i-1], ')])')
    right <- paste0('min(x[c(', dig[i], ')])')
    fn_string <- paste0(fn_string, left, ' < ', right, ifelse(i == length(dig), '', ' && '))
  }

  fn_string
}


# Checks whether the user input to k_sample is valid
.check_user_input <- function(hyp, ns, ss) {
  len <- c(length(hyp), length(ns), length(ss))

  if (len[1] < 2) {
    stop('Need at least two hypotheses to compare!')
  }

  if (len[2] != len[3]) {
    stop('Need information about sample size and observed sum of squares for all groups!')
  }

  for (i in seq(length(hyp))) {
    check <- gsub('[0-9>=, ]', '.', hyp[i])

    if (check != gsub('.', '.', hyp[i])) {
      stop('Hypothesis ', i, ' is misspecified. Only relations of the form [>,=] are allowed!')
    }
  }
}


# Checks whether hypothesis only includes ordinal and no constraints (e.g., '1<2=3<4')
.is_only_ordered_and_equal <- function(hyp) {
  s <- strsplit(hyp, '>')[[1]]

  check <- sapply(s, function(si) {
    length(strsplit(si, ',')[[1]]) == 1
  })

  all(check)
}


# Checks whether hypothesis is all equal (e.g., '1=2=3=4')
.is_allequal <- function(hyp) {
  s <- strsplit(hyp, '=')[[1]]
  all(nchar(s) < 2)
}


# Creates a list that is passed onto Stan
.prepare_standat <- function(hyp, ns, ss, a, priors_only = FALSE) {
  k <- length(ns)

  if (k != length(ss)) {
    stop('Need same amount of information for each group.')
  }

  check <- gsub('[0-9<=, ]', '.', hyp)
  rel <- strsplit(gsub('[0-9]', '', hyp), '')[[1]]
  rel[which(rel == '=')] <- 1
  rel[which(rel == '<')] <- 2
  rel[which(rel == ',')] <- 3
  rel <- as.numeric(rel)

  nr_equal <- sum(rel == 1)
  nr_ordered <- sum(rel == 2)
  strip_num <- strsplit(gsub('[<,=]', ' ', hyp), ' ')[[1]]
  strip_punct <- strsplit(gsub('[0-9]', ' ', hyp), ' +')[[1]]

  strip <- c()
  for (i in seq(length(strip_num) - 1)) {
    strip <- c(strip, strip_num[i], strip_punct[i+1])
  }
  strip <- c(strip, strip_num[length(strip_num)])

  index <- 1
  count <- 1
  index_vector <- numeric(k)

  for (i in seq(length(strip))) {
    if (strip[i] %in% c('<', ',')) {
      count <- count + 1
    }

    if (grepl('[1-9]', strip[i])) {
      index_vector[index] <- count
      index <- index + 1
    }
  }

  standat <- list(
    k = k, s2 = ss,
    N = ns, alpha = a,
    index_vector = index_vector,
    priors_only = priors_only,
    nr_equal = nr_equal, nr_ordered = nr_ordered + 1
  )

  # Otherwise Stan throws an error in the K = 2 case
  # dim(standat$relations) <- length(standat$relations)
  standat
}


# extend print function to only show parameters tau and rho
print.bfvar <- function(x) {
  pars <- c('tau', 'rho')
  print(x$fit, pars = pars)

  if (!is.null(x$logml)) {
    cat(paste0('\nLog Marginal Likelihood = ', round(x$logml, 3), '\n'))
  }
}


#' Estimates the model using Stan and computes the marginal likelihood
#'
#' @params hypothesis a string specifying the hypothesis
#' @params ss a vector containing sample sum of squares
#' @params ns a vector containing sample sizes
#' @params a a vector specifying the value of the parameters of the Dirichlet prior
#' @params compute_ml a logical specifying whether the marginal likelihood should be computes
#' @params priors_only a logical specifying whether we should only sample from the prior
#' @params precision a logical specifying whether parameterization in terms of precision or variance is used
#' @params ... arguments to rstan::sampling
#'
#' @returns an object of class 'bfvar', which is a stanfit object with a log marginal likelihood (if specified)
.create_bfvar_object <- function(hyp, ns, ss, a, compute_ml = TRUE, priors_only = FALSE, ...) {

  # translate hypothesis on variance / standard deviation into hypothesis on precision
  hyp <- gsub('>', '<', hyp)
  logml <- NULL
  standat <- .prepare_standat(hyp, ns, ss, a, priors_only = priors_only)

  if (!.is_only_ordered_and_equal(hyp) || .is_allequal(hyp)) {
    fit <- rstan::sampling(stanmodels$Mixed, data = standat, ...)

    if (compute_ml && !priors_only) {
      logml <- bridgesampling::bridge_sampler(fit)$logml
    }

  } else if (.is_only_ordered_and_equal(hyp)) {
    fit <- rstan::sampling(stanmodels$Ordered, data = standat, ...)

    if (compute_ml && !priors_only) {
      logml <- bridgesampling::bridge_sampler(fit)$logml
    }
  } else {
    fit <- rstan::sampling(stanmodels$Mixed, data = standat, ...)

    if (compute_ml && !priors_only) {
      rho <- rstan::extract(fit, 'rho')$rho
      hyp_fn <- eval(parse(text = .create_hyp_fn(hyp)))

      # Kluglist & Hoijtink (2005) trick
      BF_mixedequal_fullequal <- log(mean(apply(rho, 1, hyp_fn))) - log(.compute_prior_restr(hyp))
      logml_fullequal <- bridgesampling::bridge_sampler(fit)$logml
      logml <- BF_mixedequal_fullequal + logml_fullequal # this is logml_mixedequal
    }
  }

  res <- list('fit' = fit, 'logml' = logml)
  class(res) <- 'bfvar'
  res
}
