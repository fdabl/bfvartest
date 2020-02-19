# Computes the log marginal likelihood for 1 = 2 hypothesis
.compute_logml_restr_k1 <- function(n, s2, popsd, interval, alpha = 0.50) {
  popvar <- popsd^2
  tau0 <- 1 / popvar

  # Prior distribution, centered at \tau_0
  scaled_lbetaprime <- function(tau, tau0, alpha) {
    value <- (alpha - 1) * log(tau / tau0) - 2*alpha *
             log(1 + tau / tau0) - log(tau0) - lbeta(alpha, alpha)
    value
  }

  # Change normalizing constant when prior is restricted
  lo <- interval[1]
  hi <- interval[2]
  hiZ <- integrate(function(tau) exp(scaled_lbetaprime(tau, tau0, alpha)), 0, hi)$value
  loZ <- ifelse(lo == 0, 0, integrate(function(tau) exp(scaled_lbetaprime(tau, tau0, alpha)), 0, lo)$value)
  Z <- hiZ - loZ

  # Integrate likelihood with respect to prior
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

  # Transform \delta to \rho
  # Change normalizing constant when prior is restricted
  to_rho <- function(delta) ifelse(delta == Inf, 1, delta^2 / (1 + delta^2))
  lo <- to_rho(interval[1])
  hi <- to_rho(interval[2])
  Z <- pbeta(hi, alpha, alpha) - pbeta(lo, alpha, alpha)

  # Rmpfr provides arbitrary precision floating point arithmetic
  # Integrate likelihood with respect to prior
  value <- Rmpfr::integrateR(function(rho) {
    rho <- Rmpfr::mpfr(rho, 100)
    llh <- ((n1 - 1)/2 + alpha - 1) * log(rho) +
           ((n2 - 1)/2 + alpha - 1) * log(1 - rho) +
           ((2 - n)/2) * log(rho*n1*s1 + (1 - rho)*n2*s2)

    exp(llh - log(Z) - lbeta(alpha, alpha))
  }, lo, hi)$value

  log(value)
}


# Checks whether the intervals for the null and the alternative are correctly specified
.check_interval_input <- function(alternative_interval, null_interval) {

  alt_condition <- all(is.numeric(alternative_interval)) &&
                   alternative_interval[2] > alternative_interval[1] &&
                   all(alternative_interval >= 0)

  if (!alt_condition) {
    stop('Something is off with your specification of the alternative_interval!')
  }

  if (!is.null(null_interval)) {
    null_condition <- all(is.numeric(null_interval)) &&
                      all(null_interval >= 0) &&
                      null_interval[2] > null_interval[1] &&
                      all(null_interval >= 0)

    if (!null_condition) {
      stop('Something is off with your specification of the null_interval!')
    }
  }
}


# Computes the prior probability with which \rho follows the hypothesis constraint
# 1<2<3 -> 1 / factorial(3)
# 1<2,3 -> factorial(2) / factorial(3)
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


# Remove equality constraints from the hypothesis, e.g. '1,2<3=4,5<6' -> '1,2<3,4<5'
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

  # Add order-constraints
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


# Creates a list that is passed to the Stan model
.prepare_standat <- function(hyp, ns, ss, a, priors_only = FALSE) {
  k <- length(ns)

  if (k != length(ss)) {
    stop('Need the same amount of information for each group!')
  }

  # Split hypothesis and count number of equalities (=),
  # order constraints (<), and no constraints (,)
  check <- gsub('[0-9<=, ]', '.', hyp)
  rel <- strsplit(gsub('[0-9]', '', hyp), '')[[1]]
  rel[which(rel == '=')] <- 1
  rel[which(rel == '<')] <- 2
  rel[which(rel == ',')] <- 3
  rel <- as.numeric(rel)

  nr_equal <- sum(rel == 1)
  nr_ordered <- sum(rel == 2)

  # Get a vector whose elements give the parameter of the group
  # This is used for equality constraints; for example, 1=2<3 results in c(1, 1, 2)
  index <- 1
  count <- 1
  index_vector <- numeric(k)
  strip <- strsplit(hyp, '')[[1]]

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

  standat
}


# Extend print function to only show parameters tau and rho
print.bfvar <- function(x) {
  pars <- c('tau', 'rho')
  print(x$fit, pars = pars)

  if (!is.null(x$logml)) {
    cat(paste0('\nLog Marginal Likelihood = ', round(x$logml, 3), '\n'))
  }
}


#' Estimates the model using Stan and computes the marginal likelihood
#'
#' @param hypothesis a string specifying the hypothesis
#' @param ss a vector containing sample sum of squares
#' @param ns a vector containing sample sizes
#' @param a a vector specifying the value of the parameters of the Dirichlet prior
#' @param compute_ml a logical specifying whether the marginal likelihood should be computes
#' @param priors_only a logical specifying whether we should only sample from the prior
#' @param precision a logical specifying whether parameterization in terms of precision or variance is used
#' @param ... arguments to rstan::sampling
#'
#' @returns an object of class 'bfvar', which is a stanfit object with a log marginal likelihood (if desired)
.create_bfvar_object <- function(hyp, ns, ss, a, compute_ml = TRUE, priors_only = FALSE, ...) {

  # Translate hypothesis on variance / standard deviation into hypothesis on precision
  hyp <- gsub('>', '<', hyp)
  logml <- NULL
  standat <- .prepare_standat(hyp, ns, ss, a, priors_only = priors_only)

  # If hypothesis contrains only equalities (e.g., 1=2=3) or is of
  # mixed type (e.g., 1,2<3), then we use the Mixed.stan model
  if (!.is_only_ordered_and_equal(hyp) || .is_allequal(hyp)) {
    fit <- rstan::sampling(stanmodels$Mixed, data = standat, ...)

    if (compute_ml && !priors_only) {
      logml <- bridgesampling::bridge_sampler(fit)$logml
    }

  # If hypothesis contains only equalities and order-constraints (e.g., 1=2>3),
  # then we use the Ordered.stan model
  } else if (.is_only_ordered_and_equal(hyp)) {
    fit <- rstan::sampling(stanmodels$Ordered, data = standat, ...)

    if (compute_ml && !priors_only) {
      logml <- bridgesampling::bridge_sampler(fit)$logml
    }

  # For mixed hypothesis with no order-constraints (e.g., 1=2,3) we also use Mixed.stan
  # However, in contrast to above, we cannot rely on bridgesampling to
  # compute marginal likelihoods but have to use the Kluglist & Hoijtink (2005) trick
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
