# Computes the log marginal likelihood for H1 in the one-sample test
.compute_logml_restr_k1 <- function(n, s2, popsd, interval, alpha = 0.50) {
  popvar <- popsd^2
  tau0 <- 1 / popvar

  # Prior distribution centered at \tau_0
  scaled_lbetaprime <- function(tau, tau0, alpha) {
    value <- (alpha - 1) * log(tau / tau0) - 2 * alpha *
             log(1 + tau / tau0) - log(tau0) - lbeta(alpha, alpha)
    value
  }

  # Change normalizing constant when prior is restricted
  lo <- interval[1]
  hi <- interval[2]
  hi_Z <- stats::integrate(function(tau) exp(scaled_lbetaprime(tau, tau0, alpha)), 0, hi)$value
  lo_Z <- ifelse(lo == 0, 0, stats::integrate(function(tau) exp(scaled_lbetaprime(tau, tau0, alpha)), 0, lo)$value)
  Z <- hi_Z - lo_Z

  # Integrate likelihood with respect to prior
  value <- stats::integrate(function(tau) {
    llh <- (n - 1) / 2 * log(tau) - 0.5 * tau * n * s2
    lprior <- scaled_lbetaprime(tau, tau0, alpha)

    exp(llh + lprior - log(Z))
  }, lo, hi)$value

  log(value)
}


# Computes the log marginal likelihood for H1 in the two-sample test
.compute_logml_restr_k2 <- function(n1, n2, s1, s2, interval, alpha = 0.50) {
  n <- n1 + n2

  # Transform \delta to \rho
  # Change normalizing constant when prior is restricted
  to_rho <- function(delta) ifelse(delta == Inf, 1, delta^2 / (1 + delta^2))
  lo <- to_rho(interval[1])
  hi <- to_rho(interval[2])
  Z <- stats::pbeta(hi, alpha, alpha) - stats::pbeta(lo, alpha, alpha)

  # Rmpfr provides arbitrary precision floating point arithmetic
  # Integrate likelihood with respect to prior
  value <- Rmpfr::integrateR(function(rho) {
    rho <- Rmpfr::mpfr(rho, 100)
    llh <- ((n1 - 1) / 2 + alpha - 1) * log(rho) +
           ((n2 - 1) / 2 + alpha - 1) * log(1 - rho) +
           ((2 - n) / 2) * log(rho * n1 * s1 + (1 - rho) * n2 * s2)

    exp(llh - log(Z) - lbeta(alpha, alpha))
  }, lo, hi)$value


  log(value)
}


# Required for the posterior of \delta in the K = 2 case
.Gauss2F1 <- function(a, b, c, x){
  if(x >= 0 & x < 1){
    gsl::hyperg_2F1(a, b, c, x)
  } else {
    gsl::hyperg_2F1(c - a, b, c, 1 - 1 / (1 - x))/ (1 - x)^b
  }
}


# Checks whether the intervals for the null and the alternative are correctly specified
.check_interval_input <- function(alternative_interval, null_interval) {

  is_correct <- function(interval) {
    all(is.numeric(interval)) &&
    interval[2] > interval[1] &&
    all(interval >= 0)
  }

  if (!is_correct(alternative_interval)) {
    stop('Something is off with your specification of the alternative_interval!')
  }

  if (!is.null(null_interval)) {
    null_condition <- all(is.numeric(null_interval)) &&
                      null_interval[2] > null_interval[1] &&
                      all(null_interval >= 0)

    if (!is_correct(null_interval)) {
      stop('Something is off with your specification of the null_interval!')
    }
  }
}


# Computes the prior probability with which \rho follows the hypothesis constraint
# 1<2<3 -> 1 / factorial(3)
# 1<2,3 -> factorial(2) / factorial(3)
.compute_prior_restr <- function(hyp) {
  sblock <- strsplit(hyp, '<')[[1]]
  nblock <- length(sblock)

  ss <- strsplit(hyp, '')[[1]]
  ngroups <- as.numeric(ss[length(ss)])
  num <- factorial(ngroups)

  for (i in seq(nblock)) {
    num <- num / factorial(length(strsplit(sblock[i], ',')[[1]]))
  }

  1 / num
}


# Remove equality constraints from the hypotheses
# e.g. c('1=2', '3,4,5=6', '7,8') -> c('1', '2,3,4', '5,6')
.reduce_hyp_equals <- function(hyp) {
  n <- length(hyp)
  res <- c()
  j <- 1

  for (i in seq(n)) {
    si <- strsplit(hyp[i], ',')[[1]]
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


# Checks whether hypothesis only includes ordinal and no constraints (e.g., '1>2=3>4')
.is_only_ordered_and_equal <- function(hyp) {
  s <- strsplit(hyp, '>')[[1]]

  check <- sapply(s, function(si) {
    length(strsplit(si, ',')[[1]]) == 1
  })

  all(check) && !.is_allequal(hyp)
}


# Checks whether hypothesis only includes equality and no constraints (e.g., '1,2=3,4')
.is_only_unconstr_and_equal <- function(hyp) {
  s <- strsplit(hyp, ',')[[1]]

  check <- sapply(s, function(si) {
    length(strsplit(si, '<')[[1]]) == 1
  })

  all(check) && !.is_allequal(hyp) && !.is_allunequal(hyp)
}


# Checks whether hypothesis is all equal (e.g., '1=2=3=4')
.is_allequal <- function(hyp) {
  s <- strsplit(hyp, '=')[[1]]
  all(nchar(s) < 2)
}


# Checks whether hypothesis is all equal (e.g., '1,2,3,4')
.is_allunequal <- function(hyp) {
  s <- strsplit(hyp, ',')[[1]]
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
#' @export
print.bfvar <- function(x, ...) {
  print(x$fit, pars = 'sds')

  if (!is.null(x$logml)) {
    cat(paste0('\nLog Marginal Likelihood = ', round(x$logml, 3), '\n'))
  }
}


#' Estimates the model using Stan and computes the marginal likelihood
#'
#' @param hyp a string specifying the hypothesis
#' @param ss a vector containing sample sum of squares
#' @param ns a vector containing sample sizes
#' @param a a vector specifying the value of the parameters of the Dirichlet prior
#' @param compute_ml a logical specifying whether the marginal likelihood should be computed
#' @param priors_only a logical specifying whether we should only sample from the prior
#' @param silent a logical specifying whether to print results from sampling and bridgesampling
#' @param ... arguments to rstan::sampling
#' @returns an object of class 'bfvar', which is a stanfit object with a log marginal likelihood (if desired)
.create_bfvar_object <- function(hyp, ns, ss, a, compute_ml = TRUE, priors_only = FALSE, silent = TRUE, ...) {

  # Translate hypothesis on variance / standard deviation into hypothesis on precision
  hyp <- gsub('>', '<', hyp)
  logml <- NULL
  standat <- .prepare_standat(hyp, ns, ss, a, priors_only = priors_only)
  refresh <- ifelse(silent, 0, 200)
  standat$nr_equal <- standat$nr_equal * !.is_allequal(hyp)

  # If hypothesis contrains only equalities (e.g., 1=2=3) or
  # contains no constraints (e.g., 1,2,3) or a mix (e.g., 1=2,3), then we use the Mixed.stan model
  if (.is_allunequal(hyp) || .is_allequal(hyp) || .is_only_unconstr_and_equal(hyp)) {
    fit <- suppressWarnings(rstan::sampling(stanmodels$Mixed, data = standat, refresh = refresh, ...))

    if (compute_ml && !priors_only) {
      logml <- bridgesampling::bridge_sampler(fit, silent = silent)$logml
    }

  # If hypothesis contains only equalities and order-constraints (e.g., 1=2>3),
  # then we use the Ordered.stan model
  } else if (.is_only_ordered_and_equal(hyp)) {
    fit <- suppressWarnings(rstan::sampling(stanmodels$Ordered, data = standat, refresh = refresh, ...))

    if (compute_ml && !priors_only) {
      logml <- bridgesampling::bridge_sampler(fit, silent = silent)$logml
    }

  # For mixed hypothesis with order-constraints (e.g., 1=2,3>4) we also use Mixed.stan
  # However, in contrast to above, we cannot rely on bridgesampling to
  # compute marginal likelihoods but have to use the Kluglist & Hoijtink (2005) trick
  } else {

    # Example:
    # (1) Remove order-constraints from the hypothesis '1=2,3>4' (Hr) to yield '1=2,3,4' (H1)
    # (2) Compute the Bayes factor BFr1 in favour of '1=2,3>4' against '1=2,3=4' using Kluglist & Hoijtink (2005)
    # (3) Get the marginal likelihood of H1 using bridgesampling
    # (4) Get the marginal likelihood of Hr by multiplying BFr1 with the marginal likelihood of H1

    hyp1 <- gsub('<', ',', hyp) # (1)
    standat <- .prepare_standat(hyp1, ns, ss, a, priors_only = priors_only)
    fit <- suppressWarnings(rstan::sampling(stanmodels$Mixed, data = standat, refresh = refresh, ...))

    if (compute_ml && !priors_only) {
      rho <- rstan::extract(fit, 'rho')$rho
      rho <- t(apply(rho, 1, unique)) # only get not equal parameters
      hyp_fn <- eval(parse(text = .create_hyp_fn(hyp)))

      # Kluglist & Hoijtink (2005) trick
      BFr1 <- log(mean(apply(rho, 1, hyp_fn))) - log(.compute_prior_restr(hyp)) # (2)
      logml1 <- bridgesampling::bridge_sampler(fit, silent = silent)$logml # (3)
      logml <- BFr1 + logml1 # This is the log marginal likelihood of Hr # (4)
    }
  }

  res <- list('fit' = fit, 'logml' = logml)
  class(res) <- 'bfvar'
  res
}
