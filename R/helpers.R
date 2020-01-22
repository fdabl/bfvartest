# Computes the log marginal likelihood for K = 2 hypothesis
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


# Takes a hypothesis and creates a matrix encoding the order-constraints of the elements (used in Stan)
.create_constraint_matrix <- function(hyp) {
  check <- gsub('[0-9>=, ]', '.', hyp)
  correct <- sapply(strsplit(check, '')[[1]], function(el) el == '.')

  # Check whether input is correct
  if (!all(correct)) { stop('Hypothesis can contain only [digits > = ,]') }
  # if (length(correct) < 5) { stop('Hypothesis needs at least three group elements.')}

  # Check if hypothesis is ordered from low to high
  check <- gsub('[>=,]', '.', hyp)
  num <- as.numeric(unlist(strsplit(check, '\\.')))

  if (!all(num == seq(length(num)))) {
    stop('Hypothesis needs to be given in ascending numeric order (e.g., 1>2>3 instead of 3>2>1)')
  }

  # Split the hypothesis into ordered blocks, initialize constraints between these blocks
  ordered_blocks <- strsplit(hyp, '>')[[1]]
  constraint_mat <- matrix(NA, nrow = length(ordered_blocks), ncol = 2)

  row <- 1
  for (i in seq(1, length(ordered_blocks))) {

    if (i == 1) {

      # First element does not have a lower constraint
      prev <- -99

      # Upper constraint is the first element of the next ordered block
      nnext <- ordered_blocks[i + 1]

    } else if (i == length(ordered_blocks)) {

      # Last element does not have an upper constraint
      nnext <- -99

      # Lower constraint is the last element of the previous ordered block
      prev <- ordered_blocks[i - 1]

    } else {

      # In between blocks have both lower and upper constraints

      # Lower constraint is the last element of the previous ordered block
      prev <- ordered_blocks[i - 1]

      # Upper constraint is the last element of the previous ordered block
      nnext <- ordered_blocks[i + 1]
    }

    # since our constraints are of the form '1>2'
    constraint_mat[row, ] <- as.numeric(cbind(nnext, prev))
    row <- row + 1
  }

  # If the ordered block is of length one, we do not have any order-constraints!
  if (length(ordered_blocks) == 1) {
    constraint_mat <- matrix(c(-99, -99), 1, 2)
  }

  constraint_mat
}


# Creates a list that is passed onto Stan
.prepare_standat <- function(hyp, ns, ss, a, priors_only = FALSE) {
  k <- length(ns)

  if (k != length(ss)) {
    stop('Need same amount of information for each group.')
  }

  constraint_mat <- .create_constraint_matrix(hyp)

  check <- gsub('[0-9>=, ]', '.', hyp)
  rel <- strsplit(gsub('[0-9]', '', hyp), '')[[1]]
  rel[which(rel == '=')] <- 1
  rel[which(rel == '>')] <- 2
  rel[which(rel == ',')] <- 3
  rel <- as.numeric(rel)

  nr_equal <- sum(rel == 1)
  nr_ordered <- sum(rel == 2)
  nr_free <- sum(rel == 3)

  strip_num <- strsplit(gsub('[>,=]', ' ', hyp), ' ')[[1]]
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

    if (strip[i] %in% c('>', ',')) {
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
    constraint_mat = constraint_mat,
    relations = rel, nr_rel = length(rel), priors_only = priors_only,
    nr_equal = nr_equal, nr_ordered = nr_ordered, nr_free = nr_free
  )

  # hack so that it works with only two groups (for testing)
  if (standat$nr_rel == 1) {
    standat$nr_rel <- 2
    if (standat$relations[1] == 3) {
      standat$relations <- c(standat$relations, 3)
    } else {
      standat$relations <- c(standat$relations, -99)
    }
  }

  standat
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
#' @returns the log marginal likelihood (and samples if specified) of the hypothesis
.create_levene_object <- function(hyp, ns, ss, a, compute_ml = TRUE, priors_only = FALSE, ...) {


  if (hyp == 'ord') {
    standat <- .prepare_standat('1>2', ns, rev(ss), a, priors_only = priors_only)
    stanres <- rstan::sampling(stanmodels$BayesLeveneOrdered, data = standat, ...)
  } else {
    standat <- .prepare_standat(hyp, ns, ss, a, priors_only = priors_only)
    stanres <- rstan::sampling(stanmodels$BayesLevene, data = standat, ...)
  }

  if (!compute_ml) {
    return(list(
      'fit' = stanres
    ))
  }

  logml <- bridgesampling::bridge_sampler(stanres)$logml
  res <- list('fit' = stanres, 'logml' = logml)
  class(res) <- 'Levene'
  res
}
