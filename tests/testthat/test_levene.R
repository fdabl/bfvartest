library('rstan')
library('BayesLevene')


context('K = 1 Sample Test')
test_that('One-sample Test Gives Expected Result', {
  expect_lt(one_sample(100, 1, 1, 0.5), 0) # log BF10 < 0
  expect_gt(one_sample(100, 1, 2, 0.5), 0) # log BF10 > 0
})


context('Interval Specification')
test_that('Function stops if wrong interval input', {
  expect_error(.check_interval_input(alternative_interval = NULL, null_interval = NULL))
  expect_error(.check_interval_input(alternative_interval = c(-2, 4), null_interval = NULL))

  expect_error(.check_interval_input(alternative_interval = c(0, Inf), null_interval = c(-2, 4)))
  expect_error(.check_interval_input(alternative_interval = c(0, Inf), null_interval = c(.5, .2)))
})


test_that('Two-sample Test Gives Expected Result', {
  expect_true(two_sample(100, 100, 1, 1, 0.5) < 0)
  expect_true(two_sample(100, 100, 1, 2, 0.5) > 0)

  # 2 is larger than 1
  expect_true(two_sample(100, 100, 1, 2, 0.5, alternative_interval = c(1, Inf)) > 0)

  # 1 is larger than 2
  expect_true(two_sample(100, 100, 1, 2, 0.5, alternative_interval = c(0, 1)) < 0)
})


context('K > 2 Sample Helper Functions')
test_that('Constrained matrix errors expectedly', {
  expect_error(.create_constraint_matrix('2>1')) # needs to be ascending
  expect_error(.create_constraint_matrix('1<2')) # < not allowed
  expect_error(.create_constraint_matrix('1u2,3>4')) # u not allowed
})

test_that('Constrained matrix handles mixed hypotheses', {
  hyp <- '1,2>3=4>5>6=7=8>9'
  cmat <- .create_constraint_matrix(hyp)
  emat <- cbind(
    c(3, 5, 6, 9, -99), # lower constraints
    c(-99, 2, 4, 5, 8) # upper constraints
  )
  expect_true(all(cmat == emat))
})

test_that('Group sizes K > 9 constraints work', {
  hyp <- c('1>2>3>4>5>6>7>8>9>10>11>12')
  lower <- c(seq(2, 12), -99)
  upper <- c(-99, seq(11))
  expect_true(all(.create_constraint_matrix(hyp) == cbind(lower, upper)))
})


context('K > 2 Sample Test')
test_that('Stan order_parameters function works', {
  rstan::expose_stan_functions(stanmodels$BayesLevene)

  get_oparams <- function(hyp, k) {
    rho <- runif(k)
    beta <- runif(k, 0, 1)

    standat <- .prepare_standat(hyp, rep(100, k), rep(1, k), .50)
    cmat <- lapply(1:nrow(standat$constraint_mat), function(i) standat$constraint_mat[i, ])
    order_parameters(
      rho, beta, length(rho), nr_rel = standat$nr_rel,
      relations = standat$relations, constraint_mat = cmat
    )
  }

  r <- get_oparams('1>2>3>4', 4)
  expect_true(r[1] > r[2] && r[2] > r[3] && r[3] > r[4])

  r <- get_oparams('1=2=3=4', 4)
  expect_true(all(r == 1/4))

  r <- get_oparams('1,2,3,4', 4)
  expect_true(all(diff(r) != 0) && sum(r) == 1)

  r <- get_oparams('1,2>3=4>5>6=7=8>9', 9)
  expect_true(
    r[1] != r[2] && r[1] > r[3] && r[2] > r[3] && r[3] == r[4] &&
    r[4] > r[5] && r[5] > r[6] && r[6] == r[7] && r[7] == r[8] && r[8] > r[9]
  )
})

# rho <- t(sapply(seq(1000), function(i) get_oparams('1>2', 2)))
# rho <- t(replicate(1000, sim(.50)))


# sim <- function(a = .50) {
#   b <- runif(1)
#   rho1 <- rbeta(1, a, a)
#   rho2 <- 1 - rho1
#
#   rho <- c(b * (1 - rho2) + rho2)
#   rho <- c(rho, (1 - b) * rho[1])
#   rho / sum(rho)
# }


test_that('k_sample errors correctly', {
  ss <- c(1, 2, 3)
  ns <- c(100, 100, 100)

  expect_error(k_sample('1>2>3', ns, ss, alpha = 0.50))
  expect_error(k_sample(c('1>2>3', '1=2=3'), ns[-1], ss, alpha = 0.50))
})


test_that('k_sample gives same result as (undirected) two_sample for K = 2', {
  ss <- c(1, 2)
  ns <- c(100, 100)
  hyp <- c('1,2', '1=2')

  res <- k_sample(hyp, ns, ss, alpha = 0.50)
  bf10 <- two_sample(ns[1], ns[2], ss[1], ss[2])
  expect_true(abs(res$BF[1, 2] - bf10) < .1)
})


test_that('k_sample gives same result as (directed) two_sample for K = 2', {
  ss <- c(1, 2)
  ns <- c(100, 100)
  hyp <- c('1,2', '1>2', 'ord')

  bf10 <- two_sample(ns[1], ns[2], ss[1], ss[2])
  bfr0 <- two_sample(ns[1], ns[2], ss[1], ss[2], alternative_interval = c(1, Inf))
  bfr1 <- bfr0 - bf10

  res <- k_sample(hyp, ns, ss, alpha = 0.50, priors_only = FALSE, iter = 10000)
  res <- k_sample(hyp, ns, ss, alpha = 0.50, priors_only = TRUE)
  expect_true(abs(res$BF[2, 1] - bfr1) < .1)
})

# fit1 <- res[[1]]$fit
# rho <- extract(fit1, 'rho')$rho
# rhou <- extract(fit1, 'rho_unconstrained')$rho_unconstrained
#
# fit2 <- res[[2]]$fit
# rho2 <- extract(fit2, 'rho')$
# rhou2 <- extract(fit2, 'rho_unconstrained')$rho_unconstrained
#
# fit3 <- res[[3]]$fit
# rho3 <- extract(fit3, 'rho')$rho[, 2:1]
#
# standat1 <- .prepare_standat(hyp[1], ns, ss, .5)
# standat2 <- .prepare_standat(hyp[2], ns, ss, .5)
# get_oparams(hyp[2], 2)
# stanres1 <- rstan::sampling(stanmodels$BayesLevene, data = standat1)
# stanres2 <- rstan::sampling(stanmodels$BayesLevene, data = standat2)
# res <- k_sample(hyp, ns, ss, alpha = 0.50, priors_only = FALSE)

test_that('All cases work for K = 3', {
  ss <- c(1, 2, 3)
  ns <- c(100, 100, 100)
  hyp <- c('1>2>3', '1=2=3')

  k_sample(hyp, ns, ss, alpha = 0.50)
})


test_that('k_sample works for K > 9 groups', {
  hyp <- c('1>2>3>4>5>6>7>8>9>10>11>12', '1,2,3,4,5,6,7,8,9,10,11,12')
  ss <- seq(12)
  ns <- rep(100, 12)
  res <- k_sample(hyp, ns, ss)
})
