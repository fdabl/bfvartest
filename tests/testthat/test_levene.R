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
test_that('Parses the input correctly into a function', {
  hyp <- '1<2'
  hyp_fn <- .create_hyp_fn(hyp)
  hyp_want <- 'function(x) max(x[c(1)]) < min(x[c(2)])'
  expect_true(identical(hyp_fn, hyp_want))

  hyp <- '1<2,3'
  hyp_fn <- .create_hyp_fn(hyp)
  hyp_want <- 'function(x) max(x[c(1)]) < min(x[c(2,3)])'
  expect_true(identical(hyp_fn, hyp_want))

  hyp <- '1<2,3=4'
  hyp_fn <- .create_hyp_fn(hyp)
  hyp_want <- 'function(x) max(x[c(1)]) < min(x[c(2,3)])'
  expect_true(identical(hyp_fn, hyp_want))

  hyp <- '1,2<3,4=5<6,7=8<9'
  hyp_fn <- .create_hyp_fn(hyp)
  hyp_want <- 'function(x) max(x[c(1,2)]) < min(x[c(3,4)]) && max(x[c(3,4)]) < min(x[c(5,6)]) && max(x[c(5,6)]) < min(x[c(7)])'
  expect_true(identical(hyp_fn, hyp_want))
})


test_that('Only Ordered and equal check works', {
  expect_true(.is_only_ordered_and_equal('1<2<3'))
  expect_true(.is_only_ordered_and_equal('1<2=3'))
  expect_false(.is_only_ordered_and_equal('1<2,3'))
  expect_false(.is_only_ordered_and_equal('1<2=3,4<5'))
})


test_that('Prior probability of restriction is correct', {
  hyp <- c('1<2,3,4<5,6<7')
  ns <- rep(100, 7)
  ss <- c(1, 2, 2, 5, 6, 6, 7)
  standat <- .prepare_standat(hyp, ns, ss, 0.50, priors_only = TRUE)
  fit <- rstan::sampling(stanmodels$BayesLeveneMixed, data = standat, iter = 100000, cores = 6, chains = 6)
  rho <- extract(fit, 'rho')$rho

  hyp_fn <- eval(parse(text = .create_hyp_fn(hyp)))
  expect_true(abs(mean(apply(rho, 1, hyp_fn)) - .compute_prior_restr(hyp)) < 0.001)
})


test_that('User input check is valid', {
  ss <- c(1, 2, 3)
  ns <- c(100, 100, 100)

  expect_error(.check_user_input('1<2<3', ns, ss))
  expect_error(.check_user_input(c('1<2<3', '1=2=3'), ns[-1], ss))
  expect_error(.check_user_input(c('1>2>3', '1=2=3'), ns, ss))
})



context('K > 2 Sample Test')
test_that('k_sample gives same result as (undirected) two_sample for K = 2', {
  ss <- c(1, 1)
  ns <- c(100, 100)
  hyp <- c('1,2', '1=2')

  res <- k_sample(hyp, ns, ss, alpha = 0.50)
  bf10 <- two_sample(ns[1], ns[2], ss[1], ss[2])
  expect_true(abs(res$BF[1, 2] - bf10) < .1)
})


test_that('k_sample gives same result as (directed) two_sample for K = 2', {
  ss <- c(1, 2)
  ns <- c(100, 100)
  hyp <- c('1,2', '1<2')

  bf10 <- two_sample(ns[1], ns[2], ss[1], ss[2])
  bfr0 <- two_sample(ns[1], ns[2], ss[1], ss[2], alternative_interval = c(1, Inf))
  bfr1 <- bfr0 - bf10

  res <- k_sample(hyp, ns, ss, alpha = 0.50, priors_only = FALSE)
  expect_true(abs(res$BF[2, 1] - bfr1) < .1)
})


test_that('All cases work for K = 3', {
  ss <- c(1, 2, 3)
  ns <- c(100, 100, 100)
  hyp <- c('1<2<3', '1=2=3')

  k_sample(hyp, ns, ss, alpha = 0.50)
})


test_that('k_sample works for K > 9 groups', {
  hyp <- c('1<2<3<4<5<6<7<8<9<10<11<12', '1,2,3,4,5,6,7,8,9,10,11,12')
  ss <- seq(12)
  ns <- rep(100, 12)
  res <- k_sample(hyp, ns, ss)
})



# logml <- bridgesampling::bridge_sampler(fit)$logml
#
# ldir <- sum(lgamma(rep(a, k))) - lgamma(a*k)
# ldir2 <- sum(lgamma(rep(a, k-1))) - lgamma(a*(k-1))
# l <- log(mean(apply(post, 1, function(rho) rho[1] > rho[2] & rho[2] > rho[4])) / (1/6)) + logml
#
# log_marginal_likelihood(ss, ns, hypothesis = c('1<2=3<4,5'))
# log_marginal_likelihood(ss, ns, hypothesis = c('1,2=3,4,5'))
# lml <- log_marginal_likelihoods(ss, ns, hypotheses = c('1,2=3,4,5', '1<2=3<4,5', '1,2,3,4,5'))
# bayes_factors(lml, log.BF = TRUE)
