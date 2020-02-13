library('rstan')
library('bfvartest')


context('K = 1 Sample Test')
test_that('One-sample Test Gives Expected Result', {
  expect_lt(onesd_test(100, 1, 1, 0.5), 0) # log BF10 < 0
  expect_gt(onesd_test(100, 1, 2, 0.5), 0) # log BF10 > 0
})


context('Interval Specification')
test_that('Function stops if wrong interval input', {
  expect_error(.check_interval_input(alternative_interval = NULL, null_interval = NULL))
  expect_error(.check_interval_input(alternative_interval = c(-2, 4), null_interval = NULL))
  expect_error(.check_interval_input(alternative_interval = c(0, Inf), null_interval = c(-2, 4)))
  expect_error(.check_interval_input(alternative_interval = c(0, Inf), null_interval = c(.5, .2)))
})


test_that('Two-sample Test Gives Expected Result', {
  expect_true(twosd_test(100, 100, 1, 1, 0.5) < 0)
  expect_true(twosd_test(100, 100, 1, 2, 0.5) > 0)

  # 2 is larger than 1
  expect_true(twosd_test(100, 100, 1, 2, 0.5, alternative_interval = c(1, Inf)) > 0)

  # 1 is larger than 2
  expect_true(twosd_test(100, 100, 1, 2, 0.5, alternative_interval = c(0, 1)) < 0)
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
  hyp <- c('1>2,3,4>5,6>7')
  hyp_prec <- gsub('>', '<', hyp)

  ns <- rep(100, 7)
  sds <- c(1, 2, 2, 5, 6, 6, 7)
  bfvar <- .create_bfvar_object(hyp, ns, sds, 0.50, priors_only = TRUE, iter = 10000, cores = 6, chains = 6)
  rho <- rstan::extract(bfvar$fit, 'rho')$rho

  hyp_fn <- eval(parse(text = .create_hyp_fn(hyp_prec)))
  expect_true(abs(mean(apply(rho, 1, hyp_fn)) - .compute_prior_restr(hyp)) < 0.005)
})


test_that('User input check is valid', {
  sds <- c(1, 2, 3)
  ns <- c(100, 100, 100)

  expect_error(.check_user_input('1>2>3', ns, sds))
  expect_error(.check_user_input(c('1>2>3', '1=2=3'), ns[-1], sds))
  expect_error(.check_user_input(c('1<2<3', '1=2=3'), ns, sds))
  expect_true(is.null(.check_user_input(c('1>2>3', '1=2=3'), ns, sds)))
})



context('K > 2 Sample Test')
test_that('ksd_test gives same result as (undirected) twosd_test for K = 2', {
  sds <- c(1, 1)
  ns <- c(100, 100)
  hyp <- c('1,2', '1=2')

  res <- ksd_test(hyp, ns, sds, alpha = 0.50)
  bf10 <- twosd_test(ns[1], ns[2], sds[1], sds[2])
  expect_true(abs(res$BF[1, 2] - bf10) < .05)
  # -288.4641 '1,2'
  # -284.7882 '1=2'
})


test_that('ksd_test gives same result as (directed) twosd_test for K = 2', {
  sds <- c(2, 1)
  ns <- c(100, 100)
  hyp <- c('1,2', '1>2')

  bf10 <- twosd_test(ns[1], ns[2], sds[1], sds[2])
  bfr0 <- twosd_test(ns[1], ns[2], sds[1], sds[2], alternative_interval = c(0, 1))
  bfr1 <- bfr0 - bf10

  res <- ksd_test(hyp, ns, sds, alpha = 0.50, priors_only = FALSE)
  expect_true(abs(res$BF[2, 1] - bfr1) < .05)
})


test_that('All cases work', {
  sds <- c(3, 2, 1)
  ns <- c(100, 100, 100)
  hyp <- c('1>2>3', '1=2=3', '1,2,3', '1,2=3', '1,2>3', '1=2>3')

  res <- ksd_test(hyp, ns, sds, alpha = 0.50)
})


# Z <- function(k, nr_equal, a = 0.50) {
#   lgamma(a * (k - nr_equal)) - sum(lgamma(rep(a, k - nr_equal)))
# }

# ns <- rep(50, 3)
# sds <- c(1, 1, 1)
# hyp <- c('1,2,3', '1=2=3', '1>2>3')
# hyp_mul <- c('1,2,3', '1=2=3', '1<2<3')
# res <- ksd_test(hyp, ns, sds)
#
# lml <- log_marginal_likelihoods(sds, ns, hypotheses = hyp_mul)
# mulderBF <- bayes_factors(lml, log.BF = TRUE)$AFBF
# colnames(mulderBF) <- rownames(mulderBF) <- colnames(res$BF)
# mulderBF
# res$BF
#
# ns <- rep(50, 2)
# sds <- c(1, 1)
# hyp <- c('1,2', '1=2', '1>2')
# hyp_mul <- c('1,2', '1=2', '1<2')
# res <- ksd_test(hyp, ns, sds)
#
# lml <- log_marginal_likelihoods(sds, ns, hypotheses = hyp_mul)
# mulderBF <- bayes_factors(lml, log.BF = TRUE)$AFBF
# colnames(mulderBF) <- rownames(mulderBF) <- colnames(res$BF)
# mulderBF
# res$BF
#
# bf10 <- twosd_test(ns[1], ns[2], sds[1], sds[2])
