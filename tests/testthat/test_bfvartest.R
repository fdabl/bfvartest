context('K = 1 Sample Test')
test_that('K = 1 Test gives expected results', {
  expect_lt(onesd_test(100, 1, 1, 0.5), 0) # log BF10 < 0
  expect_gt(onesd_test(100, 1, 2, 0.5), 0) # log BF10 > 0

  # H0 true
  BF_10 <- onesd_test(100, 1, 1, 0.5)
  BF_overlapping_10 <- onesd_test(100, 1, 1, 0.5, null_interval = c(0.90, 1.10))
  BF_nonoverlapping_10 <- onesd_test(
    100, 1, 1, 0.5, null_interval = c(0.90, 1.10), nonoverlapping_interval = TRUE
  )

  # More evidence for sharp H0
  expect_gt(BF_overlapping_10, BF_10)

  # Non-overlapping H1 should do worse
  expect_gt(BF_overlapping_10, BF_nonoverlapping_10)

  # H1 true
  BF_10 <- onesd_test(100, 1.2, 1, 0.5)
  BF_overlapping_10 <- onesd_test(100, 1.2, 1, 0.5, null_interval = c(0.90, 1.10))
  BF_nonoverlapping_10 <- onesd_test(
    100, 1.2, 1, 0.5, null_interval = c(0.90, 1.10),
    nonoverlapping_interval = TRUE
  )

  # More evidence against sharp H0
  expect_gt(BF_10, BF_overlapping_10)

  # Non-overlapping H1 should do better
  expect_gt(BF_nonoverlapping_10, BF_overlapping_10)
})


context('K = 2 Sample Test')
test_that('K = 2 test errors correctly (input check)', {
  expect_error(.check_interval_input(alternative_interval = NULL, null_interval = NULL))
  expect_error(.check_interval_input(alternative_interval = c(-2, 4), null_interval = NULL))
  expect_error(.check_interval_input(alternative_interval = c(0, Inf), null_interval = c(-2, 4)))
  expect_error(.check_interval_input(alternative_interval = c(0, Inf), null_interval = c(.5, .2)))
})


test_that('K = 2 test gives expected results', {
  expect_lt(twosd_test(100, 100, 1, 1, 0.5), 0)
  expect_gt(twosd_test(100, 100, 1, 2, 0.5), 0)

  # H1 is true: 2 is larger than 1
  expect_gt(twosd_test(100, 100, 1, 2, 0.5, alternative_interval = c(1, Inf)), 0)

  # H1 is true: 1 is larger than 2
  expect_lt(twosd_test(100, 100, 1, 2, 0.5, alternative_interval = c(0, 1)), 0)

  # H0 true
  BF_10 <- twosd_test(100, 100, 1, 1, 0.5)
  BF_overlapping_10 <- twosd_test(100, 100, 1, 1, 0.5, null_interval = c(0.90, 1.10))
  BF_nonoverlapping_10 <- twosd_test(
    100, 100, 1, 1, 0.5, null_interval = c(0.90, 1.10), nonoverlapping_interval = TRUE
  )

  # More evidence for sharp H0
  expect_gt(BF_overlapping_10, BF_10)

  # Non-overlapping H1 should do worse
  expect_gt(BF_overlapping_10, BF_nonoverlapping_10)

  # H1 true
  BF_10 <- twosd_test(100, 100, 1, 1.5, 0.5)
  BF_overlapping_10 <- twosd_test(100, 100, 1, 1.5, 0.5, null_interval = c(0.90, 1.10))
  BF_nonoverlapping_10 <- twosd_test(
    100, 100, 1, 1.5, 0.5, null_interval = c(0.90, 1.10), nonoverlapping_interval = TRUE
  )

  # More evidence against sharp H0
  expect_gt(BF_10, BF_overlapping_10)

  # Non-overlapping H1 should do better
  expect_gt(BF_nonoverlapping_10, BF_overlapping_10)
})


test_that('Posterior integrates to one', {
  int1 <- integrate(function(x) Vectorize(dphi2)(x, 100, 100, 1, 1), 0, Inf)$value
  int2 <- integrate(function(x) Vectorize(dphi2)(x, 100, 100, 1, 2), 0, Inf)$value
  int3 <- integrate(function(x) Vectorize(dphi2)(x, 100, 100, 2, 1), 0, Inf)$value
  int4 <- integrate(function(x) Vectorize(dphi2)(x, 100, 200, 2, 1), 0, Inf)$value
  expect_true(all(abs(c(int1, int2, int3, int4) - 1) < .001))
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
  expect_true(.is_only_ordered_and_equal('1<2=3=4=5=6=7=8=9<10<11=12'))

  expect_false(.is_only_ordered_and_equal('1<2,3'))
  expect_false(.is_only_ordered_and_equal('1<2=3,4<5'))
  expect_false(.is_only_ordered_and_equal('1<2=3=4,5=6=7=8=9<10<11=12'))
})


test_that('Is all equal / unequal works', {
  hyp0 <- '1=2=3=4=5=6=7=8=9=10=11=12'
  hyp1 <- '1=2,3=4=5=6=7=8=9=10=11=12'
  hyp2 <- '1,2,3,4,5,6,7,8,9,10,11,12'

  expect_true(.is_allequal(hyp0))
  expect_false(.is_allequal(hyp1))
  expect_false(.is_allequal(hyp2))

  expect_true(.is_allunequal(hyp2))
  expect_false(.is_allunequal(hyp1))
  expect_false(.is_allunequal(hyp0))
})


test_that('Prior probability of restriction is correct', {
  hyp <- c('1>2,3,4>5,6>7')
  hyp_prec <- gsub('>', '<', hyp)

  ns <- rep(100, 7)
  sds <- c(1, 2, 2, 5, 6, 6, 7)
  bfvar <- .create_bfvar_object(hyp, ns, sds, 0.50, priors_only = TRUE, iter = 1000, cores = 1, chains = 6)
  rho <- rstan::extract(bfvar$fit, 'rho')$rho

  hyp_fn <- eval(parse(text = .create_hyp_fn(hyp_prec)))
  expect_equal(mean(apply(rho, 1, hyp_fn)), .compute_prior_restr(hyp), tolerance = 0.05)
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

  res <- ksd_test(hyp, ns, sds, u = 0.50, chains = 6, iter = 6000)
  bf10 <- twosd_test(ns[1], ns[2], sds[1], sds[2])
  expect_equal(res$BF[1, 2], bf10, tolerance = 0.05)
})


test_that('ksd_test gives same result as (directed) twosd_test for K = 2', {
  sds <- c(2, 1)
  ns <- c(100, 100)
  hyp <- c('1,2', '1>2')

  bf10 <- twosd_test(ns[1], ns[2], sds[1], sds[2])
  bfr0 <- twosd_test(ns[1], ns[2], sds[1], sds[2], alternative_interval = c(0, 1))
  bfr1 <- bfr0 - bf10

  res <- ksd_test(hyp, ns, sds, u = 0.50, chains = 6, iter = 6000)
  expect_equal(res$BF[2, 1], bfr1, tolerance = 0.05)
})


test_that('ksd_test works for large K', {
  hyp0 <- c('1=2=3=4=5=6=7=8=9=10=11=12=13=14=15=16=17=18=19=20=21')
  hyp1 <- c('1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21')
  hypr <- c('1,2,3,4,5,6=7=8=9=10=11=12=13,14,15,16>17>18>19=20,21')

  ns <- rep(50, 21)
  sds <- rep(10, 21)
  m <- ksd_test(c(hyp0, hyp1, hypr), ns = ns, sds = sds)

  expect_true(m$BF[1, 2] > 0 && m$BF[1, 3] > 0 && m$BF[2, 3] > 0)
})


test_that('Respects evidence bound in ordinal hypotheses', {
  sds <- c(10, 5, 1)
  ns <- c(500, 500, 500)
  hyp <- c('1,2,3', '1>2>3')
  res <- ksd_test(hyp, ns, sds, u = 0.50, chains = 6, iter = 6000)

  expect_equal(res$BF[2, 1], lfactorial(3), tolerance = 0.05)
})


test_that('Mixed equality and ordinal hypotheses make sense I', {
  sds <- c(2, 2, 1)
  ns <- c(500, 500, 500)
  hyp <- c('1=2,3', '1=2>3')
  res <- ksd_test(hyp, ns, sds, u = 0.50, chains = 6, iter = 1000, cores = 1)

  expect_equal(res$BF[2, 1], lfactorial(2), tolerance = 0.05)
})


test_that('Mixed equality and ordinal hypotheses make sense II', {
  sds <- c(4, 4, 4, 1)
  ns <- c(500, 500, 500, 500)
  hyp <- c('1=2=3,4', '1=2=3>4')
  res <- ksd_test(hyp, ns, sds, u = 0.50, chains = 6, iter = 1000, cores = 1)

  expect_equal(res$BF[2, 1], lfactorial(2), tolerance = 0.05)
})


test_that('Mixed equality and ordinal hypotheses make sense III', {
  ns <- rep(100, 5)
  sds <- c(5, 4, 4, 3, 1)
  hyp <- c('1=2=3=4=5', '1,2,3,4,5', '1>2=3>4,5')

  res <- ksd_test(hyp, ns, sds, u = 0.50, chains = 6, iter = 6000)
  expect_gt(res$BF[3, 1], 1)
  expect_gt(res$BF[3, 2], 1)
})


test_that('Evidence ordering makes sense', {
  sds <- c(2, 1.5, 1)
  ns <- c(100, 100, 100)
  hyp <- c('1>2>3', '1,2>3', '1,2,3', '1=2>3', '1=2=3')
  res <- ksd_test(hyp, ns, sds, u = 0.50)

  res$BF <- NULL
  logmls <- sapply(res, function(x) x$logml)
  expect_true(all(logmls == sort(logmls, decreasing = TRUE)))
})
