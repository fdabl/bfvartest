library('BayesLevene')


context('K = 1 Sample Test')
test_that('One-sample Test Gives Expected Result', {
  expect_lt(one_sample(100, 1, 1, 0.5), 0) # log BF10 < 0
  expect_gt(one_sample(100, 1, 2, 0.5), 0) # log BF10 > 0
})


context('K = 2 Sample Test')
test_that('Function stops if wrong interval input', {
  expect_error(two_sample(100, 100, 1, 1, 0.5, alternative_interval = NULL))
  expect_error(two_sample(100, 100, 1, 1, 0.5, alternative_interval = c(-2, 4)))

  expect_error(two_sample(100, 100, 1, 1, 0.5, null_interval = c(-2, 4)))
  expect_error(two_sample(100, 100, 1, 1, 0.5, null_interval = c(.5, .2)))
})


test_that('Two-sample Test Gives Expected Result', {
  expect_true(Rmpfr::asNumeric(two_sample(100, 100, 1, 1, 0.5)) < 0)
  expect_true(Rmpfr::asNumeric(two_sample(100, 100, 1, 2, 0.5)) > 0)

  # 2 is larger than 1
  expect_true(Rmpfr::asNumeric(two_sample(100, 100, 1, 2, 0.5, alternative_interval = c(1, Inf))) > 0)

  # 1 is larger than 2
  expect_true(Rmpfr::asNumeric(two_sample(100, 100, 1, 2, 0.5, alternative_interval = c(0, 1))) < 0)
})


context('K > 2 Sample Helper Functions')
test_that('Constrained matrix works', {
  expect_error(.create_constraint_matrix('2>1')) # needs to be ascending
  expect_error(.create_constraint_matrix('1<2')) # < not allowed
  expect_error(.create_constraint_matrix('1u2,3>4')) # u not allowed
})


test_that('Group sizes K > 9 work', {
  hyp <- c('1>2>3>4>5>6>7>8>9>10>11>12')
  lower <- c(-1, seq(11))
  upper <- c(seq(2, 12), -2)
  expect_true(all(.create_constraint_matrix(hyp) == cbind(lower, upper)))
})



context('K > 2 Sample Test')
test_that('k_sample errors correctly', {
  ss <- c(1, 2, 3)
  ns <- c(100, 100, 100)

  expect_error(k_sample('1>2>3', ns, ss, alpha = 0.50))
  expect_error(k_sample(c('1>2>3', '1=2=3'), ns[-1], ss, alpha = 0.50))
})

test_that('k_sample gives same result as (undirected) two_sample for K = 2', {
  ss <- c(1, 2)
  prec <- 1/ss
  ns <- c(100, 100)
  hyp <- c('1,2', '1=2')

  res <- k_sample(hyp, ns, prec, alpha = 0.50)
  bf10 <- Rmpfr::asNumeric(two_sample(ns[1], ns[2], ss[1], ss[2]))
  expect_true(abs(res$BF[1, 2] - bf10) < .1)
})


test_that('k_sample gives same result as (directed) two_sample for K = 2', {
  ss <- c(1, 2)
  prec <- 1/ss
  ns <- c(100, 100)
  hyp <- c('1,2', '1=2', '1>2')

  bf10 <- Rmpfr::asNumeric(two_sample(ns[1], ns[2], ss[1], ss[2]))
  bfr0 <- Rmpfr::asNumeric(two_sample(ns[1], ns[2], ss[1], ss[2], alternative_interval = c(1, Inf)))
  bfr1 <- bfr0 - bf10

  res <- k_sample(hyp, ns, prec, alpha = 0.50)
  expect_true(abs(res$BF[1, 2] - bf10) < .01)
})



  # tau <- extract(res[[1]]$fit, 'tau')$tau
  # rho <- extract(res[[1]]$fit, 'rho')$rho
  #
  # tau1 <- 2 * tau * rho[, 1]
  # sigma1 <- sqrt(1 / tau1)
  # delta <- sqrt(rho[, 1] / rho[, 2])

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
