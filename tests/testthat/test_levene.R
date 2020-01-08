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
  expect_true(two_sample(100, 100, 1, 1, 0.5) < 0)
  expect_true(two_sample(100, 100, 1, 2, 0.5) > 0)
})
