context("transformations")

test_that("tranformation functions", {
  
  ## logit function
  # Test range
  expect_equal(logit2(0), -Inf)
  expect_equal(logit2(0.5), Inf)
  # Test values
  expect_equal(logit2(0.1), -1.386294, tolerance=1e-6)
  expect_equal(logit2(0.25), 0)
  expect_equal(logit2(0.4), 1.386294, tolerance=1e-6)

  expect_true(is.na(logit2(NA)))
  expect_error(logit2("abc"))
  
  ## inv.logit2 function
  # Test range
  expect_equal(inv.logit2(-Inf), 0)
  expect_equal(inv.logit2(Inf), 0.5)
  # Test values
  expect_equal(inv.logit2(-1.386294), 0.1, tolerance=1e-7)
  expect_equal(inv.logit2(0), 0.25)
  expect_equal(inv.logit2(1.386294), 0.4, tolerance=1e-6)
  
  ## logit function
  # Test range
  expect_equal(logit(0), -Inf)
  expect_equal(logit(1), Inf)
  # Test values
  expect_equal(logit(0.2), -1.386294, tolerance=1e-6)
  expect_equal(logit(0.5), 0)
  expect_equal(logit(0.8), 1.386294, tolerance=1e-6)
  
  ## Test range
  expect_equal(inv.logit2(-Inf), 0)
  expect_equal(inv.logit2(Inf), 0.5)
  ## Test values
  expect_equal(inv.logit2(-1.386294), 0.1, tolerance=1e-7)
  expect_equal(inv.logit2(0), 0.25)
  expect_equal(inv.logit2(1.386294), 0.4, tolerance=1e-6)
  
  expect_equal(inv.logit2(NA), 0.5)
  expect_error(inv.logit2("abc"))
})
