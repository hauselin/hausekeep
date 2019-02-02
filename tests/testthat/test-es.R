context("test-es")

test_that("test d)", {
  y <- es(d = 0.5)
  expect_equal(y$r, 0.243)
})
