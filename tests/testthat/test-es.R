context("test-es")

test_that("es() returns correct values", {
  d <- 0.5
  y <- es(d = d)
  expect_equal(y$r, 0.243)
  expect_equal(round(es(r = y$r)$d, 1), d)
  expect_equal(es(f = 0.2)$d, 0.4)
  expect_equal(es(r = 0.2)$R2, 0.2^2)
  expect_equal(es(R2 = 0.25)$r, sqrt(0.25))
  expect_error(es())
  expect_error(es(c(d = 0.2, a = 0.3)))
})
