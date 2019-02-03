context("test-fit_ezddm")

test_that("fit_ezddm works as intended", {
  library(rtdists)
  data1 <- rdiffusion(n = 100, a = runif(1, 2, 4), v = runif(1, -3, 3), t0 = runif(1, 0.2, 0.5), z = 0.5 * 2) # simulate data
  data2 <- rdiffusion(n = 100, a = runif(1, 2, 4), v = runif(1, -3, 3), t0 = runif(1, 0.2, 0.5), z = 0.5 * 2) # simulate data
  dataAll <- rbind(data1, data2) # join data

  dataAll$response_char <- dataAll$response
  dataAll$response <- ifelse(dataAll$response == "upper", 1, 0) # convert responses to 1 and 0
  dataAll$subject <- rep(c(1, 2), each = 100) # assign subject id

  dataAll$cond1 <- base::sample(c("a", "b"), 200, replace = TRUE) # randomly assign conditions a/b
  dataAll$cond2 <- base::sample(c("y", "z"), 200, replace = TRUE) # randomly assign conditions y/z

  # fit model to just entire data set (assumes all data came from 1 subject)
  a <- fit_ezddm(data = dataAll, rts = "rt", responses = "response")
  expect_equal(a[, .N], 1)
  # fit model to just entire data set (assumes all data came from 1 subject)
  b <- fit_ezddm(data = dataAll, rts = "rt", responses = "response", id = "subject")
  expect_equal(b[, .N], 2)
  # fit model to each subject by cond1
  c <- fit_ezddm(data = dataAll, rts = "rt", responses = "response", id = "subject", group = "cond1")
  expect_equal(c[, .N], 4)
  # fit model to each subject by cond1,cond2
  d <- fit_ezddm(data = dataAll, rts = "rt", responses = "response", id = "subject", group = c("cond1", "cond2"))
  expect_equal(d[, .N], 8)

  dataAll$rt <- dataAll$rt * runif(nrow(dataAll), 100, 2000)
  expect_error(fit_ezddm(data = dataAll, rts = "rt", responses = "response"))

})
