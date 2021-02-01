context("templates")

test_that("macquarie perch arguments can be defined and work with simulate", {

  # simulate from a Macquarie perch object
  dyn <- macquarie_perch()
  arg <- get_args("macquarie_perch")
  sim <- simulate(dyn, args = arg)
  expect_equal(dim(sim), c(1L, 30L, 51L))

})

test_that("macquarie perch arguments can be changed correctly", {

  # simulate from a Macquarie perch object
  dyn <- macquarie_perch()
  arg <- get_args("macquarie_perch")
  sim <- simulate(dyn, args = arg)
  expect_equal(dim(sim), c(1L, 30L, 51L))

  # try adding some YOY and juveniles
  dyn <- macquarie_perch()
  arg <- get_args("macquarie_perch", n = c(5, 1, 0))
  sim <- simulate(dyn, args = arg)
  expect_equal(dim(sim), c(1L, 30L, 51L))

  # change density dependence parameters
  dyn <- macquarie_perch()
  arg <- get_args("macquarie_perch", allee_strength = 2)
  sim <- simulate(dyn, args = arg)
  expect_equal(dim(sim), c(1L, 30L, 51L))

  # and envstoch parameters
  dyn <- macquarie_perch()
  arg <- get_args("macquarie_perch",
                  contributing_min = 0.5,
                  contributing_max = 0.95,
                  recruit_failure = 0.25)
  sim <- simulate(dyn, args = arg)
  expect_equal(dim(sim), c(1L, 30L, 51L))

})

test_that("get_args errors informatively", {

  # check random species name
  expect_error(
    get_args("labradoodle"),
    "labradoodle does not have arguments defined"
  )

})
