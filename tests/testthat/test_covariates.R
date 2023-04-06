context("covariates")

# define some covariates applicable to all species
set.seed(525110129)
n <- 25
covars <- list(

  macquarie_perch = data.frame(
    river_height_change = runif(n, 0.5, 2),
    average_daily_flow = rlnorm(n),
    min_daily_flow = rlnorm(n, mean = -1),
    water_level_change = runif(n, 0.5, 2),
    spawning_variability = runif(n, 0, 2),
    spawning_flow = rlnorm(n),
    temperature_effect = runif(n, 0.7, 1.2)
  ),

  australian_bass = data.frame(
    water_temperature = rlnorm(n),
    discharge = rlnorm(n)
  ),

  barred_galaxias = data.frame(
    riparian = runif(n, 0.5, 1.5),
    bushfire = sample(c(0, 1), size = n, replace = TRUE),
    ctf = sample(c(0, 1), size = n, replace = TRUE),
    trout = sample(c(0, 1), size = n, replace = TRUE)
  ),

  common_carp = data.frame(
    floodplain_access = sample(c(0, 1), size = n, replace = TRUE),
    flow_variability = rlnorm(n)
  ),

  estuary_perch = data.frame(
    x  = runif(n, 0.2, 1)
  ),

  murray_cod = data.frame(
    minimum_daily_flow = rlnorm(n, mean = -1),
    spawning_flow_variability = rlnorm(n),
    proportional_spring_flow = rlnorm(n),
    proportional_max_antecedent = rlnorm(n),
    proportional_summer_flow = rlnorm(n, mean = -0.5),
    proportional_winter_flow = rlnorm(n),
    spawning_temperature = rlnorm(n)
  ),

  murray_rainbowfish = data.frame(
    instream_cover = runif(n, 0, 1),
    nday_gt20 = rpois(n, lambda = 50),
    nday_lt10 = rpois(n, lambda = 5),
    gambusia = sample(c(0, 1), size = n, replace = TRUE),
    redfin = sample(c(0, 1), size = n, replace = TRUE),
    spawning_flow_variability = rnorm(n),
    proportional_spring_flow = rnorm(n),
    proportional_summer_flow = rnorm(n)
  ),

  platypus = data.frame(
    proportional_spring_flow = rlnorm(n),
    proportional_summer_flow = rlnorm(n, mean = -0.5),
    proportional_winter_flow = rlnorm(n),
    proportional_maximum_flow = rlnorm(n, mean = 1),
    spawning_flow_variability = rlnorm(n),
    ctf_duration = rpois(n, lambda = 20)
  ),

  pygmy_perch = data.frame(
    habitat = runif(n, 0.3, 0.9),
    predators = sample(c(0, 1), size = n, replace = TRUE),
    dry = sample(c(0, 1), size = n, replace = TRUE)
  ),

  river_blackfish = data.frame(
    instream_cover = runif(n, 0, 1),
    nday_gt16 = rpois(n, lambda = 50),
    nday_lt5 = rpois(n, lambda = 5),
    spawning_flow_variability = rnorm(n),
    proportional_spring_flow = rnorm(n),
    proportional_summer_flow = rnorm(n),
    proportional_winter_flow = rnorm(n),
    antecedent_flow = rnorm(n)
  )

)

test_that("simulate works with covariates", {

  # simulate from a Macquarie perch object with covariates
  dyn <- macquarie_perch()
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(covariates = format_covariates(covars$macquarie_perch))
  )
  expect_equal(dim(sim), c(10L, 30L, n + 1))

  # repeat but for the river model
  dyn <- macquarie_perch(system = "river")
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(covariates = format_covariates(covars$macquarie_perch))
  )
  expect_equal(dim(sim), c(10L, 30L, n + 1))

  # simulate from an Australian bass object with covariates
  dyn <- australian_bass()
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(covariates = format_covariates(covars$australian_bass))
  )
  expect_equal(dim(sim), c(10L, 50L, n + 1))

  # simulate from a barred galaxias object with covariates
  dyn <-  barred_galaxias()
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(covariates = format_covariates(covars$barred_galaxias))
  )
  expect_equal(dim(sim), c(10L, 4L, n + 1))

  # simulate from a common carp object with covariates
  dyn <- common_carp()
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(covariates = format_covariates(covars$common_carp))
  )
  expect_equal(dim(sim), c(10L, 28L, n + 1))

  # simulate from an estuary perch object with covariates
  dyn <- estuary_perch()
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(covariates = format_covariates(covars$estuary_perch))
  )
  expect_equal(dim(sim), c(10L, 40L, n + 1))

  # simulate from a murray cod object with covariates
  dyn <- murray_cod()
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(
      covariates = c(
        format_covariates(covars$murray_cod),
        list(system = "murray_river")
      )
    )
  )
  expect_equal(dim(sim), c(10L, 50L, n + 1))

  # repeat murray cod sims for different covariate systems
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(
      covariates = c(
        format_covariates(covars$murray_cod),
        list(system = "broken_river")
      )
    )
  )
  expect_equal(dim(sim), c(10L, 50L, n + 1))
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(
      covariates = c(
        format_covariates(covars$murray_cod),
        list(system = "ovens_river")
      )
    )
  )
  expect_equal(dim(sim), c(10L, 50L, n + 1))
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(
      covariates = c(
        format_covariates(covars$murray_cod),
        list(system = "goulburn_river")
      )
    )
  )
  expect_equal(dim(sim), c(10L, 50L, n + 1))
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(
      covariates = c(
        format_covariates(covars$murray_cod),
        list(system = "king_river")
      )
    )
  )
  expect_equal(dim(sim), c(10L, 50L, n + 1))
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(
      covariates = c(
        format_covariates(covars$murray_cod),
        list(system = "broken_creek")
      )
    )
  )
  expect_equal(dim(sim), c(10L, 50L, n + 1))

  # and with manually set covar effects
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(
      covariates = c(
        format_covariates(covars$murray_cod),
        list(coefs = c(-10, 30, 50, -25, 50, 10))
      )
    )
  )
  expect_equal(dim(sim), c(10L, 50L, n + 1))

  # simulate from a rainbowfish object with covariates
  dyn <- murray_rainbowfish()
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(
      covariates = format_covariates(covars$murray_rainbowfish)
    )
  )
  expect_equal(dim(sim), c(10L, 5L, n + 1))

  # and with different coefficients for covariate effects
  dyn <- murray_rainbowfish()
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(
      covariates = c(
        format_covariates(covars$murray_rainbowfish),
        list(coefs = c(200, 100, 20))
      )
    )
  )
  expect_equal(dim(sim), c(10L, 5L, n + 1))

  # simulate from a river blackfish object with covariates
  dyn <- river_blackfish()
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(covariates = format_covariates(covars$river_blackfish))
  )
  expect_equal(dim(sim), c(10L, 11L, n + 1))

  # and with different coefficients for covariate effects
  dyn <- river_blackfish()
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(
      covariates = c(
        format_covariates(covars$river_blackfish),
        list(coefs = c(200, 100, 20, 5, 10))
      )
    )
  )
  expect_equal(dim(sim), c(10L, 11L, n + 1))

  # simulate from a platypus object with covariates
  dyn <- platypus()
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(covariates = format_covariates(covars$platypus))
  )
  expect_equal(dim(sim), c(10L, 2L, n + 1))

  # simulate from a pygmy perch object with covariates
  dyn <- pygmy_perch()
  sim <- simulate(
    dyn,
    nsim = 10,
    args = list(covariates = format_covariates(covars$pygmy_perch))
  )
  expect_equal(dim(sim), c(10L, 4L, n + 1))

})

test_that("simulate errors informatively when covariates are specified
  incorrectly", {

  # simulate from a murray cod object with covariates
  dyn <- murray_cod()
  sim <- expect_error(
    simulate(
      dyn,
      nsim = 10,
      args = list(
        covariates = c(
          format_covariates(covars$murray_cod),
          list(coefs = c(-10, 30, 50, -25, 50, 10, 5))
        )
      )
    ),
    "coefs must include six values"
  )

  # simulate from a rainbowfish object with covariates
  dyn <- murray_rainbowfish()
  sim <- expect_error(
    simulate(
      dyn,
      nsim = 10,
      args = list(
        covariates = c(
          format_covariates(covars$murray_rainbowfish),
          list(coefs = c(200, 100, 20, 100))
        )
      )
    ),
    "coefs must include three values"
  )

  # simulate from a river blackfish object with covariates
  dyn <- river_blackfish()
  sim <- expect_error(
    simulate(
      dyn,
      nsim = 10,
      args = list(
        covariates = c(
          format_covariates(covars$river_blackfish),
          list(coefs = c(200, 100, 20, 5, 10, 10))
        )
      )
    ),
    "coefs must include five values"
  )

})
