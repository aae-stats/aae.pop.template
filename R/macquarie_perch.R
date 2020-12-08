#' @rdname macquarie_perch
#'
#' @export
#'
#' @importFrom stats pnorm rnorm runif
macquarie_perch <- function(
  k = 1000,
  reproductive = 3:30,
  system = "lake",
  ...
) {
  get_template(
    sp = "macquarieperch",
    k = k,
    reproductive = reproductive,
    system = system,
    ...)
}

# internal function: define species defaults
template_macquarieperch <- function(
  k = 1000,                     # adult female carrying capacity
  reproductive = 3:30,          # reproductive age classes
  system = "lake",              # define covariate type by system
  genetic_factor = 1.0          # define change in early survival
) {

  # set default system
  if (length(system) > 1)
    system <- system[1]

  # and check it
  if (!system %in% c("lake", "river"))
    stop("system must be one of lake or river", call. = FALSE)

  # define a survival function, adding dots to soak up extra
  #    arguments within simulate
  survival_gen <- function(
    mat, mean_real, sd_real, perfect_correlation = TRUE, ...
  ) {
    rmultiunit_from_real(
      n = 1,
      mean_real = mean_real,
      sd_real = sd_real,
      perfect_correlation = perfect_correlation
    )
  }

  # define a reproduction function, being careful with argument
  #   names to avoid conflicts with any arguments in
  #   survival_gen, which would then require multiple different
  #   arguments with the same name in simulate
  reproduction_gen <- function(
    mat,
    fec_mean = c(1.68, -0.302, 2.886),
    fec_sd = c(0.3, 0.05, 0.15),
    early_mean,
    early_sd,
    recruit_failure = 0.25,
    contributing_min = 0.5,
    contributing_max = 1.0,
    ...
  ) {

    # generate stochastic values for early life
    #   survival (eggs, larvae, young-of-year)
    early_surv <- rmultiunit_from_real(n = 1, mean = early_mean, sd = early_sd)

    # otherwise draw random variates for the three model parameters
    y1 <- rnorm(n = 1, mean = fec_mean[1], sd = fec_sd[1])
    y2 <- rnorm(n = 1, mean = fec_mean[2], sd = fec_sd[2])
    y3 <- rnorm(n = 1, mean = fec_mean[3], sd = fec_sd[3])

    # generate reproduction estimates for all adult age classes, incorporating
    #   stochastic early life estimates
    y2_term <- exp(y2 %o% reproductive)
    y1_y2 <- log(
      43.15 * exp(sweep(y2_term, 1, -y1, "*"))
    )
    reprod <- exp(sweep(2.295 * y1_y2, 1, y3, "+"))

    # did recruitment fail?
    recruit_binary <- ifelse(recruit_failure >= runif(1), 0, 1)

    # make contributing a stochastic variables
    if (contributing_min > contributing_max) {
      contributing_min <- contributing_max
      warning(
        "contributing_min was greater than contributing_max; ",
        "both values have been set to contributing_max",
        call. = FALSE
      )
    }
    contributing <- runif(1, min = contributing_min, max = contributing_max)

    # add early life survival and muliply by 0.5
    #   to account for a 50:50 sex ratio
    0.5 * recruit_binary * contributing * reprod * prod(early_surv)

  }

  # define mean survival
  survival_mean <- c(
    0.25, 0.44, 0.56, 0.63, 0.69, 0.72, 0.75, 0.78, 0.79,
    0.81, 0.82, 0.83, 0.83, 0.84, 0.84, 0.84, 0.85, 0.85, 0.84,
    0.84, 0.84, 0.83, 0.82, 0.80, 0.78, 0.76, 0.71, 0.63, 0.48
  )

  # calculate reproduction at mean parameter values
  reproduction_mean <- exp(
    2.295 *
      log(43.15 * exp(- 1.68 * exp(-0.302 * reproductive))) +
      2.886) *
    0.5 *                                       # 50:50 sex ratio
    prod(genetic_factor * c(0.5, 0.013, 0.13))  # add early life survival

  # define population matrix
  nclass <- length(survival_mean) + 1
  popmat <- matrix(0, nrow = nclass, ncol = nclass)
  popmat[transition(popmat)] <- survival_mean
  popmat[reproduction(popmat, dims = reproductive)] <- reproduction_mean

  # define density dependence, only affects adult survival
  #   and reproductive stages
  density_masks <- list(
    transition(popmat, dims = 2:30),
    reproduction(popmat, dims = reproductive)
  )

  # top-down effects of competition for habitat, at
  #   carrying capacity k
  topdown_fn <- function(mat, pop, ...) {
    sum_n <- sum(pop[reproductive])
    ifelse(sum_n > k, k / sum_n, 1) * mat
  }

  # positive density dependence (Allee effect)
  allee_fn <- function(mat, pop, allee_strength = 1, allee_factor = 10, ...) {
    sum_n <- sum(pop[reproductive])
    allee <- (2 / (1 + exp(-sum_n / (allee_strength * allee_factor)))) - 1
    allee * mat
  }

  # and collate masks and functions in a single object
  dens_depend <- density_dependence(
    masks = density_masks,
    funs = list(topdown_fn, allee_fn)
  )

  # define environmental stochasticity
  envstoch <- environmental_stochasticity(
    masks = list(
      transition(popmat),
      reproduction(popmat, dims = reproductive)
    ),
    funs = list(survival_gen, reproduction_gen)
  )

  # use density_dependence_n to include stocking or
  #   translocations (removals)
  dd_n_masks <- list(
    all_classes(popmat, dim = 1),
    all_classes(popmat, dim = 2),
    all_classes(popmat, dim = reproductive)
  )
  dd_n_fns <- list(
    function(pop, n_yoy, ...)
      add_remove(pop = pop, n = n_yoy, add = TRUE),
    function(pop, n_twoplus, ...)
      add_remove(pop = pop, n = n_twoplus, add = TRUE),
    function(pop, n_adult, ...)
      add_remove(pop = pop, n = n_adult, add = TRUE)
  )
  dens_depend_n <- density_dependence_n(
    masks = dd_n_masks,
    funs = dd_n_fns
  )

  # define covariate effects on recruitment in all systems
  recruit_effects <- list(
    function(
      mat, x, spawning_param = c(-0.01, -0.05), variability_param = -0.003, ...
    ) {

      # effect of spawning-flow magnitude on recruitment
      # negative effect of Nov/Dec discharge on recruitment
      log_flow <- log(x$spawning_flow + 0.01)
      scale_factor <- exp(
        spawning_param[1] * log_flow + spawning_param[2] * (log_flow ^ 2)
      )
      scale_factor[scale_factor > 1] <- 1
      scale_factor[scale_factor < 0] <- 0
      mat <- mat * scale_factor

      # effect of spawning-flow variability on recruitment
      # negative effect of Nov/Dec discharge variability on recruitment
      #   (days with more than 100% change from previous)
      mat <- mat * exp(-variability_param * x$spawning_variability)

      # return
      mat

    }
  )

  # define covariate effects on adult survival in all systems
  survival_effects <- NULL

  # define covariate masks in all systems
  recruit_masks <- list(reproduction(popmat))
  survival_masks <- NULL

  # add system-specific covariate effects
  if (system == "lake") {

    # negative effect of rising lake level on YOY
    recruit_effects_lake <- function(
      mat, x, recruit_param = -0.5, shift = 10, ...
    ) {
      mat * (1 / (1 + exp(-recruit_param * (x$water_level_change + shift))))
    }

    # update effects and masks
    recruit_effects <- c(recruit_effects, list(recruit_effects_lake))
    recruit_masks <- c(recruit_masks, list(reproduction(popmat)))

  }
  if (system == "river") {

    # negative effect of rising river level on YOY
    recruit_effects_river <- function(
      mat, x, recruit_param = -0.01, shift = 200, ...
    ) {
      mat * (1 / (1 + exp(recruit_param * (x$river_height_change + shift))))
    }

    # negative effect of low flows on adult survival
    survival_effects_river <- function(
      mat, x, survival_param = c(0.2, -0.2), ...
    ) {

      # positive effect of flow on overall population growth rate
      #    (based on individuals 1 year or older)
      log_flow <- log(x$average_daily_flow + 0.01)
      scale_factor <- exp(
        survival_param[1] * log_flow + survival_param[2] * (log_flow ^ 2)
      )
      scale_factor[scale_factor > 1] <- 1
      scale_factor[scale_factor < 0] <- 0
      mat <- mat * scale_factor

      # return
      mat

    }

    # update effects and masks
    recruit_effects <- c(recruit_effects, list(recruit_effects_river))
    recruit_masks <- c(recruit_masks, list(reproduction(popmat)))
    survival_effects <- c(survival_effects, list(survival_effects_river))
    survival_masks <- c(
      survival_masks, list(transition(popmat, dims = 2:30))
    )

  }

  # compile covariates process
  covars <- covariates(
    masks = c(recruit_masks, survival_masks),
    funs = c(recruit_effects, survival_effects)
  )

  # nolint start
  # print a message to ensure args are included
  message('Arguments are required to simulate Macquarie perch dynamics.\n',
          'Default arguments can be accessed with get_args("macquarie_perch").')
  # nolint end

  # return
  list(
    matrix = popmat,
    covariates = covars,
    environmental_stochasticity = envstoch,
    demographic_stochasticity = NULL,
    density_dependence = dens_depend,
    density_dependence_n = dens_depend_n
  )

}
