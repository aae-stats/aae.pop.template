#' @name macquarie_perch
#' @title Parameterised Macquarie perch population model
#' @description Use a pre-defined population dynamics object to
#'   simulate population dynamics of Macquarie perch
#'   (*Macquaria australasica*). Model parameters are based
#'   on existing data sets and published studies on Macquarie
#'   perch.
NULL

#' @rdname macquarie_perch
#'
#' @importFrom stats pnorm rnorm runif
#' @importFrom Rdpack reprompt
#'
#' @export
#'
#' @param k integer carrying capacity (maximum population size)
#'   for the simulated Macquarie perch population. Defaults to 1000
#' @param reproductive integer vector defining the adult stages in the
#'   Macquarie perch population. Defaults to \code{3:30}
#' @param system one of \code{"lake"} or \code{"river"}, defining
#'   whether covariate associations are based on lake or river
#'   conditions. Covariate associations are included in simulations
#'   only if covariates are provided to \code{\link[aae.pop]{simulate}}
#'   through the \code{args} argument
#' @param genetic_factor positive real value representing proportional
#'   changes in early life (egg, larvae, young-of-year) survival,
#'   included here to allow scenarios of genetic mixing of Macquarie
#'   perch populations (e.g., \insertCite{lutz20}{aae.pop.templates}).
#'   Note that non-default values of \code{genetic_factor} must be
#'   provided when defining the population model with
#'   \code{macquarie_perch} and again when defining arguments with
#'   \code{get_args("macquarie_perch")}
#'
#' @details The \code{macquarie_perch} population model is an
#'   age-structured model with 30 age classes and includes negative
#'   and positive density dependence, environmental and demographic
#'   stochasticity, and optional associations with hydrological conditions
#'   in rivers or lakes (controlled with the \code{system} parameter)
#'   and stocking or angling effects.
#'
#'   The full population model is developed and described in
#'   \insertCite{todd_lint15}{aae.pop.templates}
#'
#' @references
#' \insertRef{todd_lint15}{aae.pop.templates}
#'
#' @examples
#' # define a basic model for Macquarie perch with
#' #   carrying capacity = 1000
#' mp <- macquarie_perch(k = 1000)
#'
#' # define required arguments
#' mp_args <- get_args("macquarie_perch")
#'
#' # simulate from this model
#' sims <- simulate(mp, nsim = 100, args = mp_args)
#'
#' # plot the simulated values
#' plot(sims)
macquarie_perch <- function(
  k = 1000,
  reproductive = 3:30,
  system = "lake",
  genetic_factor = 1.0
) {
  get_template(
    sp = "macquarie_perch",
    k = k,
    reproductive = reproductive,
    system = system)
}

# internal function: define species defaults
template_macquarie_perch <- function(
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
    function(pop, n_yoy, add_yoy, ...)
      add_remove(pop = pop, n = n_yoy, add = add_yoy),
    function(pop, n_twoplus, add_twoplus, ...)
      add_remove(pop = pop, n = n_twoplus, add = add_twoplus),
    function(pop, n_adult, add_adult, ...)
      add_remove(pop = pop, n = n_adult, add = add_adult)
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

#' @rdname macquarie_perch
#'
#' @importFrom stats pnorm rnorm runif
#'
#' @export
#'
#' @param n an integer, vector of integers, or list specifying the
#'   number of young-of-year, 2+, and adult fish to add or
#'   remove in any given year. If a single integer is provided,
#'   all age classes are assumed to have the same number of
#'   additions or removals. If a vector is provided, it must
#'   have one value for each age class. If a list is provided,
#'   it must have one element for each age class, with a
#'   value for each time step. Addition versus removal is
#'   controlled with \code{add}. Defaults to
#'   \code{c(0, 0, 0))}, which will specify no additions
#'   or removals
#' @param ntime number of time steps used in population
#'   simulation. Defaults to \code{50} and is required
#'   to ensure simulated additions or removals are defined
#'   for every simulated time step
#' @param start time step at which additions or removals start
#'   if \code{n} is an integer or vector. Defaults to
#'   \code{c(1, 1, 1)}
#' @param end time step at which additions or removals finish
#'   if \code{n} is an integer or vector. Defaults to
#'   \code{c(1, 1, 1)}
#' @param add logical indicating whether individuals are
#'   added or removed from the population. Defaults to
#'   \code{TRUE}, which defines additions (e.g., stocking). Can
#'   be specified as a single value, in which case all stages
#'   are set equal, as three values, in which case stages can
#'   differ, or as a matrix with three rows and \code{ntime}
#'   columns, in which case the effects of \code{add} can
#'   change through time
#' @param allee_strength strength of Allee effect. Defaults
#'   to \code{1}. See \insertCite{todd_lint15}{aae.pop.templates}
#'   for details
#' @param contributing_min minimum proportion of adult females
#'   contributing to reproduction at any given time step.
#'   Defaults to \code{1.0}
#' @param contributing_max maximum proportion of adult females
#'   contributing to reproduction at any given time step.
#'   Defaults to \code{1.0}
#' @param recruit_failure probability of complete recruitment
#'   failure at any given time step. Defaults to \code{0}
#' @param genetic_factor proportional change in early life
#'   survival due to genetic mixing
#'
#' @details The Macquarie perch population template requires
#'   several additional arguments and allows several optional
#'   arguments. Arguments can be defined with a call to
#'   \code{get_args("macquarie_perch")} and are described
#'   individually above.
#'
args_macquarie_perch <- function(
  n = c(0, 0, 0),
  ntime = 50,
  start = c(1, 1, 1),
  end = c(1, 1, 1),
  add = TRUE,
  allee_strength = 1,
  contributing_min = 0.75,
  contributing_max = 1.0,
  recruit_failure = 0,
  genetic_factor = 1.0
) {

  # expand n, start, end, add if required
  if (length(n) == 1)
    n <- rep(n, 3)
  if (length(start) == 1)
    start <- rep(start, 3)
  if (length(end) == 1)
    end <- rep(end, 3)
  if (length(add) == 1)
    add <- rep(add, 3)

  # check for other lengths of n, start, or end
  if (any(c(length(n), length(start), length(end)) != 3)) {
    stop("n, start, and end must be vectors with 1 or 3 elements",
         call. = FALSE)
  }

  # helper to define removals process
  define_removals <- function(start, end, n, add = add) {

    # set up a sequence of iterations at which individuals are removed
    if (!is.list(n)) {
      n <- mapply(
        zeros_and_fill,
        n,
        start,
        end,
        MoreArgs = list(len = ntime),
        SIMPLIFY = FALSE
      )
    } else {
      if (!all(sapply(n, length) == ntime)) {
        stop("if n is a list, each element must be a vector ",
             "with one value for each time step",
             call. = FALSE)
      }
    }

    # set up a sequence of add flags
    if (!is.matrix(add)) {
      add <- matrix(rep(add, ntime), nrow = 3)
    } else {
      if (nrow(add) != 3 | ncol(add) != ntime) {
        stop("if add is a matrix, it must have three rows ",
             "and ntime columns (ntime = ", ntime, ")",
             call. = FALSE)
      }
    }

    # define this as a function
    translocate <- function(obj, pop, iter) {

      # return
      list(
        n_yoy = n[[1]][iter],
        n_twoplus = n[[2]][iter],
        n_adult = n[[3]][iter],
        add_yoy = add[1, iter],
        add_twoplus = add[2, iter],
        add_adult = add[3, iter]
      )

    }

    # return
    translocate

  }

  early_surv <- c(0.5, 0.013, 0.13)

  # helper to calculate real-valued parameters for survival
  #   simulation
  transform_survival <- function(obj, pop, iter) {

    # pull out the population matrix in the current time step
    mat <- obj$matrix
    if (is.list(mat))
      mat <- mat[[iter]]

    # wrap up all survival means and SDs, including early life
    #  (this allows a single call to unit_to_real, which is slow)
    survival_mean <- c(
      genetic_factor * early_surv,  # early life survival with gene mixing
      mat[transition(mat)]    # from population matrix in current time step
    )

    survival_sd <- c(
      0.1, 0.007, 0.028,  # early life
      0.05, 0.09, 0.11, 0.10, 0.10, 0.07, 0.08, 0.08, 0.08,
      0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,
      0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.07, 0.06, 0.05
    )

    # convert unit interval to real line equivalents
    out <- unit_to_real(
      unit_mean = survival_mean,
      unit_sd = survival_sd
    )

    # separate early life from other estimates
    idx <- seq_len(nrow(out)) > 3

    # return
    list(
      mean_real = out[idx, 1],    # for survival_gen
      sd_real = out[idx, 2],      # for survival_gen
      early_mean = out[!idx, 1],  # for reproduction_gen
      early_sd = out[!idx, 2]     # for reproduction_gen
    )

  }

  # return named list of args
  list(

    # set as 1 (default) or 2
    density_dependence = list(allee_strength = allee_strength),

    # set contributing as random uniform on 0.75-1.0 by default
    # set recruit_failure at 0 by default
    # add function to pre-transform unit to real and back
    environmental_stochasticity = list(
      contributing_min = contributing_min,
      contributing_max = contributing_max,
      recruit_failure = recruit_failure,
      transform_survival
    ),

    # to include additions or removals of individuals
    density_dependence_n = list(
      define_removals(
        start = start, end = end, n = n, add = add
      )
    )

  )

}
