#' @name murray_cod
#' @title Parameterised Murray cod population model
#' @description Use a pre-defined population dynamics object to
#'   simulate population dynamics of Murray cod
#'   (*Maccullochella peelii*). Model parameters are based
#'   on existing data sets and published studies on Murray
#'   cod.
NULL

#' @rdname murray_cod
#'
#' @importFrom stats rnorm
#'
#' @export
#'
#' @param k carrying capacity
#' @param reproductive integer vector defining the adult stages in the
#'   Murray cod population. Defaults to \code{5:50}
#' @param system one of \code{"murray"}, \code{"largetrib"}, or
#'   \code{"smalltrib"}, defining whether the maximum size of
#'   individuals matches observations in the Murray River (1300 mm),
#'   large tributaries (e.g., Goulburn River, Campaspe River; 1100 mm),
#'   or small tributaries (e.g., XYZ; 900 mm). Defaults to "murray"
#'
#' @details The \code{murray_cod} population model is a mixed
#'   age- and stage-structured model with 25 classes and includes negative
#'   density dependence, environmental and demographic
#'   stochasticity, and optional associations with hydrological conditions
#'   in rivers and stocking or angling effects.
#'
#' @examples
#' # define a basic model for Murray cod
#' mc <- murray_cod()
#'
#' # define required arguments
#' mc_args <- get_args("murray_cod")
#'
#' # simulate from this model
#' sims <- simulate(mc, nsim = 100, args = mc_args)
#'
#' # plot the simulated values
#' plot(sims)
murray_cod <- function(k = 20000, system = "murray") {
  get_template(sp = "murray_cod", k = k, system = system)
}

# internal function: define species defaults
template_murray_cod <- function(
  k = 20000,
  reproductive = 5:50,          # reproductive age classes
  system = "murray"
) {

  # how many stages are we going to work with?
  nstage <- 50

  # check system
  if (!system %in% c("murray", "largetrib", "smalltrib")) {
    stop(
      "system must be one of murray, largetrib or smalltrib",
      call. = FALSE
    )
  }

  # set system-specific max lengths
  max_length_options <- c(
    "murray" = 1300,
    "largetrib" = 1100,
    "smalltrib" = 900
  )
  max_length <- max_length_options[system]

  # force evaluation of max length
  force(max_length)

  # helper function to calculate length from age (stochastically)
  estimate_length <- function(age) {

    # some length/growth parameters
    l_inf <- 1360.465
    k_est <- 0.0674 + rbeta(1, 849.86, 884.55) - 0.49
    max_lt <- rnorm(1, max_length, 40 * (max_length / l_inf))

    # return estimate lengths
    max_lt * (1 - exp(-k_est * (age + 1.535)))

  }

  # helper function to calculate weight from length (stochastically)
  estimate_weight <- function(length_est) {

    # length-to-weight parameter
    a_est <- rnorm(1, -13.52, 0.1352)

    # estimate and return weight
    (exp(a_est) / 1000) * (length_est ^ 3.36)

  }

  # define  a survival function
  survival_gen <- function(
    mat,
    mean_real,
    sd_real,
    perfect_correlation = TRUE,
    ...
  ) {

    rmultiunit_from_real(
      n = 1,
      mean_real = mean_real,
      sd_real = sd_real,
      perfect_correlation = perfect_correlation
    )
  }

  # define a reproduction function
  reproduction_gen <- function(
    mat,
    fec_params = c(389, 2723, 5344, 53.44),
    early_mean,
    early_sd,
    ...
  ) {

    # estimate weights from ages
    length_est <- estimate_length(reproductive)
    weight_est <- estimate_weight(length_est)

    # simulate parameters for fecundity from weight conversion
    x <- rnorm(n = 1, mean = fec_params[1], sd = fec_params[2])
    y <- rnorm(n = 1, mean = fec_params[3], sd = fec_params[4])

    # generate stochastic values for early life
    #   survival (eggs, larvae, young-of-year)
    early_surv <- rmultiunit_from_real(n = 1, mean = early_mean, sd = early_sd)

    # estimate baseline, per-capita fecundity
    reprod <- -x + y * weight_est - 69.5 * (weight_est ^ 2)

    # estimate and return reproduction estimates
    0.5 * reprod * prod(early_surv)

  }

  # define mean survival
  survival_params <- c(1107.608, 1.2668, -7.6079)
  mean_surv <- (survival_params[1] / (c(1:(nstage)) ^ survival_params[2])) +
    survival_params[3]
  mean_surv <- mean_surv[2:(nstage)] / mean_surv[1:(nstage - 1)]
  mean_surv[mean_surv < 0] <- 0

  # define base matrix
  popmat <- matrix(0, nrow = nstage, ncol = nstage)

  # add survival
  popmat[transition(popmat)] <- mean_surv

  # add reproduction
  early_surv <- c(0.5, 0.012, 0.38, 0.31)
  reproduction_mask <- reproduction(popmat, dims = reproductive)
  mean_lengths <- max_length * (1 - exp(-0.0674 * (reproductive + 1.535)))
  mean_weights <- (exp(-13.52) / 1000) * (mean_lengths ^ 3.36)
  reprod <- -389 + 5344 * mean_weights - 69.5 * (mean_weights ^ 2)
  popmat[reproduction_mask] <- 0.5 * prod(early_surv) * reprod

  # define basic biomass-based density dependence
  biomass_dd <- function(k, stages) {
    function(x, n) {
      sum_n <- sum(n[min(stages$dims):length(n)])
      x * ifelse(
        sum_n > (stages$scaling * k),
        (stages$scaling * k) / sum_n,
        1
      )
    }
  }
  dd_stages <- list(
    list(dims = 1, scaling = 4),
    list(dims = 2:9, scaling = 1),
    list(dims = 10:14, scaling = 1),
    list(dims = 15:24, scaling = 1),
    list(dims = 25:50, scaling = 1)
  )
  biomass_dd_list <- lapply(dd_stages, biomass_dd, k = k)
  biomass_mask_list <- lapply(
    dd_stages, transition, mat = popmat
  )
  dens_depend <- density_dependence(biomass_mask_list, biomass_dd_list)

  # covariate effects based on standardised discharge metrics
  #   - associations estimated and described in Tonkin et al. 2020 (STOTEN)
  recruitment_effects <- function(mat, x, system, ...) {

    # define system-specific coefficients
    coefs <- list(
      murray = c(-93.1, 25.2, 168.7, -34.0, 60.9, -14.9),
      ovens = c(26.5, -10.7, -69.7, 12.2, 2.7, 2.4),
      goulburn = c(92.0, 0.1, -28.7, -36.0, 4.5, -25.3),
      king = c(-21.9, -6.9, -18.5, -30.9, 14.3, -3.6),
      broken = c(-0.7, 7.1, -1.2, -7.3, 0.6, 16.1)
    )

    # switch for currently used covariates
    if (!system %in% c("murray", "goulburn", "campaspe")) {
      stop('system must be set as "murray", "goulburn", or "campaspe" ',
           'in a call to get_args("murray_cod"), with the returned ',
           'arguments passed to simulate.',
           call. = FALSE)
    }
    if (system == "goulburn" | system == "murray")
      coefs <- coefs$murray
    if (system == "campaspe")
      coefs <- coefs$broken

    # pull out relevant system
    names(coefs) <- c(
      "flow_var",
      "spring_flow",
      "max_ante",
      "summer_flow",
      "winter_flow",
      "spawning_temp"
    )

    # calculate scaling factor by year
    metrics <- c(
      x$proportional_flow_variability,
      x$proportional_spring_flow,
      x$proportional_max_antecedent,
      x$proportional_summer_flow,
      x$proportional_winter_flow,
      x$spawning_temperature
    )
    effect <- metrics * (coefs / 100)

    # fill NAs
    effect[is.na(effect)] <- 0

    # summarise over all metrics
    effect <- sum(effect)

    # calculate change in fecundity
    mat <- mat + effect * (mat / 5)
    mat[mat < 0] <- 0

    # return
    mat

  }

  # low-flow effect that reduces survival when flow falls below 5 ML/day due
  #    to increased risk of blackwater/hypoxic events
  survival_effects <- function(mat, x, ...) {
    scaling <- 0.8 + (0.2 / (1 + exp(-2 * (x$minimum_daily_flow - 5))))
    scaling * mat
  }
  cov_effects <- list(
    recruitment_effects,
    survival_effects
  )
  cov_masks <- list(
    reproduction(popmat, dims = reproductive),
    transition(popmat, dims = reproductive)
  )
  covars <- covariates(
    masks = cov_masks,
    funs = cov_effects
  )

  # define environmental stochasticity
  envstoch <- environmental_stochasticity(
    masks = list(
      transition(popmat),
      reproduction(popmat, dims = reproductive)
    ),
    funs = list(survival_gen, reproduction_gen)
  )

  # define angling effects
  go_fishing <- function(
    n, p_capture, slot, ...
  ) {

    # only needed if p_capture > 0
    if (p_capture > 0) {

      # flatten ages
      age_vector <- rep(seq_along(n), times = n)

      #   with n_available individuals within slot
      length_est <- estimate_length(age_vector)

      # which individuals are within the slot?
      catchable <- length_est >= slot[1] & length_est <= slot[2]

      # binary switch to determine which ones get caught
      caught <- rbinom(n = sum(catchable), size = 1, prob = p_capture)

      # which ages were caught?
      caught <- age_vector[caught == 1]

      # remove caught individuals from relevant age classes,
      #   checking to make sure at least one individual was
      #   caught
      if (length(caught) > 0) {
        to_remove <- table(caught)
        idx <- as.numeric(names(to_remove))
        n[idx] <- n[idx] - to_remove
      }

    }

    # return
    n

  }

  # use density_dependence_n to include stocking,
  #   translocations, or angling
  dd_n_masks <- list(
    all_classes(popmat, dim = 1),
    all_classes(popmat, dim = 2),
    all_classes(popmat, dim = reproductive),
    all_classes(popmat)
  )
  dd_n_fns <- list(
    function(pop, n_yoy, add_yoy, ...)
      add_remove(pop = pop, n = n_yoy, add = add_yoy),
    function(pop, n_twoplus, add_twoplus, ...)
      add_remove(pop = pop, n = n_twoplus, add = add_twoplus),
    function(pop, n_adult, add_adult, ...)
      add_remove(pop = pop, n = n_adult, add = add_adult),
    go_fishing
  )
  dens_depend_n <- density_dependence_n(
    masks = dd_n_masks,
    funs = dd_n_fns
  )

  # nolint start
  # print a message to ensure args are included
  message('Arguments are required to simulate Murray cod dynamics.\n',
          'Default arguments can be accessed with get_args("murray_cod").')
  # nolint end

  # return template
  list(
    matrix = popmat,
    covariates = covars,
    environmental_stochasticity = envstoch,
    demographic_stochasticity = NULL,
    density_dependence = dens_depend,
    density_dependence_n = dens_depend_n
  )

}

#' @rdname murray_cod
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
#' @param p_capture probability of capture by recreational anglers,
#'   defaults to 0, in which case recreational fishing does not
#'   occur
#' @param slot length slot within which Murray cod can legally
#'   be removed by recreational anglers. Defined in millimetres,
#'   defaults to 550-750 mm
#' @param system one of \code{"murray"}, \code{"goulburn"}, or
#'   \code{"campaspe"}, defining flow associations based on those
#'   observed in the Murray River, Goulburn River (currently identical
#'   to observations in the Murray River), or the Campaspe River
#'   (currently based on observations in the Broken River). Defaults
#'   to \code{"murray"}
#'
#' @details The Murray cod population template requires
#'   several additional arguments and allows several optional
#'   arguments. Arguments can be defined with a call to
#'   \code{get_args("murray_cod")} and are described
#'   individually above.
#'
args_murray_cod <- function(
  n = c(0, 0, 0),
  ntime = 50,
  start = c(1, 1, 1),
  end = c(1, 1, 1),
  add = TRUE,
  p_capture = 0.0,
  slot = c(550, 750),
  system = "murray"
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

  # define mean early-life survival
  early_surv <- c(0.5, 0.012, 0.38, 0.31)

  # force evaluation of early-life survival
  force(early_surv)

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
      early_surv,  # early life survival
      mat[transition(mat)]    # from population matrix in current time step
    )
    survival_sd <- c(
      0.2 * survival_mean[1:7],
      0.15 * survival_mean[8:10],
      0.1 * survival_mean[11:length(survival_mean)]
    )

    # convert unit interval to real line equivalents
    out <- unit_to_real(
      unit_mean = survival_mean,
      unit_sd = survival_sd
    )

    # separate early life from other estimates
    idx <- seq_len(nrow(out)) > 4

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

    # set contributing as random uniform on 0.75-1.0 by default
    # set recruit_failure at 0 by default
    # add function to pre-transform unit to real and back
    environmental_stochasticity = list(
      transform_survival
    ),

    # to include additions or removals of individuals
    density_dependence_n = list(
      define_removals(
        start = start, end = end, n = n, add = add
      ),
      p_capture = p_capture,
      slot = slot
    ),

    # to set covariates relative to a specific system
    covariates = list(
      system = system
    )

  )

}
