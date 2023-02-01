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
#' @importFrom stats rnorm rbinom rbeta qlogis plogis
#' @import aae.pop
#'
#' @export
#'
#' @param k carrying capacity
#' @param reproductive integer vector defining the adult stages in the
#'   Murray cod population. Defaults to \code{5:50}
#' @param system one of \code{"murray_river"}, \code{"goulburn_river"},
#'   \code{"campaspe_river"}, \code{"ovens_river"}, \code{"broken_river"}
#'   or \code{"broken_creek"}, defining whether the maximum size of
#'   individuals matches observations in the Murray River (1300 mm),
#'   large tributaries (e.g., Goulburn River, Campaspe River; 1100 mm),
#'   or small tributaries (e.g., XYZ; 900 mm). Defaults to "murray" and
#'   ignores case
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
#' @details The \code{murray_cod} population model is a mixed
#'   age- and stage-structured model with 25 classes and includes negative
#'   density dependence, environmental and demographic
#'   stochasticity, and optional associations with hydrological conditions
#'   in rivers and stocking or angling effects.
#'
#'   The Murray cod population template requires
#'   several additional arguments and allows several optional
#'   arguments. Arguments are described
#'   individually above.
#'
#' @examples
#' # define a basic model for Murray cod
#' mc <- murray_cod()
#'
#' # simulate from this model
#' sims <- simulate(mc, nsim = 100)
#'
#' # plot the simulated values
#' plot(sims)
murray_cod <- function(
  k = 20000,
  system = "murray_river",
  reproductive = 5:50,
  n = c(0, 0, 0),
  ntime = 50,
  start = c(1, 1, 1),
  end = c(1, 1, 1),
  add = TRUE,
  p_capture = 0.0,
  slot = c(550, 750)
) {
  get_template(
    sp = "murray_cod",
    k = k,
    system = system,
    reproductive = reproductive,
    n = n,
    ntime = ntime,
    start = start,
    end = end,
    add = add,
    p_capture = p_capture,
    slot = slot
  )
}

# internal function: define species defaults
template_murray_cod <- function(
  k = 20000,
  system = "murray_river",
  reproductive = 5:50,
  n = c(0, 0, 0),
  ntime = 50,
  start = c(1, 1, 1),
  end = c(1, 1, 1),
  add = TRUE,
  p_capture = 0.0,
  slot = c(550, 750)
) {

  # how many stages are we going to work with?
  nstage <- 50

  # check system
  # nolint start
  system <- tolower(system)
  if (
    !system %in% c("murray_river", "goulburn_river", "campaspe_river",
                   "ovens_river", "broken_river", "broken_creek")
  ) {
    stop('system must be one of "murray_river"',
         '"goulburn_river", "campaspe_river", "ovens_river"',
         '"broken_river", or "broken_creek" ',
         'in a call to murray_cod',
         call. = FALSE)
  }

  # set as small or large trib
  sys_set <- switch(
    system,
    "murray_river" = "murray",
    "goulburn_river" = "largetrib",
    "campaspe_river" = "largetrib",
    "ovens_river" = "largetrib",
    "broken_river" = "largetrib",
    "broken_creek" = "largetrib",
    "murray"
  )

  # set system-specific max lengths
  max_length_options <- c(
    "murray" = 1300,
    "largetrib" = 1100,
    "smalltrib" = 900
  )
  max_length <- max_length_options[sys_set]

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
  survival_gen <- function(mat, ...) {
    plogis(rnorm(length(mat), mean = qlogis(mat), sd = 0.5 * abs(qlogis(mat))))
  }

  # define a reproduction function, being careful with argument
  #   names to avoid conflicts with any arguments in
  #   survival_gen, which would then require multiple different
  #   arguments with the same name in simulate
  reproduction_gen <- function(
    mat,
    ...
  ) {

    # simulate average number of recruits based on vital rates
    reprod <- rnorm(length(mat), mean = mat, sd = 0.25 * mat)

    # catch any negative values
    reprod[reprod < 0] <- 0

    # multiply by 0.5 to account for a 50:50 sex ratio
    #   (don't need product of early survival because
    #    it's already accounted for in pop mat
    0.5 * reprod

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
    function(x, n, kdyn = NULL) {
      if (!is.null(kdyn))
        k <- kdyn
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
  recruitment_effects <- function(mat, x, system, coefs = NULL, threshold = 0, ...) {

    # define system-specific coefficients
    if (is.null(coefs)) {
      coefs <- list(
        murray_river = c(-93.1, 25.2, 168.7, -34.0, 60.9, -14.9),
        ovens_river = c(26.5, -10.7, -69.7, 12.2, 2.7, 2.4),
        goulburn_river = c(92.0, 0.1, -28.7, -36.0, 4.5, -25.3),
        king_river = c(-21.9, -6.9, -18.5, -30.9, 14.3, -3.6),
        broken_river = c(-0.7, 7.1, -1.2, -7.3, 0.6, 16.1),
        broken_creek = c(-0.7, 7.1, -1.2, -7.3, 0.6, 16.1)
      )

      # need to know the system if coefs are not provided
      if (!system %in% names(coefs)) {
        stop(
          "system must be specified if coefs are not provided, ",
          "and must be one of murray_river, ovens_river, ",
          "goulburn_river, king_river, broken_river, or broken_creek",
          call. = FALSE
        )
      }
      coefs <- coefs[[system]]

    } else {

      # need six coefficients
      if (length(coefs) != 6) {
        stop(
          "coefs must include six values, one for each of ",
          "flow variability, spring flows, max. antecedent flows, ",
          "summer flows, winter flows, and spawning temperature",
          .call = FALSE
        )
      }

    }

    # nolint end

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
      x$spawning_flow_variability,
      x$proportional_spring_flow,
      x$proportional_max_antecedent,
      x$proportional_summer_flow,
      x$proportional_winter_flow,
      x$spawning_temperature
    )
    effect <- 1 + (metrics * (coefs / 100))

    # fill NAs
    effect[is.na(effect)] <- 0

    # and set negative values to zero
    effect[effect < 0] <- 0

    # summarise over all metrics
    #   - using product, which will limit to bottlenecks
    #     (e.g., 100% mortality in any stage = 0 survival overall)
    effect <- prod(effect)

    # threshold effect to avoid really small values (recruitment always
    #   occurs for MC)
    effect[effect < threshold] <- threshold

    # divide by early life survival to avoid double counting
    #  (is included in reproduction_gen function and needs to be
    #   there so non-covar models still work)
    effect <- effect

    # calculate change in fecundity
    mat <- mat * effect

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

  # collect arguments for species if required
  arguments <- get_args(
    "murray_cod",
    n = n,
    ntime = ntime,
    start = start,
    end = end,
    add = add,
    p_capture = p_capture,
    slot = slot
  )

  # return
  list(
    dynamics = list(
      matrix = popmat,
      covariates = covars,
      environmental_stochasticity = envstoch,
      demographic_stochasticity = NULL,
      density_dependence = dens_depend,
      density_dependence_n = dens_depend_n
    ),
    arguments = arguments
  )

}

#' @rdname murray_cod
#'
#' @importFrom stats pnorm rnorm runif
#'
#' @export
#'
args_murray_cod <- function(
  n = c(0, 0, 0),
  ntime = 50,
  start = c(1, 1, 1),
  end = c(1, 1, 1),
  add = TRUE,
  p_capture = 0.0,
  slot = c(550, 750)
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

  # return named list of args
  list(

    # to include additions or removals of individuals
    density_dependence_n = list(
      define_removals(
        start = start, end = end, n = n, add = add
      ),
      p_capture = p_capture,
      slot = slot
    )

  )

}
