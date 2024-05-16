#' @name common_carp
#' @title Parameterised common carp population model
#' @description Use a pre-defined population dynamics object to
#'   simulate population dynamics of common carp (*Cyrpinus carpio*).
#'   Model parameters are based on existing data sets and
#'   published studies on common carp.
NULL

#' @rdname common_carp
#'
#' @importFrom stats rnorm qlogis plogis
#' @import aae.pop
#'
#' @export
#'
#' @param k carrying capacity, defaults to 50000 but a rule-of-thumb
#'   for setting this value is 2 adult females per linear metre of
#'   a large river (e.g. the Murray River)
#' @param reproductive integer vector defining the adult stages in the
#'   common carp population. Defaults to \code{3:28}
#' @param system one of \code{"main_channel"}, \code{"river_wetland"},
#'   \code{"ephemeral_wetland"}, \code{"permanent_wetland"},
#'   or \code{"floodplain"}, defining the early life survival.
#'   Defaults to "main_channel" and ignores case
#' @param n an integer, vector of integers, or list specifying the
#'   number of young-of-year, 1+, or older fish to add or
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
#'
#' @details The \code{common_carp} population model is a mixed
#'   age- and stage-structured model with 28 classes and includes negative
#'   density dependence, environmental and demographic
#'   stochasticity, and optional associations with hydrological conditions
#'   in rivers and removal efforts.
#'
#'   The common carp population template requires
#'   several additional arguments and allows several optional
#'   arguments. Arguments are described
#'   individually above.
#'
#' @examples
#' # define a basic model for common carp
#' carp <- common_carp()
#'
#' # simulate from this model
#' sims <- simulate(carp, nsim = 100)
#'
#' # plot the simulated values
#' plot(sims)
common_carp <- function(
  k = 50000,
  reproductive = 3:28,
  system = "main_channel",
  n = c(0, 0, 0),
  ntime = 50,
  start = c(1, 1, 1),
  end = c(1, 1, 1),
  add = TRUE
) {
  get_template(
    sp = "common_carp",
    k = k,
    system = system,
    reproductive = reproductive,
    n = n,
    ntime = ntime,
    start = start,
    end = end,
    add = add
  )
}

# internal function: define species defaults
template_common_carp <- function(
  k = 50000,
  reproductive = 3:28,
  system = "main_channel",
  n = c(0, 0, 0),
  ntime = 50,
  start = c(1, 1, 1),
  end = c(1, 1, 1),
  add = TRUE
) {

  # how many stages are we going to work with?
  nstage <- 28

  # check system
  system <- tolower(system)
  if (
    !system %in% c("main_channel",
                   "river_wetland",
                   "ephemeral_wetland",
                   "permanent_wetland",
                   "floodplain")
  ) {
    stop("system must be one of main_channel, ",
         "river_wetland, ephemeral_wetland, ",
         "permanent_wetland, or floodplain",
         "in a call to common_carp",
         call. = FALSE)
  }

  # set system-specific max lengths
  early_surv_options <- c(
    "main_channel" = 0.0245 * 0.0524 * 0.0689 * 0.11,
    "river_wetland" = 0.1207 * 0.10 * 0.2141 * 0.155,
    "ephemeral_wetland" = 0.0796 * 0.0570 * 0.1683 * 0.0796,
    "permanent_wetland" = 0.0645 * 0.0654 * 0.1484 * 0.2112,
    "floodplain" = 0.109 * 0.0815 * 0.2031 * 0.2139
  )
  early_surv <- early_surv_options[system]

  # force evaluation of max length
  force(early_surv)

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
  mean_surv <- c(0.20, 0.54, 0.67, 0.74, 0.78, 0.80, 0.82, 0.84, 0.85,
                 0.86, 0.86, 0.86, 0.87, 0.87, 0.87, 0.87, 0.86, 0.86,
                 0.85, 0.85, 0.83, 0.82, 0.80, 0.77, 0.72, 0.64, 0.48)

  # define base matrix
  popmat <- matrix(0, nrow = nstage, ncol = nstage)

  # add survival
  popmat[transition(popmat)] <- mean_surv

  # add reproduction
  reproduction_mask <- reproduction(popmat, dims = reproductive)
  reprod <- exp(0.0051 * log(reproductive) + 13.04)
  popmat[reproduction_mask] <- 0.5 * prod(early_surv) * reprod

  # bottom-up effects of density on early survival
  #   through competition for resources
  bh <- function(x, pop, theta = 0.2, ...) {
    x / (1 + theta * x * sum(pop[reproductive]) / k)
  }

  # and collate masks and functions in a single object
  dens_depend <- density_dependence(
    masks = reproduction(popmat, dims = reproductive),
    funs = bh
  )

  # covariate effects based on stable flows and floodplain
  #   connectivity
  recruit_effects <- function(mat, x, variability_param = -0.05, ...) {

      # floodplain early survival is approx. 40 times higher
      #   than in main channel when flows are very high (connected
      #   to floodplain); steep drop-off as flows reduce within bank
      early_surv <- 40 / (1 + exp(-8 * (x$floodplain_access - 0.5)))
      # early_surv <- x$floodplain_access * 40

      # and they tend to prefer stable conditions with higher summer flows,
      # so add a negative effect of flow variability on recruitment
      mat <- mat * exp(-variability_param * x$flow_variability)

      # no summer flows effect for now but TBD

      # return
      early_surv * mat

  }

  # define covars objet
  covars <- covariates(
    masks = reproduction(popmat, dims = reproductive),
    funs = recruit_effects
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
  #   removals
  dd_n_masks <- list(
    all_classes(popmat, dims = 1),
    all_classes(popmat, dims = 2),
    all_classes(popmat, dims = reproductive)
  )
  dd_n_fns <- list(
    function(pop, n_yoy, add_yoy, ...) {
      add_remove(pop = pop, n = n_yoy, add = add_yoy)
    },
    function(pop, n_twoplus, add_twoplus, ...) {
      add_remove(pop = pop, n = n_twoplus, add = add_twoplus)
    },
    function(pop, n_adult, add_adult, ...) {
      add_remove(pop = pop, n = n_adult, add = add_adult)
    }
  )
  dens_depend_n <- density_dependence_n(
    masks = dd_n_masks,
    funs = dd_n_fns
  )

  # collect arguments for species if required
  arguments <- get_args(
    "common_carp",
    n = n,
    ntime = ntime,
    start = start,
    end = end,
    add = add
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

#' @rdname common_carp
#'
#' @importFrom stats pnorm rnorm runif
#'
#' @export
#'
args_common_carp <- function(
  n = c(0, 0, 0),
  ntime = 50,
  start = c(1, 1, 1),
  end = c(1, 1, 1),
  add = TRUE
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
      if (nrow(add) != 3 || ncol(add) != ntime) {
        stop("if add is a matrix, it must have three rows ",
             "and ntime columns (ntime = ", ntime, ")",
             call. = FALSE)
      }
    }

    # force evaluation of n and add so they don't get lost below
    force(n)
    force(add)

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
      define_removals(start = start, end = end, n = n, add = add)
    )

  )

}
