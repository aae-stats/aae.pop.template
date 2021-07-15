#' @name australian_bass
#' @title Parameterised Australian bass population model
#' @description Use a pre-defined population dynamics object to
#'   simulate population dynamics of Australian bass
#'   (*Percalates novemaculeata*). Model parameters are based
#'   on existing data sets and published studies on Australian bass.
NULL

#' @rdname australian_bass
#'
#' @importFrom stats rnorm
#'
#' @export
#'
#' @param k carrying capacity
#' @param n an integer, vector of integers, or list specifying the
#'   number of young-of-year, 2+, 3+, and adult fish to add or
#'   remove in any given year. If a single integer is provided,
#'   all age classes are assumed to have the same number of
#'   additions or removals. If a vector is provided, it must
#'   have one value for each age class. If a list is provided,
#'   it must have one element for each age class, with a
#'   value for each time step. Addition versus removal is
#'   controlled with \code{add}. Defaults to
#'   \code{c(0, 0, 0, 0))}, which will specify no additions
#'   or removals
#' @param ntime number of time steps used in population
#'   simulation. Defaults to \code{50} and is required
#'   to ensure simulated additions or removals are defined
#'   for every simulated time step
#' @param start time step at which additions or removals start
#'   if \code{n} is an integer or vector. Defaults to
#'   \code{c(1, 1, 1, 1)}
#' @param end time step at which additions or removals finish
#'   if \code{n} is an integer or vector. Defaults to
#'   \code{c(1, 1, 1, 1)}
#' @param add logical indicating whether individuals are
#'   added or removed from the population. Defaults to
#'   \code{TRUE}, which defines additions (e.g., stocking). Can
#'   be specified as a single value, in which case all stages
#'   are set equal, as four values, in which case stages can
#'   differ, or as a matrix with four rows and \code{ntime}
#'   columns, in which case the effects of \code{add} can
#'   change through time
#' @param p_capture probability of capture by recreational anglers,
#'   defaults to 0, in which case recreational fishing does not
#'   occur
#'
#' @details The \code{australian_bass} population model is a
#'   stage-structured model with 50 classes and includes negative
#'   density dependence, environmental and demographic
#'   stochasticity, and optional associations with environmental
#'   covariates, including recreational angling.
#'
#'   The Australian bass population template requires
#'   several additional arguments and allows several optional
#'   arguments, described individually above.
#'
#' @examples
#' # define a basic model for Australian bass
#' p <- australian_bass()
#'
#' # simulate from this model
#' sims <- simulate(p, nsim = 100)
#'
#' # plot the simulated values
#' plot(sims)
australian_bass <- function(
  k = 30000,
  n = 0,
  ntime = 50,
  start = 1,
  end = 1,
  add = TRUE,
  p_capture = 0.0
) {
  get_template(
    sp = "australian_bass",
    k = k,
    n = n,
    ntime = ntime,
    start = start,
    end = end,
    add = add,
    p_capture = p_capture
  )
}

# internal function: define species defaults
template_australian_bass <- function(
  k = 30000,
  n = 0,
  ntime = 50,
  start = 1,
  end = 1,
  add = TRUE,
  p_capture = 0.0
) {

  # how many stages are we going to work with?
  nstage <- 50

  # and which are reproductive?
  reproductive <- 4:50

  # define  a survival function
  survival_gen <- function(mat, ...) {
    plogis(
      rnorm(
        length(mat),
        mean = qlogis(mat),
        sd = 0.5 * abs(qlogis(mat)))
    )
  }

  # define a reproduction function
  reproduction_gen <- function(
    mat,
    proportion_female = 0.5,
    proportion_breeding = 0.85,
    ...
  ) {

    # simulate realised fecundity in a given year
    n_offspring <- rpois(1, lambda = mat)

    # estimate and return reproduction estimates
    n_offspring * proportion_female * proportion_breeding

  }

  # define base matrix
  popmat <- matrix(0, nrow = nstage, ncol = nstage)

  # survival of egg, larvae and 0+ fish
  early_survival <- c(0.5, 0.0061, 0.1226)

  # add survival
  age_frequency <- (1020.23 * (seq_len(nstage) ^ -1.16)) - 13.5
  popmat[transition(popmat)] <-
    age_frequency[-1] / age_frequency[-length(age_frequency)]

  # tweak new stages that give wacky survival estimates
  popmat[transition(popmat)][41:49] <- (1 - seq(0.1, 0.9, length = 9)) * popmat[transition(popmat)][40]

  # add reproduction
  fecundity_by_age <- exp(
    -11.83 +
      4.21 * log(384.69 * (1 - exp(-0.19 * (reproductive + 2.27))))
  )
  popmat[reproduction(popmat, dims = reproductive)] <-
    prod(early_survival) * fecundity_by_age

  # define custom density dependence function
  theta_ricker <- function(x, n, theta = 4, r = 0.4) {
    x * exp(r * (1 - (sum(n[reproductive]) / k) ^ theta)) / exp(r)
  }
  dens_depend <- density_dependence(
    reproduction(popmat),
    theta_ricker
  )

  # covariate effects based on probability of
  #   spawning and larval survival
  #     due to flow conditions prior to and following spawning,
  #   Values of x calculated and returned independently of template
  #     at present
  fecundity_effects <- function(
    mat,
    x,
    alpha = 0.5,
    beta = 0,
    gamma = -0.01,
    offset = 14,
    delta = 3,
    ...
  ) {

    # include effects of water temperature
    scaling <- alpha +
      delta * (beta * (x$water_temperature - offset) +
                 gamma * (x$water_temperature - offset) ^ 2)

    # catch negative scaling values
    scaling[scaling < 0] <- 0

    # return rescaled fecundity with an additional term
    #   for discharge effects (calculated externally)
    mat * x$discharge * scaling

  }

  # combine into a covariates object
  covars <- covariates(
    masks = reproduction(popmat, dims = reproductive),
    funs = fecundity_effects
  )

  # define environmental stochasticity
  envstoch <- environmental_stochasticity(
    masks = list(
      transition(popmat),
      reproduction(popmat, dims = reproductive)
    ),
    funs = list(
      survival_gen,
      reproduction_gen
    )
  )

  # define angling effects
  go_fishing <- function(
    n, p_capture, ...
  ) {

    # only needed if p_capture > 0
    if (p_capture > 0) {

      # flatten ages
      age_vector <- rep(seq_along(n), times = n)

      # which individuals are within the slot?
      catchable <- age_vector >= 4

      # binary switch to determine which ones get caught
      caught <- rbinom(
        n = sum(catchable), size = 1, prob = p_capture
      )

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

  # use density_dependence_n to include stocking, translocations,
  #   and fishing
  dd_n_masks <- list(
    all_classes(popmat),
    all_classes(popmat, dim = 2),
    all_classes(popmat, dim = 3),
    all_classes(popmat, dim = reproductive),
    all_classes(popmat)
  )
  dd_n_fns <- list(
    function(pop, n_yoy, add_yoy, ...) {
      n_yoy <- floor(theta_ricker(n_yoy, pop))
      pop[1] <- add_remove(pop = pop[1], n = n_yoy, add = add_yoy)
      pop
    },
    function(pop, n_twoplus, add_twoplus, ...)
      add_remove(pop = pop, n = n_twoplus, add = add_twoplus),
    function(pop, n_threeplus, add_threeplus, ...)
      add_remove(pop = pop, n = n_threeplus, add = add_threeplus),
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
    "australian_bass",
    n = n,
    ntime = ntime,
    start = start,
    end = end,
    add = add,
    p_capture = p_capture
  )

  # return template
  list(
    dynamics = list(
      matrix = popmat,
      covariates = covars,
      environmental_stochasticity = envstoch,
      density_dependence = dens_depend,
      density_dependence_n = dens_depend_n
    ),
    arguments = arguments
  )

}

#' @rdname australian_bass
#'
#' @importFrom stats pnorm rnorm runif
#'
#' @export
#'
args_australian_bass <- function(
  n = c(0, 0, 0, 0),
  ntime = 50,
  start = c(1, 1, 1, 1),
  end = c(1, 1, 1, 1),
  add = TRUE,
  p_capture = 0.0
) {

  # expand n, start, end, add if required
  if (length(n) == 1)
    n <- rep(n, 4)
  if (length(start) == 1)
    start <- rep(start, 4)
  if (length(end) == 1)
    end <- rep(end, 4)
  if (length(add) == 1)
    add <- rep(add, 4)

  # check for other lengths of n, start, or end
  if (any(c(length(n), length(start), length(end)) != 4)) {
    stop("n, start, and end must be vectors with 1 or 4 elements",
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
      add <- matrix(rep(add, ntime), nrow = 4)
    } else {
      if (nrow(add) != 4 | ncol(add) != ntime) {
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
        n_threeplus = n[[3]][iter],
        n_adult = n[[4]][iter],
        add_yoy = add[1, iter],
        add_twoplus = add[2, iter],
        add_threeplus = add[3, iter],
        add_adult = add[4, iter]
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
      p_capture = p_capture
    )

  )

}
