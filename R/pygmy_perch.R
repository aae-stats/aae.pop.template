#' @name pygmy_perch
#' @title Parameterised pygmy perch population model
#' @description Use a pre-defined population dynamics object to
#'   simulate population dynamics of pygmy perch
#'   (*Nannoperca australis*). Model parameters are based
#'   on existing data sets and published studies on pygmy perch.
NULL

#' @rdname pygmy_perch
#'
#' @importFrom stats rnorm rpois
#' @import aae.pop
#'
#' @export
#'
#' @param k carrying capacity
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
#'
#' @details The \code{pygmy_perch} population model is a
#'   stage-structured model with 4 classes and includes negative
#'   density dependence, environmental and demographic
#'   stochasticity, and optional associations with environmental
#'   covariates
#'
#'   The pygmy perch population template requires
#'   several additional arguments and allows several optional
#'   arguments, described individually above.
#'
#' @examples
#' # define a basic model for pygmy perch
#' p <- pygmy_perch()
#'
#' # simulate from this model
#' sims <- simulate(p, nsim = 100)
#'
#' # plot the simulated values
#' plot(sims)
pygmy_perch <- function(
  k = 2000,
  n = 0,
  ntime = 50,
  start = 1,
  end = 1,
  add = TRUE
) {
  get_template(
    sp = "pygmy_perch",
    k = k,
    n = n,
    ntime = ntime,
    start = start,
    end = end,
    add = add
  )
}

# internal function: define species defaults
template_pygmy_perch <- function(
  k = 2000,
  n = 0,
  ntime = 50,
  start = 1,
  end = 1,
  add = TRUE
) {

  # how many stages are we going to work with?
  nstage <- 4

  # and which are reproductive?
  reproductive <- 1:4

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
    proportion_breeding = 1.0,
    ...
  ) {

    # simulate realised fecundity in a given year
    n_offspring <- rpois(length(mat), lambda = mat)

    # estimate and return reproduction estimates
    n_offspring * proportion_female * proportion_breeding

  }

  # define base matrix
  popmat <- matrix(0, nrow = nstage, ncol = nstage)

  # surivval of egg, larvae and 0+ fish
  early_survival <- 0.0055

  # add survival
  popmat[transition(popmat)] <- c(0.5270, 0.4872, 0.0526)

  # add reproduction
  fecundity_by_age <- exp(
    0.05068 * exp(0.1933 * reproductive + 3.4597) + 3.2043
  )
  # includes YOY, so not possible to use reproduction(popmat)
  popmat[1, ] <- prod(early_survival) * fecundity_by_age

  # define custom density dependence function
  theta_ricker <- function(x, n, theta = 4, r = 0.4) {
    x * exp(r * (1 - (sum(n[reproductive]) / k) ^ theta)) / exp(r)
  }
  dens_depend <- density_dependence(
    reproduction(popmat),
    theta_ricker
  )

  # leave basic covariate effects for now but will be ignored
  #   in most cases
  fecundity_effects <- function(
    mat,
    x,
    ...
  ) {

    # return as is
    mat

  }
  survival_effects <- function(
    mat,
    x,
    ...
  ) {

    # account for predators
    if (x$predators)
      mat <- 0.9 * x$habitat * mat

    # account for persistent water availability
    if (x$dry)
      mat <- 0.2 * mat

    # return directly scaled value because x is calculated
    #   as a probability
    mat

  }

  # combine into a covariates object
  covars <- covariates(
    masks = list(
      reproduction(popmat, dims = reproductive),
      transition(popmat)
    ),
    funs = list(
      fecundity_effects,
      survival_effects
    )
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

  # use density_dependence_n to include stocking, translocations,
  #   and fishing
  dd_n_masks <- list(
    all_classes(popmat, dim = 1),
    all_classes(popmat, dim = 2),
    all_classes(popmat, dim = reproductive)
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
    "pygmy_perch",
    n = n,
    ntime = ntime,
    start = start,
    end = end,
    add = add
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

#' @rdname pygmy_perch
#'
#' @importFrom stats pnorm rnorm runif
#'
#' @export
#'
args_pygmy_perch <- function(
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
      )
    )

  )

}
