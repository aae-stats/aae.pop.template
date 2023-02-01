#' @name barred_galaxias
#' @title Parameterised barred galaxias population model
#' @description Use a pre-defined population dynamics object to
#'   simulate population dynamics of barred galaxias
#'   (*Galaxias fuscus*). Model parameters are based
#'   on existing data sets and published studies on barred galaxias.
NULL

#' @rdname barred_galaxias
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
#' @details The \code{barred_galaxias} population model is a
#'   stage-structured model with 4 classes and includes negative
#'   density dependence, environmental and demographic
#'   stochasticity, and optional associations with environmental
#'   covariates
#'
#'   The barred galaxias population template requires
#'   several additional arguments and allows several optional
#'   arguments, described individually above.
#'
#' @examples
#' # define a basic model for barred galaxias
#' p <- barred_galaxias()
#'
#' # simulate from this model
#' sims <- simulate(p, nsim = 100)
#'
#' # plot the simulated values
#' plot(sims)
barred_galaxias <- function(
  k = 1000,
  n = 0,
  ntime = 50,
  start = 1,
  end = 1,
  add = TRUE
) {
  get_template(
    sp = "barred_galaxias",
    k = k,
    n = n,
    ntime = ntime,
    start = start,
    end = end,
    add = add
  )
}

# internal function: define species defaults
template_barred_galaxias <- function(
  k = 1000,
  n = 0,
  ntime = 50,
  start = 1,
  end = 1,
  add = TRUE
) {

  # how many stages are we going to work with?
  nstage <- 4

  # and which are reproductive?
  reproductive <- 3:4

  # define  a survival function
  survival_gen <- function(mat, ...) {

    # make sure values are not exactly 0 or 1
    eps <- 1e-4
    mat[mat == 1] <- 1 - eps
    mat[mat == 0] <- eps

    # return
    plogis(
      rnorm(length(mat), mean = qlogis(mat), sd = 0.5 * abs(qlogis(mat)))
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

  # define early survival
  early_survival <- 0.1

  # add survival
  popmat[transition(popmat)] <- c(0.2, 0.1, 0.1)
  popmat[survival(popmat)] <- c(0.4, 0.5, 0.6, 0.65)

  # add reproduction
  popmat[reproduction(popmat, dims = reproductive)] <-
    prod(early_survival) * c(200, 500)

  # define custom density dependence function
  theta_ricker <- function(x, n, theta = 4, r = 0.4) {
    x * exp(r * (1 - (sum(n[reproductive]) / k) ^ theta)) / exp(r)
  }
  dens_depend <- density_dependence(
    reproduction(popmat, dims = reproductive),
    theta_ricker
  )

  # define effects of covariates on juvenile survival (stages 1 and 2)
  juvenile_survival_effects <- function(mat, x, ...) {

    # account for trout
    mat <- ifelse(x$trout, 0.7 * mat, mat)

    # account for bushfire
    if (x$bushfire) {
      scale_factor <- 1 / (2 - x$riparian)
      mat <- scale_factor * mat
    }

    # account for ctf
    if (x$ctf) {
      scale_factor <- 1 / (5 - 3 * x$riparian)
      mat <- scale_factor * mat
    }

    # account for riparian
    mat <- x$riparian * mat
    mat <- ifelse(mat > 1, 1, mat)

    # return
    mat

  }

  # define effects of covariates on adult survival (stages 3 and 4)
  adult_survival_effects <- function(mat, x, ...) {

    # account for trout
    mat <- ifelse(x$trout, 0, mat)

    # account for bushfire
    if (x$bushfire) {
      scale_factor <- 3 / (4 - x$riparian)
      mat <- scale_factor * mat
    }

    # account for ctf
    if (x$ctf) {
      scale_factor <- 3 / (5 - 1.25 * x$riparian)
      mat <- scale_factor * mat
    }

    # account for riparian
    mat <- x$riparian * mat
    mat <- ifelse(mat > 1, 1, mat)

    # return
    mat

  }

  # define effects of covariates on reproduction
  reproduction_effects <- function(mat, x, ...) {

    # account for bushfire
    if (x$bushfire) {
      scale_factor <- 1 / (2 - x$riparian)
      mat <- scale_factor * mat
    }

    # account for ctf
    if (x$ctf) {
      scale_factor <- 1 / (5 - 3 * x$riparian)
      mat <- scale_factor * mat
    }

    # account for riparian
    mat <- x$riparian * mat

    # return
    mat

  }

  # collate masks and functions into lists
  covars <- covariates(
    masks = list(
      combine(survival, transition)(popmat, dims = 1:2),
      combine(survival, transition)(popmat, dims = reproductive),
      reproduction(popmat, dims = reproductive)
    ),
    funs = list(
      juvenile_survival_effects,
      adult_survival_effects,
      reproduction_effects
    )
  )

  # define environmental stochasticity
  envstoch <- environmental_stochasticity(
    masks = list(
      combine(survival, transition)(popmat),
      reproduction(popmat, dims = reproductive)
    ),
    funs = list(
      survival_gen,
      reproduction_gen
    )
  )

  # define demographic stochasticity
  demostoch <- demographic_stochasticity(
    masks = list(
      all_classes(popmat, dims = 1)
    ),
    funs = list(
      function(x) rpois(length(x), lambda = x)
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

  # collect arguments for species if required
  arguments <- get_args(
    "barred_galaxias",
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
      demographic_stochasticity = demostoch,
      environmental_stochasticity = envstoch,
      density_dependence = dens_depend,
      density_dependence_n = dens_depend_n
    ),
    arguments = arguments
  )

}

#' @rdname barred_galaxias
#'
#' @importFrom stats pnorm rnorm runif
#'
#' @export
#'
args_barred_galaxias <- function(
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
      )
    )

  )

}
