#' @name platypus
#' @title Parameterised platypus population model
#' @description Use a pre-defined population dynamics object to
#'   simulate population dynamics of platypus
#'   (*Ornithorhynchus anatinus*). Model parameters are based
#'   on existing data sets and published studies on platypus.
NULL

#' @rdname platypus
#'
#' @importFrom stats rnorm
#'
#' @export
#'
#' @param k carrying capacity. Previously estimated densities
#'   suggest 1-7 adult platypus per km in rivers, with larger values
#'   in systems that more closely resemble reference (or natural)
#'   conditions. Defaults to 400, which assumes a 100 km river
#'   reach at 4 individuals per km. The value of \code{k} is
#'   multiplied by \code{0.65} to give adult female carrying
#'   capacity, assuming a 65:35 female:male sex ratio in
#'   adult platypus
#' @param n an integer or list of integers specifying the
#'   number of adult platypus to add in any given year.
#'   If a single integer is provided,
#'   all years are assumed to have the same number of
#'   additions or removals. If a list is provided,
#'   it must have one value for each year. Defaults to
#'   \code{0}, which will specify no additions
#' @param ntime number of time steps used in population
#'   simulation. Defaults to \code{50} and is required
#'   to ensure simulated additions or removals are defined
#'   for every simulated time step
#' @param start time step at which additions or removals start
#'   if \code{n} is an integer. Defaults to
#'   \code{1}
#' @param end time step at which additions or removals finish
#'   if \code{n} is an integer. Defaults to
#'   \code{1}
#'
#' @details The \code{platypus} population model is a
#'   stage-structured model with 2 classes and includes negative
#'   density dependence, environmental and demographic
#'   stochasticity, and optional associations with hydrological
#'   conditions, instream and riparian habitat, and predators.
#'
#'   The platypus population template requires
#'   several additional arguments and allows several optional
#'   arguments, described individually above.
#'
#' @examples
#' # define a basic model for platypus
#' p <- platypus()
#'
#' # define required arguments
#' p_args <- get_args("platypus")
#'
#' # simulate from this model
#' sims <- simulate(p, nsim = 100, args = p_args)
#'
#' # plot the simulated values
#' plot(sims)
platypus <- function(k = 400, n = 0, ntime = 50, start = 1, end = 1) {
  get_template(
    sp = "platypus",
    k = k,
    n = n,
    ntime = ntime,
    start = start,
    end = end
  )
}

# internal function: define species defaults
template_platypus <- function(k = 400, n = 0, ntime = 50, start = 1, end = 1) {

  # how many stages are we going to work with?
  nstage <- 2

  # define  a survival function
  # define  a survival function
  survival_gen <- function(mat, ...) {
    plogis(rnorm(length(mat), mean = qlogis(mat), sd = 0.5 * abs(qlogis(mat))))
  }

  # define a reproduction function
  reproduction_gen <- function(
    mat,
    proportion_female = 0.5,
    proportion_breeding = 0.62,
    ...
  ) {

    # simulate realised fecundity in a given year
    n_offspring <- rpois(1, lambda = mat)

    # estimate and return reproduction estimates
    n_offspring * proportion_female * proportion_breeding

  }

  # define base matrix
  popmat <- matrix(0, nrow = nstage, ncol = nstage)

  # add survival
  popmat[2, 1] <- 0.29
  popmat[2, 2] <- 0.89

  # add reproduction
  popmat[1, 2] <- 0.5 * 0.62 * 2

  # define custom density dependence function
  k_female <- k * 0.65     # assumes 65:35 F:M sex ratio in adults
  theta_ricker <- function(x, n, theta = 4, r = 0.4) {
    x * exp(r * (1 - (n[2] / k_female) ^ theta)) / exp(r)
  }
  dens_depend <- density_dependence(
    reproduction(popmat),
    theta_ricker
  )

  # covariate effects based on CTF events, spring/summer/winter flows,
  #   and extremes (max flows and variability)
  survival_effects <- function(mat, x, ...) {

    # default is no change
    scale <- 1

    # add negative effects of low spring/summer/winter flows
    scale <- c(1 / (1 + exp(- 4 * x$proportional_spring_flow)))
    scale <- c(1 / (1 + exp(- 6 * x$proportional_summer_flow)))
    scale <- c(1 / (1 + exp(- 2 * x$proportional_winter_flow)))

    # and some negative effects of extremes
    scale <- c(1 / (1 + exp(-500 * (1 / x$proportional_maximum_flow))))
    scale <- c(1 / (1 + exp(-750 * (1 / x$spawning_flow_variability))))

    # add CTF effects
    scale <- c(1 / (1 + 100 * exp(-200 * (1 / x$ctf_duration))))

    # and return
    mat * scale

  }

  # combine into a covariates object
  covars <- covariates(
    masks = combine(survival(popmat, dims = 2), transition(popmat, dims = 1)),
    funs = survival_effects
  )

  # define environmental stochasticity
  envstoch <- environmental_stochasticity(
    masks = list(
      masks = combine(survival(popmat, dims = 2), transition(popmat, dims = 1)),
      reproduction(popmat)
    ),
    funs = list(
      survival_gen,
      reproduction_gen
    )
  )

  # define demographic stochasticity
  demostoch <- demographic_stochasticity(
    masks = list(
      all_classes(popmat)
    ),
    funs = list(
      function(x, ...) rpois(n = length(x), lambda = x)
    )
  )

  # use density_dependence_n to include translocations
  dd_n_masks <- list(
    all_classes(popmat, dims = 2)
  )
  dd_n_fns <- list(
    add_remove
  )
  dens_depend_n <- density_dependence_n(
    masks = dd_n_masks,
    funs = dd_n_fns
  )

  # collect arguments for species if required
  arguments <- get_args(
    "platypus",
    n = n,
    ntime = ntime,
    start = start,
    end = end
  )

  # return template
  list(
    dynamics = list(
      matrix = popmat,
      covariates = covars,
      environmental_stochasticity = envstoch,
      demographic_stochasticity = demostoch,
      density_dependence = dens_depend,
      density_dependence_n = dens_depend_n
    ),
    arguments = arguments
  )

}

#' @rdname platypus
#'
#' @importFrom stats pnorm rnorm runif
#'
#' @export
#'
args_platypus <- function(n = 0, ntime = 50, start = 1, end = 1) {

  # helper to define removals process
  define_translocation <- function(start, end, n) {

    # set up a sequence of iterations at which individuals are removed
    if (length(n) == 1) {
      x <- rep(0, ntime)
      x[start:end] <- n
      n <- x
    } else {
      if (length(n) != ntime) {
        stop("if n is a vector it must have one value for each time step",
             call. = FALSE)
      }
    }

    # define this as a function
    translocate <- function(obj, pop, iter) {

      # return
      list(
        n = n[iter], add = TRUE
      )

    }

    # return
    translocate

  }


  # return named list of args
  list(

    # to define additions through translocation
    density_dependence_n = list(
      define_translocation(
        start = start, end = end, n = n
      )
    )

  )

}
