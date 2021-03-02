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
#'
#' @details The \code{platypus} population model is a
#'   stage-structured model with 2 classes and includes negative
#'   density dependence, environmental and demographic
#'   stochasticity, and optional associations with hydrological
#'   conditions, instream and riparian habitat, and predators.
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
platypus <- function(k = 400) {
  get_template(sp = "platypus", k = k)
}

# internal function: define species defaults
template_platypus <- function(k = 400) {

  # how many stages are we going to work with?
  nstage <- 2

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
    proportion_female = 0.5,
    proportion_breeding = 0.62,
    ...
  ) {

    # simulate realised fecundity in a given year
    n_offspring <- rpois(1, lambda = mat[1, 2])

    # estimate and return reproduction estimates
    n_offspring * proportion_female * proportion_breeding

  }

  # define base matrix
  popmat <- matrix(0, nrow = nstage, ncol = nstage)

  # add survival
  popmat[2, 1] <- 0.29
  popmat[2, 2] <- 0.89

  # add reproduction
  popmat[1, 2] <- 0.5 * 0.62 * 1.5

  # define BH density dependence
  k_female <- k * 0.65     # assumes 65:35 F:M sex ratio in adults
  dens_depend <- density_dependence(
    combine(survival(popmat), reproduction(popmat)),
    beverton_holt(k = k_female)
  )

  # covariate effects based on six-month cumulativee discharge
  #    - associations modified from Bino et al. (2015, Sci Rep; DOI: 10.1038/srep16073)
  #    - modifications to retain same overall shape of curves but based on
  #        proportional changes in survival, assuming same relative change in
  #        juveniles and adults
  #    - x is cumulative discharge over preceding 6 months, in GL
  survival_effects <- function(mat, x, ...) {

    # set parameters
    alpha <- -7.38
    beta <- 0.0008
    gamma <- -0.0001
    delta <- 0.013
    offset <- 1.05
    ave_weight <- 680

    # calculate proportional change in survival
    term <- alpha + beta * x + gamma * (x ^ 2) + delta * ave_weight + offset
    effect <- exp(term) / (1 + exp(term))
    effect <- effect / max(effect)

    # calculate and return absolute change in survival
    mat * effect

  }

  # combine into a covariates object
  covars <- covariates(
    masks = survival(popmat),
    funs = survival_effects
  )

  # define environmental stochasticity
  envstoch <- environmental_stochasticity(
    masks = list(
      survival(popmat),
      reproduction(popmat)
    ),
    funs = list(
      survival_gen,
      reproduction_gen
    )
  )

  # optional: use density_dependence_n to include predation
  #   or translocations
  # dd_n_masks <- list(
  #   all_classes(popmat)
  # )
  # dd_n_fns <- list(
  #   add_remove
  # )
  # dens_depend_n <- density_dependence_n(
  #   masks = dd_n_masks,
  #   funs = dd_n_fns
  # )

  # nolint start
  # print a message to ensure args are included
  message('Arguments are required to simulate platypus dynamics.\n',
          'Default arguments can be accessed with get_args("platypus").')
  # nolint end

  # return template
  list(
    matrix = popmat,
    covariates = covars,
    environmental_stochasticity = envstoch,
    demographic_stochasticity = NULL,
    density_dependence = dens_depend,
    density_dependence_n = NULL
  )

}

#' @rdname platypus
#'
#' @importFrom stats pnorm rnorm runif
#'
#' @export
#'
#' @param n an integer, vector of integers, or list specifying the
#'   number of juvenile and adult platypus to add or
#'   remove in any given year. If a single integer is provided,
#'   all age classes are assumed to have the same number of
#'   additions or removals. If a vector is provided, it must
#'   have one value for each age class. If a list is provided,
#'   it must have one element for each age class, with a
#'   value for each time step. Addition versus removal is
#'   controlled with \code{add}. Defaults to
#'   \code{c(0, 0))}, which will specify no additions
#'   or removals
#' @param ntime number of time steps used in population
#'   simulation. Defaults to \code{50} and is required
#'   to ensure simulated additions or removals are defined
#'   for every simulated time step
#' @param start time step at which additions or removals start
#'   if \code{n} is an integer or vector. Defaults to
#'   \code{c(1, 1)}
#' @param end time step at which additions or removals finish
#'   if \code{n} is an integer or vector. Defaults to
#'   \code{c(1, 1)}
#' @param add logical indicating whether individuals are
#'   added or removed from the population. Defaults to
#'   \code{TRUE}, which defines additions (e.g., translocation). Can
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
#' @details The platypus population template requires
#'   several additional arguments and allows several optional
#'   arguments. Arguments can be defined with a call to
#'   \code{get_args("platypus")} and are described
#'   individually above.
#'
args_platypus <- function(...) {

  # helper to calculate real-valued parameters for survival
  #   simulation
  transform_survival <- function(obj, pop, iter) {

    # pull out the population matrix in the current time step
    mat <- obj$matrix
    if (is.list(mat))
      mat <- mat[[iter]]

    # wrap up all survival means and SDs
    survival_mean <- c(mat[2, 1], mat[2, 2])
    survival_sd <- c(0.01, 0.05)

    # convert unit interval to real line equivalents
    out <- unit_to_real(
      unit_mean = survival_mean,
      unit_sd = survival_sd
    )

    # return
    list(
      mean_real = out[, 1],
      sd_real = out[, 2]
    )

  }

  # return named list of args
  list(

    # add function to pre-transform unit to real and back
    environmental_stochasticity = list(
      transform_survival
    )

  )

}
