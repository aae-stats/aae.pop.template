#' @name murray_darling_rainbowfish
#' @title Parameterised Murray-Darling rainbowfish population model
#' @description Use a pre-defined population dynamics object to
#'   simulate population dynamics of Murray-Darling rainbowfish
#'   (*Melanotaenia fluviatilis*). Model parameters are based
#'   on existing data sets and published studies on Murray-Darling rainbowfish.
NULL

#' @rdname murray_darling_rainbowfish
#'
#' @importFrom stats rnorm rpois plogis qlogis
#' @import aae.pop
#'
#' @export
#'
#' @param k carrying capacity. Defaults to \code{10000}
#' @param ntime number of time steps used in population
#'   simulation. Defaults to \code{50}
#'
#' @details The \code{murray_darling_rainbowfish} population model is a
#'   stage-structured model with 5 classes and includes negative
#'   density dependence, environmental and demographic
#'   stochasticity, and optional associations with hydrological
#'   conditions, instream and riparian habitat, and predators.
#'
#'   The overarching model structure was based the life history of Murray-
#'   Darling rainbowfish. The primary model parameters (survival, fecundity,
#'   density dependence) were estimated from observed counts of rainbowfish
#'   in 5 length classes (0-30 mm, 30-40 mm, 40-50 mm, 50-70 mm, > 70 mm).
#'   The effects of hydrology, habitat, and predation were based on
#'   (minimal) published literature on the species and are subject to
#'   moderate levels of uncertainty.
#'
#'   The murray_darling_rainbowfish population template does not
#'   currently require additional arguments.
#'
#' @examples
#' # define a basic model for Murray-Darling rainbowfish
#' p <- murray_darling_rainbowfish()
#'
#' # simulate from this model
#' sims <- simulate(p, nsim = 100)
#'
#' # plot the simulated values
#' plot(sims)
murray_darling_rainbowfish <- function(k = 10000, ntime = 50) {
  get_template(
    sp = "murray_darling_rainbowfish",
    k = k,
    ntime = ntime
  )
}

# internal function: define species defaults
template_murray_darling_rainbowfish <- function(k = 10000, ntime = 50) {

  # assume five length classes (in mm): 0-30, 30-40, 40-50, 50-70, >70
  nstage <- 5

  # reproductive for all but the first size class
  reproductive <- c(2L:5L)

  # define  a survival function
  survival_gen <- function(mat, ...) {
    plogis(rnorm(length(mat), mean = qlogis(mat), sd = 0.5 * abs(qlogis(mat))))
  }

  # define a reproduction function
  reproduction_gen <- function(
    mat,
    ...
  ) {

    # simulate realised fecundity in a given year
    rpois(length(mat), lambda = mat)

  }

  # basic reproduction settings
  sex_ratio <- 0.5
  early_surv <- 0.3
  fecundity <- c(0, 50, 100, 150, 200)

  # define base matrix
  # fecundity of ~ 50-300, but can have multiple breeding attempts
  popmat <- rbind(
    sex_ratio * early_surv * fecundity,
    c(0.3, 0.25, 0,    0,    0),
    c(0,   0.35, 0.31, 0,    0),
    c(0,   0,    0.29, 0.35, 0),
    c(0,   0,    0,    0.38, 0.41)
  )

  # define contest competition
  k_female <- k * 0.5     # assumes 50:50 F:M sex ratio in adults
  bh <- aae.pop::beverton_holt(k = k_female)
  dens_depend <- density_dependence(
    reproduction(popmat, dims = reproductive),
    bh
  )

  # covariate effects based on low water temperatures, the
  #   presence of exotic predators, and habitat condition
  survival_effects <- function(
    mat,
    x,
    ...
  ) {

    # default to mean survival with no negative impacts of
    #   water temperatures, exotic species, or habitat condition
    scale <- 1

    # survival requires water temperatures above 10C in winter
    if (!is.null(x$nday_lt10))
      scale <- ifelse(x$nday_lt10 > 5, 0, scale)

    # and negative effects of redfin presence
    if (!is.null(x$redfin))
      scale <- ifelse(x$redfin, 0.5 * scale, scale)

    # and habitat condition for spawning
    if (!is.null(x$instream_cover))
      scale <- x$instream_cover * scale

    # and return
    mat * scale

  }

  # temperature effects on reproduction:
  #    consider number of days above 20C (breeding threshold)
  #    and assume one breeding attempt every 30 days where
  #    temperatures are suitable
  reproduction_effects <- function(
    mat,
    x,
    ...
  ) {

    # default to a single breeding attempt with no negative
    #   impacts of exotic species or habitat condition
    scale <- 1

    # the number of breeding attempts is linked to the number of
    #    days where water temperature exceeds 20C
    if (!is.null(x$nday_gt20))
      scale <- floor(x$nday_gt20 / 30)

    # and negative effects of gambusia presence
    if (!is.null(x$gambusia))
      scale <- ifelse(x$gambusia, 0.3 * scale, scale)

    # and habitat condition for spawning
    if (!is.null(x$instream_cover))
      scale <- x$instream_cover * scale

    # and return
    mat * scale

  }

  # combine into a covariates object
  covars <- covariates(
    masks = list(
      aae.pop::combine(
        aae.pop::survival(popmat, dims = 2L:5L),
        aae.pop::transition(popmat, dims = 1L:4L)
      ),
      aae.pop::reproduction(popmat, dims = 2L:5L)
    ),
    funs = list(
      survival_effects,
      reproduction_effects
    )
  )

  # define environmental stochasticity
  envstoch <- environmental_stochasticity(
    masks = list(
      aae.pop::combine(
        aae.pop::survival(popmat, dims = 2L:5L),
        aae.pop::transition(popmat, dims = 1L:4L)
      ),
      aae.pop::reproduction(popmat, dims = 2L:5L)
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

  # collect arguments for species if required
  arguments <- get_args(
    "murray_darling_rainbowfish",
    ntime = ntime
  )

  # return template
  list(
    dynamics = list(
      matrix = popmat,
      covariates = covars,
      environmental_stochasticity = envstoch,
      demographic_stochasticity = demostoch,
      density_dependence = dens_depend
    ),
    arguments = arguments
  )

}

#' @rdname murray_darling_rainbowfish
#'
#' @export
#'
args_murray_darling_rainbowfish <- function(ntime = 50) {

  # return named list of args
  list()

}
