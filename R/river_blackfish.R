#' @name river_blackfish
#' @title Parameterised river blackfish population model
#' @description Use a pre-defined population dynamics object to
#'   simulate population dynamics of river blackfish
#'   (*Gadopsis marmoratus*). Model parameters are based
#'   on existing data sets and published studies on river blackfish.
NULL

#' @rdname river_blackfish
#'
#' @importFrom stats rnorm rpois plogis qlogis
#' @import aae.pop
#'
#' @export
#'
#' @param k carrying capacity. Defaults to \code{1000}
#' @param ntime number of time steps used in population
#'   simulation. Defaults to \code{50}
#'
#' @details The \code{river_blackfish} population model is a
#'   stage-structured model with 4 classes and includes negative
#'   density dependence, environmental and demographic
#'   stochasticity, and optional associations with hydrological
#'   conditions and instream and riparian habitat.
#'
#'   The overarching model structure was based the life history of river
#'   blackfish. The primary model parameters (survival, fecundity,
#'   density dependence) were estimated from observed counts of river
#'   blackfish in 4 length classes (0-80 mm, 80-140 mm, 140-200 mm, > 200 mm).
#'   These length classes roughly correspond with 0+, 1+, 2+, and older
#'   individuals based on preliminary assessments of survey data in
#'   Victorian rivers. The effects of hydrology and habitat were based on
#'   (minimal) published literature on the species and are subject to
#'   moderate levels of uncertainty.
#'
#'   The river_blackfish population template does not
#'   currently require additional arguments.
#'
#' @examples
#' # define a basic model for river blackfish
#' p <- river_blackfish()
#'
#' # simulate from this model
#' sims <- simulate(p, nsim = 100)
#'
#' # plot the simulated values
#' plot(sims)
river_blackfish <- function(k = 1000, ntime = 50) {
  get_template(
    sp = "river_blackfish",
    k = k,
    ntime = ntime
  )
}

# internal function: define species defaults
template_river_blackfish <- function(k = 1000, ntime = 50) {

  # assume five length classes (in mm): 0-80, 80-140, 140-200, >200,
  #   reproductive for all but the first size class
  ## UPDATE: trial 11-age class model
  nstage <- 11
  reproductive <- seq_len(nstage)[-1]

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
  early_surv <- 0.1
  fecundity <- seq(150, 300, length = nstage - 1L)

  # define base matrix
  # fecundity of ~ 50-300, but can have multiple breeding attempts
  popmat <- matrix(0, nrow = nstage, ncol = nstage)
  popmat[reproduction(popmat, dims = reproductive)] <-
    sex_ratio * early_surv * fecundity
  popmat[transition(popmat)] <-
    c(0.2, 0.4, 0.5, 0.55, 0.4, 0.3, 0.2, 0.2, 0.1, 0.05)

  # define contest competition
  k_female <- k * 0.5     # assumes 50:50 F:M sex ratio in adults
  bh <- aae.pop::beverton_holt(k = k_female)
  dens_depend <- density_dependence(
    reproduction(popmat, dims = reproductive),
    bh
  )

  # covariate effects based on low water temperatures (> 5C), high
  #   average flows in recent years, and habitat condition
  #   (snags or rocky cover)
  survival_effects <- function(
    mat,
    x,
    ...
  ) {

    # default to mean survival with no negative impacts of
    #   water temperatures, exotic species, or habitat condition
    scale <- 1

    # set a small threshold so the qlogis/plogis transformations work
    eps <- 1e-5

    # survival requires water temperatures above 5C in winter
    if (!is.null(x$nday_lt5))
      scale <- ifelse(x$nday_lt5 > 5, 0.05 * scale, scale)

    # and habitat condition for shelter and nesting (instream hab first)
    if (!is.null(x$instream_cover))
      scale <- x$instream_cover * scale

    # and habitat condition for shelter and nesting (then overhang)
    if (!is.null(x$veg_overhang))
      scale <- x$veg_overhang * scale

    # make sure scale is between 0 and 1
    scale <- ifelse(scale > 1, 1 - eps, ifelse(scale < 0, eps, scale))

    # and return
    mat * scale

  }

  # temperature effects on reproduction:
  #    - require water temperatures above 16C (breeding threshold)
  #        during the spawning period (Oct-Dec) and assume a large
  #        decline in breeding output if temperatures are below this
  #        threshold prior to Dec
  #    - assume high winter, spring, and antecedent flows result in high
  #        reproductive output, along with low summer flows and limited
  #        flow variability
  #    - coldwater during summer (< 18C) has a negative effect on
  #        reproductive output
  reproduction_effects <- function(
    mat,
    x,
    coefs = NULL,
    temperature_coefficient = 0.1,
    coldwater_coefficient = 0.05,
    ...
  ) {

    # set default coefs if not specified
    if (is.null(coefs)) {

      # defaults
      coefs <- c(-100, 20, 5, 20, 20)

    } else {

      # need three coefficients
      if (length(coefs) != 5) {
        stop(
          "coefs must include five values, one for each of ",
          "flow variability, spring flows, summer flows, winter flows, ",
          "and antecedent flows",
          .call = FALSE
        )
      }

    }

    # default to mean survival with no impacts of flows or
    #   water temperature
    scale <- 1

    # set a small threshold so the qlogis/plogis transformations work
    eps <- 1e-5

    # spawning requires temperatures above 16C
    if (!is.null(x$nday_gt16))
      scale <- 1 / (1 + exp(-temperature_coefficient * (x$nday_gt16 - 30)))

    # and above 18C for later in spring/summer
    if (!is.null(x$nday_lt18))
      scale <- 1 / (1 + exp(coldwater_coefficient * (x$nday_lt18 - 100)))

    # calculate scaling factor by year
    metrics <- c(
      x$spawning_flow_variability,
      x$proportional_spring_flow,
      x$proportional_summer_flow,
      x$proportional_winter_flow,
      x$antecedent_flow
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

    # and multiply by previously calculated scaling effect
    scale <- effect * scale

    # make sure scale is between 0 and 1
    scale <- ifelse(scale > 1, 1 - eps, ifelse(scale < 0, eps, scale))

    # and return
    mat * scale

  }

  # combine into a covariates object
  covars <- covariates(
    masks = list(
      aae.pop::transition(popmat),
      aae.pop::reproduction(popmat, dims = reproductive)
    ),
    funs = list(
      survival_effects,
      reproduction_effects
    )
  )

  # define environmental stochasticity
  envstoch <- environmental_stochasticity(
    masks = list(
      aae.pop::transition(popmat),
      aae.pop::reproduction(popmat, dims = reproductive)
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
    "river_blackfish",
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

#' @rdname river_blackfish
#'
#' @export
#'
args_river_blackfish <- function(ntime = 50) {

  # return named list of args
  list()

}
