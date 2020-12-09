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
#' @param k carrying capacity
#'
#' @details The \code{murray_cod} population model is a mixed
#'   age- and stage-structured model with 25 classes and includes negative
#'   density dependence, environmental and demographic
#'   stochasticity, and optional associations with hydrological conditions
#'   in rivers and stocking or angling effects.
#'
#' @examples
#' # define a basic model for Murray cod with
#' #   carrying capacity = 25000
#' mc <- murray_cod(k = 25000)
#'
#' # simulate from this model
#' sims <- simulate(mc, nsim = 100)
#'
#' # plot the simulated values
#' plot(sims)
murray_cod <- function(k = 20000) {
  get_template(sp = "murraycod", k = k)
}

# internal function: define species defaults
template_murraycod <- function(k = 20000) {

  # how many stages are we going to work with?
  nstage <- 25

  # define base matrix
  mat <- matrix(0, nrow = nstage, ncol = nstage)
  survival_mask <- combine(
    transition(mat), survival(mat, dims = nstage)
  )
  mat[survival_mask] <- c(
    0.4790, 0.5846, 0.6552, 0.7054, 0.7431, 0.7722, 0.7954, 0.8144, 0.8301,
    0.8434, 0.8547, 0.8646, 0.8731, 0.8807, 0.8874, 0.8934, 0.8988, 0.9037,
    0.9081, 0.9121, 0.9158, 0.9192, 0.9224, 0.9253, 0.9375
  )
  yoy_surv <- 0.5 * 0.0122 * 0.1225
  reproduction_mask <- reproduction(mat, dims = c(5:nstage))
  mat[reproduction_mask] <- yoy_surv * c(
    3000, 5000, 7000, 9000, 12000, 16000,
    20000, 25000, 30000, 34000, 38000,
    41000, 43000, 45000, 47000, 48000,
    48000, 49000, 49000, 49000, 50000
  )

  # define basic BH density dependence
  biomass_dd <- function(k, dims) {
    function(x, n) {
      sum_n <- sum(n[min(dims):length(n)])
      x * ifelse(sum_n > k, k / sum_n, 1)
    }
  }
  dd_stages <- list(
    c(3:4),
    c(5:7),
    c(8:10),
    c(11:14),
    c(15:25)
  )
  biomass_dd_list <- lapply(dd_stages, biomass_dd, k = k)
  biomass_mask_list <- lapply(
    dd_stages, transition, mat = mat
  )
  biomass_mask_list[[length(biomass_mask_list)]][nstage, nstage] <- TRUE
  dd_fns <- c(
    list(beverton_holt(k = k)), biomass_dd_list
  )
  dd_masks <- c(
    list(reproduction(mat, dims = c(5:nstage))),
    biomass_mask_list
  )
  dd <- density_dependence(dd_masks, dd_fns)

  # basic single variable covariate function
  cov_funs <- function(mat, x) {
    mat * (1 / (1 + exp(-0.5 * (x + 10))))
  }
  cov_masks <- transition(mat)
  covars <- covariates(
    masks = cov_masks,
    funs = cov_funs
  )

  # define environmental stochasticity based on known standard deviations of
  #   parameters
  # survival
  survival_env <- function(x) {

    # define base parameters for SD
    sd_scale <- c(0.2, rep(0.15, 3), rep(0.1, 3), rep(0.075, 3), rep(0.05, 15))
    sd_tmp <- sd_scale * c(
      0.4790, 0.5846, 0.6552, 0.7054, 0.7431, 0.7722, 0.7954, 0.8144, 0.8301,
      0.8434, 0.8547, 0.8646, 0.8731, 0.8807, 0.8874, 0.8934, 0.8988, 0.9037,
      0.9081, 0.9121, 0.9158, 0.9192, 0.9224, 0.9253, 0.9375
    )

    # simulate normal random variates
    out <- rnorm(length(x), mean = x, sd = sd_tmp)

    # check none sit outside 0/1
    out <- ifelse(out > 1, 1, out)
    out <- ifelse(out < 0, 0, out)

    # return
    out

  }
  # reproduction
  reproduction_env <- function(x) {

    # define base parameters for SD
    yoy_surv <- 0.5 * 0.0122 * 0.1225
    sd_tmp <- yoy_surv * c(
      1500, 2400, 3200, 4000, 5300, 6900,
      8400, 10500, 12600, 13900, 15600,
      16800, 17600, 18500, 18800, 19200,
      19200, 19600, 19600, 19600, 20000
    )

    # simulate normal random variates
    out <- rnorm(length(x), mean = x, sd = sd_tmp)

    # check none are negative (but can be > 1)
    out <- ifelse(out < 0, 0, out)

    # return
    out

  }
  # collate into a enviro stochasticity object
  envstoch <- environmental_stochasticity(
    masks = list(survival_mask, reproduction_mask),
    funs = list(survival_env, reproduction_env)
  )

  # set demographic stochasticity
  # survival
  demo_fn <- function(x) {
    rpois(length(x), lambda = x)
  }
  demostoch <- demographic_stochasticity(
    masks = all_classes(mat),
    funs = demo_fn
  )

  # return template
  list(
    matrix = mat,
    covariates = covars,
    environmental_stochasticity = envstoch,
    demographic_stochasticity = demostoch,
    density_dependence = dd,
    density_dependence_n = NULL
  )

}
