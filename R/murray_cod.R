murray_cod <- function(k = 20000, ...) {
  get_template(sp = "murraycod", k = k, ...)
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

# internal function: define macquarie perch arguments
args_macquarieperch <- function(
  n = c(0, 0, 0),
  ntime = 50,
  start = c(1, 1, 1),
  end = c(1, 1, 1),
  add = TRUE,
  allee_strength = 1,
  contributing_min = 0.75,
  contributing_max = 1.0,
  recruit_failure = 0,
  genetic_factor = 1.0
) {

  # expand n, start, end if required
  if (length(n) == 1)
    n <- rep(n, 3)
  if (length(start) == 1)
    start <- rep(start, 3)
  if (length(end) == 1)
    end <- rep(end, 3)

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
      if (!all.equal(sapply(n, length), ntime)) {
        stop("if n is a list, each element must be a vector ",
             "with one value for each time step",
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
        add = add
      )

    }

    # return
    translocate

  }

  early_surv <- c(0.5, 0.013, 0.13)

  # helper to calculate real-valued parameters for survival
  #   simulation
  transform_survival <- function(obj, pop, iter) {

    # pull out the population matrix in the current time step
    mat <- obj$matrix
    if (is.list(mat))
      mat <- mat[[iter]]

    # wrap up all survival means and SDs, including early life
    #  (this allows a single call to unit_to_real, which is slow)
    survival_mean <- c(
      genetic_factor * early_surv,  # early life survival with gene mixing
      mat[transition(mat)]    # from population matrix in current time step
    )
    survival_sd <- c(
      0.1, 0.007, 0.028,  # early life
      0.05, 0.09, 0.11, 0.10, 0.10, 0.07, 0.08, 0.08, 0.08,
      0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08,
      0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.07, 0.06, 0.05
    )

    # convert unit interval to real line equivalents
    out <- unit_to_real(
      unit_mean = survival_mean,
      unit_sd = survival_sd
    )

    # separate early life from other estimates
    idx <- seq_len(nrow(out)) > 3

    # return
    list(
      mean_real = out[idx, 1],    # for survival_gen
      sd_real = out[idx, 2],      # for survival_gen
      early_mean = out[!idx, 1],  # for reproduction_gen
      early_sd = out[!idx, 2]     # for reproduction_gen
    )

  }

  # return named list of args
  list(

    # set as 1 (default) or 2
    density_dependence = list(allee_strength = allee_strength),

    # set contributing as random uniform on 0.75-1.0 by default
    # set recruit_failure at 0 by default
    # add function to pre-transform unit to real and back
    environmental_stochasticity = list(
      contributing_min = contributing_min,
      contributing_max = contributing_max,
      recruit_failure = recruit_failure,
      transform_survival
    ),

    # to include additions or removals of individuals
    density_dependence_n = list(
      define_removals(
        start = start, end = end, n = n, add = add
      )
    )

  )

}
