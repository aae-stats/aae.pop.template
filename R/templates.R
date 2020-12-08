#' @name templates
#' @title Parameterised population dynamics objects
#' @description Use pre-defined population dynamics objects to
#'   define a matrix model for a species with known parameters.
NULL

# internal method
get_template <- function(sp, ...) {

  # unpack dots
  arg_list <- list(...)

  # draw up relevant parameters based on corrected species name
  sp <- parse_species(sp)
  all_parameters <- do.call(get(paste0("template_", sp)), arg_list)

  # return collated dynamics object
  do.call(dynamics, all_parameters)

}

#' @rdname templates
#'
#' @export
#'
#' @importFrom stats rnorm
#'
#' @param k carrying capacity
#' @param reproductive integer vector specifying reproductive
#'   age classes or stages in the population matrix
#' @param system ecosystem type defining elements of population
#'   dynamics. Currently implemented for \code{macquarie_perch},
#'   which has different covariate effects in lakes and rivers
#'
#' @param \dots additional arguments passed to templates or args
#'   functions
#'
#' @details These functions (e.g. \code{murray_cod}) return a collated
#'   \code{\link{dynamics}} object parameterised with values based
#'   on existing data sets and published works. Currently implemented
#'   species are: Murray cod (*Maccullochella peelii*) and Macquarie
#'   perch (*Macquaria australasica*).
#'
#'   In some cases, additional arguments might need to be passed to
#'   \code{\link{simulate}}, and it may be useful to define these as
#'   part of the template. The \code{get_args} function
#'   handles this situation. User-defined arguments functions
#'   should be specified with \code{args_} followed by
#'   a species name or identifier (e.g \code{args_my_species}).
#'   Arguments functions should return a series of named lists
#'   for any of \code{args}, \code{args.dyn}, or \code{args.fn} (see
#'   \code{\link{simulate}} for descriptions of these terms).
#'
#'   Additional arguments are currently required and provided
#'   for Macquarie perch through a call to
#'   \code{get_args("macquarie_perch")}. This function has no
#'   visible arguments but will accept the following options:
#'
#'   - n: an integer, vector of integers, or list specifying the
#'     number of young-of-year, 2+, and adult fish to add or
#'     remove in any given year. If a single integer is provided,
#'     all age classes are assumed to have the same number of
#'     additions or removals. If a vector is provided, it must
#'     have one value for each age class. If a list is provided,
#'     it must have one element for each age class, with a
#'     value for each time step. Addition versus removal is
#'     controlled with \code{add}. Defaults to
#'     \code{c(0, 0, 0))}, which will specify no additions
#'     or removals.
#'
#'   - ntime: number of time steps over which individuals
#'     are added or removed if \code{n} is an integer or vector.
#'     Defaults to \code{50}.
#'
#'   - start: time step at which additions or removals start
#'     if \code{n} is an integer or vector. Defaults to
#'     \code{c(1, 1, 1)}.
#'
#'   - end: time step at which additions or removals finish
#'     if \code{n} is an integer or vector. Defaults to
#'     \code{c(1, 1, 1)}.
#'
#'   - add: logical specifying whether individuals are added
#'     or removed from the population. Defaults to \code{TRUE},
#'     which sets additions.
#'
#'   - allee_strength: strength of Allee effect. Defaults to
#'     \code{1}.
#'
#'   - contributing_min: minimum proportion of adult females
#'     contributing to reproduction at any given time step.
#'     Defaults to \code{0.75}.
#'
#'   - contributing_max: maximum proportion of adult females
#'     contributing to reproduction at any given time step.
#'     Defaults to \code{1.0}.
#'
#'   - recruit_failure: probability of complete recruitment
#'     failure in any given year. Defaults to \code{0}.
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


# convert species names to standardised common names
#  (could replace with sci names if ambiguous)
parse_species <- function(sp) {

  # create lookup table with common synonyms
  sp_list <- list(
    "murraycod" = "murraycod",
    "maccullochellapeelii" = "murraycod",
    "troutcod" = "troutcod",
    "maccullochellamacquariensis" = "troutcod",
    "goldenperch" = "goldenperch",
    "macquariaambigua" = "goldenperch",
    "silverperch" = "silverperch",
    "bidyanusbidyanus" = "silverperch",
    "macquarieperch" = "macquarieperch",
    "macquariaaustralasica" = "macquarieperch",
    "australiansmelt" = "australiansmelt",
    "retropinnasemoni" = "australiansmelt",
    "commongalaxias" = "commongalaxias",
    "galaxiasmaculatus" = "commongalaxias"
  )

  # clean sp to remove underscores or spaces or hyphens
  sp_clean <- gsub("_|-| ", "", sp)

  # pull out a match if there, otherwise check the template
  #   function exists and return as is (error otherwise)
  if (sp_clean %in% names(sp_list)) {
    sp <- sp_list[sp_clean]
  } else {
    if (!exists(paste0("template_", sp))) {
      stop(sp, " does not have a matching template function",
           call. = FALSE)
    }
  }

  # return appropriate species
  sp

}
