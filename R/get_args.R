#' @name get_args
#'
#' @title Define additional arguments relevant to specific population
#'   model templates
#'
#' @description Define additional arguments used to control the
#'   simulation of population dynamics in \code{\link[aae.pop]{simulate}}.
#'
#' @export
#'
#' @details The \code{get_args} function can be used to access
#'   previously defined arguments to be passed to
#'   \code{\link[aae.pop]{simulate}} for any given species.
#'   In general, arguments for a particular species can be
#'   accessed with \code{get_args("species_name")}, where
#'   \code{species_name} is replaced with the function name used
#'   to define the \code{dynamics} object for that species. For example,
#'   arguments for Macquarie perch (function \code{macquarie_perch})
#'   can be accessed with \code{get_args("macquarie_perch")}.
#'
#'   When defining templates, arguments are best contained in a
#'   separate function named \code{args_species_name}, where
#'   \code{species_name} is described above. This function must
#'   return a named list with elements corresponding to any of
#'   the processes required by \code{dynamics} or
#'   \code{multispecies}.
#'
#'   The \code{get_args} function has been deprecated, and is now
#'   called directly within a call to \code{get_template}.
#'
#' @examples
#' # define a basic model for Macquarie perch with
#' #   carrying capacity = 1000
#' mp <- macquarie_perch(k = 1000)
#'
#' # define required arguments
#' mp_args <- get_args("macquarie_perch")
#'
#' # simulate from this model
#' sims <- simulate(mp, nsim = 100, args = mp_args)
#'
#' # plot the simulated values
#' plot(sims)
get_args <- function(sp, ...) {

  # draw up relevant parameters based on corrected species name
  sp <- check_species_args(sp)

  # initialise NULL output
  out <- NULL

  # check if species has arguments available
  available <- exists(
    paste0("args_", sp),
    envir = getNamespace("aae.pop"),
    mode = "function"
  )

  # collapse dots into a list
  arg_list <- list(...)

  # get arguments if available
  if (available)
    out <- do.call(get(paste0("args_", sp)), arg_list)

  # return
  out

}

# internal function: check whether a species has corresponding
#   arguments in their population model template
check_species_args <- function(x) {

  # currently implemented species
  sp_list <- c("macquarie_perch", "murray_cod", "platypus", "estuary_perch")

  # give x a fighting chance with fuzzy matching
  if (any(agrepl(x, sp_list)))
    x <- sp_list[agrepl(x, sp_list)]

  # error if species not known
  if (!x %in% sp_list) {
    stop(x, " does not have arguments defined in the ",
         "aae.pop.templates package",
         call. = FALSE)
  }

  # return species name
  #   (only needed if using partial matching)
  x

}
