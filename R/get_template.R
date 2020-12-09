#' @name get_template
#'
#' @title Compile population dynamics objects from a template
#'
#' @description Function to compile population dynamics objects from
#   a provided template
#'
#' @export
#'
#' @details In most cases, \code{get_template} will not be required
#'   because species templates are defined with wrappers
#'   (e.g., \code{macquarie_perch}). This function is exported
#'   primarily to support user-defined species templates.
#'
#'   When defining templates, species templates must be contained
#'   in a separate function named \code{template_species_name}, where
#'   \code{species_name} can be any unique name or code used to
#'   describe a given species. This function must return a named
#'   list of processes to be passed to \code{\link[aae.pop]{dynamics}}.
#'
#'   User-defined species templates can optionally include arguments
#'   passed to \code{\link[aae.pop]{simulate}}. This processed is
#'   described in \code{\link{get_args}}.
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
#' sims <- simulate(mc, nsim = 100, args = mp_args)
#'
#' # plot the simulated values
#' plot(sims)
get_template <- function(sp, ...) {

  # unpack dots
  arg_list <- list(...)

  # draw up relevant parameters based on corrected species name
  sp <- parse_species(sp)
  all_parameters <- do.call(get(paste0("template_", sp)), arg_list)

  # return collated dynamics object
  do.call(dynamics, all_parameters)

}
