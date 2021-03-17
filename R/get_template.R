#' @name get_template
#'
#' @title Compile population dynamics objects from a template
#'
#' @description Function to compile population dynamics objects from
#'   a provided template
#'
#' @param sp character naming the species to be collected. Will be
#'   partially matched to included species
#' @param \dots additional arguments passed to the species template
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
#'   passed to \code{\link[aae.pop]{simulate}}. This process is
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
#' sims <- simulate(mp, nsim = 100, args = mp_args)
#'
#' # plot the simulated values
#' plot(sims)
get_template <- function(sp, ...) {

  # collate arguments
  args <- list(...)

  # get relevant parameters based on corrected species name
  sp <- check_species_template(sp)
  template <- do.call(get(paste0("template_", sp)), args)

  # collate dynamics and arguments objects
  template <- list(
    dynamics = do.call(dynamics, template$dynamics),
    arguments = template$arguments
  )

  # set class and return
  as_template(template)

}

# internal function: check whether a species has a corresponding
#   population model template
check_species_template <- function(x) {

  # currently implemented species
  sp_list <- c("murray_cod", "macquarie_perch", "platypus")

  # give x a fighting chance with fuzzy matching
  if (any(agrepl(x, sp_list)))
    x <- sp_list[agrepl(x, sp_list)]

  # error if species not known
  if (!x %in% sp_list) {
    stop(x, " is not a defined population model in the ",
         "aae.pop.templates package",
         call. = FALSE)
  }

  # return corrected species name
  x

}

# internal function: set template class
as_template <- function(x) {
  type <- "list"
  as_class(x, name = "template", type = type)
}
