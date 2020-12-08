#' @rdname arguments
#'
#' @export
get_args <- function(sp, ...) {

  # draw up relevant parameters based on corrected species name
  sp <- parse_species(sp)

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

