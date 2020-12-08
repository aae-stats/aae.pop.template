# internal function to initialise a vector of zeros
#   and fill a subset of values
zeros_and_fill <- function(x, start, end, len) {
  out <- rep(0, len)
  out[start:end] <- x
  out
}

# internal function to handle translocations or stocking
add_remove <- function(
  pop,
  n,
  add = TRUE,
  ...
) {

  # only add/remove individuals if required
  if (n > 0) {

    # check there are enough individuals to remove
    if (sum(pop) < n & !add) {
      n <- sum(pop)
      warning("removal required more individuals than ",
              "were available; reduced to ",
              sum(pop),
              " individuals",
              call. = FALSE)
    }

    # are we removing?
    if (!add) {

      # if so, expand n to remove from random age classes
      n_by_age <- rep(seq_along(pop), times = pop)
      idx <- sample.int(
        length(n_by_age), size = n, replace = FALSE
      )
      n_change <- table(n_by_age[idx])

    } else {

      # otherwise, add to age classes at random
      n_change <- table(
        sample(seq_along(pop), size = n, replace = TRUE)
      )

    }

    n_expanded <- rep(0, length(pop))
    names(n_expanded) <- as.character(seq_along(pop))
    n_expanded[names(n_change)] <- n_change

    # are we adding or removing?
    if (add)
      n_expanded <- -n_expanded

    # update pop abundances
    pop <- pop - n_expanded

  }

  # return
  pop

}
