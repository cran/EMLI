#' @title minimize_1D_likelihood
#'
#' @description Internal function for minimizing the one-dimensional likelihood function by golden section search.
#'
#' @noRd

minimize_1D_likelihood <- function(n, x_c, dy, o, tol) {
  solvephi <- (sqrt(5) - 1) / 2
  solvephi2 <- (3 - sqrt(5)) / 2

  a <- 0
  b <- 1

  l0 <- calculate_1D_likelihood(n, x_c, dy, 0, o)
  l1 <- calculate_1D_likelihood(n, x_c, dy, 1, o)

  h <- 1
  steps <- ceiling(log(tol) / log(solvephi))

  c <- solvephi2
  d <- solvephi

  lc <- calculate_1D_likelihood(n, x_c, dy, c, o)
  ld <- calculate_1D_likelihood(n, x_c, dy, d, o)

  for (i in 1:steps) {
    if (lc < ld) {
      b <- d
      d <- c
      ld <- lc
      h <- solvephi * h
      c <- a + solvephi2 * h
      lc <- calculate_1D_likelihood(n, x_c, dy, c, o)

    }

    else {
      a <- c
      c <- d
      lc <- ld
      h <- solvephi * h
      d <- a + solvephi * h
      ld <- calculate_1D_likelihood(n, x_c, dy, d, o)

    }

  }

  if (lc < ld) {
    s <- (a + d) / 2

  }

  else {
    s <- (c + b) / 2

  }

  ls <- calculate_1D_likelihood(n, x_c, dy, s, o)

  minimum <- min(ls, l0, l1)

  if (minimum == ls) {
    return(s)
  }
  else if (minimum == l0) {
    return(0)
  }
  else {
    return(1)
  }

}
