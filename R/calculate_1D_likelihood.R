#' @title calculate_1D_likelihood
#'
#' @description Internal function for calculating the one-dimensional likelihood function that depends on parameter s.
#'
#' @noRd

calculate_1D_likelihood <- function(n, x_c, dy, s, o) {
  tmp <- calculate_selected(n, x_c, dy, s, o)

  h <- log(tmp[[4]] / (tmp[[5]] ^ (1 - 1 / n)))
  g <- tmp[[3]]

  return (h + log(g))

}
