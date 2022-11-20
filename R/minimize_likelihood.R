#' @title minimize_likelihood
#'
#' @description Internal function for finding the parameter values that minimize the likelihood function.
#'
#' @noRd

minimize_likelihood <- function(n, x, dy, o, tol) {
  mu_x <- (t(x) %*% o) / n
  
  x_c <- x - o %*% t(mu_x)
  Sigma_x <- (t(x_c) %*% x_c) / n
  
  s <- minimize_1D_likelihood(n, x_c, dy, o, tol)
  
  tmp <- calculate_selected(n, x_c, dy, s, o)
  Sigma_s <- tmp[[1]]
  mu_eta <- tmp[[2]]
  fg <- (tmp[[4]] / tmp[[5]]) * tmp[[3]]
  
  return(list(mu_x, Sigma_x, s, Sigma_s, mu_eta, fg))
}