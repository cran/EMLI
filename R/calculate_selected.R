#' @title calculate_selected
#'
#' @description Internal function for calculating certain parameter estimates and a few complementary terms that all depend on s.
#'
#' @noRd

calculate_selected <- function(n, x_c, dy, s, o) {
  tmp1 <- calculate_recursive(n, o, o, s)[[1]]
  tmp2 <- calculate_recursive(n, x_c, o, s)[[1]]
  tmp3 <- solve(calculate_recursive(n, x_c, x_c, s)[[1]])

  Sigma_s <-
    solve(diag(ncol(x_c)) - (tmp3 %*% tmp2 %*% t(tmp2)) / tmp1) %*% tmp3 %*% calculate_recursive(n, x_c, dy - (calculate_recursive(n, dy, o, s)[[1]] / tmp1) * o, s)[[1]]

  tmp <- dy - x_c %*% Sigma_s
  mu_eta <- calculate_recursive(n, tmp, o, s)[[1]] / tmp1

  tmp <- tmp - mu_eta * o
  tmp <- calculate_recursive(n, tmp, tmp, s)

  return(list(Sigma_s, mu_eta, tmp[[1]], tmp[[2]], tmp[[3]]))

}
