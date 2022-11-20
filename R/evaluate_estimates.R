#' @title evaluate_estimates
#'
#' @description Calculates a discrepancy-function-based metric of accuracy of the statistical measure estimates for the output-differenced version of the one-dimensional cumulative structural equation model with shock-error output measurement equation and assumptions of normality and independence. Suitable when there are no contradictions in the factuals/estimates.
#'
#' @param f A list consisting of 3 elements: 1) the factual Sigma ((m + 1) x (m + 1) matrix of finite numeric elements); 2) the factual sigma_y^2 (vector of length 1, finite numeric element); 3) the factual mu ((m + 1) x 1 matrix of finite numeric elements).
#' @param e Analogous to parameter f but with estimates instead of factuals.
#' @param n The number of time moments used for obtaining parameter e (vector of length  1, finite positive integer).
#'
#' @return Calculated accuracy metric value (vector of length 1, numeric element). The lower the value, the better the accuracy, with 0 indicating perfect accuracy.
#'
#' @examples
#' set.seed(1)
#'
#' m <- 4
#' k <- 2
#'
#' L <- matrix(runif((m + 1) * k, min = -10, max = 10), nrow = m + 1)
#' sigma <- matrix(runif(m + 2, min = 0, max = 10), nrow = m + 2)
#' mu <- matrix(runif(m + 1, min = -10, max = 10), nrow = m + 1)
#'
#' n <- 100
#' data <- generate_data(n, L, sigma, mu)
#'
#' Sigma <- L %*% t(L) + diag(sigma[1:(m + 1), ] ^ 2)
#' sigma_y_squared <- sigma[m + 2, ] ^ 2
#' Sigma[m + 1, m + 1] <- Sigma[m + 1, m + 1] + 2 * sigma_y_squared
#'
#' factual_parameters <- list(Sigma, sigma_y_squared, mu)
#' estimated_parameters <- estimate_parameters(data, 0.00001)
#'
#' evaluate_estimates(factual_parameters, estimated_parameters, n)
#'
#' @export

evaluate_estimates <- function(f, e, n) {
  if (is.list(f) && is.list(e) && is.vector(n)) {
    if (length(f) == 3 && length(e) == 3 && length(n) == 1) {
      if (is.matrix(f[[1]]) &&
          is.matrix(e[[1]]) &&
          is.vector(f[[2]]) &&
          is.vector(e[[2]]) &&
          is.matrix(f[[3]]) && is.matrix(e[[3]])) {
        if (dim(f[[1]])[1] == dim(f[[1]])[2] &&
            dim(f[[1]])[1] == dim(f[[3]])[1] &&
            all.equal(dim(e[[1]]), dim(f[[1]])) &&
            dim(f[[3]])[2] == 1 &&
            all.equal(dim(f[[3]]), dim(e[[3]])) &&
            length(f[[2]]) == 1 && length(e[[2]]) == 1) {
          if (is.numeric(f[[1]]) &&
              is.numeric(f[[2]]) &&
              is.numeric(f[[3]]) &&
              is.numeric(e[[1]]) &&
              is.numeric(e[[2]]) &&
              is.numeric(e[[3]]) && is.numeric(n)) {
            if (all(is.finite(f[[1]])) &&
                is.finite(f[[2]]) &&
                all(is.finite(f[[3]])) &&
                all(is.finite(e[[1]])) &&
                is.finite(e[[2]]) &&
                all(is.finite(e[[3]])) && is.finite(n)) {
              if (n == as.integer(n) && n >= 1) {
                if (isSymmetric(f[[1]]) && isSymmetric(e[[1]])) {
                  if (min(eigen(f[[1]])$values) >= 0 &&
                      min(eigen(e[[1]])$values) >= 0) {
                    m <- nrow(f[[1]]) - 1
                    inv_f1 <- solve(f[[1]])
                    tmp <- f[[2]] * inv_f1[m + 1, m + 1]
                    e_tmp <- e[[2]] * solve(e[[1]])[m + 1, m + 1]
                    if (tmp >= 0 &&
                        tmp <= 0.5 && e_tmp >= 0 && e_tmp <= 0.5) {
                      mu_x <- f[[3]][1:m]
                      mu_dy <- f[[3]][m + 1]
                      inv_Sigma_x <- solve(f[[1]][1:m, 1:m])
                      if (tmp == 0) {
                        s <- 0
                      }
                      else {
                        s <- exp(-acosh(1 / 2 * 1 / tmp))
                      }
                      Sigma_xdy <-
                        inv_Sigma_x %*% f[[1]][1:m, m + 1]
                      s_t <- s ^ 2 + 1
                      
                      e_mu_x <- e[[3]][1:m]
                      e_mu_dy <- e[[3]][m + 1]
                      e_Sigma_x <- e[[1]][1:m, 1:m]
                      if (e_tmp == 0) {
                        e_s <- 0
                      }
                      else {
                        e_s <- exp(-acosh(1 / 2 * 1 / e_tmp))
                      }
                      
                      tmp1 <- 1
                      tmp2 <- 1
                      tmp4 <- 0
                      tmp5 <- n
                      if (n > 1) {
                        for (i in 1:(n - 1)) {
                          a <- 2 * i
                          b <- s ^ a
                          tmp1 <- tmp1 + b
                          tmp2 <- tmp2 + e_s ^ a
                          c <- (n - i) * (i + 1)
                          tmp4 <- tmp4 + s ^ (a - 2) * c
                          tmp5 <- tmp5 + b * c
                        }
                      }
                      tmp3 <- tmp1 - 1
                      a <- 2 * n
                      tmp1 <- tmp1 + s ^ a
                      tmp2 <- tmp2 + e_s ^ a
                      
                      l <-
                        n * log(det(f[[1]]) / det(e[[1]])) + n * log((e_s ^ 2 + 1) / s_t) + log(tmp1 /
                                                                                                  tmp2)
                      t <-
                        (sum(diag(
                          e_Sigma_x %*% inv_Sigma_x
                        )) * s + e[[2]] * inv_f1[m + 1, m + 1] * s_t) * s * (n * tmp3 - s_t * tmp4) /
                        tmp1 + sum(diag(e[[1]] %*% inv_f1)) * s_t * tmp5 / tmp1
                      tmp1 <- 1
                      tmp2 <- 1
                      f1 <- 1
                      f2 <- f1 + s ^ 2
                      if (n > 1) {
                        for (i in 2:n) {
                          f3 <- f2 + s ^ (2 * i)
                          tmp2 <- tmp2 * s * f1 / f2 + 1
                          tmp1 <- tmp1 + tmp2 ^ 2 * s_t * f2 / f3
                          f1 <- f2
                          f2 <- f3
                        }
                      }
                      mu_d <- e_mu_x - mu_x
                      t_mu_d <- t(mu_d)
                      q <-
                        t_mu_d %*% inv_Sigma_x %*% mu_d * n + (e_mu_dy - mu_dy - t_mu_d %*% Sigma_xdy) ^
                        2 * tmp1 * inv_f1[m + 1, m + 1]
                      return(c("Accuracy of estimates" = as.numeric((
                        l + t + q - (m + 1) * n
                      ) / n)))
                    }
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  warning("The input is not valid. Consider checking it against the requirements of the function.")
  return(NULL)
}