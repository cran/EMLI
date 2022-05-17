#' @title generate_data
#'
#' @description Generates data according to the one-dimensional cumulative structural equation model with normality assumptions with given model parameter values.
#'
#' @param n The number of time moments to generate the data for (vector of length  1, finite positive integer).
#' @param L Factor loadings ((m + 1) x k matrix of finite numeric elements: the first m rows correspond to the input measurement equation; the last row corresponds to the transition equation).
#' @param sigma Standard deviations of the error/noise terms ((m + 2) x 1 matrix of finite non-negative numeric elements: the first m rows correspond to the input measurement equation; the row before the last one corresponds to the transition equation; the last row corresponds to the output measurement equation).
#' @param mu Intercept terms ((m + 1) x 1 matrix of finite numeric elements; the first m rows correspond to the input measurement equation; the last row corresponds to the transition equation).
#'
#' @return A list consisting of 2 elements: 1) observed input data (n x m matrix of numeric elements); 2) observed output differences data (n x 1 matrix of numeric elements).
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

#' generate_data(10, L, sigma, mu)
#'
#' @export
#'
#' @import stats

generate_data <- function(n, L, sigma, mu) {
  if (is.vector(n) &&
      is.matrix(L) && is.matrix(sigma) && is.matrix(mu)) {
    if (length(n) == 1 &&
        dim(L)[1] == dim(sigma)[1] - 1 &&
        dim(L)[1] == dim(mu)[1] && dim(sigma)[2] == 1 && dim(mu)[2] == 1) {
      if (is.numeric(n) &&
          is.numeric(L) && is.numeric(sigma) && is.numeric(mu)) {
        if (is.finite(n) &&
            all(is.finite(L)) && all(is.finite(sigma)) && all(is.finite(mu))) {
          if (n == as.integer(n) && n >= 1 && all(sigma >= 0)) {
            m <- nrow(L) - 1
            k <- ncol(L)

            xi <- matrix(rnorm(n * k, 0, 1), nrow = n)

            deltazeta <- matrix(NA, nrow = n, ncol = m + 1)
            for (j in 1:(m + 1)) {
              deltazeta[, j] <- rnorm(n, 0, sigma[j,])

            }

            epsilon <-
              matrix(rnorm(n + 1, 0, sigma[m + 2,]), nrow = n + 1)

            x <- matrix(NA, nrow = n, ncol = m)
            for (j in 1:m) {
              for (i in 1:n) {
                x[i, j] <- mu[j,] + L[j,] %*% xi[i,]  + deltazeta[i, j]

              }

            }

            dy <- matrix(NA, nrow = n, ncol = 1)
            for (i in 1:n) {
              dy[i, ] <-
                mu[m + 1,] + L[m + 1,] %*% xi[i,]  + deltazeta[i, m + 1] + (epsilon[i + 1,] - epsilon[i,])

            }

            return(list(x, dy))

          }
        }

      }

    }

  }

  warning("The input is not valid. Consider checking it against the requirements of the function.")
  return(NULL)

}
