#' @title identify_model
#'
#' @description Calculates maximum likelihood estimates of the statistical parameters of the one-dimensional cumulative structural equation model with normality assumptions.
#'
#' @param x Observed input data (n x m matrix of finite numeric elements).
#' @param dy Observed output differences data (n x 1 matrix of finite numeric elements).
#' @param tol A tolerance parameter of the golden section search algorithm used for minimizing the one-dimensional likelihood function (vector of length 1, finite positive numeric element).
#'
#' @return A list consisting of 3 elements: 1) estimate of the covariance of cbind(x, dy) at lag 0 (Sigma; (m + 1) x (m + 1) matrix of numeric elements); 2) estimate of the only non-zero element of the negative covariance of cbind(x, dy) at lag 1 (sigma_y^2; vector of length 1, numeric element); 3) estimate of the mean of cbind(x, dy) (mu; (m + 1) x 1 matrix of numeric elements).
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
#' data <- generate_data(100, L, sigma, mu)
#'
#' identify_model(data[[1]], data[[2]], 0.00001)
#'
#' @export

identify_model <- function(x, dy, tol) {
  if (is.matrix(x) && is.matrix(dy) && is.vector(tol)) {
    if (dim(x)[1] == dim(dy)[1] &&
        dim(dy)[2] == 1 && length(tol) == 1) {
      if (is.numeric(x) && is.numeric(dy) && is.numeric(tol)) {
        if (all(is.finite(x)) && all(is.finite(dy)) && is.finite(tol)) {
          if (tol > 0) {
            n <- nrow(x)
            o <- matrix(1, n, 1)
            opt <- minimize_likelihood(n, x, dy, o, tol)

            Sigma_xdy <- opt[[2]] %*% opt[[4]]
            Sigma_dy <-
              opt[[6]] * (opt[[3]] ^ 2 + 1) + t(Sigma_xdy) %*% opt[[4]]

            Sigma <-
              rbind(cbind(opt[[2]], Sigma_xdy), cbind(t(Sigma_xdy), Sigma_dy))
            sigma_y_squared <- opt[[3]] * opt[[6]]
            mu <- rbind(opt[[1]], opt[[5]])

            return(list(Sigma, sigma_y_squared, mu))

          }

        }

      }

    }

  }

  warning("The input is not valid. Consider checking it against the requirements of the function.")
  return(NULL)

}
