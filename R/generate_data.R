#' @title generate_data
#'
#' @description Generates data according to the one-dimensional cumulative structural equation model with shock-error output measurement equation and assumptions of normality and independence with given model parameter values.
#'
#' @param n The number of time moments to generate the data for (vector of length  1, finite positive integer).
#' @param L Factor loadings ((m + 1) x k matrix of finite numeric elements: the first m rows correspond to the input measurement equation; the last row corresponds to the transition equation).
#' @param sigma Standard deviations of the error/noise terms ((m + 2) x 1 matrix of finite non-negative numeric elements: the first m rows correspond to the input measurement equation; the row before the last one corresponds to the transition equation; the last row corresponds to the output measurement equation).
#' @param mu Intercept terms ((m + 1) x 1 matrix of finite numeric elements; the first m rows correspond to the input measurement equation; the last row corresponds to the transition equation).
#'
#' @return An (n + 1) x (m + 1) data frame of numeric elements (except for row 1 columns 1 to m that contain NA's) containing observed input (columns 1 to m) and output (column m + 1) data.
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
        dim(L)[1] == dim(mu)[1] &&
        dim(sigma)[2] == 1 && dim(mu)[2] == 1) {
      if (is.numeric(n) &&
          is.numeric(L) && is.numeric(sigma) && is.numeric(mu)) {
        if (is.finite(n) &&
            all(is.finite(L)) &&
            all(is.finite(sigma)) && all(is.finite(mu))) {
          if (n == as.integer(n) && n >= 1 && all(sigma >= 0)) {
            m <- nrow(L) - 1
            k <- ncol(L)
            
            xi <- matrix(rnorm(n * k, 0, 1), nrow = n)
            deltazeta <-
              apply(as.matrix(sigma[1:(m + 1), ]), 1, function(x)
                rnorm(n, mean = 0, sd = x))
            
            xdeltaeta <-
              matrix(1, nrow = n) %*% t(mu) + xi %*% t(L) + deltazeta
            
            x <- xdeltaeta[, 1:m]
            x0 <- matrix(NA, nrow = n + 1, ncol = m)
            x0[-1, ] <- x
            
            deltaeta <- xdeltaeta[, m + 1]
            eta <- append(cumsum(deltaeta), 0, 0)
            
            epsilon <-
              matrix(rnorm(n + 1, 0, sigma[m + 2,]), nrow = n + 1)
            y <- eta + epsilon
            
            df <- as.data.frame(cbind(x0, y))
            colnames(df) <- c(paste0("x", 1:m), "y")
            rownames(df) <- c(0:n)
            
            return(df)
          }
        }
      }
    }
  }
  warning("The input is not valid. Consider checking it against the requirements of the function.")
  return(NULL)
}