#' @title estimate_parameters
#'
#' @description Calculates maximum likelihood estimates of the statistical measures of the output-differenced version of the one-dimensional cumulative structural equation model with shock-error output measurement equation and assumptions of normality and independence.
#'
#' @param dat An (n + 1) x (m + 1) data frame of finite numeric elements (possibly except for row 1 columns 1 to m) containing observed input (columns 1 to m) and output (column m + 1) data of the original model.
#' @param tol A tolerance parameter of the golden section search algorithm used for minimizing the one-dimensional likelihood function (vector of length 1, finite positive numeric element).
#'
#' @return A list consisting of 3 elements: 1) estimate of the covariance at lag 0 of the data that result from the output-differenced model (Sigma; (m + 1) x (m + 1) matrix of numeric elements); 2) estimate of the only non-zero element of the negative covariance at lag 1 of the data that result from the output-differenced model (sigma_y^2; vector of length 1, numeric element); 3) estimate of the mean of the data that result from the output-differenced model (mu; (m + 1) x 1 matrix of numeric elements).
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
#' estimate_parameters(data, 0.00001)
#'
#' @export

estimate_parameters <- function(dat, tol) {
  if (is.data.frame(dat) && is.vector(tol)) {
    m <- ncol(dat) - 1
    if (nrow(dat) >= 2 && m >= 1 && length(tol) == 1) {
      if (all(sapply(dat, is.numeric)) && is.numeric(tol)) {
        if (all(sapply(dat[-1, 1:m], is.finite)) &&
            all(is.finite(dat[, m + 1])) && is.finite(tol)) {
          if (tol > 0) {
            x <- as.matrix(dat)[-1, 1:m]
            dy <- diff(as.matrix(dat[, m + 1]))
            colnames(dy) <- "dy"
            
            n <- nrow(x)
            o <- matrix(1, n, 1)
            opt <- minimize_likelihood(n, x, dy, o, tol)
            
            Sigma_xdy <- opt[[2]] %*% opt[[4]]
            Sigma_dy <-
              opt[[6]] * (opt[[3]] ^ 2 + 1) + t(Sigma_xdy) %*% opt[[4]]
            
            Sigma <-
              rbind(cbind(opt[[2]], Sigma_xdy), cbind(t(Sigma_xdy), Sigma_dy))
            sigma_y_squared <-
              c("sigma_y_squared" = as.numeric(opt[[3]] * opt[[6]]))
            mu <- rbind(opt[[1]], opt[[5]])
            rownames(mu)[m + 1] <- "dy"
            colnames(mu) <- "mu"
            
            estimates <- list(Sigma, sigma_y_squared, mu)
            names(estimates) <-
              c(
                "Estimated covariance matrix at lag 0",
                "Estimated non-zero element of the negative covariance matrix at lag 1",
                "Estimated mean vector"
              )
            
            return(estimates)
          }
        }
      }
    }
  }
  warning("The input is not valid. Consider checking it against the requirements of the function.")
  return(NULL)
}