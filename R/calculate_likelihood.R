#' @title calculate_likelihood
#'
#' @description Calculates the likelihood function value for given data and statistical measure values of the output-differenced version of the one-dimensional cumulative structural equation model with shock-error output measurement equation and assumptions of normality and independence. Suitable when there are no contradictions in statistical measure values.
#' 
#' @param dat An (n + 1) x (m + 1) data frame of finite numeric elements (possibly except for row 1 columns 1 to m) containing observed input (columns 1 to m) and output (column m + 1) data of the original model.
#' @param params A list consisting of 3 elements: 1) Sigma ((m + 1) x (m + 1) matrix of finite numeric elements); 2) sigma_y^2 (vector of length 1, finite numeric element); 3) mu ((m + 1) x 1 matrix of finite numeric elements).
#'
#' @return Calculated likelihood function value (vector of length 1, numeric element).
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
#' estimated_parameters <- estimate_parameters(data, 0.00001)
#'
#' calculate_likelihood(data, estimated_parameters)
#'
#' @export

calculate_likelihood <- function(dat, params) {
  if (is.list(params) && is.data.frame(dat)) {
    m <- ncol(dat) - 1
    if (length(params) == 3 && nrow(dat) >= 2 && m >= 1) {
      if (is.matrix(params[[1]]) &&
          is.vector(params[[2]]) && is.matrix(params[[3]])) {
        if (dim(params[[1]])[1] == dim(params[[1]])[2] &&
            dim(params[[1]])[1] == dim(params[[3]])[1] &&
            dim(params[[3]])[2] == 1 && length(params[[2]]) == 1) {
          if (is.numeric(params[[1]]) &&
              is.numeric(params[[2]]) &&
              is.numeric(params[[3]]) && all(sapply(dat, is.numeric))) {
            if (all(is.finite(params[[1]])) &&
                is.finite(params[[2]]) &&
                all(is.finite(params[[3]])) &&
                all(sapply(dat[-1, 1:m], is.finite)) &&
                all(is.finite(dat[, m + 1]))) {
              if (isSymmetric(params[[1]]) &&
                  min(eigen(params[[1]])$values) >= 0) {
                inv_Sigma_dy <- solve(params[[1]])[m + 1, m + 1]
                tmp <- params[[2]] * inv_Sigma_dy
                if (tmp >= 0 && tmp <= 0.5) {
                  x <- as.matrix(dat)[-1, 1:m]
                  dy <- diff(as.matrix(dat[, m + 1]))
                  colnames(dy) <- "dy"
                  
                  n <- nrow(x)
                  o <- matrix(1, n, 1)
                  
                  mu_x <- params[[3]][1:m]
                  mu_dy <- params[[3]][m + 1]
                  Sigma_x <- params[[1]][1:m, 1:m]
                  inv_Sigma_x <- solve(Sigma_x)
                  if (tmp == 0) {
                    s <- 0
                  }
                  else {
                    s <- exp(-acosh(1 / 2 * 1 / tmp))
                  }
                  Sigma_xdy <-
                    inv_Sigma_x %*% params[[1]][1:m, m + 1]
                  
                  tmp <-
                    dy - o %*% mu_dy - (x - o %*% t(mu_x)) %*% Sigma_xdy
                  rec <- calculate_recursive(n, tmp, tmp, s)
                  
                  tmp <- x - o %*% t(mu_x)
                  likelihood <-
                    log(det(Sigma_x)) + sum(diag(tmp %*% inv_Sigma_x %*% t(tmp))) / n + log((rec[[3]]) ^ (1 / n) / ((s ^ 2 + 1) * inv_Sigma_dy)) + (s ^ 2 + 1) * inv_Sigma_dy * rec[[2]] / rec[[3]] * rec[[1]]
                  likelihood <-
                    c("Likelihood function value" = as.numeric(likelihood))
                  
                  return(likelihood)
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