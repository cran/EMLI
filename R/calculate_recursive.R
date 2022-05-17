#' @title calculate_recursive
#'
#' @description Internal function for recursively calculating the terms that involve inversion of a matrix whose dimensions depend on the sample size.
#'
#' @noRd

calculate_recursive <- function(n, u, v, s) {
  a <- u[1,]
  b <- v[1,]
  c <- a %*% t(b)
  f1 <- 1
  f2 <- f1 + s ^ 2

  if (n > 1) {
    for (i in 2:n) {
      f3 <- f2 + s ^ (2 * i)

      ab_ <- s * (f1 / f2)
      c_ <- (f1 * f3) / f2 ^ 2

      f1 <- f2
      f2 <- f3

      a <- a * ab_  + u[i,]
      b <- b * ab_ + v[i,]
      c <- c * (1 - 1 / i) * c_ + (a %*% t(b)) / i


    }

  }

  if (sum(dim(c)) == 2) {
    return(list(c[1, 1], f1, f2))

  }

  else {
    return(list(c, f1, f2))

  }

}
