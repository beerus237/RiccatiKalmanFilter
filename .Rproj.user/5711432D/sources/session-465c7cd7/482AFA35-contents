#' Solve Riccati Equation Iteratively
#'
#' This function solves the Riccati equation iteratively.
#'
#' @param A A numeric matrix. The state transition matrix.
#' @param B A numeric matrix. The input matrix.
#' @param Q A numeric matrix. The state cost matrix.
#' @param R A numeric matrix. The input cost matrix.
#' @param tol A numeric value. The tolerance for convergence. Default is 1e-8.
#' @param max_iter An integer. The maximum number of iterations. Default is 100000.
#' @param alpha A numeric value. The step size for the iterative update. Default is 1.
#' @return A numeric matrix. The solution matrix P for the Riccati equation.
#' @details
#' This function uses an iterative method to solve the discrete-time Riccati equation:
#' \deqn{P_{k+1} = A P_k A' - A P_k B (R + B' P_k B)^{-1} B' P_k A + Q}
#' The iteration stops when the change in P between iterations is smaller than the specified tolerance \code{tol},
#' or when the maximum number of iterations \code{max_iter} is reached.
#' @examples
#' A <- matrix(c(1, 0.1, 0, 1), nrow = 2, ncol = 2)
#' B <- matrix(c(0.1, 0.2), nrow = 2, ncol = 1)
#' Q <- matrix(c(0.1, 0, 0, 0.1), nrow = 2, ncol = 2)
#' R <- matrix(0.01, nrow = 1, ncol = 1)
#' solve_riccati_iterative(A, B, Q, R)
solve_riccati_iterative <- function(A, B, Q, R, tol = 1e-8, max_iter = 100000, alpha = 1) {
  n <- nrow(A)
  P <- matrix(0.01 * runif(n^2), n, n) # Random initialization of matrix P

  for (iter in 1:max_iter) {
    # Calculate the temporary matrix
    temp <- A %*% P + P %*% A - P %*% B %*% solve(R) %*% t(B) + Q
    # Update the matrix P using the iterative approach
    P_new <- P + alpha * solve(temp, Q)

    # Check the convergence criteria
    if (norm(P_new - P, "F") / norm(P_new, "F") < tol) {
      break
    }
    P <- P_new
  }
  return(P)
}

