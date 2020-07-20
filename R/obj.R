#' Maximum Likelihood Objective Function
#'
#' Calculates the maximum likelihood objective function.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family objective functions
#' @keywords optimization
#' @param Sigmatheta Model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma}\left( \boldsymbol{\theta} \right)}).
#' @param Sigma Variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma}}).
#' @import jeksterslabRmatrix
#' @export
ml <- function(Sigmatheta,
               Sigma) {
  log(det(Sigmatheta)) + tr(Sigma %*% solve(Sigmatheta)) - log(det(Sigma)) - nrow(Sigma)
}

#' Generalized Least Squares Objective Function
#'
#' Calculates the generalized least squares objective function.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family objective functions
#' @keywords optimization
#' @inheritParams ml
#' @export
gls <- function(Sigmatheta,
                Sigma) {
  (1 / 2) %*% tr((diag(nrow(Sigma)) - Sigmatheta %*% solve(Sigma))^2)
}

#' Unweighted Least Squares Objective Function
#'
#' Calculates the unweighted least squares objective function.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family objective functions
#' @keywords optimization
#' @inheritParams ml
#' @export
uls <- function(Sigmatheta,
                Sigma) {
  (1 / 2) %*% tr((Sigma - Sigmatheta))^2
}
