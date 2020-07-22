#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Reticular Action Model - Model-Implied Variance-Covariance Matrix
#'   \eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'
#' @description Derives the model-implied variance-covariance matrix
#'   \eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The model-implied variance-covariance matrix
#'   \eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'   as a function of Reticular Action Model (RAM) matrices
#'   is given by
#'
#'   \deqn{
#'     \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)
#'     =
#'     \mathbf{F} \left( \mathbf{I} - \mathbf{A} \right)^{-1} \mathbf{S}
#'     \left[ \left( \mathbf{I} - \mathbf{A} \right)^{-1} \right]^{T} \mathbf{F}^{T}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{m \times m}} represents asymmetric paths (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{S}_{m \times m}} represents symmetric paths (double-headed arrows),
#'     such as variances and covariances,
#'   - \eqn{\mathbf{F}_{k \times m}} represents the filter matrix
#'     used to select the observed variables,
#'   - \eqn{\mathbf{I}_{m \times m}} represents an identity matrix,
#'   - \eqn{k} number of observed variables,
#'   - \eqn{q} number of latent variables, and
#'   - \eqn{m} number of observed and latent variables, that is \eqn{k + q} .
#'
#' @family SEM notation functions
#' @keywords matrix ram
#' @param A `m x m` numeric matrix
#'   \eqn{\mathbf{A}_{m \times m}}.
#'   Asymmetric paths (single-headed arrows),
#'   such as regression coefficients and factor loadings.
#' @param S `m x m` numeric matrix
#'   \eqn{\mathbf{S}_{m \times m}}.
#'   Symmetric paths (double-headed arrows),
#'   such as variances and covariances.
#' @param filter `k x m` numeric matrix
#'   \eqn{\mathbf{F}_{k \times m}}.
#'   Filter matrix used to select variables.
#' @return Returns the model-implied variance-covariance matrix
#'   \eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'   derived from the `A`, `S`, and `filter` matrices.
#' @references
#'   McArdle, J. J. (2013).
#'   The development of the RAM rules for latent variable structural equation modeling.
#'   In A. Maydeu-Olivares & J. J. McArdle (Eds.),
#'   *Contemporary Psychometrics: A festschrift for Roderick P. McDonald* (pp. 225--273).
#'   Lawrence Erlbaum Associates.
#'
#'   McArdle, J. J., & McDonald, R. P. (1984).
#'   Some algebraic properties of the Reticular Action Model for moment structures.
#'   *British Journal of Mathematical and Statistical Psychology*, *37* (2), 234--251.
#' @export
ramSigmatheta <- function(A,
                          S,
                          filter) {
  invIminusA <- solve(diag(nrow(A)) - A)
  filter %*% invIminusA %*% S %*% t(invIminusA) %*% t(filter)
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Reticular Action Model - Model-Implied Mean Vector
#'   \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'
#' @description Derives the model-implied mean vector
#'   \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The model-implied mean vector
#'   \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'   as a function of Reticular Action Model (RAM) matrices
#'   is given by
#'
#'   \deqn{
#'     \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)
#'     =
#'     \mathbf{F} \left( \mathbf{I} - \mathbf{A} \right)^{-1} \mathbf{M}
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{m \times m}} represents asymmetric paths (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{I}_{m \times m}} represents an identity matrix,
#'   - \eqn{\mathbf{F}_{k \times m}} represents the filter matrix
#'     used to select the observed variables,
#'   - \eqn{\mathbf{M}_{m \times 1}} represents the mean structure,
#'     that is, a vector of means and intercepts,
#'   - \eqn{k} number of observed variables,
#'   - \eqn{q} number of latent variables, and
#'   - \eqn{m} number of observed and latent variables, that is \eqn{k + q} .
#'
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams ramSigmatheta
#' @inherit ramSigmatheta references
#' @param M `m x 1` numeric vector \eqn{\mathbf{M}_{m \times 1}}.
#'   Mean structure. Vector of means and intercepts.
#' @return Returns the model-implied mean vector
#'   \eqn{ \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'   derived from the \eqn{\mathbf{A}}, \eqn{\mathbf{F}}, \eqn{\mathbf{I}},
#'   matrices and \eqn{\mathbf{M}} vector.
#' @export
rammutheta <- function(M,
                       A,
                       filter) {
  filter %*% solve(diag(nrow(A)) - A) %*% M
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Reticular Action Model - Mean Structure Vector \eqn{\mathbf{M}}
#'
#' @description Derives the mean structure vector \eqn{\mathbf{M}}
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The mean structure vector \eqn{\mathbf{M}}
#'   as a function of Reticular Action Model (RAM) matrices is given by
#'
#'   \deqn{
#'     \mathbf{M}
#'     =
#'     \left( \mathbf{I} - \mathbf{A} \right)^{-1}
#'     \mathbf{F}^{T} \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)
#'   }
#'
#'   where
#'
#'   - \eqn{\mathbf{A}_{m \times m}} represents asymmetric paths (single-headed arrows),
#'     such as regression coefficients and factor loadings,
#'   - \eqn{\mathbf{I}_{m \times m}} represents an identity matrix,
#'   - \eqn{\mathbf{F}_{k \times m}} represents the filter matrix
#'     used to select the observed variables,
#'   - \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'     is the \eqn{k \times 1} model-implied mean vector
#'   - \eqn{k} number of observed variables,
#'   - \eqn{q} number of latent variables, and
#'   - \eqn{m} number of observed and latent variables, that is \eqn{k + q} .
#'
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams ramSigmatheta
#' @inherit ramSigmatheta references
#' @param mu `k x 1` numeric vector
#'   \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)_{k \times 1}} .
#'   Model-implied meam vector.
#' @return Returns the mean structure vector \eqn{\mathbf{M}}
#'   derived from the \eqn{\mathbf{A}}, \eqn{\mathbf{F}},
#'   \eqn{ \mathbf{I}}, matrices and
#'   \eqn{ \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'   vector.
#' @export
ramM <- function(mu,
                 A,
                 filter) {
  solve(diag(nrow(A)) - A) %*% t(filter) %*% mu
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Reticular Action Model - The \eqn{\mathbf{S}} Matrix from \eqn{\sigma^2}
#'
#' @description Derives the \eqn{\mathbf{S}} matrix
#'   using the Reticular Action Model (RAM) notation
#'   from variable variances \eqn{\sigma^2}.
#'   The off-diagonal elements of the \eqn{\mathbf{S}} matrix are assumed to be zeroes.
#'
#' @details The \eqn{\mathbf{S}} matrix is derived using the \eqn{\mathbf{A}} matrix
#'   and sigma squared \eqn{\left( \sigma^2 \right)} vector (variances).
#'   **Note that the first or last (see `start` argument) element
#'   in the `A` and `S` matrices should be an exogenous variable.**
#'
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams ramSigmatheta
#' @inherit ramSigmatheta references
#' @param sigma2 Numeric vector.
#'   Vector of variances \eqn{\sigma^2}.
#' @param start Logical.
#'   If `TRUE`, an exogenous variable is positioned as the first element in the matrices.
#'   If `FALSE`, an exogenous variable is positioned as the last element in the matrices.
#' @return Returns the  \eqn{\mathbf{S}} matrix.
#' @export
ramsigma2 <- function(sigma2,
                      A,
                      start = TRUE) {
  message(
    "The off-diagonal elements of the S matrix are assumed to be zeroes."
  )
  S <- matrix(
    data = 0,
    ncol = dim(A)[1],
    nrow = dim(A)[2]
  )
  # filter as identity matrix to retain all variables
  filter <- diag(nrow(A))
  Sigmatheta_temp <- ramSigmatheta(
    A = A,
    S = S,
    filter = filter
  )
  if (start) {
    index <- 1:nrow(A)
  } else {
    index <- nrow(A):1
  }
  for (i in index) {
    S[i, i] <- sigma2[i] - Sigmatheta_temp[i, i]
    Sigmatheta_temp <- ramSigmatheta(
      A = A,
      S = S,
      filter = filter
    )
  }
  S
}
