#' Reticular Action Model - Model-Implied Variance Covariance Matrix
#' \eqn{\left( \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right) \right)}
#'
#' Derives the model-implied variance-covariance matrix
#' \eqn{\left( \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right) \right)}
#' using the Reticular Action Model notation.
#'
#' The model-implied variance-covariance matrix
#' \eqn{\left( \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right) \right)}
#' as a function of Reticular Action Model matrices
#' is given by
#'   \deqn{
#'     \boldsymbol{\Sigma}
#'     \left(
#'       \boldsymbol{\theta}
#'     \right)
#'     =
#'     \mathbf{F}
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'     \mathbf{S}
#'     \left[
#'       \left(
#'         \mathbf{I} - \mathbf{A}
#'       \right)^{-1}
#'     \right]^{\prime}
#'     \mathbf{F}^{\prime}
#'   }
#' where
#' \eqn{\mathbf{A}}
#' represents asymmetric paths (single-headed arrows),
#' such as regression coefficients and factor loadings.
#' \eqn{\mathbf{S}}
#' represents symmetric paths (double-headed arrows),
#' such as variances and covariances.
#' \eqn{\mathbf{F}}
#' represents the filter matrix
#' used to select the observed variables.
#' \eqn{\mathbf{I}}
#' represents the identity matrix.
#'
#' When all the variables
#' in the \eqn{\mathbf{A}} and \eqn{\mathbf{S}} matrices
#' are observed,
#'   \deqn{
#'     \boldsymbol{\Sigma}
#'     \left(
#'       \boldsymbol{\theta}
#'     \right)
#'     =
#'     \mathbf{F}
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'     \mathbf{S}
#'     \left[
#'       \left(
#'         \mathbf{I} - \mathbf{A}
#'       \right)^{-1}
#'     \right]^{\prime}
#'     \mathbf{F}^{\prime}
#'     =
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'     \mathbf{S}
#'     \left[
#'       \left(
#'         \mathbf{I} - \mathbf{A}
#'       \right)^{-1}
#'     \right]^{\prime} , \\
#'     \text{ when }
#'     \mathbf{F} = \mathbf{I} .
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param A Matrix.
#' Asymmetric paths (single-headed arrows),
#' such as regression coefficients and factor loadings.
#' @param S Matrix.
#' Symmetric paths (double-headed arrows),
#' such as variances and covariances.
#' @param F Matrix.
#' Filter matrix
#' used to select the observed variables.
#' @param I Matrix.
#' Identity matrix.
#' @param filter Logical.
#' If `TRUE`, the `F` matrix is used.
#' If `FALSE`, the `F` matrix is omitted and the argument `F` is ignored.
#' See `Value` and `Details` below for more information.
#' @return Returns the model-implied variance-covariance matrix
#' \eqn{\left( \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right) \right)}
#' derived from the \eqn{\mathbf{A}},
#' \eqn{\mathbf{S}},
#' \eqn{\mathbf{F}}, and
#' \eqn{\mathbf{I}} matrices.
#' If `filter = TRUE`, the `F` matrix is used as usual.
#'   \deqn{
#'     \mathbf{F}
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'     \mathbf{S}
#'     \left[
#'       \left(
#'         \mathbf{I} - \mathbf{A}
#'       \right)^{-1}
#'     \right]^{\prime}
#'     \mathbf{F}^{\prime} .
#'   }
#' If `filter = FALSE`, the `F` matrix is omitted.
#' \deqn{
#'   \left(
#'     \mathbf{I} - \mathbf{A}
#'   \right)^{-1}
#'   \mathbf{S}
#'   \left[
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'   \right]^{\prime}
#' } .
#' @references
#' McArdle, J. J. (2013).
#' The development of the RAM rules for latent variable structural equation modeling.
#' In A. Maydeu-Olivares & J. J. McArdle (Eds.),
#' *Contemporary Psychometrics: A festschrift for Roderick P. McDonald* (pp. 225--273).
#' Lawrence Erlbaum Associates.
#'
#' McArdle, J. J., & McDonald, R. P. (1984).
#' Some algebraic properties of the Reticular Action Model for moment structures.
#' *British Journal of Mathematical and Statistical Psychology*, *37* (2), 234--251.
#' @family SEM notation functions
#' @keywords matrix ram
#' @examples
#' #################################################
#' # Simple mediation model with observed variables.
#' #################################################
#' # parameters
#' alpha <- 0.5
#' beta <- 0.4
#' tauprime <- 0.25
#' sigma2X <- 1.00
#' sigma2epsilonM <- 0.75
#' sigma2epsilonY <- 0.6775
#' # A and S matrices only contain observed variables X, M, and Y
#' A <- matrix(data = 0, nrow = 3, ncol = 3)
#' A[2, 1] <- alpha
#' A[3, 1] <- beta
#' A[3, 2] <- tauprime
#' S <- diag(c(sigma2X, sigma2epsilonM, sigma2epsilonY))
#' F <- I <- diag(3)
#' ram_Sigmatheta(A = A, S = S, F = F, I = I)
#' # Equivalent
#' ram_Sigmatheta(A = A, S = S, I = I, filter = FALSE)
#' # A and S matrices also contain latent variables sigma2epsilonM and sigma2epsilonY
#' A <- matrix(data = 0, nrow = 5, ncol = 5)
#' A[2, 1] <- alpha
#' A[2, 4] <- 1
#' A[3, 1] <- beta
#' A[3, 2] <- tauprime
#' A[3, 5] <- 1
#' S <- diag(c(sigma2X, 0, 0, sigma2epsilonM, sigma2epsilonY))
#' I <- diag(nrow(A))
#' F <- matrix(data = 0, nrow = 3, ncol = 5)
#' F[1, 1] <- 1
#' F[2, 2] <- 1
#' F[3, 3] <- 1
#' # In this case since F is not equal to I, we CANNOT use filter = FALSE
#' ram_Sigmatheta(A = A, S = S, F = F, I = I)
#' #######################
#' # One-factor CFA model.
#' #######################
#' A <- matrix(data = 0, ncol = 6, nrow = 6)
#' for (i in 2:6) {
#'   A[i, 1] <- 0.5
#' }
#' S <- diag(c(1, 0.75, 0.75, 0.75, 0.75, 0.75))
#' F <- I <- diag(nrow(A))
#' F <- diag(nrow(A) - 1)
#' F <- cbind(0, F)
#' ram_Sigmatheta(
#'   A = A,
#'   S = S,
#'   F = F,
#'   I = I
#' )
#' @export
ram_Sigmatheta <- function(A,
                           S,
                           F,
                           I,
                           filter = TRUE) {
  #  x <- solve(I - A)
  #  F %*% x %*% S %*% t(x) %*% t(F)
  inverse <- solve(I - A)
  full <- inverse %*% S %*% t(inverse)
  if (filter) {
    return(F %*% full %*% t(F))
  } else {
    return(full)
  }
}

#' Reticular Action Model - Model-Implied Mean Vector
#' \eqn{\left( \boldsymbol{\mu} \right)}
#'
#' Derives the model-implied
#' \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)} vector
#' using the Reticular Action Model notation.
#'
#' \deqn{
#'   \boldsymbol{\mu}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   =
#'   \mathbf{F}
#'   \left(
#'     \mathbf{I}
#'     -
#'     \mathbf{A}
#'   \right)^{-1}
#'   \mathbf{M}
#' }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams ram_Sigmatheta
#' @inherit ram_Sigmatheta references
#' @param M Numeric vector.
#' Mean structure.
#' Vector of means and intercepts.
#' @return Returns the means
#' \eqn{\left( \boldsymbol{\mu} \left( \boldsymbol{\theta} \right) \right)}
#' derived from the \eqn{\mathbf{A}},
#' \eqn{\mathbf{F}},
#' \eqn{\mathbf{I}}, matrices and
#' \eqn{\mathbf{M}} vector.
#' @examples
#' A <- matrix(
#'   data = c(
#'     0, 0.26^(1 / 2), 0,
#'     0, 0, 0.26^(1 / 2),
#'     0, 0, 0
#'   ),
#'   ncol = 3
#' )
#' F <- I <- diag(3)
#' M <- c(100, 49.0098049, 49.0098049)
#' ram_mutheta(A = A, F = F, I = I, M = M)
#' @export
ram_mutheta <- function(A,
                        F,
                        I,
                        M) {
  F %*% solve(I - A) %*% M
}

#' Reticular Action Model - Mean Structure Vector
#' \eqn{\left( \mathbf{M} \right)}
#'
#' Mean structure vector
#' \eqn{\left( \mathbf{M} \right)}
#' using the Reticular Action Model notation.
#'
#' \deqn{\mathbf{M}
#'   =
#'   \mathbf{F}
#'   \left(
#'     \mathbf{I}
#'     -
#'     \mathbf{A}
#'   \right)^{-1}
#'   \boldsymbol{\mu}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#' }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams ram_Sigmatheta
#' @inherit ram_Sigmatheta references
#' @param mu Numeric vector.
#' Vector of means
#' \eqn{\left( \boldsymbol{\mu} \left( \boldsymbol{\theta} \right) \right)}.
#' @return Returns the mean structure vector \eqn{\left( \mathbf{M} \right)}
#' derived from the \eqn{\mathbf{A}},
#' \eqn{\mathbf{F}},
#' \eqn{\mathbf{I}}, matrices and
#' \eqn{\left( \boldsymbol{\mu} \left( \boldsymbol{\theta} \right) \right)} vector.
#' @examples
#' A <- matrix(
#'   data = c(
#'     0, 0.26^(1 / 2), 0,
#'     0, 0, 0.26^(1 / 2),
#'     0, 0, 0
#'   ),
#'   ncol = 3
#' )
#' F <- I <- diag(3)
#' mu <- c(100, 100, 100)
#' ram_M(A = A, F = F, I = I, mu = mu)
#' @export
ram_M <- function(A,
                  F,
                  I,
                  mu) {
  F %*% (I - A) %*% mu
}

#' Reticular Action Model - Model-Implied Variance Covariance Matrix or S Matrix
#' \eqn{\left( \boldsymbol{\Sigma}} or \eqn{\mathbf{S} \right)}
#'
#' Model-implied variance-covariance matrix
#' \eqn{\left( \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right) \right)}
#' or \eqn{\mathbf{S}} Matrix
#' using the Reticular Action Model notation.
#'
#' Derives the
#' model-implied variance-covariance matrix
#' \eqn{\left( \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right) \right)}
#' or the
#' \eqn{\mathbf{S}}
#' matrix
#' from the
#' \eqn{\mathbf{A}}
#' matrix
#' and sigma squared
#' \eqn{\left( \sigma^2 \right)}
#' vector (variances).
#' **Note that the first element in the A ans S matrices should be an exogenous variable.**
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams ram_Sigmatheta
#' @inherit ram_Sigmatheta references
#' @param sigma2 Numeric vector.
#' Vector of variances \eqn{\left( \sigma^2 \right)}.
#' **The first element should be the variance of an exogenous variable**
#' @param SigmaMatrix Logical.
#' If `TRUE`,
#' returns
#' the model-implied variance-covariance matrix
#' \eqn{\left( \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right) \right)}.
#' If `FALSE`,
#' returns the
#' \eqn{\mathbf{S}}
#' matrix.
#' @return Returns
#' the model-implied variance-covariance matrix
#' \eqn{\left( \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right) \right)}
#' or the
#' \eqn{\mathbf{S}}
#' matrix
#' derived from the
#' \eqn{\mathbf{A}}
#' matrix and
#' \eqn{\sigma^2}
#' vector.
#' @examples
#' A <- matrix(
#'   data = c(
#'     0, 0.26^(1 / 2), 0,
#'     0, 0, 0.26^(1 / 2),
#'     0, 0, 0
#'   ),
#'   ncol = 3
#' )
#' sigma2 <- c(15^2, 15^2, 15^2)
#' F <- I <- diag(3)
#' # Returns the model-implied variance-covariance matrix
#' ram_S(A = A, sigma2 = sigma2, F = F, I = I, SigmaMatrix = TRUE)
#' # Returns the model-implied S matrix
#' ram_S(A = A, sigma2 = sigma2, F = F, I = I, SigmaMatrix = FALSE)
#' @export
ram_S <- function(A,
                  sigma2,
                  F,
                  I,
                  SigmaMatrix = TRUE) {
  #  S <- matrix(
  #    data = 0,
  #    ncol = dim(A)[1],
  #    nrow = dim(A)[2]
  #  )
  #  Sigmatheta <- ram_Sigmatheta(
  #    A = A,
  #    S = S,
  #    F = F,
  #    I = I
  #  )
  #  for (i in 1:nrow(A)) {
  #    S[i, i] <- sigma2[i] - Sigmatheta[i, i]
  #    Sigmatheta <- ram_Sigmatheta(
  #      A = A,
  #      S = S,
  #      F = F,
  #      I = I
  #    )
  #  }
  #  if (SigmaMatrix) {
  #    return(Sigmatheta)
  #  } else {
  #    return(S)
  #  }
  S <- matrix(
    data = 0,
    ncol = dim(A)[1],
    nrow = dim(A)[2]
  )
  Sigmatheta_temp <- ram_Sigmatheta(
    A = A,
    S = S,
    F = F,
    I = I,
    filter = FALSE
  )
  for (i in 1:nrow(A)) {
    S[i, i] <- sigma2[i] - Sigmatheta_temp[i, i]
    # filter = FALSE makes it flexible to handle latent variables
    Sigmatheta_temp <- ram_Sigmatheta(
      A = A,
      S = S,
      F = F,
      I = I,
      filter = FALSE
    )
  }
  if (SigmaMatrix) {
    # filter = TRUE filter out the latent variables
    Sigmatheta <- ram_Sigmatheta(
      A = A,
      S = S,
      F = F,
      I = I,
      filter = TRUE
    )
    return(Sigmatheta)
  } else {
    return(S)
  }
}

#' Reticular Action Model - Residuals
#' \eqn{\left( \hat{\Sigma} - \hat{\Sigma} \left( \hat{\theta} \right) \right)}
#'
#' \deqn{
#'   \hat{\Sigma}
#'   -
#'   \hat{\Sigma}
#'   \left(
#'     \hat{\theta}
#'   \right)
#' }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams ram_Sigmatheta
#' @param Sigmahat Matrix.
#' Estimated variance-covariance matrix
#' \eqn{\left( \hat{\boldsymbol{\Sigma}} \right)}
#' @param Sigmathetahat Matrix.
#' Model-implied variance-covariance matrix
#' as a function of estimated parameters
#' \eqn{\left( \boldsymbol{\Sigma}\left( \hat{\boldsymbol{\theta}} \right) \right)}.
#' @param Ahat Matrix.
#' Estimated asymmetric paths (single-headed arrows),
#' such as regression coefficients and factor loadings.
#' @param Shat Matrix.
#' Estimated symmetric paths (double-headed arrows),
#' representing variances and covariances.
#' @export
ram_e <- function(Sigmahat,
                  Sigmathetahat = NULL,
                  Ahat,
                  Shat,
                  F,
                  I) {
  if (is.null(Sigmathetahat)) {
    Sigmathetahat <- ram_Sigmatheta(
      A = Ahat,
      S = Shat,
      F = F,
      I = I,
      filter = TRUE
    )
  }
  Sigmahat - Sigmathetahat
}
