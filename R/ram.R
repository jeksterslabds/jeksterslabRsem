#' Reticular Action Model - Model-Implied Variance Covariance Matrix
#' \eqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-Sigmatheta)
#' }
#'
#' @description Derives the model-implied variance-covariance matrix
#' \eqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-Sigmatheta)
#' }
#' using the Reticular Action Model (RAM) notation.
#'
#' @details The model-implied variance-covariance matrix
#' \eqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-Sigmatheta)
#' }
#' as a function of Reticular Action Model (RAM) matrices
#' is given by
#' \deqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   =
#'   \mathbf{F}
#'   \left(
#'     \mathbf{I} - \mathbf{A}
#'   \right)^{-1}
#'   \mathbf{S}
#'   \left[
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'   \right]^{\prime}
#'   \mathbf{F}^{\prime}
#'   %(\#eq:sem-ram-Sigmatheta)
#' }
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
#' \deqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   =
#'   \mathbf{F}
#'   \left(
#'     \mathbf{I} - \mathbf{A}
#'   \right)^{-1}
#'   \mathbf{S}
#'   \left[
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'   \right]^{\prime}
#'   \mathbf{F}^{\prime}
#'   %(\#eq:sem-ram-Sigmatheta)
#' }
#' is equl to
#' \deqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   =
#'   \left(
#'     \mathbf{I} - \mathbf{A}
#'   \right)^{-1}
#'   \mathbf{S}
#'   \left[
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'   \right]^{\prime} ,
#'   \quad
#'   \mathrm{when}
#'   \quad
#'   \mathbf{F}
#'   =
#'   \mathbf{I} .
#'   %(\#eq:sem-ram-Sigmatheta-observed)
#' }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family SEM notation functions
#' @keywords matrix ram
#' @param A \eqn{\mathbf{A}} matrix.
#' Asymmetric paths (single-headed arrows),
#' such as regression coefficients and factor loadings.
#' @param S \eqn{\mathbf{S}} matrix.
#' Symmetric paths (double-headed arrows),
#' such as variances and covariances.
#' @param F \eqn{\mathbf{F}} matrix.
#' Filter matrix
#' used to select the observed variables.
#' @param I \eqn{\mathbf{I}} matrix.
#' Identity matrix.
#' @param filter Logical.
#' If `TRUE`, the `F` matrix is used.
#' If `FALSE`, the `F` matrix is omitted and the argument `F` is ignored.
#' See Value and Details below for more information.
#' @return Returns the model-implied variance-covariance matrix
#' \eqn{\left( \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right) \right)}
#' derived from the \eqn{\mathbf{A}},
#' \eqn{\mathbf{S}},
#' \eqn{\mathbf{F}}, and
#' \eqn{\mathbf{I}} matrices.
#' If `filter = TRUE`, the `F` matrix is used as usual.
#' \deqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   =
#'   \mathbf{F}
#'   \left(
#'     \mathbf{I} - \mathbf{A}
#'   \right)^{-1}
#'   \mathbf{S}
#'   \left[
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'   \right]^{\prime}
#'   \mathbf{F}^{\prime}
#'   %(\#eq:sem-ram-Sigmatheta)
#' }
#' If `filter = FALSE`, the `F` matrix is omitted.
#' \deqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   =
#'   \left(
#'     \mathbf{I} - \mathbf{A}
#'   \right)^{-1}
#'   \mathbf{S}
#'   \left[
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'   \right]^{\prime} .
#'   %(\#eq:sem-ram-Sigmatheta-filter-false)
#' }
#' @examples
#' # Simple mediation model with observed variables------------------------------
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
#' # One-factor CFA model--------------------------------------------------------
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
ram_Sigmatheta <- function(A,
                           S,
                           F,
                           I,
                           filter = TRUE) {
  inverse <- solve(I - A)
  full <- inverse %*% S %*% t(inverse)
  if (filter) {
    return(F %*% full %*% t(F))
  } else {
    return(full)
  }
}

#' Reticular Action Model - Model-Implied Mean Vector
#' \eqn{
#'   \boldsymbol{\mu}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-mutheta)
#' }
#'
#' @description Derives the model-implied mean vector
#' \eqn{
#'   \boldsymbol{\mu}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-mutheta)
#' }
#' using the Reticular Action Model (RAM) notation.
#'
#' @details The model-implied mean vector
#' \eqn{
#'   \boldsymbol{\mu}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-mutheta)
#' }
#' as a function of Reticular Action Model (RAM) matrices
#' is given by
#' \deqn{
#'   \boldsymbol{\mu}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   =
#'   \mathbf{F}
#'   \left(
#'     \mathbf{I} - \mathbf{A}
#'   \right)^{-1}
#'   \mathbf{M}
#'   %(\#eq:sem-ram-mutheta)
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
#' \eqn{
#'   \boldsymbol{\mu}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-mutheta)
#' }
#' derived from the
#' \eqn{\mathbf{A}},
#' \eqn{\mathbf{F}},
#' \eqn{\mathbf{I}},
#' matrices and
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
                        M,
                        filter = TRUE) {
  # F %*% solve(I - A) %*% M
  inverse <- solve(I - A)
  full <- solve(I - A) %*% M
  if (filter) {
    return(F %*% full)
  } else {
    return(full)
  }
}

#' Reticular Action Model - Mean Structure Vector
#' \eqn{
#'   \mathbf{M}
#' }
#'
#' @description Derives the mean structure vector
#' \eqn{
#'   \mathbf{M}
#' }
#' using the Reticular Action Model (RAM) notation.
#'
#' @details The mean structure vector
#' \eqn{
#'   \mathbf{M}
#' }
#' as a function of Reticular Action Model (RAM) matrices
#' is given by
#' \deqn{
#'   \mathbf{M}
#'   =
#'   \mathbf{F}
#'   \left(
#'     \mathbf{I} - \mathbf{A}
#'   \right)^{-1}
#'   \boldsymbol{\mu}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-ram-M)
#' }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams ram_Sigmatheta
#' @inherit ram_Sigmatheta references
#' @param mu Numeric vector.
#' Vector of means
#' \eqn{
#'   \boldsymbol{\mu}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-mutheta)
#' } .
#' @return Returns the mean structure vector
#' \eqn{
#'   \mathbf{M}
#' }
#' derived from the
#' \eqn{\mathbf{A}},
#' \eqn{\mathbf{F}},
#' \eqn{\mathbf{I}},
#' matrices and
#' \eqn{
#'   \boldsymbol{\mu}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-mutheta)
#' }
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
#' F <- I <- diag(3)
#' mu <- c(100, 100, 100)
#' ram_M(A = A, F = F, I = I, mu = mu)
#' @export
ram_M <- function(A,
                  F,
                  I,
                  mu,
                  filter = TRUE) {
  # F %*% (I - A) %*% mu
  inverse <- solve(I - A)
  full <- F %*% (I - A) %*% mu
  if (filter) {
    return(F %*% full)
  } else {
    return(full)
  }
}

#' Reticular Action Model -
#' The
#' \eqn{
#'   \mathbf{S}
#' }
#' Matrix
#' and
#' The Model-Implied Variance Covariance Matrix
#' \eqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-Sigmatheta)
#' }
#'
#' @description Derives the
#' \eqn{
#'   \mathbf{S}
#' }
#' matrix
#' and
#' the model-implied variance-covariance matrix
#' \eqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-Sigmatheta)
#' }
#' using the Reticular Action Model (RAM) notation.
#'
#' @details
#' The
#' \eqn{
#'   \mathbf{S}
#' }
#' matrix
#' and
#' the model-implied variance-covariance matrix
#' \eqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-Sigmatheta)
#' }
#' are derived using
#' the \eqn{\mathbf{A}}
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
#' Vector of variances \eqn{\sigma^2}.
#' **The first element should be the variance of an exogenous variable.**
#' @param SigmaMatrix Logical.
#' If `TRUE`,
#' returns
#' the model-implied variance-covariance matrix
#' \eqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-Sigmatheta)
#' } .
#' If `FALSE`,
#' returns the
#' \eqn{
#'   \mathbf{S}
#' }
#' matrix.
#' @param both Logical.
#' If `TRUE`, returns both
#' the
#' \eqn{
#'   \mathbf{S}
#' }
#' marix
#' and
#' the model-implied variance-covariance matrix
#' \eqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-Sigmatheta)
#' }
#' in a list.
#' @return If `both = TRUE`,
#' returns both
#' the
#' \eqn{
#'   \mathbf{S}
#' }
#' matrix
#' and
#' the model-implied variance-covariance matrix
#' \eqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-Sigmatheta)
#' }
#' derived from the
#' \eqn{\mathbf{A}} matrix and
#' \eqn{\sigma^2}
#' vector.
#' in a list.
#' Otherwise,
#' If `both = FALSE`,
#' the model-implied variance-covariance matrix
#' \eqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#'   %(\#eq:sem-Sigmatheta)
#' }
#' is returned if `SigmaMatrix = TRUE`.
#' If `both = FALSE`,
#' \eqn{
#'   \mathbf{S}
#' }
#' matrix
#' is returned if `SigmaMatrix = FALSE`.
#'
#' NOTE:
#'   SigmaMatrix - option to filter out latent variables
#'     If `filter = TRUE`, variance-covariance matrix of observed variables.
#'     If `filter = FALSE`, combined variance-covariance matrix of observed and latent variables.
#'     S - always the full S matrix with the same dimensions as the input A matrix
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
#' # One-factor CFA model--------------------------------------------------------
#' A <- matrix(data = 0, ncol = 6, nrow = 6)
#' for (i in 2:6) {
#'   A[i, 1] <- 0.5
#' }
#' sigma2 <- c(1, 1, 1, 1, 1, 1)
#' F <- I <- diag(nrow(A))
#' F <- diag(nrow(A) - 1)
#' F <- cbind(0, F)
#' ram_S(
#'   A = A,
#'   sigma2 = sigma2,
#'   F = F,
#'   I = I
#' )
#' ram_S(
#'   A = A,
#'   sigma2 = sigma2,
#'   F = F,
#'   I = I,
#'   SigmaMatrix = FALSE
#' )
#' @export
ram_S <- function(A,
                  sigma2,
                  F,
                  I,
                  SigmaMatrix = TRUE,
                  filter = TRUE,
                  both = FALSE) {
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
  # apply the filter
  if (filter) {
    Sigmatheta <- ram_Sigmatheta(
      A = A,
      S = S,
      F = F,
      I = I,
      filter = TRUE
    )
  }
  if (SigmaMatrix) {
    return(Sigmatheta)
  } else {
    return(S)
  }
  if (both) {
    return(
      list(
        S = S,
        Sigmatheta = Sigmatheta
      )
    )
  }
}

#' Reticular Action Model - Residuals
#' \eqn{
#'   \boldsymbol{
#'     \hat{\Sigma}
#'   }
#'   -
#'   \boldsymbol{
#'     \hat{\Sigma}
#'   }
#'   \left(
#'     \boldsymbol{
#'       \hat{\theta}
#'     }
#'   \right)
#'   %(\#eq:sem-residuals)
#' }
#'
#' \deqn{
#'   \boldsymbol{
#'     \hat{\Sigma}
#'   }
#'   -
#'   \boldsymbol{
#'     \hat{\Sigma}
#'   }
#'   \left(
#'     \boldsymbol{
#'       \hat{\theta}
#'     }
#'   \right)
#'   %(\#eq:sem-residuals)
#' }
#' where
#' \deqn{
#'   \boldsymbol{
#'     \hat{\Sigma}
#'   }
#'   \left(
#'     \hat{
#'       \boldsymbol{\theta}
#'     }
#'   \right)
#'   =
#'   \mathbf{F}
#'   \left(
#'     \mathbf{I} - \mathbf{\hat{A}}
#'   \right)^{-1}
#'   \mathbf{\hat{S}}
#'   \left[
#'     \left(
#'       \mathbf{I} - \mathbf{\hat{A}}
#'     \right)^{-1}
#'   \right]^{\prime}
#'   \mathbf{F}^{\prime} .
#'   %(\#eq:sem-ram-Sigmahatthetahat)
#' }
#' Items with a hat (^) are estimates using the sample data.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams ram_Sigmatheta
#' @param Sigmahat Matrix.
#' Estimated variance-covariance matrix
#' \eqn{
#'   \boldsymbol{
#'     \hat{\Sigma}
#'   }
#' }
#' @param Sigmathetahat Matrix.
#' Model-implied variance-covariance matrix
#' as a function of estimated parameters
#' \eqn{
#'   \boldsymbol{
#'     \hat{\Sigma}
#'   }
#'   \left(
#'     \boldsymbol{
#'       \hat{\theta}
#'     }
#'   \right)
#' } .
#' @param Ahat \eqn{mathbf{\hat{A}}} matrix.
#' Estimated asymmetric paths (single-headed arrows),
#' such as regression coefficients and factor loadings
#' @param Shat \eqn{mathbf{\hat{S}}} matrix.
#' Estimated symmetric paths (double-headed arrows),
#' representing variances and covariances.
#' @export
ram_residuals <- function(Sigmahat,
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
