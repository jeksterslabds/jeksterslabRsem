#' Reticular Action Model (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'
#' Model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   using the Reticular Action Model notation.
#'
#' The model-implied variance-covariance matrix
#' (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#' as a function of Reticular Action Model matrices
#' is given by:
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
#' represent asymmetric paths,
#' such as regression coefficients and factor loadings.
#' \eqn{\mathbf{S}}
#' represent symmetric matrix
#'   representing variances and covariances.
#' \eqn{\mathbf{F}}
#' represent filter matrix
#'   used to select the observed variables.
#' \eqn{\mathbf{I}}
#' represent identity matrix.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param A Asymmetric paths,
#'   such as regression coefficients and factor loadings.
#' @param S Symmetric matrix
#'   representing variances and covariances.
#' @param F Filter matrix
#'   used to select the observed variables.
#' @param I Identity matrix.
#' @return Returns the model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   derived from the \eqn{\mathbf{A}},
#'   \eqn{\mathbf{S}},
#'   \eqn{\mathbf{F}}, and
#'   \eqn{\mathbf{I}} matrices.
#' @references
#'   McArdle, J. J. (2013).
#'     The development of the RAM rules for latent variable structural equation modeling.
#'     In A. Maydeu-Olivares & J. J. McArdle (Eds.),
#'     \emph{Contemporary Psychometrics: A festschrift for Roderick P. McDonald} (pp. 225--273).
#'     Lawrence Erlbaum Associates.
#'
#'   McArdle, J. J., & McDonald, R. P. (1984).
#'     Some algebraic properties of the Reticular Action Model for moment structures.
#'     \emph{British Journal of Mathematical and Statistical Psychology, 37}(2), 234--251.
#' @family SEM notation functions
#' @keywords matrix ram
#' @examples
#' A <- matrix(
#'   data = c(
#'     0, 0.26^(1 / 2), 0,
#'     0, 0, 0.26^(1 / 2),
#'     0, 0, 0
#'   ),
#'   ncol = 3
#' )
#' S <- F <- I <- diag(3)
#' S[1, 1] <- 225
#' S[2, 2] <- 166.5
#' S[3, 3] <- 166.5
#' ram_Sigmatheta(A = A, S = S, F = F, I = I)
#' @export
ram_Sigmatheta <- function(A,
                           S,
                           F,
                           I) {
  x <- solve(I - A)
  F %*% x %*% S %*% t(x) %*% t(F)
}

#' Reticular Action Model (\eqn{\boldsymbol{\mu}})
#'
#' Model-implied \eqn{\boldsymbol{\mu}} vector
#'   using the Reticular Action Model notation.
#'
#' @details \deqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)
#'   =
#'   \mathbf{F}
#'   \left( \mathbf{I} - \mathbf{A} \right)^{-1}
#'   \mathbf{M}
#' }
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams ram_Sigmatheta
#' @param M Vector of means and intercepts.
#' @return Returns the expected values (\eqn{\boldsymbol{\mu}})
#'   derived from the \eqn{\mathbf{A}},
#'   \eqn{\mathbf{F}},
#'   \eqn{\mathbf{I}}, matrices and
#'   \eqn{\mathbf{M}} vector.
#' @inherit ram_Sigmatheta references
#' @family SEM notation functions
#' @keywords matrix ram
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

#' Reticular Action Model (\eqn{\mathbf{M}})
#'
#' Mean Structure
#'   (\eqn{\mathbf{M}} vector)
#'   using the Reticular Action Model notation.
#'
#' @details \deqn{\mathbf{M}
#'   =
#'   \mathbf{F}
#'   \left( \mathbf{I} - \mathbf{A} \right)^{-1}
#'   \boldsymbol{\mu}
#' }
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams ram_Sigmatheta
#' @param mu Vector of expected values (\eqn{\boldsymbol{\mu}}).
#' @return Returns the mean structure (\eqn{\mathbf{M}} vector)
#'   derived from the \eqn{\mathbf{A}},
#'   \eqn{\mathbf{F}},
#'   \eqn{\mathbf{I}}, matrices and
#'   \eqn{\boldsymbol{\mu}} vector.
#' @inherit ram_Sigmatheta references
#' @family SEM notation functions
#' @keywords matrix ram
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

#' Reticular Action Model (\eqn{\boldsymbol{\Sigma}} or \eqn{\mathbf{S}})
#'
#' Model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   or \eqn{\mathbf{S}} Matrix
#'   using the Reticular Action Model notation.
#'
#' @details Derives the
#'   model-implied variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'   or the
#'   \eqn{\mathbf{S}}
#'   matrix
#'   from the
#'   \eqn{\mathbf{A}}
#'   matrix
#'   and sigma squared
#'   (\eqn{\sigma^2})
#'   vector (variances).
#'   \strong{Note that the first element in the matrix should be an exogenous variable.}
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams ram_Sigmatheta
#' @param sigma2 Vector of variances (\eqn{\sigma^2}).
#' @param SigmaMatrix Logical.
#'   If \code{TRUE},
#'     returns
#'     the model-implied variance-covariance matrix
#'     (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}).
#'   If \code{FALSE},
#'     returns the
#'     \eqn{\mathbf{S}}
#'     matrix.
#' @return Returns
#'     the model-implied variance-covariance matrix
#'     (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#'     or the
#'     \eqn{\mathbf{S}}
#'     matrix
#'     derived from the
#'     \eqn{\mathbf{A}}
#'     matrix and
#'     \eqn{\sigma^2}
#'     vector.
#' @inherit ram_Sigmatheta references
#' @family SEM notation functions
#' @keywords matrix ram
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
  S <- matrix(
    data = 0,
    ncol = dim(A)[1],
    nrow = dim(A)[2]
  )
  Sigmatheta <- ram_Sigmatheta(
    A = A,
    S = S,
    F = F,
    I = I
  )
  for (i in 1:nrow(A)) {
    S[i, i] <- sigma2[i] - Sigmatheta[i, i]
    Sigmatheta <- ram_Sigmatheta(
      A = A,
      S = S,
      F = F,
      I = I
    )
  }
  if (SigmaMatrix) {
    return(Sigmatheta)
  } else {
    return(S)
  }
}

#' Reticular Action Model
#' (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#' from Parameters of a \eqn{k}-Variable Linear Regression Model)
#'
#' Model-implied variance-covariance matrix
#' (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)})
#' using the Reticular Action Model notation
#' from parameters of a \eqn{k}-variable linear regression model.
#'
#' The parameters of a linear regression model
#' for the covariance structure
#' are
#' the slopes
#' (\eqn{\boldsymbol{\beta}_{\mathrm{slopes}}}),
#' the variance of the error term \eqn{\boldsymbol{\varepsilon}}
#' (\eqn{\sigma^2}),
#' and
#' the covariances of
#'   \eqn{
#'     {X}_{2},
#'     {X}_{3},
#'     \dots,
#'     {X}_{k}
#'   }
#' (\eqn{\boldsymbol{\Sigma}_{\mathbf{X}}}).
#' The parameters form the following matrices,
#' \deqn{
#'   \mathbf{A}_{k \times k}
#'   =
#'   \begin{bmatrix}
#'     0 & \beta_2 & \beta_3 & \cdots & \beta_k \\
#'     0 & 0       & 0       & \cdots & 0       \\
#'     0 & 0       & 0       & \cdots & 0       \\
#'     0 & 0       & 0       & \cdots & 0       \\
#'     0 & 0       & 0       & \cdots & 0
#'   \end{bmatrix} ,
#' }
#' \deqn{
#'   \mathbf{S}_{k \times k}
#'   =
#'   \begin{bmatrix}
#'     \sigma^2 & 0                         & 0                         & \cdots & 0   \\
#'     0        & \mathrm{Var}_{X_2}        & 0                         & \cdots & 0   \\
#'     0        & \mathrm{Cov}_{X_{3}X_{2}} & \mathrm{Var}_{X_3}        & \cdots & 0   \\
#'     0        & \vdots                    & \vdots                    & \ddots & 0   \\
#'     0        & \mathrm{Cov}_{X_{k}X_{2}} & \mathrm{Cov}_{X_{k}X_{3}} & \cdots & \mathrm{Var}_{X_k}
#'   \end{bmatrix} ,
#' }
#' \deqn{
#'   \mathbf{I}_{k \times k}
#'   =
#'   \begin{bmatrix}
#'     1 & 0   & 0   & \cdots & 0   \\
#'     0 & 1   & 0   & \cdots & 0   \\
#'     0 & 0   & 1   & \cdots & 0   \\
#'     0 & 0   & 0   & 1      & 0   \\
#'     0 & 0   & 0   & \cdots & 1
#'   \end{bmatrix} ,
#' }
#' and
#' \deqn{
#'   \mathbf{F}_{k \times k}
#'   =
#'   \begin{bmatrix}
#'     1 & 0   & 0   & \cdots & 0   \\
#'     0 & 1   & 0   & \cdots & 0   \\
#'     0 & 0   & 1   & \cdots & 0   \\
#'     0 & 0   & 0   & 1      & 0   \\
#'     0 & 0   & 0   & \cdots & 1
#'   \end{bmatrix} .
#' }
#' The model-implied variance-covariance matrix
#' (\eqn{
#'   \boldsymbol{\Sigma}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#' })
#' is calculated using
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
#'     \mathbf{F}^{\prime} .
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param beta_slopes Vector or `k - 1` by `1` matrix.
#'   \eqn{\left( k - 1 \right) \times 1}
#'   vector of regression slopes.
#' @param sigma2 Numeric.
#'   Variance of the error term
#'   \eqn{\boldsymbol{\varepsilon}}
#'   (\eqn{\sigma^2}).
#' @param SigmaX Matrix.
#'   Covariances of
#'   \eqn{
#'     {X}_{2},
#'     {X}_{3},
#'     \dots,
#'     {X}_{k}
#'   }
#'   (\eqn{\boldsymbol{\Sigma}_{\mathbf{X}}}).
#' @return
#'   Returns the model-implied
#'   variance-covariance matrix
#'   (\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}).
#' @export
ram_Sigmatheta_linreg <- function(beta_slopes,
                                  sigma2,
                                  SigmaX) {
  beta_slopes <- as.vector(beta_slopes)
  p <- length(beta_slopes)
  k <- p + 1
  A <- matrix(
    data = 0,
    nrow = k,
    ncol = k
  )
  A[1, ] <- c(
    0,
    beta_slopes
  )
  I <- F <- diag(k)
  SigmaX[upper.tri(SigmaX)] <- 0
  S <- SigmaX
  S <- cbind(
    0,
    S
  )
  S <- rbind(
    0,
    S
  )
  S[1, 1] <- sigma2
  ram_Sigmatheta(
    A = A,
    S = S,
    F = F,
    I = I
  )
}

#' Reticular Action Model
#' (\eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#' from Parameters of a \eqn{k}-Variable Linear Regression Model)
#'
#' Model-implied mean vector
#' (\eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)})
#' using the Reticular Action Model notation
#' from parameters of a \eqn{k}-variable linear regression model.
#'
#' The parameters of a linear regression model
#' for the mean structure
#' are
#' the regression coefficients
#' (\eqn{\boldsymbol{\beta}}),
#' and
#' the means of
#' the regressors
#'   \eqn{
#'     {X}_{2},
#'     {X}_{3},
#'     \dots,
#'     {X}_{k}
#'   }
#' excluding the mean for the vector of constants
#' (\eqn{\boldsymbol{\mu}_{\mathbf{X}}}).
#' The parameters form the following matrices,
#' \deqn{
#'   \mathbf{A}_{k \times k}
#'   =
#'   \begin{bmatrix}
#'     0 & \beta_2 & \beta_3 & \cdots & \beta_k \\
#'     0 & 0       & 0       & \cdots & 0       \\
#'     0 & 0       & 0       & \cdots & 0       \\
#'     0 & 0       & 0       & \cdots & 0       \\
#'     0 & 0       & 0       & \cdots & 0
#'   \end{bmatrix} ,
#' }
#' \deqn{
#'   \mathbf{M}_{k \times 1}
#'   =
#'   \begin{bmatrix}
#'     \beta_1   \\
#'     \mu_{X_2} \\
#'     \mu_{X_3} \\
#'     \vdots    \\
#'     \mu_{X_k}
#'   \end{bmatrix} ,
#' }
#' \deqn{
#'   \mathbf{I}_{k \times k}
#'   =
#'   \begin{bmatrix}
#'     1 & 0   & 0   & \cdots & 0   \\
#'     0 & 1   & 0   & \cdots & 0   \\
#'     0 & 0   & 1   & \cdots & 0   \\
#'     0 & 0   & 0   & 1      & 0   \\
#'     0 & 0   & 0   & \cdots & 1
#'   \end{bmatrix} ,
#' }
#' and
#' \deqn{
#'   \mathbf{F}_{k \times k}
#'   =
#'   \begin{bmatrix}
#'     1 & 0   & 0   & \cdots & 0   \\
#'     0 & 1   & 0   & \cdots & 0   \\
#'     0 & 0   & 1   & \cdots & 0   \\
#'     0 & 0   & 0   & 1      & 0   \\
#'     0 & 0   & 0   & \cdots & 1
#'   \end{bmatrix} .
#' }
#' The model-implied mean vector
#' (\eqn{
#'   \boldsymbol{\mu}
#'   \left(
#'     \boldsymbol{\theta}
#'   \right)
#' })
#' is calculated using
#'   \deqn{
#'     \boldsymbol{\mu}
#'     \left(
#'       \boldsymbol{\theta}
#'     \right)
#'     =
#'     \mathbf{F}
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'     \mathbf{M} .
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param beta Vector or `k` by `1` matrix.
#'   The vector
#'   \eqn{\boldsymbol{\beta}}
#'   is a \eqn{k \times 1} vector
#'   of \eqn{k} regression coefficients.
#' @param muX Vector or `k - 1` by `1` matrix.
#'   \eqn{\left( k - 1 \right) \times 1} vector
#'   of means of
#'   the regressors
#'     \eqn{
#'       {X}_{2},
#'       {X}_{3},
#'       \dots,
#'       {X}_{k}
#'     }
#'   excluding the mean for the vector of constants
#'   (\eqn{\boldsymbol{\mu}_{\mathbf{X}}}).
#' @export
ram_mutheta_linreg <- function(beta,
                               muX) {
  beta <- as.vector(beta)
  k <- length(beta)
  beta_intercept <- beta[1]
  beta_slopes <- beta[-1]
  A <- matrix(
    data = 0,
    nrow = k,
    ncol = k
  )
  A[1, ] <- c(
    0,
    beta_slopes
  )
  M <- c(
    beta_intercept,
    as.vector(muX)
  )
  I <- F <- diag(k)
  ram_mutheta(
    A = A,
    F = F,
    I = I,
    M = M
  )
}

#' Reticular Action Model
#' (Mean Structure \eqn{\mathbf{M}}
#' from Regression Slopes
#' and Means of Observed Variables
#' of a \eqn{k}-Variable Linear Regression Model)
#'
#' Mean structure column vector
#' (\eqn{\mathbf{M}})
#' using the Reticular Action Model notation
#' from regression slopes
#' and means of observed variables
#' of a \eqn{k}-variable Linear Regression Model.
#'
#' The parameters of a linear regression model
#' for the mean structure
#' are
#' the slopes
#' (\eqn{\boldsymbol{\beta}_{\mathrm{slopes}}}),
#' the mean of
#' the regressand
#' \eqn{\mu_\mathbf{y}}
#' and
#' the means of
#' the regressors
#'   \eqn{
#'     {X}_{2},
#'     {X}_{3},
#'     \dots,
#'     {X}_{k}
#'   }
#' excluding the mean for the vector of constants
#' (\eqn{\boldsymbol{\mu}_{\mathbf{X}}}).
#' The parameters form the following matrices,
#' \deqn{
#'   \mathbf{A}_{k \times k}
#'   =
#'   \begin{bmatrix}
#'     0 & \beta_2 & \beta_3 & \cdots & \beta_k \\
#'     0 & 0       & 0       & \cdots & 0       \\
#'     0 & 0       & 0       & \cdots & 0       \\
#'     0 & 0       & 0       & \cdots & 0       \\
#'     0 & 0       & 0       & \cdots & 0
#'   \end{bmatrix} ,
#' }
#' \deqn{
#'   \boldsymbol{\mu}_{k \times 1}
#'   =
#'   \begin{bmatrix}
#'     \mu_{\mathbf{y}} \\
#'     \mu_{X_2} \\
#'     \mu_{X_3} \\
#'     \vdots    \\
#'     \mu_{X_k}
#'   \end{bmatrix} ,
#' }
#' \deqn{
#'   \mathbf{I}_{k \times k}
#'   =
#'   \begin{bmatrix}
#'     1 & 0   & 0   & \cdots & 0   \\
#'     0 & 1   & 0   & \cdots & 0   \\
#'     0 & 0   & 1   & \cdots & 0   \\
#'     0 & 0   & 0   & 1      & 0   \\
#'     0 & 0   & 0   & \cdots & 1
#'   \end{bmatrix} ,
#' }
#' and
#' \deqn{
#'   \mathbf{F}_{k \times k}
#'   =
#'   \begin{bmatrix}
#'     1 & 0   & 0   & \cdots & 0   \\
#'     0 & 1   & 0   & \cdots & 0   \\
#'     0 & 0   & 1   & \cdots & 0   \\
#'     0 & 0   & 0   & 1      & 0   \\
#'     0 & 0   & 0   & \cdots & 1
#'   \end{bmatrix} .
#' }
#' The mean structure column vector
#' (\eqn{\mathbf{M}})
#' is calculated using
#'   \deqn{
#'     \mathbf{M}
#'     =
#'     \mathbf{F}
#'     \left(
#'       \mathbf{I} - \mathbf{A}
#'     \right)^{-1}
#'     \boldsymbol{\mu} .
#'   }
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @inheritParams ram_Sigmatheta_linreg
#' @inheritParams ram_mutheta_linreg
#' @param muy Numeric.
#'   Mean of the regressand variable
#'   \eqn{\mu_{\mathbf{y}}}.
#' @export
ram_M_linreg <- function(beta_slopes,
                         muy,
                         muX) {
  beta_slopes <- as.vector(beta_slopes)
  p <- length(beta_slopes)
  k <- p + 1
  A <- matrix(
    data = 0,
    nrow = k,
    ncol = k
  )
  A[1, ] <- c(
    0,
    beta_slopes
  )
  I <- F <- diag(k)
  mu <- c(
    muy,
    as.vector(muX)
  )
  ram_M(
    A = A,
    F = F,
    I = I,
    mu = mu
  )
}

#' Reticular Action Model (Residuals)
#'
#' @param Sigmahat Matrix.
#'   Estimated variance-covariance matrix
#'   (\eqn{\hat{\boldsymbol{\Sigma}}})
#' @param Sigmathetahat Matrix.
#'   Model-implied variance-covariance matrix
#'   as a function of estimated parameters
#'   (\eqn{\boldsymbol{\Sigma}\left( \hat{\boldsymbol{\theta}} \right)}).
#' @param Ahat Estimated asymmetric paths,
#'   such as regression coefficients and factor loadings.
#' @param Shat Estimated symmetric matrix
#'   representing variances and covariances.
#' @inheritParams ram_Sigmatheta
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
      I = I
    )
  }
  Sigmahat - Sigmathetahat
}
