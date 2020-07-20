#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Reticular Action Model - Model-Implied Variance-Covariance Matrix
#'   \eqn{ \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'
#' @description Derives the model-implied variance-covariance matrix
#'   \eqn{ \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)
#'     %(\#eq:sem-Sigmatheta)
#'   }
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The model-implied variance-covariance matrix
#'   \eqn{ \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'   as a function of Reticular Action Model (RAM) matrices
#'   is given by
#'   \deqn{ \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)
#'     =
#'     \mathbf{F} \left( \mathbf{I} - \mathbf{A} \right)^{-1} \mathbf{S}
#'     \left[ \left( \mathbf{I} - \mathbf{A} \right)^{-1} \right]^{\prime}
#'     \mathbf{F}^{\prime}}
#'   where
#'   - \eqn{\mathbf{A}} represents asymmetric paths (single-headed arrows),
#'     such as regression coefficients and factor loadings.
#'   - \eqn{\mathbf{S}} represents symmetric paths (double-headed arrows),
#'     such as variances and covariances.
#'   - \eqn{\mathbf{F}} represents the filter matrix
#'     used to select the observed variables.
#'   - \eqn{\mathbf{I}} represents the identity matrix.
#'
#'   The dimensions of the matrices are as follows:
#'
#'   \deqn{
#'     \mathbf{A}_{p \times p}, \\
#'     \mathbf{S}_{p \times p}, \\
#'     \mathbf{I}_{p \times p}, \\
#'     \mathbf{F}_{k \times p}, \\
#'     \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)_{k \times k}
#'   }
#'
#'   where
#'
#'   - \eqn{k} number of observed variables
#'   - \eqn{q} number of latent variables
#'   - \eqn{p} number of observed and latent variables, that is \eqn{k + q}
#'
#'   \eqn{\mathbf{I}} is an identity matrix,
#'   meaning the diagonal elements are ones,
#'   and the off-diagonal elements are zeroes.
#'   Since \eqn{\mathbf{I}} has the same dimensions
#'   as \eqn{\mathbf{A}_{p \times p}}, and \eqn{\mathbf{S}_{p \times p}},
#'   the function derives \eqn{\mathbf{I}} automatically.
#'
#'   \eqn{\mathbf{F}} is a filter matrix.
#'   The number of columns \eqn{p}
#'   represent the number of observed and latent variables.
#'   The number of rows represent the number of observed variables.
#'   The off-diagonal elements are set to zero.
#'   The diagonal elements are set to zero or one.
#'   Setting the \eqn{k}-th diagonal element to **one filters in**
#'   the variable that corresponds to that position.
#'
#' @family SEM notation functions
#' @keywords matrix ram
#' @param A `p x p` numeric matrix
#'   \eqn{\mathbf{A}_{p \times p}}.
#'   Asymmetric paths (single-headed arrows),
#'   such as regression coefficients and factor loadings.
#' @param S `p x p` numeric matrix
#'   \eqn{\mathbf{S}_{p \times p}}.
#'   Symmetric paths (double-headed arrows),
#'   such as variances and covariances.
#' @param filter `k x p` numeric matrix
#'   \eqn{\mathbf{F}_{k \times p}}.
#'   Filter matrix used to select variables.
#' @return Returns the model-implied variance-covariance matrix
#'   \eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'   derived from the `A`, `S`, and `filter` matrices.
#' @examples
#' library(lavaan)
#' library(semPlot)
#' #-----------------------------------------------------------------------------
#' # Simple mediation model with observed variables------------------------------
#' # The variables in the model are X, M, and Y.
#' # All three are observed variables.
#' #-----------------------------------------------------------------------------
#' alpha <- 0.3385926
#' beta <- 0.4510391
#' tauprime <- 0.2076475
#' sigma2X <- 1.2934694
#' sigma2epsilonM <- 0.9296694
#' sigma2epsilonY <- 0.9310601
#' A <- matrix(data = 0, nrow = 3, ncol = 3)
#' A[2, 1] <- alpha
#' A[3, 1] <- beta
#' A[3, 2] <- tauprime
#' S <- diag(c(sigma2X, sigma2epsilonM, sigma2epsilonY))
#' filter <- diag(3)
#' colnames(A) <- c("X", "M", "Y")
#' rownames(A) <- c("X", "M", "Y")
#' colnames(S) <- c("X", "M", "Y")
#' rownames(S) <- c("X", "M", "Y")
#' colnames(filter) <- c("X", "M", "Y")
#' rownames(filter) <- c("X", "M", "Y")
#' Sigma <- ram_Sigmatheta(A = A, S = S, filter = filter)
#' Sigma
#'
#' #-----------------------------------------------------------------------------
#' # Simple mediation model with observed variables------------------------------
#' # and latent variable error terms---------------------------------------------
#' # The same mediation model can be parameterized by explicitly
#' # including the error terms epsilonM and epsilonY in the matrices.
#' # The variables in the model are X, M, Y, epsilonM, and epsilonY.
#' # X, M, and Y are observed variables.
#' # epsilonM, and epsilonY are latent variables.
#' #-----------------------------------------------------------------------------
#' A <- matrix(data = 0, nrow = 5, ncol = 5)
#' A[2, 1] <- alpha
#' A[2, 4] <- 1
#' A[3, 1] <- beta
#' A[3, 2] <- tauprime
#' A[3, 5] <- 1
#' S <- diag(c(sigma2X, 0, 0, sigma2epsilonM, sigma2epsilonY))
#' filter <- matrix(data = 0, nrow = 3, ncol = 5)
#' filter[1, 1] <- 1
#' filter[2, 2] <- 1
#' filter[3, 3] <- 1
#' colnames(A) <- c("X", "M", "Y", "epsilonM", "epsilonY")
#' rownames(A) <- c("X", "M", "Y", "epsilonM", "epsilonY")
#' colnames(S) <- c("X", "M", "Y", "epsilonM", "epsilonY")
#' rownames(S) <- c("X", "M", "Y", "epsilonM", "epsilonY")
#' colnames(filter) <- c("X", "M", "Y", "epsilonM", "epsilonY")
#' rownames(filter) <- c("X", "M", "Y")
#' Sigma <- ram_Sigmatheta(A = A, S = S, filter = filter)
#' Sigma
#' # model specification
#' model <- "
#'   M ~ X
#'   Y ~ X + M
#'   X ~~ X
#'   M ~~ M
#'   Y ~~ Y
#' "
#' # model fitting
#' fit <- sem(
#'   model,
#'   sample.cov = Sigma,
#'   sample.nobs = 50
#' )
#' # results
#' summary(fit, fit.measures = TRUE, standardized = TRUE)
#' semPaths(fit, what = "path", whatLabels = "est", style = "ram")
#'
#' # One-factor CFA model--------------------------------------------------------
#' # The variables in the model are eta, y1, y2, y3, y4, y5.
#' # eta is a latent variable. y1 to y5 are observed variables.
#' #-----------------------------------------------------------------------------
#' A <- matrix(data = 0, ncol = 6, nrow = 6)
#' loadings <- 0.5 # tau equivalence
#' for (i in 2:6) {
#'   A[i, 1] <- 0.5
#' }
#' S <- diag(c(1, 0.75, 0.75, 0.75, 0.75, 0.75))
#' filter <- diag(nrow(A) - 1)
#' filter <- cbind(0, filter)
#' colnames(A) <- c("eta", "y1", "y2", "y3", "y4", "y5")
#' rownames(A) <- c("eta", "y1", "y2", "y3", "y4", "y5")
#' colnames(S) <- c("eta", "y1", "y2", "y3", "y4", "y5")
#' rownames(S) <- c("eta", "y1", "y2", "y3", "y4", "y5")
#' colnames(filter) <- c("eta", "y1", "y2", "y3", "y4", "y5")
#' rownames(filter) <- c("y1", "y2", "y3", "y4", "y5")
#' Sigma <- ram_Sigmatheta(A = A, S = S, filter = filter)
#' model <- "
#'   eta =~ NA * y1 + y2 + y3 + y4 + y5
#'   eta ~~ 1 * eta
#' "
#' # model fitting
#' fit <- sem(
#'   model,
#'   sample.cov = Sigma,
#'   sample.nobs = 100
#' )
#' # results
#' summary(fit, fit.measures = TRUE, standardized = TRUE)
#' semPaths(fit, what = "path", whatLabels = "est", style = "ram")
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
                           filter) {
  invIminusA <- solve(diag(nrow(A)) - A)
  filter %*% invIminusA %*% S %*% t(invIminusA) %*% t(filter)
}

#' Reticular Action Model - Model-Implied Mean Vector
#'   \eqn{ \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'
#' @description Derives the model-implied mean vector
#'   \eqn{ \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The model-implied mean vector
#'   \eqn{ \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'   as a function of Reticular Action Model (RAM) matrices
#'   is given by
#'   \deqn{
#'     \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)
#'     = \mathbf{F} \left( \mathbf{I} - \mathbf{A} \right)^{-1} \mathbf{M}}
#'   where \eqn{\mathbf{M}} is the \eqn{p \times 1} mean structure,
#'   that is, a vector of means and intercepts.
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams ram_Sigmatheta
#' @inherit ram_Sigmatheta references
#' @param M `p x 1` numeric vector \eqn{\mathbf{M}_{p \times 1}}.
#'   Mean structure. Vector of means and intercepts.
#' @return Returns the model-implied mean vector
#'   \eqn{ \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'   derived from the \eqn{\mathbf{A}}, \eqn{\mathbf{F}}, \eqn{\mathbf{I}},
#'   matrices and \eqn{\mathbf{M}} vector.
#' @examples
#' #-----------------------------------------------------------------------------
#' # Simple mediation model with observed variables------------------------------
#' # The variables in the model are X, M, and Y.
#' # All three are observed variables.
#' #-----------------------------------------------------------------------------
#' alpha <- 0.3385926
#' beta <- 0.4510391
#' tauprime <- 0.2076475
#' muX <- 70.18000
#' deltaM <- -20.70243
#' deltaY <- -12.71288
#' A <- matrix(data = 0, nrow = 3, ncol = 3)
#' A[2, 1] <- alpha
#' A[3, 1] <- beta
#' A[3, 2] <- tauprime
#' filter <- diag(3)
#' M <- c(X = muX, M = deltaM, Y = deltaY)
#' colnames(A) <- c("X", "M", "Y")
#' rownames(A) <- c("X", "M", "Y")
#' colnames(filter) <- c("X", "M", "Y")
#' rownames(filter) <- c("X", "M", "Y")
#' ram_mutheta(M = M, A = A, filter = filter)
#'
#' #-----------------------------------------------------------------------------
#' # Simple mediation model with observed variables------------------------------
#' # and latent variable error terms---------------------------------------------
#' # The same mediation model can be parameterized by explicitly
#' # including the error terms epsilonM and epsilonY in the matrices.
#' # The variables in the model are X, M, Y, epsilonM, and epsilonY.
#' # X, M, and Y are observed variables.
#' # epsilonM, and epsilonY are latent variables.
#' #-----------------------------------------------------------------------------
#' A <- matrix(data = 0, nrow = 5, ncol = 5)
#' A[2, 1] <- alpha
#' A[2, 4] <- 1
#' A[3, 1] <- beta
#' A[3, 2] <- tauprime
#' A[3, 5] <- 1
#' filter <- matrix(data = 0, nrow = 3, ncol = 5)
#' filter[1, 1] <- 1
#' filter[2, 2] <- 1
#' filter[3, 3] <- 1
#' M <- c(X = muX, M = deltaM, Y = deltaY, epsilonM = 0, epsilonY = 0)
#' colnames(A) <- c("X", "M", "Y", "epsilonM", "epsilonY")
#' rownames(A) <- c("X", "M", "Y", "epsilonM", "epsilonY")
#' colnames(filter) <- c("X", "M", "Y", "epsilonM", "epsilonY")
#' rownames(filter) <- c("X", "M", "Y")
#' ram_mutheta(M = M, A = A, filter = filter)
#'
#' # One-factor CFA model--------------------------------------------------------
#' # The variables in the model are eta, y1, y2, y3, y4, y5.
#' # eta is a latent variable. y1 to y5 are observed variables.
#' #-----------------------------------------------------------------------------
#' A <- matrix(data = 0, ncol = 6, nrow = 6)
#' loadings <- 0.5 # tau equivalence
#' for (i in 2:6) {
#'   A[i, 1] <- 0.5
#' }
#' filter <- diag(nrow(A) - 1)
#' filter <- cbind(0, filter)
#' M <- c(eta = 0, y1 = 0, y2 = 0, y3 = 0, y4 = 0, y5 = 0)
#' colnames(A) <- c("eta", "y1", "y2", "y3", "y4", "y5")
#' rownames(A) <- c("eta", "y1", "y2", "y3", "y4", "y5")
#' colnames(filter) <- c("eta", "y1", "y2", "y3", "y4", "y5")
#' rownames(filter) <- c("y1", "y2", "y3", "y4", "y5")
#' ram_mutheta(M = M, A = A, filter = filter)
#' @export
ram_mutheta <- function(M,
                        A,
                        filter) {
  filter %*% solve(diag(nrow(A)) - A) %*% M
}

#' Reticular Action Model - Mean Structure Vector \eqn{\mathbf{M}}
#'
#' @description Derives the mean structure vector \eqn{\mathbf{M}}
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The mean structure vector \eqn{\mathbf{M}}
#'   as a function of Reticular Action Model (RAM) matrices is given by
#'   \deqn{ \mathbf{M} = \left( \mathbf{I} - \mathbf{A} \right)^{-1}
#'     \mathbf{F}^{\prime} \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'   where \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'   is the \eqn{k \times 1} model-implied mean vector.
#' @author Ivan Jacob Agaloos Pesigan
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams ram_Sigmatheta
#' @inherit ram_Sigmatheta references
#' @param mu `k x 1` numeric vector
#'   \eqn{\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)_{k \times 1}} .
#'   Model-implied meam vector.
#' @return Returns the mean structure vector \eqn{\mathbf{M}}
#'   derived from the \eqn{\mathbf{A}}, \eqn{\mathbf{F}},
#'   \eqn{ \mathbf{I}}, matrices and
#'   \eqn{ \boldsymbol{\mu} \left( \boldsymbol{\theta} \right)}
#'   vector.
#' @examples
#' #-----------------------------------------------------------------------------
#' # Simple mediation model with observed variables------------------------------
#' # The variables in the model are X, M, and Y.
#' # All three are observed variables.
#' #-----------------------------------------------------------------------------
#' alpha <- 0.3385926
#' beta <- 0.4510391
#' tauprime <- 0.2076475
#' muX <- 70.18000
#' muM <- 3.339346
#' muY <- 3.515903
#' A <- matrix(data = 0, nrow = 3, ncol = 3)
#' A[2, 1] <- alpha
#' A[3, 1] <- beta
#' A[3, 2] <- tauprime
#' filter <- diag(3)
#' mu <- c(X = muX, M = muM, Y = muY)
#' colnames(A) <- c("X", "M", "Y")
#' rownames(A) <- c("X", "M", "Y")
#' colnames(filter) <- c("X", "M", "Y")
#' rownames(filter) <- c("X", "M", "Y")
#' ram_M(mu = mu, A = A, filter = filter)
#'
#' #-----------------------------------------------------------------------------
#' # Simple mediation model with observed variables------------------------------
#' # and latent variable error terms---------------------------------------------
#' # The same mediation model can be parameterized by explicitly
#' # including the error terms epsilonM and epsilonY in the matrices.
#' # The variables in the model are X, M, Y, epsilonM, and epsilonY.
#' # X, M, and Y are observed variables.
#' # epsilonM, and epsilonY are latent variables.
#' #-----------------------------------------------------------------------------
#' A <- matrix(data = 0, nrow = 5, ncol = 5)
#' A[2, 1] <- alpha
#' A[2, 4] <- 1
#' A[3, 1] <- beta
#' A[3, 2] <- tauprime
#' A[3, 5] <- 1
#' filter <- matrix(data = 0, nrow = 3, ncol = 5)
#' filter[1, 1] <- 1
#' filter[2, 2] <- 1
#' filter[3, 3] <- 1
#' mu <- c(X = muX, M = muM, Y = muY)
#' colnames(A) <- c("X", "M", "Y", "epsilonM", "epsilonY")
#' rownames(A) <- c("X", "M", "Y", "epsilonM", "epsilonY")
#' colnames(filter) <- c("X", "M", "Y", "epsilonM", "epsilonY")
#' rownames(filter) <- c("X", "M", "Y")
#' ram_M(mu = mu, A = A, filter = filter)
#' #'
#' # One-factor CFA model--------------------------------------------------------
#' # The variables in the model are eta, y1, y2, y3, y4, y5.
#' # eta is a latent variable. y1 to y5 are observed variables.
#' #-----------------------------------------------------------------------------
#' A <- matrix(data = 0, ncol = 6, nrow = 6)
#' loadings <- 0.5 # tau equivalence
#' for (i in 2:6) {
#'   A[i, 1] <- 0.5
#' }
#' filter <- diag(nrow(A) - 1)
#' filter <- cbind(0, filter)
#' mu <- c(y1 = 0, y2 = 0, y3 = 0, y4 = 0, y5 = 0)
#' colnames(A) <- c("eta", "y1", "y2", "y3", "y4", "y5")
#' rownames(A) <- c("eta", "y1", "y2", "y3", "y4", "y5")
#' colnames(filter) <- c("eta", "y1", "y2", "y3", "y4", "y5")
#' rownames(filter) <- c("y1", "y2", "y3", "y4", "y5")
#' ram_M(mu = mu, A = A, filter = filter)
#' @export
ram_M <- function(mu,
                  A,
                  filter) {
  solve(diag(nrow(A)) - A) %*% t(filter) %*% mu
}

#' Reticular Action Model - The \eqn{ \mathbf{S} } Matrix
#'
#' @description Derives the \eqn{ \mathbf{S} } matrix
#'   using the Reticular Action Model (RAM) notation
#'   from variable variances \eqn{\sigma^2}.
#'   The off-diagonal elements of the \eqn{ \mathbf{S} } matrix are assumed to be zeroes.
#'
#' @details The \eqn{ \mathbf{S} } matrix
#'   is derived using
#'   the \eqn{\mathbf{A}}
#'   matrix
#'   and sigma squared
#'   \eqn{\left( \sigma^2 \right)}
#'   vector (variances).
#'   **Note that the first or last (see `start` argument) element
#'   in the `A` and `S` matrices should be an exogenous variable.**
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams ram_Sigmatheta
#' @inherit ram_Sigmatheta references
#' @param sigma2 Numeric vector.
#'   Vector of variances \eqn{\sigma^2}.
#' @param start Logical.
#'   If `TRUE`, an exogenous variable is positioned in the first position.
#'   If `FALSE`, an exogenous variable is positioned in the last position.
#' @return If `both = TRUE`, returns both the \eqn{ \mathbf{S} } matrix
#' is derived from the \eqn{\mathbf{A}} matrix and \eqn{\sigma^2} vector.
#' @examples
#' #-----------------------------------------------------------------------------
#' # Simple mediation model with observed variables------------------------------
#' # The variables in the model are X, M, and Y.
#' # All three are observed variables.
#' #-----------------------------------------------------------------------------
#' alpha <- 0.3385926
#' beta <- 0.4510391
#' tauprime <- 0.2076475
#' sigma2X <- 1.2934694
#' sigma2M <- 1.0779592
#' sigma2Y <- 1.2881633
#' A <- matrix(data = 0, nrow = 3, ncol = 3)
#' A[2, 1] <- alpha
#' A[3, 1] <- beta
#' A[3, 2] <- tauprime
#' sigma2 <- c(X = sigma2X, M = sigma2M, Y = sigma2Y)
#' colnames(A) <- c("X", "M", "Y")
#' rownames(A) <- c("X", "M", "Y")
#' ram_sigma2(sigma2 = sigma2, A = A)
#'
#' #-----------------------------------------------------------------------------
#' # Simple mediation model with observed variables------------------------------
#' # and latent variable error terms---------------------------------------------
#' # The same mediation model can be parameterized by explicitly
#' # including the error terms epsilonM and epsilonY in the matrices.
#' # The variables in the model are X, M, Y, epsilonM, and epsilonY.
#' # X, M, and Y are observed variables.
#' # epsilonM, and epsilonY are latent variables.
#' #-----------------------------------------------------------------------------
#' A <- matrix(data = 0, nrow = 5, ncol = 5)
#' A[2, 1] <- alpha
#' A[2, 4] <- 1
#' A[3, 1] <- beta
#' A[3, 2] <- tauprime
#' A[3, 5] <- 1
#' sigma2 <- c(X = sigma2X, M = sigma2M, Y = sigma2Y)
#' colnames(A) <- c("X", "M", "Y", "epsilonM", "epsilonY")
#' rownames(A) <- c("X", "M", "Y", "epsilonM", "epsilonY")
#' ram_sigma2(sigma2 = sigma2, A = A)
#'
#' # One-factor CFA model--------------------------------------------------------
#' # The variables in the model are eta, y1, y2, y3, y4, y5.
#' # eta is a latent variable. y1 to y5 are observed variables.
#' #-----------------------------------------------------------------------------
#' A <- matrix(data = 0, ncol = 6, nrow = 6)
#' loadings <- 0.5 # tau equivalence
#' for (i in 2:6) {
#'   A[i, 1] <- 0.5
#' }
#' sigma2 <- c(eta = 1, y1 = 0.75, y2 = 0.75, y3 = 0.75, y4 = 0.75, y5 = 0.75)
#' colnames(A) <- c("eta", "y1", "y2", "y3", "y4", "y5")
#' rownames(A) <- c("eta", "y1", "y2", "y3", "y4", "y5")
#' ram_sigma2(sigma2 = sigma2, A = A)
#'
#' # One-factor CFA model--------------------------------------------------------
#' # The variables in the model are eta, y1, y2, y3, y4, y5.
#' # eta is a latent variable. y1 to y5 are observed variables.
#' #-----------------------------------------------------------------------------
#' A <- matrix(data = 0, ncol = 6, nrow = 6)
#' loadings <- 0.5 # tau equivalence
#' for (i in 1:5) {
#'   A[i, 6] <- 0.5
#' }
#' sigma2 <- c(y1 = 0.75, y2 = 0.75, y3 = 0.75, y4 = 0.75, y5 = 0.75, eta = 1)
#' colnames(A) <- c("y1", "y2", "y3", "y4", "y5", "eta")
#' rownames(A) <- c("y1", "y2", "y3", "y4", "y5", "eta")
#' ram_sigma2(sigma2 = sigma2, A = A, start = FALSE)
#' @export
ram_sigma2 <- function(sigma2,
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
  Sigmatheta_temp <- ram_Sigmatheta(
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
    Sigmatheta_temp <- ram_Sigmatheta(
      A = A,
      S = S,
      filter = filter
    )
  }
  S
}

#' @author Ivan Jacob Agaloos Pesigan
#'
#' @title Reticular Action Model - The \eqn{ \mathbf{S} } Matrix
#'
#' @description Derives the \eqn{ \mathbf{S} } matrix
#'   using the Reticular Action Model (RAM) notation.
#'
#' @details The \eqn{ \mathbf{S} } matrix
#' is derived using
#' the \eqn{\mathbf{A}}
#' matrix
#' and the model-implied variance-covariance matrix
#' \eqn{ \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @family SEM notation functions
#' @keywords matrix ram
#' @inheritParams ram_Sigmatheta
#' @param Sigmatheta Matrix.
#'   Model-implied variance-covariance matrix
#'   \eqn{ \boldsymbol{ \Sigma } \left( \boldsymbol{ \theta } \right) } .
#' @export
ram_S <- function(Sigmatheta,
                  A,
                  filter) {
  invIminusA <- solve(I - A)
  filter %*% invIminusA %*% Sigmatheta %*% t(invIminusA) %*% t(filter)
}

#' Reticular Action Model - Residuals
#' \eqn{ \boldsymbol{ \hat{\Sigma} }
#'   -
#'   \boldsymbol{
#'     \hat{\Sigma}
#'   }
#'   \left(
#'     \boldsymbol{
#'       \hat{\theta}
#'     }
#'   \right)
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
#'   Estimated asymmetric paths (single-headed arrows),
#'   such as regression coefficients and factor loadings
#' @param Shat \eqn{mathbf{\hat{S}}} matrix.
#'   Estimated symmetric paths (double-headed arrows),
#'   representing variances and covariances.
#' @export
ram_residuals <- function(Sigmahat,
                          Sigmathetahat = NULL,
                          Ahat,
                          Shat,
                          filter) {
  if (is.null(Sigmathetahat)) {
    Sigmathetahat <- ram_Sigmatheta(
      A = Ahat,
      S = Shat,
      filter = filter
    )
  }
  Sigmahat - Sigmathetahat
}
