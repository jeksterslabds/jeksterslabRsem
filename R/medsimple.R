#' Simple Mediation - Estimator
#'
#' @author Ivan Jacob Agaloos Pesigan
#' @param theta Numeric vector.
#'   Vector of parameters to estimate.
#'   `theta[1] =` \eqn{\alpha},
#'   `theta[2] =` \eqn{\tau^{\prime}},
#'   `theta[3] =` \eqn{\beta},
#'   `theta[4] =` \eqn{\sigma^{2}_{X}},
#'   `theta[5] =` \eqn{\sigma^{2}_{\varepsilon_{M}}},
#'   `theta[6] =` \eqn{\sigma^{2}_{\varepsilon_{Y}}}.
#' @param Sigmahat Matrix.
#'   \eqn{3 \times 3} sample variance-covariance matrix.
#'   `Sigmahat[1,1] =` \eqn{\mathrm{Var} \left( X \right)},
#'   `Sigmahat[2,2] =` \eqn{\mathrm{Var} \left( M \right)},
#'   `Sigmahat[3,3] =` \eqn{\mathrm{Var} \left( Y \right)}.
#' @param obj Function.
#'   Objectove function.
#'   `fml` for [`fml()`],
#'   `fgls` for [`fgls()`],
#'   `fuls` for [`fuls()`].
#' @examples
#' X <- jeksterslabRdatarepo::thirst
#' Sigmahat <- cov(X)
#' jeksterslabRdist::opt(
#'   FUN = medsimpleobj,
#'   start_values = rep(x = 0.20, times = 6),
#'   optim = TRUE,
#'   Sigmahat = Sigmahat
#' )
#' @export
medsimpleobj <- function(theta,
                         Sigmahat,
                         obj = fml) {
  # Ensure that variance and residual variances are positive
  vars <- c(
    theta[4],
    theta[5],
    theta[6]
  )
  vars <- ifelse(
    test = vars < 0,
    yes = TRUE,
    no = FALSE
  )
  if (any(vars)) {
    return(NA)
  }
  A <- matrix(
    data = 0,
    ncol = 3,
    nrow = 3
  )
  A[2, 1] <- theta[1]
  A[3, 1] <- theta[2]
  A[3, 2] <- theta[3]
  S <- matrix(
    data = 0,
    ncol = 3,
    nrow = 3
  )
  S[1, 1] <- theta[4]
  S[2, 2] <- theta[5]
  S[3, 3] <- theta[6]
  I <- F <- diag(3)
  # return NA if there are errors in computing Sigmatheta
  Sigmatheta <- tryCatch(
    {
      ram_Sigmatheta(
        A = A,
        S = S,
        F = F,
        I = I
      )
    },
    error = function(e) {
      return(NA)
    }
  )
  # return NA if there are errors in fml
  out <- tryCatch(
    {
      obj(
        Sigmatheta = Sigmatheta,
        Sigma = Sigmahat
      )
    },
    error = function(e) {
      return(NA)
    }
  )
  out
}
