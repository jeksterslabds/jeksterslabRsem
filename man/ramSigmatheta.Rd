% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ram.R
\name{ramSigmatheta}
\alias{ramSigmatheta}
\title{Reticular Action Model - Model-Implied Variance-Covariance Matrix
\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}}
\usage{
ramSigmatheta(A, S, filter)
}
\arguments{
\item{A}{\verb{m x m} numeric matrix
\eqn{\mathbf{A}_{m \times m}}.
Asymmetric paths (single-headed arrows),
such as regression coefficients and factor loadings.}

\item{S}{\verb{m x m} numeric matrix
\eqn{\mathbf{S}_{m \times m}}.
Symmetric paths (double-headed arrows),
such as variances and covariances.}

\item{filter}{\verb{k x m} numeric matrix
\eqn{\mathbf{F}_{k \times m}}.
Filter matrix used to select variables.}
}
\value{
Returns the model-implied variance-covariance matrix
\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
derived from the \code{A}, \code{S}, and \code{filter} matrices.
}
\description{
Derives the model-implied variance-covariance matrix
\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
using the Reticular Action Model (RAM) notation.
}
\details{
The model-implied variance-covariance matrix
\eqn{\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)}
as a function of Reticular Action Model (RAM) matrices
is given by

\deqn{
    \boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)
    =
    \mathbf{F} \left( \mathbf{I} - \mathbf{A} \right)^{-1} \mathbf{S}
    \left[ \left( \mathbf{I} - \mathbf{A} \right)^{-1} \right]^{T} \mathbf{F}^{T}
  }

where
\itemize{
\item \eqn{\mathbf{A}_{m \times m}} represents asymmetric paths (single-headed arrows),
such as regression coefficients and factor loadings,
\item \eqn{\mathbf{S}_{m \times m}} represents symmetric paths (double-headed arrows),
such as variances and covariances,
\item \eqn{\mathbf{F}_{k \times m}} represents the filter matrix
used to select the observed variables,
\item \eqn{\mathbf{I}_{m \times m}} represents an identity matrix,
\item \eqn{k} number of observed variables,
\item \eqn{q} number of latent variables, and
\item \eqn{m} number of observed and latent variables, that is \eqn{k + q} .
}
}
\references{
McArdle, J. J. (2013).
The development of the RAM rules for latent variable structural equation modeling.
In A. Maydeu-Olivares & J. J. McArdle (Eds.),
\emph{Contemporary Psychometrics: A festschrift for Roderick P. McDonald} (pp. 225--273).
Lawrence Erlbaum Associates.

McArdle, J. J., & McDonald, R. P. (1984).
Some algebraic properties of the Reticular Action Model for moment structures.
\emph{British Journal of Mathematical and Statistical Psychology}, \emph{37} (2), 234--251.
}
\seealso{
Other SEM notation functions: 
\code{\link{ramM}()},
\code{\link{rammutheta}()},
\code{\link{ramsigma2}()}
}
\author{
Ivan Jacob Agaloos Pesigan
}
\concept{SEM notation functions}
\keyword{matrix}
\keyword{ram}
