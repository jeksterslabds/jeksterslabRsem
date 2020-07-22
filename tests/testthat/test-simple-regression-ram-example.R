#' ---
#' title: "Tests: The Simple Linear Regression Model (RAM)"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Tests: The Simple Linear Regression Model (RAM)}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
#'
#+ include = FALSE
knitr::opts_chunk$set(
  error = TRUE,
  collapse = TRUE,
  comment = "#>",
  out.width = "100%"
)
#'
#'
# Structural Equation Modeling - Simple Linear Regression RAM Example {#simple-regression-ram-example}
#'
#+ echo = FALSE, message = FALSE, warning = FALSE
library(testthat)
library(jeksterslabRsem)
#'
#'
#+ echo = FALSE
varnames <- c("y", "x")
wages <- jeksterslabRdatarepo::wages
x <- wages$education
y <- wages$wages
n <- length(x)
obj <- lm(y ~ x)
beta <- unname(coef(obj))
names(beta) <- c("beta1", "beta2")
# covariance structure
beta1 <- beta[1]
beta2 <- beta[2]
sigma2x <- var(x)
sigma2y <- var(y)
# sigma^2 has some discrepancy
# sigma2epsilon <- summary(obj)$sigma^2
sigma2epsilon <- sigma2y - (beta2^2 * sigma2x)
# mean structure
mux <- mean(x)
#'
#'
#+ echo = FALSE
x <- rnorm(n = n, mean = mux, sd = sqrt(sigma2x))
epsilon <- rnorm(n = n, mean = 0, sd = sqrt(sigma2epsilon))
y <- beta1 + beta2 * x + epsilon
obj <- lm(y ~ x)
beta <- unname(coef(obj))
names(beta) <- c("beta1", "beta2")
# covariance structure
beta1 <- beta[1]
beta2 <- beta[2]
sigma2x <- var(x)
sigma2y <- var(y)
sigmayx <- cov(y, x)
ryx <- cor(y, x)
Sigmatheta <- matrix(
  data = c(
    sigma2y,
    sigmayx,
    sigmayx,
    sigma2x
  ),
  nrow = 2,
  ncol = 2
)
# sigma^2 has some discrepancy
# sigma2epsilon <- summary(obj)$sigma^2
sigma2epsilon <- sigma2y - (beta2^2 * sigma2x)
# mean structure
mux <- mean(x)
muy <- mean(y)
mutheta <- matrix(
  data = c(muy, mux),
  ncol = 1
)
#'
#'
#' In this example, we show how the Reticular Action Model (RAM) Matrices
#' are specified using the simple linear regression model.
#'
#' We also demonstrate how to use
#' RAM notation functions in the `jeksterslabRsem` package namely
#'
#' - `jeksterslabRsem::rammutheta()` - used to derive model-implied mean vector
#'    $\boldsymbol{\mu} \left( \boldsymbol{\theta} \right)$
#' - `jeksterslabRsem::ramSigmatheta()` - used to derive model-implied variance-covariance matrix
#'   $\boldsymbol{\Sigma} \left( \boldsymbol{\theta} \right)$
#' - `jeksterslabRsem::ramM()` - used to derive the $\mathbf{M}$ matrix
#' - `jeksterslabRsem::ramsigma2()` - used to derive the $\mathbf{S}$ matrix when the diagonals are all zeroes
#'    from variances of the variables $\boldsymbol{\sigma}^{2}$
#'
#' In this hypothetical example,
#' we assume that we have population data
#' and we are interested in the association between wages and education.
#' The regressor variable is years of education.
#' The regressand variable is hourly wage in US dollars.
#' The intercept is the predicted wage of an employee with 0 years of education.
#' The slope is the increase in hourly wage in US dollars
#' for one year increase in education.
#'
#' The the following vectors represent the parameters of the simple linear regression model.
#'
#' \begin{equation}
#'   \boldsymbol{\theta}_{\text{mean structure}}
#'   =
#'   \begin{bmatrix}
#'     \beta_1 \\
#'     \mu_x
#'   \end{bmatrix}
#' \end{equation}
#'
#' \begin{equation}
#'   \boldsymbol{\theta}_{\text{covariance structure}}
#'   =
#'   \begin{bmatrix}
#'     \beta_2 \\
#'     \sigma_{\varepsilon}^{2} \\
#'     \sigma_{x}^{2}
#'   \end{bmatrix}
#' \end{equation}
#'
#+ echo = FALSE
Variable <- c(
  "`beta1`",
  "`mux`"
)
Description <- c(
  "Intercept",
  "Mean of $x$"
)
Notation <- c(
  "$\\beta_1$",
  "$\\mu_x$"
)
Value <- c(
  beta1,
  mux
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Notation,
    Value
  ),
  caption = "$\\boldsymbol{\\theta}_{\\text{mean structure}}$",
  row.names = FALSE
)
#'
#'
#+ echo = FALSE
Variable <- c(
  "`beta2`",
  "`sigma2epsilon`",
  "`sigma2x`"
)
Description <- c(
  "Slope",
  "Variance of $\\varepsilon$",
  "Variance of $x$"
)
Notation <- c(
  "$\\beta_2$",
  "$\\sigma_{\\varepsilon}^{2}$",
  "$\\sigma_{x}^{2}$"
)
Value <- c(
  beta2,
  sigma2epsilon,
  sigma2x
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Notation,
    Value
  ),
  caption = "$\\boldsymbol{\\theta}_{\\text{covariance structure}}$",
  row.names = FALSE
)
#'
#'
#' ## Parameterization Including the Error Variance
#'
#' Variables are ordered as follows
#'
#' - $y$,
#' - $x$, and
#' - $\varepsilon$
#'
#' ### M Matrix
#'
#+
varnames <- c("y", "x", "epsilon")
k <- 2
q <- 1
p <- k + q
M <- matrix(
  data = c(
    beta1,
    mux,
    0
  ),
  ncol = 1
)
rownames(M) <- varnames
colnames(M) <- "M"
#'
#'
#+ echo = FALSE
knitr::kable(
  x = M
)
#'
#'
#' ### A Matrix
#'
#+
A <- matrix(
  data = c(
    0,
    0,
    0,
    beta2,
    0,
    0,
    1,
    0,
    0
  ),
  nrow = p
)
rownames(A) <- varnames
colnames(A) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = A
)
#'
#'
#' ### S Matrix
#'
#+
S <- matrix(
  data = c(
    0,
    0,
    0,
    0,
    sigma2x,
    0,
    0,
    0,
    sigma2epsilon
  ),
  nrow = p
)
rownames(S) <- varnames
colnames(S) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = S
)
#'
#'
#' ### I Matrix
#'
#+
I <- diag(nrow(A))
rownames(I) <- varnames
colnames(I) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = I
)
#'
#'
#' ### F Matrix
#'
#+
filter <- matrix(
  data = c(
    1,
    0,
    0,
    1,
    0,
    0
  ),
  nrow = 2,
  ncol = 3
)
rownames(filter) <- c("y", "x")
colnames(filter) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = filter
)
#'
#'
#' ### Model-Implied Mean Vector
#'
#+
result01_mutheta <- jeksterslabRsem::rammutheta(
  M = M,
  A = A,
  filter = filter
)
#'
#'
#+ echo = FALSE
knitr::kable(
  x = result01_mutheta
)
#'
#'
#' ### Model-Implied Variance-Covariance Matrix
#'
#+
result01_Sigmatheta <- jeksterslabRsem::ramSigmatheta(
  A = A,
  S = S,
  filter = filter
)
#'
#'
#+ echo = FALSE
knitr::kable(
  x = result01_Sigmatheta
)
#'
#'
#' ### M Matrix from mutheta and A Matrix
#'
#+
result01_M <- jeksterslabRsem::ramM(
  mu = mutheta,
  A = A,
  filter = filter
)
rownames(result01_M) <- varnames
colnames(result01_M) <- "M"
#'
#'
#+ echo = FALSE
knitr::kable(
  x = result01_M
)
#'
#'
#' ### S Matrix from sigma2 and A Matrix
#'
#' This should only be used when the off-diagonal elements
#' of the $\mathbf{S}$ matrix are all zeroes.
#'
#+
result01_S <- jeksterslabRsem::ramsigma2(
  sigma2 = c(sigma2y, sigma2x, sigma2epsilon),
  A = A,
  start = FALSE
)
rownames(result01_S) <- varnames
colnames(result01_S) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = result01_S
)
#'
#'
#' ## Parameterization Excluding the Error Variance
#'
#' Variables are ordered as follows
#'
#' - $y$, and
#' - $x$
#'
#' ### M Matrix
#'
#+
varnames <- c("y", "x")
k <- 2
q <- 0
p <- k + q
M <- matrix(
  data = c(
    beta1,
    mux
  ),
  ncol = 1
)
rownames(M) <- varnames
colnames(M) <- "M"
#'
#'
#+ echo = FALSE
knitr::kable(
  x = M
)
#'
#'
#' ### A Matrix
#'
#+
A <- matrix(
  data = c(
    0,
    0,
    beta2,
    0
  ),
  nrow = p
)
rownames(A) <- varnames
colnames(A) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = A
)
#'
#'
#' ### S Matrix
#'
#+
S <- matrix(
  data = c(
    sigma2epsilon,
    0,
    0,
    sigma2x
  ),
  nrow = p
)
rownames(S) <- varnames
colnames(S) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = S
)
#'
#'
#' ### I Matrix
#'
#+
I <- diag(nrow(A))
rownames(I) <- varnames
colnames(I) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = I
)
#'
#'
#' ### F Matrix
#'
#+
filter <- diag(nrow(A))
rownames(filter) <- c("y", "x")
colnames(filter) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = filter
)
#'
#'
#' ### Model-Implied Mean Vector
#'
#+
result02_mutheta <- jeksterslabRsem::rammutheta(
  M = M,
  A = A,
  filter = filter
)
#'
#'
#+ echo = FALSE
knitr::kable(
  x = result02_mutheta
)
#'
#'
#' ### Model-Implied Variance-Covariance Matrix
#'
#+
result02_Sigmatheta <- jeksterslabRsem::ramSigmatheta(
  A = A,
  S = S,
  filter = filter
)
#'
#'
#+ echo = FALSE
knitr::kable(
  x = result02_Sigmatheta
)
#'
#'
#' ### M Matrix from mutheta and A Matrix
#'
#+
result02_M <- jeksterslabRsem::ramM(
  mu = mutheta,
  A = A,
  filter = filter
)
rownames(result02_M) <- varnames
colnames(result02_M) <- "M"
#'
#'
#+ echo = FALSE
knitr::kable(
  x = result02_M
)
#'
#'
#' ### S Matrix from sigma2 and A Matrix
#'
#' This should only be used when the off-diagonal elements
#' of the $\mathbf{S}$ matrix are all zeroes.
#'
#+
result02_S <- jeksterslabRsem::ramsigma2(
  sigma2 = c(sigma2y, sigma2x),
  A = A,
  start = FALSE
)
rownames(result02_S) <- varnames
colnames(result02_S) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = result02_S
)
#'
#'
#' ## Parameterization Excluding the Error Variance
#'
#' Variables are ordered as follows
#'
#' - $x$, and
#' - $y$
#'
#' ### M Matrix
#'
#+
varnames <- c("x", "y")
k <- 2
q <- 0
p <- k + q
M <- matrix(
  data = c(
    mux,
    beta1
  ),
  ncol = 1
)
rownames(M) <- varnames
colnames(M) <- "M"
#'
#'
#+ echo = FALSE
knitr::kable(
  x = M
)
#'
#'
#' ### A Matrix
#'
#+
A <- matrix(
  data = c(
    0,
    beta2,
    0,
    0
  ),
  nrow = p
)
rownames(A) <- varnames
colnames(A) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = A
)
#'
#'
#' ### S Matrix
#'
#+
S <- matrix(
  data = c(
    sigma2x,
    0,
    0,
    sigma2epsilon
  ),
  nrow = p
)
rownames(S) <- varnames
colnames(S) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = S
)
#'
#'
#' ### I Matrix
#'
#+
I <- diag(nrow(A))
rownames(I) <- varnames
colnames(I) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = I
)
#'
#'
#' ### F Matrix
#'
#+
filter <- diag(nrow(A))
rownames(filter) <- c("x", "y")
colnames(filter) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = filter
)
#'
#'
#' ### Model-Implied Mean Vector
#'
#+
result03_mutheta <- jeksterslabRsem::rammutheta(
  M = M,
  A = A,
  filter = filter
)
#'
#'
#+ echo = FALSE
knitr::kable(
  x = result03_mutheta
)
#'
#'
#' ### Model-Implied Variance-Covariance Matrix
#'
#+
result03_Sigmatheta <- jeksterslabRsem::ramSigmatheta(
  A = A,
  S = S,
  filter = filter
)
#'
#'
#+ echo = FALSE
knitr::kable(
  x = result03_Sigmatheta
)
#'
#'
#' ### M Matrix from mutheta and A Matrix
#'
#+
result03_M <- jeksterslabRsem::ramM(
  mu = mutheta,
  A = A,
  filter = filter
)
rownames(result03_M) <- varnames
colnames(result03_M) <- "M"
#'
#'
#+ echo = FALSE
knitr::kable(
  x = result03_M
)
#'
#'
#' ### S Matrix from sigma2 and A Matrix
#'
#' This should only be used when the off-diagonal elements
#' of the $\mathbf{S}$ matrix are all zeroes.
#'
#+
result03_S <- jeksterslabRsem::ramsigma2(
  sigma2 = c(sigma2x, sigma2y),
  A = A,
  start = TRUE
)
rownames(result03_S) <- varnames
colnames(result03_S) <- varnames
#'
#'
#+ echo = FALSE
knitr::kable(
  x = result03_S
)
#'
#'
#+ echo = FALSE
context("Test simple-regression-ram")
test_that("muy", {
  expect_equivalent(
    muy,
    result01_mutheta["y", 1],
    result02_mutheta["y", 1],
    result03_mutheta["y", 1]
  )
})
test_that("mux", {
  expect_equivalent(
    mux,
    result01_mutheta["x", 1],
    result02_mutheta["x", 1],
    result03_mutheta["x", 1]
  )
})
test_that("M - beta1", {
  expect_equivalent(
    beta1,
    result01_M["y", 1],
    result02_M["y", 1],
    result03_M["y", 1]
  )
})
test_that("M - mux", {
  expect_equivalent(
    mux,
    result01_M["x", 1],
    result02_M["x", 1],
    result03_M["x", 1]
  )
})
test_that("M - epsilon", {
  expect_equivalent(
    0,
    result01_M["epsilon", 1]
  )
})
test_that("sigma2x", {
  expect_equivalent(
    sigma2x,
    result01_Sigmatheta["x", "x"],
    result02_Sigmatheta["x", "x"],
    result03_Sigmatheta["x", "x"]
  )
})
test_that("sigma2y", {
  expect_equivalent(
    sigma2y,
    result01_Sigmatheta["y", "y"],
    result02_Sigmatheta["y", "y"],
    result03_Sigmatheta["y", "y"]
  )
})
test_that("sigmayx", {
  expect_equivalent(
    sigmayx,
    result01_Sigmatheta["x", "y"],
    result02_Sigmatheta["x", "y"],
    result03_Sigmatheta["x", "y"],
    result01_Sigmatheta["y", "x"],
    result02_Sigmatheta["y", "x"],
    result03_Sigmatheta["y", "x"]
  )
})
test_that("S - xy", {
  expect_equivalent(
    0,
    result01_S["x", "y"],
    result02_S["x", "y"],
    result03_S["x", "y"],
    result01_S["y", "x"],
    result02_S["y", "x"],
    result03_S["y", "x"]
  )
})
test_that("S - y", {
  expect_equivalent(
    sigma2epsilon,
    result01_S["y", "y"],
    result02_S["y", "y"],
    result03_S["y", "y"]
  )
})
test_that("S - x", {
  expect_equivalent(
    sigma2x,
    result01_S["x", "x"],
    result02_S["x", "x"],
    result03_S["x", "x"]
  )
})
#'
