#' ---
#' title: "Test: Reticular Action Model"
#' author: "Ivan Jacob Agaloos Pesigan"
#' date: "`r Sys.Date()`"
#' output: rmarkdown::html_vignette
#' vignette: >
#'   %\VignetteIndexEntry{Test: Reticular Action Model}
#'   %\VignetteEngine{knitr::rmarkdown}
#'   %\VignetteEncoding{UTF-8}
#' ---
#'
#+ knitr_options, include=FALSE, cache=FALSE
knitr::opts_chunk$set(
  error = TRUE,
  collapse = TRUE,
  comment = "#>",
  out.width = "100%"
)
#'
#+ setup, echo = FALSE, message = FALSE
library(testthat)
library(lavaan)
library(MASS)
library(jeksterslabRsem)
context("Test Reticular Action Model.")
#'
#' ## Set test parameters
#'
#' \begin{equation}
#'   Y_i = \delta_Y + \tau^{\prime} X_i + \beta M_i + \varepsilon_{Y_{i}}
#' \end{equation}
#'
#' \begin{equation}
#'   M_i = \delta_M + \alpha X_i + \varepsilon_{M_{i}}
#' \end{equation}
#'
#' ```{tikz, simple_mediation, echo = FALSE, fig.cap = "The Simple Mediation Model", fig.ext = "png", cache = TRUE}
#' \usetikzlibrary{er, arrows, positioning}
#' \begin{tikzpicture}[
#'   auto,
#'   node distance = 20mm,
#'   manifest/.style = {
#'     rectangle,
#'     draw,
#'     thick,
#'     inner sep = 0pt,
#'     minimum width = 15mm,
#'     minimum height = 10mm
#'   },
#'   inv/.style = {
#'     rectangle,
#'     draw=none,
#'     fill=none,
#'     inner sep = 0pt,
#'     minimum width = 15mm,
#'     minimum height = 10mm
#'   },
#'   error/.style = {
#'     ellipse,
#'     draw,
#'     thick,
#'     inner sep = 0pt,
#'     minimum size = 7mm,
#'     align = center
#'   },
#'   mean/.style={
#'     regular polygon,
#'     regular polygon sides = 3,
#'     draw,
#'     thick,
#'     inner sep = 0pt,
#'     minimum width = 7mm,
#'     minimum height = 7mm
#'   },
#'   path/.style = {
#'     ->,
#'     thick,
#'     >=stealth'
#'   },
#'   cov/.style = {
#'     <->,
#'     thick,
#'     >=stealth'
#'   },
#' ]
#'   \node[manifest] (X) {$X$};
#'   \node[mean] (1) [above right = of X] {$1$};
#'   \node[manifest] (M) [above = of 1] {$M$};
#'   \node[manifest] (Y) [below right = of 1] {$Y$};
#'   \node[error] (epsilon_M) [right = of M] {$\epsilon_M$};
#'   \node[error] (epsilon_Y) [right = of Y] {$\epsilon_Y$};
#'   \draw [path] (X) to node {$\tau^{\prime}$} (Y);
#'   \draw [path] (X) to node {$\alpha$} (M);
#'   \draw [path] (M) to node {$\beta$} (Y);
#'   \draw [path] (epsilon_M) to node {$1$} (M);
#'   \draw [path] (epsilon_Y) to node {$1$} (Y);
#'   \draw [path] (1) to node {$\mu_X$} (X);
#'   \draw [path] (1) to node {$\delta_M$} (M);
#'   \draw [path] (1) to node {$\delta_Y$} (Y);
#'   \draw [cov] (X) to[out=170,in=190,looseness=5] node[left] {$\sigma^{2}_{X}$} (X);
#'   \draw [cov] (epsilon_M) to[out=70,in=110,looseness=5] node[above] {$\sigma^{2}_{\epsilon_{M}}$} (epsilon_M);
#'   \draw [cov] (epsilon_Y) to[out=70,in=110,looseness=5] node[above] {$\sigma^{2}_{\epsilon_{Y}}$} (epsilon_Y);
#'   \draw [cov] (1) to[out=-60,in=-120,looseness=7] node[below] {1} (1);
#' \end{tikzpicture}
#' ```
#'
#+ parameters, echo = FALSE
var_names <- c("X", "M", "Y")
slopes <- runif(
  n = 3,
  min = .10,
  max = .50
)
tau_prime <- slopes[1]
beta <- slopes[2]
alpha <- slopes[3]
sigma2 <- runif(
  n = 3,
  min = 10^2,
  max = 15^2
)
sigma2X <- sigma2[1]
sigma2M <- sigma2[2]
sigma2Y <- sigma2[3]
mu <- runif(
  n = 3,
  min = 5^2,
  max = 10^2
)
muX <- mu[1]
muM <- mu[2]
muY <- mu[3]
sigma2epsilonM <- sigma2M - alpha^2 * sigma2X
sigma2epsilonY <- sigma2Y - (beta^2 * alpha^2 * sigma2X) - (beta^2 * sigma2epsilonM) - (2 * alpha * beta * tau_prime * sigma2X) - (tau_prime^2 * sigma2X)
delta_M <- muM - alpha * muX
delta_Y <- muY - tau_prime * muX - beta * muM
A <- matrix(
  data = c(
    0,
    alpha,
    tau_prime,
    0,
    0,
    beta,
    0,
    0,
    0
  ),
  ncol = 3
)
colnames(A) <- var_names
rownames(A) <- var_names
S <- matrix(
  data = c(
    sigma2X,
    0,
    0,
    0,
    sigma2epsilonM,
    0,
    0,
    0,
    sigma2epsilonY
  ),
  ncol = 3
)
colnames(S) <- var_names
rownames(S) <- var_names
F <- diag(3)
colnames(F) <- var_names
rownames(F) <- var_names
I <- diag(3)
colnames(I) <- var_names
rownames(I) <- var_names
theta <- c(
  tau_prime = tau_prime,
  beta = beta,
  alpha = alpha,
  sigma2X = sigma2X,
  sigma2epsilonM = sigma2epsilonM,
  sigma2epsilonY = sigma2epsilonY,
  muX = muX,
  delta_M = delta_M,
  delta_Y = delta_Y
)
Variable <- c(
  "`tau_prime`",
  "`beta`",
  "`alpha`",
  "`sigma2X`",
  "`sigma2epsilonM`",
  "`sigma2epsilonY`",
  "`muX`",
  "`delta_M`",
  "`delta_Y`"
)
Description <- c(
  "Path from $X$ to $Y$ ($\\tau^{\\prime}$).",
  "Path from $M$ to $Y$ ($\\beta$).",
  "Path from $X$ to $M$ ($\\alpha$).",
  "Variance of $X$ ($\\sigma^{2}_{X}$)",
  "Residual variance of $\\epsilon_{M}$ ($\\sigma^{2}_{\\varepsilon_{M}}$)",
  "Residual variance of $\\epsilon_{Y}$ ($\\sigma^{2}_{\\varepsilon_{Y}}$)",
  "Mean of $X$ ($\\mu_X$)",
  "Intercept for model predicting $M$ ($\\delta_M$)",
  "Intercept for model predicting $Y$ ($\\delta_Y$)"
)
Value <- c(
  tau_prime,
  beta,
  alpha,
  sigma2X,
  sigma2epsilonM,
  sigma2epsilonY,
  muX,
  delta_M,
  delta_Y
)
knitr::kable(
  x = data.frame(
    Variable,
    Description,
    Value
  ),
  row.names = FALSE,
  caption = "Parameters"
)
knitr::kable(
  x = A,
  row.names = TRUE,
  caption = "A Matrix"
)
knitr::kable(
  x = S,
  row.names = TRUE,
  caption = "S Matrix"
)
knitr::kable(
  x = F,
  row.names = TRUE,
  caption = "F Matrix"
)
knitr::kable(
  x = I,
  row.names = TRUE,
  caption = "I Matrix"
)
#'
#' ## Model-Implied Matrices
#'
#+ ram
Sigmatheta <- ram_Sigmatheta(
  A = A,
  S = S,
  F = F,
  I = I
)
Sigmatheta2 <- ram_S(
  A = A,
  sigma2 = sigma2,
  F = F,
  I = I,
  SigmaMatrix = TRUE
)
S_ver2 <- ram_S(
  A = A,
  sigma2 = sigma2,
  F = F,
  I = I,
  SigmaMatrix = FALSE
)
M <- ram_M(
  A = A,
  F = F,
  I = I,
  mu = mu
)
mutheta <- ram_mutheta(
  A = A,
  F = F,
  I = I,
  M = M
)
M <- as.vector(M)
names(M) <- var_names
#'
#' ## Generate Data with $\hat{\boldsymbol{\Sigma}}$ equal to $\boldsymbol{\Sigma}$
#'
#+ data
data <- mvrnorm(
  n = 1000,
  mu = mu,
  Sigma = Sigmatheta,
  empirical = TRUE
)
#'
#' ## Estimate Parameters from Data
#'
#+ estimate
model <- "
  Y ~ tau_prime*X + beta*M
  M ~ alpha*X
  X ~~ sigma2X*X
  M ~~ sigma2epsilonM*M
  Y ~~ sigma2epsilonY*Y
  X ~ muX*1
  M ~ delta_M*1
  Y ~ delta_Y*1
"
fit <- sem(
  model = model,
  data = data,
  likelihood = "wishart"
)
thetahat <- as.vector(coef(fit))
muhat <- as.vector(colMeans(data))
#'
#' ## Results
#'
#+ results, echo = FALSE
knitr::kable(
  x = Sigmatheta,
  row.names = TRUE,
  caption = "Model-implied variance-covariance matrix"
)
knitr::kable(
  x = mutheta,
  row.names = TRUE,
  col.names = "$\\boldsymbol{\\mu}$",
  caption = "Model-implied mean vector"
)
knitr::kable(
  x = data.frame(
    Item = c(
      Description,
      "Mean of $X$ ($\\mu_X$)",
      "Mean of $M$ ($\\mu_M$)",
      "Mean of $Y$ ($\\mu_Y$)"
    ),
    Parameters = c(
      theta,
      mu
    ),
    Estimates = c(
      thetahat,
      muhat
    )
  ),
  row.names = FALSE,
  caption = "Results"
)
#'
#' ## testthat
#'
#+ testthat_01
test_that("theta and thetahat are equivalent", {
  expect_equivalent(
    round(
      x = theta,
      digits = 2
    ),
    round(
      x = thetahat,
      digits = 2
    )
  )
})
#'
#+ testthat_02
test_that("mu, mutheta, and muhat are equivalent", {
  expect_equivalent(
    round(
      x = mu,
      digits = 2
    ),
    round(
      x = mutheta,
      digits = 2
    ),
    round(
      x = muhat,
      digits = 2
    )
  )
})
#'
#+ testthat_03
test_that("ram_M returns the correct values", {
  expect_equivalent(
    round(
      x = theta[7:9],
      digits = 2
    ),
    round(
      x = thetahat[7:9],
      digits = 2
    ),
    round(
      x = M,
      digits = 2
    )
  )
})
#'
#+ testthat_04
test_that("muX", {
  expect_equivalent(
    round(
      x = muX,
      digits = 2
    ),
    round(
      x = M[1],
      digits = 2
    )
  )
})
#'
#+ testthat_05
test_that("muM", {
  expect_equivalent(
    round(
      x = muM,
      digits = 2
    ),
    round(
      x = M[2] + alpha * muX,
      digits = 2
    )
  )
})
#'
#+ testthat_06
test_that("muY", {
  expect_equivalent(
    round(
      x = muY,
      digits = 2
    ),
    round(
      x = M[3] + tau_prime * muX + beta * muM,
      digits = 2
    )
  )
})
#'
#+ testthat_07
test_that("S", {
  expect_equivalent(
    round(
      x = S,
      digits = 2
    ),
    round(
      x = S_ver2,
      digits = 2
    )
  )
})
#'
#+ testthat_08
test_that("Sigma", {
  expect_equivalent(
    round(
      x = Sigmatheta,
      digits = 2
    ),
    round(
      x = Sigmatheta2,
      digits = 2
    ),
    round(
      x = cov(data),
      digits = 2
    )
  )
})
