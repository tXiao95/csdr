#' Example data for causal sufficient dimension reduction
#'
#' A small simulated dataset for running the package examples. It contains a
#' continuous outcome, two continuous exposures, and three continuous
#' covariates. The exposure and covariate matrices have informative column and
#' row names.
#'
#' @format A named list with three components:
#' \describe{
#'   \item{Y}{A numeric outcome vector with 60 observations.}
#'   \item{A}{A numeric 60 by 2 exposure matrix.}
#'   \item{C}{A numeric 60 by 3 covariate matrix.}
#' }
#' @source Simulated for the `csdr` package with seed 2026.
"csdr_example"
