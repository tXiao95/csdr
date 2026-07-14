# csdr

`csdr` estimates causal sufficient dimension reduction targets for continuous
multivariate exposures. The main user-facing entry point is `csdr()`.

## Installation

```r
# install.packages("remotes")
remotes::install_github("tXiao95/csdr")
library(csdr)
```

## Quick Start

The package includes a small example dataset containing an outcome vector `Y`,
an exposure matrix `A`, and a covariate matrix `C`.

```r
data(csdr_example)

fit <- csdr(
  Y = csdr_example$Y,
  A = csdr_example$A,
  C = csdr_example$C,
  variants = "DR",
  d = 1,
  L = 5,
  seed = 1,
  learners = csdr_learners(sl_library = "SL.glm"),
  verbose = FALSE
)

fit
summary(fit)
```

Extract the estimated dimension-reduction directions, reduced exposure scores,
and generated target data:

```r
beta_hat <- coef(fit, variant = "DR")
Z_hat <- scores(fit, variant = "DR")
target_data <- targets(fit, variant = "DR")
```

The quick start uses `SL.glm`, which is included with SuperLearner and keeps the
example lightweight. The default SuperLearner library can also use `glmnet`,
`xgboost`, and `earth` when those optional packages are installed.

## Which Variant Should I Use?

| Variant | Nuisance models | Main modeling reliance | Relative cost |
|---|---|---|---|
| `RA` | Outcome regression | Correct outcome regression | Low |
| `DR` | Outcome regression and GPS | Either nuisance component provides robustness under the method's assumptions | Medium |
| `PO` | Outcome regression and GPS | Pseudo-outcome construction and marginalization | High |
| `RP` | `E[Y \| C]` and each `E[A_j \| C]` | Residualization models | Medium |

`DR` is the package default and a reasonable starting point when both the
outcome regression and generalized propensity score can be estimated. Use the
other variants when their assumptions or target construction better match the
analysis. See `?csdr` and `?csdr_target` for parameter-level details.

## Multiple Target Variants

`csdr()` can construct and fit multiple target variants in one call.

```r
fit_all <- csdr(
  Y = Y,
  A = A,
  C = C,
  variants = c("RA", "DR", "PO", "RP"),
  L = 5,
  seed = 1
)
```

For RA, DR, and PO, the generated exposure is the original `A`. For RP, the
generated exposure is the residualized exposure.

## Default Learners

Default nuisance learners:

- Outcome regression `E[Y | A, C]`: SuperLearner with `SL.glm`, `SL.glmnet`, `SL.xgboost`, `SL.earth`
- GPS `f(A | C)`: MVN GPS with linear conditional mean
- RP regression `E[Y | C]`: SuperLearner with the default library
- RP regression `E[A_j | C]`: SuperLearner with the default library

If optional learner packages are not installed, the SuperLearner fitting engine
uses the available requested wrappers.

## Modify the SuperLearner Library

Use one SuperLearner library for all default regression learner roles:

```r
fit_sl <- csdr(
  Y = Y,
  A = A,
  C = C,
  learners = csdr_learners(
    sl_library = c("SL.glm", "SL.ranger", "SL.xgboost")
  )
)
```

Use role-specific learner choices:

```r
fit_custom_sl <- csdr(
  Y = Y,
  A = A,
  C = C,
  learners = csdr_learners(
    outcome = sl_regression(SL.library = c("SL.glm", "SL.ranger")),
    gps = mvn_gps(method = "linear"),
    rp_y = sl_regression(SL.library = c("SL.glm")),
    rp_a = sl_regression(SL.library = c("SL.glm", "SL.ranger"))
  )
)
```

Use MVN GPS with a SuperLearner conditional mean:

```r
fit_gps_sl <- csdr(
  Y = Y,
  A = A,
  C = C,
  learners = csdr_learners(
    gps = mvn_gps(
      method = "SuperLearner",
      SL.library = c("SL.glm", "SL.ranger")
    )
  )
)
```

## Custom Learners

Custom regression learners should follow the contract `fitter(Y, W, ...)` and
return an object with a working `predict(object, newdata = ...)` method.

```r
my_fitter <- function(Y, W, ...) {
  fit <- stats::lm(Y ~ ., data = data.frame(Y = Y, W))
  structure(list(fit = fit), class = "my_regression")
}

predict.my_regression <- function(object, newdata, ...) {
  as.numeric(stats::predict(object$fit, newdata = as.data.frame(newdata)))
}

fit_custom <- csdr(
  Y = Y,
  A = A,
  C = C,
  learners = csdr_learners(
    outcome = custom_regression(my_fitter, label = "custom lm")
  )
)
```

Custom GPS learners should follow the contract `fitter(A, C, ...)` and return a
fitted density object whose `predict(object, newdata = ...)` method returns
numeric conditional density estimates.
