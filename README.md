# csdr

`csdr` estimates causal sufficient dimension reduction targets for continuous
multivariate exposures. The main user-facing entry point is `csdr()`.

## Citation

The methods implemented in `csdr` are introduced in:

> Hsiao, T. W., Chang, H. H., and Nabi, R. (2026). *Causal Sufficient
> Dimension Reduction for Multiple Continuous Exposures with an Application to
> Environmental Mixtures*. arXiv:2606.14840.
> https://doi.org/10.48550/arXiv.2606.14840

If you use `csdr` in your work, please cite this preprint:

```bibtex
@misc{hsiao2026causal,
  title = {Causal Sufficient Dimension Reduction for Multiple Continuous
           Exposures with an Application to Environmental Mixtures},
  author = {Hsiao, Thomas W. and Chang, Howard H. and Nabi, Razieh},
  year = {2026},
  eprint = {2606.14840},
  archivePrefix = {arXiv},
  primaryClass = {stat.ME},
  doi = {10.48550/arXiv.2606.14840},
  url = {https://arxiv.org/abs/2606.14840}
}
```

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

Inspect the original MAVE fit and, when dimension selection was run, its raw
selection object:

```r
mave <- mave_fits(fit, variant = "DR")
mave$fit
mave$dimension_selection
```

MAVE objects are retained by default. Nuisance fits can be substantially larger
and are retained only when requested:

```r
fit_with_nuisance <- csdr(
  Y = csdr_example$Y,
  A = csdr_example$A,
  C = csdr_example$C,
  d = 1,
  L = 2,
  learners = csdr_learners(sl_library = "SL.glm"),
  keep_nuisance = TRUE,
  verbose = FALSE
)

outcome_folds <- nuisance_fits(fit_with_nuisance, role = "outcome")
gps_fold_1 <- nuisance_fits(fit_with_nuisance, role = "gps", fold = 1)
```

Set `keep_mave = FALSE`, `keep_targets = FALSE`, or leave
`keep_nuisance = FALSE` to reduce the size of saved fit objects when those raw
components are not needed.

The quick start uses `SL.glm`, which is included with SuperLearner and keeps the
example lightweight. The default SuperLearner library can also use `glmnet`,
`xgboost`, and `earth` when those optional packages are installed.

## Which Variant Should I Use?

| Variant | Nuisance models | Main modeling reliance | Relative cost |
|---|---|---|---|
| `RA` | Outcome regression | Correct outcome regression | Low |
| `DR` | Outcome regression and GPS | One of the outcome regression or GPS is correct | Medium |
| `PO` | Outcome regression and GPS | Same as DR | Medium |
| `RP` | `E[Y \| C]` and each `E[A_j \| C]` | Residualization models | High |

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

## Outcome Families

Gaussian is the default outcome family and is used for continuous and count
outcomes. For a binary outcome, configure the package-provided SuperLearner
nuisance regressions with `outcome_family = "binomial"`:

```r
binary_learners <- csdr_learners(
  outcome_family = "binomial",
  sl_library = "SL.glm"
)

binary_fit <- csdr(
  Y = binary_Y,
  A = A,
  C = C,
  learners = binary_learners
)
```

The selected family applies to `E[Y | A, C]` for RA, DR, and PO and to
`E[Y | C]` for RP. Regressions of `A` on `C` remain Gaussian, and the GPS is
unchanged. Explicit `outcome` or `rp_y` learner specifications override this
shortcut with a warning. Custom learners remain responsible for their own
response-family behavior. Poisson-specific nuisance modeling and automatic
family inference are not currently provided.

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

## Troubleshooting

- **Input errors before fitting:** Check the argument and column names reported
  by the error. `Y`, `A`, and `C` must be finite numeric data with compatible
  rows; `A` and `C` need unique, nonoverlapping names and nonconstant columns.
- **Nuisance failure in a fold:** The error reports the learner role and fold.
  Start with `SL.glm`, verify the role-specific learner contract, and inspect
  sparse or constant training-fold inputs.
- **MAVE fitting failure:** Inspect the generated target with `targets()` and
  review `mave_control`. If automatic dimension selection fails, consider a
  scientifically justified explicit `d`.
- **Raw object unavailable:** Refit with `keep_mave = TRUE` or
  `keep_nuisance = TRUE`, as identified by the accessor error. Retaining
  nuisance fits can substantially increase saved-object size.

Package errors have structured classes. Validation errors inherit from
`csdr_validation_error`; fitting errors inherit from `csdr_fit_error` and retain
the original error in `condition$parent`. When available, `stage`, `variant`,
`fold`, and `role` identify where the failure occurred.
