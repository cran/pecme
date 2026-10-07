# pecme

**Penalized ECME Estimation for Censored Linear Mixed Models**

`pecme` fits Gaussian random-intercept linear mixed models

<p align="center">

<em>Y\*<sub>ij</sub> = x<sub>ij</sub>'β + b<sub>j</sub> + e<sub>ij</sub>,   e<sub>ij</sub> \~ N(0, σ<sup>2</sup>),   b<sub>j</sub> \~ N(0, σ<sub>b</sub><sup>2</sup>)</em>

</p>

where the response `Y*` may be left-, right-, or interval-censored, and/or completely missing, using the **ECME algorithm** (Liu & Rubin, 1994), with **simultaneous penalized variable selection** for `beta`: `none`, `lasso`, `adaptive` (Adaptive Lasso), `scad`, `mcp`, `elastic` (Elastic Net), `ridge`.

- Random intercept integrated out **numerically** by fixed-order Gauss-Hermite quadrature (a numerical approximation, not an exact evaluation of the integral -- use `pecme_quadrature_check()` to check sensitivity to the number of nodes).
- `beta` updated by penalized coordinate descent (a true CM-step on the expected complete-data objective).
- `(sigma, sigma_b)` updated, by default, by directly maximizing the **actual observed-data marginal likelihood** (not merely its EM lower bound) -- this is the "either" in ECME, and what distinguishes it from plain EM/ECM.
- A single fit `pecme()` or a full penalty-parameter grid search ([`grid_pecme()`](R/grid-search.R)), sequential or parallel (`workers = <n_cores>`), with a user-supplied or data-dependent lambda grid and selection by **AIC/BIC/EBIC/GCV or cross-validation** (never by raw log-likelihood).
- Missing **response** values are estimated via the same E-step used for censoring (`missing_y = TRUE`); missing **covariates** have documented, practical handling (`missing_x`), with a Rubin's-rules combiner ([`pecme_mi_combine()`](R/missing-data.R)) for rigorous multiple imputation workflows.

## Installation

``` r
# from a local checkout / the unzipped source tree:
install.packages("devtools")
devtools::document("path/to/pecme")   # regenerates NAMESPACE/man from
                                       # the roxygen comments already in R/
devtools::install("path/to/pecme")
```

The package should be checked locally with `R CMD check`, and the test suite can be run with `testthat::test_dir("tests/testthat")` before use in a final analysis.

## Quick start

``` r
library(pecme)

sim <- simulate_pecme_data(
  n_groups = 40, n_per_group = 10, p = 12, n_active = 4,
  beta_active = 1.8, sigma = 0.6, sigma_b = 0.5,
  censor_type = "left", censor_prob = 0.2, seed = 1
)

fixed <- reformulate(paste0("x", 1:12), response = "y")

# one fit at a fixed lambda
fit <- pecme(fixed, random_var = "group", data = sim$data,
             censor = "censor", left = "left", right = "right",
             penalty = "lasso", lambda = 0.3)
summary(fit)

# grid search, BIC-selected, parallel across lambda on a 32-core machine
g <- grid_pecme(fixed, random_var = "group", data = sim$data,
                 censor = "censor", left = "left", right = "right",
                 penalty = "scad", nlambda = 40, criterion = "bic",
                 workers = 32)
plot(g)
```

`response` is not a separate argument: `pecme()`/`grid_pecme()` read it from the left-hand side of `fixed_formula` (here, `y`).

See `vignette("pecme-intro")` for a full walk-through (single fits, grid search, missing response, missing covariates).

## Diagnostics and valid inference for penalized fits

- `pecme_quadrature_check(fit)` -- Gauss-Hermite quadrature is a numerical approximation, not exact integration; this reports how sensitive the fit is to the number of nodes.
- `pecme_effective_df(fit)` -- a closed-form effective degrees-of- freedom for Ridge (the naive nonzero-coefficient count used for AIC/BIC/GCV is particularly misleading there, since Ridge rarely zeroes anything); `NA` with an explanation for Elastic Net/SCAD/MCP, where no closed form is implemented.
- `pecme_bootstrap(fit)` -- a cluster (group-level) bootstrap: refits the same penalized model with groups resampled with replacement, for a genuinely valid source of SEs/CIs where naive ones do not apply.
- `fit$refit` -- an automatic "relaxed refit" (unpenalized, on the selected predictors only) for approximate SEs/p-values; reduces shrinkage bias but is **not** valid unconditional post-selection inference (the selection step itself is not accounted for).
- `pecme_mi_combine(fits, use_refit = TRUE)` -- pools relaxed refits across multiply-imputed data sets via Rubin's rules; requires `use_refit = TRUE` for penalized fits (naive SEs are `NA` by design and cannot be pooled), with the same post-selection-inference caveat.

## Thesis-reporting metrics

`pecme_metrics()` collects everything you'd typically report for a mixed-model fit in one call: MSE, RMSE, MAE, MAPE, SMAPE, a predictive (OLS-style) R-squared *and* the marginal/conditional pseudo-R-squared that are the standard for mixed models (Nakagawa & Schielzeth, 2013), logLik, AIC, BIC, EBIC, and GCV.

``` r
pecme_metrics(fit)                                   # in-sample
pecme_metrics(fit, newdata = test_df)                # out-of-sample (response column auto-detected)
pecme_metrics(fit, loo = TRUE, loo_folds = 10)        # + leave-one-group-out CV
```

**On LOO and WAIC**: `pecme` is a frequentist penalized-likelihood estimator -- it does not produce posterior draws, so a literal WAIC (or DIC) is not a well-defined quantity here, and `pecme_metrics()` deliberately reports `$WAIC` as `NA` with an explanatory note rather than fabricate one. What it *does* provide, via `loo = TRUE`, is a genuine out-of-sample criterion in the same spirit: a leave-one-group- out (or leave-one-fold-out) cross-validated predictive log-likelihood, computed by actually refitting the model with each group held out -- this is the frequentist counterpart of what WAIC and PSIS-LOO both aim to approximate, at the cost of being more expensive to compute.

## Review-sensitive implementation notes

The current implementation explicitly addresses the review points that affect statistical validity and reproducibility: censoring thresholds are tied to the requested `censor_prob`; the outer penalized objective uses the correct minimization monotonicity direction; failed variance optimization cannot create artificial zero likelihood/objective values; the `"random"` initialization changes values actually passed to the optimizer; transformed responses are taken from the evaluated model frame; refit utilities preserve the original analysis data and fitting arguments; EBIC uses the number of selected penalized predictors; and convergence/inner-optimizer diagnostics are retained in the top-level fitted object.

`pecme_effective_df()` is deliberately described as a working Ridge-type/conditional approximation where a full effective-df result has not been derived for the complete censored mixed-effects estimator. Likewise, `pecme_lambda_max()` is a data-dependent reference scale, not a formally established exact sparsity threshold.

## Package layout

```         
R/
  utils-math.R        Gauss-Hermite quadrature (cached), log-sum-exp
  truncation.R         Truncated-normal moments & log-lik (incl. "missing")
  estep.R               E-step, sequential and cluster-parallel over groups
  loglik.R              Observed-data marginal log-likelihood
  penalty.R              Penalty values + closed-form coordinate updates
  coordinate-descent.R    Penalized CM-step for beta
  objective.R              Full penalized objective / Q-function
  variance-update.R         EM and true-ECME updates for sigma, sigma_b
  standardize.R              Design-matrix standardization
  random-effects-init.R       lmer()-based starting values
  missing-data.R                Missing y / missing X handling, MI combiner
  validate.R                     Input validation
  prepare.R                       Formula/data -> matrices (no listwise deletion)
  engine.R                         Core ECME loop for one (penalty, lambda)
  lambda-grid.R                     Data-dependent lambda_max / grid
  parallel-utils.R                   Cluster setup/teardown
  fit.R                               pecme_setup()/pecme_run()/pecme()
  grid-search.R                       grid_pecme()
  metrics.R                            AIC/BIC/EBIC/GCV, prediction metrics, effective df
  thesis-metrics.R                      pecme_metrics(), pecme_loo()
  quadrature-check.R                    pecme_quadrature_check()
  pecme-bootstrap.R                     pecme_bootstrap() (cluster bootstrap)
  methods.R                             print/summary/coef/predict/plot/...
  simulate.R                             simulate_pecme_data()
tests/testthat/        Unit tests
vignettes/pecme-intro.Rmd  Worked example
```

## License

MIT
