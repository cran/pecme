# pecme 0.1.1

* Hardened grid search against failed/nonfinite candidates and restricted EBIC recomputation to `criterion = "ebic"`.
* Fixed LOO handling of all-failed folds, required converged fold refits, and recomputed Adaptive Lasso pilots/weights within each training fold.
* Preserved `missing_x` behavior in bootstrap refits and repaired complete-case group remapping and regression-imputation intercept handling.
* Stored the full fixed-effect covariance matrix for valid Hessian-based fits and used \(x^T V_beta x\) for prediction standard errors.
* Information criteria now use the number of informative response rows when completely missing responses are present.
* Added explicit handling/documentation for transformed responses under censoring and for model-matrix missing-X imputation.
* `pecme_mi_combine(use_refit = TRUE)` now requires converged refits with valid finite standard errors.
