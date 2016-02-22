# MDMR 0.3.2

* Added two arguments to delta() that allow users to specify which effect sizes they are interested in measuring with respect to both predictors (x.inds) and outcomes (y.inds). This was primarily added for situations where there are many predictors, not all of which are found to be statistically significant, and the user is only interested in effect sizes of the significant predictors. Note, however, that the conditional effect sizes reported by delta() are still conditioned on all variables comprising the matrix of predictors even if only a subset of effect sizes are requested to be computed and reported.

# MDMR 0.3.1

* Substantial improvement to computation time
* Updated delta() with an option to plot results in grayscale

# MDMR 0.3.0

* Added an option to compute permutation-based p-values to mdmr()
* Fixed an error resulting from passing a univariate predictor to delta()

# MDMR 0.2.1

* Reduced computation time of examples to comply with CRAN guidelines.

# MDMR 0.2.0

* First stable release
