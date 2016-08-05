# MDMR 0.4.1

* Updated documentation of the package to reflect the updated status of the manuscript upon which this package is based, as well as the addition of a new coauthor to the paper. At the request of the editor, we have conducted a second round of minor revisions, and the mew version is currently under review. 

* Due to the fact that many of the mathematical details of this package rely on access to the not-yet-published paper, we encourage interested users to email us with questions regarding the methodology underlying this package.

# MDMR 0.4.0

* Substantially changed how categorical predictors are treated. In the new version, omnibus effects are tested rather than individually testing dummy/contrast codes.

# MDMR 0.3.3

* Fixed a bug that occurred with missingness on X
* Added significance indicators to summary.mdmr() that mirror those used by lm()

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
