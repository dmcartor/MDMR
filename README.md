# MDMR
Conduct multivariate distance matrix regression (Anderson, 2001; McArdle & Anderson, 2001) using the analytical p-values described by McArtor, Lubke, \& Bergeman (second revision under review) and compute permutation-based measures of effect size proposed by McArtor et al.

The most recent version can be installed from github using the devtools package:

    devtools::install_github("dmcartor/mdmr", build_vignettes = TRUE)
    library(MDMR)

or directly from CRAN:

    install.packages(MDMR)
    library(MDMR)


## Usage

There are two primary functions that comprise this package: mdmr(), which regresses a distance matrix onto a set of predictors, and delta(), which computes measures of univariate effect size in the context of multivariate distance matrix regression. The help files of both functions provide more general information than the package vignette. 

For a complete illustration of the package, see the package vignette by running the following line in the R console:

    vignette('mdmr-vignette')

## Example

    # Import data
    data(mdmrdata)
    
    # Compute distance matrix
    D <- dist(Y.mdmr, method = 'euclidean')
    
    # Conduct MDMR
    mdmr.res <- mdmr(X = X.mdmr, D = D)
    
    # Check results
    summary(mdmr.res)
    
    # Study univariate effect sizes
    delta.res <- delta(X.mdmr, Y = Y.mdmr, dtype = 'euclidean', 
                      niter = 10, plot.res = T)
