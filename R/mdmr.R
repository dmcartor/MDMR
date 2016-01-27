#' Compute Gower's centered similarity matrix from a distance matrix
#'
#' Compute Gower's centered similarity matrix \eqn{G}, which is the matrix
#' decomposed by the MDMR test statistic.
#'
#' @param d.mat Symmetric distance matrix (or R distance object) computed from
#'  the outcome data to be used in MDMR.
#'
#' @return G Gower's centered dissimilarity matrix computed from D.
#'
#' @author \packageAuthor{MDMR}
#'
#' @references Gower, J. C. (1966). Some distance properties of latent root and
#' vector methods used in multivariate analysis. Biometrika, 53(3-4), 325-338.
#'
#' @export
gower <- function(d.mat){
  # Convert distance object to matrix form
  d.mat <- as.matrix(d.mat)

  # Dimensionality of distance matrix
  n <- nrow(d.mat)

  # Create Gower's symmetry matrix (Gower, 1966)
  A <- -0.5 * d.mat^2

  # Subtract column means As = (I - 1/n * 11')A
  As <- A - rep(colMeans(A), rep.int(n, n))

  # Subtract row means G = As(I- 1/n * 11')
  # the transpose is just to deal with how "matrix minus vector" operations work
  # in R
  return(t(As) - rep(rowMeans(As), rep.int(n, n)))
}





#' Conduct MDMR with analytical p-values
#'
#' \code{mdmr} (multivariate distance matrix regression) is used to regress a
#' distance matrix onto a set of predictors. It returns the test statistic,
#' pseudo R-square statistic, and analytical p-values for all predictors
#' jointly and for each predictor individually, conditioned on the rest.
#'
#' This function is the fastest approach to conducting MDMR. It uses the
#' fastest known computational strategy to compute the MDMR test statistic (see
#' Appendix A of McArtor & Lubke, 2015), and it uses fast, analytical p-values.
#'
#' The slowest part of conducting MDMR is now the necessary eigendecomposition
#' of the \code{G} matrix, whose computation time is a function of
#' \eqn{n^3}. If MDMR is to be conducted multiple times on the same
#' distance matrix, it is recommended to compute eigenvalues of \code{G} in
#' advance and pass them to the function rather than computing them every
#' time \code{mdmr} is called, as is the case if the argument \code{lambda}
#' is left \code{NULL}.
#'
#' The distance matrix \code{D} can be passed to \code{mdmr} as either a
#' distance object or a symmetric matrix.
#'
#' @param X A \eqn{n x p} matrix or data frame of predictors. Factors are
#' allowed to be included, and will be dummy-coded prior to analysis.
#' @param D Distance matrix computed on the outcome data. Can be either a
#' matrix or an R \code{\link{dist}} object. Either \code{D} or \code{G}
#' must be passed to \code{mdmr()}.
#' @param G Gower's centered similarity matrix computed from \code{D}.
#' Either \code{D} or \code{G} must be passed to \code{mdmr}.
#' @param lambda Optional argument: Eigenvalues of \code{G}.
#' Eigendecomposition of large \code{G} matrices can be somewhat time
#' consuming, and the theoretical p-values require the eigenvalues of
#' \code{G}. If MDMR is to be conducted multiple times on one distance
#' matrix, it is advised to conduct the eigendecomposition once and pass the
#' eigenvalues to \code{mdmr()} directly each time.
#' @param return.lambda Logical; indicates whether or not the eigenvalues of
#' \code{G} should be returned, if calculated. Default is \code{FALSE}.
#' @param start.acc Starting accuracy of the Davies (1980) algorithm
#' implemented in the \code{\link{davies}} function in the \code{CompQuadForm}
#' package (Duchesne &  De Micheaux, 2010) that \code{mdmr()} uses to compute
#' MDMR p-values.
#' @param ncores Integer; if \code{ncores} > 1, the \code{\link{parallel}}
#' package is used to speed computation. Note: Windows users must set
#' \code{ncores = 1} because the \code{parallel} pacakge relies on forking. See
#' \code{mc.cores} in the \code{\link{mclapply}} function in the
#' \code{parallel} pacakge for more details.
#'
#' @return An object with five elements and a summary function. Calling
#' \code{summary(mdmr.res)} produces a data frame comprised of:
#' \item{Statistic}{Value of the corresponding MDMR test statistic}
#' \item{Pseudo R2}{Size of the corresponding effect on the
#' distance matrix}
#' \item{p-value}{Analytical p-value}
#' In addition to the information in the three columns comprising
#' \code{summary(res)}, the \code{res} object also contains  a vector reporting
#' the precision of each p-value as reported by the \code{davies} function in
#' \code{CompQuadForm} and a vector of the eigenvalues of \code{G} (if
#' \code{return.lambda = T}).
#'
#' @author \packageAuthor{MDMR}
#'
#' @references Davies, R. B. (1980). The Distribution of a Linear Combination of
#'  chi-square Random Variables. Journal of the Royal Statistical Society.
#'  Series C (Applied Statistics), 29(3), 323-333.
#'
#'  Duchesne, P., & De Micheaux, P.L. (2010). Computing the distribution of
#'  quadratic forms: Further comparisons between the Liu-Tang-Zhang
#'  approximation and exact methods. Computational Statistics and Data
#'  Analysis, 54(4), 858-862.
#'
#'  McArtor, D.B. & Lubke, G.H. (submitted). Extending multivariate distance
#'  matrix regression with an effect size measure and the distribution of the
#'  test statistic.
#'
#' @examples
#'# --- The following two approaches yield equivalent results --- #
#'# Approach 1
#'data(mdmrdata)
#'D <- dist(Y.mdmr, method = 'euclidean')
#'mdmr(X = X.mdmr, D = D)
#'
#'# Approach 2
#'data(mdmrdata)
#'D <- dist(Y.mdmr, method = 'euclidean')
#'G <- gower(D)
#'mdmr(X = X.mdmr, G = G)
#'
#' @importFrom CompQuadForm davies
#' @importFrom parallel mclapply
#' @export
mdmr <- function(X, D = NULL, G = NULL, lambda = NULL, return.lambda = F,
                 start.acc = 1e-20, ncores = 1){
  # Make sure "D" is not interpreted as the D function
  if(is.function(D)){
    stop(paste0('Please provide either a distance matrix or a ',
                'Gower-centered distance matrix.'))
  }

  # Make sure either D or G was provided
  if(is.null(D) & is.null(G)){
    stop(paste0('Please provide either a distance matrix ',
                'or a Gower-centered distance matrix.'))
  }

  # Make sure D is a distance matrix
  if(!is.null(D)){
    if(nrow(as.matrix(D)) != ncol(as.matrix(D))){
      stop(paste0('Please provide either a distance matrix ',
                  'or a Gower-centered distance matrix.'))
    }
  }

  # If G was not provided, compute it from D
  if(is.null(G)){
    G <- gower(D)
  }

  # Handle potential factors, no-named variables
  X <- stats::model.matrix(~ . , data = as.data.frame(X))
  xnames <- colnames(X)

  # Record the number of items and sample size
  p <- ncol(X)-1
  n <- nrow(X)


  # ===================== FUNCTION TO COMPUTE P-VALUES ======================= #
  pmdmr <- function(q, lambda, k, p, n = length(lambda),
                    lim = 50000, acc = start.acc){
    # Use the eigenvalues of G and the test statistic q to make a vector
    # corresponding to all of the weights for the composite chi-square variables
    # See Equation 12 in McArtor & Lubke (2015)
    gamma <- c(lambda,  -q * lambda)

    # Aggregate the degrees of freedom for each composite chi-square variable
    # See Equation 12 in McArtor & Lubke (2015)
    nu <- c(rep(k, length(lambda)), rep(n-p-1, length(lambda)))

    # Call the Davies function at zero using the given weights and df, along
    # with the starting values of the metaparameters of the algorithm
    pv <- CompQuadForm::davies(0, lambda = gamma, h = nu, lim = lim, acc = acc)

    # Check error status. If there was an error, return entire davies object
    if(pv$ifault != 0){
      #     warning('Error in davies() procedure, please check results.')
      return(pv)
    }

    # If the p-value is below zero, interpret as an error and break
    if(pv$Qq < 0){
      #     warning('Error in davies() procedure, please check results.')
      return(pv)
    }

    # If there was no error, return only the p-value
    if(pv$ifault == 0){
      return(pv$Qq)
    }
  }


  # ======================= Omnibus Test Statistic =========================== #
  # Compute hat matrix
  H <- tcrossprod(tcrossprod(X, solve(crossprod(X))), X)

  # Computational trick: H is idempotent, so H = HH. tr(ABC) = tr(CAB), so
  # tr(HGH) = tr(HHG) = tr(HG). Also, tr(AB) = vec(A)'vec(B), so
  vh <- matrix(H, nrow = 1)
  vg <- matrix(G, ncol = 1)

  # -- Do the computation in chunks to not overload memory -- #
  # (This approach is also faster than doing the whole multiplication at once)
  len <- ncol(vh)
  # minimum chunk size is 3000^2 or n^2. Even my personal laptop has no problems
  # with this size
  chunksize <- min(len, 3000^2)
  # Set number of chunks, and make sure the last chunk isn't too big or small
  nchunk <- ceiling(len/chunksize)
  chunkEnd <- seq(chunksize, nchunk*chunksize, chunksize)
  chunkEnd[nchunk] <- len

  # Start the computation then iterate over all chunks
  numer <- vh[,1:chunksize] %*% vg[1:chunksize,]
  if(nchunk>1){
    for(i in 2:nchunk){
      ind <- (chunkEnd[i-1]+1):chunkEnd[i]
      numer <- numer + vh[,ind] %*% vg[ind,]
    }
  }

  # Numerical trick: tr((I-H)G(I-H)) = tr(G) - tr(HGH), so
  trG <- sum(diag(G))
  denom <- trG - numer

  # Save omnibus F
  pr2.omni <- as.numeric(numer / trG)
  f.omni <- as.numeric(numer / denom)

  # Omnibus p-value
  # Compute the eigenvalues of G if they were not provided
  if(is.null(lambda)){
    lambda <- eigen(G, only.values = T)$values
  }

  # Initiailze accuracy of the Davies algorithm
  acc.omni <- start.acc
  pv.omni <- pmdmr(f.omni, lambda, p, p, n, acc = acc.omni)

  # If the davies procedure threw an error, decrease the accuracy
  while(length(pv.omni) > 1){
    acc.omni <- acc.omni * 10
    pv.omni <- pmdmr(q = f.omni, lambda = lambda, k = p, p = p, n = n,
                     acc = acc.omni)
  }

  # ====================== Tests of Each Predictor =========================== #

  # --- CASE 1: NO PARALLELIZATION --- #
  if(ncores == 1){

    # Get vectorized hat matrices for each conditional effect
    Hs <- lapply(2:(p+1), function(k){
      Xs <- X[,-k]
      Hs <- tcrossprod(tcrossprod(Xs, solve(crossprod(Xs))), Xs)
      matrix(H - Hs, nrow = 1)
    })

    # Compute SSD due to conditional effect
    numer.x <- unlist(lapply(Hs, function(vhs){
      numer <- vhs[,1:chunksize] %*% vg[1:chunksize,]
      if(nchunk>1){
        for(i in 2:nchunk){
          ind <- (chunkEnd[i-1]+1):chunkEnd[i]
          numer <- numer + vhs[,ind] %*% vg[ind,]
        }
      }
      numer
    }))

    # Rescale to get either test statistic or pseudo r-square
    f.x <- numer.x / denom
    pr2.x <- numer.x / trG

    # Get p-values
    p.res <- lapply(1:p, function(k){
      item.acc <- start.acc
      pv <- pmdmr(q = f.x[k], lambda = lambda, k = 1,
                  p = p, n = n, acc = item.acc)

      # If the davies procedure threw an error, decrease the accuracy
      while(length(pv) > 1){
        item.acc <- item.acc * 10
        pv <- pmdmr(q = f.x[k], lambda = lambda, k = 1,
                    p = p, n = n, acc = item.acc)
      }
      c(pv.x = pv, acc = item.acc)
    })
    pv.x <- unlist(lapply(p.res, function(p){p[[1]]}))
    acc.x <- unlist(lapply(p.res, function(p){p[[2]]}))

  }

  # --- CASE 2: WITH PARALLELIZATION --- #
  if(ncores > 1){

    # Get vectorized hat matrices for each conditional effect
    Hs <- parallel::mclapply(2:(p+1), function(k){
      Xs <- X[,-k]
      Hs <- tcrossprod(tcrossprod(Xs, solve(crossprod(Xs))), Xs)
      matrix(H - Hs, nrow = 1)
    },
    mc.preschedule = TRUE, mc.set.seed = TRUE,
    mc.silent = FALSE, mc.cores = ncores,
    mc.cleanup = TRUE, mc.allow.recursive = TRUE)

    # Compute SSD due to conditional effect
    numer.x <- unlist(parallel::mclapply(Hs, function(vhs){
      numer <- vhs[,1:chunksize] %*% vg[1:chunksize,]
      if(nchunk>1){
        for(i in 2:nchunk){
          ind <- (chunkEnd[i-1]+1):chunkEnd[i]
          numer <- numer + vhs[,ind] %*% vg[ind,]
        }
      }
      numer
    },
    mc.preschedule = TRUE, mc.set.seed = TRUE,
    mc.silent = FALSE, mc.cores = ncores,
    mc.cleanup = TRUE, mc.allow.recursive = TRUE))

    # Rescale to get either test statistic or pseudo r-square
    f.x <- numer.x / denom
    pr2.x <- numer.x / trG

    # Get p-values
    p.res <- parallel::mclapply(1:p, function(k){
      item.acc <- start.acc
      pv <- pmdmr(q = f.x[k], lambda = lambda, k = 1,
                  p = p, n = n, acc = item.acc)

      # If the davies procedure threw an error, decrease the accuracy
      while(length(pv) > 1){
        item.acc <- item.acc * 10
        pv <- pmdmr(q = f.x[k], lambda = lambda, k = 1,
                    p = p, n = n, acc = item.acc)
      }
      c(pv.x = pv, acc = item.acc)
    },
    mc.preschedule = TRUE, mc.set.seed = TRUE,
    mc.silent = FALSE, mc.cores = ncores,
    mc.cleanup = TRUE, mc.allow.recursive = TRUE)

    pv.x <- unlist(lapply(p.res, function(p){p[[1]]}))
    acc.x <- unlist(lapply(p.res, function(p){p[[2]]}))
  }




  # ====================== Combine Output =========================== #
  stat <- c(f.omni, f.x)
  pr2 <- c(pr2.omni, pr2.x)
  pv <- c(pv.omni, pv.x)
  pv.acc <- c(acc.omni, acc.x)

  names(stat) <- names(pr2) <- names(pv) <- names(pv.acc) <-
    c('Omnibus', xnames[-1])

  if(return.lambda){
    out <- list('stat' = stat, 'pr.sq' = pr2, 'pv' = pv, 'p.prec' = pv.acc,
                lambda = lambda)
  }
  if(!return.lambda){
    out <- list('stat' = stat, 'pr.sq' = pr2, 'pv' = pv, 'p.prec' = pv.acc,
                lambda = NULL)
  }

  class(out) <- c('mdmr', class(out))

  return(out)
}


#' Print MDMR Object
#'
#' \code{print} method for class \code{mdmr}
#'
#' @param x Output from \code{mdmr}
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' \item{p-value}{Analytical p-values for the omnibus test and each predictor}
#'
#' @author \packageAuthor{MDMR}
#'
#'
#' @export
print.mdmr <- function(x, ...){
  cat('p-values for each predictor and the omnibus test:', fill = T)
  hold <- data.frame(format.pval(x$pv), row.names = names(x$pv))
  names(hold) <- ''
  print(hold)
}

#' Summarizing MDMR Results
#'
#' \code{summary} method for class \code{mdmr}
#'
#' @param object Output from \code{mdmr}
#' @param ... Further arguments passed to or from other methods.
#'
#' @return
#' \item{Statistic}{Value of the corresponding MDMR test statistic}
#' \item{Pseudo R2}{Effect size of the corresponding effect on the
#' distance matrix}
#' \item{p-value}{Analytical p-value}
#'
#' @author \packageAuthor{MDMR}
#'
#' @references Davies, R. B. (1980). The Distribution of a Linear Combination of
#'  chi-square Random Variables. Journal of the Royal Statistical Society.
#'  Series C (Applied Statistics), 29(3), 323-333.
#'
#'  Duchesne, P., & De Micheaux, P.L. (2010). Computing the distribution of
#'  quadratic forms: Further comparisons between the Liu-Tang-Zhang
#'  approximation and exact methods. Computational Statistics and Data
#'  Analysis, 54(4), 858-862.
#'
#'  McArtor, D.B. & Lubke, G.H. (submitted). Extending multivariate distance
#'  matrix regression with an effect size measure and the distribution of the
#'  test statistic.
#'
#' @examples
#'# --- The following two approaches yield equivalent results --- #
#'# Approach 1
#'data(mdmrdata)
#'D <- dist(Y.mdmr, method = 'euclidean')
#'mdmr.res <- mdmr(X = X.mdmr, D = D)
#'summary(mdmr.res)
#'
#' @export
summary.mdmr <- function(object, ...){
  print(data.frame('Statistic' =
                     format(object$stat, digits = 3),
                   'Pseudo R2' = format.pval(object$pr.sq, digits = 3),
                   'p-value' = format.pval(object$pv)))
  res <- data.frame('Statistic' = object$stat,
                    'Pseudo R2' = object$pr.sq,
                    'p-value' = object$pv)
  invisible(res)
}






#' Compute univariate MDMR effect sizes
#'
#' \code{delta} computes permutation-based effect sizes on individual items
#' comprising the distance matrix outcome used in multivariate distance matrix
#' regression. It returns the omnibus estimates of delta (i.e. effect size of
#' the entire design matrix on each outcome) as well as estimates of each
#' pair-wise effect size (i.e. the effect of each predictor on each outcome
#' variable, conditional on the rest of the predictors).
#'
#' See McArtor & Lubke (submitted) for a detailed description of how delta is
#' computed. Note that it is a relative measure of effect, quantifying which
#' effects are strong (high values of delta) and weak (low values of delta)
#' within a single analysis, but estimates of delta cannot be directly compared
#' across different datasets.
#'
#' There are two options for using this function. The first option is to
#' specify the predictor matrix \code{X}, the outcome matrix \code{Y}, the
#' distance type \code{dtype} (supported by "dist" in R), and number of
#' iterations \code{niter}. This option conducts the permutation of each Y-item
#' \code{niter} times (to average out random association in each permutation)
#' and reports the median estimates of delta over the \code{niter} reps.
#'
#' The second option is to specify \code{X}, \code{G}, and \code{G.list}, a
#' list of G matrices where the permutation has already been done for each item
#' comprising Y. The names of the elements in \code{G.list} should correspond
#' to the names of the variables that were permuted. This option is implemented
#' so that delta can be computed when MDMR is being used in conjunction with
#' distance metrics not supported by \code{dist}.
#'
#' @param X A \eqn{n x p} matrix or data frame of predictors. Factors are
#' allowed to be included, and will be dummy-coded prior to analysis.
#' @param Y Outcome data: \eqn{n x q} matrix of scores along the
#' dependent variables.
#' @param dtype Measure of dissimilarity that will be used by \code{\link{dist}}
#' to compute the distance matrix based on \code{Y}. As is the case when calling
#' \code{dist} directly, this must be one of \code{"euclidean"},
#' \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"} or
#' \code{"minkowski"}, and unambiguous substring can be given.
#' @param niter Number of times to permute each outcome item in the procedure
#' to compute delta. The final result is the average of all \code{niter}
#' iterations. Higher values of \code{niter} require more computation time, but
#' result in more precise estimates.
#' @param G Gower's centered similarity matrix computed from \code{D}.
#' Either \code{D} or \code{G} must be passed to \code{mdmr()}.
#' @param G.list List of length \eqn{q} where the \eqn{i^{th}}
#' element contains the \code{G} matrix computed from distance a matrix that
#' was computed on a version of \code{Y} where the \eqn{i^{th}}
#' column has been randomly permuted.
#' @param plot.res Logical; Indicates whether or not a heat-map of the results
#' should be plotted.
#' @param ncores Integer; if \code{ncores} > 1, the \code{\link{parallel}}
#' package is used to speed computation. Note: Windows users must set
#' \code{ncores = 1} because the \code{parallel} pacakge relies on forking. See
#' \code{mc.cores} in the \code{\link{mclapply}} function in the
#' \code{parallel} pacakge for more details.
#' @param seed Integer; sets seed for the permutations of each variable
#' comprising Y so that results can be replicated.
#'
#' @return A data frame whose rows correspond to the omnibus effects and the
#' effect of each individual predictor (conditional on the rest), and whose
#' columns correspond to each outcome variable whose effect sizes are being
#' quantified. If \code{plot.res = TRUE}, a heat-map is plotted of this data
#' frame to easily identify the strongest effects. Note that the heatmap is
#' partitioned into the omnibus effect (first row) and pair-wise effects
#' (remaining rows), because otherwise the omnibus effect would dominate the
#' heatmap.
#'
#' @author \packageAuthor{MDMR}
#'
#' @references  McArtor, D.B. & Lubke, G.H. (submitted). Extending
#' multivariate distance matrix regression with an effect size measure and the
#' distribution of the test statistic.
#'
#' @examples
#' data(mdmrdata)
#' # --- Method 1 --- #
#' delta(X.mdmr, Y = Y.mdmr, dtype = 'euclidean', niter = 1, seed = 12345)
#'
#' # --- Method 2 --- #
#' D <- dist(Y.mdmr, method = 'euclidean')
#' G <- gower(D)
#' q <- ncol(Y.mdmr)
#' G.list <- vector(mode = 'list', length = q)
#' names(G.list) <- names(Y.mdmr)
#' for(i in 1:q){
#'    Y.shuf <- Y.mdmr
#'    Y.shuf[,i] <- sample(Y.shuf[,i])
#'    G.list[[i]] <- gower(dist(Y.shuf, method = 'euclidean'))
#' }
#' delta(X.mdmr, G = G, G.list = G.list)
#'
#' @export
#' @importFrom parallel mclapply
delta <- function(X, Y = NULL, dtype = NULL, niter = 10,
                  G = NULL, G.list = NULL, plot.res = F,
                  ncores = 1, seed = NULL){
  # ============================================================================
  # Step 1: Check input type
  # ============================================================================
  # If a list of distance matrices is not provided....
  if(is.null(G.list)){
    # Make sure raw data and distance METRIC is provided
    if(any(is.null(Y), is.null(dtype))){
      stop(paste0('Input either Y and the distance type ',
                  'or a list of distance matrices'))
    }
    # Set the method to raw (see below)
    method <- 'raw'
  }
  # If a list of distance matrices is provided...
  if(!is.null(G.list)){
    # Set the method to list (see below)
    method <- 'list'
  }

  # Handle potential factors, no-named variables
  X <- stats::model.matrix(~ . , data = as.data.frame(X))[,-1]
  xnames <- colnames(X)
  p <- ncol(X)
  n <- nrow(X)


  # ============================================================================
  # Step 2: Function to compute pseudo R-square (it's done a lot here with the
  # same X and H: only G changes)
  # ============================================================================

  # --- Do all the computations that are required every time pseudo.r2 is called
  # Overall hat matrix
  X.hat <- cbind(1, X)
  H <- tcrossprod(tcrossprod(X.hat, solve(crossprod(X.hat))), X.hat)
  vh <- matrix(H, nrow = 1)

  # Do the computation in chunks to not overload memory
  len <- ncol(vh)
  chunksize <- min(len, 3000^2)
  nchunk <- ceiling(len/chunksize)
  chunkEnd <- seq(chunksize, nchunk*chunksize, chunksize)
  chunkEnd[nchunk] <- len

  # Hat matrices for each conditional effect
  Xs <- lapply(2:(p+1), function(i){
    X.hat[,-i]
  })
  Hs <- lapply(Xs, function(x){
    H - tcrossprod(tcrossprod(x, solve(crossprod(x))), x)
  })
  vhs <- lapply(Hs, function(h){
    matrix(h, nrow = 1)
  })

  # ----- FUNCTION TO COMPUTE PSEUDO-R-SQUARE WITHIN THIS RUN OF DELTA ----- #

  pseudo.r2 <- function(G){

    # =========================== Omnibus Test =============================== #
    # Computational trick: H is idempotent, so H = HH. tr(ABC) = tr(CAB), so
    # tr(HGH) = tr(HHG) = tr(HG). Also, tr(AB) = vec(A)'vec(B), so
    vg <- matrix(G, ncol = 1)
    # Start the computation then iterate over all chunks
    numer <- vh[,1:chunksize] %*% vg[1:chunksize,]
    if(nchunk>1){
      for(i in 2:nchunk){
        ind <- (chunkEnd[i-1]+1):chunkEnd[i]
        numer <- numer + vh[,ind] %*% vg[ind,]
      }
    }

    # pseudo-R2 is defined as tr(HGH)/tr(G), so
    denom <- sum(diag(G))

    # Omnibus pseudo R-square
    res <- numer/denom

    # ===================== Tests of Each Predictor ========================== #
    numer.x <- unlist(lapply(1:p, function(k){
      numer <- vhs[[k]][,1:chunksize] %*% vg[1:chunksize,]
      if(nchunk>1){
        for(i in 2:nchunk){
          ind <- (chunkEnd[i-1]+1):chunkEnd[i]
          numer <- numer + vhs[[k]][,ind] %*% vg[ind,]
        }
      }
      numer
    }))

    r2.x <- numer.x / denom

    # ====================== Return Results =========================== #
    res <- c(res, r2.x)
    names(res) <- c('Omnibus', xnames[-1])
    return(res)
  }

  # ============================================================================
  # Step 3: Computation if raw data are provided
  # ============================================================================
  if(method == 'raw'){
    # ----- Manage Input ----- #
    Y <- as.matrix(Y)
    q <- ncol(Y)
    ynames <- colnames(data.frame(Y))
    if(all(ynames == paste0('X', 1:q))){
      ynames <- paste0('Y', 1:q)
    }

    G <- gower(stats::dist(Y, method = dtype))

    n <- nrow(X)

    # ----- Populate delta matrices using jackknife procedure ----- #
    # Get the "real" pseudo R-square
    pr2 <- pseudo.r2(G)

    # IF THE USER IS NOT USING PARALLELIZATION
    if(ncores == 1){
      # If a seed is provided by the function, use it
      if(!is.null(seed)){
        set.seed(seed)
      }
      # Compute pseudo R-square with each item jackknifed "niter" times
      jack.pr2 <-
        lapply(1:niter,
               function(i){
                 jackknifed.y <- lapply(1:q, function(k){
                   y.jack <- Y
                   y.jack[,k] <- sample(y.jack[,k], size = n, replace = F)
                   y.jack
                 })
                 res <- lapply(jackknifed.y, function(yy){
                   GG <- gower(stats::dist(yy, method = dtype))
                   pseudo.r2(GG)
                 })
                 res <- matrix(unlist(res), nrow = p+1, ncol = q)
                 dimnames(res) <- list(c('Omnibus', xnames),
                                       paste0(ynames, '.jack'))
                 res
               })
      # Subtract the jackknifed pseudo R-squares from the real pseudo R-squares
      # to get the delta statistics for each rep
      jack.pr2 <- lapply(jack.pr2, function(jack){
        apply(jack, 2, function(x){pr2 - x})
      })
    }

    # IF THE USER IS USING PARALLELIZATION
    if(ncores > 1){

      # If no seed is specified by the user, generate a random one - mclapply
      # requires one, which is why "seed" is an argument rather than using the
      # standard approach of just setting a seed prior to running the function
      if(is.null(seed)){
        seed <- round(stats::runif(1,0,1) * 1e5)
      }

      # Compute pseudo R-square with each item jackknifed "niter" times
      jack.pr2 <-
        parallel::mclapply(1:niter,
                           function(i){
                             set.seed(seed+i)
                             jackknifed.y <- lapply(1:q, function(k){
                               y.jack <- Y
                               y.jack[,k] <- sample(y.jack[,k], size = n,
                                                    replace = F)
                               y.jack
                             })
                             res <- lapply(jackknifed.y, function(yy){
                               GG <- gower(stats::dist(yy, method = dtype))
                               pseudo.r2(GG)
                             })
                             res <- matrix(unlist(res), nrow = p+1, ncol = q)
                             dimnames(res) <- list(c('Omnibus', xnames),
                                                   paste0(ynames, '.jack'))
                             res
                           },
                           mc.preschedule = TRUE, mc.set.seed = TRUE,
                           mc.silent = FALSE, mc.cores = ncores,
                           mc.cleanup = TRUE, mc.allow.recursive = TRUE)

      # Subtract the jackknifed pseudo R-squares from the real pseudo R-squares
      # to get the delta statistics for each rep
      jack.pr2 <- parallel::mclapply(jack.pr2, function(jack){
        apply(jack, 2, function(x){pr2 - x})
      }, mc.preschedule = TRUE, mc.set.seed = TRUE,
      mc.silent = FALSE, mc.cores = ncores,
      mc.cleanup = TRUE, mc.allow.recursive = TRUE)
    }


    # --- Compute Median of Delta --- #
    delta.med <- matrix(0, nrow = p + 1, ncol = q)
    dimnames(delta.med) <-  list(c('Omnibus', xnames),
                                 paste0(ynames, '.jack'))
    for(i in 1:nrow(delta.med)){
      for(j in 1:ncol(delta.med)){
        delta.med[i,j] <- stats::median(
          unlist(lapply(jack.pr2, function(x){x[i,j]})))
      }
    }

  }

  # ============================================================================
  # Step 4: Computation if list of distance matrices is provided
  # ============================================================================
  if(method == 'list'){

    # ----- Manage Input ----- #
    X <- as.matrix(X)
    xnames <- colnames(data.frame(X))
    p <- ncol(X)
    q <- length(G.list)

    ynames <- names(G.list)
    if(is.null(ynames)){
      ynames <- paste0('Y', 1:q)
    }


    n <- nrow(X)

    # ----- Populate delta matrices using jackknife procedure ----- #
    # Get the "real" pseudo R-square
    pr2 <- pseudo.r2(G)
    # Get each permuted pseudo R-square
    # IF THE USER IS NOT USING PARALLELIZATION
    if(ncores == 1){
      # Compute pseudo R-square with each item jackknifed "niter" times
      jack.pr2 <-
        matrix(unlist(lapply(G.list, function(GG){pseudo.r2(GG)})),
               nrow = p+1, ncol = q)
      # Subtract the jackknifed pseudo R-squares from the real pseudo R-squares
      # to get the delta statistics for each rep
      jack.pr2 <- apply(jack.pr2, 2, function(x){pr2 - x})
      dimnames(jack.pr2) <- list(c('Omnibus', xnames),
                                 paste0(ynames, '.jack'))
    }

    # IF THE USER IS USING PARALLELIZATION
    if(ncores > 1){

      # If no seed is specified by the user, generate a random one - mclapply
      # requires one, which is why "seed" is an argument rather than using the
      # standard approach of just setting a seed prior to running the function
      if(is.null(seed)){
        seed <- round(stats::runif(1,0,1) * 1e5)
      }

      # Compute pseudo R-square with each item jackknifed "niter" times
      jack.pr2 <-
        matrix(unlist(
          parallel::mclapply(G.list, function(GG){pseudo.r2(GG)},
                             mc.preschedule = TRUE, mc.set.seed = TRUE,
                             mc.silent = FALSE, mc.cores = ncores,
                             mc.cleanup = TRUE, mc.allow.recursive = TRUE)),
          nrow = p+1, ncol = q)

      # Subtract the jackknifed pseudo R-squares from the real pseudo R-squares
      # to get the delta statistics for each rep
      jack.pr2 <- apply(jack.pr2, 2, function(x){pr2 - x})
      dimnames(jack.pr2) <- list(c('Omnibus', xnames),
                                 paste0(ynames, '.jack'))
    }

    # --- In this case, there's only one rep, so the median is the single rep
    delta.med <- jack.pr2
  }


  # ============================================================================
  # Step 5: Plot
  # ============================================================================
  if(plot.res){
    graphics::plot(NA, xlim = c(0.5, q+0.5), ylim = c(0.5,p+0.5+1), xaxt = 'n',
         yaxt = 'n',  xlab = '', ylab = '', bty = 'n',
         main = 'MDMR Effect Sizes')
    graphics::axis(1, at = 1:q, labels = c(ynames), las = 2)
    graphics::axis(2, at = (p+1):1, labels = c('Omnibus', xnames), las = 1)

    # Convert to z scores for shading
    z.scores <- matrix(scale(c(delta.med)), nrow = p+1, ncol = q)
    z.scores[delta.med < 0] <- -9999
    omni.cols <- stats::pnorm(z.scores[1,])
    pairwise.cols <- stats::pnorm(z.scores[-1,])

    for(i in 1:nrow(delta.med)){
      for(j in 1:ncol(delta.med)){
        x.low <- j-0.5
        y.low <- p-i+1-0.5+1
        x.up <- j+0.5
        y.up <- p-i+1+0.5+1

        if(i == 1){
          # Y importances
          graphics::rect(x.low, y.low, x.up, y.up,
               col = grDevices::rgb(0, 0, 1, omni.cols[j]))
        }
        if(i > 1){
          # XY Importances
          graphics::rect(x.low, y.low, x.up, y.up,
               col = grDevices::rgb(0,0.75,0, pairwise.cols[i-1,j]))
        }

        # Effect Size text
        graphics::text(x = j, y = p - i + 1 + 1, col = 'white',
             labels = formatC(delta.med[i,j], format = 'g', digits = 2),
             cex = 0.75)
      }
    }
  }


  # ============================================================================
  # Step 6: Output
  # ============================================================================
  return(delta.med)
}




#' Simulated predictor data to illustrate the MDMR package.
#'
#' See package vignette by calling \code{vignette('mdmr-vignette')}.
"X.mdmr"



#' Simulated outcome data to illustrate the MDMR package.
#'
#' See package vignette by calling \code{vignette('mdmr-vignette')}.
"Y.mdmr"
