# This function computes the correlation between any of the N rows of X and any 
# of the M rows of Y, according to the correlation function specified in input 'string'.
# Output is an NxM matrix.
# Each row of X/Y is a k-dimensional vector. Correlation lengths to use for 
# computation of the correlation function are specified in d. 
#
#
# INPUTS:
# X: matrix or data-frame, dimension (N,k)
# Y: matrix or data frame, dimension (M,k)
# d: 1-dim or k-dim vector of positive correlation lengths. 
#    If 1-dim, it is extended to a vector (d,d,...,d) of length k.
# string: one of the following, to specify which correlation function to use:
#         'exp2'     (squared exponential)
#         'matern32' (matern 3/2, regularity 3/2)
#         'matern52' (matern 5/2, regularity 5/2)
#         'abs_exp'  (absolute exponential, not C^1).
#
# OUTPUT:
# C: matrix of dimension (N,M): C[i,j] is the correlation between X[i,] and Y[j,]

# COMPUTATIONAL NOTE: 
# The function uses 3D tensor multiplication to substantially increase speed wrt nested 
# for loops in i&j: approx 1 order of magnitude of gain in speed is 
# observed for any two orders of magnitude increase in NxM. 

Corr.fun <- function(X, Y, d, string){

### DEFINE VARIABLES AND CHECK CORRECT INPUT DIMENSIONS ###
  
  #Convert dataframes into matrices if needed  
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    d <- as.numeric(d)
  
  # Store sizes
    N <- dim(X)[1]
    M <- dim(Y)[1]
    k <- dim(X)[2]
    
  # Stop if dimensions of inputs don't match
    if (dim(X)[2] != dim(Y)[2])
      stop('Inputs X and Y must have the same number of columns')
  
  # Stop if wrong number of correlation lengths is provided (or extend if just one is provided)  
    if (length(d)==1)
      d <- rep(d,k)
    else if (length(d) != k)
      stop('Length of input d must be equal to the number of columns of X and Y, or 1')
    
  
### CREATE 3D ARRAY Diff, WITH Diff[i,j,] = X[i,]-Y[j,] for i=1:N, j=1:M ###
  
  # Stack up copies of X and Y into 3D arrays
    X3 <- array(X, dim=c(dim(X),M) )  # dim = (N,k,M)
    Y3 <- array(Y, dim=c(dim(Y),N))   # (M,k,N)
  
  # Permute dimensions so that:
  # X3[i,j,:] is the k-dim vector X[i,:]
  # Y3[i,j,:] is the k-dim vector Y[j,:].
    X3 <- aperm(X3, c(1,3,2))         # (N,M,k)
    Y3 <- aperm(Y3, c(3,1,2))         # (N,M,k)  
  
  # Define Diff=X3-Y3, Diff[i,j,:] = X[i,]-Y[j,:]
    Diff <- X3-Y3                     # (N,M,k)
  # and rescale the k dimensions by correlation lengths
    for (i in 1:k)
      Diff[,,i] <- Diff[,,i]/d[i] 
  
### RETURN NxM MATRIX C OF CORRELATIONS ###
    
# Computations differ according to type of correlation function specified in 'string'. 
    A <- rowSums(Diff^2, dims=2) # A[i,j] = sum of the k elements Diff[i,j,:]^2 
    
    if (identical(string,'exp2'))
      C <- exp(-A)
    else if (identical(string,'abs_exp'))
      C <- exp(-sqrt(A))
    else if (identical(string,'matern32')){
      A1 <- sqrt(3*A)
      C <- (1+ A1)*exp(-A1)
    }
    else if (identical(string,'matern52')){
      A1 <- sqrt(5*A)
      A2 <- 5*A/3
      C <- (1+ A1 + A2)*exp(-A1)
    }
    else
      stop('Input \'string\' must be one of: \'exp2\', \'matern32\', \'matern52\', \'abs_exp\'.')
      
    return(C)
}





Cross_Val <- function(Design.points.full, Design.regr.full, y.full, beta, Cov.beta, d, nu, sigma2, string){
  n <- dim(Design.points.full)[1]

  # PRIOR MEAN OF DESIGN POINTS
  prior.mean.Design.full <- Design.regr.full %*% beta              # n-dim vector

  # PRIOR COVARIANCE OF DESIGN POINTS
  a <- Design.regr.full %*% Cov.beta %*% t(Design.regr.full)       # nxn, prior covariance of regression part
  b <- sigma2*corr.fun(Design.points.full, Design.points.full, d, string)  # nxn, prior GP cov
  c <- nu*diag(n)                                                  # variance of nugget term
  A.full <- a + b + c  
  
  # PRIOR ERROR (REAL VALUE MINUS PREDICTION)
  e.full <- y - prior.mean.Design.full
  
  # VARIABLES TO STORE POSTERIOR CROSS-VALIDATED MEAN AND VARIANCE
  M <- array(dim=c(n,1))
  V <- array(dim=c(n,1))
  
  ind.full <- 1:n
  
  for (i in 1:n){
    
    ind <- ind.full[-i]
    
    # Quantities needed in emulation computation
    y <- y.full[ind]
    A <- A.full[ind, ind]
    e <- e.full[ind]
    tx <- A.full[ind,i, drop=F]
    
    # POSTERIOR MEAN AND COVARIANCE FOR INPUT PARAMETERS
    K <- solve(A, tx)                                  
    post.mean <- prior.mean.Design.full[i] + t(e)%*%K  
    post.var <-   A.full[i,i]   - t(tx)%*%K
    
    M[i] <- post.mean
    V[i] <- post.var
    
    if (i%%10 == 0){
      cat(sprintf("Iteration number %d out of %d.\n", i, n))
    }
    
  }
  
  return( list("Mean"=M, "Variance"=V) )
}



###########################################################################################
# The following function creates a list with element names as in "element.names".
# By default (ie, no variable "indices" provided) each element of the list is an Nx2 matrix 
# with all entries equal to "val", and column names "Mean", "Var".
# If "indices" is provided as a vector of integers, only the list element corresponding
# to those indices will be populated by the matrix as above, the remaining ones will be empty.

Create.my.list <- function(N, element.names, indices=c(1:length(element.names)), val=0){
  
  my.list <- sapply(element.names, function(x) NULL)
  
  mat <- matrix(data = val, nrow = N, ncol = 2)
  colnames(mat) <- c("Mean", "Var")
  
  for (i in indices){
    my.list[[i]] <- mat
  }
  return(my.list)
  
}



