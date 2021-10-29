################################################################################
### BRIEF DESCRIPTION
#
# The function BL.Emul trains a Bayes Linear emulator of a computer simulator f, 
# based on observed outputs y=(y1,...,yn) at n design points. 
# Emulation (ie, output prediction) is performed at N test points.
#
# Specification of inputs and outputs below, with full details available at: 
# https://github.com/dario-domi/Bayes-Linear-Emulation/blob/main/Documentation.pdf


# OVERVIEW OF EMULATOR'S MODEL:
#
# The following statistical model of f is assumed to build the emulator:
#
# f(x) = regression_term(x_A) + eta(x_A) + nugget_term(x),
#
# where x_A are the "active inputs", explaining most of the variability in f. 
# The three terms are as follows:
#
# 1) The regression term is a linear combination of basis functions g_j(x_A) (regressors):
#    regression_term(x_A) = Sum_j(beta_j * g_j(x_A)).
# 2) The zero-mean process eta has stationary covariance function specified by the input 'kernel', 
#    with constant variance sigma2.
# 3) The nugget accounts for the residual variability of f due to inactive inputs.
#    This is modeled via zero-mean uncorrelated errors with variance nu2.
# Prior mean and covariance for beta can be specified, otherwise a constant 0 regression term is considered.


### INPUTS (with default choices)
#
# ActInp.Design:  nxp matrix/data.frame, with values of the p active inputs at the n design points
#                    at which the simulator has been run.
# ActInp.Test:    Nxp matrix/data.frame, with values of the p active inputs at the N test points 
#                    at which the emulator predictions will be computed.
# y:              n-dim vector of simulator responses at the design points.
# Regress.Design: nxq matrix/data.frame, with values of the q regressors g_j() at the n design points.
#                    Default: ActiveInp.Design.
# Regress.Test:   Nxq matrix/data.frame, with values of the q regressors g_j() at the N test points.
#                    Default: ActiveInp.Test.
# beta:           q-dim vector of prior mean coefficients of regression term (1).
#                    Default: [0,...,0].
# Cov.beta:       qxq prior covariance matrix of regression term coefficients 
#                    Default: qxq 0-matrix.
# sigma2:         constant prior variance of zero-mean process eta(x_A).
#                    Default: var(y).
# kernel:         correlation kernel to use as prior correlation function of eta. 
#                    One of: 'exp2', 'abs_exp', 'matern32', 'matern52': see 'Corr.fun()' for more details.
#                    Default: 'exp2'.
# d:              p-dim vector of correlation lengths for correlation function of eta. See 'Corr.fun()'.
#                    Default: for each dimension, a third of the range covered by the design points.
# nu2:            variance of nugget process. Default: 0.
#
#
### OUTPUT
#
# An Nx2 matrix with named columns 'Mean' and 'Variance', containing the mean and variance of the emulator's
# predicted response at the N test points.


### NOTE
# The computation of mean and variance predictions is divided in blocks, to prevent memory problems
# when the number N of provided test points is large (order 10^6 or more)

################################################################################

BL.Emul <- function(ActInp.Design,
                    ActInp.Test,
                    y,
                    Regress.Design = ActInp.Design,
                    Regress.Test = ActInp.Test,
                    beta = rep(0, ncol(Regress.Design)),
                    Cov.beta = matrix(0, ncol(Regress.Design), ncol(Regress.Design)),
                    sigma2 = var(y),
                    kernel = 'exp2',
                    d = (apply(ActInp.Design, 2, max) - apply(ActInp.Design, 2, min)) /3,
                    nu2 = 0){
  
## PRELIMINARY DATA MANIPULATION AND VARIABLE DEFINITIONS ##
  
# Convert data frames into matrices if necessary
  ActInp.Test    <- as.matrix(ActInp.Test)
  Regress.Test   <- as.matrix(Regress.Test)
  ActInp.Design  <- as.matrix(ActInp.Design)
  Regress.Design <- as.matrix(Regress.Design)
  
# Define size variables
  N <- nrow(ActInp.Test)        # Size of Test Set
  n <- nrow(ActInp.Design)      # Size of Design Set

# Vectors which will contain posterior means and variances at test points
  Post.mean <- array(dim = c(N,1))          
  Post.var <-  array(dim = c(N,1))
  
## PRIOR MEAN AT DESIGN POINTS ##
  Prior.mean.Design <- Regress.Design %*% beta      # n-dim vector

## PRIOR COVARIANCE AT DESIGN POINTS ##
    
  # Prior covariance of regression part (the mean)
    a <- Regress.Design %*% Cov.beta %*% t(Regress.Design)         # nxn
  # Prior covariance of zero-mean process eta
    b <- sigma2*Corr.fun(ActInp.Design, ActInp.Design, d, kernel)  # nxn 
  # Prior covariance of nugget term (independent noise)
    c <- nu2*diag(n)                                               # nxn 
  # Total prior covariance at design points
    A <- a + b + c                                                 # nxn 
    
# THE FOLLOWING COMPUTATIONS ARE DIVIDED IN BLOCKS TO AVOID MEMORY PROBLEMS
# IF NUMBER OF TEST POINTS IS VERY LARGE 
  
  Block_size <- 2.e4                 # Size of each block of test points
  N_blocks <- ceiling(N/Block_size)  # Number of blocks to emulate
  
for (block in 1:N_blocks){
      
  # Start & end indices for the block to be considered, and length of the block
    ind_start <- (block-1)*Block_size +1
    ind_end   <- min( block*Block_size, N)
    ind <- ind_start:ind_end
    L   <- length(ind)
      
  # Select active inputs and regressors for test points in the current block
    ActiveInp.TestBlock <- ActInp.Test[ind, , drop=F]   # Lxp
    Regress.TestBlock   <- Regress.Test[ind, , drop=F]  # Lxq
    
  ## PRIOR MEAN AT TEST POINTS ## 
    Prior.mean.TestBlock <- Regress.TestBlock %*% beta  # Lx1
    
  ## PRIOR VARIANCE AT TEST POINTS ##
      
  # First compute the prior variance "Var.Reg" of the regression part:
  # Var.Reg = diag(Regress.SubsetTest %*% Cov.beta %*% Regress.SubsetTest).
  # To save time and memory, only diagonal terms are computed, not the whole matrix.
    K <- Regress.TestBlock %*% Cov.beta               # Lxq
    Var.Reg <- array(dim=c(L,1))                              
    for (j in 1:L){
      Var.Reg[j] <- K[j,] %*% Regress.TestBlock[j,]   # (1xq) x (qx1)
    }                           
  # Prior variance at test points is sum of "Var.Reg" and prior variance of eta
    Prior.var.SubsetTest <- Var.Reg + sigma2          # Lx1
    
  ## PRIOR COVARIANCE BETWEEN TEST AND DESIGN POINTS ##
      
  # Prior covariance due to regression term
    a <- Regress.TestBlock %*% Cov.beta %*% t(Regress.Design)           # Lxn = Lxq x qxq x qxn
  # Prior covariance due to eta term
    b <- sigma2*Corr.fun(ActiveInp.TestBlock, ActInp.Design, d, kernel) # Lxn
  # Total prior covariance t(x)
    tx <- a + b                                                         # Lxn
    
  ## FINAL COMPUTATIONS FOR EMULATOR: POSTERIOR MEAN AND VARIANCE AT TEST POINTS ##
    
  # Posterior mean at test points
    e <- y - Prior.mean.Design                        # n-dim vector
    K <- t( solve(A, t(tx)) )                         # Lxn: tx*(A^-1)
    Post.mean[ind,] <- Prior.mean.TestBlock + K%*%e   # Lx1
    
  # Posterior variance at test points
  # "a" is the "correction" term to be subtracted from the prior variance
    a <- array(dim=c(L,1))              # Lx1
    for (j in 1:L){
      a[j] <- K[j,] %*% tx[j,]          # (1xn) x (nx1)
    }
  # Posterior variance, add nu2 to account for variability due to inactive inputs
    Post.var[ind, ] <- Prior.var.SubsetTest  - a + nu2   # Nx1 matrix
  
  } # end for block=1:Nblocks
  
  return_matrix <- cbind("Mean" = c(Post.mean), "Var" = c(Post.var))
  return(return_matrix)
}


###############################################################################

# The function Corr.fun is used inside BL.Emul to compute prior correlation between inputs.

# X is an Nxp matrix, Y is an Mxp matrix. Correlation is computed between any row of X 
# and any row of Y, using the specified kernel. Output is an NxM matrix.
# The p correlation lengths used to compute the correlation function are specified in d. 

# INPUTS:
# X:      Nxp matrix or data-frame
# Y:      Mxp matrix or data frame
# d:      scalar or p-dim vector of positive correlation lengths. 
#         If scalar, it is extended to a vector (d,d,...,d) of length p.
# kernel: one of the following strings, to specify kernel to use:
#         'exp2'     (squared exponential)
#         'matern32' (matern 3/2, regularity: higher than C^1)
#         'matern52' (matern 5/2, regularity: higher than C^2)
#         'abs_exp'  (absolute exponential, continuoius but not C^1).
#
# OUTPUT:
# C: matrix of dimension NxM: C[i,j] is the correlation between X[i,] and Y[j,]

# COMPUTATIONAL NOTE: 
# The function uses 3D tensor multiplication to substantially increase speed, 
# rather than nested for loops in i&j: approx 1 order of magnitude of gain in speed 
# is observed for any two orders of magnitude increase in N*M. 

Corr.fun <- function(X, Y, d, kernel){
  
### DEFINE VARIABLES AND CHECK CORRECT INPUT DIMENSIONS ###
  
  #Convert dataframes into matrices if needed  
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  d <- as.numeric(d)
  
  # Store sizes
  N <- dim(X)[1]
  M <- dim(Y)[1]
  p <- dim(X)[2]
  
  # Stop if dimensions of inputs don't match
  if (dim(X)[2] != dim(Y)[2])
    stop('Inputs X and Y must have the same number of columns')
  
  # Stop if wrong number of correlation lengths is provided (or extend if just one is provided)  
  if (length(d)==1)
    d <- rep(d,p)
  else if (length(d) != p)
    stop('Length of input d must be equal to the number of columns of X and Y, or 1')

  
### CREATE 3D ARRAY Diff, with Diff[i,j,] = X[i,]-Y[j,] for i=1:N, j=1:M ###
  
  # Stack up copies of X and Y into 3D arrays
  X3 <- array(X, dim=c(dim(X),M) )  # dim = (N,p,M)
  Y3 <- array(Y, dim=c(dim(Y),N))   # (M,p,N)
  
  # Permute dimensions so that:
  # X3[i,j,:] is the p-dim vector X[i,:]
  # Y3[i,j,:] is the p-dim vector Y[j,:].
  X3 <- aperm(X3, c(1,3,2))         # (N,M,p)
  Y3 <- aperm(Y3, c(3,1,2))         # (N,M,p)  
  
  # Define Diff=X3-Y3, so that Diff[i,j,:] = X[i,]-Y[j,:]
  Diff <- X3-Y3                     # (N,M,p)
  # and rescale the k dimensions by correlation lengths
  for (i in 1:p)
    Diff[,,i] <- Diff[,,i]/d[i] 
  
  
### RETURN NxM MATRIX C OF CORRELATIONS, ACCORDING TO SPECIFIED KERNEL ###

  A <- rowSums(Diff^2, dims=2) # A[i,j] = sum of the k elements Diff[i,j,:]^2 
  
  if (identical(kernel,'exp2'))                 
    C <- exp(-A)
  else if (identical(kernel,'abs_exp'))         
    C <- exp(-sqrt(A))
  else if (identical(kernel,'matern32')){      
    A1 <- sqrt(3*A)
    C <- (1+ A1)*exp(-A1)
  }
  else if (identical(kernel,'matern52')){       
    A1 <- sqrt(5*A)
    A2 <- 5*A/3
    C <- (1+ A1 + A2)*exp(-A1)
  }
  else
    stop('Input \'string\' must be one of: \'exp2\', \'matern32\', \'matern52\', \'abs_exp\'.')
  
  return(C)
}
