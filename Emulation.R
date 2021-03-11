### DESCRIPTION
#
# This function trains a Bayes Linear emulator of a computer simulator f, based on observed outputs y=(y1,...,yn) at n design points. 
# Emulation (ie, output prediction) is carried out at N Test points. 
# Full details on the meaning of the function's arguments available at: 
# https://github.com/dario-domi/Bayes-Linear-Emulation/blob/main/Documentation.pdf
#
# The emulator form is based on the following model:
#
# f(x) = regression_term(x_A) + eta(x_A) + nugget_term(x),
#
# where x_A are the "active inputs", explaining most of the variability in f. 
# The three terms are as follows:
#
# 1) The regression term is a linear combination of basis functions g_j(x_A) (regressors):
#    regression_term(x_A) = Sum_i(beta_j * g_j(x_A)).
# 2) The zero-mean process eta has stationary covariance function specified in input 'kernel', with constant variance sigma2.
# 3) The nugget accounts for the residual variability of f due to inactive inputs.
#    This is modeled via zero-mean uncorrelated errors with variance nu.
# Prior mean and covariance for beta can be specified, otherwise a constant 0 regression term is considered.
#
#
### INPUTS  
#
# ActInp.Design:  nxp matrix/data.frame, with values of the active inputs at the n design points
#                    (number of active inputs = p).
# ActInp.Test:    Nxp matrix/data.frame, with values of the active inputs at the N Test points 
#                    at which the emulator is evaluated.
# y:                 n-dim vector of simulator responses at the design points.
# Regress.Design:    nxq matrix/data.frame, with values of the q regressors g_j()
#                    evaluated at each of the n rows in ActiveInp.Design
#                    Default=ActiveInp.Design.
# Regress.Test:      Nxq matrix/data.frame, with values of the q regressors g_j()
#                    evaluated at each of the N rows in ActiveInp.Test 
#                    Default=ActiveInp.Test.
# beta:              q-dim vector of prior mean coefficients of regression term 
#                    Default: [0,...,0].
# Cov.beta:          qxq prior covariance matrix of regression term coefficients 
#                    Default: qxq 0-matrix.
# sigma2:            prior variance of zero-mean process eta(x_A) at all x_A
#                    Default: var(y).
# kernel:            correlation kernel to use for eta. One of: 'exp2', 'abs_exp', 'matern32', 'matern52'.
#                    See 'Corr.fun()' for more details. Default: 'exp2'.
# d:                 p-dim vector of correlation lengths for the covariance function of eta. See 'Corr.fun()'.
# nu:                variance of nugget process. Default: 0.
#
#
### OUTPUT
#
# An Nx2 matrix, with named columns 'Mean' and 'Variance'. These contain the mean and variance of the emulator's
# predicted response at the N Test points.
#
#
# NOTE: The computation of mean and variance predictions is divided in blocks, to prevent memory problems
#       when the number of Test inputs is large (order 10^6 or more)


BL.Emul <- function(ActInp.Design,
                    ActInp.Test,
                    y,
                    Regress.Design = ActInp.Design,
                    Regress.Test = ActInp.Test,
                    beta = rep(0, ncol(Regress.Design)),
                    Cov.beta = matrix(0, ncol(Regress.Design) , ncol(Regress.Design)),
                    sigma2 = var(y),
                    kernel = 'exp2',
                    d = (apply(ActInp.Design, 2, max) - apply(ActInp.Design, 2, min)) /3,
                    nu = 0){
  
## PRELIMINARY DATA MANIPULATION AND VARIABLE DEFINITIONS ##
  
  # Convert data frames into matrices if necessary
    ActInp.Test <- as.matrix(ActInp.Test)
    Regress.Test <- as.matrix(Regress.Test)
    ActInp.Design <- as.matrix(ActInp.Design)
    Regress.Design <- as.matrix(Regress.Design)
  
  # Define size variables
    N <- nrow(ActInp.Test)        # Size of Test Set
    n <- nrow(ActInp.Design)      # Size of Design Set

  # Vectors which will contain posterior means and variances at Test points
    Post.mean <- array(dim = c(N,1))          
    Post.var <-  array(dim = c(N,1))
  
## PRIOR MEAN AT DESIGN POINTS ##
    
    Prior.mean.Design <- Regress.Design %*% beta                         # n-dim vector
  
## PRIOR COVARIANCE AT DESIGN POINTS ##
    
  # Prior covariance of regression part (the mean)
    a <- Regress.Design %*% Cov.beta %*% t(Regress.Design)               # nxn
  # Prior covariance of zero-mean process eta
    b <- sigma2*Corr.fun(ActInp.Design, ActInp.Design, d, kernel)  # nxn 
  # Prior covariance of nugget term (independent noise)
    c <- nu*diag(n)                                                      # nxn 
  # Total prior covariance at design points
    A <- a + b + c                                                       # nxn 
    
# TO AVOID MEMORY PROBLEMS IF NUMBER OF TEST POINTS IS VERY LARGE, 
# THE FOLLOWING COMPUTATIONS ARE DIVIDED IN BLOCKS
  
    N_block <- 2.e4                 # Size of each block of Test points
    N_loops <- ceiling(N/N_block)   # Number of iterations consequently needed
  
    for (i in 1:N_loops){
      
      # Start & end indices for the block to be considered, and length of the block
        ind_start <- (i-1)*N_block +1
        ind_end   <- min( i*N_block, N)
        ind <- ind_start:ind_end
        L   <- length(ind)
      
      # Select active inputs and regressors for Test points in the current block
        ActiveInp.SubsetTest <- ActInp.Test[ind, , drop=F]
        Regress.SubsetTest   <- Regress.Test[ind, , drop=F]
    
    ## PRIOR MEAN AT TEST POINTS ## 
        Prior.mean.SubsetTest <- Regress.SubsetTest %*% beta  # Lx1
    
    ## PRIOR VARIANCE AT TEST POINTS ##
      
      # First compute the prior variance "Var.Reg" of the regression part:
      # Var.Reg = diag(Regress.SubsetTest %*% Cov.beta %*% Regress.SubsetTest).
      # To save time and memory, we only compute those diagonal terms, not the whole matrix.
        K <- Regress.SubsetTest %*% Cov.beta                  # Lxq
        Var.Reg <- array(dim=c(L,1))                              
        for (j in 1:L){
          Var.Reg[j] <- K[j,] %*% Regress.SubsetTest[j,]      # (1xq) x (qx1)
        }                           
      # Prior variance at Test points is sum of "Var.Reg" and prior variance of eta
        Prior.var.SubsetTest <- Var.Reg + sigma2              # Lx1
    
    ## PRIOR COVARIANCE BETWEEN TEST AND DESIGN POINTS ##
      
      # Prior covariance due to regression term
        a <- Regress.SubsetTest %*% Cov.beta %*% t(Regress.Design)               # Lxn = Lxq x qxq x qxn
      # Prior covariance due to eta term
        b <- sigma2*Corr.fun(ActiveInp.SubsetTest, ActInp.Design, d, kernel)  # Lxn
      # Total prior covariance t(x)
        tx <- a + b                                                              # Lxn
    
    ## FINAL COMPUTATIONS FOR EMULATOR: POSTERIOR MEAN AND VARIANCE AT TEST POINTS ##
    
      # Posterior mean at test points
        e <- y - Prior.mean.Design                          # n-dim vector
        K <- t( solve(A, t(tx)) )                           # Lxn: tx*(A^-1)
        Post.mean[ind,] <- Prior.mean.SubsetTest + K%*%e    # Lx1
    
      # Posterior variance at test points
      # "a" is the "correction" term to be subtracted from the prior variance
        a <- array(dim=c(L,1))                              # Lx1
        for (j in 1:L){
          a[j] <- K[j,] %*% tx[j,]                          # (1xn) x (nx1)
        }
      # Posterior variance, add nu to account for variability due to inactive inputs
        Post.var[ind, ] <- Prior.var.SubsetTest  - a + nu   # Nx1 matrix
    }
  
  return_matrix <- cbind("Mean" = c(Post.mean), "Var" = c(Post.var))
  return(return_matrix)
}
