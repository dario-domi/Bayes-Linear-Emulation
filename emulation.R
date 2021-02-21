### DESCRIPTION
#
# This function trains a Bayes Linear emulator of a computer simulator f, based on observed outputs y=(y1,...,yn) at n design points. 
# Emulation (ie, outputs prediction) is carried out at N Test points. The emulator form is based on the following model:
#
# f(x) = regression_term(x_A) + eta(x_A) + nugget_term(x),
#
# where x_A are the "active inputs", explaining most of the variability in f. 
# The three terms are as follows:
#
# 1) The regression term is a linear combination of basis functions g_i(x_A) (regressors):
#    regression_term(x_A) = Sum_i(beta_i * g_i(x_A)).
# 2) The zero-mean process eta has specified stationary covariance function with constant variance sigma2.
# 3) The nugget accounts for the residual variability due to inactive inputs.
#    This is modeled via zero-mean uncorrelated errors with variance nu.
# Prior mean and covariance for beta can be specified, otherwise a constant mean(y) regression term is considered.
#
#
### INPUTS  
#
# ActiveInp.Test:    Nxp matrix/data.frame, with values of the active inputs at the N Test points 
#                    at which the emulator is evaluated (number of active inputs = p)
# Regress.Test:      Nxq matrix/data.frame of regressors associated with the active inputs in Test.ActiveInp
#                    (number of regressors = q: g_1(x_A), ..., g_q(x_A))
# ActiveInp.Design:  nxp matrix/data.frame of active inputs at the n design points.
# Regress.Design:    nxq matrix/data.frame of regressors associated with ActiveInp.Design
# y:                 n-dim vector of simulator responses at the design points
# beta:              q-dim vector of prior mean coefficients of regression term
# Cov.beta:          qxq prior covariance matrix of regression term coefficients
# sigma2:            prior variance of zero-mean process eta(x_A) at all x_A
# corr.string:       correlation function to use for eta. One of: 'exp2', 'abs_exp', 'matern32', 'matern52'
#                    See 'Corr.fun' for more details
# d:                 p-dim vector of correlation lengths for the covariance function of eta. See Corr.fun
# nu:                variance of nugget process
#
#
### OUTPUT
#
# An Nx2 matrix, with named columns 'Mean' and 'Variance'. These contain the mean and variance of the emulator
# predicted response at the N Test points.
#
#
# NOTE: The computation of mean and variance predictions is divided in blocks, to prevent memory problems
#       when the number of Test inputs is large (order 10^6 or more)


BL.Emul <- function(ActiveInp.Test, 
                 ActiveInp.Design, 
                 y, 
                 Regress.Test = ActiveInp.Test, 
                 Regress.Design = ActiveInp.Design, 
                 beta = rep(0, ncol(Regress.Design)), 
                 Cov.beta = matrix(0, ncol(Regress.Design) , ncol(Regress.Design)), 
                 sigma2 = var(y), 
                 corr.string = 'exp2', 
                 d = (apply(ActiveInp.Design, 2, max) - apply(ActiveInp.Design, 2, min))/3, 
                 nu = 0){
  
## PRELIMINARY DATA MANIPULATION AND VARIABLE DEFINITIONS ##
  
  # Convert data frames into matrices if necessary
    ActiveInp.Test <- as.matrix(ActiveInp.Test)
    Regress.Test <- as.matrix(Regress.Test)
    ActiveInp.Design <- as.matrix(ActiveInp.Design)
    Regress.Design <- as.matrix(Regress.Design)
  
  # Define size variables
    N <- nrow(ActiveInp.Test)        # Size of Test Set
    n <- nrow(ActiveInp.Design)      # Size of Design Set

  # Vectors which will contain posterior means and variances at Test points
    Post.mean <- array(dim = c(N,1))          
    Post.var <-  array(dim = c(N,1))
  
## PRIOR MEAN AT DESIGN POINTS ##
    
    Prior.mean.Design <- Regress.Design %*% beta                         # n-dim vector
  
## PRIOR COVARIANCE AT DESIGN POINTS ##
    
  # Prior covariance of regression part (the mean)
    a <- Regress.Design %*% Cov.beta %*% t(Regress.Design)               # nxn
  # Prior covariance of zero-mean process eta
    b <- sigma2*Corr.fun(ActiveInp.Design, ActiveInp.Design, d, corr.string)  # nxn 
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
        ActiveInp.SubsetTest <- ActiveInp.Test[ind, , drop=F]
        Regress.SubsetTest   <- Regress.Test[ind, , drop=F]
    
    ## PRIOR MEAN AT TEST POINTS ## 
        Prior.mean.SubsetTest <- Regress.SubsetTest %*% beta        # Lx1
    
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
        b <- sigma2*Corr.fun(ActiveInp.SubsetTest, ActiveInp.Design, d, corr.string)  # Lxn
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
