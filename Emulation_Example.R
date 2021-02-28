source('Corr_Fun.R')
source("Emulation.R")
library("plot3D")
library("randtoolbox")


############################################################
########     Definition of f and plots      ################
############################################################

## Let's define a toy (non-linear) function f over a 2D square. 
## We'll pretend this is costly and slow to evaluate, so we can only observe its 
## outputs y at a few points of the 2D square. We will then emulate it. 

# Define f
f <- function(x) cos(pi*x[1]) + 2*x[1]*exp(x[2]) + x[2]^2
# Sample n random points in the square [-1,1]x[-1,1]
n <- 50
Design <- 2*sobol(n, dim=2, scrambling = 3) -1
# Convert into a data.frame (though not necessary to run BL.Emul())
Design <- data.frame('x1'=Design[,1], 'x2'=Design[,2] )
# Evaluate g at the points
y <- apply(Design, 1, f)
# Show scatter plot of the values of y
scatter2D(Design[,1], Design[,2], colvar=y, pch=20, cex=1)

## If you want to see how the actual function looks like, run the following
L = 200
x = seq(-1, 1, length=L)
xlist <- expand.grid(x, x)
ylist <- apply(xlist, MARGIN = 1, FUN = f)
scatter2D(xlist[,1], xlist[,2], colvar = ylist, pch=20, cex=1.5)


####################################################################
########        Emulation of f with default options     ############
####################################################################

## We now emulate f, leaving all optional arguments of the function BL.Emul()
## to their defaults

# Sample N points at random in the square [-1,1]x[-1,1]
N <- 10000
Test <- 2*matrix(runif(N*2), N, 2) -1
Test <- data.frame('x1'=Test[,1], 'x2'=Test[,2])
# Emulate f at these points
Em.def <- BL.Emul(ActInp.Design = Design, ActInp.Test = Test, y = y)



####################################################################
########      Emulation of f with customised options     ###########
####################################################################

# We will now customise a bit our emulator. To guess which regressors may play a role, 
# let's run a linear model with interations.
mod.full <- lm(y ~ x1*x2, data = Design)
summary(mod.full)

## For this simple illustration, given the above result, we decide to:
## 1) select x[1] and x[1]*x[2] as regressors, estimating the associated beta and Cov.beta from linear regression;
## 2) compute the variability in y once the contribution from the regressors has been taken into account,
##    and assign it as sigma2 (the prior variance of the process eta)
## 3) given the smoothness of f, increase the correlation lengths to 0.75.
## Note that, in this case, we know y is *only* a function of the active inputs x[1] and x[2],
## hence we won't include a nugget term in the emulator (ie, leave nu=0).


# Select x1 and x1*x2 as regressors
mod <- lm(y ~ x1 + x1:x2 -1, data = Design)
summary(mod)

# Formula for the built linear model
forml <- formula(mod)[c(1,3)]

# Generate regressors for Design and Test points according to above formula
Regr.Design <- model.matrix(forml, data = Design)
Regr.Test   <- model.matrix(forml, data = Test)

# Store estimates of beta and variance of residuals, to be fitted through the eta component of the emulator
b <- mod$coefficients
Cov.b <- vcov(mod)
s2 <- var(mod$residuals)

# Evaluate emulator
Em.cust <- BL.Emul(ActInp.Design = Design, ActInp.Test = Test,
                   Regress.Design = Regr.Design, Regress.Test = Regr.Test,
                   y = y, beta = b, Cov.beta = Cov.b,
                   sigma2 = s2, d=c(0.75, 0.75))


#############################################################################
#################      Goodness of each emulator         ####################
#############################################################################

# Let's compare the predictions of each emulator 
# with the actual values of f at the Test points.
true.resp <- apply(Test, 1, f)

# Do the following once with Emul <- Em.def, once with Emul <- Emul.cust
Emul <- Em.def
Emul <- Em.cust
# Raw emulator errors
Raw.Errs <- Emul[,1]-true.resp
# Errors standardised by the prediction standard deviation
Std.Errs <- Raw.Errs/sqrt(Emul[,2])

# Histograms of errors
hist(Raw.Errs, 200, main = 'Raw Errors', xlim = c(-0.2, 0.2))
hist(Std.Errs, 50, main = 'Standardised Errors')

cat('The average standard deviation of predictions across the domain is',
    mean(sqrt(Emul[,2])))

cat('95% of raw errors in emulator predictions are between',
    quantile(Raw.Errs, 0.025), 'and', quantile(Raw.Errs, 0.975))

cat('95% of standardised errors in emulator predictions are between',
    quantile(Std.Errs, 0.025), 'and', quantile(Std.Errs, 0.975), '.',
    '\nGiven Pukelseim theorem, we expect these to be contained in the interval [-3,3].')


# IMPORTANT: 
# Observe how areas of low uncertainty in the emulator prediction 
# correspond to areas where the value of f has been observed
scatter2D(Test[,1], Test[,2], colvar=log(Emul[,2]), pch=20)
scatter2D(Design[,1], Design[,2], col = 'black', pch=20, add=T)



##########################################################
################      CONCLUSIONS        #################
##########################################################


# Notice how both emulators, even the one with no custom choices, are very precise
# in approximating the (highly non-linear) function f.
# However, you should see that the second emulator is notably more accurate (lower Raw.Errs)
# and more certain (lower standard deviations) than the first emulator, as it is to be expected.
