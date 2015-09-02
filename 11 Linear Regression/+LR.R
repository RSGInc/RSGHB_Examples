# --------------------------------------------------------
# Code for estimating a multiple linear regression model
# 
# 
# Jeff Keller & Jeff Dumont
#
# --------------------------------------------------------

library(RSGHB)


# ------------------
# DATA PREPARATION
# ------------------

# Synthesize some data
set.seed(1987)
N <- 1000
x1 <- rnorm(N)
x2 <- rnorm(N)
#    B0  Bx1        Bx2
y <- 2 + 0.5 * x1 - 1.5 * x2 + rnorm(N)

# Assumes that observations are identified with an ID column
# In this case there is only one observation per unit/individual
ID <- data.frame(ID = seq_len(N))

# Specify any variables here that you'd like to use in the
# utility equations in the likelihood function below
# These can be any variables within the data or transformations of
# those variables


# ------------------------------------
# ESTIMATION CONTROL
# Setting control list for estimation
# ?doHB for more estimation options
# ------------------------------------

modelname <- "Linear_Regression"		# used for output

# Names for the random variables
gVarNamesFixed <- c("B0", "Bx1", "Bx2")

# STARTING VALUES
FC <- c(0, 0, 0)  # for the fixed coefficients
# The selection of the mean here is important when working with non-normal distributions

# ITERATION SETTINGS
gNCREP    <- 5000  # Number of iterations to use prior to convergence
gNEREP    <- 5000  # Number of iterations to keep for averaging after convergence has been reached
gNSKIP    <- 1		 # Number of iterations to do in between retaining draws for averaging
gINFOSKIP <- 100   # How frequently to print info about the iteration process


# CONTROL LIST TO PASS TO doHB
control <- list(
     modelname = modelname,
     gVarNamesFixed = gVarNamesFixed,
     FC = FC,
     gNCREP = gNCREP,
     gNEREP = gNEREP,
     gNSKIP = gNSKIP,
     gINFOSKIP = gINFOSKIP,
     gSeed = 1987
)

# ---------------------------------------------------------------------------------------
# likelihood
# USE:     Calculates the likelihood of y | x
#          Returns likelihood values for each observation
# NOTES:   This is where the bulk of the computation resides so coding this efficiently
#          is essential to reducing run time.
# ---------------------------------------------------------------------------------------
likelihood <- function(fc, b) {
     
     cc <- 1
     B0 <- fc[cc]; cc <- cc + 1
     Bx1 <- fc[cc]; cc <- cc + 1
     Bx2 <- fc[cc]; cc <- cc + 1
     
     # Number of parameters
     k <- length(fc)
     
     y_hat <- B0 + Bx1 * x1 + Bx2 * x2
     
     p <- dnorm(y - y_hat, mean = 0, sd = sqrt(sum((y - y_hat)^2)/(N - k)))
     
     return(p)
}

# Estimate the model
model <- doHB(likelihood, ID, control)

# Plot model statistics
plot(model)
plot(model, type = "F")

# Parameter estimates
colMeans(model[["F"]][, -1])
apply(model[["F"]][, -1], MARGIN = 2, FUN = sd)

# Compare to frequentist regression
summary(lm(y ~ x1 + x2))

# Save model object
save(model, file = paste0(model$modelname, ".RData"))

# Save in CSV format (Sawtooth-esque)
writeModel(model)
