# --------------------------------------------------------
# Code for a Latent Class Model
# 
# 
# Jeff Dumont & Jeff Keller
#
# --------------------------------------------------------

library(RSGHB)

# ------------------
# DATA PREPARATION
# ------------------

# assumes that respondents are identified with a ID column
# also assumes that the data is sorted by respondent then experiment
choicedata <- read.table("lc_simulated.csv", sep = ",", header = TRUE)

# Specify any variables here that you'd like to use in the
# utility equations in the likelihood function below
# These can be any variables within the data or transformations of
# those variables
attach(choicedata)
brand_alt1_level1 <- ( brand_alt1 == 1 ) * 1 
brand_alt2_level1 <- ( brand_alt2 == 1 ) * 1 
brand_alt1_level2 <- ( brand_alt1 == 2 ) * 1 
brand_alt2_level2 <- ( brand_alt2 == 2 ) * 1 
brand_alt1_level3 <- ( brand_alt1 == 3 ) * 1 
brand_alt2_level3 <- ( brand_alt2 == 3 ) * 1 

# The choice vectors
# Dummying coding the choice vector allows for easier coding of the 
# the likelihood calculations. So we will have one column for each 
# alternative in the design
Choice <-  choice
Choice1 <- ( choice == 1 )
Choice2 <- ( choice == 2 )
Choice3 <- ( choice == 3 )

# ------------------------------------
# ESTIMATION CONTROL
# Setting control list for estimation
# ?doHB for more estimation options
# ------------------------------------

modelname <- "Latent_Class"	# used for output

# Names for the Fixed Variables
gVarNamesFixed <- c("BBrand_2_class1", "BBrand_3_class1", "Bfeature_class1", "Bprice*scale_class1", "Bnone_class1",
                    "BBrand_2_class2", "BBrand_3_class2", "Bfeature_class2", "Bprice*scale_class2", "Bnone_class2", "delta2")

# STARTING VALUES
FC <- rep(0, length(gVarNamesFixed))  # for the fixed coefficients
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

likelihood <- function(fc,b) {
     
     cc <- 1 
     BBrand_2_class1     <- fc[cc]; cc <- cc + 1
     BBrand_3_class1     <- fc[cc]; cc <- cc + 1
     Bfeature_class1     <- fc[cc]; cc <- cc + 1
     Bprice_scale_class1 <- fc[cc]; cc <- cc + 1
     Bnone_class1        <- fc[cc]; cc <- cc + 1

     BBrand_2_class2     <- fc[cc]; cc <- cc + 1
     BBrand_3_class2     <- fc[cc]; cc <- cc + 1
     Bfeature_class2     <- fc[cc]; cc <- cc + 1
     Bprice_scale_class2 <- fc[cc]; cc <- cc + 1
     Bnone_class2        <- fc[cc]; cc <- cc + 1
     
     delta2 <- fc[cc]; cc <- cc + 1
          
     # Alternative utilities for Class 1
     V1.class1 <- Bprice_scale_class1 * ( BBrand_2_class1 * brand_alt1_level2 + BBrand_3_class1 * brand_alt1_level3 + Bfeature_class1 * feature_alt1 + price_alt1 )
     V2.class1 <- Bprice_scale_class1 * ( BBrand_2_class1 * brand_alt2_level2 + BBrand_3_class1 * brand_alt2_level3 + Bfeature_class1 * feature_alt2 + price_alt2 )
     V3.class1 <- Bprice_scale_class1 * Bnone_class1 
     
     # Alternative utilities for Class 2
     V1.class2 <- Bprice_scale_class2 * ( BBrand_2_class2 * brand_alt1_level2 + BBrand_3_class2 * brand_alt1_level3 + Bfeature_class2 * feature_alt1 + price_alt1 )
     V2.class2 <- Bprice_scale_class2 * ( BBrand_2_class2 * brand_alt2_level2 + BBrand_3_class2 * brand_alt2_level3 + Bfeature_class2 * feature_alt2 + price_alt2 )
     V3.class2 <- Bprice_scale_class2 * Bnone_class2 
     
     # Choice probabilities given class
     P.class1 <- ( Choice1 * exp(V1.class1) + Choice2 * exp(V2.class1) + Choice3 * exp(V3.class1) ) / ( exp(V1.class1) + exp(V2.class1) + exp(V3.class1) )
     P.class2 <- ( Choice1 * exp(V1.class2) + Choice2 * exp(V2.class2) + Choice3 * exp(V3.class2) ) / ( exp(V1.class2) + exp(V2.class2) + exp(V3.class2) )

     # aggregate() is very slow. This likelihood function could be sped up by
     # using something like data.table for the class probability aggregation by ID (not shown here)
     P.class1 <- aggregate(P.class1, by = list(ID), prod)[, 2]
     P.class2 <- aggregate(P.class2, by = list(ID), prod)[, 2]
     
     # Class membership - could be a function of additional socio-demographic variables (not shown here)
     # and not just the constant delta2
     P.membership1 <- 1/(1 + exp(delta2))
     P.membership2 <- exp(delta2)/(1 + exp(delta2))
     
     P <- P.class1 * P.membership1 + P.class2 * P.membership2
     
     P <- P[ID]^(1/9)
     
     return(P)
     
}

# Estimate the model
model <- doHB(likelihood, choicedata, control)

# Plot model statistics
plot(model)
plot(model, type = "F")

# Save model object
save(model, file = paste0(model$modelname, ".RData"))

# Save in CSV format (Sawtooth-esque)
writeModel(model)
