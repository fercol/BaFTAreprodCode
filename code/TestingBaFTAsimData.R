# ============================== CODE METADATA =============================== #
# AUTHOR: Fernando Colchero
# DATE CREATED: 2024-12-01
# DESCRIPTION: Testing BaFTA on simulated data.
# NOTES: 
# 1.- Data types: 
#   a) 'Aggr': aggregated data
#   b) 'Simp': Individual data with discrete ages (seasonal)
#   c) 'Ext': Individual data with variable interbirth interval (IBI)
# 2.- Species type:
#   a) 'Fast': low survival, high reproductive output
#   b) 'Slow': high survival, low reproductive output
# 3.- Uncertainty types:
#   a) 'Reg': no uncertainty
#   b) 'UnkAge': some individuals have unknown age
#   c) 'UnkObs': some parents were not followed constantly and the number of
#                offspring can be a fraction of the total number of offspring.
# ================================ CODE START ================================ #
# ==================== #
# ==== LOAD DATA: ====
# ==================== #
# Installed packages:
instPacks <- installed.packages()[, 1]

# Libraries:
if (!"snowfall" %in% instPacks) {
  install.packages("snowfall")
}

if (!"RColorBrewer" %in% instPacks) {
  install.packages("RColorBrewer")
}

if (!"BaFTA" %in% instPacks) {
  if (!"devtools" %in% instPacks) {
    install.packages("devtools")
  }
  library(devtools)
  
  install_git("https://github.com/fercol/BaFTA", subdir = "pkg/")
}

# Load libraries:
library(snowfall)
library(BaFTA)
library(RColorBrewer)

# Set working directory:
setwd("~/FERNANDO/PROJECTS/1.ACTIVE/FertilityEstimation/analysis/ReprodCode/")

# Source additional functions:
source("code/TestingBaFTAextraFunctions.R")

# Logical to save plots:
saveResults <- FALSE

# ============================ #
# ==== SELECT DATA TYPES: ====
# ============================ #
# Data type (options: 'Simp', 'Ext', 'Aggr'):
dType <- "Ext"

# species type (options: 'Slow', 'Fast'):
spType <- "Slow"

# Uncertainty type (options: 'Reg', 'UnkAge', 'UnkObs'):
uType <- "Reg"

# Number of individuals:
ndat <- 100

# Prop unk ages:
pUnkAge <- 0.5

# Proportion of uncertain births:
pUnkOffs <- 0.25

# ===================== #
# ==== PREP. DATA: ====
# ===================== #
# Full individual data:
if (dType %in% c("Simp", "Aggr")) {
  # Load full simulate data:
  load(sprintf("data/simDataSimp%sSps.RData", spType))
  
  # Simulation parameter values:
  simPars <- c(simSettings$beta, gamma = simSettings$gamma)
  
} else {
  load(sprintf("data/simDataExt%sSps.RData", spType))
  
  # Simulation parameter values:
  simPars <- c(simSettings$beta, gamma = simSettings$gamma, 
               eta = simSettings$eta, kappa = simSettings$kappa, 
               vSd = simSettings$vSd)
}

# Number of simulation parameters:
nthe <- length(simPars)

# Minimum and maximum ages:
minAge <- simSettings$alpha
maxAge <- simSettings$omega

# Total number of records:
nr <- nrow(reprdat)

# Individual IDs:
indID <- unique(reprdat$indID)

# total number of individuals:
nind <- length(indID)

# Individual data:
indsamp <- indID[sample(1:nind, size = ndat, replace = FALSE)]
idind <- which(reprdat$indID %in% indsamp)
anDat <- reprdat[idind, ]
anDat$indID <- as.character(anDat$indID)

# Aggregated data or unknow ages/offspring:
if (dType == "Aggr") {
  # ---------------- #
  # Aggregated data:
  # ---------------- #
  # Find individuals in longevity data:
  idlong <- which(longdat$indID %in% anDat$indID)
  
  # Create vector of discrete ages:
  ages <- minAge:c(min(c(maxAge, floor(max(anDat$Age)))))
  
  # length of age vector:
  nages <- length(ages)
  
  # Temporal empty table for aggregated data:
  tempDat <- data.frame(Age = ages, nParents = rep(0, nages), 
                         nOffspring = rep(0, nages), Fertility = rep(0, nages))
  
  # Populate the table:
  for (ia in 1:nages) {
    # Find available parents at age:
    idpar <- which(longdat$ageEntry[idlong] <= ages[ia] + 1 & 
                     longdat$ageDepart[idlong] >= ages[ia])
    
    # Number of available parents at that age:
    npars <- length(idpar)
    
    # Number of offspring produced at that age:
    idoff <- which(floor(anDat$Age) == ages[ia])
    noffs <- sum(anDat$nOffspring[idoff])
    
    # Fill up table:
    tempDat$nParents[ia] <- npars
    tempDat$nOffspring[ia] <- noffs
    tempDat$Fertility[ia] <- noffs / npars
  }
  # Avoid NAs:
  tempDat$Fertility[which(is.na(tempDat$Fertility))] <- 0
  
  # Update analysis data table:
  anDat <- tempDat
} else if (uType == "UnkAge") {
  # ---------------------------------- #
  # Individual data with unknown ages:
  # ---------------------------------- #
  # Number of individuals with unknown age:
  nUnkAge <- round(ndat * pUnkAge)
  
  # Select individuals with unknown age:
  iduage <- sort(sample(1:ndat, size = nUnkAge))
  
  # Extract unique individual IDs in sampled data:
  sunID <- unique(anDat$indID)
  
  # Real ages:
  realAges <- rep(NA, length(iduage))
  
  # Use built-in function to calculate min-max ages for those individuals:
  anDat <- CalcMinMaxAge(object = anDat, minAge = minAge, 
                         maxAge = maxAge, 
                         unkAgeID = sunID[iduage])
  
  # Randomize assigned age by the user:
  for (ii in 1:nUnkAge) {
    # Find all records for the ii individual:
    idi <- which(anDat$indID == sunID[iduage[ii]])
    
    # Fill first real age:
    realAges[ii] <- min(anDat$Age[idi])
    
    # Find range of delta value based on min max age for that individual:
    updRan <- c(anDat$MinAge[idi[1]], anDat$MaxAge[idi[1]]) -
      anDat$Age[idi[1]]
    
    # Draw a random delta value for that individual:
    if (dType == "Ext") {
      deltaAge <- rtnorm(n = 1, mean = 0, sd = 1, lower = updRan[1],
                          upper = updRan[2])
    } else {
      deltaAge <- rtnorm(n = 1, mean = 0, sd = 1, lower = updRan[1],
                          upper = updRan[2])
    }
    
    # Update assigned age for that individual:
    anDat$Age[idi] <- anDat$Age[idi] + deltaAge
  }
  
} else if (uType == "UnkObs" & dType == "Simp") {
  # --------------------------------------------------- #
  # Individual data with uncertain number of offspring:
  # --------------------------------------------------- #
  # number of individual records in analysis data:
  nsub <- nrow(anDat)
  
  # Number of partial observations in analysis data:
  npartial <- floor(nsub * pUnkOffs)
  
  # Select records with partial observations:
  idpartial <- sort(sample(1:nsub, size = npartial, replace = FALSE))
  
  # Real offspring:
  realOffs <- anDat$nOffspring[idpartial]
  
  # Create vector of observation proportions (1 for fully observed, < 1 otherw.)
  obsProp <- rep(1, nsub)
  
  # Fill up vector of observation proportions for partially observed inds.:
  obsProp[idpartial] <- sample(c(0.25, 0.5, 0.75), size = length(idpartial),
                               replace = TRUE)
  
  # Collate new observation proportion to analysis data frame:
  anDat <- data.frame(anDat, obsProp = obsProp)
  
  # Update number of observed offspring for partial observations:
  anDat$nOffspring[idpartial] <- 
    rbinom(n = npartial, size = anDat$nOffspring[idpartial], 
           prob = obsProp[idpartial])
  
} else if (uType == "UnkObs" & dType == "Ext") {
  stop("Unknown observations and uncertain offspring only implemented for
       dType == 'Simp'")
}

# ============================= #
# ==== RUN BaFTA ANALYSES: ====
# ============================= #
# data type for main BaFTA function:
if (dType == "Ext") {
  dataType <- "indivExtended"
} else if (dType == "Simp") {
  dataType <- "indivSimple"
} else if (dType == "Aggr") {
  dataType <- "aggregated"
}

# -------------------------- #
# ---- Quadratic model: ----
# -------------------------- #
# run analysis:
out0 <- bafta(object = anDat, dataType = dataType, model = "quadratic",
              niter = 20000, burnin = 1001, nsim = 4, ncpus = 4)

# print output to the screen:
out0

# plot traces:
plot(out0)

# Plot parameter posterior densities:
plot(out0, type = "density")

# plot fertility:
plot(out0, type = "fertility")

# plot predicted number of offspring:
plot(out0, type = "predictive")

# ---------------------- #
# ---- Gamma model: ----
# ---------------------- #
# Run analysis:
out1 <- bafta(object = anDat, dataType = dataType, model = "gamma",
              niter = 20000, burnin = 1001, nsim = 4, ncpus = 4)

# print output to the screen:
out1

# plot traces:
plot(out1)

# Plot parameter posterior densities:
plot(out1, type = "density")

# plot fertility:
plot(out1, type = "fertility")

# plot predicted number of offspring:
plot(out1, type = "predictive")

# Combined plot of posterior densities and traces:
PlotDensTrace(out = out1, simPars = simPars)
