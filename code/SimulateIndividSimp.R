# ============================== CODE METADATA =============================== #
# AUTHOR: Fernando Colchero
# DATE CREATED: 2024-12-19
# DESCRIPTION: Code to simulate individual level fertility data on discrete 
#              ages to test BaFTA performance.
# ================================ CODE START ================================ #
# ======================== #
# ==== GENERAL SETUP: ====
# ======================== #
# Installed packages:
instPacks <- installed.packages()[, 1]

# Libraries:
if (!"paramDemo" %in% instPacks) {
  if (!"devtools" %in% instPacks) {
    install.packages("devtools")
  }
  library(devtools)
  
  install_git("https://github.com/fercol/paramDemo/", subdir = "pkg/")
}

# Load libraries:
library(paramDemo)

# Set working directory:
setwd("~/FERNANDO/PROJECTS/1.ACTIVE/FertilityEstimation/analysis/ReprodCode/")

# ==================== #
# ==== FUNCTIONS: ====
# ==================== #
# Function to construct Leslie matrix:
FillMatr <- function(p, f, n) {
  idcol <- 1:n
  idrow <- c(2:n, n)
  idpa <- (idcol - 1) * n + idrow
  Aa <- matrix(0, n, n)
  Aa[1, ] <- f
  Aa[idpa] <- p
  return(Aa)
}

# ============================== #
# ==== PARAMETRIC SETTINGS: ====
# ============================== #
# Choose species type (i.e., Fast vs Slow):
spType <- "Slow"

# Initial population size:
N <- 2000

# Number of years to simulate:
Nyears <- 44
Tmax <- Nyears

if (spType == "Slow") {
  # Siler mortality parameters:
  theta <- c(a0 = -1, a1 = 2, c = 0.0001, b0 = -6, b1 = 0.15)
  
  # Gamma fertility parameters:
  beta <- c(b0 = 4, b1 = 2, b2 = 0.15)
  
  # Parameter for individual effects on fertility:
  gamma <- 10
  
  # Eta inter-birth interval parameter:
  eta <- 2
  
  # Gamma parameter for age at first birth:
  kappa <- 2
  
  # Standard error of IBI individual effects:
  vSd <- 0.2
  
  # Age at maturity:
  alpha <- 4
  
  # Gestation and weaning period:
  gest <- 0.5
  wean <- 0.5
  tau <- gest + wean
} else {
  # Siler mortality parameters:
  theta <- c(a0 = 0, a1 = 4, c = 0.5, b0 = -2, b1 = 0.1)
  
  # Gamma fertility parameters:
  beta <- c(b0 = 40, b1 = 2, b2 = 0.5)
  
  # Parameter for individual effects on fertility:
  gamma <- 20
  
  # Eta inter-birth interval parameter:
  eta <- 10
  
  # Gamma parameter for age at first birth:
  kappa <- 10
  
  # Standard error of IBI individual effects:
  vSd <- 0.2
  
  # Age at maturity:
  alpha <- 1
  
  # Gestation and weaning period:
  gest <- 0.06
  wean <- 0.07
  tau <- gest + wean
}

# Vector of ages:
dx <- 0.001
x <- seq(0, 100, dx)

# ============================ #
# ==== VISUAL INSPECTION: ====
# ============================ #
# Calculate mortality and survival:
Sx <- CalcSurv(theta = theta, x = x, model = "GO", shape = "bathtub")
idsx <- which(Sx >= 0.001)
xv <- x[idsx]
omega <- floor(max(xv))
Sx <- Sx[idsx]
mux <- CalcMort(theta = theta, x = xv, model = "GO", shape = "bathtub")

# Calculate age-specific fertility:
gx <- CalcFert(beta = beta, x = xv[which(xv >= alpha)] - alpha, 
               modelFert = "gamma")

# Age at minimum mortality:
xmin <- xv[which(mux == min(mux))]

# Plot vital rates:
par(mfrow = c(3, 1), mar = c(4, 4, 1, 1))
plot(xv, mux, type = 'l', col = 'dark red', lwd = 4, xlab = "", 
     ylab = "Mortality")
abline(v = alpha, col = 'orange')
abline(v = xmin, col = 'orange', lty = 2)
plot(xv, Sx, type = 'l', col = 'dark red', lwd = 4, xlab = "", 
     ylab = "Survival")
abline(v = alpha, col = 'orange')
abline(v = xmin, col = 'orange', lty = 2)
plot(xv[which(xv >= alpha)], gx, type = 'l', col = 'dark red', lwd = 4, 
     xlab = "Age", ylab = "Fertility", xlim = c(0, max(xv)), 
     ylim = c(0, max(gx)*qgamma(0.99, shape = gamma, rate = gamma)))
for (ii in c(-1, 1)) {
  lines(xv[which(xv >= alpha)], 
        gx * qgamma(c(0.025, 0.975)[ii], shape = gamma, rate = gamma), 
        col = 'dark red')
}
abline(v = alpha, col = 'orange')
abline(v = xmin, col = 'orange', lty = 2)


# ===================================== #
# ==== SIMULATE DISCRETE AGE DATA: ====
# ===================================== #
# Individual IDs:
indID <- sapply(1:N, function(n) {
  Nch <- as.character(n)
  nch <- nchar(Nch)
  NumChar <- paste(paste(rep(0, 4 - nch), collapse = ""), Nch, sep = "")
  return(NumChar)
})
indID <- sprintf("I%s", indID)

# Ages at death:
X <- SampleRandAge(n = N, theta = theta, dx = 1, model = "GO", 
                            shape = "bathtub")

# Create longevity data:
longdat <- data.frame(indID = indID, ageEntry = rep(0, N), 
                      ageDepart = X)

# Find maximum age:
omega <- max(X)

# Age vector:
xv <- 0:omega

# Age specific fertility:
gx <- xv * 0
gx[which(xv >= alpha)] <- CalcFert(beta = beta, 
                                   x = xv[which(xv >= alpha)] - alpha + 0.5, 
                                   modelFert = "gamma")


# Individual effects:
uini <- rgamma(n = N, shape = gamma, rate = gamma)

# Individuals alive after sexual maturity:
idSim <- which(X > alpha)

# Simulate births:
for (ii in idSim) {
  iages <- alpha:X[ii]
  idages <- which(xv %in% iages)
  nages <- length(idages)
  gxi <- gx[idages]
  noffsi <- rpois(n = nages, lambda = gxi * uini[ii])
  temp <- data.frame(indID = rep(indID[ii], nages), Age = iages, 
                     nOffspring = noffsi,
                     ageDeath = rep(X[ii], nages), stringsAsFactors = FALSE)
  if (ii == idSim[1]) {
    reprdat <- temp
  } else {
    reprdat <- rbind(reprdat, temp)
  }
}


# List of settings:
simSettings <- list(N = N, Nyears = Nyears, theta = theta, beta = beta, 
                    gamma = gamma, alpha = alpha, omega = omega, dx = dx, 
                    x = x)

# Save tables:
save(list = c("simSettings", "longdat", "reprdat"),
     file = sprintf("data/simDataSimp%sSps.RData", spType))
