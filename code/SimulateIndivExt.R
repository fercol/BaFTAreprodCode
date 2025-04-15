# ============================== CODE METADATA =============================== #
# AUTHOR: Fernando Colchero
# DATE CREATED: 2024-12-19
# DESCRIPTION: Code to simulate different levels of age-specific fertility 
#              datasets to test with package.
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
spType <- "Fast"

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
  theta <- c(a0 = 0, a1 = 1, c = 0.5, b0 = -2, b1 = 0.1)
  
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
  alpha <- 0.25
  
  # Gestation and weaning period:
  gest <- 0.06
  wean <- 0.07
  tau <- gest + wean
}

# Vector of ages:
dx <- 0.001
x <- seq(0, 100, dx)
z <- seq(0, 5, dx)
w <- seq(0, 5, dx)

# ============================ #
# ==== VISUAL INSPECTION: ====
# ============================ #
# Calculate mortality and survival:
Sx <- CalcSurv(theta = theta, x = x, model = "GO", shape = "bathtub")
idsx <- which(Sx >= 1 / N)
xv <- x[idsx]
Sx <- Sx[idsx]
mux <- CalcMort(theta = theta, x = xv, model = "GO", shape = "bathtub")

# Calculate age-specific fertility:
gx <- CalcFert(beta = beta, x = xv[which(xv >= alpha)] - alpha, 
               modelFert = "gamma")

# Calculate interbirth interval (after gestation and weaning):
fz <- dexp(z, rate = eta)

# Calculate time since age at maturity to first birth:
fw <- dexp(w, rate = kappa)

# Age at minimum mortality:
xmin <- xv[which(mux == min(mux))]

# Plot vital rates:
par(mfrow = c(3, 2), mar = c(4, 4, 1, 1))
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

plot(w + alpha, fw, type = 'l', col = 'dark red', lwd = 4, 
     xlim = c(0, max(w) + alpha), xlab = "Age at first birth", 
     ylab = "PDF age first birth")
abline(v = alpha, col = 'orange')
abline(v = xmin, col = 'orange', lty = 2)

plot(z + tau, fz, type = 'l', col = 'dark red', lwd = 4, 
     xlim = c(0, max(z) + tau), xlab = "Interbirth interval", ylab = "PDF IBI")
abline(v = gest, col = 'red')
abline(v = tau, col = 'red')

# ---------------------------- #
# Simulate single life course:
# ---------------------------- #
# Individual IDs:
indID <- sapply(1:N, function(n) {
  Nch <- as.character(n)
  nch <- nchar(Nch)
  NumChar <- paste(paste(rep(0, 4 - nch), collapse = ""), Nch, sep = "")
  return(NumChar)
})
indID <- sprintf("I%s", indID)

# index of ages after sexual maturity:
idgx <- which(xv >= alpha)

# Age specific fertility:
gx <- xv * 0
gx[idgx] <- CalcFert(beta = beta, x = xv[idgx] - alpha, 
                                   modelFert = "gamma")

# Ages at death:
X <- SampleRandAge(n = N, theta = theta, model = "GO", 
                      shape = "bathtub")

# Birth date:
B <- rep(0, N)

# Death date:
D <- B + X

# Individual effects:
u <- rgamma(n = N, shape = gamma, rate = gamma)
v <- rnorm(n = N, mean = 0, sd = vSd)

# Longevity table:
longdat <- data.frame(indID = indID, birthDate = B, entryDate = 0, 
                      departDate = D, departType = rep("D", N), 
                      ageEntry = rep(0, N), ageDeath = X, 
                      parent = rep(NA, N), u = u, v = v)

# Start individual counter:
icount <- 0

# Logical to continue with loop:
Start <- Sys.time()
for (i in 1:N) {
  # Extract dates and ages for ind. i:
  Ai <- longdat$ageEntry[i]
  if (Ai > alpha) {
    ti <- longdat$entryDate[i]
    xit <- Ai
  } else {
    ti <- longdat$entryDate[i] + alpha - Ai
    xit <- alpha
  }
  Bi <- longdat$birthDate[i]
  Di <- longdat$departDate[i]
  Xi <- longdat$ageDeath[i]
  ui <- longdat$u[i]
  vi <- longdat$v[i]
  ibi <- 0
  
  if (ti <= Tmax) {
    # Start reproductive events counter:
    irep <- 0
    
    while (xit < Xi & ti <= Tmax) {
      # Age at time ti:
      if (xit == alpha) {
        wi <- rexp(n = 1, rate = kappa)
        xit <- xit + wi
        ti <- ti + wi
        r1i <- 1
      } else {
        ibi <- tau + rexp(n = 1, rate = eta * exp(vi))
        xit <- xit + ibi
        ti <- ti + ibi
        r1i <- 0
      }
      
      if (xit < Xi & ti <= Tmax) {
        # Update ind. counter:
        icount <- icount + 1
        
        # Birth event:
        idgxi <- findInterval(xit, xv)
        yit <- rpois(n = 1, lambda = gx[idgxi] * ui)
        reprtemp <- data.frame(indID = indID[i], Age = xit, nOffspring = yit, 
                               IBI = ibi, First = r1i)
        irep <- irep + 1
        if (icount == 1 & irep == 1) {
          reprdat <- reprtemp
        } else {
          reprdat <- rbind(reprdat, reprtemp)
        }
      }
    }
  }
}
End <- Sys.time()

print(End - Start)

# Find maximum age:
omega <- ceiling(max(max(reprdat$Age)))

# List of settings:
simSettings <- list(N = N, Nyears = Nyears, theta = theta, beta = beta, 
                    gamma = gamma, eta = eta, kappa = kappa, vSd = vSd, 
                    alpha = alpha, omega = omega, gest = gest, wean = wean, 
                    tau = tau, dx = dx, x = x, z = z, w = w)

# Save tables:
save(list = c("simSettings", "longdat", "reprdat"),
     file = sprintf("data/simDataExt%sSps.RData", spType))

