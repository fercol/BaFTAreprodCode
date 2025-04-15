# ============================== CODE METADATA =============================== #
# AUTHOR: Fernando Colchero
# DATE CREATED: 2024-12-01
# DESCRIPTION: Addtional functions for testing BaFTA on simulated data.
# NOTES: 
# ================================ CODE START ================================ #
# Truncated normal:
rtnorm <- function(n, mean, sd, lower = -Inf, upper = Inf) {
  Flow <- pnorm(lower, mean, sd)
  Fup <- pnorm(upper, mean, sd)
  ru <- runif(n, Flow, Fup)
  rx <- qnorm(ru, mean, sd)
  return(rx)
}

dtnorm <- function(x, mean, sd, lower = -Inf, upper = Inf, log = FALSE) {
  Flow <- pnorm(lower, mean, sd)
  Fup <- pnorm(upper, mean, sd)
  densx <- dnorm(x, mean, sd) / (Fup - Flow)
  if (log) densx <- log(densx)
  return(densx)
}

ptnorm <- function(q, mean, sd, lower = -Inf, upper = Inf, log = FALSE) {
  p <- (pnorm(q, mean, sd) - pnorm(lower, mean, sd)) / 
    (pnorm(upper, mean, sd) - pnorm(lower, mean, sd))
  if (log) {
    p <- log(p)
  }
  return(p)
}

qtnorm <- function (p, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  p2 <- (p) * (pnorm(upper, mean, sd) - pnorm(lower, mean, sd)) + 
    pnorm(lower, mean, sd)
  q <- qnorm(p2, mean, sd)
  return(q)
}

# Empty plot:
plot0 <- function(x = c(0, 1), y = c(0, 1), ...) {
  plot(x, y, col = NA, xlab = "", ylab = "", axes = FALSE)
}

PlotDensTrace <- function(out, plotDir, figName, saveResults = FALSE,
                          simPars = NULL) {
  # Number of MCMC chains:
  nsim <- out$settings$nsim
  
  # Number of parameters:
  p <- nrow(out$coefficients)
  
  # Parameter names:
  pnames <- rownames(out$coefficients)
  
  # Fancy names:
  labpnames <- c(b0 = expression(theta[0]), b1 = expression(theta[1]), 
                 b2 = expression(theta[2]), b3 = expression(theta[3]),
                 b4 = expression(theta[4]), b5 = expression(theta[5]), 
                 b6 = expression(theta[6]), gamma = expression(gamma),
                 eta = expression(eta), kappa = expression(kappa), 
                 vSd = expression(italic(sigma[v])))
  labpnames <- labpnames[which(names(labpnames) %in% pnames)]
  
  # Converged parameters:
  idSamp <- out$params$idSamp
  if (out$settings$dataType == "indivExtended") {
    idSamp <- c(idSamp, max(idSamp) + 1)
  }
  
  # Converged parameter matrix:
  parsMat <- out$theta[, pnames]
  
  # data for trace plots:
  keep <- seq(1, out$settings$niter, out$settings$thinning)
  nkeep <- length(keep)
  fullPars <- array(0, dim = c(nkeep, p, nsim))
  for (ii in 1:nsim) {
    fullPars[,, ii] <- out$runs[[ii]]$theta[keep, pnames]
  }
  
  # Y limits for traces:
  ylimtr <- sapply(1:p, function(ip) {
    range(fullPars[, ip, ])
  })
  dimnames(ylimtr) <- list(c("low", "upp"), pnames)
  
  # x limits for densities:
  xlimden <- apply(parsMat, 2, function(pv) {
    quantile(pv, c(0.001, 0.999))
  })
  dimnames(xlimden) <- list(c("low", "upp"), pnames)
  
  # Prepare layout:
  ppmat <- matrix(1:(p * 2) + 4, p, 2)
  ppmat <- cbind(ppmat[, 1], rep(2, p), ppmat[, 2])
  laymat <- cbind(c(rep(1, p), 0), rbind(ppmat, c(3, 0, 4)))
  widths <- c(0.2, 1, 0.2, 1)
  heights <- c(rep(0.5, p), 0.2)
  
  # PDF settings:
  whr <- sum(widths) / sum(heights)
  pdfw <- 6
  pdfh <- pdfw / whr
  
  # Margins:
  mar <- c(2, 2, 1, 1)
  
  # Trace colors:
  coltr <- brewer.pal(9, "Set1")[-6][1:nsim]
  
  # Produce plot:
  if (saveResults) {
    pdf(file = sprintf("%s%s.pdf", plotDir, figName), width = pdfw, 
        height = pdfh)
  }
  
  layout(mat = laymat, widths = widths, heights = heights)
  
  # Y label 1:
  par(mar = mar * c(1, 0, 1, 0))
  plot0()
  text(0.5, 0.5, "Parameter posterior density", srt = 90, cex = 2)
  
  # Y label 2:
  par(mar = mar * c(1, 0, 1, 0))
  plot0()
  text(0.5, 0.5, "Parameter traces", srt = 90, cex = 2)
  
  # X label 1:
  par(mar = mar * c(0, 1, 0, 1))
  plot0()
  text(0.5, 0.5, "Parameter value", cex = 2)
  
  # Y label 2:
  par(mar = mar * c(0, 1, 0, 1))
  plot0()
  text(0.5, 0.5, "MCMC step", cex = 2)
  
  # Density plots
  for (ip in 1:p) {
    # Posterior density:
    thev <- parsMat[, ip]
    idth <- which(thev >= xlimden[1, ip] & thev <= xlimden[2, ip])
    thev <- thev[idth]
    thest <- out$coefficients[ip, c(1:4)]
    dthe <- density(thev)
    ylim <- c(0, max(dthe$y))
    idin <- which(dthe$x >= thest["Lower"] & dthe$x <= thest["Upper"])
    idmean <- which(abs(dthe$x - thest["Mean"]) == 
                      min(abs(dthe$x - thest["Mean"])))[1]
    par(mar = mar, las = 1)
    plot(xlimden[, ip], ylim, col = NA, lwd = 2, xlab = "", ylab = "")
    polygon(dthe$x[c(idin, rev(idin))], c(dthe$y[idin], rep(0, length(idin))),
            col = "orange", border = NA)
    lines(dthe$x[rep(idmean, 2)], c(0, dthe$y[idmean]), col = 'white', lwd = 2)
    lines(dthe$x, dthe$y, col = 'orange', lwd = 2)
    if (!is.null(simPars)) {
      lines(x = rep(simPars[ip], 2), y = ylim, col = 'dark red', 
            lwd = 2)
    }
    text(xlimden[1, ip] + diff(xlimden[, ip]) * 0.025, 
         ylim[2] - diff(ylim) * 0.2,
         labpnames[ip], cex = 2, adj = 0)
  }
  
  # Trace plots:
  xlim <- range(keep)
  for (ip in 1:p) {
    par(mar = mar, las = 1)
    plot(xlim, ylimtr[, ip], col = NA, xlab = "", ylab = "")
    for (isim in 1:nsim) {
      lines(keep, fullPars[, ip, isim], col = coltr[isim])
    }
  }
  if (saveResults) dev.off()
  
}
