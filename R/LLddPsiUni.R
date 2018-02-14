#-------------------------------------------------------------------------------
#
#          LL
#
#-------------------------------------------------------------------------------

#' loglikelihood for model
#'
#' evaluates to the loglikelihood for MCMC sampling
#'
#' @param hr hour of day
#' @param thetaX1 the starting value for the X-variable, which is the underlying pup population
#' @param thetaY1 the starting value for the Y-variable, which is the underlying nonpup population
#' @param betaXdy day-of-year linear regression coefficient for pups
#' @param betaXdy2 day-of-year quadratic regression coefficient for pups
#' @param betaXhr hour-of-day linear regression coefficient for pups
#' @param betaXhr2 hour-of-day quadratic regression coefficient for pups
#' @param betaYdy day-of-year linear regression coefficient for nonpups
#' @param betaYdy2 day-of-year quadratic regression coefficient for nonpups
#' @param betaYhr hour-of-day linear regression coefficient for nonpups
#' @param betaYhr2 hour-of-day quadratic regression coefficient for nonpups
#' @param delta annual survivorship of nonpups 
#' @param phi annual fecundity of nonpups to time of survey
#' @param kappa1 annual survivorship of pups from time to survey to yearling
#' @param sigX measurement error for pups
#' @param sigY measurement error for nonpups
#' @param d2 dataset
#' @param H harvest dataset
#' @param adelta beta parameter a for prior on delta
#' @param bdelta beta parameter b for prior on delta
#' @param aphi beta parameter a for prior on phi
#' @param bphi beta parameter b for prior on phi
#' @param akappa beta parameter a for prior on kappa
#' @param bkappa beta parameter b for prior on kappa
#' @param psi proportion of animals that are out of the water
#' @param rho density dependent fecundity parameter
#'
#' @return the loglikelihood
#'
#' @author Jay Ver Hoef

LLddPsiUni <- function(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
  betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
  adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
{
  thetaY <- rep(NA, times = 30)
  thetaX <- thetaY
  thetaY[1] <- thetaY1
  thetaX[1] <- thetaX1
  for(i in 2:30) {
    thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1] - H[i]
    thetaX[i] <- exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1]
  }  
  lamX <- psi*thetaX[d2$yr] + betaXdy*d2$dy + betaXdy2*d2$dy^2 + 
    betaXhr*hr + betaXhr2*hr^2
  lamY <- psi*thetaY[d2$yr] + betaYdy*d2$dy + betaYdy2*d2$dy^2 + 
    betaYhr*hr + betaYhr2*hr^2
  sum(dnorm(d2[!is.na(d2$pupcount),"pupcount"], lamX[!is.na(d2$pupcount)],
      sigX, log = TRUE)) +
    sum(dnorm(d2[!is.na(d2$aducount),"aducount"], lamY[!is.na(d2$aducount)],
      sigY, log = TRUE)) +
    sum(dnorm(d2[!is.na(d2$totcount),"totcount"], lamX[is.na(d2$pupcount)] + 
      lamY[is.na(d2$aducount)], sqrt(sigX^2 + sigY^2), log = TRUE)) + 
    sum(dbeta(delta,adelta,bdelta,log = TRUE)) + sum(dbeta(phi,aphi,bphi, log = TRUE)) +
    sum(dbeta(kappa1,akappa,bkappa, log = TRUE))
}
