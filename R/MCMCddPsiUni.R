#-------------------------------------------------------------------------------
#
#          MCMC
#
#-------------------------------------------------------------------------------

#' MCMC sampling
#'
#' MCMC sampling
#'
#' @param burnin the number of MCMC for burnin
#' @param niter the number of MCMC iterations
#' @param thin thinning of the MCMC chain
#' @param adelta beta parameter a for prior on delta
#' @param bdelta beta parameter b for prior on delta
#' @param aphi beta parameter a for prior on phi
#' @param bphi beta parameter b for prior on phi
#' @param akappa beta parameter a for prior on kappa
#' @param bkappa beta parameter b for prior on kappa
#' @param harvestData harvest dataset
#' @param data  survey dataset
#' @param psi.lower lower bound on uniform prior for psi
#' @param psi.upper upper bound on uniform prior for psi
#' @param maxSight the maximum sightability across all surveys
#'
#' @return a list of the MCMC chains values
#'
#' @author Jay Ver Hoef
#' @export

MCMCddPsiUni <- function(burnin = 10000, niter = 100000, thin = 100, 
  adelta, bdelta, aphi, bphi, akappa, bkappa, harvestData, data, 
	psi.lower, psi.upper, xi, maxSight)
{
  d2 = data
  #hard-wired starting values
  delta = runif(29)*.01 + .82
  phi = runif(29)*.01 + .29
  kappa1 = runif(29)*.01 + .74
  betaXdy = 1
  betaXdy2 = -2
  betaXhr = 1
  betaXhr2 = -2
  betaYdy = 1
  betaYdy2 = -2
  betaYhr = 1
  betaYhr2 = -2
  thetaX1 = 200
  thetaY1 = 500
  sigX = 20
  sigY = 90
  thetaY <- rep(NA, times = 30)
  thetaX <- thetaY
  thetaY[1] <- thetaY1
  thetaX[1] <- thetaX1
  psi = .5
	rho = 0.03
  harvDat <- harvestData[5:10,]
  indxObs <- (harvDat$Year - 1984)[(harvDat$Year - 1984) > 0]  
  H <- rep(0, times = 30)
  for(i in 2:30) {
    thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1] - H[i]
    thetaX[i] <- phi[i-1]*thetaY[i-1]
  }  
  nMissHr <- sum(is.na(d2$hr))
  hr <- d2$hr
  hr[is.na(hr)] <- sample(d2[!is.na(d2$hr),"hr"],nMissHr, replace = TRUE)


  testLamMu <- function(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2) 
  {
    thetaY <- rep(NA, times = 30)
    thetaX <- thetaY
    thetaY[1] <- thetaY1
    thetaX[1] <- thetaX1
    for(i in 2:30) {
      thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1] - H[i]
      thetaX[i] <- exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1]
    }  
    lamX <- thetaX[d2$yr] + betaXdy*d2$dy + betaXdy2*d2$dy^2 + 
      betaXhr*hr + betaXhr2*hr^2
    lamY <- thetaY[d2$yr] + betaYdy*d2$dy + betaYdy2*d2$dy^2 + 
      betaYhr*hr + betaYhr2*hr^2
    !(any(
         (lamX[!is.na(d2$pupcount)] < 
            d2[!is.na(d2$pupcount),"pupcount"]/maxSight) |
         (lamY[!is.na(d2$aducount)] < 
            d2[!is.na(d2$aducount),"aducount"]/maxSight)
       ) |
      any((lamX[is.na(d2$pupcount)] +  lamY[is.na(d2$aducount)]) < 
          d2[!is.na(d2$totcount),"totcount"]/maxSight
        )
      )
  }

harvSample <- function(H) {
    H = harvestData
    H$year = H$year - 1984
    H.Ig.Obs = H[H$community == 'Igiugig',]
    H.Ig = rep(0, times = 30)
    H.Ig[H.Ig.Obs$year] = H.Ig.Obs$harvest
    H.Ig[H.Ig.Obs$year]
    H.Ig[!((1:30) %in% H.Ig.Obs$year)] = sample(H.Ig.Obs$harvest, 
      30 - length(H.Ig.Obs[,1]), replace = TRUE)
    H.Il.Obs = H[H$community == 'Iliamna',]
    H.Il = rep(0, times = 30)
    H.Il[H.Il.Obs$year] = H.Il.Obs$harvest
    H.Il[H.Il.Obs$year]
    H.Il[!((1:30) %in% H.Il.Obs$year)] = sample(H.Il.Obs$harvest, 
      30 - length(H.Il.Obs[,1]), replace = TRUE)
    H.Ko.Obs = H[H$community == 'Kokhanok',]
    H.Ko = rep(0, times = 30)
    H.Ko[H.Ko.Obs$year] = H.Ko.Obs$harvest
    H.Ko[H.Ko.Obs$year]
    H.Ko[!((1:30) %in% H.Ko.Obs$year)] = sample(H.Ko.Obs$harvest, 
      30 - length(H.Ko.Obs[,1]), replace = TRUE)
    H.Le.Obs = H[H$community == 'Levelock',]
    H.Le = rep(0, times = 30)
    H.Le[H.Le.Obs$year] = H.Le.Obs$harvest
    H.Le[H.Le.Obs$year]
    H.Le[!((1:30) %in% H.Le.Obs$year)] = sample(H.Le.Obs$harvest, 
      30 - length(H.Le.Obs[,1]), replace = TRUE)
    H.Ne.Obs = H[H$community == 'Newhalen',]
    H.Ne = rep(0, times = 30)
    H.Ne[H.Ne.Obs$year] = H.Ne.Obs$harvest
    H.Ne[H.Ne.Obs$year]
    H.Ne[!((1:30) %in% H.Ne.Obs$year)] = sample(H.Ne.Obs$harvest, 
      30 - length(H.Ne.Obs[,1]), replace = TRUE)
    H.Pe.Obs = H[H$community == 'Pedro Bay',]
    H.Pe = rep(0, times = 30)
    H.Pe[H.Pe.Obs$year] = H.Pe.Obs$harvest
    H.Pe[H.Pe.Obs$year]
    H.Pe[!((1:30) %in% H.Pe.Obs$year)] = sample(H.Pe.Obs$harvest, 
      30 - length(H.Pe.Obs[,1]), replace = TRUE)
    H.Ig + H.Il + H.Ko + H.Le + H.Ne + H.Pe
  }


    for(k in 1:burnin) {

    H.try <- harvSample(harvestData)

    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H.try,
       betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
      H <- H.try
    }  
    hr.old = hr  
    hr <- d2$hr
    hr[is.na(hr)] <- sample(d2[!is.na(d2$hr),"hr"], 
      nMissHr, replace = TRUE)
    if(!testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H.try,
       betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
      hr <- hr.old
    }
    U <- log(runif(1))
    thetaX1.try <- max(rnorm(1, thetaX1, 2),0)
    if(testLamMu(thetaX1.try, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1.try, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) thetaX1 <- thetaX1.try
    }
    U <- log(runif(1))
    thetaY1.try <- max(rnorm(1, thetaY1, 2),0)
    if(testLamMu(thetaX1, thetaY1.try, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1.try, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) thetaY1 <- thetaY1.try
    }
    U <- log(runif(1))
    betaXdy.try <- rnorm(1, betaXdy, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy.try, H,
      betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
     LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy.try, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaXdy <- betaXdy.try
    }
    U <- log(runif(1))
    betaXdy2.try <- rnorm(1, betaXdy2, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2.try, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
     LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2.try, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaXdy2 <- betaXdy2.try
    }
    U <- log(runif(1))
    betaXhr.try <- rnorm(1, betaXhr, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr.try, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr.try, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaXhr <- betaXhr.try
    }
    U <- log(runif(1))
    betaXhr2.try <- rnorm(1, betaXhr2, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr2.try, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2.try,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaXhr2 <- betaXhr2.try
    }
    U <- log(runif(1))
    betaYdy.try <- rnorm(1, betaYdy, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr, betaYdy.try, betaYdy2, betaYhr, betaYhr2, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy.try, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaYdy <- betaYdy.try
    }
    U <- log(runif(1))
    betaYdy2.try <- rnorm(1, betaYdy2, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2.try, betaYhr, betaYhr2, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2.try, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaYdy2 <- betaYdy2.try
    }
    U <- log(runif(1))
    betaYhr.try <- rnorm(1, betaYhr, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr.try, betaYhr2, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr.try, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaYhr <- betaYhr.try
    }
    U <- log(runif(1))
    betaYhr2.try <- rnorm(1, betaYhr2, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2.try, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2.try, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaYhr2 <- betaYhr2.try
    }
    for(i in 1:length(delta)) {
      delta.try <- delta
      U <- log(runif(1))
      delta.try[i] <- min(max(delta[i] + runif(1)*0.04 - 0.02, 0),1)
      if(testLamMu(thetaX1, thetaY1, delta.try, kappa1, phi, rho, betaXdy, H,
        betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
      LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta.try, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
        LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
          betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
          adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
      if(LLdif > U) delta <- delta.try
      }
    }
    for(i in 1:length(phi)) {
      phi.try <- phi
      U <- log(runif(1))
      phi.try[i] <- min(max(phi[i] + runif(1)*0.04 - 0.02, 0),1)
      if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi.try, rho, betaXdy, H,
        betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
      LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi.try, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
        LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
          betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
          adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
      if(LLdif > U) phi <- phi.try
      }
    }
    for(i in 1:length(kappa1)) {
      kappa1.try <- kappa1
      U <- log(runif(1))
      kappa1.try[i] <- min(max(kappa1[i] + runif(1)*0.04 - 0.02, 0),1)
      if(testLamMu(thetaX1, thetaY1, delta, kappa1.try, phi, rho, betaXdy, H,
        betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
      LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1.try, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
        LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
          betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
          adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
      if(LLdif > U) kappa1 <- kappa1.try
      }
    }    
    U <- log(runif(1))
    sigX.try <- min(max(sigX + runif(1)*4 - 2, 0),30)
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX.try, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) sigX <- sigX.try
    U <- log(runif(1))
    sigY.try <- min(max(sigY + runif(1)*4 - 2, 0),100)
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY.try, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) sigY <- sigY.try
    U <- log(runif(1))
    psi.try <- min(max(psi + runif(1)*0.04 - 0.02, psi.lower),psi.upper)
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi.try, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) psi <- psi.try
		U <- log(runif(1))
    rho.try <- min(max(rho + runif(1)*0.06 - 0.03, 0),1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho.try, betaXdy, H,
        betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
		LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho.try) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) rho <- rho.try
    }
    if(k%%thin == 0) {
      cat("\r", "Burnin Number: ", k)
    }
  }
  W <- vector("list", 18)
  for(k in 1:niter) {

    H.try <- harvSample(harvestData)

    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H.try,
       betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
      H <- H.try
    }  
    hr.old = hr  
    hr <- d2$hr
    hr[is.na(hr)] <- sample(d2[!is.na(d2$hr),"hr"], 
      nMissHr, replace = TRUE)
    if(!testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H.try,
       betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
      hr <- hr.old
    }
    U <- log(runif(1))
    thetaX1.try <- max(rnorm(1, thetaX1, 2),0)
    if(testLamMu(thetaX1.try, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1.try, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) thetaX1 <- thetaX1.try
    }
    U <- log(runif(1))
    thetaY1.try <- max(rnorm(1, thetaY1, 2),0)
    if(testLamMu(thetaX1, thetaY1.try, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1.try, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) thetaY1 <- thetaY1.try
    }
    U <- log(runif(1))
    betaXdy.try <- rnorm(1, betaXdy, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy.try, H,
      betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
     LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy.try, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaXdy <- betaXdy.try
    }
    U <- log(runif(1))
    betaXdy2.try <- rnorm(1, betaXdy2, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2.try, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
     LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2.try, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaXdy2 <- betaXdy2.try
    }
    U <- log(runif(1))
    betaXhr.try <- rnorm(1, betaXhr, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr.try, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr.try, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaXhr <- betaXhr.try
    }
    U <- log(runif(1))
    betaXhr2.try <- rnorm(1, betaXhr2, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr2.try, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2.try,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaXhr2 <- betaXhr2.try
    }
    U <- log(runif(1))
    betaYdy.try <- rnorm(1, betaYdy, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr, betaYdy.try, betaYdy2, betaYhr, betaYhr2, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy.try, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaYdy <- betaYdy.try
    }
    U <- log(runif(1))
    betaYdy2.try <- rnorm(1, betaYdy2, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2.try, betaYhr, betaYhr2, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2.try, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaYdy2 <- betaYdy2.try
    }
    U <- log(runif(1))
    betaYhr.try <- rnorm(1, betaYhr, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr.try, betaYhr2, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr.try, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaYhr <- betaYhr.try
    }
    U <- log(runif(1))
    betaYhr2.try <- rnorm(1, betaYhr2, .1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho, betaXdy, H,
      betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2.try, d2)) {
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2.try, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) betaYhr2 <- betaYhr2.try
    }
    for(i in 1:length(delta)) {
      delta.try <- delta
      U <- log(runif(1))
      delta.try[i] <- min(max(delta[i] + runif(1)*0.04 - 0.02, 0),1)
      if(testLamMu(thetaX1, thetaY1, delta.try, kappa1, phi, rho, betaXdy, H,
        betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
      LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta.try, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
        LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
          betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
          adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
      if(LLdif > U) delta <- delta.try
      }
    }
    for(i in 1:length(phi)) {
      phi.try <- phi
      U <- log(runif(1))
      phi.try[i] <- min(max(phi[i] + runif(1)*0.04 - 0.02, 0),1)
      if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi.try, rho, betaXdy, H,
        betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
      LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi.try, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
        LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
          betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
          adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
      if(LLdif > U) phi <- phi.try
      }
    }
    for(i in 1:length(kappa1)) {
      kappa1.try <- kappa1
      U <- log(runif(1))
      kappa1.try[i] <- min(max(kappa1[i] + runif(1)*0.04 - 0.02, 0),1)
      if(testLamMu(thetaX1, thetaY1, delta, kappa1.try, phi, rho, betaXdy, H,
        betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
      LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1.try, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
        LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
          betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
          adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
      if(LLdif > U) kappa1 <- kappa1.try
      }
    }    
    U <- log(runif(1))
    sigX.try <- min(max(sigX + runif(1)*4 - 2, 0),30)
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX.try, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) sigX <- sigX.try
    U <- log(runif(1))
    sigY.try <- min(max(sigY + runif(1)*4 - 2, 0),100)
    LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY.try, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) sigY <- sigY.try
    U <- log(runif(1))
    psi.try <- min(max(psi + runif(1)*0.04 - 0.02, psi.lower),psi.upper)
		LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi.try, xi, rho) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) psi <- psi.try
		U <- log(runif(1))
    rho.try <- min(max(rho + runif(1)*0.06 - 0.03, 0),1)
    if(testLamMu(thetaX1, thetaY1, delta, kappa1, phi, rho.try, betaXdy, H,
        betaXdy2, betaXhr, betaXhr2, betaYdy, betaYdy2, betaYhr, betaYhr2, d2)) {
		LLdif <- LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho.try) -
      LLddPsiUni(hr, thetaX1, thetaY1, betaXdy, betaXdy2, betaXhr, betaXhr2,
        betaYdy, betaYdy2, betaYhr, betaYhr2, delta, phi, kappa1, sigX, sigY, d2, H,
        adelta, bdelta, aphi, bphi, akappa, bkappa, psi, xi, rho)
    if(LLdif > U) rho <- rho.try
		}
    if(k%%thin == 0) {
      cat("\r", "Iteration Number: ", k)
      W[[1]] <- c(W[[1]],thetaX1)
      W[[2]] <- c(W[[2]],thetaY1)
      W[[3]] <- c(W[[3]],betaXdy)
      W[[4]] <- c(W[[4]],betaXdy2)
      W[[5]] <- c(W[[5]],betaXhr)
      W[[6]] <- c(W[[6]],betaXhr2)
      W[[7]] <- c(W[[7]],betaYdy)
      W[[8]] <- c(W[[8]],betaYdy2)
      W[[9]] <- c(W[[9]],betaYhr)
      W[[10]] <- c(W[[10]],betaYhr2)
      W[[11]] <- cbind(W[[11]],delta)
      W[[12]] <- cbind(W[[12]],phi)
      W[[13]] <- cbind(W[[13]],kappa1)
      W[[14]] <- c(W[[14]],sigX)
      W[[15]] <- c(W[[15]],sigY)
      W[[16]] <- cbind(W[[16]],H)
      W[[17]] <- cbind(W[[17]],psi)
			W[[18]] <- c(W[[18]],rho)
    }
  }
  cat("\n")
  names(W) <- c("thetaX1", "thetaY1", "betaXdy", "betaXdy2", "betaXhr", 
    "betaXhr2", "betaYdy", "betaYdy2", "betaYhr", "betaYhr2", "delta", 
    "phi", "kappa", "sigX", "sigY", "H", "psi", "rho")
  W
}

