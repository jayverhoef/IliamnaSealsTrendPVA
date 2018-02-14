library(IliamnaSealsTrendPVA)
data(Iliamna_counts_revised)
data(harvest)
d2 <- Iliamna_counts_revised

ind <- is.na(d2$puptotal)

d2[,"pupcount"] <- d2[,"puptotal"]
d2[,"aducount"] <- d2[,"adulttotal"]
d2[,"totcount"] <- d2[,"adulttotal"]
d2[!ind,"totcount"] <- NA
d2[ind,"aducount"] <- NA
d2[,"year"] <- as.POSIXlt(d2[,"datetime"])$year + 1900
d2[,"yday"] <- as.POSIXlt(d2[,"datetime"])$yday
d2[,"hour"] <- as.POSIXlt(d2[,"datetime"])$hour
d2[d2$hour == 0,"hour"] <- NA
d2[,"yr"] <- as.POSIXlt(d2[,"datetime"])$year - 83
d2[,"dy"] <- (d2$yday - mean(d2$yday))/sqrt(var(d2$yday))
d2[,"hr"] <- d2$hour
d2[!is.na(d2$hr),"hr"] <- (d2[!is.na(d2$hr),"hr"] - 
  mean(d2[!is.na(d2$hr),"hr"]))/sqrt(var(d2[!is.na(d2$hr),"hr"]))

W1 <- MCMCdd(burnin = 10000, niter = 100000, thin = 10,
  adelta = 272.63,
  bdelta = 33.52,
  aphi = 36.86,
  bphi = 125.43,
  akappa = 40.34,
  bkappa = 11.85,
  xi = 100,
  maxSight = .9,
  harvestData = harvestData, data = d2)
save(W1, file = "/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/data/W1.rda")


W2 <- MCMCdd(burnin = 10000, niter = 100000, thin = 10,
  adelta = 38.25,
  bdelta = 7.14,
  aphi = 3.53,
  bphi = 12.91,
  akappa = 4.11,
  bkappa = 1.51,
  xi = 100,
  maxSight = .9,
  harvestData = harvestData, data = d2)
save(W2, file = "/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/data/W2.rda")


W3 = MCMCddUniPrior(burnin = 1000, niter = 10000, thin = 10,
  delta.lower = 0.5,
  delta.upper = 1,
  phi.lower = 0,
  phi.upper = 0.5,
  kappa.lower = 0.5,
  kappa.upper = 1,
  xi = 100,
  maxSight = .9,
  harvestData = harvestData, data = d2)
save(W3, file = "/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/data/W3.rda")


W4 <- MCMCdd(burnin = 10000, niter = 100000, thin = 10,
  adelta = 38.25,
  bdelta = 7.14,
  aphi = 3.53,
  bphi = 12.91,
  akappa = 4.11,
  bkappa = 1.51,
  xi = 50,
  maxSight = .9,
  harvestData = harvestData, data = d2)
save(W4, file = "/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/data/W4.rda")

W5 <- MCMCdd(burnin = 10000, niter = 100000, thin = 10,
  adelta = 38.25,
  bdelta = 7.14,
  aphi = 3.53,
  bphi = 12.91,
  akappa = 4.11,
  bkappa = 1.51,
  xi = 10,
  maxSight = .9,
  harvestData = harvestData, data = d2)
save(W5, file = "/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/data/W5.rda")

W6 = MCMCddPsiUni(burnin = 10000, niter = 100000, thin = 10,
  adelta = 38.25,
  bdelta = 7.14,
  aphi = 3.53,
  bphi = 12.91,
  akappa = 4.11,
  bkappa = 1.51,
	psi.lower = .1,
	psi.upper = .9,
  maxSight = .9,
  harvestData = harvestData, data = d2)
save(W6, file = "/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/data/W6.rda")


save(W1, file = "/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/data/W1.rda")
save(W2, file = "/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/data/W2.rda")
save(W3, file = "/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/data/W3.rda")
save(W4, file = "/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/data/W4.rda")
save(W5, file = "/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/data/W5.rda")
save(W6, file = "/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/data/W6.rda")

mu = .5
phi = 100
a = mu*phi
b = (1 - mu)*phi

plot(dbeta((1:99)/100,a,b))
