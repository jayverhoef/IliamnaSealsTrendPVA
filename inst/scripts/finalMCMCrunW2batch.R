# can run this in command window, run this as R CMD BATCH '/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/inst/scripts/finalMCMCrunW2batch.R'
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


W2 <- MCMCdd(burnin = 50000, niter = 100000, thin = 10,
  adelta = 38.25,
  bdelta = 7.14,
  aphi = 3.53,
  bphi = 12.91,
  akappa = 4.11,
  bkappa = 1.51,
	mupsi = 0.5,
  xi = 10,
  maxSight = .9,
  harvestData = harvestData, data = d2)
save(W2, file = "/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/data/W2.rda")


#plot((1:99)/100, dbeta((1:99)/100, .5*50, (1 - .5)*50), type = 'l', ylim = c(0,15))
#hist(W2$psi, freq = FALSE, add = TRUE)


