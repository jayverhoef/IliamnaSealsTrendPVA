# can run this in command window, run this as R CMD BATCH '/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/inst/scripts/finalMCMCrunW3batch.R'
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


W3 = MCMCddU(burnin = 50000, niter = 100000, thin = 10,
  delta.lower = 0.5,
  delta.upper = 1,
  phi.lower = 0,
  phi.upper = 0.5,
  kappa.lower = 0.5,
  kappa.upper = 1,  
  mupsi = .5,
  xi = 10,
  maxSight = .9,
  harvestData = harvestData, data = d2)
save(W3, file = "/mnt/Hitachi2GB/00NMML/activePapers/IliamnaSealsTrendPVA/package/IliamnaSealsTrendPVA/data/W3.rda")

#plot((1:99)/100, (dbeta((1:99)/100, .5*10, (1 - .5)*10))^2, type = 'l', ylim = c(0,40))
#hist(W3$psi, freq = FALSE, add = TRUE)

