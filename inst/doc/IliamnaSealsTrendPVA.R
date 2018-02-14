## ----setup, include=FALSE------------------------------------------------
dig <- function(x,d) {formatC(x, format="f", digits=d)}

## ----sum-table, echo=FALSE, include = FALSE, message=FALSE---------------
library(IliamnaSealsTrendPVA)
library(dplyr)
library(lubridate)

data(Iliamna_counts_revised)
counts_df <- tbl_df(Iliamna_counts_revised)

counts_df$datetime <- ymd_hms(counts_df$datetime,tz="US/Alaska")

by_agency_year <- group_by(counts_df,data_source = agency,year = as.factor(year(datetime)))
by_agency_year <- arrange(by_agency_year,agency,year(datetime))
sum_tbl <- summarize(by_agency_year,
          first_survey = as.character(min(as.Date(datetime))),
          last_survey = as.character(max(as.Date(datetime))),
          n_surveys = n())
sum_tbl = as.data.frame(sum_tbl)
names(sum_tbl) <- c("Data Source","Year","First Survey","Last Survey","No. of Surveys")
sum_tbl[,1] = as.character(sum_tbl[,1])
sum_tbl[sum_tbl['Data Source'] == 'adfg','Data Source'] = 'ADFG'
sum_tbl[sum_tbl['Data Source'] == 'mathisen_kline','Data Source'] = 'MK92'
sum_tbl[sum_tbl['Data Source'] == 'newhalen','Data Source'] = 'Newhalen'
sum_tbl[sum_tbl['Data Source'] == 'nmml','Data Source'] = 'NMML'
sum_tbl[sum_tbl['Data Source'] == 'pebble_abr','Data Source'] = 'ABR/Pebble'

sum_tbl = sum_tbl[order(sum_tbl$Year,sum_tbl[,"Data Source"]),]
library(xtable)
survTable = xtable(sum_tbl)

## ----echo = FALSE, include = FALSE---------------------------------------
  data(harvest)
  d2 <- Iliamna_counts_revised

  ind <- is.na(d2$puptotal)

  d2[,"pupcount"] <- d2[,"puptotal"]
  d2[,"aducount"] <- d2[,"adulttotal"]
  d2[,"totcount"] <- d2[,"adulttotal"]
  d2[!ind,"totcount"] <- NA
  d2[ind,"aducount"] <- NA
  #cbind(d2$pupcount, d2$aducount, d2$totcount)
  d2[,"year"] <- as.POSIXlt(d2[,"datetime"])$year + 1900
  d2[,"yday"] <- as.POSIXlt(d2[,"datetime"])$yday
  d2[,"hour"] <- as.POSIXlt(d2[,"datetime"])$hour
  d2[d2$hour == 0,"hour"] <- NA
  d2[,"yr"] <- as.POSIXlt(d2[,"datetime"])$year - 83
  d2[,"dy"] <- (d2$yday - mean(d2$yday))/sqrt(var(d2$yday))
  d2[,"hr"] <- d2$hour
  d2[!is.na(d2$hr),"hr"] <- (d2[!is.na(d2$hr),"hr"] - 
    mean(d2[!is.na(d2$hr),"hr"]))/sqrt(var(d2[!is.na(d2$hr),"hr"]))

  harvDat <- harvestData[5:10,]
  indxObs <- (harvDat$Year - 1984)[(harvDat$Year - 1984) > 0]
  H <- rep(NA, times = 30)
  H[indxObs] <- harvDat[,2]
   
  data(W1)
  data(W2)
  data(W3)
  data(W4)
  data(W5)
	data(W6)

  LesMat = matrix(c(0, .224, .773, .891), ncol = 2, byrow = TRUE)
  alpha4LesMat = .948


## ----CovEffects, echo = FALSE, include = FALSE, cache = TRUE-------------
  W = W2
  # A Pups and Day-of-Year
  ydaylim <- c(-2:2)*sqrt(var(d2$yday)) + mean(d2$yday)
  fits <- NULL
  for(k in 1:10000) {
     fit <- W[[3]][k]*(-20:20)/10 + W[[4]][k]*((-20:20)/10)^2
     fits <- cbind(fits, fit)
  }
  meanPupFit <- apply(fits,1,mean)
  ind = as.POSIXlt(seq(as.Date("2000/6/1"), by = "day", length.out = 120))$yday == round(
    ((-20:20)/10)[meanPupFit == max(meanPupFit)]*sqrt(var(d2$yday))+mean(d2$yday))
  optDayPups <- as.POSIXlt(seq(as.Date("2000/6/1"), by = "day", length.out = 120))[ind]
  optDayPupsC <- format(ISOdate(year = 2000, month = optDayPups$mon + 1, 
    day = optDayPups$mday), "%d %B")
  hrlim <- (-2:2)*sqrt(var(d2[!is.na(d2$hour),"hour"])) + mean(d2[!is.na(d2$hour),"hour"])
  hrlim <- round(hrlim*10)/10
  fits <- NULL
  for(k in 1:10000) {
    fit <- W[[5]][k]*(-20:20)/10 + W[[6]][k]*((-20:20)/10)^2
    fits <- cbind(fits, fit)
  }
  meanPupHourFit = apply(fits,1,mean)
  optHourPups = round(((-20:20)/10)[meanPupHourFit == max(meanPupHourFit)]*
    sqrt(var(d2[!is.na(d2$hour),"hour"])) + mean(d2[!is.na(d2$hour),"hour"]))*100
    ydaylim <- c(-2:2)*sqrt(var(d2$yday))+mean(d2$yday)
  fits <- NULL
  for(k in 1:10000) {
     fit <- W[[7]][k]*(-200:200)/100 + W[[8]][k]*((-200:200)/100)^2
     fits <- cbind(fits, fit)
  }
  meanNonPupFit <- apply(fits,1,mean)
  ind = as.POSIXlt(seq(as.Date("2000/6/1"), by = "day", length.out = 120))$yday == round(
    ((-200:200)/100)[meanNonPupFit == max(meanNonPupFit)]*sqrt(var(d2$yday))+mean(d2$yday))
  optDayNonPups <- as.POSIXlt(seq(as.Date("2000/6/1"), by = "day", length.out = 120))[ind]
  optDayNonPupsC <- format(ISOdate(year = 2000, month = optDayNonPups$mon + 1, 
    day = optDayNonPups$mday), "%d %B")
  fits <- NULL
  for(k in 1:10000) {
    fit <- W[[9]][k]*(-20:20)/10 + W[[10]][k]*((-20:20)/10)^2
    fits <- cbind(fits, fit)
  }
  meanNonpupHourFit = apply(fits,1,mean)
  optHourNonpups = round(((-20:20)/10)[meanNonpupHourFit == max(meanNonpupHourFit)]*
    sqrt(var(d2[!is.na(d2$hour),"hour"])) + mean(d2[!is.na(d2$hour),"hour"]))*100
    ydaylim <- c(-2:2)*sqrt(var(d2$yday))+mean(d2$yday)

## ----plot-histtrend, fig.width=12, fig.height=9, echo= FALSE, include = FALSE, cache = TRUE----
  delta1 = W1[[11]]
  phi1 = W1[[12]]
  kappa1 = W1[[13]]
  delta2 = W2[[11]]
  phi2 = W2[[12]]
  kappa2 = W2[[13]]
  delta3 = W3[[11]]
  phi3 = W3[[12]]
  kappa3 = W3[[13]]
  EV1 = NULL
  EV2 = NULL
  EV3 = NULL
  for(i in 1:length(delta1)) {
    EV1 = c(EV1, eigen(matrix(c(0,kappa1[i],phi1[i],delta1[i]), nrow = 2))$values[1])
    EV2 = c(EV2, eigen(matrix(c(0,kappa2[i],phi2[i],delta2[i]), nrow = 2))$values[1])
    EV3 = c(EV3, eigen(matrix(c(0,kappa3[i],phi3[i],delta3[i]), nrow = 2))$values[1])
  }
  # optimal day of year for pup counts
  optval <- function(W, Windex, k)
  {
    ind <- W[[Windex]][k]*(-20:20)/10 + W[[Windex+1]][k]*((-20:20)/10)^2 == 
      max(W[[Windex]][k]*(-20:20)/10 + W[[Windex+1]][k]*((-20:20)/10)^2)
    ((-20:20)/10)[ind]
  }
  optXdy <- rep(NA, times = 10000)
  for(j in 1:10000) optXdy[j] <- optval(W1, 3,j)
  optXhr <- rep(NA, times = 10000)
  for(j in 1:10000) optXhr[j] <- optval(W1, 5,j)
  optYdy <- rep(NA, times = 10000)
  for(j in 1:10000) optYdy[j] <- optval(W1, 7,j)
  optYhr <- rep(NA, times = 10000)
  for(j in 1:10000) optYhr[j] <- optval(W1, 9,j)
  Xdy <- mean(optXdy)
  Xhr <- mean(optXhr)
  Ydy <- mean(optYdy)
  Yhr <- mean(optYhr)
  trend5yr1 <- NULL
  trend10yr1 <- NULL
  trend15yr1 <- NULL
  fits <- NULL
  k <- 1
  for(k in 1:10000) {
    delta <- W1[[11]][,k]
    phi <- W1[[12]][,k]
    kappa1 <- W1[[13]][,k]
    thetaY <- rep(NA, times = 30)
    thetaX <- thetaY
    thetaY[1] <- W1[[2]][k]
    thetaX[1] <- W1[[1]][k]
    Hsamp <- W1[[16]][,k]
	  rho = W1[[18]][k]
    for(i in 2:30) {
      thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1]- Hsamp[i]
      thetaX[i] <- exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1]
    }
    fit <- thetaX + W1[[3]][k]*optXdy[k] + W1[[4]][k]*optXdy[k]^2 + 
       W1[[5]][k]*optXhr[k] + W1[[6]][k]*optXhr[k]^2 +
       thetaY + W1[[7]][k]*optYdy[k] + W1[[8]][k]*optYdy[k]^2 +
       W1[[9]][k]*optYhr[k] + W1[[10]][k]*optYhr[k]^2
    fits <- cbind(fits, fit)
    trend5yr1 <- c(trend5yr1,coef(lm(y ~ x, 
      data = data.frame(x = 1:5, y = fit[26:30])))[2])
    trend10yr1 <- c(trend10yr1, coef(lm(y ~ x, 
      data = data.frame(x = 1:10, y = fit[21:30])))[2])
    trend15yr1 <- c(trend15yr1, coef(lm(y ~ x, 
      data = data.frame(x = 1:15, y = fit[16:30])))[2])
  }
  optXdy <- rep(NA, times = 10000)
  for(j in 1:10000) optXdy[j] <- optval(W2, 3,j)
  optXhr <- rep(NA, times = 10000)
  for(j in 1:10000) optXhr[j] <- optval(W2, 5,j)
  optYdy <- rep(NA, times = 10000)
  for(j in 1:10000) optYdy[j] <- optval(W2, 7,j)
  optYhr <- rep(NA, times = 10000)
  for(j in 1:10000) optYhr[j] <- optval(W2, 9,j)
  Xdy <- mean(optXdy)
  Xhr <- mean(optXhr)
  Ydy <- mean(optYdy)
  Yhr <- mean(optYhr)  
  trend5yr2 <- NULL
  trend10yr2 <- NULL
  trend15yr2 <- NULL
  fits <- NULL
  k <- 1
  for(k in 1:10000) {
    delta <- W2[[11]][,k]
    phi <- W2[[12]][,k]
    kappa1 <- W2[[13]][,k]
    thetaY <- rep(NA, times = 30)
    thetaX <- thetaY
    thetaY[1] <- W2[[2]][k]
    thetaX[1] <- W2[[1]][k]
    Hsamp <- W2[[16]][,k]
	  rho = W2[[18]][k]
    for(i in 2:30) {
      thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1]- Hsamp[i]
      thetaX[i] <- exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1]
	  }
    fit <- thetaX + W2[[3]][k]*optXdy[k] + W2[[4]][k]*optXdy[k]^2 + 
       W2[[5]][k]*optXhr[k] + W2[[6]][k]*optXhr[k]^2 +
       thetaY + W2[[7]][k]*optYdy[k] + W2[[8]][k]*optYdy[k]^2 +
       W2[[9]][k]*optYhr[k] + W2[[10]][k]*optYhr[k]^2
    fits <- cbind(fits, fit)
    trend5yr2 <- c(trend5yr2,coef(lm(y ~ x, 
      data = data.frame(x = 1:5, y = fit[26:30])))[2])
    trend10yr2 <- c(trend10yr2, coef(lm(y ~ x, 
      data = data.frame(x = 1:10, y = fit[21:30])))[2])
    trend15yr2 <- c(trend15yr2, coef(lm(y ~ x, 
      data = data.frame(x = 1:15, y = fit[16:30])))[2])
  }
  optXdy <- rep(NA, times = 10000)
  for(j in 1:10000) optXdy[j] <- optval(W3, 3,j)
  optXhr <- rep(NA, times = 10000)
  for(j in 1:10000) optXhr[j] <- optval(W3, 5,j)
  optYdy <- rep(NA, times = 10000)
  for(j in 1:10000) optYdy[j] <- optval(W3, 7,j)
  optYhr <- rep(NA, times = 10000)
  for(j in 1:10000) optYhr[j] <- optval(W3, 9,j)
  Xdy <- mean(optXdy)
  Xhr <- mean(optXhr)
  Ydy <- mean(optYdy)
  Yhr <- mean(optYhr)  
  trend5yr3 <- NULL
  trend10yr3 <- NULL
  trend15yr3 <- NULL
  fits <- NULL
  k <- 1
  for(k in 1:10000) {
    delta <- W3[[11]][,k]
    phi <- W3[[12]][,k]
    kappa1 <- W3[[13]][,k]
    thetaY <- rep(NA, times = 30)
    thetaX <- thetaY
    thetaY[1] <- W3[[2]][k]
    thetaX[1] <- W3[[1]][k]
    Hsamp <- W3[[16]][,k]
	  rho = W3[[18]][k]
    for(i in 2:30) {
      thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1]- Hsamp[i]
      thetaX[i] <- exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1]
		}
    fit <- thetaX + W3[[3]][k]*optXdy[k] + W3[[4]][k]*optXdy[k]^2 + 
       W3[[5]][k]*optXhr[k] + W3[[6]][k]*optXhr[k]^2 +
       thetaY + W3[[7]][k]*optYdy[k] + W3[[8]][k]*optYdy[k]^2 +
       W3[[9]][k]*optYhr[k] + W3[[10]][k]*optYhr[k]^2
    fits <- cbind(fits, fit)
    trend5yr3 <- c(trend5yr3,coef(lm(y ~ x, 
      data = data.frame(x = 1:5, y = fit[26:30])))[2])
    trend10yr3 <- c(trend10yr3, coef(lm(y ~ x, 
      data = data.frame(x = 1:10, y = fit[21:30])))[2])
    trend15yr3 <- c(trend15yr3, coef(lm(y ~ x, 
      data = data.frame(x = 1:15, y = fit[16:30])))[2])
  }

  layout(matrix(1:12, nrow = 3, ncol = 4, byrow = TRUE), heights = c(1.19, 1, 1))
  postDensEV1 = density(EV1, bw = .01)
  postDensEV2 = density(EV2, bw = .01)
  postDensEV3 = density(EV3, bw = .01)
  postDens15.1 = density(trend15yr1, bw = 3)
  postDens10.1 = density(trend10yr1, bw = 3)
  postDens5.1 = density(trend5yr1, bw = 3)
  postDens15.2 = density(trend15yr2, bw = 3)
  postDens10.2 = density(trend10yr2, bw = 3)
  postDens5.2 = density(trend5yr2, bw = 3)  
  postDens15.3 = density(trend15yr3, bw = 3)
  postDens10.3 = density(trend10yr3, bw = 3)
  postDens5.3 = density(trend5yr3, bw = 3)
  xMin = min(postDens15.1$x, postDens10.1$x, postDens5.1$x,
             postDens15.2$x, postDens10.2$x, postDens5.2$x,
             postDens15.3$x, postDens10.3$x, postDens5.3$x)
  xMax = max(postDens15.1$x, postDens10.1$x, postDens5.1$x,
             postDens15.2$x, postDens10.2$x, postDens5.2$x,
             postDens15.3$x, postDens10.3$x, postDens5.3$x)
  xMinEV = min(postDensEV1$x, postDensEV2$x, postDensEV3$x)
  xMaxEV = max(postDensEV1$x, postDensEV2$x, postDensEV3$x)
  par(mar = c(5,5,5,0))
  plot(postDensEV1$x, postDensEV1$y, type = 'l', lwd = 3, xlab = '',
    ylab = '', cex.lab = 2.5, cex.axis = 1.5,
    main = 'All Years', cex.main = 2.5, xlim = c(xMinEV, xMaxEV))
  plot(postDens15.1$x, postDens15.1$y, type = 'l', lwd = 3, xlab = '', 
    ylab = '', cex.lab = 2.5, cex.axis = 1.5,
    xlim = c(-60, 100), main = '15-year', cex.main = 2.5)
  plot(postDens10.1$x, postDens10.1$y, type = 'l', lwd = 3, xlab = '', 
    ylab = '', cex.lab = 2.5, cex.axis = 1.5, xlim = c(-60, 100),
    main = '10-year', cex.main = 2.5)
  plot(postDens5.1$x, postDens5.1$y, type = 'l', lwd = 3, xlab = '', 
    ylab = '', cex.lab = 2.5, cex.axis = 1.5, xlim = c(-60, 100),
    main = '5-year', cex.main = 2.5)
  par(mar = c(5,5,1,0))
  plot(postDensEV2$x, postDensEV2$y, type = 'l', lwd = 3, xlab = '',
    ylab = 'Posterior Density', cex.lab = 2.8, cex.axis = 1.5,
    main = '', cex.main = 2.5, xlim = c(xMinEV, xMaxEV))
  plot(postDens15.2$x, postDens15.2$y, type = 'l', lwd = 3, xlab = '', 
    ylab = '', cex.lab = 2.5, cex.axis = 1.5,
    xlim = c(-60, 100), main = '', cex.main = 2.5)
  plot(postDens10.2$x, postDens10.2$y, type = 'l', lwd = 3, xlab = '', 
    ylab = '', cex.lab = 2.5, cex.axis = 1.5, xlim = c(-60, 100),
    main = '', cex.main = 2.5)
  plot(postDens5.2$x, postDens5.2$y, type = 'l', lwd = 3, xlab = '', 
    ylab = '', cex.lab = 2.5, cex.axis = 1.5, xlim = c(-60, 100),
    main = '', cex.main = 2.5)
  par(mar = c(5,5,1,0))
  plot(postDensEV3$x, postDensEV3$y, type = 'l', lwd = 3, xlab = 'Growth Rate',
    ylab = '', cex.lab = 2.5, cex.axis = 1.5,
    main = '', cex.main = 2.5, xlim = c(xMinEV, xMaxEV))
  plot(postDens15.3$x, postDens15.3$y, type = 'l', lwd = 3, xlab = 'Trend', 
    ylab = '', cex.lab = 2.5, cex.axis = 1.5,
    xlim = c(-60, 100), main = '', cex.main = 2.5)
  plot(postDens10.3$x, postDens10.3$y, type = 'l', lwd = 3, xlab = 'Trend', 
    ylab = '', cex.lab = 2.5, cex.axis = 1.5, xlim = c(-60, 100),
    main = '', cex.main = 2.5)
  plot(postDens5.3$x, postDens5.3$y, type = 'l', lwd = 3, xlab = 'Trend', 
    ylab = '', cex.lab = 2.5, cex.axis = 1.5, xlim = c(-60, 100),
    main = '', cex.main = 2.5)

## ----quants, echo = FALSE, include = FALSE, cache = TRUE-----------------
optval <- function(W, Windex, k)
{
  ind <- W[[Windex]][k]*(-20:20)/10 + W[[Windex+1]][k]*((-20:20)/10)^2 == 
    max(W[[Windex]][k]*(-20:20)/10 + W[[Windex+1]][k]*((-20:20)/10)^2)
  ((-20:20)/10)[ind]
}
optXdy <- rep(NA, times = 10000)
for(j in 1:10000) optXdy[j] <- optval(W2, 3,j)
optXhr <- rep(NA, times = 10000)
for(j in 1:10000) optXhr[j] <- optval(W2, 5,j)
optYdy <- rep(NA, times = 10000)
for(j in 1:10000) optYdy[j] <- optval(W2, 7,j)
optYhr <- rep(NA, times = 10000)
for(j in 1:10000) optYhr[j] <- optval(W2, 9,j)
Xdy <- mean(optXdy)
Xhr <- mean(optXhr)
Ydy <- mean(optYdy)
Yhr <- mean(optYhr)
fits = NULL
totsum <- rep(0, times = 30)
k <- 1
for(k in 1:10000) {
  delta <- W2[[11]][,k]
  phi <- W2[[12]][,k]
  kappa1 <- W2[[13]][,k]
  thetaY <- rep(NA, times = 30)
  thetaX <- thetaY
  thetaY[1] <- W2[[2]][k]
  thetaX[1] <- W2[[1]][k]
  Hsamp <- W2[[16]][,k]
	rho = W2[[18]][k]
  for(i in 2:30) {
    thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1]- Hsamp[i]
    thetaX[i] <- exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1]
	}
  fit <- thetaX + W2[[3]][k]*optXdy[k] + W2[[4]][k]*optXdy[k]^2 + 
     W2[[5]][k]*optXhr[k] + W2[[6]][k]*optXhr[k]^2 +
     thetaY + W2[[7]][k]*optYdy[k] + W2[[8]][k]*optYdy[k]^2 +
     W2[[9]][k]*optYhr[k] + W2[[10]][k]*optYhr[k]^2
  fits <- cbind(fits, fit)
  totsum <- totsum + fit
}
W2E2013 = (totsum/10000)[30]
W2l2.5 = apply(fits,1,quantile, p = .025)[30]
W2u97.5 = apply(fits,1,quantile, p = .975)[30]
W2E2008 = (totsum/10000)[25]
W22008l2.5 = apply(fits,1,quantile, p = .025)[25]
W22008u97.5 = apply(fits,1,quantile, p = .975)[25]

optXdy <- rep(NA, times = 10000)
for(j in 1:10000) optXdy[j] <- optval(W3, 3,j)
optXhr <- rep(NA, times = 10000)
for(j in 1:10000) optXhr[j] <- optval(W3, 5,j)
optYdy <- rep(NA, times = 10000)
for(j in 1:10000) optYdy[j] <- optval(W3, 7,j)
optYhr <- rep(NA, times = 10000)
for(j in 1:10000) optYhr[j] <- optval(W3, 9,j)
Xdy <- mean(optXdy)
Xhr <- mean(optXhr)
Ydy <- mean(optYdy)
Yhr <- mean(optYhr)
fits = NULL
totsum <- rep(0, times = 30)
k <- 1
for(k in 1:10000) {
  delta <- W3[[11]][,k]
  phi <- W3[[12]][,k]
  kappa1 <- W3[[13]][,k]
  thetaY <- rep(NA, times = 30)
  thetaX <- thetaY
  thetaY[1] <- W3[[2]][k]
  thetaX[1] <- W3[[1]][k]
  Hsamp <- W3[[16]][,k]
	rho = W3[[18]][k]
  for(i in 2:30) {
    thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1]- Hsamp[i]
    thetaX[i] <- exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1]
	}
  fit <- thetaX + W3[[3]][k]*optXdy[k] + W3[[4]][k]*optXdy[k]^2 + 
     W3[[5]][k]*optXhr[k] + W3[[6]][k]*optXhr[k]^2 +
     thetaY + W3[[7]][k]*optYdy[k] + W3[[8]][k]*optYdy[k]^2 +
     W3[[9]][k]*optYhr[k] + W3[[10]][k]*optYhr[k]^2
  fits <- cbind(fits, fit)
  totsum <- totsum + fit
}
W3E2013 = (totsum/10000)[30]
W3l2.5 = apply(fits,1,quantile, p = .025)[30]
W3u97.5 = apply(fits,1,quantile, p = .975)[30]

optXdy <- rep(NA, times = 10000)
for(j in 1:10000) optXdy[j] <- optval(W4, 3,j)
optXhr <- rep(NA, times = 10000)
for(j in 1:10000) optXhr[j] <- optval(W4, 5,j)
optYdy <- rep(NA, times = 10000)
for(j in 1:10000) optYdy[j] <- optval(W4, 7,j)
optYhr <- rep(NA, times = 10000)
for(j in 1:10000) optYhr[j] <- optval(W4, 9,j)
Xdy <- mean(optXdy)
Xhr <- mean(optXhr)
Ydy <- mean(optYdy)
Yhr <- mean(optYhr)
fits = NULL
totsum <- rep(0, times = 30)
k <- 1
for(k in 1:10000) {
  delta <- W4[[11]][,k]
  phi <- W4[[12]][,k]
  kappa1 <- W4[[13]][,k]
  thetaY <- rep(NA, times = 30)
  thetaX <- thetaY
  thetaY[1] <- W4[[2]][k]
  thetaX[1] <- W4[[1]][k]
  Hsamp <- W4[[16]][,k]
	rho = W4[[18]][k]
  for(i in 2:30) {
    thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1]- Hsamp[i]
    thetaX[i] <- exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1]
	}
	fit <- thetaX + W4[[3]][k]*optXdy[k] + W4[[4]][k]*optXdy[k]^2 + 
     W4[[5]][k]*optXhr[k] + W4[[6]][k]*optXhr[k]^2 +
     thetaY + W4[[7]][k]*optYdy[k] + W4[[8]][k]*optYdy[k]^2 +
     W4[[9]][k]*optYhr[k] + W4[[10]][k]*optYhr[k]^2
  fits <- cbind(fits, fit)
  totsum <- totsum + fit
}
W4E2013 = (totsum/10000)[30]
W4l2.5 = apply(fits,1,quantile, p = .025)[30]
W4u97.5 = apply(fits,1,quantile, p = .975)[30]

## ----echo = FALSE, results = 'asis'--------------------------------------
print(survTable, type="latex",include.rownames=FALSE, caption.placement = 'top',
     hline.after = c(-1,-1,0,0,nrow(survTable)), only.contents = TRUE)

## ----plot-demparms, fig.width=12, fig.height=14, echo=FALSE, include = FALSE, cache = TRUE----
  layout(matrix(1:12, nrow = 4, ncol = 3))

  #nonpup survivorship prior

  a1 = 138.45
  b1 = 18.38
  ymax = max(dbeta((0:10000)/10000,a1,b1))
  postDens = density(W1[[11]], bw = .008)
  ymax = max(ymax,postDens$y)
  par(mar = c(6,6,1,1))
  plot((0:10000)/10000, dbeta((0:10000)/10000,a1,b1), type = "l", lwd = 4,
    xlab = expression(delta), ylab = "Density", cex.lab = 3.2, cex.axis = 1.8,
    ylim = c(0,ymax), lty = 2)
  lines(postDens$x,postDens$y, lwd = 4)


  # nonpup fecundity and survival to time of surveys

  a1 = 15.17
  b1 = 53.09
  ymax = max(dbeta((0:10000)/10000,a1,b1))
  postDens = density(W1[[12]], bw = .008)
  ymax = max(ymax,postDens$y)
  par(mar = c(6,6,1,1))
  plot((0:10000)/10000, dbeta((0:10000)/10000,a1,b1), type = "l", lwd = 4,
    xlab = expression(phi), ylab = "Density", cex.lab = 3.2, cex.axis = 1.8,
    ylim = c(0, ymax), lty = 2)
  lines(postDens$x,postDens$y, lwd = 4)

  # nonpup density dependent fecundity and survival to time of surveys

  postDens1 = density(exp(-W1[[18]]*.3), bw = .004)
	postDens2 = density(exp(-W1[[18]]*.6), bw = .004)
	postDens3 = density(exp(-W1[[18]]*1.2), bw = .004)
  ymax = max(postDens1$y)
  par(mar = c(6,6,1,1))
  plot(postDens1$x, postDens1$y, type = "l", lwd = 4,
    xlab = expression(italic(c)), ylab = "Density", cex.lab = 3.2, cex.axis = 1.8,
    ylim = c(0, ymax), xlim = c(0.5, 1))
  lines(postDens2$x,postDens2$y, lwd = 4, lty = 2)
  lines(postDens3$x,postDens3$y, lwd = 4, lty = 3)
	legend(0.5,ymax, legend = c(expression(italic(N) == 300), 
			expression(italic(N) == 600), expression(italic(N) == 1200)),
     lty = 1:3, cex = 2.5, lwd = c(3,3,4))


  # pup survival from surveys to nonpup status

  a1 = 21.6
  b1 = 7.13
  ymax = max(dbeta((0:10000)/10000,a1,b1))
  postDens = density(W1[[13]], bw = .01)
  ymax = max(ymax,postDens$y)
  par(mar = c(6,6,1,1))
  plot((0:10000)/10000, dbeta((0:10000)/10000,a1,b1), type = "l", lwd = 4,
    xlab = expression(kappa), ylab = "Density", cex.lab = 3.2, cex.axis = 1.8,
    ylim = c(0, ymax), lty = 2)
  lines(postDens$x,postDens$y, lwd = 4)

  #nonpup survivorship prior

  a1 = 40.93
  b1 = 7.36
  ymax = max(dbeta((0:10000)/10000,a1,b1))
  postDens = density(W2[[11]], bw = .008)
  ymax = max(ymax,postDens$y)
  par(mar = c(6,6,1,1))
  plot((0:10000)/10000, dbeta((0:10000)/10000,a1,b1), type = "l", lwd = 4,
    xlab = expression(delta), ylab = "Density", cex.lab = 3.2, cex.axis = 1.8,
    ylim = c(0,ymax), lty = 2)
  lines(postDens$x,postDens$y, lwd = 4)


  # nonpup fecundity and survival to time of surveys

  a1 = 3.37
  b1 = 12.44
  ymax = max(dbeta((0:10000)/10000,a1,b1))
  postDens = density(W2[[12]], bw = .008)
  ymax = max(ymax,postDens$y)
  par(mar = c(6,6,1,1))
  plot((0:10000)/10000, dbeta((0:10000)/10000,a1,b1), type = "l", lwd = 4,
    xlab = expression(phi), ylab = "Density", cex.lab = 3.2, cex.axis = 1.8,
    ylim = c(0, ymax), lty = 2)
  lines(postDens$x,postDens$y, lwd = 4)

  # nonpup density dependent fecundity and survival to time of surveys

  postDens1 = density(exp(-W2[[18]]*.3), bw = .004)
	postDens2 = density(exp(-W2[[18]]*.6), bw = .004)
	postDens3 = density(exp(-W2[[18]]*1.2), bw = .004)
  ymax = max(postDens1$y)
  par(mar = c(6,6,1,1))
  plot(postDens1$x, postDens1$y, type = "l", lwd = 4,
    xlab = expression(italic(c)), ylab = "Density", cex.lab = 3.2, cex.axis = 1.8,
    ylim = c(0, ymax), xlim = c(0.5, 1))
  lines(postDens2$x,postDens2$y, lwd = 4, lty = 2)
  lines(postDens3$x,postDens3$y, lwd = 4, lty = 3)
	legend(0.5,ymax, legend = c(expression(italic(N) == 300), 
			expression(italic(N) == 600), expression(italic(N) == 1200)),
     lty = 1:3, cex = 2.5, lwd = c(3,3,4))


  # pup survival from surveys to nonpup status

  a1 = 4.9
  b1 = 1.93
  ymax = max(dbeta((0:10000)/10000,a1,b1))
  postDens = density(W2[[13]], bw = .01)
  ymax = max(ymax,postDens$y)
  par(mar = c(6,6,1,1))
  plot((0:10000)/10000, dbeta((0:10000)/10000,a1,b1), type = "l", lwd = 4,
    xlab = expression(kappa), ylab = "Density", cex.lab = 3.2, cex.axis = 1.8,
    ylim = c(0, ymax), lty = 2)
  lines(postDens$x,postDens$y, lwd = 4)

  #nonpup survivorship prior

  ymax = .5
  postDens = density(W3[[11]], bw = .008)
  ymax = max(ymax,postDens$y)
  par(mar = c(6,6,1,1))
  plot(c(.5,1), c(2,2), type = "l", lwd = 4,
    xlab = expression(delta), ylab = "Density", cex.lab = 3.2, cex.axis = 1.8,
    ylim = c(0,ymax), xlim = c(0,1), lty = 2)
  lines(c(0,.5),c(0,0), lty = 2, lwd = 4)
  lines(c(.5,.5),c(0,2), lty = 2, lwd = 4)
  lines(postDens$x,postDens$y, lwd = 4)

  # nonpup fecundity and survival to time of surveys

  ymax = .5
  postDens = density(W3[[12]], bw = .008)
  ymax = max(ymax,postDens$y)
  par(mar = c(6,6,1,1))
  plot(c(0,.5), c(2,2), type = "l", lwd = 4,
    xlab = expression(phi), ylab = "Density", cex.lab = 3.2, cex.axis = 1.8,
    ylim = c(0, ymax), xlim = c(0,1), lty = 2)
  lines(c(.5,1),c(0,0), lty = 2, lwd = 4)
  lines(c(.5,.5),c(0,2), lty = 2, lwd = 4)
  lines(postDens$x,postDens$y, lwd = 4)

  # nonpup density dependent fecundity and survival to time of surveys

  postDens1 = density(exp(-W3[[18]]*.3), bw = .004)
	postDens2 = density(exp(-W3[[18]]*.6), bw = .004)
	postDens3 = density(exp(-W3[[18]]*1.2), bw = .004)
  ymax = max(postDens1$y)
  par(mar = c(6,6,1,1))
  plot(postDens1$x, postDens1$y, type = "l", lwd = 4,
    xlab = expression(italic(c)), ylab = "Density", cex.lab = 3.2, cex.axis = 1.8,
    ylim = c(0, ymax), xlim = c(0.5, 1))
  lines(postDens2$x,postDens2$y, lwd = 4, lty = 2)
  lines(postDens3$x,postDens3$y, lwd = 4, lty = 3)
	legend(0.5,ymax, legend = c(expression(italic(N) == 300), 
			expression(italic(N) == 600), expression(italic(N) == 1200)),
     lty = 1:3, cex = 2.5, lwd = c(3,3,4))

  # pup survival from surveys to nonpup status

  ymax = .5
  postDens = density(W3[[13]], bw = .01)
  ymax = max(ymax,postDens$y)
  par(mar = c(6,6,1,1))
  plot(c(.5,1), c(2,2), type = "l", lwd = 4,
    xlab = expression(kappa), ylab = "Density", cex.lab = 3.2, cex.axis = 1.8,
    ylim = c(0, ymax), xlim = c(0,1), lty = 2)
  lines(c(0,.5),c(0,0), lty = 2, lwd = 4)
  lines(c(.5,.5),c(0,2), lty = 2, lwd = 4)
  lines(postDens$x,postDens$y, lwd = 4)

## ----plot-covariates, fig.width=8, fig.height=10, echo= FALSE, include = FALSE, cache = TRUE, dev='jpeg'----
  layout(matrix(1:4, nrow = 2, ncol = 2, byrow = TRUE))
  W = W2
  # A Pups and Day-of-Year
  ydaylim <- c(-2:2)*sqrt(var(d2$yday))+mean(d2$yday)
	ydaylim = 110 + (0:4)*48
  par(mar = c(5,5,3,1))
  plot(c(-2,2),c(-20,20), type = "n",  ylab = "Pup Counts", xaxt = "n", xlab = "Day of Year",
    cex.lab = 2, cex.axis = 1.5)
  fits <- NULL
  for(k in 1:10000) {
     fit <- W[[3]][k]*(-20:20)/10 + W[[4]][k]*((-20:20)/10)^2
     lines((-20:20)/10,fit, col = rgb(.5,.5,.5,.05))
     fits <- cbind(fits, fit)
  }
  lines((-20:20)/10, apply(fits,1,mean), lwd = 3)
  axis(1, at = c(-2:2), labels = as.character(round(ydaylim)), cex.axis = 1.5)
  mtext('A', side = 3, line = 1, adj = -.25, cex = 2.5)

  meanPupFit <- apply(fits,1,mean)
  ind = as.POSIXlt(seq(as.Date("2000/6/1"), by = "day", length.out = 120))$yday == round(
    ((-20:20)/10)[meanPupFit == max(meanPupFit)]*sqrt(var(d2$yday))+mean(d2$yday))
  optDayPups <- as.POSIXlt(seq(as.Date("2000/6/1"), by = "day", length.out = 120))[ind]
  optDayPupsC <- format(ISOdate(year = 2000, month = optDayPups$mon + 1, 
    day = optDayPups$mday), "%d %B")

  # B Pups and Hour-of-Day
  hrlim <- (-2:2)*sqrt(var(d2[!is.na(d2$hour),"hour"])) + mean(d2[!is.na(d2$hour),"hour"])
  hrlim <- round(hrlim*10)/10
	hrlim = 9 + (0:4)*2.5
  par(mar = c(5,5,3,1))
  plot(c(-2,2),c(-40,40), type = "n", xaxt = "n", xlab = "Hour of Day", ylab = "Pup Counts",
    cex.lab = 2, cex.axis = 1.5)
  fits <- NULL
  for(k in 1:10000) {
    fit <- W[[5]][k]*(-20:20)/10 + W[[6]][k]*((-20:20)/10)^2
    lines((-20:20)/10, fit, col = rgb(.5,.5,.5,.05))
    fits <- cbind(fits, fit)
  }
  lines((-20:20)/10, apply(fits,1,mean), lwd = 3)
  axis(1, at = c(-2:2), labels = as.character(round(hrlim,1)), cex.axis = 1.5)
  mtext('B', side = 3, line = 1, adj = -.25, cex = 2.5)

  # C Nonpups and Day-of-Year
  par(mar = c(5,5,3,1))
  plot(c(-2,2),c(-40,20), type = "n",  ylab = "Nonpup Counts", xaxt = "n", xlab = "Day of Year",
    cex.lab = 2, cex.axis = 1.5)
  fits <- NULL
  for(k in 1:10000) {
     fit <- W[[7]][k]*(-200:200)/100 + W[[8]][k]*((-200:200)/100)^2
     lines((-200:200)/100,fit, col = rgb(.5,.5,.5,.05))
     fits <- cbind(fits, fit)
  }
  lines((-200:200)/100, apply(fits,1,mean), lwd = 3)
  axis(1, at = c(-2:2), labels = as.character(round(ydaylim)), cex.axis = 1.5)
  mtext('C', side = 3, line = 1, adj = -.25, cex = 2.5)

  meanNonPupFit <- apply(fits,1,mean)
  ind = as.POSIXlt(seq(as.Date("2000/6/1"), by = "day", length.out = 120))$yday == round(
    ((-200:200)/100)[meanNonPupFit == max(meanNonPupFit)]*sqrt(var(d2$yday))+mean(d2$yday))
  optDayNonPups <- as.POSIXlt(seq(as.Date("2000/6/1"), by = "day", length.out = 120))[ind]
  optDayNonPupsC <- format(ISOdate(year = 2000, month = optDayNonPups$mon + 1, 
    day = optDayNonPups$mday), "%d %B")

  # D Nonpups and Hour-of-Day
  par(mar = c(5,5,3,1))
  plot(c(-2,2),c(-40,40), type = "n", xaxt = "n", xlab = "Hour of Day", ylab = "Nonpup Counts",
    cex.lab = 2, cex.axis = 1.5)
  fits <- NULL
  for(k in 1:10000) {
    fit <- W[[9]][k]*(-20:20)/10 + W[[10]][k]*((-20:20)/10)^2
    lines((-20:20)/10, fit, col = rgb(.5,.5,.5,.05))
    fits <- cbind(fits, fit)
  }
  lines((-20:20)/10, apply(fits,1,mean), lwd = 3)
  axis(1, at = c(-2:2), labels = as.character(round(hrlim,.5)), cex.axis = 1.5)
  mtext('D', side = 3, line = 1, adj = -.25, cex = 2.5)

## ----plot-trendlines, fig.width=12, fig.height=12, echo= FALSE, include = FALSE, cache = TRUE, dev='jpeg'----
layout(matrix(1:6, nrow = 3, ncol = 2, byrow = TRUE), width = c(2,1))

# absolute time-line of totals
total <- apply(cbind(d2$pupcount+ d2$aducount,d2$totcount),1,sum, na.rm= TRUE) 
fyr <- d2$yr + d2$yday/365

# optimal day of year for pup counts
optval <- function(W, Windex, k)
{
  ind <- W[[Windex]][k]*(-20:20)/10 + W[[Windex+1]][k]*((-20:20)/10)^2 == 
    max(W[[Windex]][k]*(-20:20)/10 + W[[Windex+1]][k]*((-20:20)/10)^2)
  ((-20:20)/10)[ind]
}
optXdy <- rep(NA, times = 10000)
for(j in 1:10000) optXdy[j] <- optval(W1, 3,j)
optXhr <- rep(NA, times = 10000)
for(j in 1:10000) optXhr[j] <- optval(W1, 5,j)
optYdy <- rep(NA, times = 10000)
for(j in 1:10000) optYdy[j] <- optval(W1, 7,j)
optYhr <- rep(NA, times = 10000)
for(j in 1:10000) optYhr[j] <- optval(W1, 9,j)
Xdy <- mean(optXdy)
Xhr <- mean(optXhr)
Ydy <- mean(optYdy)
Yhr <- mean(optYhr)

par(mar = c(6,6,4,2))
plot(c(1985,2014),c(0,800), type = "n", xlab = "Year", ylab = "Seal Counts",
  cex.main = 2, cex.axis = 1.5, cex.lab = 2.5)
totsum <- rep(0, times = 30)
trend5yr1 <- NULL
trend10yr1 <- NULL
trend15yr1 <- NULL
fits <- NULL
k <- 1
for(k in 1:10000) {
  delta <- W1[[11]][,k]
  phi <- W1[[12]][,k]
  kappa1 <- W1[[13]][,k]
  thetaY <- rep(NA, times = 30)
  thetaX <- thetaY
  thetaY[1] <- W1[[2]][k]
  thetaX[1] <- W1[[1]][k]
  Hsamp <- W1[[16]][,k]
	rho <- W1[[18]][k]
  for(i in 2:30) {
    thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1]- Hsamp[i]
    thetaX[i] <-  exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1]
  }
  fit <- thetaX + W1[[3]][k]*optXdy[k] + W1[[4]][k]*optXdy[k]^2 + 
     W1[[5]][k]*optXhr[k] + W1[[6]][k]*optXhr[k]^2 +
     thetaY + W1[[7]][k]*optYdy[k] + W1[[8]][k]*optYdy[k]^2 +
     W1[[9]][k]*optYhr[k] + W1[[10]][k]*optYhr[k]^2
  fits <- cbind(fits, fit)
  lines((1984:2013) + .56, fit, col = rgb(.5,.5,.5,.05))
  totsum <- totsum + fit
  trend5yr1 <- c(trend5yr1,coef(lm(y ~ x, 
    data = data.frame(x = 1:5, y = fit[26:30])))[2])
  trend10yr1 <- c(trend10yr1, coef(lm(y ~ x, 
    data = data.frame(x = 1:10, y = fit[21:30])))[2])
  trend15yr1 <- c(trend15yr1, coef(lm(y ~ x, 
    data = data.frame(x = 1:15, y = fit[16:30])))[2])
}
points(1983 + fyr,total, pch = 19)
avefit <- totsum/10000
lines((1984:2013) + .56, avefit, lwd = 4)
lines((1984:2013) + .56, apply(fits,1,quantile, p = .025), lwd = 3, lty = 2)
lines((1984:2013) + .56, apply(fits,1,quantile, p = .975), lwd = 3, lty = 2)
lines((1984:2013) + .56, avefit*mean(W1[['psi']]), lwd = 2)
mtext('A', side = 3, line = 1, adj = -.1, cex = 2.5)

priorVals =  exp((0:100))/2.68812e+41
postDens = density(as.vector(W1[["psi"]]), bw = .01)
maxY = max(priorVals, postDens$y)
par(mar = c(6,6,4,1))
plot(c(0,1), c(0,30), type = "n", ylab = "Density",
   xlab = expression(psi), cex.lab = 2.5, cex.axis = 1.5)
lines((1:99)/100, dbeta((1:99)/100, .5*50, (1 - .5)*50), lty = 2, lwd = 4)
lines(postDens$x, postDens$y, type = "l", lwd = 4)
mtext('B', side = 3, line = 1, adj = -.25, cex = 2.5)

optXdy <- rep(NA, times = 10000)
for(j in 1:10000) optXdy[j] <- optval(W2, 3,j)
optXhr <- rep(NA, times = 10000)
for(j in 1:10000) optXhr[j] <- optval(W2, 5,j)
optYdy <- rep(NA, times = 10000)
for(j in 1:10000) optYdy[j] <- optval(W2, 7,j)
optYhr <- rep(NA, times = 10000)
for(j in 1:10000) optYhr[j] <- optval(W2, 9,j)
Xdy <- mean(optXdy)
Xhr <- mean(optXhr)
Ydy <- mean(optYdy)
Yhr <- mean(optYhr)

par(mar = c(6,6,4,2))
plot(c(1985,2014),c(0,800), type = "n", xlab = "Year", ylab = "Seal Counts",
  cex.main = 2, cex.axis = 1.5, cex.lab = 2.5)
totsum <- rep(0, times = 30)
trend5yr2 <- NULL
trend10yr2 <- NULL
trend15yr2 <- NULL
fits <- NULL
k <- 1
for(k in 1:10000) {
  delta <- W2[[11]][,k]
  phi <- W2[[12]][,k]
  kappa1 <- W2[[13]][,k]
  thetaY <- rep(NA, times = 30)
  thetaX <- thetaY
  thetaY[1] <- W2[[2]][k]
  thetaX[1] <- W2[[1]][k]
  Hsamp <- W2[[16]][,k]
	rho = W2[[18]][k]
  for(i in 2:30) {
    thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1]- Hsamp[i]
    thetaX[i] <-  exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1]
  }
  fit <- thetaX + W2[[3]][k]*optXdy[k] + W2[[4]][k]*optXdy[k]^2 + 
     W2[[5]][k]*optXhr[k] + W2[[6]][k]*optXhr[k]^2 +
     thetaY + W2[[7]][k]*optYdy[k] + W2[[8]][k]*optYdy[k]^2 +
     W2[[9]][k]*optYhr[k] + W2[[10]][k]*optYhr[k]^2
  fits <- cbind(fits, fit)
  lines((1984:2013) + .56, fit, col = rgb(.5,.5,.5,.05))
  totsum <- totsum + fit
  trend5yr2 <- c(trend5yr2,coef(lm(y ~ x, 
    data = data.frame(x = 1:5, y = fit[26:30])))[2])
  trend10yr2 <- c(trend10yr2, coef(lm(y ~ x, 
    data = data.frame(x = 1:10, y = fit[21:30])))[2])
  trend15yr2 <- c(trend15yr2, coef(lm(y ~ x, 
    data = data.frame(x = 1:15, y = fit[16:30])))[2])
}
points(1983 + fyr,total, pch = 19)
avefit <- totsum/10000
lines((1984:2013) + .56, avefit, lwd = 4)
lines((1984:2013) + .56, apply(fits,1,quantile, p = .025), lwd = 3, lty = 2)
lines((1984:2013) + .56, apply(fits,1,quantile, p = .975), lwd = 3, lty = 2)
lines((1984:2013) + .56, avefit*mean(W2[['psi']]), lwd = 2)
mtext('C', side = 3, line = 1, adj = -.1, cex = 2.5)


postDens = density(as.vector(W2[["psi"]]), bw = .01)
par(mar = c(6,6,4,1))
plot(c(0,1), c(0,30), type = "n", ylab = "Density",
   xlab = expression(psi), cex.lab = 2.5, cex.axis = 1.5)
lines((1:99)/100, dbeta((1:99)/100, .5*50, (1 - .5)*50), lty = 2, lwd = 4)
lines(postDens$x, postDens$y, type = "l", lwd = 4)
mtext('D', side = 3, line = 1, adj = -.25, cex = 2.5)

optXdy <- rep(NA, times = 10000)
for(j in 1:10000) optXdy[j] <- optval(W3, 3,j)
optXhr <- rep(NA, times = 10000)
for(j in 1:10000) optXhr[j] <- optval(W3, 5,j)
optYdy <- rep(NA, times = 10000)
for(j in 1:10000) optYdy[j] <- optval(W3, 7,j)
optYhr <- rep(NA, times = 10000)
for(j in 1:10000) optYhr[j] <- optval(W3, 9,j)
Xdy <- mean(optXdy)
Xhr <- mean(optXhr)
Ydy <- mean(optYdy)
Yhr <- mean(optYhr)

par(mar = c(6,6,4,2))
plot(c(1985,2014),c(0,800), type = "n", xlab = "Year", ylab = "Seal Counts",
  cex.main = 2, cex.axis = 1.5, cex.lab = 2.5)
totsum <- rep(0, times = 30)
trend5yr3 <- NULL
trend10yr3 <- NULL
trend15yr3 <- NULL
fits <- NULL
k <- 1
for(k in 1:10000) {
  delta <- W3[[11]][,k]
  phi <- W3[[12]][,k]
  kappa1 <- W3[[13]][,k]
  thetaY <- rep(NA, times = 30)
  thetaX <- thetaY
  thetaY[1] <- W3[[2]][k]
  thetaX[1] <- W3[[1]][k]
  Hsamp <- W3[[16]][,k]
	rho = W3[[18]][k]
  for(i in 2:30) {
    thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1]- Hsamp[i]
    thetaX[i] <-  exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1]
  }
  fit <- thetaX + W3[[3]][k]*optXdy[k] + W3[[4]][k]*optXdy[k]^2 + 
     W3[[5]][k]*optXhr[k] + W3[[6]][k]*optXhr[k]^2 +
     thetaY + W3[[7]][k]*optYdy[k] + W3[[8]][k]*optYdy[k]^2 +
     W3[[9]][k]*optYhr[k] + W3[[10]][k]*optYhr[k]^2
  fits <- cbind(fits, fit)
  lines((1984:2013) + .56, fit, col = rgb(.5,.5,.5,.05))
  totsum <- totsum + fit
  trend5yr3 <- c(trend5yr3,coef(lm(y ~ x, 
    data = data.frame(x = 1:5, y = fit[26:30])))[2])
  trend10yr3 <- c(trend10yr3, coef(lm(y ~ x, 
    data = data.frame(x = 1:10, y = fit[21:30])))[2])
  trend15yr3 <- c(trend15yr3, coef(lm(y ~ x, 
    data = data.frame(x = 1:15, y = fit[16:30])))[2])
}
points(1983 + fyr,total, pch = 19)
avefit <- totsum/10000
lines((1984:2013) + .56, avefit, lwd = 4)
lines((1984:2013) + .56, apply(fits,1,quantile, p = .025), lwd = 3, lty = 2)
lines((1984:2013) + .56, apply(fits,1,quantile, p = .975), lwd = 3, lty = 2)
lines((1984:2013) + .56, avefit*mean(W3[['psi']]), lwd = 2)
mtext('E', side = 3, line = 1, adj = -.1, cex = 2.5)


postDens = density(as.vector(W3[['psi']]), bw = .01)
par(mar = c(6,6,4,1))
plot(c(0,1), c(0,30), type = "n", ylab = "Density",
   xlab = expression(psi), cex.lab = 2.5, cex.axis = 1.5)
lines((1:99)/100, dbeta((1:99)/100, .5*50, (1 - .5)*50), lty = 2, lwd = 4)
lines(postDens$x, postDens$y, type = "l", lwd = 4)
mtext('F', side = 3, line = 1, adj = -.25, cex = 2.5)

## ----plot-xieffect, fig.width=12, fig.height=12, echo= FALSE, include = FALSE, cache = TRUE, dev='jpeg'----
layout(matrix(1:6, nrow = 3, ncol = 2, byrow = TRUE), width = c(2,1))

# absolute time-line of totals
total <- apply(cbind(d2$pupcount+ d2$aducount,d2$totcount),1,sum, na.rm= TRUE) 
fyr <- d2$yr + d2$yday/365

# optimal day of year for pup counts
optval <- function(W, Windex, k)
{
  ind <- W[[Windex]][k]*(-20:20)/10 + W[[Windex+1]][k]*((-20:20)/10)^2 == 
    max(W[[Windex]][k]*(-20:20)/10 + W[[Windex+1]][k]*((-20:20)/10)^2)
  ((-20:20)/10)[ind]
}

# -------------------------- A ---------------------

optXdy <- rep(NA, times = 10000)
for(j in 1:10000) optXdy[j] <- optval(W6, 3,j)
optXhr <- rep(NA, times = 10000)
for(j in 1:10000) optXhr[j] <- optval(W6, 5,j)
optYdy <- rep(NA, times = 10000)
for(j in 1:10000) optYdy[j] <- optval(W6, 7,j)
optYhr <- rep(NA, times = 10000)
for(j in 1:10000) optYhr[j] <- optval(W6, 9,j)
Xdy <- mean(optXdy)
Xhr <- mean(optXhr)
Ydy <- mean(optYdy)
Yhr <- mean(optYhr)

par(mar = c(6,6,4,2))
plot(c(1985,2014),c(0,2500), type = "n", xlab = "Year", ylab = "Seal Counts",
  cex.main = 2, cex.axis = 1.5, cex.lab = 2.5)
totsum <- rep(0, times = 30)
fits <- NULL
k <- 1
for(k in 1:10000) {
  delta <- W6[[11]][,k]
  phi <- W6[[12]][,k]
  kappa1 <- W6[[13]][,k]
  thetaY <- rep(NA, times = 30)
  thetaX <- thetaY
  thetaY[1] <- W6[[2]][k]
  thetaX[1] <- W6[[1]][k]
  Hsamp <- W6[[16]][,k]
	rho = W6[[18]][k]
  for(i in 2:30) {
    thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1]- Hsamp[i]
    thetaX[i] <-  exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1]
  }
  fit <- thetaX + W6[[3]][k]*optXdy[k] + W6[[4]][k]*optXdy[k]^2 + 
     W6[[5]][k]*optXhr[k] + W6[[6]][k]*optXhr[k]^2 +
     thetaY + W6[[7]][k]*optYdy[k] + W6[[8]][k]*optYdy[k]^2 +
     W6[[9]][k]*optYhr[k] + W6[[10]][k]*optYhr[k]^2
  fits <- cbind(fits, fit)
  lines((1984:2013) + .56, fit, col = rgb(.5,.5,.5,.05))
  totsum <- totsum + fit
}
points(1983 + fyr,total, pch = 19)
avefit <- totsum/10000
lines((1984:2013) + .56, avefit, lwd = 4)
lines((1984:2013) + .56, apply(fits,1,quantile, p = .025), lwd = 3, lty = 2)
lines((1984:2013) + .56, apply(fits,1,quantile, p = .975), lwd = 3, lty = 2)
lines((1984:2013) + .56, avefit*mean(W6[['psi']]), lwd = 2)
mtext('A', side = 3, line = 1, adj = -.1, cex = 2.5)

# -------------------------- B ---------------------

postDens = density(as.vector(W6[["psi"]]), bw = .005)
par(mar = c(6,6,4,1))
plot(c(0,1), c(0,30), type = "n", ylab = "Density",
   xlab = expression(psi), cex.lab = 2.5, cex.axis = 1.5)
lines(c(0.1,0.1),c(0,1.25), lty = 2, lwd = 4)
lines(c(0.9,0.9),c(0,1.25), lty = 2, lwd = 4)
lines(c(0.1,0.9),c(1.25,1.25), lty = 2, lwd = 4)
lines(postDens$x, postDens$y, type = "l", lwd = 4)
mtext('B', side = 3, line = 1, adj = -.25, cex = 2.5)

# -------------------------- C ---------------------

optXdy <- rep(NA, times = 10000)
for(j in 1:10000) optXdy[j] <- optval(W5, 3,j)
optXhr <- rep(NA, times = 10000)
for(j in 1:10000) optXhr[j] <- optval(W5, 5,j)
optYdy <- rep(NA, times = 10000)
for(j in 1:10000) optYdy[j] <- optval(W5, 7,j)
optYhr <- rep(NA, times = 10000)
for(j in 1:10000) optYhr[j] <- optval(W5, 9,j)
Xdy <- mean(optXdy)
Xhr <- mean(optXhr)
Ydy <- mean(optYdy)
Yhr <- mean(optYhr)

par(mar = c(6,6,4,2))
plot(c(1985,2014),c(0,2500), type = "n", xlab = "Year", ylab = "Seal Counts",
  cex.main = 2, cex.axis = 1.5, cex.lab = 2.5)
totsum <- rep(0, times = 30)
fits <- NULL
k <- 1
for(k in 1:10000) {
  delta <- W5[[11]][,k]
  phi <- W5[[12]][,k]
  kappa1 <- W5[[13]][,k]
  thetaY <- rep(NA, times = 30)
  thetaX <- thetaY
  thetaY[1] <- W5[[2]][k]
  thetaX[1] <- W5[[1]][k]
  Hsamp <- W5[[16]][,k]
	rho = W5[[18]][k]
  for(i in 2:30) {
    thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1]- Hsamp[i]
    thetaX[i] <-  exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1]
  }
  fit <- thetaX + W5[[3]][k]*optXdy[k] + W5[[4]][k]*optXdy[k]^2 + 
     W5[[5]][k]*optXhr[k] + W5[[6]][k]*optXhr[k]^2 +
     thetaY + W5[[7]][k]*optYdy[k] + W5[[8]][k]*optYdy[k]^2 +
     W5[[9]][k]*optYhr[k] + W5[[10]][k]*optYhr[k]^2
  fits <- cbind(fits, fit)
  lines((1984:2013) + .56, fit, col = rgb(.5,.5,.5,.05))
  totsum <- totsum + fit
}
points(1983 + fyr,total, pch = 19)
avefit <- totsum/10000
lines((1984:2013) + .56, avefit, lwd = 4)
lines((1984:2013) + .56, apply(fits,1,quantile, p = .025), lwd = 3, lty = 2)
lines((1984:2013) + .56, apply(fits,1,quantile, p = .975), lwd = 3, lty = 2)
lines((1984:2013) + .56, avefit*mean(W5[['psi']]), lwd = 2)
mtext('C', side = 3, line = 1, adj = -.1, cex = 2.5)

# -------------------------- D ---------------------

postDens = density(as.vector(W5[["psi"]]), bw = .01)
par(mar = c(6,6,4,1))
plot(c(0,1), c(0,30), type = "n", ylab = "Density",
   xlab = expression(psi), cex.lab = 2.5, cex.axis = 1.5)
lines((0:100)/100, exp((0:100)/10)/2202.55, lty = 2, lwd = 4)
lines(postDens$x, postDens$y, type = "l", lwd = 4)
mtext('D', side = 3, line = 1, adj = -.25, cex = 2.5)

# -------------------------- E ---------------------

optXdy <- rep(NA, times = 10000)
for(j in 1:10000) optXdy[j] <- optval(W4, 3,j)
optXhr <- rep(NA, times = 10000)
for(j in 1:10000) optXhr[j] <- optval(W4, 5,j)
optYdy <- rep(NA, times = 10000)
for(j in 1:10000) optYdy[j] <- optval(W4, 7,j)
optYhr <- rep(NA, times = 10000)
for(j in 1:10000) optYhr[j] <- optval(W4, 9,j)
Xdy <- mean(optXdy)
Xhr <- mean(optXhr)
Ydy <- mean(optYdy)
Yhr <- mean(optYhr)

par(mar = c(6,6,4,2))
plot(c(1985,2014),c(0,2500), type = "n", xlab = "Year", ylab = "Seal Counts",
  cex.main = 2, cex.axis = 1.5, cex.lab = 2.5)
totsum <- rep(0, times = 30)
fits <- NULL
k <- 1
for(k in 1:10000) {
  delta <- W4[[11]][,k]
  phi <- W4[[12]][,k]
  kappa1 <- W4[[13]][,k]
  thetaY <- rep(NA, times = 30)
  thetaX <- thetaY
  thetaY[1] <- W4[[2]][k]
  thetaX[1] <- W4[[1]][k]
  Hsamp <- W4[[16]][,k]
	rho = W4[[18]][k]
  for(i in 2:30) {
    thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1]- Hsamp[i]
    thetaX[i] <-  exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1]
  }
  fit <- thetaX + W4[[3]][k]*optXdy[k] + W4[[4]][k]*optXdy[k]^2 + 
     W4[[5]][k]*optXhr[k] + W4[[6]][k]*optXhr[k]^2 +
     thetaY + W4[[7]][k]*optYdy[k] + W4[[8]][k]*optYdy[k]^2 +
     W4[[9]][k]*optYhr[k] + W4[[10]][k]*optYhr[k]^2
  fits <- cbind(fits, fit)
  lines((1984:2013) + .56, fit, col = rgb(.5,.5,.5,.05))
  totsum <- totsum + fit
}
points(1983 + fyr,total, pch = 19)
avefit <- totsum/10000
lines((1984:2013) + .56, avefit, lwd = 4)
lines((1984:2013) + .56, apply(fits,1,quantile, p = .025), lwd = 3, lty = 2)
lines((1984:2013) + .56, apply(fits,1,quantile, p = .975), lwd = 3, lty = 2)
lines((1984:2013) + .56, avefit*mean(W4[['psi']]), lwd = 2)
mtext('E', side = 3, line = 1, adj = -.1, cex = 2.5)

# -------------------------- F ---------------------

postDens = density(as.vector(W4[["psi"]]), bw = .02)
par(mar = c(6,6,4,1))
plot(c(0,1), c(0,30), type = "n", ylab = "Density",
   xlab = expression(psi), cex.lab = 2.5, cex.axis = 1.5)
lines((0:100)/100, exp((0:100)/2)/1.03694e+20, lty = 2, lwd = 4)
lines(postDens$x, postDens$y, type = "l", lwd = 4)
mtext('F', side = 3, line = 1, adj = -.25, cex = 2.5)



## ----plot-PVA, fig.width=12, fig.height=12, echo=FALSE, include=FALSE, cache=TRUE, dev='jpeg'----
  set.seed(1001)
  layout(matrix(1:6,nrow = 3, byrow = TRUE))
  par(mar = c(5,5,4,1))
  plot(c(1,130),c(-1,8), type = "n", xaxt = "n", yaxt = 'n', xlab = "Year", 
    ylab = "Abundance", cex.lab = 2, cex.axis = 1.5)
  axis(1, at = (0:13)*10, labels = as.character(round((0:13)*10+1984)), cex.axis = 1.5)
	axis(2, at = (0:7), labels = as.character(round(exp((0:7)))), cex.axis = 1.5, las = 1)
  extinct10 <- NULL
  extinct50 <- NULL
  extinct100 <- NULL
  for(k in 1:10000) {
    delta <- W1[[11]][,k]
    phi <- W1[[12]][,k]
    kappa1 <- W1[[13]][,k]
    thetaY <- rep(NA, times = 100)
    thetaX <- thetaY
    thetaY[1] <- W1[[2]][k]
    thetaX[1] <- W1[[1]][k]
    Hsamp <- W1[[16]][,k]
		rho = W1[[18]][k]
    for(i in 2:30) {
      thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1] - Hsamp[i]
      thetaX[i] <- exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1] 
    }
      for(i in 31:130) {
				yrindx = sample(1:29,1)
				MCMCindx = sample(1:10000,1)
        deltai <- W1[[11]][yrindx,MCMCindx]
        phii <- W1[[12]][yrindx,MCMCindx]
        kappai <- W1[[13]][yrindx,MCMCindx]
				rhoi <- W1[[18]][MCMCindx]
        Hi <- sample(Hsamp,1)
        thetaY[i] <- deltai*thetaY[i-1] + kappai*thetaX[i-1] - Hi*thetaY[i-1]/350
        thetaX[i] <- exp(-rhoi*(thetaY[i-1]+thetaX[i-1])/1000)*phii*thetaY[i-1]
      }
      fit <- thetaX + thetaY
      ind <- fit[31:130] <= 10
      if(any(ind) ) extinct10 <- c(extinct10, min((1:100)[ind]))
      ind <- fit[31:130] <= 50
      if(any(ind) ) extinct50 <- c(extinct50, min((1:100)[ind]))
      ind <- fit[31:130] <= 100
      if(any(ind) ) extinct100 <- c(extinct100, min((1:100)[ind]))
      lines(1:130, log(fit), col = rgb(.5,.5,.5,.04))
    
  }
  lines(c(1,130), c(log(10),log(10)), lty = 3, lwd = 3)
  lines(c(1,130), c(log(50),log(50)), lty = 2, lwd = 3)
  lines(c(1,130), c(log(100),log(100)), lty = 1, lwd = 3)
  legend(1,log(6), legend = c("100", "50", "10"),
     lty = 1:3, cex = 1.5, lwd = rep(2, times = 3))
	mtext('A', side = 3, line = 1, adj = -.1, cex = 3)
  #
  #                                   B
  #
  v10 <- sort(extinct10)
  v50 <- sort(extinct50)
  v100 <- sort(extinct100)
  pExtYr10 <- matrix(0, nrow = 100, ncol = 2)
  pExtYr50 <- matrix(0, nrow = 100, ncol = 2)
  pExtYr100 <- matrix(0, nrow = 100, ncol = 2)
  for(i in 1:100) {
    pExtYr10[i,1] <- i
    pExtYr10[i,2] <- sum(v10 <= i)/10000
    pExtYr50[i,1] <- i
    pExtYr50[i,2] <- sum(v50 <= i)/10000
    pExtYr100[i,1] <- i
    pExtYr100[i,2] <- sum(v100 <= i)/10000
 
  }
  plot(pExtYr100, type = "l", lwd = 3, xlab = "Years from 2013", 
       ylab = "Probability(<Threshold)",
       cex.lab = 2, cex.axis = 1.5, ylim = c(0,.2))
  lines(pExtYr50, lwd = 3, lty = 2)
  lines(pExtYr10, lwd = 3, lty = 3)
  legend(0,.2, legend = c("<100", "<50", "<10"),
     lty = 1:3, 
     cex = 2, lwd = rep(2, times = 3))
  mtext('B', side = 3, line = 1, adj = -.1, cex = 3)
  #
  #                                                    C
  #
  par(mar = c(5,5,4,1))
  plot(c(1,130),c(-1,8), type = "n", xaxt = "n", yaxt = 'n', xlab = "Year", 
    ylab = "Abundance", cex.lab = 2, cex.axis = 1.5)
  axis(1, at = (0:13)*10, labels = as.character(round((0:13)*10+1984)), cex.axis = 1.5)
	axis(2, at = (0:7), labels = as.character(round(exp((0:7)))), cex.axis = 1.5, las = 1)
  extinct10 <- NULL
  extinct50 <- NULL
  extinct100 <- NULL
  for(k in 1:10000) {
    delta <- W2[[11]][,k]
    phi <- W2[[12]][,k]
    kappa1 <- W2[[13]][,k]
    thetaY <- rep(NA, times = 100)
    thetaX <- thetaY
    thetaY[1] <- W2[[2]][k]
    thetaX[1] <- W2[[1]][k]
    Hsamp <- W2[[16]][,k]
		rho = W2[[18]][k]
    for(i in 2:30) {
      thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1] - Hsamp[i]
      thetaX[i] <- exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1] 
    }
      for(i in 31:130) {
				yrindx = sample(1:29,1)
				MCMCindx = sample(1:10000,1)
        deltai <- W2[[11]][yrindx,MCMCindx]
        phii <- W2[[12]][yrindx,MCMCindx]
        kappai <- W2[[13]][yrindx,MCMCindx]
				rhoi <- W2[[18]][MCMCindx]
        Hi <- sample(Hsamp,1)
        thetaY[i] <- deltai*thetaY[i-1] + kappai*thetaX[i-1] - Hi*thetaY[i-1]/350
        thetaX[i] <- exp(-rhoi*(thetaY[i-1]+thetaX[i-1])/1000)*phii*thetaY[i-1]
      }
      fit <- thetaX + thetaY
      ind <- fit[31:130] <= 10
      if(any(ind) ) extinct10 <- c(extinct10, min((1:100)[ind]))
      ind <- fit[31:130] <= 50
      if(any(ind) ) extinct50 <- c(extinct50, min((1:100)[ind]))
      ind <- fit[31:130] <= 100
      if(any(ind) ) extinct100 <- c(extinct100, min((1:100)[ind]))
      lines(1:130, log(fit), col = rgb(.5,.5,.5,.04))
  }
  lines(c(1,130), c(log(10),log(10)), lty = 3, lwd = 3)
  lines(c(1,130), c(log(50),log(50)), lty = 2, lwd = 3)
  lines(c(1,130), c(log(100),log(100)), lty = 1, lwd = 3)
  legend(1,log(6), legend = c("100", "50", "10"),
     lty = 1:3, cex = 1.5, lwd = rep(2, times = 3))
  mtext('C', side = 3, line = 1, adj = -.1, cex = 3)
  #
  #                                                    D
  #
  v10 <- sort(extinct10)
  v50 <- sort(extinct50)
  v100 <- sort(extinct100)
  pExtYr10 <- matrix(0, nrow = 100, ncol = 2)
  pExtYr50 <- matrix(0, nrow = 100, ncol = 2)
  pExtYr100 <- matrix(0, nrow = 100, ncol = 2)
  for(i in 1:100) {
    pExtYr10[i,1] <- i
    pExtYr10[i,2] <- sum(v10 <= i)/10000
    pExtYr50[i,1] <- i
    pExtYr50[i,2] <- sum(v50 <= i)/10000
    pExtYr100[i,1] <- i
    pExtYr100[i,2] <- sum(v100 <= i)/10000
 
  }
  plot(pExtYr100, type = "l", lwd = 3, xlab = "Years from 2013", 
       ylab = "Probability(<Threshold)",
       cex.lab = 2, cex.axis = 1.5, ylim = c(0,.2))
  lines(pExtYr50, lwd = 3, lty = 2)
  lines(pExtYr10, lwd = 3, lty = 3)
  legend(0,.2, legend = c("<100", "<50", "<10"),
     lty = 1:3, 
     cex = 2, lwd = rep(2, times = 3))
  mtext('D', side = 3, line = 1, adj = -.1, cex = 3)
  #
  #                                                    E
  #
  par(mar = c(5,5,4,1))
  plot(c(1,130),c(-1,8), type = "n", xaxt = "n", yaxt = 'n', xlab = "Year", 
    ylab = "Abundance", cex.lab = 2, cex.axis = 1.5)
  axis(1, at = (0:13)*10, labels = as.character(round((0:13)*10+1984)), cex.axis = 1.5)
	axis(2, at = (0:7), labels = as.character(round(exp((0:7)))), cex.axis = 1.5, las = 1)
  extinct10 <- NULL
  extinct50 <- NULL
  extinct100 <- NULL
  for(k in 1:10000) {
    delta <- W3[[11]][,k]
    phi <- W3[[12]][,k]
    kappa1 <- W3[[13]][,k]
    thetaY <- rep(NA, times = 100)
    thetaX <- thetaY
    thetaY[1] <- W3[[2]][k]
    thetaX[1] <- W3[[1]][k]
    Hsamp <- W3[[16]][,k]
		rho = W3[[18]][k]
    for(i in 2:30) {
      thetaY[i] <- delta[i-1]*thetaY[i-1] + kappa1[i-1]*thetaX[i-1] - Hsamp[i]
      thetaX[i] <- exp(-rho*(thetaY[i-1]+thetaX[i-1])/1000)*phi[i-1]*thetaY[i-1] 
    }
      for(i in 31:130) {
				yrindx = sample(1:29,1)
				MCMCindx = sample(1:10000,1)
        deltai <- W3[[11]][yrindx,MCMCindx]
        phii <- W3[[12]][yrindx,MCMCindx]
        kappai <- W3[[13]][yrindx,MCMCindx]
				rhoi <- W3[[18]][MCMCindx]
        Hi <- sample(Hsamp,1)
        thetaY[i] <- deltai*thetaY[i-1] + kappai*thetaX[i-1] - Hi*thetaY[i-1]/350
        thetaX[i] <- exp(-rhoi*(thetaY[i-1]+thetaX[i-1])/1000)*phii*thetaY[i-1]
      }
      fit <- thetaX + thetaY
      ind <- fit[31:130] <= 10
      if(any(ind) ) extinct10 <- c(extinct10, min((1:100)[ind]))
      ind <- fit[31:130] <= 50
      if(any(ind) ) extinct50 <- c(extinct50, min((1:100)[ind]))
      ind <- fit[31:130] <= 100
      if(any(ind) ) extinct100 <- c(extinct100, min((1:100)[ind]))
      lines(1:130, log(fit), col = rgb(0.5,0.5,0.5,.04))
  }
  lines(c(1,130), c(log(10),log(10)), lty = 3, lwd = 3)
  lines(c(1,130), c(log(50),log(50)), lty = 2, lwd = 3)
  lines(c(1,130), c(log(100),log(100)), lty = 1, lwd = 3)
  legend(1,log(6), legend = c("100", "50", "10"),
     lty = 1:3, cex = 1.5, lwd = rep(2, times = 3))
  mtext('E', side = 3, line = 1, adj = -.1, cex = 3)
  #
  #                                                  F
  #
  v10 <- sort(extinct10)
  v50 <- sort(extinct50)
  v100 <- sort(extinct100)
  pExtYr10 <- matrix(0, nrow = 100, ncol = 2)
  pExtYr50 <- matrix(0, nrow = 100, ncol = 2)
  pExtYr100 <- matrix(0, nrow = 100, ncol = 2)
  for(i in 1:100) {
    pExtYr10[i,1] <- i
    pExtYr10[i,2] <- sum(v10 <= i)/10000
    pExtYr50[i,1] <- i
    pExtYr50[i,2] <- sum(v50 <= i)/10000
    pExtYr100[i,1] <- i
    pExtYr100[i,2] <- sum(v100 <= i)/10000
 
  }
  plot(pExtYr100, type = "l", lwd = 3, xlab = "Years from 2013", 
       ylab = "Probability(<Threshold)",
       cex.lab = 2, cex.axis = 1.5, ylim = c(0,.2))
  lines(pExtYr50, lwd = 3, lty = 2)
  lines(pExtYr10, lwd = 3, lty = 3)
  legend(0,.2, legend = c("<100", "<50", "<10"),
     lty = 1:3, 
     cex = 2, lwd = rep(2, times = 3))
  mtext('F', side = 3, line = 1, adj = -.1, cex = 3)

## ----echo=FALSE, include = FALSE-----------------------------------------
library(xtable)
ntot = 165
novu = 143
nimp = 132
nfail = 14
ovurate = novu/ntot
pregfail = nfail/nimp

pp1 = (ovurate - pregfail)
pp2 = (1 - 0.036)
pp = pp1*pp2
vpp1 = (ovurate*(1-ovurate)/ntot + pregfail*(1-pregfail)/nimp)
vpp2 = 0.026^2

sepp = sqrt(pp1^2*vpp2 + pp2^2*vpp1 - vpp1*vpp2)

Stab = data.frame(
  parm = c("Survival Small Male Pups", 
           "Survival Large Male Pups", 
           "Survival Male 1-3", 
           "Survival Male 3+",
           "Survival Small Female Pups", 
           "Survival Large Female Pups", 
           "Survival Female 1-3", 
           "Survival Female 3+",
           "Fecundity Female 3+"),
#part b  est = c(.293, .686, .773, .866, .424, .796, .858, .920, pp),
#partb  se = c(.109, .121, .051, .050, .127, .098, .038, .037, sqrt(pp*(1-pp)/nn))
#part a
est = c(.405, .717, .782, .879, .549, .820, .865, .929, pp),
se = c(.089, .088, .035, .038, .093, .069, .027, .026, sepp)
)
SurvTable = Stab
colnames(SurvTable) = c("Parameter", "Estimate", "Standard Error")

## ----results = 'asis', echo = FALSE--------------------------------------
print(xtable(SurvTable, digits = c(1,3,3,3),
  caption = "Demographic parameters used to develop informative priors. \\label{tab:demparms}"), 
  include.rownames = FALSE, caption.placement = 'top')

## ----echo=FALSE, include = FALSE-----------------------------------------
MalePupBirth2Survey = Stab[Stab$parm == "Survival Small Male Pups","est"]/
  Stab[Stab$parm == "Survival Large Male Pups","est"]
FemalePupBirth2Survey = Stab[Stab$parm == "Survival Small Female Pups","est"]/
  Stab[Stab$parm == "Survival Large Female Pups","est"]
M = matrix(0, nrow = 8, ncol = 8)
M[2,1] = Stab[Stab$parm == "Survival Small Male Pups","est"]
M[3,2] = Stab[Stab$parm == "Survival Male 1-3","est"]
M[4,3] = Stab[Stab$parm == "Survival Male 1-3","est"]
M[4,4] = Stab[Stab$parm == "Survival Male 3+","est"]
M[6,5] = Stab[Stab$parm == "Survival Small Female Pups","est"]
M[7,6] = Stab[Stab$parm == "Survival Female 1-3","est"]
M[8,7] = Stab[Stab$parm == "Survival Female 1-3","est"]
M[8,8] = Stab[Stab$parm == "Survival Female 3+","est"]
M[1,8] = Stab[Stab$parm == "Fecundity Female 3+","est"]/2
M[5,8] = Stab[Stab$parm == "Fecundity Female 3+","est"]/2
digitMat = matrix(0, nrow = 8, ncol = 9)
digitMat[2,2] = digitMat[3,3] = digitMat[4,4] = digitMat[4,5] = digitMat[6,6] = 
  digitMat[7,7] = digitMat[8,8] = digitMat[8,9] = digitMat[1,9] = 
  digitMat[5,9] = 3
EVM0 = Re(eigen(M)$values[1])

## ----results = 'asis', echo = FALSE--------------------------------------
  print(xtable(M, digits = digitMat, align = rep("l",times = 9)), floating=FALSE, 
    hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE, 
    tabular.environment = "array", only.contents = TRUE)

## ----echo=FALSE, include = FALSE-----------------------------------------
MalePupBirth2Survey = Stab[Stab$parm == "Survival Small Male Pups","est"]/
  Stab[Stab$parm == "Survival Large Male Pups","est"]
FemalePupBirth2Survey = Stab[Stab$parm == "Survival Small Female Pups","est"]/
  Stab[Stab$parm == "Survival Large Female Pups","est"]
M = matrix(0, nrow = 8, ncol = 8)
M[2,1] = Stab[Stab$parm == "Survival Large Male Pups","est"]
M[3,2] = Stab[Stab$parm == "Survival Male 1-3","est"]
M[4,3] = Stab[Stab$parm == "Survival Male 1-3","est"]
M[4,4] = Stab[Stab$parm == "Survival Male 3+","est"]
M[6,5] = Stab[Stab$parm == "Survival Large Female Pups","est"]
M[7,6] = Stab[Stab$parm == "Survival Female 1-3","est"]
M[8,7] = Stab[Stab$parm == "Survival Female 1-3","est"]
M[8,8] = Stab[Stab$parm == "Survival Female 3+","est"]
M[1,8] = Stab[Stab$parm == "Fecundity Female 3+","est"]/2*MalePupBirth2Survey
M[5,8] = Stab[Stab$parm == "Fecundity Female 3+","est"]/2*FemalePupBirth2Survey
EVM1 = Re(eigen(M)$values[1])
eVecM1 = as.matrix(Re((eigen(M)$vectors)[,1])/sum(Re((eigen(M)$vectors)[,1])))

## ----results = 'asis', echo = FALSE--------------------------------------
  print(xtable(M, digits = digitMat, align = rep("l",times = 9)), floating=FALSE, 
    hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE, 
    tabular.environment = "array", only.contents = TRUE)

## ----results = 'asis', echo = FALSE--------------------------------------
  print(xtable(eVecM1, digits = c(0,3)), floating=FALSE, 
    hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE, 
    tabular.environment = "array", only.contents = TRUE)

## ----echo = FALSE, include = FALSE---------------------------------------
  Mvals = rep(NA, times = 8)
  Mvals[1] = Stab[Stab$parm == "Survival Large Male Pups","est"]
  Mvals[2] = Stab[Stab$parm == "Survival Male 1-3","est"]
  Mvals[3] = Stab[Stab$parm == "Survival Male 3+","est"]
  Mvals[4] = Stab[Stab$parm == "Survival Large Female Pups","est"]
  Mvals[5] = Stab[Stab$parm == "Survival Female 1-3","est"]
  Mvals[6] = Stab[Stab$parm == "Survival Female 3+","est"]
  Mvals[7] = Stab[Stab$parm == "Fecundity Female 3+","est"]/2*MalePupBirth2Survey
  Mvals[8] = Stab[Stab$parm == "Fecundity Female 3+","est"]/2*FemalePupBirth2Survey 
  deltaFunc = function(Mentries) {
    M = matrix(0, nrow = 8, ncol = 8)
    M[2,1] = Mentries[1]
    M[3,2] = Mentries[2]
    M[4,3] = Mentries[2]
    M[4,4] = Mentries[3]
    M[6,5] = Mentries[4]
    M[7,6] = Mentries[5]
    M[8,7] = Mentries[5]
    M[8,8] = Mentries[6]
    M[1,8] = Mentries[7]
    M[5,8] = Mentries[8]
    stN = Re(eigen(M)$vectors[,1])/sum(Re(eigen(M)$vectors[,1]))
    SN2N = sum(stN[c(2:4,6:8)]*c(Mentries[2],Mentries[2],Mentries[3],
      Mentries[5],Mentries[5],Mentries[6]))/
      sum(stN[c(2:4,6:8)])
    return(SN2N)
  }
  kappaFunc = function(Mentries) {
    M = matrix(0, nrow = 8, ncol = 8)
    M[2,1] = Mentries[1]
    M[3,2] = Mentries[2]
    M[4,3] = Mentries[2]
    M[4,4] = Mentries[3]
    M[6,5] = Mentries[4]
    M[7,6] = Mentries[5]
    M[8,7] = Mentries[5]
    M[8,8] = Mentries[6]
    M[1,8] = Mentries[7]
    M[5,8] = Mentries[8]
    stN = Re(eigen(M)$vectors[,1])/sum(Re(eigen(M)$vectors[,1]))
    SP2N = sum(stN[c(1,5)]*c(Mentries[1],Mentries[4]))/
      sum(stN[c(1,5)])
    return(SP2N)
  }
  phiFunc = function(Mentries) {
    M = matrix(0, nrow = 8, ncol = 8)
    M[2,1] = Mentries[1]
    M[3,2] = Mentries[2]
    M[4,3] = Mentries[2]
    M[4,4] = Mentries[3]
    M[6,5] = Mentries[4]
    M[7,6] = Mentries[5]
    M[8,7] = Mentries[5]
    M[8,8] = Mentries[6]
    M[1,8] = Mentries[7]
    M[5,8] = Mentries[8]
    stN = Re(eigen(M)$vectors[,1])/sum(Re(eigen(M)$vectors[,1]))
    FF = sum(stN[8]*(Mentries[7] + Mentries[8]))/
       sum(stN[c(2:4,6:8)])
    return(FF)
  }
  delta = deltaFunc(Mvals)
  kappa1 = kappaFunc(Mvals)
  phi = phiFunc(Mvals)
  M2x2 = matrix(c(0,phi,kappa1,delta), ncol = 2, nrow = 2, byrow = TRUE)
  M2x2eVal = Re(eigen(M2x2)$value[1])
  M2x2eVec = Re(eigen(M2x2)$vec[,1])/sum(Re(eigen(M2x2)$vec[,1]))

## ----results = 'asis', echo = FALSE--------------------------------------
  print(xtable(M2x2, digits = c(0,3,3)), floating=FALSE, 
    hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE, 
    tabular.environment = "array", only.contents = TRUE)

## ----echo = FALSE, include = FALSE---------------------------------------
  SigMentries = matrix(0, nrow = 8, ncol = 8)
  SigMentries[1,1] = Stab[Stab$parm == "Survival Large Male Pups","se"]^2
  SigMentries[2,2] = Stab[Stab$parm == "Survival Male 1-3","se"]^2
  SigMentries[3,3] = Stab[Stab$parm == "Survival Male 3+","se"]^2
  SigMentries[4,4] = Stab[Stab$parm == "Survival Large Female Pups","se"]^2
  SigMentries[5,5] = Stab[Stab$parm == "Survival Female 1-3","se"]^2
  SigMentries[6,6] = Stab[Stab$parm == "Survival Female 3+","se"]^2
  SigMentries
  varcol8 = diag(c(
    Stab[Stab$parm == "Fecundity Female 3+","se"]^2/4,
    Stab[Stab$parm == "Survival Small Male Pups","se"]^2,
    Stab[Stab$parm == "Survival Large Male Pups","se"]^2,
    Stab[Stab$parm == "Survival Small Female Pups","se"]^2,
    Stab[Stab$parm == "Survival Large Female Pups","se"]^2
    ))
  df1 = c(MalePupBirth2Survey,
    Stab[Stab$parm == "Fecundity Female 3+","est"]/2/
       Stab[Stab$parm == "Survival Large Male Pups","est"],
    Stab[Stab$parm == "Fecundity Female 3+","est"]/2*
       Stab[Stab$parm == "Survival Small Male Pups","est"]/
       Stab[Stab$parm == "Survival Large Male Pups","est"]^2,
    0,0)
  df2 = c(FemalePupBirth2Survey,  0, 0,
    Stab[Stab$parm == "Fecundity Female 3+","est"]/2/
       Stab[Stab$parm == "Survival Large Female Pups","est"],
    Stab[Stab$parm == "Fecundity Female 3+","est"]/2*
       Stab[Stab$parm == "Survival Small Female Pups","est"]/
       Stab[Stab$parm == "Survival Large Female Pups","est"]^2)
  DD = cbind(df1, df2)

  varMentries78 = t(DD) %*% varcol8 %*% DD
  SigMentries[7:8,7:8] = varMentries78
  df11 = df1[1:3]
  df12 = c(0, 0, 1)
  DD = cbind(df11, df12)
  covMentries18 = t(DD) %*% varcol8[1:3,1:3] %*% DD
  SigMentries[1,8] = SigMentries[8,1] = covMentries18[1,2]
  df21 = df2[c(1,4,5)]
  df22 = c(0, 0, 1)
  DD = cbind(df21, df22)
  covMentries58 = t(DD) %*% varcol8[c(1,4,5),c(1,4,5)] %*% DD
  SigMentries[5,8] = SigMentries[8,5] = covMentries58[1,2]
  digitMat = matrix(0, nrow = 8, ncol = 9)
  digitMat[1,2] = digitMat[2,3] = digitMat[3,4] = 
    digitMat[4,5] = digitMat[5,6] = digitMat[6,7] = 
    digitMat[7,8] = digitMat[8,9] = digitMat[1,9] = 
    digitMat[5,9] = digitMat[8,6] = digitMat[1,9] = 
    digitMat[8,2] = digitMat[7,9] = digitMat[8,8] = 5

## ----results = 'asis', echo = FALSE--------------------------------------
  print(xtable(SigMentries, digits = digitMat, align = rep("l",times = 9)), floating=FALSE, 
    hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE, 
    tabular.environment = "array", only.contents = TRUE)

## ----echo = FALSE, include = FALSE---------------------------------------
  h = 0.000001
  Ddelta = rep(0, times = 8)
  for(i in 1:8) {
    Mentriesh = Mvals
    Mentriesh[i] =  Mvals[i] + h
    Ddelta[i] = (deltaFunc(Mentriesh) - deltaFunc(Mvals))/h
  }
  Dkappa = rep(0, times = 8)
  for(i in 1:8) {
    Mentriesh = Mvals
    Mentriesh[i] =  Mvals[i] + h
    Dkappa[i] = (kappaFunc(Mentriesh) - kappaFunc(Mvals))/h
  }
  Dphi = rep(0, times = 8)
  for(i in 1:8) {
    Mentriesh = Mvals
    Mentriesh[i] =  Mvals[i] + h
    Dphi[i] = (phiFunc(Mentriesh) - phiFunc(Mvals))/h
  }
  DD = cbind(Ddelta, Dphi, Dkappa)
  covDelPhiKap = t(DD) %*% SigMentries %*% DD
  varDelta = diag(covDelPhiKap)[1]
  varPhi = diag(covDelPhiKap)[2]
  varKappa = diag(covDelPhiKap)[3]

## ----results = 'asis', echo = FALSE--------------------------------------
  print(xtable(covDelPhiKap, digits = 6), floating=FALSE, 
    hline.after=NULL, include.rownames=FALSE, include.colnames=FALSE, 
    tabular.environment = "array", only.contents = TRUE)

## ----echo = FALSE, include = FALSE---------------------------------------
  mu = delta
  sigma2 = varDelta
  a1 = (mu^2*(1 - mu) - mu*sigma2)/sigma2
  b1 = a1*(1 - mu)/mu
  aDelta = a1
  bDelta = b1
  mu = kappa1
  sigma2 = varKappa
  a1 = (mu^2*(1 - mu) - mu*sigma2)/sigma2
  b1 = a1*(1 - mu)/mu
  aKappa = a1
  bKappa = b1
  mu = phi
  sigma2 = varPhi
  a1 = (mu^2*(1 - mu) - mu*sigma2)/sigma2
  b1 = a1*(1 - mu)/mu
  aPhi = a1
  bPhi = b1

## ----echo = FALSE, include = FALSE---------------------------------------
  alpha = 1/eigen(M2x2)$values[1]
  eigValAlphaM2x2 = eigen(alpha*M2x2)$values[1]
  M2x2eigVal1 = alpha*M2x2
  delta2 = M2x2eigVal1[2,2]
  phi2 = M2x2eigVal1[1,2]
  kappa2 = M2x2eigVal1[2,1]
  varDelta2 = 9*varDelta
  varPhi2 = 9*varPhi
  varKappa2 = 9*varKappa
  mu = delta2
  sigma2 = varDelta2
  a1 = (mu^2*(1 - mu) - mu*sigma2)/sigma2
  b1 = a1*(1 - mu)/mu
  aDelta2 = a1
  bDelta2 = b1
  mu = kappa2
  sigma2 = varKappa2
  a1 = (mu^2*(1 - mu) - mu*sigma2)/sigma2
  b1 = a1*(1 - mu)/mu
  aKappa2 = a1
  bKappa2 = b1
  mu = phi2
  sigma2 = varPhi2
  a1 = (mu^2*(1 - mu) - mu*sigma2)/sigma2
  b1 = a1*(1 - mu)/mu
  aPhi2 = a1
  bPhi2 = b1
  #var uniform 0 to .5
  varUnif = .5^2/12
  scen3Mat = matrix(c(0, 0.25, 0.75, 0.75), nrow = 2, ncol = 2, byrow = TRUE)
  eigValScen3 = Re(eigen(scen3Mat)$values[1])

