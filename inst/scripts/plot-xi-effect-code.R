<<plot-xieffect, fig.width=12, fig.height=12, echo= FALSE, include = FALSE, cache = TRUE, dev='jpeg'>>=
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


@
  