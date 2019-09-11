
rm(list=ls())
setwd("C://Users/mkarnauskas/Desktop/red_snapper/Mar2016")

library(ncdf4)
library(maps)             # libraries
library(matlab)
library(grDevices)
library(sp)

nc <- nc_open('GEBCO_2014_2D_-84_24_-74_36.nc')        # open netcdf file from GEBCO
v1 <- nc$var[[1]]
z <-ncvar_get(nc, v1)                                
x <- v1$dim[[1]]$vals
y <- v1$dim[[2]]$vals
nc_close(nc)

z[which(z>0)] <- 0
image(x,y,z)

c15<- contourLines(x,y,z, levels=c(-15))                  # isobaths for 0, 15, 64m
  for (i in 1:1) {  lines(c15[[i]]$x, c15[[i]]$y)  }
c64 <- contourLines(x,y,z, levels=c(-64))
  for (i in 1:1) {  lines(c64[[i]]$x, c64[[i]]$y)  }
c0 <- contourLines(x,y,z, levels=c(0))
  for (i in 1:1) {  lines(c0[[i]]$x, c0[[i]]$y)  }
  
##### define shallow and deep polygon boundaries
sh <- cbind(c15[[1]]$x[seq(1,5650,5)], c15[[1]]$y[seq(1,5650,5)])
dp <- cbind(c64[[1]]$x[seq(1,4500,3)], c64[[1]]$y[seq(1,4500,3)])

########## connect with Easternmost polygon from Gulf grid
d <- read.table("C://Users/mkarnauskas/Desktop/red_snapper_CMS/redsnapperRadialGrid_STATE.xyz", sep="\t")
d1 <- d[which(d$V3==75),]
polygon(d1$V1, d1$V2, col=2)

#image(x,y,z, xlim=c(-83.5, -80), ylim=c(24, 26))
#points(sh[,1], sh[,2], pch=20, cex=0.2)
#points(dp[,1], dp[,2], pch=20, cex=0.2)

polygon(d1$V1, d1$V2, col=2)
points(d1[which.min(d1$V2),1:2], pch=20, col=5)    # 34
d1[which.min(d1$V2),1:2]
points(-82.55262, 24.54242, pch=20, col=5)         # 26
a <- which.min(abs(sh[,1] - (-82.55262)) + abs(sh[,2] - (24.54242)))
points(sh[a,1], sh[a,2], col=4, pch=20)

sh <- sh[1:a,]
dp <- dp[which(dp[,1] > (-83.21364)),]

image(x,y,z)
points(sh, cex=0.2)
points(dp, cex=0.2)
polygon(d1$V1, d1$V2, col=2)

##########################   start making file   ###############################

pols <- data.frame()   ##  polygon holder

##########  first polygon -- connect with boundaries of #75 from Gulf

points(d1[26:34,], pch=20)
a <- which(sh[,1] < (-81.7))
b <- which(dp[,1] < (-81.7))
b <- b[-c(1:3)]     # otherwise concave polygon

p <- data.frame(rbind(sh[a,], d1[26:34,1:2], dp[b,]))
polygon(p, col=1)
plot(p)
polygon(p, col=1)
p$id <- 76            # one more than polygon in Gulf                                                
image(x,y,z, xlim=c(-83.5, -80), ylim=c(24, 26))  # check
polygon(d1$V1, d1$V2, col=2)
polygon(p, col=1)

pols <- rbind(pols, p)    # add to holder variable

map("usa", xlim=c(-84, -74), ylim=c(24, 36))
#image(x,y,z, add=T)
polygon(d1$V1, d1$V2, col=2)

########## polygon # 77
a <- which(sh[,1] >= (-81.7) & sh[,1] < (-80.5) & sh[,2] < 25)
b <- which(dp[,1] >= (-81.7) & dp[,1] < (-80.5) & dp[,2] < 25)
p <- as.data.frame(rbind(sh[a,], dp[b,]), col.names=c("V1", "V2"))
polygon(p, col=3)
p$id <- 77
pols <- rbind(pols, p)

########## polygon # 78

a1 <- which.max(sh[,2] < 26)
b1 <- which.min(dp[,2] < 26)
p <- as.data.frame(rbind(sh[min(a):a1,], dp[b1:max(b),]), col.names=c("V1", "V2"))
polygon(p, col=4)
p$id <- 78
pols <- rbind(pols, p)

######## polygons 79 - 87                                                                       
latcuts <- c(26, 27, 28, seq(28.5, 30.0, 0.25))
pct <- 78

for (i in 2:length(latcuts))   {
  a <- which.max(sh[,2] < latcuts[i-1])
  b <- which.min(dp[,2] <= latcuts[i-1])
  a1 <- which.max(sh[,2] <= latcuts[i])
  b1 <- which.min(dp[,2] < latcuts[i])
  p <- as.data.frame(rbind(sh[a:a1,], dp[b1:b,]), col.names=c("V1", "V2"))
  polygon(p, col=i)
  text(mean(p$V1), mean(p$V2), pct, col="white", cex=0.8) 
  p$id <- pct + 1
  pct <- pct + 1
  pols <- rbind(pols, p)  }

####  polygon 88 - 107 -- corner pieces                                     
  a <- which.max(sh[,1] < (-78))
  b <- which.min(dp[,1] < (-78))
  a1 <- which.max(sh[,2] < latcuts[i])
  b1 <- which.min(dp[,2] < latcuts[i])
  
  aa <- (a:a1)[round(seq(length(a:a1), 1, length.out = 21))]
  bb <- (b:b1)[round(seq(length(b:b1), 1, length.out = 21))]
  aa[3] <- 545
  aa[13] <- 360 
  
  for (k in 1:20) { 
    p <- as.data.frame(rbind(sh[aa[k+1]:aa[k],], dp[bb[k]:bb[k+1],]), col.names=c("V1", "V2"))
    polygon(p, col=k+1)
    p$id <- pct + 1                               
    pct <- pct + 1
    text(mean(p$V1), mean(p$V2), pct, col="white", cex=0.8) 
    pols <- rbind(pols, p)   
     }  
  
#  windows()
#  f <- pols[which(pols$id==89),]
#  f <- pols[which(pols$id==99),]
#  plot(f$V1, f$V2, col=0)
#  polygon(f)
#  text(mean(f$V1), mean(f$V2), mean(f$id)) 

### polygons 108 ++
loncuts <- c(-78, -77.825, seq(-77.5, -75.75, 0.25))

for (i in 2:length(loncuts))   {
  a <- which.max(sh[,1] < loncuts[i-1])
  b <- which.min(dp[,1] < loncuts[i-1])
  a1 <- which.max(sh[,1] < loncuts[i])
  b1 <- which.min(dp[,1] < loncuts[i])
#  if (i==10) { a1 <- 38}
  p <- as.data.frame(rbind(sh[a:a1,], dp[b1:b,]), col.names=c("V1", "V2"))
  polygon(p, col=i)
  p$id <- pct + 1
  pct <- pct + 1
  pols <- rbind(pols, p)  }
  
########## last polygon
a <- which(sh[,1] >= (-75.75) & sh[,1] < (-75.3) & sh[,2] < 35.18)
b <- which(dp[,1] >= (-75.75) & dp[,1] < (-75.3) & dp[,2] < 35.18)
p <- as.data.frame(rbind(sh[a,], dp[b,]), col.names=c("V1", "V2"))
polygon(p, col=3)
p$id <- pct + 1
pct <- pct + 1
pols <- rbind(pols, p) 

  
#####################  end of polygon-making   ###############
map("usa", xlim=c(-84, -74), ylim=c(24, 36))
# image(x,y,z)  # check
polygon(d1$V1, d1$V2, col=2)
for (i in 76:117)  {
  f <- pols[which(pols$id==i),]
  polygon(f, col=i)  
  text(mean(f$V1), mean(f$V2), mean(f$id), cex=0.5, col="white")
  }

c0 <- contourLines(x,y,z, levels=c(0))
for (i in 1:1) {  lines(c0[[i]]$x, c0[[i]]$y)  }

windows()
for (i in 106:125)  {
  f <- pols[which(pols$id==i),]
  plot(f$V1, f$V2, col=0)
  polygon(f)
  text(mean(f$V1), mean(f$V2), i)
  Sys.sleep(2) }


############  concatenate with Gulf   ###########                  

names(pols) <- names(d)
pols <- rbind(d, pols)

write.table(pols, file="redSnapperSett_GOM_ATL_hires.txt", sep="\t", col.names=F, row.names=F)

########## check file
rm(list=ls())
d <- read.table("C://Users/mkarnauskas/Desktop/red_snapper/Mar2016/redSnapperSett_GOM_ATL_hires.txt", sep="\t")
map("state", ylim=c(24, 36), xlim=c(-100, -75))  # check
for (i in 1:117)  {
  f <- d[which(d$V3==i),]
  polygon(f, col=i)  }


     
              