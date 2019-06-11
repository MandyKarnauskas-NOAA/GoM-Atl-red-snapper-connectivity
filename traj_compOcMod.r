

rm(list=ls())

#############################  get bathymetry  #################################

nc <- nc_open("C:/Users/mkarnauskas/Desktop/RS_FATEproject/GEBCO_GOM_SATL.nc")        # open netcdf file from GEBCO
v1 <- nc$var[[1]]
z <-ncvar_get(nc, v1)
x <- v1$dim[[1]]$vals
y <- v1$dim[[2]]$vals
nc_close(nc)
z[z>0] <- 100
z[z< (-800)] <- -800

cols <- rainbow(100, start=0.55, end=0.62)[100:1]
cols[100] <- "white"

#########################   GET HYCOM TRAJECTORIES  ############################

setwd("C://Users/mkarnauskas/Desktop/RS_FATEproject/compOcMod/HYCOM")

lona <- c()
lata <- c()
depa <- c()
coda <- c()
dat1 <- c()
for (m in 1:7) {
  nam <- paste("traj_file_", m, ".nc", sep="")
nc <- nc_open(nam, write=FALSE, readunlim=TRUE, verbose=FALSE)
  v2 <- nc$var[[2]]
  lon <- ncvar_get(nc, v2)
  lon[which(lon>10000)] <- NA
  v3 <- nc$var[[3]]
  lat <- ncvar_get(nc, v3)
  lat[which(lat>10000)] <- NA
  v4 <- nc$var[[4]]
  dep <- ncvar_get(nc, v4)
  dep[which(dep>10000)] <- NA
  v6 <- nc$var[[6]]
  cod <- ncvar_get(nc, v6)
  cod[which(cod>10000)] <- NA    # 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data
  v7 <- nc$var[[7]]
  dat <- ncvar_get(nc, v7)
    lona <- cbind(lona, lon)
    lata <- cbind(lata, lat)
    depa <- cbind(depa, dep)
    coda <- c(coda, cod)
    dat1 <- c(dat1, dat)
nc_close(nc)                   }
brks <- seq(2453000-365, 2453000+365*12, 365)
plot(dat1)
dcut1 <- as.numeric(cut(dat1, brks))

##########################   GET SABGOM TRAJECTORIES  ##########################
setwd("C://Users/mkarnauskas/Desktop/RS_FATEproject/compOcMod/SABGOM")
setwd("C://Users/mkarnauskas/Desktop/RS_FATEproject/initial_run2")

lonb <- c()
latb <- c()
depb <- c()
codb <- c()
dat2 <- c()
intb <- c()
for (m in 1:6) {
  nam <- paste("traj_file_", m, ".nc", sep="")
nc <- nc_open(nam, write=FALSE, readunlim=TRUE, verbose=FALSE)
  v1 <- nc$var[[1]]
  int <- ncvar_get(nc, v1)
  v2 <- nc$var[[2]]
  lon <- ncvar_get(nc, v2)
  lon[which(lon>10000)] <- NA
  v3 <- nc$var[[3]]
  lat <- ncvar_get(nc, v3)
  lat[which(lat>10000)] <- NA
  v4 <- nc$var[[4]]
  dep <- ncvar_get(nc, v4)
  dep[which(dep>10000)] <- NA
  v6 <- nc$var[[6]]
  cod <- ncvar_get(nc, v6)
  cod[which(cod>10000)] <- NA    # 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data
  v7 <- nc$var[[7]]
  dat <- ncvar_get(nc, v7)
    lonb <- cbind(lonb, lon)
    latb <- cbind(latb, lat)
    depb <- cbind(depb, dep)
    codb <- c(codb, cod)
    dat2 <- c(dat2, dat)
    intb <- c(intb, int)
  nc_close(nc)                   }
plot(dat2)
dcut2 <- as.numeric(cut(dat2, brks))


############################   PLOT TRAJECTORIES   #############################

cols <- gray(seq(0.6, 0.999, length.out=100))
cols[100] <- "black"
#cols[1:99] <- "white"
cols2 <- rainbow(31, start=0.01, end=0.8, alpha=0.20)

par(mfrow=c(3, 2), mex=0.4)

yrs <- 2003:2010
m <- 2
for (m in 2:7) {

map('usa', xlim=c(-98, -85), ylim=c(23.8,31), col=0)
image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1, cex=1.2); axis(2, las=2, cex=1.2)
legend(-98.7, 31.2, paste(yrs[m]), cex=1.5, text.font=2, bty="n", text.col="white")

  i <- which(dcut2 == m) #  & coda != (-4))
    k <- i[seq(1, length(i), length.out=200)]
  for (j in k) {
    lines(lonb[,j]-360, latb[,j], col="#FF000040")        # SABGOM
    points(lona[1,j]-360, lata[1,j], col="#00000060", pch=19, cex=0.8)   }

  i <- which(dcut1 == m) #  & coda != (-4))
    k <- i[seq(1, length(i), length.out=200)]
  for (j in k) {
    lines(lona[,j]-360, lata[,j], col="#0000FF15")       # HYCOM
    points(lona[1,j]-360, lata[1,j], col="#00000060", pch=19, cex=0.8)  }
    legend("bottomright", c("HYCOM", "SABGOM"), col=c("#0000FF80", "#FF000080"), lty=1, lwd=2, bty="n")
    }
    
cols2 <- rainbow(31, start=0.01, end=0.8, alpha=0.20)
xloc <- seq(-92,-86, length.out=31)
cols3 <- rainbow(31, start=0.01, end=0.8, alpha=0.80)

######################  PLOT TRAJECTORIES BY MODEL   ###########################
### SABGOM
par(mfrow=c(3, 2), mex=0.4)
for (m in 2:7) {
  map('usa', xlim=c(-98, -85), ylim=c(23.8,31), col=0)
  image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1, cex=1.2); axis(2, las=2, cex=1.2)
  legend(-98.7, 31.2, paste("SABGOM", yrs[m]), cex=1.2, text.font=2, bty="n", text.col="white")
    i <- which(dcut2 == m) #  & coda != (-4))
      k <- i[seq(1, length(i), length.out=500)]
  for (j in k) {
    points(lonb[1,j]-360, latb[1,j], col="#00000060", pch=19, cex=0.7)
    for (f in 2:length(lona[,1])) {
    lines(lonb[(f-1):f,j]-360, latb[(f-1):f,j], col=cols2[f]) }    }
 for (j in 1:31) {   polygon(c(xloc[j], xloc[j+1],xloc[j+1], xloc[j]), c(24.35, 24.35, 24.7, 24.7), col=cols3[j], border=NA) }
 text(-89, 25.1, "days after release", cex=1)
 text(xloc[seq(4,30,4)], 24.1, seq(4,30,4), cex=1)
    }

### HYCOM
par(mfrow=c(3, 2), mex=0.4)
for (m in 2:7) {
  map('usa', xlim=c(-98, -85), ylim=c(23.8,31), col=0)
  image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1, cex=1.2); axis(2, las=2, cex=1.2)
  legend(-98.7, 31.2, paste("HYCOM", yrs[m]), cex=1.2, text.font=2, bty="n", text.col="white")
    i <- which(dcut1 == m) #  & coda != (-4))
      k <- i[seq(1, length(i), length.out=500)]
  for (j in k) {
    points(lona[1,j]-360, lata[1,j], col="#00000060", pch=19, cex=0.7)
    for (f in 2:length(lona[,1])) {
    lines(lona[(f-1):f,j]-360, lata[(f-1):f,j], col=cols2[f]) }    }
 for (j in 1:31) {   polygon(c(xloc[j], xloc[j+1],xloc[j+1], xloc[j]), c(24.35, 24.35, 24.7, 24.7), col=cols3[j], border=NA) }
 text(-89, 25.1, "days after release", cex=1)
 text(xloc[seq(4,30,4)], 24.1, seq(4,30,4), cex=1)
    }

################################################################################

rm(list=ls())
setwd("C://Users/mkarnauskas/Desktop/RS_FATEproject/compOcMod")

nc <- nc_open("C:/Users/mkarnauskas/Desktop/RS_FATEproject/GEBCO_GOM_SATL.nc")        # open netcdf file from GEBCO
v1 <- nc$var[[1]]
z <-ncvar_get(nc, v1)
x <- v1$dim[[1]]$vals
y <- v1$dim[[2]]$vals
nc_close(nc)
z[z>0] <- 100
z[z< (-800)] <- -800
cols <- rainbow(100, start=0.55, end=0.62)[100:1]
cols[100] <- "white"

nam <- "traj_file_diff5.nc"
nc <- nc_open(nam, write=FALSE, readunlim=TRUE, verbose=FALSE)
  v1 <- nc$var[[1]]
  int <- ncvar_get(nc, v1)
  v2 <- nc$var[[2]]
  lon <- ncvar_get(nc, v2)
  lon[which(lon>10000)] <- NA
  v3 <- nc$var[[3]]
  lat <- ncvar_get(nc, v3)
  lat[which(lat>10000)] <- NA
  v4 <- nc$var[[4]]
  dep <- ncvar_get(nc, v4)
  dep[which(dep>10000)] <- NA
  v6 <- nc$var[[6]]
  cod <- ncvar_get(nc, v6)
  cod[which(cod>10000)] <- NA    # 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data
  v7 <- nc$var[[7]]
  dat <- ncvar_get(nc, v7)
  nc_close(nc)              
brks <- seq(2453000-365, 2453000+365*12, 365)
plot(dat)
dcut1 <- as.numeric(cut(dat, brks))

nam <- "traj_file_diff25.nc"
nc <- nc_open(nam, write=FALSE, readunlim=TRUE, verbose=FALSE)
  v1 <- nc$var[[1]]
  inta <- ncvar_get(nc, v1)
  v2 <- nc$var[[2]]
  lona <- ncvar_get(nc, v2)
  lona[which(lona>10000)] <- NA
  v3 <- nc$var[[3]]
  lata <- ncvar_get(nc, v3)
  lata[which(lata>10000)] <- NA
  v4 <- nc$var[[4]]
  depa <- ncvar_get(nc, v4)
  depa[which(depa>10000)] <- NA
  v6 <- nc$var[[6]]
  coda <- ncvar_get(nc, v6)
  coda[which(coda>10000)] <- NA    # 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data
  v7 <- nc$var[[7]]
  dat2 <- ncvar_get(nc, v7)
  nc_close(nc)                  
plot(dat2)
dcut2 <- as.numeric(cut(dat2, brks))

cols <- gray(seq(0.6, 0.999, length.out=100))
cols[100] <- "black"
cols2 <- rainbow(31, start=0.01, end=0.8, alpha=0.20)
xloc <- seq(-92,-86, length.out=31)
cols3 <- rainbow(31, start=0.01, end=0.8, alpha=0.80)
yrs <- 2003:2010

### diff 25
par(mfrow=c(5, 2), mex=0.4)
for (m in 2:6) {
  map('usa', xlim=c(-98, -85), ylim=c(23.8,31), col=0)
  image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1, cex=1.2); axis(2, las=2, cex=1.2)
  legend(-98.7, 31.2, paste(yrs[m], "diffusivity = 5"), cex=1.2, text.font=2, bty="n", text.col="white")
    i <- which(dcut1 == m) #  & cod != (-4))
      k <- i[seq(1, length(i), length.out=500)]
  for (j in k) {
    points(lon[1,j]-360, lat[1,j], col="#00000060", pch=19, cex=0.7)
    for (f in 2:length(lon[,1])) {
    lines(lon[(f-1):f,j]-360, lat[(f-1):f,j], col=cols2[f]) }    }
 for (j in 1:31) {   polygon(c(xloc[j], xloc[j+1],xloc[j+1], xloc[j]), c(24.35, 24.35, 24.7, 24.7), col=cols3[j], border=NA) }
 text(-89, 25.1, "days after release", cex=1)
 text(xloc[seq(4,30,4)], 24.1, seq(4,30,4), cex=1)

  map('usa', xlim=c(-98, -85), ylim=c(23.8,31), col=0)
  image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1, cex=1.2); axis(2, las=2, cex=1.2)
  legend(-98.7, 31.2, paste(yrs[m], "diffusivity = 25"), cex=1.2, text.font=2, bty="n", text.col="white")
    i <- which(dcut2 == m) #  & coda != (-4))
      k <- i[seq(1, length(i), length.out=500)]
  for (j in k) {
    points(lona[1,j]-360, lata[1,j], col="#00000060", pch=19, cex=0.7)
    for (f in 2:length(lona[,1])) {
    lines(lona[(f-1):f,j]-360, lata[(f-1):f,j], col=cols2[f]) }    }
 for (j in 1:31) {   polygon(c(xloc[j], xloc[j+1],xloc[j+1], xloc[j]), c(24.35, 24.35, 24.7, 24.7), col=cols3[j], border=NA) }
 text(-89, 25.1, "days after release", cex=1)
 text(xloc[seq(4,30,4)], 24.1, seq(4,30,4), cex=1)
    }










setwd("C://Users/mkarnauskas/Desktop/RS_FATEproject/GoM_Atl_conn_lores_byyear")

lonb <- c()
latb <- c()
depb <- c()
codb <- c()
dat2 <- c()
intb <- c()

for (m in 1:7) {
  nam <- paste("traj_file_", m, ".nc", sep="")
nc <- nc_open(nam, write=FALSE, readunlim=TRUE, verbose=FALSE)
  v1 <- nc$var[[1]]
  int <- ncvar_get(nc, v1)
  v2 <- nc$var[[2]]
  lon <- ncvar_get(nc, v2)
  lon[which(lon>10000)] <- NA
  v3 <- nc$var[[3]]
  lat <- ncvar_get(nc, v3)
  lat[which(lat>10000)] <- NA
  v4 <- nc$var[[4]]
  dep <- ncvar_get(nc, v4)
  dep[which(dep>10000)] <- NA
  v6 <- nc$var[[6]]
  cod <- ncvar_get(nc, v6)
  cod[which(cod>10000)] <- NA    # 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data
  v7 <- nc$var[[7]]
  dat <- ncvar_get(nc, v7)
    lonb <- cbind(lonb, lon)
    latb <- cbind(latb, lat)
    depb <- cbind(depb, dep)
    codb <- c(codb, cod)
    dat2 <- c(dat2, dat)
    intb <- c(intb, int)
  nc_close(nc)                   }

plot(dat2)
dcut2 <- as.numeric(cut(dat2, brks))


map('usa', xlim=c(-82, -78), ylim=c(27,33), col=0)
  image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1, cex=1.2); axis(2, las=2, cex=1.2)
  legend(-98.7, 31.2, paste("SABGOM", yrs[m]), cex=1.2, text.font=2, bty="n", text.col="white")
    i <- which(dcut2 == 2) #  & coda != (-4))
      k <- i[seq(1, length(i), length.out=500)]
  for (j in k) {
  if (lonb[1,j]-360 > (-81.5)) {
    points(lonb[1,j]-360, latb[1,j], col="#00000060", pch=19, cex=0.7)
    for (f in 2:length(lona[,1])) {
    lines(lonb[(f-1):f,j]-360, latb[(f-1):f,j], col=cols2[f]) }    } } 
    
xloc <- seq(-81,-79, length.out=31)
 for (j in 1:31) {   polygon(c(xloc[j], xloc[j+1],xloc[j+1], xloc[j]), c(27.2, 27.2, 27.3, 27.3), col=cols3[j], border=NA) }
 text(-80, 27.4, "days after release", cex=1)
 text(xloc[seq(4,30,4)], 27.1, seq(4,30,4), cex=1)
 
 mtext(side=3, line=2, "SABGOM trajectories", cex=1.1, font=2)
