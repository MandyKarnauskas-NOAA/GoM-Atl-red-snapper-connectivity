<<<<<<< HEAD

rm(list=ls())
library(ncdf4)
library(maps)

setwd("C://Users/mkarnauskas/Desktop/RS_FATEproject/Hatteras_issue/Matthieu_HYCOM20082009")
#setwd("C://Users/mkarnauskas/Desktop/RS_FATEproject/Hatteras_issue/CMS_inputs_outputs")
setwd("C://Users/mkarnauskas/Desktop/RS_FATEproject/Tortugas_issue/")

##############  concatenate files  ###################
#
#dat <- read.table("con_file_1_hatterasSABGOM")
dat1 <- read.table("con_file_1")
dat2 <- read.table("con_file_2")
dat <- rbind(dat1, dat2)
colnames(dat) <- c("rel_poly","ret_poly","ret_yr","ret_mo","ret_d","age","ret_dep","rel_yr","rel_mo","rel_d")
##########################################################

###########  input release file
#mat <- read.table("RS_ATL_releaseHatteras2010.txt", sep="\t", header=F)
#mat <- read.table("TortugasRel_RS.txt", sep="\t", header=F)
mat <- read.table("releaseHatteras200809.txt", sep="\t", header=F)
tapply(mat$V5, mat$V6, sum)

nrow(dat) / sum(mat$V5) * 100     #  0.7987663   < 1% of particles returned

head(dat)
table(dat$rel_poly)
table(dat$ret_poly)
table(dat$ret_yr)
table(dat$ret_mo)
table(dat$ret_dep)

########################  plot trajectories  #############################

#nc <- nc_open("traj_file_1_GLBhycom.nc", write=FALSE, readunlim=TRUE, verbose=FALSE)
#nc <- nc_open("traj_file_1_hatterasSABGOM.nc", write=FALSE, readunlim=TRUE, verbose=FALSE)                                                                        
nc <- nc_open("traj_file_1.nc", write=FALSE, readunlim=TRUE, verbose=FALSE)

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
dates <- ncvar_get(nc, v7)

v8 <- nc$var[[8]]
pol <- ncvar_get(nc, v8)

v9 <- nc$var[[9]]
dens <- ncvar_get(nc, v9)

v10 <- nc$var[[10]]
diam <- ncvar_get(nc, v10)

nc_close(nc)

hist(dens)
hist(diam)

table(cod)                  # 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data

hist(dep[2,which(cod==(-1))])
hist(dens[which(cod==(-1))])
hist(lon[1,which(cod==(-1))])
hist(lat[1,which(cod==(-1))])

table(dep[1,which(cod==(-1))])

matplot(-dep, type="l", col="#FF00FF10")
matplot(-dep[,seq(1,9000,100)], type="l", col="#FF00FF10")


############  get depth  ###############

b <- nc_open("C:/Users/mkarnauskas/Desktop/RS_FATEproject/srtm15plus_58ff_d12b_7306.nc", write=FALSE, readunlim=TRUE, verbose=FALSE)

v2 <- b$var[[1]]
z <- ncvar_get(b, v2)
x <- b$var[[1]]$dim[[1]]$vals
y <- b$var[[1]]$dim[[2]]$vals
z[which(z>0)] <- NA
image(x,y,z, col=rainbow(100, start=0.5, end=0.6)[100:1])

nc_close(b)

# 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data
 
#windows()
#image(x,y,z, col=cols, ylim=c(23,31), xlim=c(-93,-80), axes=F, xlab="", ylab=""); box()
#map('usa', ylim=c(32.5, 38.5), xlim=c(-78,-72.2), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
#image(x,y,z, col=rainbow(100, start=0.5, end=0.6)[100:1], add=T); box()
#i <- which(cod != (-2))
#i <- seq(1, length(cod), length.out=500)  # & cod == (-4))
#  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }                            
#  for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }   


windows()
#nn <- c("settled", "close to land", "left area", "can still move"); m <- 1                   
nn <- c("settled", "no oceanographic data", "left area", "can still move"); m <- 1   
par(mfrow=c(2,2), mex=0.5)
#for (k in c(-4, -2, -1, 0))  { 
for (k in c(-4, -5, -1, 0))  { 
map('usa', ylim=c(32.5, 38.5), xlim=c(-78,-72.8), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
#map('usa', ylim=c(32.5, 38.5), xlim=c(-78,-72.2), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
image(x,y,z, col=rainbow(100, start=0.5, end=0.6)[100:1], add=T); box()
mtext(side=3, line=1, nn[m]); m <- m+1
i <- which(cod == (k))
text(-77.1, 32.8, paste("n =", length(i)))
for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }  } 
                        
                                                                               
windows()
#image(x,y,z, col=cols, ylim=c(23,31), xlim=c(-93,-80), axes=F, xlab="", ylab=""); box()
map('usa', ylim=c(34, 36), xlim=c(-77.0,-74.5), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
image(x,y,z, col=rainbow(100, start=0.5, end=0.6)[100:1], add=T); box()
i <- which(cod == (-4))
mtext(side=3, line=2.2, "successfully settled", cex=1.2)
mtext(side=3, line=1.0, "HYCOM 1/50 deg - 2008", cex=1.2)
  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
  for (j in i) {  points(lon[1,i]-360, lat[1,i], col="#FFFF0005", pch=19, cex=0.65)   }
  for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }  
  legend("topleft", c("start", "end"), pch=19, cex=1, col=c(7,2), bg="gray")


map('usa', ylim=c(32.5, 38.5), xlim=c(-78,-72.8), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
#map('usa', ylim=c(32.5, 38.5), xlim=c(-78,-72.2), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
image(x,y,z, col=rainbow(100, start=0.5, end=0.6)[100:1], add=T); box()
mtext(side=3, line=2.2, "can still move", cex=1.2)
mtext(side=3, line=1.0, "HYCOM 1/50 deg - 2008", cex=1.2)
i <- which(cod == 0)
text(-77.1, 32.8, paste("n =", length(i)))
for (j in i[seq(1, length(i), length.out=300)]) {  points(lon[1,i]-360, lat[1,i], col="#FFFF0005", pch=19, cex=0.35)      }
for (j in i[seq(1, length(i), length.out=300)]) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
for (j in i[seq(1, length(i), length.out=300)]) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }   



map('usa', ylim=c(34, 36), xlim=c(-77.0,-74.5), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
image(x,y,z, col=rainbow(100, start=0.5, end=0.6)[100:1], add=T); box()
i <- which(cod == (-1))
mtext(side=3, line=2, "left area", cex=1.2)
#  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
  for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.3)   }  










###########  Tortugas issue  ###########

# 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data

##########  bathymetry  ##############

nc <- nc_open("C:/Users/mkarnauskas/Desktop/RS_FATEproject/GEBCO_GOM_SATL.nc")        # open netcdf file from GEBCO
v1 <- nc$var[[1]]
z <-ncvar_get(nc, v1)
x <- v1$dim[[1]]$vals
y <- v1$dim[[2]]$vals
nc_close(nc)
z[z>0] <- 100
z[z< (-800)] <- -800

cols <- rainbow(100, start=0.50, end=0.58)[100:1]
cols[100] <- "white"

windows()
par(mfrow=c(2, 2), mex=0.5)
map('usa', xlim=c(-85.0, -75), ylim=c(23.85,35.7), col=0)
image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1); axis(2, las=2)
i <- which(cod == (-4) & dates < 2454800)
mtext(side=3, line=3, "successfully settled", cex=1)
mtext(side=3, line=1.0, "HYCOM 1/50 deg - 2008", cex=1)
  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
  for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }  
  for (j in i) {  points(lon[1,i]-360, lat[1,i], col="#FFFF0005", pch=19, cex=0.65)   }
  legend("topleft", c("start", "end"), pch=19, cex=1, col=c(7,2), bg="gray")

map('usa', xlim=c(-85.0, -75), ylim=c(23.85,35.7), col=0)
image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1); axis(2, las=2)
i <- which(cod == (-4) & dates > 2454800)
mtext(side=3, line=3, "successfully settled", cex=1)
mtext(side=3, line=1.0, "HYCOM 1/50 deg - 2009", cex=1)
  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
  for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }  
  for (j in i) {  points(lon[1,i]-360, lat[1,i], col="#FFFF0005", pch=19, cex=0.65)   }
  legend("topleft", c("start", "end"), pch=19, cex=1, col=c(7,2), bg="gray")



map('usa', xlim=c(-85.0, -75), ylim=c(23.85,35.7), col=0)
image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1); axis(2, las=2)
i <- which(cod == (0) & dates < 2454800)
mtext(side=3, line=3, "can still move", cex=1)
mtext(side=3, line=1.0, "HYCOM 1/50 deg - 2008", cex=1)
  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
  for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }  
  for (j in i) {  points(lon[1,i]-360, lat[1,i], col="#FFFF0005", pch=19, cex=0.65)   }
  legend("topleft", c("start", "end"), pch=19, cex=1, col=c(7,2), bg="gray")

map('usa', xlim=c(-85.0, -75), ylim=c(23.85,35.7), col=0)
image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1); axis(2, las=2)
i <- which(cod == (0) & dates > 2454800)
mtext(side=3, line=3, "can still move", cex=1)
mtext(side=3, line=1.0, "HYCOM 1/50 deg - 2009", cex=1)
  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
  for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }  
  for (j in i) {  points(lon[1,i]-360, lat[1,i], col="#FFFF0005", pch=19, cex=0.65)   }
  legend("topleft", c("start", "end"), pch=19, cex=1, col=c(7,2), bg="gray")


=======

rm(list=ls())
library(ncdf4)
library(maps)

setwd("C://Users/mkarnauskas/Desktop/RS_FATEproject/Hatteras_issue/Matthieu_HYCOM20082009")
#setwd("C://Users/mkarnauskas/Desktop/RS_FATEproject/Hatteras_issue/CMS_inputs_outputs")
setwd("C://Users/mkarnauskas/Desktop/RS_FATEproject/Tortugas_issue/")

##############  concatenate files  ###################
#
#dat <- read.table("con_file_1_hatterasSABGOM")
dat1 <- read.table("con_file_1")
dat2 <- read.table("con_file_2")
dat <- rbind(dat1, dat2)
colnames(dat) <- c("rel_poly","ret_poly","ret_yr","ret_mo","ret_d","age","ret_dep","rel_yr","rel_mo","rel_d")
##########################################################

###########  input release file
#mat <- read.table("RS_ATL_releaseHatteras2010.txt", sep="\t", header=F)
#mat <- read.table("TortugasRel_RS.txt", sep="\t", header=F)
mat <- read.table("releaseHatteras200809.txt", sep="\t", header=F)
tapply(mat$V5, mat$V6, sum)

nrow(dat) / sum(mat$V5) * 100     #  0.7987663   < 1% of particles returned

head(dat)
table(dat$rel_poly)
table(dat$ret_poly)
table(dat$ret_yr)
table(dat$ret_mo)
table(dat$ret_dep)

########################  plot trajectories  #############################

#nc <- nc_open("traj_file_1_GLBhycom.nc", write=FALSE, readunlim=TRUE, verbose=FALSE)
#nc <- nc_open("traj_file_1_hatterasSABGOM.nc", write=FALSE, readunlim=TRUE, verbose=FALSE)                                                                        
nc <- nc_open("traj_file_1.nc", write=FALSE, readunlim=TRUE, verbose=FALSE)

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
dates <- ncvar_get(nc, v7)

v8 <- nc$var[[8]]
pol <- ncvar_get(nc, v8)

v9 <- nc$var[[9]]
dens <- ncvar_get(nc, v9)

v10 <- nc$var[[10]]
diam <- ncvar_get(nc, v10)

nc_close(nc)

hist(dens)
hist(diam)

table(cod)                  # 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data

hist(dep[2,which(cod==(-1))])
hist(dens[which(cod==(-1))])
hist(lon[1,which(cod==(-1))])
hist(lat[1,which(cod==(-1))])

table(dep[1,which(cod==(-1))])

matplot(-dep, type="l", col="#FF00FF10")
matplot(-dep[,seq(1,9000,100)], type="l", col="#FF00FF10")


############  get depth  ###############

b <- nc_open("C:/Users/mkarnauskas/Desktop/RS_FATEproject/srtm15plus_58ff_d12b_7306.nc", write=FALSE, readunlim=TRUE, verbose=FALSE)

v2 <- b$var[[1]]
z <- ncvar_get(b, v2)
x <- b$var[[1]]$dim[[1]]$vals
y <- b$var[[1]]$dim[[2]]$vals
z[which(z>0)] <- NA
image(x,y,z, col=rainbow(100, start=0.5, end=0.6)[100:1])

nc_close(b)

# 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data
 
#windows()
#image(x,y,z, col=cols, ylim=c(23,31), xlim=c(-93,-80), axes=F, xlab="", ylab=""); box()
#map('usa', ylim=c(32.5, 38.5), xlim=c(-78,-72.2), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
#image(x,y,z, col=rainbow(100, start=0.5, end=0.6)[100:1], add=T); box()
#i <- which(cod != (-2))
#i <- seq(1, length(cod), length.out=500)  # & cod == (-4))
#  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }                            
#  for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }   


windows()
#nn <- c("settled", "close to land", "left area", "can still move"); m <- 1                   
nn <- c("settled", "no oceanographic data", "left area", "can still move"); m <- 1   
par(mfrow=c(2,2), mex=0.5)
#for (k in c(-4, -2, -1, 0))  { 
for (k in c(-4, -5, -1, 0))  { 
map('usa', ylim=c(32.5, 38.5), xlim=c(-78,-72.8), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
#map('usa', ylim=c(32.5, 38.5), xlim=c(-78,-72.2), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
image(x,y,z, col=rainbow(100, start=0.5, end=0.6)[100:1], add=T); box()
mtext(side=3, line=1, nn[m]); m <- m+1
i <- which(cod == (k))
text(-77.1, 32.8, paste("n =", length(i)))
for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }  } 
                        
                                                                               
windows()
#image(x,y,z, col=cols, ylim=c(23,31), xlim=c(-93,-80), axes=F, xlab="", ylab=""); box()
map('usa', ylim=c(34, 36), xlim=c(-77.0,-74.5), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
image(x,y,z, col=rainbow(100, start=0.5, end=0.6)[100:1], add=T); box()
i <- which(cod == (-4))
mtext(side=3, line=2.2, "successfully settled", cex=1.2)
mtext(side=3, line=1.0, "HYCOM 1/50 deg - 2008", cex=1.2)
  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
  for (j in i) {  points(lon[1,i]-360, lat[1,i], col="#FFFF0005", pch=19, cex=0.65)   }
  for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }  
  legend("topleft", c("start", "end"), pch=19, cex=1, col=c(7,2), bg="gray")


map('usa', ylim=c(32.5, 38.5), xlim=c(-78,-72.8), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
#map('usa', ylim=c(32.5, 38.5), xlim=c(-78,-72.2), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
image(x,y,z, col=rainbow(100, start=0.5, end=0.6)[100:1], add=T); box()
mtext(side=3, line=2.2, "can still move", cex=1.2)
mtext(side=3, line=1.0, "HYCOM 1/50 deg - 2008", cex=1.2)
i <- which(cod == 0)
text(-77.1, 32.8, paste("n =", length(i)))
for (j in i[seq(1, length(i), length.out=300)]) {  points(lon[1,i]-360, lat[1,i], col="#FFFF0005", pch=19, cex=0.35)      }
for (j in i[seq(1, length(i), length.out=300)]) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
for (j in i[seq(1, length(i), length.out=300)]) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }   



map('usa', ylim=c(34, 36), xlim=c(-77.0,-74.5), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
image(x,y,z, col=rainbow(100, start=0.5, end=0.6)[100:1], add=T); box()
i <- which(cod == (-1))
mtext(side=3, line=2, "left area", cex=1.2)
#  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
  for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.3)   }  










###########  Tortugas issue  ###########

# 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data

##########  bathymetry  ##############

nc <- nc_open("C:/Users/mkarnauskas/Desktop/RS_FATEproject/GEBCO_GOM_SATL.nc")        # open netcdf file from GEBCO
v1 <- nc$var[[1]]
z <-ncvar_get(nc, v1)
x <- v1$dim[[1]]$vals
y <- v1$dim[[2]]$vals
nc_close(nc)
z[z>0] <- 100
z[z< (-800)] <- -800

cols <- rainbow(100, start=0.50, end=0.58)[100:1]
cols[100] <- "white"

windows()
par(mfrow=c(2, 2), mex=0.5)
map('usa', xlim=c(-85.0, -75), ylim=c(23.85,35.7), col=0)
image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1); axis(2, las=2)
i <- which(cod == (-4) & dates < 2454800)
mtext(side=3, line=3, "successfully settled", cex=1)
mtext(side=3, line=1.0, "HYCOM 1/50 deg - 2008", cex=1)
  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
  for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }  
  for (j in i) {  points(lon[1,i]-360, lat[1,i], col="#FFFF0005", pch=19, cex=0.65)   }
  legend("topleft", c("start", "end"), pch=19, cex=1, col=c(7,2), bg="gray")

map('usa', xlim=c(-85.0, -75), ylim=c(23.85,35.7), col=0)
image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1); axis(2, las=2)
i <- which(cod == (-4) & dates > 2454800)
mtext(side=3, line=3, "successfully settled", cex=1)
mtext(side=3, line=1.0, "HYCOM 1/50 deg - 2009", cex=1)
  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
  for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }  
  for (j in i) {  points(lon[1,i]-360, lat[1,i], col="#FFFF0005", pch=19, cex=0.65)   }
  legend("topleft", c("start", "end"), pch=19, cex=1, col=c(7,2), bg="gray")



map('usa', xlim=c(-85.0, -75), ylim=c(23.85,35.7), col=0)
image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1); axis(2, las=2)
i <- which(cod == (0) & dates < 2454800)
mtext(side=3, line=3, "can still move", cex=1)
mtext(side=3, line=1.0, "HYCOM 1/50 deg - 2008", cex=1)
  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
  for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }  
  for (j in i) {  points(lon[1,i]-360, lat[1,i], col="#FFFF0005", pch=19, cex=0.65)   }
  legend("topleft", c("start", "end"), pch=19, cex=1, col=c(7,2), bg="gray")

map('usa', xlim=c(-85.0, -75), ylim=c(23.85,35.7), col=0)
image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1); axis(2, las=2)
i <- which(cod == (0) & dates > 2454800)
mtext(side=3, line=3, "can still move", cex=1)
mtext(side=3, line=1.0, "HYCOM 1/50 deg - 2009", cex=1)
  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   }
  for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }  
  for (j in i) {  points(lon[1,i]-360, lat[1,i], col="#FFFF0005", pch=19, cex=0.65)   }
  legend("topleft", c("start", "end"), pch=19, cex=1, col=c(7,2), bg="gray")


>>>>>>> 9f03132205832f02064f087874eb662c85e0458b
