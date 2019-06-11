<<<<<<< HEAD
################################################################################
###########################   red snapper release file for CMS  ################
###########################   M. Karnauskas Jan 8, 2015         ################
#
#  code takes output from GLM model and krigging (C:/Karnauskas/Desktop/RS_mapping_with_platforms)
#  includes final fecundity map from published version
################################################################################

rm(list=ls())

if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4", repos='http://cran.us.r-project.org')
if (!"sp" %in% installed.packages()) install.packages("sp", repos='http://cran.us.r-project.org')
if (!"maps" %in% installed.packages()) install.packages("maps", repos='http://cran.us.r-project.org')
library(ncdf4)
library(sp)
library(maps)

source("C:/Users/mkarnauskas/Desktop/RS_FATEproject/plotting.r")                                                                                                  
load("C:/Users/mkarnauskas/Desktop/RS_FATEproject/FINAL_MAPPING_RESULTS.RData")        #  Fwith_platform is index to be used
head(mat)

#plotonmap(mat$Fwith_resid, mat$statlon, mat$statlat, 15, 0.9, addplus=T); text(-90.5,26.25, "index of fecundity")
plotonmap(mat$Ftotal, mat$statlon, mat$statlat, 15, 0.9, addplus=T); text(-90.5,26.25, "index of fecundity")
#plotonmaphires(mat$Ftotal, mat$statlon, mat$statlat, 15, 0.9, addplus=T, adjlab=0.4)      

dat <- data.frame(cbind(mat$statlon, mat$statlat, mat$Ftotal))
names(dat) <- c("lon", "lat", "N")

################################################################################                                                    
co <- read.table("C:/Users/mkarnauskas/Desktop/RS_FATEproject/MASTER_codes/CMS_input_files/redSnapperSett_GOM_ATL.xyz", header=F)

################  plot original recruitment habitat grid cells  ################
plot(1, xlim=c(-98,-81), ylim=c(24,31))
for (j in unique(co[,3]))  {
  m <- co[which(co[,3]==j),]; polygon(m, col=j)
  text(mean(m[,1]), mean(m[,2]), m[1,3], cex=0.6)      }
points(dat$lon, dat$lat, pch=19, cex=0.5)             # release locations

##########  get polygon numbers for release locations  ##############
pts = SpatialPoints(cbind(dat$lon, dat$lat))
m <- co[which(co[,3]==1), ]
m <- rbind(m, m[1,])
oldpol <- Polygons(list(Polygon(cbind(m[,1],m[,2]))), ID=1)
for (j in 2:max(co$V3))  {
  m <- co[which(co[,3]==j), ]
  m <- rbind(m, m[1,])
  newpol <- Polygons(list(Polygon(cbind(m[,1],m[,2]))), ID=j)
  oldpol <- c(oldpol, newpol)           }
pol = SpatialPolygons(oldpol)
polynams <- over(pts, pol)
rel <- cbind(polynams, dat)
rel
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams/rel$polynams+1)

############ find center point of each polygon  ######
ctrslat <- tapply(co$V2, co$V3, mean)
ctrslon <- tapply(co$V1, co$V3, mean)
points(ctrslon, ctrslat, pch=20, cex=2)

##############  label release site with closest polygon number  ############
for (i in 1:nrow(rel)) {
  if (is.na(rel$polynams[i])) {
     rel$polynams[i] <- which.min(abs(rel$lon[i]-ctrslon) + abs(rel$lat[i]-ctrslat))  } }
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams/(rel$polynams+1))

##########  plot new grid cells  #########################
plot(1, xlim=c(-98,-81), ylim=c(24,31))
for (j in unique(co[,3]))  {
  m <- co[which(co[,3]==j),]; polygon(m, lwd=1)
  text(m[1,1], m[1,2], m[1,3], cex=0.6)         }
map('usa', add=T)  
cols <- rainbow(102)  
points(rel$lon, rel$lat, pch=19, cex=0.5, col=cols[rel$polynams])    #  check calcs  
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams)

###############################     ADD DEPTH      ##################################################
###  !!!! NOTE !!!! 
###  Important to use nest from simulation to be run to avoid particles being trapped below surface 
#
nc <- nc_open("C:/Users/mkarnauskas/Desktop/RS_FATEproject/nest_1_20080501000000_HYCOM150.nc")    
#nc <- nc_open("C:/Users/mkarnauskas/Desktop/RS_FATEproject/nest_1_20070701000000.nc")          # open netcdf and get temp variable

v1 <- nc$var[[1]]
u <- ncvar_get(nc, v1)
dim(u)
v1 <- nc$var[[2]]
v <- ncvar_get(nc, v1)
dim(v)       
nc_close(nc)
lon <- nc$var[[1]]$dim[[1]]$vals - 360
lat <- nc$var[[1]]$dim[[2]]$vals
dep <- nc$var[[1]]$dim[[3]]$vals
cur <- sqrt(u^2 + v^2)

image(lon, lat, cur[,,1])
image(lon, lat, cur[,,1], xlim=c(-100, -80), ylim=c(24,31))

rel$depest <- NA
for (i in 1:nrow(rel)) {  rel$depest[i] <- dep[max(which(!is.na(u[which.min(abs(lon - rel$lon[i])), which.min(abs(lat - rel$lat[i])),])))] }
head(rel)

points(rel$lon, rel$lat, pch=19, col=is.na(rel$depest))
rel$depest[is.na(rel$depest)] <- 0

table(rel$depest)
hist(rel$depest)

length(which(rel$depest<45))
length(which(rel$depest>=45))

rel$spawndep <- rel$depest - 5
rel$spawndep[which(rel$spawndep > 45)] <- 45

minspawndep <- 13

length(which(rel$depest<=minspawndep))
dim(rel)
tapply(rel$N, rel$depest<=minspawndep, sum)
rel <- rel[which(rel$depest>minspawndep),]
dim(rel)

plot(rel$depest, rel$spawndep)
table(rel$spawndep < rel$depest)
table(rel$spawndep)
table(rel$depest, rel$spawndep)

plotonmap(rel$spawndep, rel$lon, rel$lat, 15, 0.6)
cols <- rainbow(102)  
points(rel$lon, rel$lat, pch=19, cex=1.5, col=cols[rel$spawndep])

################################################################################
###########   only run this section after initial Gulf-wide run   ##############
###########   cut out release areas which do not supply S Atl     ##############

# polygon 62 is first to have successful recruitment to S Atl

#rel <- rel[which(rel$polynams > 55 & rel$lat < 28.3),]        # uncomment here to apply cutoff 

plotonmap(rel$spawndep, rel$lon, rel$lat, 15, 0.6)

plot(1, xlim=c(-98,-81), ylim=c(24,31))
map('usa', add=T)  
cols <- rainbow(102)  
points(rel$lon, rel$lat, pch=19, cex=0.5, col=cols[rel$polynams])
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams)
              
###################  MAKE RELEASE FILE  #####################################
######################### input temporal information ###########################

days <- read.table("C:/Users/mkarnauskas/Desktop/RS_FATEproject/RELEASEFILE_INDEXnew.csv", sep=",", header=T)     # this inputs a file which references lunar phases and peak spawning times with all dates from 2004 - 2008
#days$red_snapper[which(days$mo==4)] <- 1
sp <- which(days$red_snapper==1 & days$yr==2008)          # adjust year as necessary to line up with simulation year
sp <- sp[seq(1, length(sp)-1, 6)]                         # every 6 days or every 3 days
#sp <- sp[seq(1, length(sp)-1, 3)]
lis <- days[sp, 1:3]
dim(lis)
lis <- lis[-1,]                                          # so that even number in each year
lis <- lis[which(lis$mo<9),]

table(lis$yr)
table(lis$yr, lis$mo)
dim(lis)
lis2 <- lis

for (i in 2009:2009) {                                  # adjust for multiple years as necessary to match simulation years
   lis2$yr <- i
   lis <- rbind(lis, lis2) }

dim(lis)
table(lis$yr, lis$mo)
table(lis$yr)
mean(table(lis$yr))

################  spawning seasonality from Porch et al. 2015  #################
lis$doy <- NA
dinmon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
for (i in 1:nrow(lis))  {  lis$doy[i] <- (sum(dinmon[1:lis$mo[i]]) + lis$da[i]) / 365  }
lis$spawnact <- (lis$doy/0.536)^(0.536/0.024) * exp((0.536-lis$doy)/0.024)
plot(lis$doy, lis$spawnact)
plot(lis$spawnact)
plot(lis$spawnact[1:70])

lis <- lis[-c(4)]

###############  input spatial site information (from above)  ###################

rel1 <- rel[-c(5)]
head(rel1)
m <- rel1       # 'm' is list of release sites with columns: polygon, lon, lat, number of releases
head(m)

m$N <- m$N * 2 #  - SCALE AS NECESSARY           
mean(m$N); min(m$N); max(m$N)
which(m$N==0)

prod(nrow(lis), nrow(m))
###################### now, making the release file ############################

mat <- as.data.frame(matrix(data=NA, nrow=nrow(lis)*nrow(m), ncol=9))   # empty matrix to be filled
mat[,1] <- rep(m[,1], nrow(lis))                                        # column 1: release polygon number
mat[,2] <- rep(m[,2], nrow(lis))                                        # column 2: release longitude 
mat[,3] <- rep(m[,3], nrow(lis))                                        # column 3: release latitude 
mat[,4] <- rep(m[,5], nrow(lis))                                        # column 4: release depth
mat[,5] <- rep(m[,4], nrow(lis))                                        # column 5: number of particles per release
mat <- mat[order(mat[,1], mat[,5], mat[,3], mat[,2]), ] # !!!!  CHECK THIS WITH NEW DATA   !!!!                                # resort matrix
# mat
mat[,6] <- rep(lis[,1], nrow(m))                                        # column 6: release year
mat[,7] <- rep(lis[,2], nrow(m))                                        # column 7: release month
mat[,8] <- rep(lis[,3], nrow(m))                                        # column 8: release day
mat[,9] <- 0                                                            # column 9: release hour

mat[,10] <- rep(lis[,4], nrow(m))                                       # column 10: scale by spawning activity!!!!
sum(mat[,5])

mat <- mat[order(mat[,6], mat[,7], mat[,8]), ] # !!!!  REORDER SO DATES ARE TOGETHER   !!!!
head(mat)

mat$V5 <- round(mat$V5 * mat$V10)
sum(mat[,5])
head(mat)

matfin <- mat[-c(10)]
head(matfin)
dim(matfin)
table((matfin$V5 >0))
table((matfin$V5 >0))[1] / nrow(matfin)     # this is % lost to rounding down to zero
matfin <- matfin[which(matfin$V5 >0),]
head(matfin)
dim(matfin)
mean(matfin$V5); min(matfin$V5); max(matfin$V5)
sum(matfin$V5)

getwd()
dim(matfin)
dim(matfin)/8

matfinGOM <- matfin
save(matfinGOM, file="C:/Users/mkarnauskas/Desktop/RS_FATEproject/fullGoM_release_HYCOM150.RData")

# double check
table(matfin$V6, matfin$V7)
table(matfin$V6)
matplot(table(matfin$V7, matfin$V6), type="l")
diff(table(matfin$V6))

tapply(matfin$V5, list(matfin$V6, matfin$V7), sum)
matplot(tapply(matfin$V5, list(matfin$V7, matfin$V6), sum), type="l")

f <- which(matfin$V6==2008 & matfin$V7 == 8 & matfin$V8 == 29); length(f)
plotonmap(matfin$V5[f], matfin$V2[f], matfin$V3[f], cexnum=0.6, pchnum=15)

f <- which(matfin$V6==2008 & matfin$V7 == 6 & matfin$V8 == 24); length(f)
plotonmap(matfin$V5[f], matfin$V2[f], matfin$V3[f], cexnum=0.6, pchnum=15)

####################   END GOM RELEASE LOCATIONS   #############################
################################################################################
=======
################################################################################
###########################   red snapper release file for CMS  ################
###########################   M. Karnauskas Jan 8, 2015         ################
#
#  code takes output from GLM model and krigging (C:/Karnauskas/Desktop/RS_mapping_with_platforms)
#  includes final fecundity map from published version
################################################################################

rm(list=ls())

if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4", repos='http://cran.us.r-project.org')
if (!"sp" %in% installed.packages()) install.packages("sp", repos='http://cran.us.r-project.org')
if (!"maps" %in% installed.packages()) install.packages("maps", repos='http://cran.us.r-project.org')
library(ncdf4)
library(sp)
library(maps)

source("C:/Users/mkarnauskas/Desktop/RS_FATEproject/plotting.r")                                                                                                  
load("C:/Users/mkarnauskas/Desktop/RS_FATEproject/FINAL_MAPPING_RESULTS.RData")        #  Fwith_platform is index to be used
head(mat)

#plotonmap(mat$Fwith_resid, mat$statlon, mat$statlat, 15, 0.9, addplus=T); text(-90.5,26.25, "index of fecundity")
plotonmap(mat$Ftotal, mat$statlon, mat$statlat, 15, 0.9, addplus=T); text(-90.5,26.25, "index of fecundity")
#plotonmaphires(mat$Ftotal, mat$statlon, mat$statlat, 15, 0.9, addplus=T, adjlab=0.4)      

dat <- data.frame(cbind(mat$statlon, mat$statlat, mat$Ftotal))
names(dat) <- c("lon", "lat", "N")

################################################################################                                                    
co <- read.table("C:/Users/mkarnauskas/Desktop/RS_FATEproject/MASTER_codes/CMS_input_files/redSnapperSett_GOM_ATL.xyz", header=F)

################  plot original recruitment habitat grid cells  ################
plot(1, xlim=c(-98,-81), ylim=c(24,31))
for (j in unique(co[,3]))  {
  m <- co[which(co[,3]==j),]; polygon(m, col=j)
  text(mean(m[,1]), mean(m[,2]), m[1,3], cex=0.6)      }
points(dat$lon, dat$lat, pch=19, cex=0.5)             # release locations

##########  get polygon numbers for release locations  ##############
pts = SpatialPoints(cbind(dat$lon, dat$lat))
m <- co[which(co[,3]==1), ]
m <- rbind(m, m[1,])
oldpol <- Polygons(list(Polygon(cbind(m[,1],m[,2]))), ID=1)
for (j in 2:max(co$V3))  {
  m <- co[which(co[,3]==j), ]
  m <- rbind(m, m[1,])
  newpol <- Polygons(list(Polygon(cbind(m[,1],m[,2]))), ID=j)
  oldpol <- c(oldpol, newpol)           }
pol = SpatialPolygons(oldpol)
polynams <- over(pts, pol)
rel <- cbind(polynams, dat)
rel
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams/rel$polynams+1)

############ find center point of each polygon  ######
ctrslat <- tapply(co$V2, co$V3, mean)
ctrslon <- tapply(co$V1, co$V3, mean)
points(ctrslon, ctrslat, pch=20, cex=2)

##############  label release site with closest polygon number  ############
for (i in 1:nrow(rel)) {
  if (is.na(rel$polynams[i])) {
     rel$polynams[i] <- which.min(abs(rel$lon[i]-ctrslon) + abs(rel$lat[i]-ctrslat))  } }
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams/(rel$polynams+1))

##########  plot new grid cells  #########################
plot(1, xlim=c(-98,-81), ylim=c(24,31))
for (j in unique(co[,3]))  {
  m <- co[which(co[,3]==j),]; polygon(m, lwd=1)
  text(m[1,1], m[1,2], m[1,3], cex=0.6)         }
map('usa', add=T)  
cols <- rainbow(102)  
points(rel$lon, rel$lat, pch=19, cex=0.5, col=cols[rel$polynams])    #  check calcs  
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams)

###############################     ADD DEPTH      ##################################################
###  !!!! NOTE !!!! 
###  Important to use nest from simulation to be run to avoid particles being trapped below surface 
#
nc <- nc_open("C:/Users/mkarnauskas/Desktop/RS_FATEproject/nest_1_20080501000000_HYCOM150.nc")    
#nc <- nc_open("C:/Users/mkarnauskas/Desktop/RS_FATEproject/nest_1_20070701000000.nc")          # open netcdf and get temp variable

v1 <- nc$var[[1]]
u <- ncvar_get(nc, v1)
dim(u)
v1 <- nc$var[[2]]
v <- ncvar_get(nc, v1)
dim(v)       
nc_close(nc)
lon <- nc$var[[1]]$dim[[1]]$vals - 360
lat <- nc$var[[1]]$dim[[2]]$vals
dep <- nc$var[[1]]$dim[[3]]$vals
cur <- sqrt(u^2 + v^2)

image(lon, lat, cur[,,1])
image(lon, lat, cur[,,1], xlim=c(-100, -80), ylim=c(24,31))

rel$depest <- NA
for (i in 1:nrow(rel)) {  rel$depest[i] <- dep[max(which(!is.na(u[which.min(abs(lon - rel$lon[i])), which.min(abs(lat - rel$lat[i])),])))] }
head(rel)

points(rel$lon, rel$lat, pch=19, col=is.na(rel$depest))
rel$depest[is.na(rel$depest)] <- 0

table(rel$depest)
hist(rel$depest)

length(which(rel$depest<45))
length(which(rel$depest>=45))

rel$spawndep <- rel$depest - 5
rel$spawndep[which(rel$spawndep > 45)] <- 45

minspawndep <- 13

length(which(rel$depest<=minspawndep))
dim(rel)
tapply(rel$N, rel$depest<=minspawndep, sum)
rel <- rel[which(rel$depest>minspawndep),]
dim(rel)

plot(rel$depest, rel$spawndep)
table(rel$spawndep < rel$depest)
table(rel$spawndep)
table(rel$depest, rel$spawndep)

plotonmap(rel$spawndep, rel$lon, rel$lat, 15, 0.6)
cols <- rainbow(102)  
points(rel$lon, rel$lat, pch=19, cex=1.5, col=cols[rel$spawndep])

################################################################################
###########   only run this section after initial Gulf-wide run   ##############
###########   cut out release areas which do not supply S Atl     ##############

# polygon 62 is first to have successful recruitment to S Atl

#rel <- rel[which(rel$polynams > 55 & rel$lat < 28.3),]        # uncomment here to apply cutoff 

plotonmap(rel$spawndep, rel$lon, rel$lat, 15, 0.6)

plot(1, xlim=c(-98,-81), ylim=c(24,31))
map('usa', add=T)  
cols <- rainbow(102)  
points(rel$lon, rel$lat, pch=19, cex=0.5, col=cols[rel$polynams])
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams)
              
###################  MAKE RELEASE FILE  #####################################
######################### input temporal information ###########################

days <- read.table("C:/Users/mkarnauskas/Desktop/RS_FATEproject/RELEASEFILE_INDEXnew.csv", sep=",", header=T)     # this inputs a file which references lunar phases and peak spawning times with all dates from 2004 - 2008
#days$red_snapper[which(days$mo==4)] <- 1
sp <- which(days$red_snapper==1 & days$yr==2008)          # adjust year as necessary to line up with simulation year
sp <- sp[seq(1, length(sp)-1, 6)]                         # every 6 days or every 3 days
#sp <- sp[seq(1, length(sp)-1, 3)]
lis <- days[sp, 1:3]
dim(lis)
lis <- lis[-1,]                                          # so that even number in each year
lis <- lis[which(lis$mo<9),]

table(lis$yr)
table(lis$yr, lis$mo)
dim(lis)
lis2 <- lis

for (i in 2009:2009) {                                  # adjust for multiple years as necessary to match simulation years
   lis2$yr <- i
   lis <- rbind(lis, lis2) }

dim(lis)
table(lis$yr, lis$mo)
table(lis$yr)
mean(table(lis$yr))

################  spawning seasonality from Porch et al. 2015  #################
lis$doy <- NA
dinmon <- c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
for (i in 1:nrow(lis))  {  lis$doy[i] <- (sum(dinmon[1:lis$mo[i]]) + lis$da[i]) / 365  }
lis$spawnact <- (lis$doy/0.536)^(0.536/0.024) * exp((0.536-lis$doy)/0.024)
plot(lis$doy, lis$spawnact)
plot(lis$spawnact)
plot(lis$spawnact[1:70])

lis <- lis[-c(4)]

###############  input spatial site information (from above)  ###################

rel1 <- rel[-c(5)]
head(rel1)
m <- rel1       # 'm' is list of release sites with columns: polygon, lon, lat, number of releases
head(m)

m$N <- m$N * 2 #  - SCALE AS NECESSARY           
mean(m$N); min(m$N); max(m$N)
which(m$N==0)

prod(nrow(lis), nrow(m))
###################### now, making the release file ############################

mat <- as.data.frame(matrix(data=NA, nrow=nrow(lis)*nrow(m), ncol=9))   # empty matrix to be filled
mat[,1] <- rep(m[,1], nrow(lis))                                        # column 1: release polygon number
mat[,2] <- rep(m[,2], nrow(lis))                                        # column 2: release longitude 
mat[,3] <- rep(m[,3], nrow(lis))                                        # column 3: release latitude 
mat[,4] <- rep(m[,5], nrow(lis))                                        # column 4: release depth
mat[,5] <- rep(m[,4], nrow(lis))                                        # column 5: number of particles per release
mat <- mat[order(mat[,1], mat[,5], mat[,3], mat[,2]), ] # !!!!  CHECK THIS WITH NEW DATA   !!!!                                # resort matrix
# mat
mat[,6] <- rep(lis[,1], nrow(m))                                        # column 6: release year
mat[,7] <- rep(lis[,2], nrow(m))                                        # column 7: release month
mat[,8] <- rep(lis[,3], nrow(m))                                        # column 8: release day
mat[,9] <- 0                                                            # column 9: release hour

mat[,10] <- rep(lis[,4], nrow(m))                                       # column 10: scale by spawning activity!!!!
sum(mat[,5])

mat <- mat[order(mat[,6], mat[,7], mat[,8]), ] # !!!!  REORDER SO DATES ARE TOGETHER   !!!!
head(mat)

mat$V5 <- round(mat$V5 * mat$V10)
sum(mat[,5])
head(mat)

matfin <- mat[-c(10)]
head(matfin)
dim(matfin)
table((matfin$V5 >0))
table((matfin$V5 >0))[1] / nrow(matfin)     # this is % lost to rounding down to zero
matfin <- matfin[which(matfin$V5 >0),]
head(matfin)
dim(matfin)
mean(matfin$V5); min(matfin$V5); max(matfin$V5)
sum(matfin$V5)

getwd()
dim(matfin)
dim(matfin)/8

matfinGOM <- matfin
save(matfinGOM, file="C:/Users/mkarnauskas/Desktop/RS_FATEproject/fullGoM_release_HYCOM150.RData")

# double check
table(matfin$V6, matfin$V7)
table(matfin$V6)
matplot(table(matfin$V7, matfin$V6), type="l")
diff(table(matfin$V6))

tapply(matfin$V5, list(matfin$V6, matfin$V7), sum)
matplot(tapply(matfin$V5, list(matfin$V7, matfin$V6), sum), type="l")

f <- which(matfin$V6==2008 & matfin$V7 == 8 & matfin$V8 == 29); length(f)
plotonmap(matfin$V5[f], matfin$V2[f], matfin$V3[f], cexnum=0.6, pchnum=15)

f <- which(matfin$V6==2008 & matfin$V7 == 6 & matfin$V8 == 24); length(f)
plotonmap(matfin$V5[f], matfin$V2[f], matfin$V3[f], cexnum=0.6, pchnum=15)

####################   END GOM RELEASE LOCATIONS   #############################
################################################################################
>>>>>>> 9f03132205832f02064f087874eb662c85e0458b
