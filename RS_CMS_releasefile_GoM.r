################################################################################
###########################   red snapper release file for CMS  ################
###########################   M. Karnauskas Jan 8, 2015         ################
#
#   Creates release file for input to larval transport simulation
#   code takes output from GLM model and krigging 
#   includes final fecundity map from published version
#
#   Mar2020 Updated work flow: 
#    - build file for separate regions using this code and RS_CMS_releasefile_Atl.r
#    - extract date range of 95% spawning activity
#    - combine using scale_releases_1to1.r
################################################################################

rm(list=ls())

if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4", repos='http://cran.us.r-project.org')
if (!"sp" %in% installed.packages()) install.packages("sp", repos='http://cran.us.r-project.org')
if (!"maps" %in% installed.packages()) install.packages("maps", repos='http://cran.us.r-project.org')
library(ncdf4)
library(sp)
library(maps)

# load data --------------------------------------------------
source("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/plotting.r")                                                                                                  
load("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINAL_MAPPING_RESULTS_ext.RData")        #  Fwith_platform is index to be used
head(mat)

#plotonmap(mat$Fwith_resid, mat$statlon, mat$statlat, 15, 0.9, addplus=T); text(-90.5,26.25, "index of fecundity")
plotonmap(mat$Ftotal, mat$statlon, mat$statlat, 15, 0.9, addplus=T); text(-90.5,26.25, "index of fecundity")
#plotonmaphires(mat$Ftotal, mat$statlon, mat$statlat, 15, 0.9, addplus=T, adjlab=0.4)      

dat <- data.frame(cbind(mat$statlon, mat$statlat, mat$Ftotal))
names(dat) <- c("lon", "lat", "N")

# assign polygon numbers -------------------------------------
co <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/CMS_input_files/redSnapperSett_GOM_ATL_hires.xyz", header=F)

# plot original recruitment habitat grid cells
plot(1, xlim=c(-98,-81), ylim=c(24,31))
for (j in unique(co[,3]))  {
  m <- co[which(co[,3]==j),]; polygon(m, col=j)
  text(mean(m[,1]), mean(m[,2]), m[1,3], cex=0.6)      }
points(dat$lon, dat$lat, pch=19, cex=0.5)             # release locations

# get polygon numbers for release locations 
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

# find center point of each polygon 
ctrslat <- tapply(co$V2, co$V3, mean)
ctrslon <- tapply(co$V1, co$V3, mean)
points(ctrslon, ctrslat, pch=20, cex=2)

# label release site with closest polygon number
for (i in 1:nrow(rel)) {
  if (is.na(rel$polynams[i])) {
     rel$polynams[i] <- which.min(abs(rel$lon[i]-ctrslon) + abs(rel$lat[i]-ctrslat))  } }
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams/(rel$polynams+1))

# plot new grid cells  
plot(1, xlim=c(-98,-81), ylim=c(24,31))
for (j in unique(co[,3]))  {
  m <- co[which(co[,3]==j),]; polygon(m, lwd=1)
  text(m[1,1], m[1,2], m[1,3], cex=0.6)         }
map('usa', add=T)  
cols <- rainbow(102)  
points(rel$lon, rel$lat, pch=19, cex=0.5, col=cols[rel$polynams])    #  check calcs  
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams)

# Assign depths ------------------------------------------------

###################  add depths from oceanographic models  ######################

###  !!!! NOTE !!!! 
###  Important to use nest from simulation to be run to avoid particles being trapped below surface 
#
nc1 <- nc_open("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/nest_1_20080501000000_HYCOM150.nc")  # Matthieu model  
nc2 <- nc_open("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/nest_1_20070701000000.nc")          # SABGOM model
#nc <- nc_open("C:/Users/mkarnauskas/Desktop/RS_FATEproject/nest_2_20070701000000.nc")          # open netcdf and get temp variable

var1 <- nc1$var[[1]]
u1 <- ncvar_get(nc1, var1)
dim(u1)
var2 <- nc1$var[[2]]
v1 <- ncvar_get(nc1, var2)
dim(v1)       
nc_close(nc1)

var1 <- nc2$var[[1]]
u2 <- ncvar_get(nc2, var1)
dim(u2)
var2 <- nc2$var[[2]]
v2 <- ncvar_get(nc2, var2)
dim(v2)       
nc_close(nc2)

lon1 <- nc1$var[[1]]$dim[[1]]$vals - 360
lat1 <- nc1$var[[1]]$dim[[2]]$vals
dep1 <- nc1$var[[1]]$dim[[3]]$vals
cur1 <- sqrt(u1^2 + v1^2)

lon2 <- nc2$var[[1]]$dim[[1]]$vals - 360
lat2 <- nc2$var[[1]]$dim[[2]]$vals
dep2 <- nc2$var[[1]]$dim[[3]]$vals
cur2 <- sqrt(u2^2 + v2^2)

image(lon1, lat1, cur1[,,1])
image(lon2, lat2, cur2[,,1]) #, xlim=c(-85, -75), ylim=c(24,31))

rel$depest <- NA
for (i in 1:nrow(rel)) {  
  depth1 <- dep1[max(which(!is.na(u1[which.min(abs(lon1 - rel$lon[i])), which.min(abs(lat1 - rel$lat[i])),])))]
  depth2 <- dep2[max(which(!is.na(u2[which.min(abs(lon2 - rel$lon[i])), which.min(abs(lat2 - rel$lat[i])),])))]  
  rel$depest[i] <- min(depth1, depth2)
}

# modify depths
points(rel$lon, rel$lat, pch=19, col=is.na(rel$depest))
rel$depest[is.na(rel$depest)] <- 0

table(rel$depest)
hist(rel$depest)

length(which(rel$depest<45))
length(which(rel$depest>=45))

rel$spawndep <- rel$depest - 5       # set spawning 5m above sea floor
rel$spawndep[which(rel$spawndep > 45)] <- 45

# ASSIGN MIN SPAWNING DEPTH HERE
minspawndep <- 13                    #  set minimum spawning depth here 

length(which(rel$depest<=minspawndep))
dim(rel)
tapply(rel$N, rel$depest<=minspawndep, sum)   # number of particles lost by setting min spawning depth
rel <- rel[which(rel$depest>minspawndep),]
dim(rel)

plot(rel$depest, rel$spawndep)
table(rel$spawndep < rel$depest)
table(rel$spawndep)
table(rel$depest, rel$spawndep)

plotonmap(rel$spawndep, rel$lon, rel$lat, 15, 0.6)
cols <- rainbow(102)  
points(rel$lon, rel$lat, pch=19, cex=1.0, col=cols[rel$spawndep])


#  STOP AND ASSESS!!!! -------------------------------

#####   only run this section after initial Gulf-wide run   #####
#####   cut out release areas which do not supply S Atl     #####

# polygon 62 is first to have successful recruitment to S Atl
# polygon 46 is first FL polygon -- cutoff for scaling abundance E/W coast

rel <- rel[which(rel$polynams >= 46),]        # uncomment here to apply cutoff     

plotonmap(rel$spawndep, rel$lon, rel$lat, 15, 0.6)

plot(1, xlim=c(-98,-81), ylim=c(24,31))
map('usa', add=T)  
cols <- rainbow(102)  
points(rel$lon, rel$lat, pch=19, cex=0.5, col=cols[rel$polynams])
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams)
abline(v=(-89))
              
# TEMPORAL INFORMATION ----------------------------------------

days <- as.Date(format("2013-01-01"))+0:364        #  modify starting year here 
days

sp <- data.frame(substr(days, 1, 4))
names(sp) <- "yr"
sp$yr <- as.numeric(as.character(sp$yr))
sp$mo <- as.numeric(as.character(substr(days, 6, 7)))
sp$da <- as.numeric(as.character(substr(days, 9, 10)))
sp$doy <- (1:365)/365   
sp              

lis <- sp[seq(4, length(days), 6), ]               # releases every 6 days 
head(lis)                                          # make sure this lines up with Atl file - first release date of May 16
lis[which(lis$mo == 5 & lis$da == 22),]

# add spawning activity by doy ------------------------------

lis$spawnact <- (lis$doy/0.536)^(0.536/0.024) * exp((0.536-lis$doy)/0.024)   # spawning relationship from Porch et al. 2015
plot(lis$doy, lis$spawnact, xlab="day of year", ylab="spawning activity")

cutoff <- 0.25; abline(h=cutoff, col=2)                   # determine cutoff here - aim for 90%  
sum(lis$spawnact[which(lis$spawnact > cutoff)]) / sum(lis$spawnact)   # calculate percentile spawning activity included with given cutoff                    

lis <- lis[which(lis$spawnact > cutoff), ]
points(lis$doy, lis$spawnact, col=2, cex=0.8, pch=20)      # check cutoff
lis

table(lis$yr)
table(lis$yr, lis$mo)
dim(lis)    

lis2 <- lis
#for (i in 2005:2010) {                             #  for multiple years
#   lis2$yr <- i
#   lis <- rbind(lis, lis2) }

dim(lis)
table(lis$yr, lis$mo)
table(lis$yr)
mean(table(lis$yr))

plot(lis$doy, lis$spawnact)
plot(lis$spawnact)

lis <- lis[-c(4)]
head(lis) 

# input spatial site information (from above) ---------------------

rel1 <- rel[-c(5)]
head(rel1)
m <- rel1       # 'm' is list of release sites with columns: polygon, lon, lat, number of releases
head(m)

m$N <- m$N * 10  #  - SCALE AS NECESSARY  - scaling up dramatically so that we can reduce to achieve desired ratio with Atl       
mean(m$N); min(m$N); max(m$N)
which(m$N==0)

prod(nrow(lis), nrow(m))

# now, making the release file ------------------------------------

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
sum(matfin$V5)                              # currently at 3 million particles -- will reduce when combined

# double check release file ----------------------------
table(matfin$V6, matfin$V7)
table(matfin$V6)
matplot(table(matfin$V7, matfin$V6), type="l")
diff(table(matfin$V6))

tapply(matfin$V5, list(matfin$V6, matfin$V7), sum)
matplot(tapply(matfin$V5, list(matfin$V7, matfin$V6), sum), type="l")

f <- which(matfin$V6==2013 & matfin$V7 == 5 & matfin$V8 == 16); length(f)   # non-peak spawning 
plotonmap(matfin$V5[f], matfin$V2[f], matfin$V3[f], cexnum=0.6, pchnum=15)

f <- which(matfin$V6==2013 & matfin$V7 == 6 & matfin$V8 == 27); length(f)   # peak spawning
plotonmap(matfin$V5[f], matfin$V2[f], matfin$V3[f], cexnum=0.6, pchnum=15)

# save output -----------------------------------------
getwd()
matfinGOM <- matfin

save(matfinGOM, file="C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/GOMreleaseForScaling.RData")

# END -------------------------------------------------

