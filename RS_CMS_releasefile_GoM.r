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
co <- read.table("C:/Users/mkarnauskas/Desktop/RS_FATEproject/MASTER_codes/CMS_input_files/redSnapperSett_GOM_ATL_hires.xyz", header=F)

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
#nc <- nc_open("C:/Users/mkarnauskas/Desktop/RS_FATEproject/nest_1_20080501000000_HYCOM150.nc")    
nc <- nc_open("C:/Users/mkarnauskas/Desktop/RS_FATEproject/nest_1_20070701000000.nc")          # this is SABGOM nest_1

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

rel$spawndep <- rel$depest - 5       # set spawning 5m above sea floor
rel$spawndep[which(rel$spawndep > 45)] <- 45

minspawndep <- 13                    #  set minimum spawning depth here -- need to check 

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
points(rel$lon, rel$lat, pch=19, cex=1.5, col=cols[rel$spawndep])

################################################################################
###########   only run this section after initial Gulf-wide run   ##############
###########   cut out release areas which do not supply S Atl     ##############

# polygon 62 is first to have successful recruitment to S Atl
# polygon 46 is first FL polygon -- cutoff for scaling abundance E/W coast

rel <- rel[which(rel$polynams >= 46),]        # uncomment here to apply cutoff     

plotonmap(rel$spawndep, rel$lon, rel$lat, 15, 0.6)

plot(1, xlim=c(-98,-81), ylim=c(24,31))
map('usa', add=T)  
cols <- rainbow(102)  
points(rel$lon, rel$lat, pch=19, cex=0.5, col=cols[rel$polynams])
points(rel$lon, rel$lat, pch=19, cex=0.5, col=rel$polynams)
              
########################     TEMPORAL INFORMATION      #########################

days <- as.Date(format("2004-01-01"))+0:364        #  modify starting year here 
days

sp <- data.frame(substr(days, 1, 4))
names(sp) <- "yr"
sp$yr <- as.numeric(as.character(sp$yr))
sp$mo <- as.numeric(as.character(substr(days, 6, 7)))
sp$da <- as.numeric(as.character(substr(days, 9, 10)))
sp$doy <- (1:365)/365   
sp              

lis <- sp[seq(1, length(days), 6), ]               # releases every 6 days 
head(lis)

#  spawning relationship from Porch et al. 2015
lis$spawnact <- (lis$doy/0.536)^(0.536/0.024) * exp((0.536-lis$doy)/0.024)
plot(lis$doy, lis$spawnact)
cutoff <- 0.15; abline(h=cutoff)                   # determine cutoff here 
sum(lis$spawnact[which(lis$spawnact > cutoff)]) / sum(lis$spawnact)   # calculate percentile spawning activity included with given cutoff                    
lis <- lis[which(lis$spawnact > cutoff), ]
points(lis$doy, lis$spawnact, col=2, cex=0.8)      # check cutoff

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

###############  input spatial site information (from above)  ###################

rel1 <- rel[-c(5)]
head(rel1)
m <- rel1       # 'm' is list of release sites with columns: polygon, lon, lat, number of releases
head(m)

m$N <- m$N * 100 #  - SCALE AS NECESSARY  - scaling up dramatically so that we can reduce to achieve desired ratio with Atl       
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
sum(matfin$V5)                              # currently at 3 million particles -- will reduce when combined

# double check
table(matfin$V6, matfin$V7)
table(matfin$V6)
matplot(table(matfin$V7, matfin$V6), type="l")
diff(table(matfin$V6))

tapply(matfin$V5, list(matfin$V6, matfin$V7), sum)
matplot(tapply(matfin$V5, list(matfin$V7, matfin$V6), sum), type="l")

f <- which(matfin$V6==2004 & matfin$V7 == 9 & matfin$V8 == 27); length(f)   # non-peak spawning 
plotonmap(matfin$V5[f]/10, matfin$V2[f], matfin$V3[f], cexnum=0.6, pchnum=15)

f <- which(matfin$V6==2004 & matfin$V7 == 6 & matfin$V8 == 23); length(f)   # peak spawning
plotonmap(matfin$V5[f]/10, matfin$V2[f], matfin$V3[f], cexnum=0.6, pchnum=15)

#############################  save output  ####################################
getwd()
matfinGOM <- matfin
save(matfinGOM, file="C:/Users/mkarnauskas/Desktop/RS_FATEproject/MASTER_codes/GOMreleaseForScaling_SABGOM.RData")

####################   END GOM RELEASE LOCATIONS   #############################
################################################################################
