
rm(list=ls())
d <- read.table("C://Users/mkarnauskas/Desktop/rfvm_sm_sta.csv", sep=",", header=T)

source("C:/Users/mkarnauskas/Desktop/completed_manuscripts/RS_mapping_paper/plotting.r")

d$LUTJANUS_CAMPECHANUS[is.na(d$LUTJANUS_CAMPECHANUS)]<- 0
plotonmap(d$LUTJANUS_CAMPECHANUS, d$sta_lon, d$sta_lat, 15, 0.4)

d1 <- d[which(d$LUTJANUS_CAMPECHANUS>0),]
plot(d1$sta_lon, d1$sta_lat)

map('usa', ylim=c(22, 31), xlim=c(-86,-80), xlab="", ylab="", col=1); box(); axis(1); axis(2, las=2)
points(d1$sta_lon, d1$sta_lat)
points(d1$end_lon, d1$end_lat, col=2)
abline(h=24.95821)

d2 <- d1[which(d1$sta_lat < 24.958),]

dat <- data.frame(cbind(d2$sta_lon, d2$sta_lat, 10))
names(dat) <- c("lon", "lat", "N")
dat

#################################
co <- read.table("C://Users/mkarnauskas/Desktop/red_snapper/Mar2016/redSnapperSett_GOM_ATL.txt", sep="\t", header=F)

################  plot original recruitment habitat grid cells  ################
plot(1, xlim=c(-98,-81), ylim=c(24,31))
for (j in unique(co[,3]))  {
  m <- co[which(co[,3]==j),]; polygon(m, col=j)
  text(m[1,1], m[1,2], m[1,3], cex=0.6)      }
points(dat$lon, dat$lat, pch=19, cex=0.5)             # release locations

##########  get polygon numbers for release locations  ##############

polynams <- 74
rel <- cbind(polynams, dat)
rel

###############################     ADD DEPTH      ##################################################

library(ncdf4)
nc <- nc_open("C:/Users/mkarnauskas/Desktop/RS_FATEproject/Hatteras_issue/Matthieu_HYCOM20082009/nest_1_20080501000000.nc")          # open netcdf and get temp variable
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

#d <- matrix(NA, 976, 694)
#for (i in 1:341) {
#  for (j in 1:257) {
#     d[i,j] <- dep[max(which(!is.na(u[i,j,])))]  }}
#image(lon, lat, d)

rel$depest <- NA
for (i in 1:nrow(rel)) {  rel$depest[i] <- dep[max(which(!is.na(u[which.min(abs(lon - rel$lon[i])), which.min(abs(lat - rel$lat[i])),])))]   }
head(rel)

length(which(rel$depest<45))
length(which(rel$depest>=45))

rel$spawndep <- rel$depest - 2
rel$spawndep[which(rel$spawndep > 45)] <- 45

minspawndep <- 13

length(which(rel$depest<=minspawndep))
dim(rel)
tapply(rel$N, rel$depest<=minspawndep, sum)
rel <- rel[which(rel$depest>minspawndep),]
rel

###################  MAKE RELEASE FILE  #####################################
######################### input temporal information ###########################

days <- read.table("C:/Users/mkarnauskas/Desktop/red_snapper/Jan2015/RELEASEFILE_INDEXnew.csv", sep=",", header=T)     # this inputs a file which references lunar phases and peak spawning times with all dates from 2004 - 2008
#days$red_snapper[which(days$mo==4)] <- 1
sp <- which(days$red_snapper==1 & days$yr==2008)
sp <- sp[seq(1, length(sp)-1, 6)]
#sp <- sp[seq(1, length(sp)-1, 3)]
lis <- days[sp, 1:3]
dim(lis)
lis <- lis[-1,]
lis <- lis[which(lis$mo<9),]

table(lis$yr)
table(lis$yr, lis$mo)
dim(lis)
lis2 <- lis

for (i in 2009) {
   lis2$yr <- i
   lis <- rbind(lis, lis2) }

dim(lis)
table(lis$yr, lis$mo)
table(lis$yr)
mean(table(lis$yr))

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

m$N <- m$N  # * 5
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
table((matfin$V5 >0))[1] / nrow(matfin)
matfin <- matfin[which(matfin$V5 >0),]
head(matfin)
dim(matfin)
mean(matfin$V5); min(matfin$V5); max(matfin$V5)
sum(matfin$V5)

getwd()
dim(matfin)
dim(matfin)/8

write.table(matfin, file="C:/Users/mkarnauskas/Desktop/RS_FATEproject/TortugasRel_RS.txt", sep="\t", col.names=F, row.names=F)

# double check
table(matfin$V6, matfin$V7)
table(matfin$V6)
matplot(table(matfin$V7, matfin$V6), type="l")
diff(table(matfin$V6))

tapply(matfin$V5, list(matfin$V6, matfin$V7), sum)
matplot(tapply(matfin$V5, list(matfin$V7, matfin$V6), sum), type="l")


