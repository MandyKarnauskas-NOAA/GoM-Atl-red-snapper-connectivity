
rm(list=ls())
library(ncdf4)
library(maps)


model <- "SABGOM"  # Mercator, SABGOM or HYCOM

if (model == "Mercator") {
  setwd("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/Hatteras_Mercator")  
  d <- read.table("C://Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/Hatteras_20132017.txt")
}
if (model == "HYCOM") {
  setwd("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/Hatteras_MattHYCOM")  
  d <- read.table("C://Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/Hatteras_20082009.txt")
}
if (model == "SABGOM") {
  setwd("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/Hatteras_SABGOM")  
  d <- read.table("C://Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/Hatteras_2017MayJun.txt")}

# concatenate connectivity files -------------------                                  
filelist <- list.files(path = ".", pattern="con_file")     #   find files
dat <- read.table(filelist[1])
filelist <- filelist[-1]
for (i in filelist)  { newdat <- read.table(i); dat <- rbind(dat, newdat) }
colnames(dat) <- c("rel_poly","ret_poly","ret_yr","ret_mo","ret_d","age","ret_dep","rel_yr","rel_mo","rel_d")

nrow(dat) / sum(d$V5) * 100     #  Mercator - 2.3% of particles returned
                                #  HYCOM - 14.1% of particles returned
                                #  SABGOM - 9.1%
head(dat)
table(dat$rel_poly)
table(dat$ret_poly)
table(dat$ret_yr)
table(dat$ret_mo)
table(dat$ret_dep)
table(dat$age)

# plot trajectories ----------------------------------------
lon <- c()
lat <- c()
dep <- c()
cod <- c()

numlarv <- 500

for (m in 1:length(grep("traj", dir()))) {
  nam <- paste("traj_file_", m, ".nc", sep="")
  nc <- nc_open(nam, write=FALSE, readunlim=TRUE, verbose=FALSE)
  v2 <- nc$var[[2]]
  lon1 <- ncvar_get(nc, v2)
  lon1[which(lon1>10000)] <- NA
  v3 <- nc$var[[3]]
  lat1 <- ncvar_get(nc, v3)
  lat1[which(lat1>10000)] <- NA
  v4 <- nc$var[[4]]
  dep1 <- ncvar_get(nc, v4)
  dep1[which(dep1>10000)] <- NA
  v6 <- nc$var[[6]]
  cod1 <- ncvar_get(nc, v6)  # 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data
  lon <- cbind(lon, lon1)
  lat <- cbind(lat, lat1)
  dep <- cbind(dep, dep1)
  cod <- c(cod, cod1)
  nc_close(nc)
}

# check!  should have only -4, -1, and 0
table(cod)                                   # 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data
length(which(cod == (-2))) / length(cod) * 100   # percent of particles messed up
hist(dep[2,which(cod==(-1))])
hist(dep[2,])

bad <- data.frame(cbind(lon[1, which(cod == (-2))]-360, lat[1, which(cod == (-2))]))
bad$V3 <- paste(bad[,1], bad[,2])
matrix(unique(bad$V3))

#  plot depth patterns by code number 
par(mfrow=c(4,1), mar = c(4,4,2,1), mgp=c(2.25,1,0))
matplot(-dep[,which(cod == (0))], col="#FF00FF05", pch=19, type="l", ylim = c(-60, 0),
        xlab = "days after spawning", ylab = "depth of particle", main = "can still move")
axis(1, at=1:31, lab=rep("", 31), tck= -0.01)
matplot(-dep[,which(cod == (-1))], col="#FF00FF30", pch=19, type="l", ylim = c(-60, 0), 
        xlab = "days after spawning", ylab = "depth of particle", main = "left area")
axis(1, at=1:31, lab=rep("", 31), tck= -0.01)
matplot(-dep[,which(cod == (-2))], col="#FF00FF30", pch=19, type="l", ylim = c(-60, 0), 
      xlab = "days after spawning", ylab = "depth of particle", main = "close to land")
axis(1, at=1:31, lab=rep("", 31), tck= -0.01)
matplot(-dep[,which(cod == (-4))], col="#FF00FF05", pch=19, type="l", ylim = c(-60, 0), 
        xlab = "days after spawning", ylab = "depth of particle", main = "settled")
axis(1, at=1:31, lab=rep("", 31), tck= -0.01)


# plot all larval fates --------------------
par(mfrow=c(2,2), mex=0.5)
nn <- c("settled", "close to land", "left area", "can still move"); m <- 1   

for (k in c(-4, -2, -1, 0))  { 
map('state', ylim=c(31, 38.5), xlim=c(-82,-72.8), xlab="", ylab="", col=1); box(); axis(1); axis(2, las=2)
# image(x,y,z, col=rainbow(100, start=0.5, end=0.6)[100:1], add=T); box()
mtext(side=3, line=1, nn[m]); m <- m+1
mm <- which(cod == (k))
text(-81, 38.2, paste("n =", length(mm)))
i <- mm[seq(1, length(mm), length.out = 500)]
for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000005")  
                points(lon[1, j] - 360, lat[1, j], col = "#00FF0010", pch = 19, cex = 1)}
for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, 
                       lat[max(which(!is.na(lat[,j]))),j], col="#FF000010", pch=19, cex=0.75) }
} 
 

# close-up of settled larvae -------------------
dev.off()
co <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/CMS_input_files/redSnapperSett_GOM_ATL_hires.xyz", sep="\t")

map('state', ylim=c(32, 36), xlim=c(-80,-74), xlab="", ylab="", col=1); box(); axis(1); axis(2, las=2)
for (i in 1:117)  {
  f <- co[which(co$V3==i),]
  polygon(f, border=1) 
  text(mean(f$V1), mean(f$V2), i, cex = 0.5)}                  

k <- which(cod == (-4))
i <- k[seq(1, length(k), length.out = 200)]
mtext(side=3, line=1, paste0(model, " - successfully settled"), cex=1.2)
  for (j in i) {  lines(lon[,j]-360, lat[,j], col="#00000010")   
                  points(lon[1,i]-360, lat[1,i], col="#00FF0001", pch=19, cex=0.65)   
                  points(lon[max(which(!is.na(lon[,j]))),j]-360, 
                         lat[max(which(!is.na(lat[,j]))),j], col="#FF000030", pch=19, cex=0.75)   }
  legend("topleft", c("start", "end"), pch=19, cex=1, col=c(3,2))

k <- which(cod == (-4))
endpts <- rep(NA, length(k))
i <- 1
  for (j in k) {   endpts[i] <- lat[max(which(!is.na(lat[,j]))),j] 
                    i <- i + 1} 
min(endpts)                             # SABGOM: 33.96841   Mercator: 33.57697  HYCOM: 33.33318
abline(h = min(endpts), col = 2)

hist(endpts)
table(endpts < 33.9) / length(endpts)
  
table(dat$rel_poly)  
table(dat$ret_poly)

min(lat[1,], na.rm = T)
min(lat, na.rm = T)

length(which(cod == -4)) / length(cod) * 100  # should match number above


# final plot -------------------

folder <- "C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/plots/"  # folder for final plots
load("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_CODES/bathy_Hatteras.RData")        # stored GEBCO data and color scheme

cols <- colorRampPalette(brewer.pal(9, 'Blues'))(180)[110:1]
cols[100] <- "white"

#png(filename = paste0(folder, "Hatteras_", model, ".png"), units="in", width=5, height=6, pointsize=12, res=72*20)

map('usa', ylim=c(32.5, 36.9), xlim=c(-77.9,-73.1), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
image(x, y, -log(-z), col= cols, add=T); box()
mtext(side=3, line=1, "ROMS", cex=1.2)

k <- which(cod == (-1) | cod == (0))
i <- k[seq(1, length(k), length.out = 300)] 
for (j in i) {  points(lon[1,i]-360, lat[1,i], col="#00b30020", pch=19, cex=0.85) }
  for (j in i) {    lines(lon[,j]-360, lat[,j], col="#00000015")  }

k <- which(cod == (-4))
i <- k[seq(1, length(k), length.out = 300)]
for (j in i) {  lines(lon[,j]-360, lat[,j], col="#FFFF0015")  } 
#  points(lon[1,i]-360, lat[1,i], col="#00FF0002", pch=19, cex=0.65)   
for (j in i) {  points(lon[max(which(!is.na(lon[,j]))),j]-360, 
         lat[max(which(!is.na(lat[,j]))),j], col="#FF000040", pch=19, cex=0.85)   }

#legend("bottomright", c("release locations", "settlement locations", "all trajectories", "trajectories of settled larvae"), 
 #      pch = c(19, 19, NA, NA), lwd = c(NA, NA, 2, 2), col = c(3, 2, 1, 7))

dev.off()





################# old code #####################


load("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_CODES/bathy_Hatteras.RData")        # stored GEBCO data and color scheme

map('usa', ylim=c(32.5, 38), xlim=c(-78,-72.8), xlab="", ylab="", col=0); box(); axis(1); axis(2, las=2)
image(x, y, -log(-z), col= cols, add=T); box()
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

