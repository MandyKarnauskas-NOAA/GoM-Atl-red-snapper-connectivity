
rm(list=ls())
library(ncdf4)
library(maps)

load("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_CODES/bathy_Hatteras.RData")        # stored GEBCO data and color scheme

cols <- colorRampPalette(brewer.pal(9, 'Blues'))(180)[110:1]
cols[100] <- "white"

png(filename = "C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/plots/Hatteras.png", units="in", width=9/1.1, height=3/1.1, pointsize=12, res=72*20)
par(mfrow = c(1, 3), mex = 0.9, mar = c(3, 6, 2, 0))

modellis <- c("SABGOM", "HYCOM", "Mercator")

for (model in modellis) { 

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

# final plot -------------------

map('usa', ylim=c(32.8, 36.9), xlim=c(-77.9,-73.3), xlab="", ylab="", col=0); box()
degs = seq(77, 74, -1)
a = sapply(degs, function(x) bquote(.(x)*degree ~ W))
axis(1, at = -degs, lab=do.call(expression, a))
degs = seq(33, 37, 1)
a = sapply(degs, function(x) bquote(.(x)*degree ~ N))
axis(2, at = degs, lab=do.call(expression, a), las = 2)
image(x, y, -log(-z), col= cols, add=T); box()
mtext(side=3, line=1, model, cex=1.0)

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
}

dev.off()

