 
rm(list=ls())

# load libraries --------------------------------------

if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4", repos='http://cran.us.r-project.org')
if (!"sp" %in% installed.packages()) install.packages("sp", repos='http://cran.us.r-project.org')
if (!"maps" %in% installed.packages()) install.packages("maps", repos='http://cran.us.r-project.org')
library(maps)
library(ncdf4)
library(sp)

# specify folder for processing ------------------------

model <- "HYCOM"

if (model == "Mercator") {
  setwd("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/Mercator_red")  
  d <- read.table("C://Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/scaledGOMATLrel20132017.txt")
  }
if (model == "HYCOM") {
  setwd("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/GOMHYCOM_small")  
  d <- read.table("C://Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/scaledGOMATLrel20152016.txt")
}
if (model == "SABGOM") {
  setwd("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/SABGOM_large")  
  d <- read.table("C://Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/scaledGOMATLrel20062010_large.txt")
}
if (model == "Atl-HYCOM") {
  setwd("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/MattHYCOM_red")  
  d <- read.table("C://Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/scaledGOMATLrel200809.txt")
}

# concatenate connectivity files -------------------                                  
filelist <- list.files(path = ".", pattern="con_file")     #   find files
dat <- read.table(filelist[1])
filelist <- filelist[-1]
for (i in filelist)  { newdat <- read.table(i); dat <- rbind(dat, newdat) }
colnames(dat) <- c("rel_poly","ret_poly","ret_yr","ret_mo","ret_d","age","ret_dep","rel_yr","rel_mo","rel_d")

# input polygon file -------------------------------

co <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/CMS_input_files/redSnapperSett_GOM_ATL_hires.xyz", sep="\t")

dat$ret_reg <- "GOM"
dat$ret_reg[which(dat$ret_poly >= 77)] <- "ATL"
dat$rel_reg <- "GOM"
dat$rel_reg[which(dat$rel_poly >= 77)] <- "ATL"

# process release file ------------------------------

d$reg <- "GOM"
d$reg[which(d$V1 >= 77)] <- "ATL"
tapply(d$V5, d$reg, sum)  
                           
table(dat$rel_reg, dat$ret_reg)
tab <- table(dat$rel_reg, dat$ret_reg)

aa <- tab[1, 1];  aa    # 0.8394142    = spawned in ATL, recruited in ATL
ga <- tab[2, 1];  ga    # 0.01136878   = spawned in GOM, recruited in ATL
gg <- tab[2, 2];  gg    # 0.4308186    = spawned in GOM, recruited in GOM
ag <- tab[1, 2];  ag    # 0.0003637645 = spawned in ATL, recruited in GOM

rel <- tapply(d$V5, list(d$reg), sum)
rel

gg / rel[2] * 100  # % Gulf-spawned eggs settling in Gulf
ga / rel[2] * 100  # % Gulf-spawned eggs settling in Atl
aa / rel[1] * 100  # % Atl-spawned eggs settling in Atl

aa / (aa + ga) * 100                       # of all recruits in Atlantic, what % came from Atl
ga / (aa + ga) * 100                      # of all recruits in Atlantic, what % came from Gulf
aa/ length(which(dat$ret_poly >= 77))  # of all recruits in Atlantic, what % came from Atl

surv <- table(dat$rel_reg)
surv/rel

#dat <- dat[which(dat$rel_poly>64),]

barplot(table(dat$rel_reg, dat$ret_yr), beside=T)
table(dat$rel_reg, dat$rel_poly)         # polygon 68 is first to have successful recruitment to S Atl
(table(dat$rel_reg) / nrow(dat)) * 100

# check particle behaviors ---------------------

rm(list = ls()[which(ls()!="model")])
getwd()

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

dev.off()
par(mfrow=c(6,5), mex=0.6)
for (i in 1:30) { hist(dep[i,], breaks=seq(0,100,5)) }

table(dep[2,]==0)/ncol(dep)
table(dep[2,]<10)/ncol(dep)
table(dep[2,]>10 & dep[2,]<20)/ncol(dep)
table(dep[2,]>20 & dep[2,]<30)/ncol(dep)

dev.off()
matplot(-dep[,seq(1, ncol(dep), length.out=500)], col="#FF00FF30", pch=19, type="l", 
        xlab = "days after spawning", ylab = "depth of particle")
axis(1, at=1:31, lab=rep("", 31), tck= -0.01)

#matplot(-dep, col="#FF00FF30", pch=19, type="l", 
#        xlab = "days after spawning", ylab = "depth of particle")
#axis(1, at=1:31, lab=rep("", 31), tck= -0.01)

#  plot by code number 
#par(mfrow=c(4,1), mar = c(4,4,2,1), mgp=c(2.25,1,0))
#matplot(-dep[,which(cod == (0))], col="#FF00FF05", pch=19, type="l", ylim = c(-90, 0),
#        xlab = "days after spawning", ylab = "depth of particle", main = "can still move")
#axis(1, at=1:31, lab=rep("", 31), tck= -0.01)

#matplot(-dep[,which(cod == (-1))], col="#FF00FF30", pch=19, type="l", ylim = c(-90, 0), 
#        xlab = "days after spawning", ylab = "depth of particle", main = "left area")
#axis(1, at=1:31, lab=rep("", 31), tck= -0.01)

matplot(-dep[,which(cod == (-2))], col="#FF00FF30", pch=19, type="l", ylim = c(-90, 0), 
        xlab = "days after spawning", ylab = "depth of particle", main = "close to land")
axis(1, at=1:31, lab=rep("", 31), tck= -0.01)

#matplot(-dep[,which(cod == (-4))], col="#FF00FF05", pch=19, type="l", ylim = c(-90, 0), 
#        xlab = "days after spawning", ylab = "depth of particle", main = "settled")
#axis(1, at=1:31, lab=rep("", 31), tck= -0.01)

# end plot -------

#plot(dens, -dep[2,], col="#FF00FF30")
#plot(diam, -dep[2,], col="#FF00FF30")
#min(dens)
#max(dens)

##########################  PLOT TRAJECTORIES  ###############################

load("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/GEBCO_bathy.RData")        # stored GEBCO data and color scheme

numlarv <- 500

#  i <- which(cod == (-4))
#  la <- rep(NA, length(i))
#      for (j in i) {   la[j] <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]  } 
#  maxlat[m] <- mean(la, na.rm=T)     }

# 4-panel plot of particle states
dev.off()
par(mfrow = c(2, 2), mex = 0.8, mgp = c(1,1,0))
map('usa', xlim=c(-90, -76), ylim=c(23, 35), col=1)
axis(1); axis(2); box(); mtext(side = 3, line = 1, "too close to land")
  i <- which(cod == (-2))
  k <- i #[seq(1, length(i), length.out=numlarv)]
  for (j in k) {  lines(lon[,j]-360, lat[,j], col="#FF00FF20")
    lo <- lon[!is.na(lon[,j]),j][length(lon[!is.na(lon[,j]),j])] - 360
    la <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]
    points(lo, la, col="#FF000010", pch=19, cex=1)     } 
  points(lon[1, i]-360, lat[1, i], col = "#00000020", pch = 19)

map('usa', xlim=c(-90, -76), ylim=c(23, 35), col=1)
axis(1); axis(2); box(); mtext(side = 3, line = 1, "left area") 
i <- which(cod == (-1))                                             #  -1 left area
k <- i[seq(1, length(i), length.out=numlarv)]  
for (j in k) {
  points(lon[1,j]-360, lat[1,j], col="#00000050", pch=19, cex=1)
  lines(lon[,j]-360, lat[,j], col="#FF00FF20")  }      

map('usa', xlim=c(-90, -76), ylim=c(23, 35), col=1)
axis(1); axis(2); box(); mtext(side = 3, line = 1, "can still move") 
i <- which(cod == (0))                                              # 0 can still move
k <- i[seq(1, length(i), length.out=numlarv)]  
for (j in k) {
  lo <- lon[!is.na(lon[,j]),j][length(lon[!is.na(lon[,j]),j])] - 360
  la <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]
  points(lon[1,j]-360, lat[1,j], col="#00000050", pch=19, cex=1)
  lines(lon[,j]-360, lat[,j], col="#FF00FF20") 
  if (lon[1,j]-360 < (-81.7))  {  points(lo, la, col="#00FF0090", pch=19, cex=1.2)  }  else { 
    points(lo, la, col="#FF000050", pch=19, cex=1)  } }      

map('usa', xlim=c(-90, -76), ylim=c(23, 35), col=1)
axis(1); axis(2); box(); mtext(side = 3, line = 1, "settled") 
i <- which(cod == (-4))                                           # -4 settled
k <- i[seq(1, length(i), length.out=numlarv)]  
for (j in k) {
  lo <- lon[!is.na(lon[,j]),j][length(lon[!is.na(lon[,j]),j])] - 360
  la <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]
  points(lon[1,j]-360, lat[1,j], col="#00000025", pch=19, cex=1)
  lines(lon[,j]-360, lat[,j], col="#FF00FF20") 
  if (lon[1, j]-360 < (-81.7))  {  points(lo, la, col="#00FF0090", pch=19, cex=1)  }  else { 
                                   points(lo, la, col="#FF000050", pch=19, cex=1)    }  }      


# plot all particles spawned in Gulf that settled in Atl ------------------
dev.off()
map('usa', xlim=c(-90, -76), ylim=c(23,35), col=1)
axis(1); axis(2); box(); mtext(side = 3, line = 1, "all GoM particles to Atl") 
k <- which(cod == (-4))
k <- k[seq(1, length(k), length.out = 500)]
for (j in k) {
  lo <- lon[!is.na(lon[,j]),j][length(lon[!is.na(lon[,j]),j])] - 360
  la <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]
  if (lon[1,j]-360 < (-81.7) & lo > (-81.7))  { 
    points(lon[1,j]-360, lat[1,j], col="#00000050", pch=19, cex=1)
    lines(lon[,j]-360, lat[,j], col="#FF00FF20") 
    points(lo, la, col="#00FF0090", pch=19, cex=1.2)  }  
}      
abline(v = -81.7)

#le <- length(which(cod == (-4)))
#endlo <- rep(NA, le)
#stalo <- rep(NA, le)
#for (i in 125270:le)  {
#  j <- which(cod == (-4))[i]
#endlo[i] <- lon[!is.na(lon[,j]),j][length(lon[!is.na(lon[,j]),j])] - 360
#stalo[i] <- lon[1, j] - 360  }
#gg <- length(which(stalo < (-81.7) & endlo < (-81.7))); gg
#ga <- length(which(stalo < (-81.7) & endlo >= (-81.7))); ga
#aa <- length(which(stalo >= (-81.7) & endlo >= (-81.7))); aa
#ag <- length(which(stalo >= (-81.7) & endlo < (-81.7))); ag
#aa / (aa + ga) 

# all trajectories on one map ---------------------------------

map('usa', xlim=c(-88, -76), ylim=c(23.85,35), col=1)
#image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1); axis(2, las=2)

i <- 1:length(cod)   
k <- i[seq(1, length(i), length.out = 300)]
  for (j in k) {  lines(lon[,j]-360, lat[,j], col="#FFFF0010") }
  i <- which(cod == (-4))
  k <- i[seq(1, length(i), length.out=200)]  
    for (j in k) {
      lo <- lon[!is.na(lon[,j]),j][length(lon[!is.na(lon[,j]),j])] - 360
      la <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]
          lines(lon[,j]-360, lat[,j], col="#00000020")
          points(lon[1,j]-360, lat[1,j], col="#00000050", pch=19, cex=1)
          if (lon[1,j]-360 < (-81.7))  {  points(lo, la, col="#00FF0070", pch=19, cex=1.2)  }  else { 
                                          points(lo, la, col="#FF000050", pch=19, cex=1.2)  }
         }
    
legend("topleft", c("all trajectories", "trajectories of successful recruits", "release locations", "settlement locations of \nlarvae released from GoM", "settlement locations of \nlarvae released from Atl"), 
 lty=c(1,1,0,0,0), lwd=c(2,2,0,0,0), col=c("yellow", 1, 1, 3, 2), pch=c(-1, -1, 19, 19, 19), cex=1.1, y.intersp=1.2, bty="n")  
#legend("topleft", c("all trajectories", "trajectories of successful recruits", "release locations", "settlement locations of \nlarvae released from Atl"), 
# lty=c(1,1,0,0), lwd=c(2,2,0,0), col=c("yellow", 1, 1, 2), pch=c(-1, -1, 19, 19), cex=1.0, y.intersp=1.0, bty="n")  

# jun 12
### all trajectories on one map  -- full GoM   

par(mar=c(7,7,1,1))
map('usa', xlim=c(-98, -76), ylim=c(23.84,35), col=0, xlab="longitude", ylab="latitude")
image(x,y,z, col=cols, axes=T, add=T, xlab="longitude", ylab="latitude"); box(); axis(1); axis(2, las=2)
mtext(side=1, line=2, "longitude"); mtext(side=2, line=2.5, "latitude")
mtext(side=3, line=1, model, font=2)

i <- 1:length(cod)   #  which(cod == (-4))        # 
k <- i[seq(1, length(i), length.out=400)]
  for (j in k) {  lines(lon[,j]-360, lat[,j], col="#FFFF0010") } 
i <- which(cod == (-4))  
k <- i[seq(1, length(i), length.out=200)]
  for (j in k) { #which(cod == (-4))  ) {
      lo <- lon[max(which(!is.na(lon[,j]))),j] - 360
      la <- lat[max(which(!is.na(lat[,j]))),j]
    if (lo > (-81.7))  {  
          lines(lon[,j]-360, lat[,j], col="#00000010")
          points(lon[1,j]-360, lat[1,j], col="#00FF0050", pch=19, cex=1)   
          points(lo, la, col="#FF000050", pch=19, cex=1) 
             } }
legend("topleft", c("all larval trajectories", "trajectories of larvae successfully recruiting to Atlantic", "source locations of larvae recruiting to Atlantic", "settlement locations of larvae released from Gulf"), 
 lty=c(1,1,0,0), lwd=c(2,2,0,0), col=c("yellow", 1, "green", 2), pch=c(-1, -1, 19, 19), cex=1.0, y.intersp=1.2, bty="n")  


######################   CONNECTIVITY MATRIX   ##########################

rm(list = ls()[which(ls()!="model")])
library(maps)
library(matlab)
getwd()

##############  concatenate files  ###################
filelist <- list.files(path = ".", pattern="con_file")     #   find files
dat <- read.table(filelist[1])
filelist <- filelist[-1]
for (i in filelist)  { newdat <- read.table(i); dat <- rbind(dat, newdat) }
colnames(dat) <- c("rel_poly","ret_poly","ret_yr","ret_mo","ret_d","age","ret_dep","rel_yr","rel_mo","rel_d")
#########################################################

d <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/CMS_input_files/redSnapperSett_GOM_ATL_hires.xyz", sep="\t")
names(d) <- c("lon", "lat", "pol")

##########  plot new grid cells  #########################

pols <- unique(d$pol)
tklocs <- c(89.5, 76.5, 95.5, 117.5, 106.5, 41.5) 
lablocs <- c(21.0, 59.0, 83.0, 92.5, 101.0, 112.0) 

barplot(table(dat$ret_poly))
barplot(table(dat$rel_poly))
barplot(table(dat$ret_yr))
barplot(table(dat$ret_mo))
barplot(table(dat$ret_dep))
barplot(table(dat$age))

# plot for all years together ----------------

con <- dat
colnames(con) <- c("source","sink","yr","mon","day","tim","ret_dep","rel_yr","rel_mo","rel_d")
pollis <- 1:max(d$pol)    # list of polygon names
limit <- length(pollis)

nf <- layout(matrix(c(1:2), 2, 1), c(10,10), c(10, 2))   # plot layout -- all years
layout.show(nf)

par(mar=c(5,5,3,2), xpd=F)
con1 <- con                     # to look at all years together
settle <- matrix(NA, length(pollis), length(pollis))
for (i in 1:length(pollis))     {
    a <- which(con1$source==pollis[i])
      if(length(a)>0) {
        for (j in 1:length(pollis)) {   settle[i,j] <- length(which(con1$sink[a]==j))    }
    } else { settle[i,] <- 0  }  }
settle <- t(settle)
nn <- 300
cl = jet.colors(nn)
image(pollis[42:limit], pollis[42:limit], log(settle[42:limit,42:limit]+1), col=cl, axes=F,  xlab="", ylab="", 
main= paste("GoM - Atl connectivity --", model))

axis(1, at=tklocs, lab=rep("", 6))
axis(2, at=tklocs, lab=rep("", 6))
axis(2, at=lablocs, lab=c("W GOM", "E GOM", "E FL", "GA", "SC", "NC"), tick=F, cex=1, las=2)
axis(1, at=lablocs, lab=c("W GOM", "E GOM", "E FL", "GA", "SC", "NC"), tick=F, cex=1, las=1)
box(); abline(0,1, lty=2)
abline(v=tklocs[2], lty=1); abline(h=tklocs[2], lty=1)

mtext(side=2, line=3.5, "Source Node", font=2)
mtext(side=1, line=3.5, "Receiving Node", font=2)

####  color bar plot
par(mar=c(0,2,0,0))   
plot(c(-0.5,1.5), c(1.05,1.18), axes=F, xlab="", ylab="", col="white")
bra <- (seq(-0.5, ceiling(log(max(settle))), (ceiling(log(max(settle)))+0.5)/nn))
i = seq(0,1,1/(nn-1))
rect(i, 1.1, 1.04, 1.12, col = cl, lwd=0, border=cl)
laa <- c(1,10,100, 1000)   # laa <- seq(0, 60000, 10000)
levs <- (0+i[1:nn]+1/(nn*2))
re <- lm(bra[1:(length(bra)-1)]~levs)
text((log(laa) - re$coef[1])/re$coef[2], 1.09, laa, pos=4, cex=1.2)
text(0.5,1.12, cex=1.2, "number of successful recruits", pos=3)



#########  plot layout -- individual years

nf <- layout(matrix(c(1:4, 1, 5:7, rep(8,4)), 3, 4, byrow=TRUE), c(1,rep(7,3)), c(7,7,1.5))
layout.show(nf)

#nf <- layout(matrix(c(1:4, 1, 5:7, 1, 8, 9, 9, rep(10,4)), 4, 4, byrow=TRUE), c(1,rep(7,3)), c(7,7,7,1.5))
#layout.show(nf)
#nf <- layout(matrix(c(1:5, 14, 1, 6:9, 14, 1, 10:13, 14, 1, rep(15,5)), 4, 6, byrow=TRUE), c(1,rep(7,4), 3), c(7,7,7,1.5))
#layout.show(nf)

####  global y label
par(mar=c(0,0,0,0))
plot.new()
mtext(side=2, line=-5, "Source Node", font=2)

#### plot connectivity matrices by year
par(mar=c(3,3,3,2)-2, xpd=F)
resmat <- matrix(NA, nrow=13, ncol=5)
limit <- length(pollis)

for (k in unique(con$yr))  {
  con1 <- con[which(con$yr == k),]

par(mar=c(5,5,3,2), xpd=F)
settle <- matrix(NA, length(pollis), length(pollis))
for (i in 1:length(pollis))     {
  a <- which(con1$source==pollis[i])
  if(length(a)>0) {
    for (j in 1:length(pollis)) {   settle[i,j] <- length(which(con1$sink[a]==j))    }
  } else { settle[i,] <- 0  }  }
settle <- t(settle)
nn <- 300
cl = jet.colors(nn)
image(pollis[42:limit], pollis[42:limit], log(settle[42:limit,42:limit]+1), col=cl, axes=F,  xlab="", ylab="", 
      main= paste(model, "model -", k))

axis(1, at=tklocs, lab=rep("", 6))
axis(2, at=tklocs, lab=rep("", 6))
axis(2, at=lablocs, lab=c("W GOM", "E GOM", "E FL", "GA", "SC", "NC"), tick=F, cex=1, las=2)
axis(1, at=lablocs, lab=c("W GOM", "E GOM", "E FL", "GA", "SC", "NC"), tick=F, cex=1, las=2)
box(); abline(0,1, lty=2)
abline(v=tklocs[2], lty=1); abline(h=tklocs[2], lty=1)
}

####  color bar plot
par(mar=c(0,2,0,0))   
plot(c(1.05,1.18), c(-0.5,1.5),  axes=F, xlab="", ylab="", col="white")
bra <- (seq(-0.5, ceiling(log(max(settle))), (ceiling(log(max(settle)))+0.5)/nn))
i = seq(0,1,1/(nn-1))
#rect(1.05, i, 1.1, 1, col = cl, lwd=0, border=cl)
rect(1.11, i, 1.12, 1.3, col = cl, lwd=0, border=cl)
laa <- c(1,10,100,1000)   # laa <- seq(0, 60000, 10000)
levs <- (0+i[1:nn]+1/(nn*2))
re <- lm(bra[1:(length(bra)-1)]~levs)
poss <- (log(laa) - re$coef[1])/re$coef[2]
text(1.12, poss, laa, pos=4, cex=1.2)
text(1.115, -0.15, cex=1.2, "number of successful recruits", pos=3)
#legend("bottom", "self-recruitment", lwd=1, lty=2, bty="n", cex=1.1)

###  global x label
par(mar=c(0,0,0,0))
plot.new()
mtext(side=1, line=-5, "Receiving Node", font=2)
#





# old code 

##########  plot new grid cells  #########################

pols <- unique(d$pol)
are <- rep(NA, length(unique(d$pol)))

plot(0, ylim=c(24, 36), xlim=c(-100, -75))  # check
for (i in 1:117)  {
  f <- d[which(d$pol==i),]
  polygon(f, border=1)  }
map('state', add=T)

text(tapply(d$lon, d$pol, mean), tapply(d$lat, d$pol, mean), unique(d$pol), cex=0.5)

are[1:41] <- "W GoM"
are[42:76] <- "E GoM"
are[77:89] <- "E FL"
are[90:95] <- "GA"
are[96:106] <- "SC"
are[107:117] <- "NC"

are <- as.factor(are)
tab <- table(pols, are)

map('state', fill = 1, interior=F, col = gray(0.95), ylim=c(23,36), xlim=c(-99,-75))
for (j in unique(d[,3]))  {
  m <- d[which(d[,3]==j),]; polygon(m, lwd=3, col=as.numeric(are)[j]+1, border = 1, lwd=0.5)  }
axis(1); axis(2); box()
for (j in unique(are)) {  text( mean(d$lon[d$pol %in% which(are==j)]),  mean(d$lat[d$pol %in% which(are==j)]), j, col=1) }
abline(v=-89.1, lty=2); text(-87.5, 24, "longitude = -89.1")

mids <- tapply(pols, are, mean)
maxs <- tapply(pols, are, max)
mins <- tapply(pols, are, min)

tklocs <- maxs + 0.5
lablocs <- sort((tklocs-mins-0.5)/2 + mins)       # used in plotting later

tklocs <- c(89.5, 76.5, 95.5, 117.5, 106.5, 41.5) 
lablocs <- c(21.0, 59.0, 83.0, 92.5, 101.0, 112.0) 

efl <- which(tab[,1]==1)
egom <- which(tab[,2]==1)
ga <- which(tab[,3]==1)
nc <- which(tab[,4]==1)
sc <- which(tab[,5]==1)
wgom <- which(tab[,6]==1)

dat$state <- NA
dat$state[which(dat$ret_poly<=max(egom))] <- "E GOM"
dat$state[which(dat$ret_poly>=min(wgom) & dat$ret_poly<=max(wgom))] <- "W GOM"
dat$state[which(dat$ret_poly>=min(efl) & dat$ret_poly<=max(efl))] <- "E FL"
dat$state[which(dat$ret_poly>=min(ga) & dat$ret_poly<=max(ga))] <- "GA"
dat$state[which(dat$ret_poly>=min(sc) & dat$ret_poly<=max(sc))] <- "SC"
dat$state[which(dat$ret_poly>=min(nc))] <- "NC"


##################################################################################
