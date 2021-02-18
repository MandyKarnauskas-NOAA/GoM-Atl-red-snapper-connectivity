 
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
spp <- "mutton"

fold <- paste0("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/final_for_subsample/", model, "_", spp)
setwd(fold) 

if (model == "Mercator")  { d <- read.table("C://Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/scaledGOMATLrel20132017.txt") }
if (model == "HYCOM")     { nd <- read.table("C://Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/scaledGOMATLrel2019.txt")     
                                 d <- rbind(nd, nd, nd, nd, nd);   d$V6 <- sort(rep(2015:2019, nrow(nd))) }    
if (model == "SABGOM")    { d <- read.table("C://Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/scaledGOMATLrel20062010.txt") }
if (model == "Atl-HYCOM") { d <- read.table("C://Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/scaledGOMATLrel200809.txt")   }

# concatenate connectivity files -------------------                                  
filelist <- list.files(path = ".", pattern="con_file")     #   find files
dat <- read.table(filelist[1])
filelist <- filelist[-1]
for (i in filelist)  { newdat <- read.table(i); dat <- rbind(dat, newdat) }
colnames(dat) <- c("rel_poly","ret_poly","ret_yr","ret_mo","ret_d","age","ret_dep","rel_yr","rel_mo","rel_d")

dat$ret_reg <- "GOM"
dat$ret_reg[which(dat$ret_poly >= 77)] <- "ATL"
dat$rel_reg <- "GOM"
dat$rel_reg[which(dat$rel_poly >= 77)] <- "ATL"

# process con file ------------------------------

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

ga/ length(which(dat$ret_poly >= 77))*100  # RECORD THIS NUMBER
gomatl <- c(model, spp, ga/ length(which(dat$ret_poly >= 77))*100)

write.table(matrix(gomatl, 1, 3), file = "C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/outputs.txt", append = T, col.names = F)

surv <- table(dat$rel_reg)
surv/rel

#dat <- dat[which(dat$rel_poly>64),]

barplot(table(dat$rel_reg, dat$ret_yr), beside=T)
table(dat$rel_reg, dat$rel_poly)         # polygon 68 is first to have successful recruitment to S Atl
(table(dat$rel_reg) / nrow(dat)) * 100

# check particle behaviors ---------------------

rm(list = ls()[-which(ls() %in% c("model", "spp", "fold"))])
getwd()
folder <- "C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/plots/"  # folder for final plots

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
table(cod)           # 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data
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
#table(dep[2,]<10)/ncol(dep)
#table(dep[2,]>10 & dep[2,]<20)/ncol(dep)
#table(dep[2,]>20 & dep[2,]<30)/ncol(dep)

dev.off()
matplot(-dep[,seq(1, ncol(dep), length.out=100)], col="#FF00FF30", pch=19, type="l", 
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

load("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/GEBCO_bathy_hires.RData")        # stored GEBCO data and color scheme

numlarv <- 500

#  i <- which(cod == (-4))
#  la <- rep(NA, length(i))
#      for (j in i) {   la[j] <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]  } 
#  maxlat[m] <- mean(la, na.rm=T)     }

# 4-panel plot of particle states
dev.off()

png(filename = paste0(folder, "statusCheck-", model, "-", spp, ".png"), units="in", width=8, height=8, pointsize=12, res=72*20)

par(mfrow = c(2, 2), mex = 0.5, mgp = c(1,1,0))

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

png(filename = paste0(folder, "GoMtoAtl-", model, "-", spp,".png"), units="in", width=4, height=6, pointsize=12, res=72*20)
par(mar = c(3, 3, 2, 2))

map('usa', xlim=c(-87.75, -75.25), ylim=c(23, 35.75), col=1)
image(x[seq(1, 3360, 5)], y[seq(1, 3360, 5)], z[seq(1, 3360, 5), seq(1, 3360, 5)], col=cols, axes=T, xlab="", ylab="", add=T)
axis(1); axis(2, las = 2); box()
mtext(side = 3, line = 1, paste(model, "-", spp), font = 2, cex = 1.2) 

k <- which(cod == (-4))
k <- k[seq(1, length(k), 20)]  # change to 2 or 3
for (j in k) {
  lo <- lon[!is.na(lon[,j]),j][length(lon[!is.na(lon[,j]),j])] - 360
  la <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]
  if (lon[1,j]-360 < (-81.7) & lo > (-81.7))  { 
    points(lon[1,j]-360, lat[1,j], col="#00FF0070", pch=19, cex=0.5)
    lines(lon[,j]-360, lat[,j], col="#FFFF0030") 
    points(lo, la, col="#FF000090", pch=19, cex=0.5)  }  
}      

dev.off()

#abline(v = -81.7)

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

png(filename = paste0(folder, "allTrajectories-", model, "-", spp, ".png"), units="in", width=6, height=9, pointsize=12, res=72*20)

par(mar = c(3, 3, 2, 2))

map('usa', xlim=c(-87.75, -75.25), ylim=c(23, 35.75), col=1)
image(x[seq(1, 3360, 5)], y[seq(1, 3360, 5)], z[seq(1, 3360, 5), seq(1, 3360, 5)], col=cols, axes=T, xlab="", ylab="", add=T)
axis(1); axis(2, las = 2); box()
mtext(side = 3, line = 1, paste(model, "-", spp), font = 2, cex = 1.2) 

i <- 1:length(cod)   
k <- i[seq(1, length(i), length.out = 100)] # change to 300
  for (j in k) {  lines(lon[,j]-360, lat[,j], col="#00000030") }

i <- which(cod == (-4))
k <- i[seq(1, length(i), length.out=100)]  
    for (j in k) {
      lo <- lon[!is.na(lon[,j]),j][length(lon[!is.na(lon[,j]),j])] - 360
      la <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]
          lines(lon[,j]-360, lat[,j], col = "#FF00FF30")
          points(lon[1,j]-360, lat[1,j], col="#00000050", pch=19, cex=1)
          if (lon[1,j]-360 < (-81.7))  {  points(lo, la, col="#FF000050", pch=19, cex=1.2)  }  else { 
                                          points(lo, la, col="#FFFF0050", pch=19, cex=1.2)  }
         }
legend("topleft", c("release locations", "settlement locations of GoM-spawned larvae", "settlement locations of Atl-spawned larvae", 
                    "all trajectories", "trajectories of successful recruits"), 
 lty=c(0,0,0,1,1), lwd=c(0,0,0,2,2), col=c(1, 2, 7, 1, 6), pch=c(19, 19, 19, -1, -1), cex=1, y.intersp=1.0, bty="n")  

dev.off()


######################   CONNECTIVITY MATRIX   ##########################

rm(list = ls()[-which(ls() %in% c("model", "spp", "folder"))])
library(maps)
library(matlab)

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

dev.off()
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

png(filename = paste0(folder, "conmatAllYears-", model, "-", spp, ".png"), units="in", width=5, height=5.7, pointsize=12, res=72*20)

nf <- layout(matrix(c(1:2), 2, 1), c(10), c(10, 1.4))   # plot layout -- all years
layout.show(nf)

par(mar=c(4,4,2,2), xpd=F)
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
image(pollis[42:limit], pollis[42:limit], log(settle[42:limit,42:limit]+1), col=cl, axes=F, 
      xlab="", ylab="", main = paste(model, "-", spp))

axis(1, at=tklocs, lab=rep("", 6))
axis(2, at=tklocs, lab=rep("", 6))
axis(2, at=lablocs, lab=c("W GOM", "E GOM", "E FL", "GA", "SC", "NC"), tick=F, cex=1, las=2)
axis(1, at=lablocs, lab=c("W GOM", "E GOM", "E FL", "GA", "SC", "NC"), tick=F, cex=1, las=1)
box(); abline(0,1, lty=2)
abline(v=tklocs[2], lty=1); abline(h=tklocs[2], lty=1)

mtext(side=2, line=3.0, "Source Node", font=2)
mtext(side=1, line=2.25, "Receiving Node", font=2)

####  color bar plot
par(mar=c(0,2,0,0))   
plot(c(-0.5,1.5), c(1.05,1.18), axes=F, xlab="", ylab="", col="white")
bra <- (seq(-0.5, ceiling(log(max(settle))), (ceiling(log(max(settle)))+0.5)/nn))
i = seq(0,1,1/(nn-1))
rect(i, 1.1, 1.04, 1.12, col = cl, lwd=0, border=cl)
laa <- c(1,10,100)   # laa <- seq(0, 60000, 10000)
levs <- (0+i[1:nn]+1/(nn*2))
re <- lm(bra[1:(length(bra)-1)]~levs)
text((log(laa) - re$coef[1])/re$coef[2], 1.07, laa, pos=4, cex=1.2)
text(0.5,1.12, cex=1.2, "number of successful recruits", pos=3)

dev.off()


#########  plot layout -- individual years

png(filename = paste0(folder, "conmatByYear-", model, "-", spp, ".png"), units="in", width=13, height=3, pointsize=12, res=72*20)

nf <- layout(matrix(c(1:7, 0, rep(8, 5), 0), 2, 7, byrow=TRUE), c(2, rep(6,5), 3), c(6, 2))
layout.show(nf)

#nf <- layout(matrix(c(1:4, 1, 5:7, 1, 8, 9, 9, rep(10,4)), 4, 4, byrow=TRUE), c(1,rep(7,3)), c(7,7,7,1.5))
#layout.show(nf)
#nf <- layout(matrix(c(1:5, 14, 1, 6:9, 14, 1, 10:13, 14, 1, rep(15,5)), 4, 6, byrow=TRUE), c(1,rep(7,4), 3), c(7,7,7,1.5))
#layout.show(nf)

####  global y label
par(mar=c(0,0,0,3))
plot.new()
mtext(side=2, line=-4, "Source Node", font=2)

#### plot connectivity matrices by year
resmat <- matrix(NA, nrow=13, ncol=5)
limit <- length(pollis)

for (k in unique(con$yr))  {
  con1 <- con[which(con$yr == k),]

par(mar=c(3, 3, 2, 2), xpd=F)
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
      main= paste(model, "-", spp, "-", k))

axis(1, at=tklocs, lab=rep("", 6))
axis(2, at=tklocs, lab=rep("", 6))
axis(2, at=lablocs, lab=c("W GOM", "E GOM", "E FL", "GA", "SC", "NC"), tick=F, cex=1, las=2)
axis(1, at=lablocs, lab=c("W GOM", "E GOM", "E FL", "GA", "SC", "NC"), tick=F, cex=1, las=2)
box(); abline(0,1, lty=2)/
abline(v=tklocs[2], lty=1); abline(h=tklocs[2], lty=1)
}

if (length(unique(con$yr)) <= 4)  { plot.new() }
if (length(unique(con$yr)) <= 3)  { plot.new() }
if (length(unique(con$yr)) <= 2)  { plot.new() }

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
text(1.11, poss, laa, pos=4, cex=1.2)
text(1.115, -0.5, cex=1.1, "number of\nsuccessful\nrecruits", pos=3)
#legend("bottom", "self-recruitment", lwd=1, lty=2, bty="n", cex=1.1)

###  global x label
par(mar=c(0,0,3,0))
plot.new()
mtext(side=1, line=-4, "Receiving Node", font=2)
#
dev.off()


# CONNECTIVITY MATRIX - regional summary ---------------------

rm(list = ls()[-which(ls() %in% c("model", "spp", "folder"))])

##############  concatenate files  ###################
filelist <- list.files(path = ".", pattern="con_file")     #   find files
dat <- read.table(filelist[1])
filelist <- filelist[-1]
for (i in filelist)  { newdat <- read.table(i); dat <- rbind(dat, newdat) }
colnames(dat) <- c("rel_poly","ret_poly","ret_yr","ret_mo","ret_d","age","ret_dep","rel_yr","rel_mo","rel_d")
#########################################################

d <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/CMS_input_files/redSnapperSett_GOM_ATL_hires.xyz", sep="\t")
names(d) <- c("lon", "lat", "pol")

for (i in 1) {
pa <- 46:52
bb <- 53:60
ta <- 61:66
sw <- 67:71
dt <- 72:76
sf <- 77:79
ef <- 80:89
ga <- 90:95
sc <- 96:106
nc <- 107:117
labs <- c("Panhandle", "Big Bend", "Central W FL", "SW FL", "Dry Tortugas", "S Florida", "E Florida", "GA", "SC", "NC")

dat$rel <- NA
dat$ret <- NA
dat$rel[which(dat$rel_poly %in% pa)] <- 1
dat$rel[which(dat$rel_poly %in% bb)] <- 2
dat$rel[which(dat$rel_poly %in% ta)] <- 3
dat$rel[which(dat$rel_poly %in% sw)] <- 4
dat$rel[which(dat$rel_poly %in% dt)] <- 5
dat$rel[which(dat$rel_poly %in% sf)] <- 6
dat$rel[which(dat$rel_poly %in% ef)] <- 7
dat$rel[which(dat$rel_poly %in% ga)] <- 8
dat$rel[which(dat$rel_poly %in% sc)] <- 9
dat$rel[which(dat$rel_poly %in% nc)] <- 10

dat$ret[which(dat$ret_poly %in% pa)] <- 1
dat$ret[which(dat$ret_poly %in% bb)] <- 2
dat$ret[which(dat$ret_poly %in% ta)] <- 3
dat$ret[which(dat$ret_poly %in% sw)] <- 4
dat$ret[which(dat$ret_poly %in% dt)] <- 5
dat$ret[which(dat$ret_poly %in% sf)] <- 6
dat$ret[which(dat$ret_poly %in% ef)] <- 7
dat$ret[which(dat$ret_poly %in% ga)] <- 8
dat$ret[which(dat$ret_poly %in% sc)] <- 9
dat$ret[which(dat$ret_poly %in% nc)] <- 10
}

##########  plot new grid cells  #########################

pols <- 1:10

dev.off()
barplot(table(dat$ret))
barplot(table(dat$rel))
barplot(table(dat$ret_yr))
barplot(table(dat$ret_mo))
barplot(table(dat$ret_dep))
barplot(table(dat$age))

# plot for all years together ----------------

png(filename = paste0(folder, "conmatLoRes-", model, "-", spp, ".png"), units="in", width=6, height=6, pointsize=12, res=72*20)

con <- dat
colnames(con) <- c("rel", "ret", "yr","mon","day","tim","ret_dep","rel_yr","rel_mo","rel_d", "source", "sink")
pollis <- 1:10    # list of polygon names
limit <- length(pollis)

#dev.off()
par(mar=c(8, 7, 2, 1), xpd=F)
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
image(pollis, pollis[1:9], (settle[, 1:9]), col=cl, axes=F,  xlab="", ylab="", 
      main = paste(model, "-", spp))

axis(1, at=1:10, lab=labs, las = 2)
axis(2, at=1:10, lab=labs, las = 2)
box()
#abline(0,1, lty=2)
abline(v = 5.5, lty=1, col = 0)
abline(h = 5.5, lty=1, col = 0)

txt <- round(settle / sum(sum(settle)) * 100, 1)
for (i in 1:10) {
  for (j in 1:10)  {
    if (txt[i,j] != 0) {
    text(i, j, txt[i,j], col = 0) }}}

mtext(side=2, line=6.2, "Source Node", font=2)
mtext(side=1, line=6.2, "Receiving Node", font=2)

dev.off()

# double check calcs

ga <- sum(as.vector(txt[6:10, 1:5]))
aa <- sum(as.vector(txt[6:10, 6:10]))

ga / (ga + aa)  # should equal as above


