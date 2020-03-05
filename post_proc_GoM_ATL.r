rm(list=ls())

setwd("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINAL_RUNS/SABGOM1")

##############  concatenate files  ###################                                   
filelist <- list.files(path = ".", pattern="con_file")     #   find files
dat <- read.table(filelist[1])
filelist <- filelist[-1]
for (i in filelist)  { newdat <- read.table(i); dat <- rbind(dat, newdat) }
colnames(dat) <- c("rel_poly","ret_poly","ret_yr","ret_mo","ret_d","age","ret_dep","rel_yr","rel_mo","rel_d")
##########################################################

co <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/CMS_input_files/redSnapperSett_GOM_ATL_hires.xyz", sep="\t")
map("state", ylim=c(24, 36), xlim=c(-100, -75))  # check
for (i in 1:117)  {
  f <- co[which(co$V3==i),]
  polygon(f, col=i)  
  text(mean(f$V1), mean(f$V2), mean(f$V3))}

dat$reg <- "GOM"
dat$reg[which(dat$ret_poly > 77)] <- "ATL"

dat$rel_reg <- "GOM"
dat$rel_reg[which(dat$rel_poly > 77)] <- "ATL"

###########################  specify release file  #############################

d <- read.table("C://Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/scaledGOMATLrel20042010.txt")

d$reg <- "GOM"
d$reg[which(d$V1 > 77)] <- "ATL"

tapply(d$V5, d$reg, sum)     #   ATL    GOM 
                            # 615783 677733
                            
length(which(dat$rel_poly > 77 & dat$ret_poly  > 77))    # 517049 = spawned in ATL, recruited in ATL
length(which(dat$rel_poly <= 77 & dat$ret_poly > 77))    # 8080   = spawned in GOM, recruited in ATL
length(which(dat$rel_poly <= 77 & dat$ret_poly <= 77))   # 291605 = spawned in GOM, recruited in GOM
length(which(dat$rel_poly > 77 & dat$ret_poly  <= 77))   # 72     = spawned in ATL, recruited in GOM

length(which(dat$rel_poly > 77 & dat$ret_poly > 77))/sum(d$V5)     # 0.8394142    = spawned in ATL, recruited in ATL
length(which(dat$rel_poly <= 77 & dat$ret_poly > 77))/sum(d$V5)    # 0.01136878   = spawned in GOM, recruited in ATL
length(which(dat$rel_poly <= 77 & dat$ret_poly <= 77))/sum(d$V5)   # 0.4308186    = spawned in GOM, recruited in GOM
length(which(dat$rel_poly > 77 & dat$ret_poly <= 77))/sum(d$V5)    # 0.0003637645 = spawned in ATL, recruited in GOM

0.8394142/(0.8394142 + 0.01136878)  # = 0.9866373        98.66 % 

rel <- tapply(d$V5, list(d$V6, d$reg), sum)

surv <- table(dat$ret_yr, dat$reg)

surv/rel

#dat <- dat[which(dat$rel_poly>64),]

barplot(table(dat$reg, dat$ret_yr), beside=T)

table(dat$reg, dat$rel_poly)         # polygon 62 is first to have successful recruitment to S Atl

(table(dat$reg) / nrow(dat)) * 100

a <- table(dat$reg, dat$ret_yr)
a

########################  END BASIC COMPARISONS  ###############################

########################  CHECK PARTICLE BEHAVIORS  ############################
rm(list=ls())

library(ncdf4)

getwd()

nc <- nc_open("traj_file_2.nc")

v1 <- nc$var[[1]]
int <- ncvar_get(nc, v1)
v2 <- nc$var[[2]]
lon <- ncvar_get(nc, v2)
v3 <- nc$var[[3]]
lat <- ncvar_get(nc, v3)
v4 <- nc$var[[4]]
dep <- ncvar_get(nc, v4)
v5 <- nc$var[[5]]
dis <- ncvar_get(nc, v5)
v6 <- nc$var[[6]]
cod <- ncvar_get(nc, v6)
v7 <- nc$var[[7]]
dates <- ncvar_get(nc, v7)
v8 <- nc$var[[8]]
pol <- ncvar_get(nc, v8)
v9 <- nc$var[[9]]
dens <- ncvar_get(nc, v9)
v10 <- nc$var[[10]]
diam <- ncvar_get(nc, v10)
nc_close(nc)

# check!  should have only -4, -1, and 0
table(cod)                                   # 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data
dates[which(cod==(-1))]
lon[which(cod==(-1))]
dep[2,which(cod==(-1))]
hist(dep[2,which(cod==(-1))])
hist(dep[2,])

windows()
par(mfrow=c(6,5), mex=0.6)
for (i in 1:30) { hist(dep[i,], breaks=seq(0,100,5)) }

table(dep[2,]==0)/ncol(dep)
table(dep[2,]<10)/ncol(dep)
table(dep[2,]>10 & dep[2,]<20)/ncol(dep)
table(dep[2,]>20 & dep[2,]<30)/ncol(dep)

windows()
matplot(-dep[,seq(1, ncol(dep), length.out=100)], col="#FF00FF30", pch=19, type="l")
matplot(-dep[,seq(1, ncol(dep), length.out=100)], col="#FF00FF30", pch=19, type="p")

plot(dens, -dep[2,], col="#FF00FF30")
plot(diam, -dep[2,], col="#FF00FF30")
 
cc <- table(cut(dep[2,], breaks=c(-1,0.5, 10, 20, 30, 40)))
cc
cc/sum(cc) 
barplot(cc)
 
min(dens)
max(dens)
 
mean(dep[2,])
min(dep[2,])
max(dep[2,])

##########################  PLOT TRAJECTORIES  ###############################

rm(list=ls())

getwd()

library(ncdf4)
library(maps)

##########  bathymetry  ##############

#save(x,y,z,cols, file="GEBCO_bathy.RData")
load("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/GEBCO_bathy.RData")        # stored GEBCO data and color scheme

###############  plots by year  
yrs <- c(2008)

par(mfrow=c(3,1), mex=0.6)

maxlat <- rep(NA, 7)

numlarv <- 500

for (m in 1:1) {
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
  cod <- ncvar_get(nc, v6)  # 0 can still move; -1 left area; -2 close to land; -3 dead; -4 settled; -5 no oceanographic data
nc_close(nc)
}

#  i <- which(cod == (-4))
#  la <- rep(NA, length(i))
#      for (j in i) {   la[j] <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]  } 
#  maxlat[m] <- mean(la, na.rm=T)     }

map('usa', xlim=c(-98, -76), ylim=c(24,35), col=1)
#map('usa', xlim=c(-85.5, -75), ylim=c(23.85,35), col=0)
#image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1); axis(2, las=2)
legend("topleft", paste(yrs[m]), cex=1.2, text.font=2, bty="n")
  i <- which(cod == (-4))
  k <- i[seq(1, length(i), length.out=numlarv*2)]
  for (j in k) {  lines(lon[,j]-360, lat[,j], col="#FF00FF20") } 
  
  
  
  i <- which(cod == (-4))
  k <- i[seq(1, length(i), length.out=numlarv)]  
    for (j in k) {
      lo <- lon[!is.na(lon[,j]),j][length(lon[!is.na(lon[,j]),j])] - 360
      la <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]
          lines(lon[,j]-360, lat[,j], col="#00000040")
          points(lon[1,j]-360, lat[1,j], col="#00000050", pch=19, cex=1)
          if (lon[1,j]-360 < (-81.7))  {  points(lo, la, col="#00FF0090", pch=19, cex=1.2)  }  else { 
          points(lo, la, col="#FF000050", pch=19, cex=1)  }
         }      
         }
#  i <- which(cod == (-4))    
#    for (j in i) {
#      lo <- lon[!is.na(lon[,j]),j][length(lon[!is.na(lon[,j]),j])] - 360
#      la <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]
#        if (lo > (-81.5))  {
#          lines(lon[,j]-360, lat[,j], col="#00000010")
#          points(lon[1,j]-360, lat[1,j], col="#00000050", pch=19, cex=1)
#          points(lo, la, col="#FF000070", pch=19, cex=1.2)
#       }  }
#    }
plot(1,1, axes=F, col=0, xlab="", ylab="")
#legend("center", c("all trajectories", "trajectories of successful \nrecruits to S ATL", "release locations of \nsuccessful S ATL recruits", "settlement locations \nof successful recruits"), 
#lty=c(1,1,0,0), lwd=c(2,2,0,0), col=c("yellow", 1, 1, 2), pch=c(-1, -1, 19, 19), cex=1.2, y.intersp=1.75, bty="n")  

#legend("center", c("all trajectories", "trajectories of \nsuccessful recruits", "release locations", "settlement locations of \nlarvae released from GoM", "settlement locations of \nlarvae released from Atl"), 
#lty=c(1,1,0,0,0), lwd=c(2,2,0,0,0), col=c("yellow", 1, "#00000050", "#00FF0090", "#FF000050"), pch=c(-1, -1, 19, 19, 19), cex=1.2, y.intersp=1.75, bty="n", pt.cex=2)  

legend("center", c("all trajectories", "trajectories of \nsuccessful recruits", "release locations", "settlement locations of \nlarvae released from Atl"), 
lty=c(1,1,0,0), lwd=c(2,2,0,0), col=c("yellow", 1, "#00000050", "#FF000050"), pch=c(-1, -1, 19, 19), cex=1.2, y.intersp=1.75, bty="n", pt.cex=2)  



    
### all trajectories on one map    
map('usa', xlim=c(-84, -76), ylim=c(23.85,35), col=0)
image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1); axis(2, las=2)

for (m in 1:3) {
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
nc_close(nc)

  i <- 1:length(cod)   #which(cod == (-4))
  k <- i[seq(1, length(i), length.out=2000)]
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
    }
#legend("topleft", c("all trajectories", "trajectories of successful recruits", "release locations", "settlement locations of \nlarvae released from GoM", "settlement locations of \nlarvae released from Atl"), 
# lty=c(1,1,0,0,0), lwd=c(2,2,0,0,0), col=c("yellow", 1, 1, 3, 2), pch=c(-1, -1, 19, 19, 19), cex=1.2, y.intersp=1.5, bty="n")  
legend("topleft", c("all trajectories", "trajectories of successful recruits", "release locations", "settlement locations of \nlarvae released from Atl"), 
 lty=c(1,1,0,0), lwd=c(2,2,0,0), col=c("yellow", 1, 1, 2), pch=c(-1, -1, 19, 19), cex=1.0, y.intersp=1.0, bty="n")  


        # jun 12
 ### all trajectories on one map  -- full GoM   
par(mar=c(7,7,1,1))
map('usa', xlim=c(-98, -76), ylim=c(23.84,35), col=0, xlab="longitude", ylab="latitude")
image(x,y,z, col=cols, axes=T, add=T, xlab="longitude", ylab="latitude"); box(); axis(1); axis(2, las=2)
#mtext(side=1, line=2, "longitude"); mtext(side=2, line=2.5, "latitude")
mtext(side=3, line=1, "HYCOM 1/50 degree model -- 2008", font=2)
#mtext(side=3, line=1, "SABGOM model -- 2004 - 2010", font=2)


for (m in 1:3) {
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
nc_close(nc)

  i <- 1:length(cod)   #  which(cod == (-4))        # 
  k <- i[seq(1, length(i), length.out=6000)]
  for (j in k) {  lines(lon[,j]-360, lat[,j], col="#FFFF0010") } 
                            i <- which(cod == (-4))  
                            k <- i[seq(1, length(i), length.out=6000)]
  for (j in k) { #which(cod == (-4))  ) {
      lo <- lon[max(which(!is.na(lon[,j]))),j] - 360
      la <- lat[max(which(!is.na(lat[,j]))),j]
    if (lo > (-81.7))  {  
          lines(lon[,j]-360, lat[,j], col="#00000010")
          points(lon[1,j]-360, lat[1,j], col="#00FF0050", pch=19, cex=1)   
          points(lo, la, col="#FF000050", pch=19, cex=1) 
             } }
    }
legend("topleft", c("all larval trajectories", "trajectories of larvae successfully recruiting to Atlantic", "source locations of larvae recruiting to Atlantic", "settlement locations of larvae released from Gulf"), 
 lty=c(1,1,0,0), lwd=c(2,2,0,0), col=c("yellow", 1, "green", 2), pch=c(-1, -1, 19, 19), cex=1.0, y.intersp=1.2, bty="n")  



#########################

map('usa', xlim=c(-83, -76), ylim=c(25,36), col=0)
image(x,y,z, col=cols, axes=T, xlab="", ylab="", add=T); box(); axis(1); axis(2, las=2)
i <- which(cod == (-4))
k <- i[seq(1, length(i), length.out=500)]
     for (j in k) {
      lo <- lon[!is.na(lon[,j]),j][length(lon[!is.na(lon[,j]),j])] - 360
      la <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]
          lines(lon[,j]-360, lat[,j], col="#00000010")
          points(lon[1,j]-360, lat[1,j], col="#00000010", pch=19, cex=0.5)
          points(lo, la, col="#FF000070", pch=19, cex=0.5)
       }  
 

 
######################   CONNECTIVITY MATRIX   ##########################

rm(list=ls())
getwd()

##############  concatenate files  ###################
filelist <- list.files(path = ".", pattern="con_file")     #   find files
dat <- read.table(filelist[1])
filelist <- filelist[-1]
for (i in filelist)  { newdat <- read.table(i); dat <- rbind(dat, newdat) }
colnames(dat) <- c("rel_poly","ret_poly","ret_yr","ret_mo","ret_d","age","ret_dep","rel_yr","rel_mo","rel_d")
##########################################################

# dat <- dat[which(dat$rel_poly>61),]

library(maps)
##############################################################################################

d <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/CMS_input_files/redSnapperSett_GOM_ATL_hires.xyz", sep="\t")
names(d) <- c("lon", "lat", "pol")

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

#pdf("13JunAtl_conmatSummary.pdf")
barplot(table(dat$ret_poly))
barplot(table(dat$rel_poly))
barplot(table(dat$ret_yr))
barplot(table(dat$ret_mo))
barplot(table(dat$ret_dep))
barplot(table(dat$age))
barplot(table(dat$state))
#dev.off()


#########################################################################################
#####################       connectivity matrix          ################################
#########################################################################################
# plot is for all years separately, see k loop to do all years together

library(matlab)

##########################

## connectivity file
con <- dat
colnames(con) <- c("source","sink","yr","mon","day","tim","ret_dep","rel_yr","rel_mo","rel_d", "state")

############## list of polygon names
pollis <- 1:max(d$pol)
#########  plot layout

nf <- layout(matrix(c(1:2), 2, 1), c(10,10), c(10, 2))
layout.show(nf)

#nf <- layout(matrix(c(1:4, 1, 5:7, 1, 8, 9, 9, rep(10,4)), 4, 4, byrow=TRUE), c(1,rep(7,3)), c(7,7,7,1.5))
#layout.show(nf)

#nf <- layout(matrix(c(1:5, 14, 1, 6:9, 14, 1, 10:13, 14, 1, rep(15,5)), 4, 6, byrow=TRUE), c(1,rep(7,4), 3), c(7,7,7,1.5))
#layout.show(nf)

####  global y label
#par(mar=c(0,0,0,0))
#plot.new()
#mtext(side=2, line=-1.5, "Source Node", font=2)
#### plot connectivity matrices by year
par(mar=c(3,3,3,2), xpd=F)
resmat <- matrix(NA, nrow=13, ncol=5)
limit <- length(pollis)

#for (k in 2004:2010)  {
#  con1 <- con[which(con$yr==k),]
#  con1 <- con[which(con$yr>=2009),]
   par(mar=c(5,5,3,2), xpd=F); con1 <- con                     # to look at all years together
  settle <- matrix(NA, length(pollis), length(pollis))
for (i in 1:length(pollis))     {
    a <- which(con1$source==pollis[i])
      if(length(a)>0) {
        for (j in 1:length(pollis)) {   settle[i,j] <- length(which(con1$sink[a]==j))    }
    } else { settle[i,] <- 0  }  }
settle <- t(settle)
nn <- 300
cl = jet.colors(nn)
#image(pollis[1:limit], pollis[1:limit], log((settle[1:limit,1:limit]+1)/(max(settle+2))), col=cl, axes=F,  xlab="", ylab="", main=k,
image(pollis[1:limit], pollis[1:limit], log(settle[1:limit,1:limit]+1), col=cl, axes=F,  xlab="", ylab="", 
main="GoM - Atl connectivity HYCOM 1/50 deg model - 2008")
#, xlim=c(57, 96), ylim=c(57,96)
#, main="GoM - Atl SABGOM connectivity 2004 - 2010")

axis(1, at=tklocs, lab=rep("", 6))
axis(2, at=tklocs, lab=rep("", 6))
#lablocs[2] <- 74
axis(2, at=lablocs, lab=c("W GOM", "E GOM", "E FL", "GA", "SC", "NC"), tick=F, cex=1, las=2)
axis(1, at=lablocs, lab=c("W GOM", "E GOM", "E FL", "GA", "SC", "NC"), tick=F, cex=1, las=2)
box(); abline(0,1, lty=2)
abline(v=tklocs[2], lty=1); abline(h=tklocs[2], lty=1)

mtext(side=2, line=3.5, "Source Node", font=2)
mtext(side=1, line=3.5, "Receiving Node", font=2)

#  text(79, 95, "N to S", font=2, col="white")
#  text(94, 78, "S to N", font=2, col="white")

#  text(limit-8,57, "GOM to ATL", font=2, col="white")
#  text(limit-8,limit-2, "ATL to ATL", font=2, col="white")
#  text(64,57, "GOM to GOM", font=2, col="white")
#  text(64,limit-2, "ATL to GOM", font=2, col="white")
  
#  text(limit-10,5, "GOM to ATL", font=2, col="white")
#  text(limit-10,limit-5, "ATL to ATL", font=2, col="white")
#  text(15,5, "GOM to GOM", font=2, col="white")
#  text(15,limit-5, "ATL to GOM", font=2, col="white")

# mtext(side=1, line=3.75, "Receiving node                ", font=2)
# mtext(side=2, line=3.75, "Source node                ", font=2)

#  resmat[k-2002,1] <- k
#  resmat[k-2002,2] <- sum(sum(settle[1:max(egom),1:max(egom)]))
#  resmat[k-2002,3] <- sum(sum(settle[1:max(egom),min(wfl):max(nc)]))
#  resmat[k-2002,4] <- sum(sum(settle[min(wfl):max(nc),1:max(egom)]))
#  resmat[k-2002,5] <- sum(sum(settle[min(wfl):max(nc),min(wfl):max(nc)]))   }

####  color bar plot
par(mar=c(0,2,0,0))   
plot(c(-0.5,1.5), c(1.05,1.18), axes=F, xlab="", ylab="", col="white")
bra <- (seq(-0.5, ceiling(log(max(settle))), (ceiling(log(max(settle)))+0.5)/nn))
i = seq(0,1,1/(nn-1))
#rect(1.05, i, 1.1, 1, col = cl, lwd=0, border=cl)
rect(i, 1.1, 1.04, 1.12, col = cl, lwd=0, border=cl)
laa <- c(1,10,100,1000)   # laa <- seq(0, 60000, 10000)
levs <- (0+i[1:nn]+1/(nn*2))
re <- lm(bra[1:(length(bra)-1)]~levs)
poss <- (log(laa) - re$coef[1])/re$coef[2]
text(poss, 1.09, laa, pos=4, cex=1.2)
text(0.5,1.12, cex=1.2, "number of successful recruits", pos=3)
#legend("bottom", "self-recruitment", lwd=1, lty=2, bty="n", cex=1.1)


###  global x label
par(mar=c(0,0,0,0))
plot.new()
mtext(side=1, line=-1.5, "Receiving Node", font=2)
#
