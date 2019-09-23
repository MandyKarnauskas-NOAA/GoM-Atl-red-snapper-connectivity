################################################################################
############   GoM-ATL scaled red snapper release file for CMS  ################
############   M. Karnauskas Sep 12, 2019                       ################
#
#  code takes output from independent regional maps for GoM and Atl
#  scales release files according to ratios reported in surveys done in both basins
#  for use in full GoM-Atl simulation
################################################################################

rm(list=ls())
library(maps)

###  Step 1: calculate the fecundity of "average fish" in GoM vs Atl

yr <- 2010    # set reference year
ages <- 1:20
BF <- 1.732*(1-exp(-0.29*ages)) ^ 6.047
plot(ages, BF)

#####  Atlantic
NAAatl <- read.table("C:/Users/mkarnauskas/Desktop/RS_FATEproject/SA.Nage.csv", sep=",", header=T)
NAAatl[which(NAAatl$year == yr),1]                  # check ref year
vec <- NAAatl[which(NAAatl$year == yr),2:21]        # extract NAA vector
BFatl <- sum(vec*BF)/sum(vec)                       # average fecundity per fish in Atl

#####  Gulf
NAAgom <- read.table("C:/Users/mkarnauskas/Desktop/RS_FATEproject/Gulf_NAA_SEDAR52.csv", sep=",", header=T, skip=1)
NAAgom[which(NAAgom$Time == yr & NAAgom$Area==1),9]                    # check ref year
vec <- NAAgom[which(NAAgom$Time == yr & NAAgom$Area==1), 13:32]        # extract NAA vector; caution with 2-area model
BFgom <- sum(vec*BF)/sum(vec)                                          # average fecundity per fish in eastern GoM

###  Step 2: calculate relative fecundity present in Gulf, Keys, Atl

# input order: GULF, KEYS, ATL
ratio <-      c(3, 0.4, 1)        # input ratio of abundance from independent fisheries survyes
                                  # REMEMBER: this is on a per-area basis!  
# these measures should be considered as representing per unit-area basis
fecGOM <- BFgom * ratio[1]        # relative fecundity per unit area GOM
fecKEY <- BFgom * ratio[2]        # relative fecundity per unit area FLK
fecATL <- BFatl * ratio[3]        # relative fecundity per unit area ATL

barplot(c(fecGOM, fecKEY, fecATL), names.arg=c("GoM", "Keys", "Atl"), ylab="relative total population fecundity"); abline(0,0)

### Step 3: import independent maps and scale them to known ratios calculated above

load("C:/Users/mkarnauskas/Desktop/RS_FATEproject/MASTER_codes/KEYSreleaseForScaling_SABGOM.RData")
KEY <- matfin
load("C:/Users/mkarnauskas/Desktop/RS_FATEproject/MASTER_codes/GOMreleaseForScaling_SABGOM.RData")
GOM <- matfinGOM
load("C:/Users/mkarnauskas/Desktop/RS_FATEproject/MASTER_codes/ATLreleaseForScaling_SABGOM.RData")
ATL <- matS

# check how many millions of particles present in each file
sum(GOM$V5)/10^6
sum(KEY$V5)/10^6
sum(ATL$V5)/10^6

areaGOM <- length(unique(paste(GOM$V2, GOM$V3)))
areaKEY <- length(unique(paste(KEY$V2, KEY$V3)))          

# Atlantic is not entire map -- subset grid cells within sampling
ATL2 <- ATL[which(ATL$V3 > 28 & ATL$V3 < 30.5), ]
dim(ATL); dim(ATL2)
areaATL <- length(unique(paste(ATL2$V2, ATL2$V3)))  

# calculate constants for scaling
const_scaler <- 10^8    # set so that final output has reasonable number of particles

constGOM <- areaGOM * fecGOM / sum(GOM$V5) * const_scaler      # check this calculation later!!!!!
constKEY <- areaKEY * fecKEY / sum(KEY$V5) * const_scaler
constATL <- areaATL * fecATL / sum(ATL2$V5) * const_scaler

constGOM
constKEY
constATL

unscaled <- rbind(GOM, KEY, ATL)
dim(unscaled)

GOMsc <- GOM
KEYsc <- KEY
ATLsc <- ATL

GOMsc$V5 <- round(GOMsc$V5 * constGOM)
KEYsc$V5 <- round(KEYsc$V5 * constKEY)
ATLsc$V5 <- round(ATLsc$V5 * constATL)

scaled <- rbind(GOMsc, KEYsc, ATLsc)
dim(scaled)
which(scaled$V5==0)

par(mfrow=c(4,1), mex=0.5)
barplot(c(fecGOM, fecKEY, fecATL), names.arg=c("GoM", "Keys", "Atl"), main="specified ratio", ylab="relative total population fecundity"); abline(0,0)
barplot(c(sum(GOMsc$V5)/areaGOM, sum(KEYsc$V5)/areaKEY, sum(ATLsc$V5)/areaATL), names.arg=c("GoM", "Keys", "Atl"), main="scaled by area", ylab="relative total population fecundity"); abline(0,0)
barplot(c(sum(GOMsc$V5), sum(KEYsc$V5), sum(ATLsc$V5)), names.arg=c("GoM", "Keys", "Atl"), main="scaled", ylab="relative total population fecundity"); abline(0,0)
barplot(c(sum(GOM$V5), sum(KEY$V5), sum(ATL$V5)), names.arg=c("GoM", "Keys", "Atl"), main="unscaled", ylab="relative total population fecundity"); abline(0,0)

# check ratios

ratio
ratio[1]/ratio[3]
(sum(GOMsc$V5)/areaGOM) / (sum(ATLsc$V5)/areaATL)

ratio[1]/ratio[2]
(sum(GOMsc$V5)) / (sum(KEYsc$V5))
(sum(GOMsc$V5)/areaGOM) / (sum(KEYsc$V5)/areaKEY)




#  plot results to check
sc2 <- scaled[which(scaled$V6==2004 & scaled$V7==5 & scaled$V8==12),]
x <- sc2$V5
pos <- c(0.005,  0.05, 0.1, 0.2, 0.5, 1, 2, 4, 5, 10, 50, 100, 500, 1000)
a <- floor(min(x))
b <- max((x)-a)*1.03
pind <- round((x-a)/b*100+1); print(min(pind)); print(max(pind))
cols <- c(rainbow(30, start=0.82, end=0.99), rainbow(70, start=0.01, end=0.17))[100:1]
map('state', fill = 1, interior=F, col = gray(0.95), ylim=c(23.4, 35), xlim=c(-88,-76))
points(sc2$V2, sc2$V3, col=cols[pind], pch=15, cex=0.8)
box(); axis(1); axis(2)
xloc <- seq(-86, -78, length.out=100)
for (j in 1:100) {   polygon(c(xloc[j], xloc[j+1],xloc[j+1], xloc[j]), c(24,24,24.3,24.3), col=cols[j], border=NA) }
w <- which.min(abs(((max(x)-min(x))/6) - pos))
if(-pos[w]<min(x)) { xx <- seq(0, max(x), pos[w]); xx <- xx[xx>min(x)] } else {  xx <- c(seq(-pos[w], min(x), -pos[w]), seq(0, max(x), pos[w])) }
text(xloc[round((xx-a)/b*100+1)], y=23.6, xx, pos=2)

#  plot results to check
sc3 <- unscaled[which(unscaled$V6==2004 & unscaled$V7==5 & unscaled$V8==12),]
x <- log(sc3$V5)
pos <- c(0.005,  0.05, 0.1, 0.2, 0.5, 1, 2, 4, 5, 10, 50, 100, 500, 1000)
a <- floor(min(x))
b <- max((x)-a)*1.03
pind <- round((x-a)/b*100+1); print(min(pind)); print(max(pind))
cols <- c(rainbow(30, start=0.82, end=0.99), rainbow(70, start=0.01, end=0.17))[100:1]
map('state', fill = 1, interior=F, col = gray(0.95), ylim=c(23.4, 35), xlim=c(-88,-76))
points(sc2$V2, sc2$V3, col=cols[pind], pch=19, cex=0.8)
box(); axis(1); axis(2)
xloc <- seq(-86, -78, length.out=100)
for (j in 1:100) {   polygon(c(xloc[j], xloc[j+1],xloc[j+1], xloc[j]), c(24,24,24.3,24.3), col=cols[j], border=NA) }
w <- which.min(abs(((max(x)-min(x))/6) - pos))
if(-pos[w]<min(x)) { xx <- seq(0, max(x), pos[w]); xx <- xx[xx>min(x)] } else {  xx <- c(seq(-pos[w], min(x), -pos[w]), seq(0, max(x), pos[w])) }
text(xloc[round((xx-a)/b*100+1)], y=23.6, xx, pos=2)







