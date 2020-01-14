################################################################################
############   GoM-ATL scaled red snapper release file for CMS  ################
############   M. Karnauskas Sep 12, 2019                       ################
#
#  code takes output from independent regional maps for GoM and Atl
#  scales release files according to ratios reported in surveys done in both basins
#  for use in full GoM-Atl simulation
################################################################################

rm(list=ls())
if (!"maps" %in% installed.packages()) install.packages("maps", repos='http://cran.us.r-project.org')
library(maps)

###  Step 1: calculate the fecundity of "average fish" in GoM vs Atl

yr <- 2015    # set reference year
ages <- 1:20
BF <- 1.732*(1-exp(-0.29*ages)) ^ 6.047
plot(ages, BF)

#####  Atlantic
NAAatl <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/SA.Nage.csv", sep=",", header=T)
NAAatl[which(NAAatl$year == yr),1]                  # check ref year
vecATL <- NAAatl[which(NAAatl$year == yr),2:21]        # extract NAA vector
BFatl <- sum(vecATL*BF)/sum(vecATL)                       # average fecundity per fish in Atl

#####  Gulf
NAAgom <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/Gulf_NAA_SEDAR52.csv", sep=",", header=T, skip=1)
NAAgom[which(NAAgom$Time == yr & NAAgom$Area==1),9]                    # check ref year
vecGOM <- NAAgom[which(NAAgom$Time == yr & NAAgom$Area==1), 13:32]        # extract NAA vector; caution with 2-area model
BFgom <- sum(vecGOM*BF)/sum(vecGOM)                                          # average fecundity per fish in eastern GoM

vec <- data.frame(t(rbind(vecATL/sum(vecATL), vecGOM/sum(vecGOM))))
barplot(t(vec), beside=T, legend=c("GOM", "ATL"), col=c(2,3), args.legend=list(x="top"), xlab="age classes", ylab="relative abundance")

BFgom     # weighted average - relative fecundity of a fish in the GOM
BFatl     # weighted average - relative fecundity of a fish in the ATL

###  Step 2: calculate relative fecundity present in Gulf, Atl

# input order: GULF, ATL           # REMEMBER: this is on a per-area basis! 
ratio <-      c(1.26, 1.83)        # input ratio of abundance from independent fisheries surveys
                                  
# these measures should be considered as representing per unit-area basis
fecGOM <- BFgom * ratio[1]        # relative fecundity per unit area GOM
fecATL <- BFatl * ratio[2]        # relative fecundity per unit area ATL

par(mfrow=c(2,1))
barplot(ratio, names.arg=c("GoM", "Atl"), main="relative population abundance per unit area\n(ratios reported from Ted)"); abline(0,0)
barplot(c(fecGOM, fecATL), names.arg=c("GoM", "Atl"), main="relative fecundity per unit area"); abline(0,0)

### Step 3: import independent maps and scale them to known ratios calculated above

load("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/GOMreleaseForScaling_SABGOM.RData")
GOM <- matfinGOM
load("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/ATLreleaseForScaling_SABGOM.RData")
ATL <- matS

# check how many millions of particles present in each file
sum(GOM$V5)/10^6
sum(ATL$V5)/10^6

areaGOM <- length(unique(paste(GOM$V2, GOM$V3)))

# Atlantic is not entire map -- subset grid cells within sampling
ATL2 <- ATL[which(ATL$V3 >= 28.0069 & ATL$V3 <= 30.5296), ]
dim(ATL); dim(ATL2)
areaATL <- length(unique(paste(ATL2$V2, ATL2$V3)))  

areaGOM; areaATL   # areas in terms of number of 10^km cells

# calculate constants for scaling 
# solve for constant based on relative fecundity per unit area, relative fecundity metric within regional map, and area of regional map
const_scaler <- 10^5    # set so that final output has reasonable number of particles

constGOM <- areaGOM * fecGOM / sum(GOM$V5) * const_scaler      # check this calculation later!!!!!
constATL <- areaATL * fecATL / sum(ATL2$V5) * const_scaler

constGOM
constATL

unscaled <- rbind(GOM, ATL)
dim(unscaled)

GOMsc <- GOM
ATLsc <- ATL

GOMsc$V5 <- round(GOMsc$V5 * constGOM)
ATLsc$V5 <- round(ATLsc$V5 * constATL)

scaled <- rbind(GOMsc, ATLsc)
dim(scaled)
which(scaled$V5==0)

par(mfrow=c(1,3), mex=0.8)
barplot(c(sum(GOM$V5), sum(ATL$V5)), names.arg=c("GoM", "Atl"), main="unscaled", ylab=""); abline(0,0)
barplot(c(fecGOM,  fecATL), names.arg=c("GoM",  "Atl"), main="specified ratio", ylab=""); abline(0,0)
#barplot(c(sum(GOMsc$V5)/areaGOM, sum(KEYsc$V5)/areaKEY, sum(ATLsc$V5)/areaATL), names.arg=c("GoM", "Keys", "Atl"), main="scaled by area", ylab="relative total population fecundity"); abline(0,0)
barplot(c(sum(GOMsc$V5), sum(ATLsc$V5)), names.arg=c("GoM", "Atl"), main="scaled", ylab=""); abline(0,0)

# check ratios

ratio                # ratio reported from Ted
ratio[1]/ratio[2]
(sum(GOMsc$V5)) / (sum(ATLsc$V5))
(sum(GOMsc$V5)/areaGOM) / (sum(ATLsc$V5)/areaATL)


#  plot results to check
windows()
sc2 <- scaled[which(scaled$V6==2004 & scaled$V7==5 & scaled$V8==12),]
x <- sc2$V5
pos <- c(0.005,  0.05, 0.1, 0.2, 0.5, 1, 2, 4, 5, 10, 50, 100, 200, 500, 1000)
a <- floor(min(x))
b <- max((x)-a)*1.03
pind <- round((x-a)/b*100+1); print(min(pind)); print(max(pind))
cols <- c(rainbow(30, start=0.82, end=0.99), rainbow(70, start=0.01, end=0.17))[100:1]
map('state', fill = 1, interior=F, col = gray(0.95), ylim=c(23.0, 35), xlim=c(-88,-76))
points(sc2$V2, sc2$V3, col=cols[pind], pch=15, cex=0.8)
box(); axis(1); axis(2)
xloc <- seq(-86, -78, length.out=100)
for (j in 1:100) {   polygon(c(xloc[j], xloc[j+1],xloc[j+1], xloc[j]), c(23.5,23.5,24.0,24.0), col=cols[j], border=NA) }
w <- which.min(abs(((max(x)-min(x))/6) - pos))
if(-pos[w]<min(x)) { xx <- seq(0, max(x), pos[w]); xx <- xx[xx>min(x)] } else {  xx <- c(seq(-pos[w], min(x), -pos[w]), seq(0, max(x), pos[w])) }
text(xloc[round((xx-a)/b*100+1)], y=23.2, xx, pos=2)
mtext(side=3, line=1.5, "scaled fecundity map", cex=1.3, font=2)
lines(x=c(-81.5, -79.5), rep(28.0069,2))
lines(x=c(-81.5, -79.5), rep(30.5296,2))
lines(x=rep(-81.5, 2), c(28.0069, 30.5296))
lines(x=rep(-79.5, 2), c(28.0069, 30.5296))


#  plot results to check
windows()
sc3 <- unscaled[which(unscaled$V6==2004 & unscaled$V7==5 & unscaled$V8==12),]
x <- sc3$V5
pos <- c(0.005,  0.05, 0.1, 0.2, 0.5, 1, 2, 4, 5, 10, 50, 100, 500, 1000)
a <- floor(min(x))
b <- max((x)-a)*1.03
pind <- round((x-a)/b*100+1); print(min(pind)); print(max(pind))
cols <- c(rainbow(30, start=0.82, end=0.99), rainbow(70, start=0.01, end=0.17))[100:1]
map('state', fill = 1, interior=F, col = gray(0.95), ylim=c(23.0, 35), xlim=c(-88,-76))
points(sc2$V2, sc2$V3, col=cols[pind], pch=15, cex=0.8)
box(); axis(1); axis(2)
xloc <- seq(-86, -78, length.out=100)
for (j in 1:100) {   polygon(c(xloc[j], xloc[j+1],xloc[j+1], xloc[j]), c(23.5,23.5,24.0,24.0), col=cols[j], border=NA) }
w <- which.min(abs(((max(x)-min(x))/6) - pos))
if(-pos[w]<min(x)) { xx <- seq(0, max(x), pos[w]); xx <- xx[xx>min(x)] } else {  xx <- c(seq(-pos[w], min(x), -pos[w]), seq(0, max(x), pos[w])) }
text(xloc[round((xx-a)/b*100+1)], y=23.2, xx, pos=2)
mtext(side=3, line=1.5, "unscaled fecundity map", cex=1.3, font=2)
lines(x=c(-81.5, -79.5), rep(28.0069,2))
lines(x=c(-81.5, -79.5), rep(30.5296,2))
lines(x=rep(-81.5, 2), c(28.0069, 30.5296))
lines(x=rep(-79.5, 2), c(28.0069, 30.5296))







