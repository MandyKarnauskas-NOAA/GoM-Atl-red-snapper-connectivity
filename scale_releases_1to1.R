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

### Step 1: import independent maps and scale them to known ratios calculated above

load("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/GOMreleaseForScaling_SABGOM.RData")
GOM <- matfinGOM
load("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/ATLreleaseForScaling_SABGOM.RData")
ATL <- matS

# check how many millions of particles present in each file
sum(GOM$V5)/10^6
sum(ATL$V5)/10^6

areaGOM <- length(unique(paste(GOM$V2, GOM$V3)))
areaATL <- length(unique(paste(ATL$V2, ATL$V3)))  

arearat <- areaGOM/areaATL

# calculate constants for scaling 
const_scaler <- sum(GOM$V5)/sum(ATL$V5)
const_scaler 

unscaled <- rbind(GOM, ATL)
dim(unscaled)

GOMsc <- GOM
ATLsc <- ATL

GOMsc$V5 <- round(GOMsc$V5 * 10)
ATLsc$V5 <- round(ATLsc$V5 * 10 * const_scaler / arearat)

scaled <- rbind(GOMsc, ATLsc)
dim(scaled)
which(scaled$V5==0)

par(mfrow=c(1,2), mex=0.8)
barplot(c(sum(GOM$V5), sum(ATL$V5)), names.arg=c("GoM", "Atl"), main="unscaled", ylab=""); abline(0,0)
#barplot(c(fecGOM,  fecATL), names.arg=c("GoM",  "Atl"), main="specified ratio", ylab=""); abline(0,0)
#barplot(c(sum(GOMsc$V5)/areaGOM, sum(KEYsc$V5)/areaKEY, sum(ATLsc$V5)/areaATL), names.arg=c("GoM", "Keys", "Atl"), main="scaled by area", ylab="relative total population fecundity"); abline(0,0)
barplot(c(sum(GOMsc$V5), sum(ATLsc$V5)), names.arg=c("GoM", "Atl"), main="scaled", ylab=""); abline(0,0)

# check ratios
(sum(GOMsc$V5)) / (sum(ATLsc$V5))
arearat
(sum(GOM$V5)) / (sum(ATL$V5))

#  plot results to check
par(mfrow=c(1,2))
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

#########################  reduce number of particles  ###########################
min(scaled$V5); max(scaled$V5); mean(scaled$V5); sum(scaled$V5)
dim(scaled)

scaled$V5 <- round(scaled$V5/50)

min(scaled$V5); max(scaled$V5); mean(scaled$V5); sum(scaled$V5)
table(scaled$V5==0)
scaled <- scaled[(scaled$V5>0),]

dim(scaled)
min(scaled$V5); max(scaled$V5); mean(scaled$V5); sum(scaled$V5)

############################  check particle reduction  ###########################
sc2 <- scaled[which(scaled$V6==2004 & scaled$V7==6 & scaled$V8==17),]
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

scaledfin <- scaled
dim(scaled)

for (y in 2005:2010)  {
  scaled2 <- scaled
  scaled2$V6 <- y
  scaledfin <- rbind(scaledfin, scaled2)   }

dim(scaledfin)/7
dim(scaled)
sum(scaledfin$V5)


###########################  save output  ####################################

write.table(scaledfin, file="scaledGOMATLrel20042010.txt", sep="\t", col.names=F, row.names=F)

################################  END  #######################################
