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

source("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/findMinDepth.R")

### Step 1: import independent maps and scale them to known ratios calculated above

# load data ----------------------------------------------------------

load("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/GOMreleaseForScaling.RData")
GOM <- matfinGOM
load("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/ATLreleaseForScaling.RData")
ATL <- matS

# check how many millions of particles present in each file ----------

sum(GOM$V5)/10^5
sum(ATL$V5)/10^5

areaGOM <- length(unique(paste(GOM$V2, GOM$V3)))
areaATL <- length(unique(paste(ATL$V2, ATL$V3)))  
areaGOM
areaATL

arearat <- areaGOM/areaATL
arearat

# calculate constants for scaling -------------------------------------
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
barplot(c(sum(GOM$V5),   sum(ATL$V5)),   names.arg=c("GoM", "Atl"), main="unscaled", ylab=""); abline(0,0)
barplot(c(sum(GOMsc$V5), sum(ATLsc$V5)), names.arg=c("GoM", "Atl"), main="scaled", ylab=""); abline(0,0)

# check ratios -----------------------------------------
(sum(GOMsc$V5)) / (sum(ATLsc$V5))
arearat
(sum(GOM$V5)) / (sum(ATL$V5))

# plot results to check --------------------------------
#par(mfrow=c(1,2))

dev.off()

par(mar = c(5, 7, 1, 1))
sc2 <- scaled[which(scaled$V6==2013 & scaled$V7==6 & scaled$V8==27),]
x <- sc2$V5 / 130
pos <- c(0.005,  0.05, 0.1, 0.2, 0.5, 1, 2, 4, 5, 10, 20, 50, 100, 1000)
a <- floor(min(x))
b <- max((x)-a)*1.03
pind <- round((x-a)/b*100+1); print(min(pind)); print(max(pind))
cols <- c(rainbow(30, start=0.82, end=0.99), rainbow(70, start=0.01, end=0.17))[100:1]
map('state', fill = 1, interior=F, col = gray(0.85), ylim=c(22.5, 35), xlim=c(-88,-76))
#mtext(side = 1, line = 2.5, "longitude")
#mtext(side = 2, line = 2.5, "latitude")
points(sc2$V2, sc2$V3, col=cols[pind], pch=15, cex=0.5)
#box(); axis(1); axis(2, las = 2)
xloc <- seq(-86, -78, length.out=100)
for (j in 1:100) {   polygon(c(xloc[j], xloc[j+1],xloc[j+1], xloc[j]), c(23.0,23.0,23.4,23.4), col=cols[j], border=NA) }
w <- which.min(abs(((max(x)-min(x))/6) - pos))
if(-pos[w]<min(x)) { xx <- seq(0, max(x), pos[w]); xx <- xx[xx>min(x)] } else {  xx <- c(seq(-pos[w], min(x), -pos[w]), seq(0, max(x), pos[w])) }
text(xloc[round((xx-a)/b*100+1)], y=22.75, xx, pos=2)
text(-82, 24.25, "relative index of spawning output", pos = 1)

degs = seq(88, 76, -2)
a = sapply(degs, function(x) bquote(.(x)*degree ~ W))
axis(1, at = -degs, lab=do.call(expression, a))
degs = seq(24, 34, 2)
a = sapply(degs, function(x) bquote(.(x)*degree ~ N))
axis(2, at = degs, lab=do.call(expression, a), las = 2)
box()

#mtext(side=3, line=1.5, "scaled fecundity map", cex=1.3, font=2)


sc3 <- unscaled[which(unscaled$V6==2013 & unscaled$V7==6 & unscaled$V8==27),]
x <- sc3$V5 / 100
pos <- c(0.005,  0.05, 0.1, 0.2, 0.5, 1, 2, 4, 5, 10, 20, 50, 100, 1000)
a <- floor(min(x))
b <- max((x)-a)*1.03
pind <- round((x-a)/b*100+1); print(min(pind)); print(max(pind))
cols <- c(rainbow(30, start=0.82, end=0.99), rainbow(70, start=0.01, end=0.17))[100:1]
map('state', fill = 1, interior=F, col = gray(0.95), ylim=c(23.0, 35), xlim=c(-88,-76))
points(sc3$V2, sc3$V3, col=cols[pind], pch=15, cex=0.8)
box(); axis(1); axis(2)
xloc <- seq(-86, -78, length.out=100)
for (j in 1:100) {   polygon(c(xloc[j], xloc[j+1],xloc[j+1], xloc[j]), c(23.5,23.5,24.0,24.0), col=cols[j], border=NA) }
w <- which.min(abs(((max(x)-min(x))/6) - pos))
if(-pos[w]<min(x)) { xx <- seq(0, max(x), pos[w]); xx <- xx[xx>min(x)] } else {  xx <- c(seq(-pos[w], min(x), -pos[w]), seq(0, max(x), pos[w])) }
text(xloc[round((xx-a)/b*100+1)], y=23.2, xx, pos=2)
mtext(side=3, line=1.5, "unscaled fecundity map", cex=1.3, font=2)


# reduce number of particles ----------------------------------
min(scaled$V5); max(scaled$V5); mean(scaled$V5); sum(scaled$V5)
dim(scaled)

scaled$V5 <- round(scaled$V5/1000)  # 200 for large, 1000 for small

min(scaled$V5); max(scaled$V5); mean(scaled$V5); sum(scaled$V5)
table(scaled$V5==0)
scaled <- scaled[(scaled$V5>0),]

dim(scaled)
min(scaled$V5); max(scaled$V5); mean(scaled$V5); sum(scaled$V5)

# check particle reduction -------------------------------------
sc2 <- scaled[which(scaled$V6==2013 & scaled$V7==7 & scaled$V8==27),]
x <- sc2$V5
pos <- c(0.005,  0.05, 0.1, 0.2, 0.5, 1, 2, 4, 5, 20, 50, 100, 200, 500, 1000)
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

#scaledfin <- scaled
dim(scaled)

setwd("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes")

d <- read.table("bad.csv", sep = ",", header = F)
d$V4 <- d$V3
d$V3 <- paste0(round(d$V1, 4), "_", round(d$V2, 4))
scaled$V10 <- paste0(round(scaled$V2, 4), "_", round(scaled$V3, 4))
prob <- which(scaled$V10 %in% d$V3)
prob2 <- which(scaled$V10 == d$V3[3])

scaledfin <- c()

for (y in 2019)  {
  scaled2 <- scaled
  scaled2$V6 <- y
  scaledfin <- rbind(scaledfin, scaled2)   }

dim(scaledfin)/5
dim(scaled)
sum(scaledfin$V5)
head(scaledfin)

scaledfin <- scaledfin[-10]
head(scaledfin)

table(scaledfin$V1)
table(scaledfin$V2)
table(scaledfin$V3)
table(scaledfin$V4)
table(scaledfin$V5)
table(scaledfin$V6)
table(scaledfin$V7)
table(scaledfin$V8)
table(scaledfin$V9)
table(scaledfin$V10)

setwd("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/")

nests <- c("nest_1_SABGOM.nc", "nest_2_SABGOM.nc", 
           "nest_1_AtlMercator.nc",  
           "nest_1_20080501000000_HYCOM150.nc", 
           "nest_1_hycomGOM90pt1.nc", "nest_1_hycomGOM32pt5.nc", "nest_2_GLBHYCOMsmall.nc")

d2 <- findMinDepth(scaledfin$V2, scaledfin$V3, nests)
plot(-d2, -scaledfin$V4)
abline(1,1)
table(scaledfin$V4 < d2)

#plot(-d2[1:100], type = "l")
#lines(-scaledfin$V4[1:100], col = 2)

scaledfin$V4 <- round(d2 - 10)
table(scaledfin$V4)
scaledfin$V4[which(scaledfin$V4 > 50)] <- 50
table(scaledfin$V4)

table(scaledfin$V4[prob])
scaledfin$V4[prob] <- scaledfin$V4[prob] - 10
table(scaledfin$V4[prob])

table(scaledfin$V4 < d2)

table(scaledfin$V4[prob2])
scaledfin$V4[prob2] <- scaledfin$V4[prob2] - 10
table(scaledfin$V4[prob2])

scaledfin$V4[which(scaledfin$V4 <= 0)] <- 1
table(scaledfin$V4)

plot(-d2, -scaledfin$V4, pch = 19, cex = 2, col = "#FF000002")
abline(1,1)

# final check ----------------------------------

tapply(scaledfin$V5, scaledfin$V1 > 77, sum)
ga <- tapply(scaledfin$V5, scaledfin$V1 > 77, sum)
ga[1] / ga[2]   # 2.079171 for small; 2.10428 for large
arearat
sum(scaledfin$V5)   # 43404 per year for small; 224604 for large

table(scaledfin$V6) # 22938 per year for small; 34058 for large 
head(scaledfin)

# save output -------------------------------------------

write.table(scaledfin, file="scaledGOMATLrel2019.txt", sep="\t", col.names=F, row.names=F)

# the end -----------------------


# plot release file by date --------------------------------------
datelist <- unique(paste(scaled$V7, "_", scaled$V8, sep=""))
datelist1 <- sort((unique(paste(scaled$V7, ".", sprintf("%02d", scaled$V8), sep=""))))
datelist <- datelist[order(as.numeric(datelist1))]

dev.off()
par(mfrow=c(4, 6), mex=0.5)

x <- scaled$V5
pos <- c(0.005,  0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 40, 100, 500, 1000)
a <- floor(min(x))
b <- max((x)-a)*1.03
cols <- c(rainbow(30, start=0.82, end=0.99), rainbow(70, start=0.01, end=0.17))[100:1]

for (i in 1:length(datelist))  {
sc2 <- scaled[which(scaled$V6==2017 & scaled$V7==strsplit(datelist[i], "_")[[1]][1] & scaled$V8==strsplit(datelist[i], "_")[[1]][2]),]
x1 <- sc2$V5 
pind <- round((x1-a)/b*100+1); print(min(pind)); print(max(pind))
map('state', fill = 1, interior=F, col = gray(0.95), ylim=c(23.0, 35), xlim=c(-89,-77))
points(sc2$V2, sc2$V3, col=cols[pind], pch=15, cex=0.8)
box(); axis(1); axis(2)
xloc <- seq(-86, -78, length.out=100)
for (j in 1:100) {   polygon(c(xloc[j], xloc[j+1],xloc[j+1], xloc[j]), c(23.5,23.5,24.0,24.0), col=cols[j], border=NA) }
w <- 11 # w <- which.min(abs(((max(x)-min(x))/6) - pos))
if(-pos[w]<min(x)) { xx <- seq(0, max(x), pos[w]); xx <- xx[xx>min(x)] } else {  xx <- c(seq(-pos[w], min(x), -pos[w]), seq(0, max(x), pos[w])) }
#text(xloc[round((xx-a)/b*100+1)], y=23.2, xx, pos=2)
mtext(side=3, line=1.5, paste(month.abb[as.numeric(strsplit(datelist[i], "_")[[1]][1])], strsplit(datelist[i], "_")[[1]][2]), cex=1, font=2)
}

################################  END  #######################################
