
rm(list = ls())

#library(sf)
library(maps)
library(shapefiles)
library(shape)
library(yarrr)

setwd("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes")

p <- read.shapefile("C://Users/mandy.karnauskas/Desktop/KevinCraig_shrimpPaper/SA_EEZ/comm_permits_SA_SG") 
p <- as.matrix(cbind(p$shp$shp[[1]]$points))
p2 <- read.shapefile("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/EEZ/MSA_FMC_GOM_FedWaters")
p2 <- as.matrix(cbind(p2$shp$shp[[1]]$points))

d <- read.table("CMS_input_files/redSnapperSett_GOM_ATL_hires.xyz", sep="\t")

map("state", ylim=c(22, 28), xlim=c(-85, -79))  # check
for (i in 1:117)  {
  f <- d[which(d$V3==i),]
  polygon(f, col=i)  
  text(mean(f$V1), mean(f$V2), f$V3)}
polygon(p, border=1, lwd=2, lty = 1)
polygon(p2, border=8, lwd=2, lty = 2)
axis(1); axis(2)
abline(h=27)
abline(v = (-82))

load("GEBCO_GoM_ATL.RData")
z[which(z > 0)] <- 0

cols <- rainbow(100, start = 0.5, end = 0.6)

labcol <- rep("Western Gulf", 117)
labcol[46:52] <- "Panhandle"
labcol[53:60] <- "Big Bend"
labcol[61:66] <- "Central W FL"
labcol[67:71] <- "SW FL"
labcol[72:75] <- "FL Keys Shelf"
labcol[76] <- " S FL"
labcol[77:79] <- "S FL"
labcol[80:89] <- "E FL"
labcol[90:95] <- "GA"
labcol[96:106] <- "SC"
labcol[107:117] <- "NC"


png(filename = "C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/plots/mapfigure.png", 
    units="in", width=6*1.5, height=5*1.2, pointsize=12, res=72*20)

map("world", xlim = c(-98, -74.2), ylim = c(23.6, 36), col = 0)
mtext(side = 1, line = 2, "Longitude")
mtext(side = 2, line = 3, "Latitude")

ss <- 1
image(x[seq(1, 6000, ss)], y[seq(1, 3240, ss)], 
      log(abs(z[seq(1, 6000, ss), seq(1, 3240, ss)])), col = cols, add = T)
axis(1, at = seq(-95, -75, 5), lab = paste0(seq(-95, -75, 5), "W"))
axis(2, las = 2, at = seq(25, 35, 5), lab = paste0(seq(25, 35, 5), "N"))
box()

col2 <- rep(c(gray(level = 0.55, alpha = 0.6), gray(level = 0.7, alpha = 0.6)), 11)
col2[8] <- transparent(6, trans.val = 0.5)

m <- 1
for (i in 2:117)  {
  f <- d[which(d$V3==i),]
  if (labcol[i] != labcol[i - 1])  { m <- m + 1}
  polygon(f, col = col2[m], border = NA)  }

map("state", add = T, boundary = F, interior = T, col = gray(0.7))

sss <- 0.7

text(-93.4, 28.7, labcol[1], cex = sss)
text(-86.6, 29.7, labcol[46], cex = sss)
text(-84.4, 29, labcol[53], cex = sss)
text(-84.1, 27.5, labcol[61], cex = sss)
text(-83, 26.2, labcol[67], cex = sss)
text(-82.7, 25.12, labcol[72], cex = sss)
text(-82.2, 24.2, labcol[76], cex = sss)
#lines(c(-83.75, -82.75), c(24.75, 24.44074))
#text(-81.0, 24.5, labcol[77], cex = sss)
text(-80.5, 29.6, labcol[80], cex = sss)
text(-80.4, 31.1, labcol[90], cex = sss)
text(-79.2, 32.6, labcol[96], cex = sss)
text(-77, 34, labcol[107], cex = sss)


polygon(p, border = 1, lty = 2, lwd = 1)
polygon(p2, border = 1, lty = 2, lwd = 1)
text(-90, 27.1, "U.S. Gulf of Mexico")
text(-78.4, 31.1, "South\nAtlantic" )

for (i in 77:79)  {
  f <- d[which(d$V3==i),]
  polygon(f, col = col2[8], border = col2[8], lwd = 3)  }


Arrows(-97, 34.8, x1 = -97, y1 = 35.4, arr.length = 0.6, arr.width = 0.4, arr.type = "triangle")
text(-97, 34.4, "N", font = 2, cex = 1.5)

dev.off()


#states <- st_as_sf(map("state", plot = FALSE, fill = TRUE))
#states <- cbind(states, st_coordinates(st_centroid(states)))
#states$ID <- toTitleCase(as.character(states$ID))
#text(states$X, states$Y, (states$ID), cex = 0.8)

