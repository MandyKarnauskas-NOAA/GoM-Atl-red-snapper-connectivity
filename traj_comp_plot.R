
rm(list = ls())
library(ncdf4)
library(splancs)
library(maps)
library(RColorBrewer)

models <- c("SABGOM", "HYCOM", "Mercator")
spps   <- c("lane", "mutton", "grey")
s1 <- c("OVM 1", "OVM 2", "OVM 3")

png(filename = "C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/plots/comp_traj_newbath1.png", units="in", width=5.9, height=7.5, pointsize=12, res=72*20)

par(mfrow = c(3, 3), mar = c(2, 4, 0, 0), xpd = F, mex = 0.9)

for (m in 1: 3)  { 
  for (s in 1:3)  {
    
    model <- models[m]
    spp <- spps[s]
    
    fold <- paste0("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/final_for_subsample/", model, "_", spp)
    setwd(fold) 
    
    lon <- c()
    lat <- c()
    dep <- c()
    cod <- c()
    
    tlis <- dir()[grep("traj", dir())]
    for (f in 1:length(tlis)) {
      nam <- tlis[1]
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
    
    load("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/GEBCO_bathy_hires.RData")        # stored GEBCO data and color scheme
    cols <- colorRampPalette(brewer.pal(9, 'Blues'))(150)[110:11]
    cols[100] <- "white"
    
    p <- cbind(c(-83, -70, -70, -81.5, -81.5, -83), c(23, 23, 40, 40, 24.58, 24.58))
    
    map('usa', xlim=c(-86.25, -75.25), ylim=c(23, 35.75), col=1)
#   map('usa', xlim=c(-85, -80), ylim=c(23, 26), col=0)
    ss <- 1
    image(x[seq(1, 3360, ss)], y[seq(1, 3360, ss)], z[seq(1, 3360, ss), seq(1, 3360, ss)], 
      xlim=c(-86.25, -75.25), ylim=c(23, 35.75), col=cols, axes=T, xlab="", ylab="", add=T)
    
    degs = seq(86, 78, -4)
    a = sapply(degs, function(x) bquote(.(x)*degree ~ W))
    axis(1, at = -degs, lab=do.call(expression, a))
    degs = seq(24, 32, 4)
    a = sapply(degs, function(x) bquote(.(x)*degree ~ N))
    axis(2, at = degs, lab=do.call(expression, a))
    box()
      #    axis(1); axis(2, las = 2); box()
    
    k <- which(cod == (-4))
    ga <- rep(NA, length(k))
    lon2 <- lon[, k]
    lat2 <- lat[, k]
    last <- function(x) { x[!is.na(x)][length(x[!is.na(x)])] }
    lo <- apply(lon2, 2, last) - 360
    la <- apply(lat2, 2, last)
    stg <- inout(cbind(lon2[1, ]-360, lat2[1, ]), p)
    eng <- inout(cbind(lo, la), p)
    
    ga[abs(stg - 1) + eng*1 == 2] <- 1
#    cat(c(model, "   "))
#    cat(table(stg, eng))
#    cat("\n")
    
   n <- which(ga == 1)    
    n1 <- n[seq(1, length(n), length.out = 1000)]
    for (j in n1) {
      lines(lon2[, j]-360, lat2[, j], col="#FFFF0010")  } 
    for (j in n1) {
      lo <- lon2[!is.na(lon2[,j]),j][length(lon2[!is.na(lon2[,j]),j])] - 360
      la <- lat2[!is.na(lat2[,j]),j][length(lat2[!is.na(lat2[,j]),j])]
        points(lon2[1,j]-360, lat2[1,j], col="#00b30020", pch=19, cex=0.4)
        points(lo, la, col="#FF000020", pch=19, cex=0.4)  }         
    # polygon(p)

    text(-86, 35, model, font = 2, cex = 1.2, pos = 4)
    text(-86, 34, s1[s], font = 2, cex = 1.2, pos = 4)
    }
}
table(stg, eng)

dev.off()

