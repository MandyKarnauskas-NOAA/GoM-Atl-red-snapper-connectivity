
rm(list = ls())

models <- c("SABGOM", "HYCOM", "Mercator")
spps   <- c("lane", "mutton", "grey")
s1 <- c("OVM 1", "OVM 2", "OVM 3")

png(filename = "C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/plots/comp_traj.png", units="in", width=5.9, height=7.25, pointsize=12, res=72*20)

par(mfrow = c(3, 3), mar = c(1, 1, 0, 0), xpd = F, mex = 0.5)

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
    
    map('usa', xlim=c(-86.25, -75.25), ylim=c(23, 35.75), col=0)
   image(x[seq(1, 3360, 1)], y[seq(1, 3360, 1)], z[seq(1, 3360, 1), seq(1, 3360, 1)], 
      xlim=c(-86.25, -75.25), ylim=c(23, 35.75), col=cols, axes=T, xlab="", ylab="", add=T)
    axis(1); axis(2, las = 2); box()
    
    k <- which(cod == (-4))
    ga <- rep(NA, length(k))
    for (j in k) {
      lo <- lon[!is.na(lon[,j]),j][length(lon[!is.na(lon[,j]),j])] - 360
      la <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]
      if (lon[1,j]-360 < (-81.7) & lo > (-81.7))  { ga[j] <- 1  }  }
    
    n <- which(ga == 1)    
    n1 <- n[seq(1, length(n), length.out = 1000)]
    for (j in n1) {
      lines(lon[,j]-360, lat[,j], col="#FFFF0010")  } 
    for (j in n1) {
      lo <- lon[!is.na(lon[,j]),j][length(lon[!is.na(lon[,j]),j])] - 360
      la <- lat[!is.na(lat[,j]),j][length(lat[!is.na(lat[,j]),j])]
        points(lon[1,j]-360, lat[1,j], col="#00FF0020", pch=19, cex=0.4)
        points(lo, la, col="#FF000020", pch=19, cex=0.4)  }         

    text(-86, 35, model, font = 2, cex = 1.2, pos = 4)
    text(-86, 34, s1[s], font = 2, cex = 1.2, pos = 4)
    }
}
dev.off()

