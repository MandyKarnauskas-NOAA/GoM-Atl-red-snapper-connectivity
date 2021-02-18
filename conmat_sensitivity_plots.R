
rm(list = ls())
library(maps)
library(matlab)

models <- c("SABGOM", "HYCOM", "Mercator")
spps   <- c("lane", "mutton", "grey")
s1 <- c("OVM 1", "OVM 2", "OVM 3")

for (m in 1: 3)  { 

  model <- models[m]  

  png(filename = paste0("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/plots/conmatALL-", model, ".png"), units="in", width=9.9, height=6.7, pointsize=12, res=72*20)
  
  nf <- layout(matrix(c(1:6, 0, 1, 7:11, 18, 1, 12:16, 0, 0, rep(17, 5), 0), 4, 7, byrow=TRUE), c(2, rep(6,5), 3), c(6,6,6,2))
  layout.show(nf)

  ####  global y label
  par(mar=c(0,0,0,3))
  plot.new()
  mtext(side=2, line=-3, "Source Node", font=2)
    
  for (s in 1:3)  {
    
  spp <- spps[s]
    
  fold <- paste0("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/final_for_subsample/", model, "_", spp)
  setwd(fold) 
    
  filelist <- list.files(path = ".", pattern="con_file")     #   find files
  dat <- read.table(filelist[1])
  filelist <- filelist[-1]
  for (i in filelist)  { newdat <- read.table(i); dat <- rbind(dat, newdat) }
  colnames(dat) <- c("rel_poly","ret_poly","ret_yr","ret_mo","ret_d","age","ret_dep","rel_yr","rel_mo","rel_d")

  d <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/MASTER_codes/CMS_input_files/redSnapperSett_GOM_ATL_hires.xyz", sep="\t")
  names(d) <- c("lon", "lat", "pol")
  pols <- unique(d$pol)
  tklocs <- c(89.5, 76.5, 95.5, 117.5, 106.5, 41.5) 
  lablocs <- c(21.0, 59.0, 83.0, 92.5, 101.0, 112.0) 

  con <- dat
  colnames(con) <- c("source","sink","yr","mon","day","tim","ret_dep","rel_yr","rel_mo","rel_d")
  pollis <- 1:max(d$pol)    # list of polygon names
  limit <- length(pollis)

for (k in unique(con$yr))  {
  con1 <- con[which(con$yr == k),]
  
  par(mar=c(3, 3, 3, 1), xpd=F)
  settle <- matrix(NA, length(pollis), length(pollis))
  for (i in 1:length(pollis))     {
    a <- which(con1$source==pollis[i])
    if(length(a)>0) {
      for (j in 1:length(pollis)) {   settle[i,j] <- length(which(con1$sink[a]==j))    }
    } else { settle[i,] <- 0  }  }
  settle <- t(settle)
  nn <- 300
  cl = jet.colors(nn)
  image(pollis[42:limit], pollis[42:limit], log(settle[42:limit,42:limit]+1), col=cl, axes=F,  xlab="", ylab="")
  mtext(side = 3, line = 0.1, paste(s1[s], "-", k), cex = 0.8, font = 2)
  
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

  }
  
###  global x label
par(mar=c(0,0,3,0))
plot.new()
mtext(side=1, line=-3, "Receiving Node", font=2)

####  color bar plot
par(mar=c(1,0,1,0))   
plot(c(1.05,1.18), c(-0.5,1.5),  axes=F, xlab="", ylab="", col="white")
bra <- (seq(-0.5, ceiling(log(max(settle))), (ceiling(log(max(settle)))+0.5)/nn))
i = seq(0,1,1/(nn-1))
#rect(1.05, i, 1.1, 1, col = cl, lwd=0, border=cl)
rect(1.11, i, 1.12, 1.3, col = cl, lwd=0, border=cl)
laa <- c(1,10,100,1000)   # laa <- seq(0, 60000, 10000)
levs <- (0+i[1:nn]+1/(nn*2))
re <- lm(bra[1:(length(bra)-1)]~levs)
poss <- (log(laa) - re$coef[1])/re$coef[2]
text(1.115, poss, laa, pos=4, cex=1.2)
text(1.115, -0.65, cex=1.1, "number of\nsuccessful\nrecruits", pos=3)
#legend("bottom", "self-recruitment", lwd=1, lty=2, bty="n", cex=1.1)


dev.off()

}
