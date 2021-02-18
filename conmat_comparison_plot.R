
rm(list = ls())
library(viridis)

models <- c("SABGOM", "HYCOM", "Mercator")
spps   <- c("lane", "mutton", "grey")
s1 <- c("OVM 1", "OVM 2", "OVM 3")

png(filename = "C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/plots/comp_conmat.png", units="in", width=8.6, height=8, pointsize=12, res=72*20)

nf <- layout(matrix(c(11, 1:3, 11, 4:6, 11, 7:9, 0, rep(10, 3)), 4, 4, byrow=TRUE), c(2, rep(10, 3)), c(rep(10, 3), 2))
layout.show(nf)

for (m in 1: 3)  { 
  for (s in 1:3)  {
 
  model <- models[m]
  spp <- spps[s]

  fold <- paste0("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/final_for_subsample/", model, "_", spp)
  setwd(fold) 
   
  filelist <- list.files(path = ".", pattern="con_file")     #   find files
  dat <- read.table(filelist[1])
  filelist <- filelist[-1]
    for (f in filelist)  { newdat <- read.table(f); dat <- rbind(dat, newdat) }
  colnames(dat) <- c("rel_poly","ret_poly","ret_yr","ret_mo","ret_d","age","ret_dep","rel_yr","rel_mo","rel_d")

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

  pols <- 1:10

# plot for all years together ----------------

colnames(dat) <- c("rel", "ret", "yr","mon","day","tim","ret_dep","rel_yr","rel_mo","rel_d", "source", "sink")
pollis <- 1:10    # list of polygon names
limit <- length(pollis)

par(mar=c(6, 6, 0.5, 0.5), xpd=F)

con1 <- dat                     # to look at all years together
settle <- matrix(NA, length(pollis), length(pollis))
for (i in 1:length(pollis))     {
  a <- which(con1$source==pollis[i])
  if(length(a)>0) {
    for (j in 1:length(pollis)) {   settle[i,j] <- length(which(con1$sink[a]==j))    }
  } else { settle[i,] <- 0  }  }

settle <- t(settle)
nn <- 300
cl = viridis(nn, end = 0.85)[nn:1]
image(pollis, pollis[1:9], (settle[, 1:9]), col=cl, breaks = c(0, seq(1, 20000, length.out = nn-1), 40000), 
      axes=F,  xlab="", ylab="")

text(0.5, 8.75, model, font = 2, cex = 1.2, pos = 4)
text(0.5, 7.75, s1[s], font = 2, cex = 1.2, pos = 4)

cat(max(settle), na.rm = T)
cat("\n")

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
      text(i, j, txt[i,j], col = 0, cex = 0.8) }}}
  }
}

par(mar = c(0, 0, 0, 0))
plot.new()  
mtext(side=1, line=-2, "               Receiving Node", font=2)  

par(mar = c(0, 0, 0, 0))
plot.new()
mtext(side=2, line= (-2), "               Source Node", font=2)

dev.off()

