
rm(list=ls())

gr <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/CMS_inputfiles/vert_matrix_gray_nobuoy", 
                 sep = "\ ", header = F, skip = 4)
la <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/CMS_inputfiles/vert_matrix_lane_nobuoy", 
                 sep = "\ ", header = F, skip = 4)
mu <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/CMS_inputfiles/vert_matrix_red_sd0_nobuoy", 
                 sep = "\ ", header = F, skip = 4)

gr1 <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/CMS_inputfiles/vert_matrix_gray_nobuoy", 
                 sep = "\t", header = F, skip = 2)
la1 <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/CMS_inputfiles/vert_matrix_lane_nobuoy", 
                 sep = "\t", header = F, skip = 2)
mu1 <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/CMS_inputfiles/vert_matrix_red_sd0_nobuoy", 
                 sep = "\t", header = F, skip = 2)

gr1[1,] == la1[1,]
gr1[1,] == mu1[1,]
gr1[2,] == la1[2,]
gr1[2,] == mu1[2,]

mu <- mu/100
gr <- gr/100
la <- la/100

tim <- as.numeric(unlist(strsplit( as.character(gr1[2,]), " ")))/24/3600
times <- rep(0, 6)
for (i in 1:5) { times[i+1] <- sum(tim[1:i])}
times
#times <- c(0, 2, 13, 17, 26, 30)                                    # fill in by hand
brks <- as.numeric(unlist(strsplit( as.character(gr1[1,]), " ")))

##### day zero, # days to preflex, # days to flex, # days to postflex, min PLD, max PLD
tm <- length(diff(times))
labs <- NULL; for (i in 1:tm) { labs[i] <- paste(times[i], "-", times[i+1]) }
x <- seq(0.5, (tm + 1), length.out = length(brks))
y <- -brks

plot(x,y, col = 0, axes = F, xlab = "days after spawning", ylab = "depth (m)", ylim = c(-70, 0))
axis(1, at=c(1:tm), lab=labs)
axis(3, line=-0.7, at=c(1:tm), lab = c("hatching","preflexion","flexion","postflexion", "settlement"), tick=F)
box(); axis(2, las=2)
abline(v=1, col=8); abline(v=2, col=8); abline(v=3, col=8); abline(v=4, col=8); abline(v=5, col=8)

lines(1 + la$V1, y, lwd=3, col = 1, lty = 1)
lines(2 + la$V2, y, lwd=3, col = 1, lty = 1)
lines(3 + la$V3, y, lwd=3, col = 1, lty = 1)
lines(4 + la$V4, y, lwd=3, col = 1, lty = 1)
lines(5 + la$V5, y, lwd=3, col = 1, lty = 1)

lines(1 + mu$V1, y, lwd=3, col = "#00000050", lty = 3)
lines(2 + mu$V2, y, lwd=3, col = "#00000050", lty = 3)
lines(3 + mu$V3, y, lwd=3, col = "#00000050", lty = 3)
lines(4 + mu$V4, y, lwd=3, col = "#00000050", lty = 3)
lines(5 + mu$V5, y, lwd=3, col = "#00000050", lty = 3)

lines(1 + gr$V1, y, lwd=3, col = "#00000090", lty = 2)
lines(2 + gr$V2, y, lwd=3, col = "#00000090", lty = 2)
lines(3 + gr$V3, y, lwd=3, col = "#00000090", lty = 2)
lines(4 + gr$V4, y, lwd=3, col = "#00000090", lty = 2)
lines(5 + gr$V5, y, lwd=3, col = "#00000090", lty = 2)

legend("bottomright", c("OVM 1", "OVM 2", "OVM 3"), col = c(1, "#00000050", "#00000090"), lty = c(1, 3, 2), lwd = 3, bty = "n")


