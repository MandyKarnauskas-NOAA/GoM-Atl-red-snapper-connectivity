
rm(list=ls())

setwd("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/final_for_subsample")

foldersel <- dir()
numsel <- 100
finres <- c(NA, NA, NA, NA)
par(mfrow = c(5, 2), mar = c(3, 3, 0, 0))

for (b in 1:10) {
  mod <- rep(NA, numsel)
  num <- rep(NA, numsel)
  
for (i in 1:numsel) {
    
  x1 <- sample(1:length(foldersel), 1)
  mod[i] <- foldersel[x1]

  filelist <- list.files(path = foldersel[x1], pattern="con_file")     #   find files
  x2 <- sample(1:length(filelist), 1)

  dat <- read.table(paste0(foldersel[x1], "/", filelist[x2]))
  colnames(dat) <- c("rel_poly","ret_poly","ret_yr","ret_mo","ret_d","age","ret_dep","rel_yr","rel_mo","rel_d")

  dat$ret_reg <- "GOM"
  dat$ret_reg[which(dat$ret_poly >= 77)] <- "ATL"
  dat$rel_reg <- "GOM"
  dat$rel_reg[which(dat$rel_poly >= 77)] <- "ATL"
  
  gomdat <- dat[which(dat$rel_reg == "GOM"),]
  atldat <- dat[which(dat$rel_reg == "ATL"),]
    
  biomfact <- b
  
  samp <- sample(1:nrow(gomdat), nrow(gomdat) * biomfact, replace = TRUE)
  length(samp)/nrow(gomdat)   # check
  gomdat2 <- gomdat[samp,]
  dat2 <- rbind(atldat, gomdat2)
  
  tab <- table(dat2$rel_reg, dat2$ret_reg)
  num[i] <- tab[2, 1] / (tab[1, 1] + tab[2, 1])   # of all recruits in Atlantic, what % came from Gulf
  }
  
  boxplot(num ~ mod, xlab = "model", ylab = "percent of recruits in Atl\nspawned from Gulf")
  
  res <- c(biomfact, quantile(num, mean(num), probs = c(0.025, 0.50, 0.975)))
  finres <- rbind(finres, res)
}

#mean(num)
#max(num)
#min(num)
#hist(num)
#tapply(num, mod, mean)

finres

