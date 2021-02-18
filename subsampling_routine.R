
rm(list=ls())

setwd("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/final_for_subsample")

foldersel <- dir()
write.table(c(), file = "C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/outputs2.txt", 
            quote = FALSE, append = F, col.names = F)


for (i in 1:length(foldersel)) {

    filelist <- list.files(path = foldersel[i], pattern="con_file")     #   find files
    
    for (j in 1:length(filelist))  {
      
    dat <- read.table(paste0(foldersel[i], "/", filelist[j]))
  
    colnames(dat) <- c("rel_poly","ret_poly","ret_yr","ret_mo","ret_d","age","ret_dep","rel_yr","rel_mo","rel_d")
    
    dat$ret_reg <- "GOM"
    dat$ret_reg[which(dat$ret_poly >= 77)] <- "ATL"
    dat$rel_reg <- "GOM"
    dat$rel_reg[which(dat$rel_poly >= 77)] <- "ATL"

    tab <- table(dat$rel_reg, dat$ret_reg)
    num <- tab[2, 1] / (tab[1, 1] + tab[2, 1])   # of all recruits in Atlantic, what % came from Gulf
    
    res <- c(foldersel[i], j, num, as.vector(tab))
    
    write.table(matrix(res, 1, 7), file = "C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/outputs2.txt", 
                quote = FALSE, append = T, col.names = F)
    }
}  
  
rm(list = ls())

d <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/outputs2.txt", sep = " ")
d <- d[,2:4]
names(d) <- c("fold", "year", "ga")
d$model <- unlist(strsplit(as.character(d$fold), "_"))[seq(1, nrow(d)*2, 2)]
d$spp   <- unlist(strsplit(as.character(d$fold), "_"))[seq(2, nrow(d)*2, 2)]

head(d)
table(as.numeric(as.factor(d$model))+1, d$model)
table(as.numeric(as.factor(d$spp)), d$spp)

plot(jitter(d$year), d$ga, col = as.numeric(as.factor(d$model))+1, pch = as.numeric(as.factor(d$spp)), cex = 2,
     xlab = "", ylab = "percentage of Atlantic recruits spawned from Gulf", axes = F)
axis(2, las = 2); box()
legend("topright", c(unique(d$model), unique(d$spp)), pch = c(19, 19, 19, 1:3), col = c(2:4, 1, 1, 1))


png(filename="summary_boxplot.png", units="in", width=3.5, height=3.5, pointsize=12, res=72*10)

ggplot(d, aes(y = ga*100, fill = spp), axes = F) +
  geom_boxplot(width = 3) + facet_wrap(~model) +
  scale_y_continuous(breaks = pretty(c(0, 35), n = 5)) +
  scale_x_discrete(name="", breaks = 1, labels = "") + 
  ylab("% of Atlantic recruits spawned in Gulf") +
  scale_fill_discrete(name = "  vertical\nmigration\n behavior") +
  theme(legend.position = c(0.85, 0.7))

dev.off()

tapply(d$ga, d$model, mean)
tapply(d$ga, d$model, sd)
tapply(d$ga, d$model, sd) / tapply(d$ga, d$model, mean)

d$model <- as.factor(d$model)
d$spp <- as.factor(d$spp)
d$year <- as.factor(d$year)

out <- lm(d$ga ~ d$model + d$spp + d$year)
summary(out)
anova(out)


rm(list = ls())
d <- read.table("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/outputs2.txt", sep = " ")
d <- d[,2:8]
names(d) <- c("fold", "year", "gaper", "aa", "ga", "ag", "gg")
d$model <- unlist(strsplit(as.character(d$fold), "_"))[seq(1, nrow(d)*2, 2)]
d$spp   <- unlist(strsplit(as.character(d$fold), "_"))[seq(2, nrow(d)*2, 2)]

d$GoM2Atl <- d$ga / (d$aa + d$ga)
round(d$gaper, 6) - round(d$GoM2Atl, 6)

d2 <- cbind(d$aa, d$ga, d$ag, d$gg)
barplot(t(d2), beside = F)
rowSums(d2)

d$model <- as.factor(d$model)

tapply(rowSums(d2), d$model, mean)
boxplot(rowSums(d2) ~ d$model)
anova(lm(rowSums(d2) ~ d$model))  # SABGOM and Mercator similar total survival rates; HYCOM slightly lower

d3 <- d2 / rowSums(d2)
barplot(t(d3), beside = F)

boxplot(d$aa ~ d$model)     # SABGOM much stronger self-recruitment in Atl
anova(lm(d$aa ~ d$model))

boxplot(d$ga ~ d$model)     # SABGOM slightly lower transport Gulf to Atl
anova(lm(d$ga ~ d$model))

boxplot(d$ag ~ d$model)
anova(lm(d$ag ~ d$model))

boxplot(d$gg ~ d$model)
anova(lm(d$gg ~ d$model))



# subsampling routine ----------------------

rm(list=ls())

setwd("C:/Users/mandy.karnauskas/Desktop/RS_FATEproject/FINALRUNS_DEC2020/final_for_subsample")

foldersel <- dir()
numsel <- 1000
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













# old code 
##########  plot new grid cells  #########################

pols <- unique(d$pol)
are <- rep(NA, length(unique(d$pol)))

plot(0, ylim=c(24, 36), xlim=c(-100, -75))  # check
for (i in 1:117)  {
  f <- d[which(d$pol==i),]
  polygon(f, border=1)  }
map('state', add=T)

text(tapply(d$lon, d$pol, mean), tapply(d$lat, d$pol, mean), unique(d$pol), cex=0.5)

are[1:41] <- "W GoM"
are[42:76] <- "E GoM"
are[77:89] <- "E FL"
are[90:95] <- "GA"
are[96:106] <- "SC"
are[107:117] <- "NC"

are <- as.factor(are)
tab <- table(pols, are)

map('state', fill = 1, interior=F, col = gray(0.95), ylim=c(23,36), xlim=c(-99,-75))
for (j in unique(d[,3]))  {
  m <- d[which(d[,3]==j),]; polygon(m, lwd=3, col=as.numeric(are)[j]+1, border = 1, lwd=0.5)  }
axis(1); axis(2); box()
for (j in unique(are)) {  text( mean(d$lon[d$pol %in% which(are==j)]),  mean(d$lat[d$pol %in% which(are==j)]), j, col=1) }
abline(v=-89.1, lty=2); text(-87.5, 24, "longitude = -89.1")

mids <- tapply(pols, are, mean)
maxs <- tapply(pols, are, max)
mins <- tapply(pols, are, min)

tklocs <- maxs + 0.5
lablocs <- sort((tklocs-mins-0.5)/2 + mins)       # used in plotting later

tklocs <- c(89.5, 76.5, 95.5, 117.5, 106.5, 41.5) 
lablocs <- c(21.0, 59.0, 83.0, 92.5, 101.0, 112.0) 

efl <- which(tab[,1]==1)
egom <- which(tab[,2]==1)
ga <- which(tab[,3]==1)
nc <- which(tab[,4]==1)
sc <- which(tab[,5]==1)
wgom <- which(tab[,6]==1)

dat$state <- NA
dat$state[which(dat$ret_poly<=max(egom))] <- "E GOM"
dat$state[which(dat$ret_poly>=min(wgom) & dat$ret_poly<=max(wgom))] <- "W GOM"
dat$state[which(dat$ret_poly>=min(efl) & dat$ret_poly<=max(efl))] <- "E FL"
dat$state[which(dat$ret_poly>=min(ga) & dat$ret_poly<=max(ga))] <- "GA"
dat$state[which(dat$ret_poly>=min(sc) & dat$ret_poly<=max(sc))] <- "SC"
dat$state[which(dat$ret_poly>=min(nc))] <- "NC"


##################################################################################
