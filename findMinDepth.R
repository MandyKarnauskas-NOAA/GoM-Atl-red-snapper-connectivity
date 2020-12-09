# M. Karnauskas 12/9/2020
# function for finding minimum depth at a list of locations for multiple oceanographic nest files
# input list of coordinates and character list of nest file names

findMinDepth <-  function(longitude, latitude, nestList)  {

if (!"ncdf4" %in% installed.packages()) install.packages("ncdf4", repos='http://cran.us.r-project.org')
library(ncdf4)  
  
for (i in 1:length(nestList))  { 
  
  nc1 <- nc_open(nestList[i])               # open netcdf and set temp variable

  u1 <- ncvar_get(nc1, nc1$var[[1]])        # get current direction
  lon1 <- nc1$var[[1]]$dim[[1]]$vals        # extract longitudes
  lat1 <- nc1$var[[1]]$dim[[2]]$vals        # extract latitudes
  dep1 <- nc1$var[[1]]$dim[[3]]$vals        # extract depths
  if (mean(lon1) > 0)  { lon1 <- lon1 - 360  }  # convert to -180 to 180 format if necessary
    nc_close(nc1)
  
  if (i == 1)  {                            # for first nest, initialize variables
    lons <- list(lon1)                      # list of longitudes
    lats <- list(lat1)                      # list of latitudes
    deps <- list(dep1)                      # list of depths
    us <- list(u1)                          # list of current velocities
    minlons <- list(min(lon1))              # 4 boundaries of nest file
    maxlons <- list(max(lon1))
    minlats <- list(min(lat1))
    maxlats <- list(max(lat1))
  } else {                                  # for subsequent nest files, append to list
    lons <- c(lons, lon1 = list(lon1)) 
    lats <- c(lats, lat1 = list(lat1))
    deps <- c(deps, dep1 = list(dep1))
    us <- c(us, u1 = list(u1))
    minlons <- append(minlons, min(lon1))
    maxlons <- append(maxlons, max(lon1))
    minlats <- append(minlats, min(lat1))
    maxlats <- append(maxlats, max(lat1))
  }
}                                           # end reading nest file info
  
#par(mfrow = c(3,1), mex = 0.2)             # for testing purposes only
#image(lons[[1]], lats[[1]], us[[1]][,,1])
#points(lon, lat, pch = 19, cex = 10)
#image(lons[[2]], lats[[2]], us[[2]][,,1])
#points(lon, lat, pch = 19, cex = 10)
#image(lons[[3]], lats[[3]], us[[3]][,,1])
#points(lon, lat, pch = 19, cex = 10)

if (length(latitude) == length(longitude)) {           # check that lat and lon are same length
  outdep <- rep(NA, length(latitude))            # placeholder for minimum depths
  
for (j in 1:length(latitude))  {                 # loop through coordinates
  depn <- rep(NA, length(nestList))         # placeholder for depths from each nest
  
  for (k in 1:length(nestList)) {           # loop through nest data
                                            # ensure that point falls within boundaries
    if (minlons[[k]] < longitude[j] & maxlons[[k]] > longitude[j] & minlats[[k]] < latitude[j] & maxlats[[k]] > latitude[j])  { 
                                            # if point is in boundary, find the maximum depth at which current data is not NA
    depn[k] <- deps[[k]][max(which(!is.na(us[[k]][which.min(abs(lons[[k]] - longitude[j])), which.min(abs(lats[[k]] - latitude[j])), ])))]
    }  
  }
  outdep[j] <- min(depn, na.rm = TRUE)    # minimum depth among nests is assigned as output
}
  }  else { print("error: lat and lon are not equal lengths") }

return(outdep)                              # return a list of depths that are the minimum across nests for each coordinate input

}
