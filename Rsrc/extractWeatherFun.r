# coordPlots = coordinates of plots; must have columns long and lat
# sY = start year, sM = start month, sD = start day, eY = end year, eM = end month, eD = end day
getWDpoints <- function(coordPlots, sY=1951, sM=1, sD=1, eY=2100, eM=12, eD=31) {
  # here we store the weather data
  wDs <- list()
  
  # this is the number of coordinate points
  nDFs <- nrow(coordPlots)
  
  for (i in 1:nDFs) {
    lon <- coordPlots$long[i]
    lat <- coordPlots$lat[i]
    
    # getting the weather data from clipick
    result <- getWeatherData(lon, lat, sY, sM, sD, eY, eM, eD)
    
    # weather datas as data.table
    wD <- data.table(t(sapply(result,c)))
    
    # fix the column names
    names(wD) <- as.character(wD[1,])
    wD <- wD[-1,]
    
    # this is the number of the day
    wD$rday <- 1:nrow(wD)
    
    wDs[[i]] <- wD
  }
  
  return(wDs)
}

