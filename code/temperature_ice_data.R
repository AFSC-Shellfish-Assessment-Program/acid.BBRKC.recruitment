#ERA 5 SST and ice cover processing

#PURPOSE: To load and process sst and ice cover data from ERA 5 monthly averaged data
  # DATA DOWNLOAD/QUERY STEPS:
  # 1) Navigate here (will need to login): https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels-monthly-means?tab=overview
  # 2) Click on "Download data" tab
  # 3) Click on the product type you’d like. We have been using “Monthly averaged reanalysis”
  # 4) The ERA5 reanalysis has a many different data types. For SST, click on the “Temperature and pressure” 
  #    drop down, and then check the “Sea surface temperature” box. For ice cover, click on the “Other” drop down, 
  #    and then check the “Sea-ice cover” box. 
  # 5) Select years you’d like the data to cover. Sometimes including the most recent year with earlier years 
  #    in your selection results in an error, which can be resolved by downloading recent year data separately.
  # 6) Select months you’d like the data to cover. 
  # 7) Select geographical area for the data. Subsetting to a specific region speeds up the data querying/download.
  #    For the data below, the region has been North 60°, West -180°, South 52°, and East -155°. 
  # 8) Select NetCDF(experimental) as the data format, and submit form to query/download data. 

### LOAD PACKAGES -------------------------------------------------------------------------------------------------------

library(tidyverse)
library(ncdf4)
library(zoo)
library(maps)
library(mapdata)
library(chron)
library(fields)
library(oce)

### PROCESS ICE DATA ------------------------------------------------------------------------------------------------------
  # 1970-1974 data
    ice_nc3 <- nc_open("./data/ERA5_sstice_1970-1974.nc")
    
    ice3 <- ncvar_get(ice_nc3, "siconc", verbose = F) # select ice variable
    ice3 <- ifelse(ice3<0.0005, 0, ice3) # set threshold for very low ice cover
    dim(ice3) #101 lon, 33 lat, 60 months
    
    #Process
    h <- (ncvar_get(ice_nc3, "time")/24) # extract hours
    d3 <- dates(h, origin = c(1, 1, 1900)) # generate dates from hours
    m3 <- months(d3) # extract months from dates
    yr3 <- chron::years(d3) # extract years from dates
    
    x3 <- ncvar_get(ice_nc3, "longitude") # extract lon
    y3 <- ncvar_get(ice_nc3, "latitude") # extract lat
    
    # Keep track of corresponding latitudes and longitudes of each column:
    lat <- rep(y3, length(x3))
    lon <- rep(x3, each = length(y3))
    
    # Create data matrix
    ice3 <- aperm(ice3, 3:1)
    mat_ice3 <- t(matrix(ice3, nrow = dim(ice3)[1], ncol = prod(dim(ice3)[2:3])))
    
    data.frame(lon = lon, lat = lat,  mat_ice3) %>%
      pivot_longer(cols = c(3:62), names_to = "month", values_to = "ice") %>% # select ice variable colums
      mutate(month = rep(m3, 3333), year = rep(yr3, 3333), ice = ice) %>% # need to repeat months and years by the nrow of prev step divided by the # months (so 199980/60 in this case) 
      na.omit()-> ice_latlon_70.74
  
  # 1975-2021 data
    ice_nc4 <- nc_open("./data/ERA5_ice_1975-2021.nc")
    
    ice4 <- ncvar_get(ice_nc4, "siconc", verbose = F) # repeat steps from above
    ice4 <- ifelse(ice4<0.0005, 0, ice4)
    dim(ice4) #101 lon, 33 lat, 564 months
    
    #Process
    h <- (ncvar_get(ice_nc4, "time")/24)
    d4 <- dates(h, origin = c(1, 1, 1900))  
    m4 <- months(d4)
    yr4 <- chron::years(d4)
    
    x4 <- ncvar_get(ice_nc4, "longitude")
    y4 <- ncvar_get(ice_nc4, "latitude")
    
    # Keep track of corresponding latitudes and longitudes of each column:
    lat <- rep(y4, length(x4))
    lon <- rep(x4, each = length(y4))
    
    #Create ice data matrix
    ice4 <- aperm(ice4, 3:1)
    mat_ice4 <- t(matrix(ice4, nrow = dim(ice4)[1], ncol = prod(dim(ice4)[2:3])))
    
    data.frame(lon = lon, lat = lat,  mat_ice4) %>%
      pivot_longer(cols = c(3:566), names_to = "month", values_to = "ice") %>%
      mutate(month = rep(m4, 3333), year = rep(yr4, 3333), ice = ice) %>%
      na.omit()-> ice_latlon
    
  #2022 data
    nc4_22 <- nc_open("./data/ERA5_ice_sst_2022.nc")
    
    ice22 <- ncvar_get(nc4_22, "siconc", verbose = F)
    
    ice22 <- ifelse(ice22<0.0005, 0, ice22)
    dim(ice22) #101 lon, 33 lat, 12 months
    
    #Process
    h <- (ncvar_get(nc4_22, "time")/24)
    d4 <- dates(h, origin = c(1, 1, 1900))  
    m4 <- months(d4)
    yr4 <- chron::years(d4)
    
    x4 <- ncvar_get(nc4_22, "longitude")
    y4 <- ncvar_get(nc4_22, "latitude")
    
    # Keep track of corresponding latitudes and longitudes of each column:
    lat <- rep(y4, length(x4))
    lon <- rep(x4, each = length(y4))
  
  
    #Create ice data matrix
    ice22 <- aperm(ice22, 3:1)
    mat_ice22 <- t(matrix(ice22, nrow = dim(ice22)[1], ncol = prod(dim(ice22)[2:3])))
    
    data.frame(lon = lon, lat = lat,  mat_ice22) %>%
      pivot_longer(cols = c(3:14), names_to = "month", values_to = "ice") %>%
      mutate(month = rep(m4, 3333), year = rep(yr4, 3333), ice = ice) %>%
      na.omit()-> ice_latlon22
  
  #bind with earlier ice data, filter by Jan-Apr and subset to BB extent, calculate mean by month
  rbind(ice_latlon_70.74, ice_latlon, ice_latlon22) %>%
    filter(lon >= -168 & lon <= -158, lat >= 54 & lat <= 60, month %in% c("Jan", "Feb", "Mar", "Apr")) %>%
    group_by(year, month) %>%
    reframe(mean.ice = mean(ice)) -> ice_JanApr_1970.2022

  #write csv
  write.csv(ice_JanApr_1970.2022, "./output/ERA5ice_JanApr_1970.2022.csv")

### PROCESS SST DATA ------------------------------------------------------------------------------------------------------
  # 1970-1974 data
    sst_nc3 <- nc_open("./data/ERA5_sstice_1970-1974.nc")
    
    sst3 <- ncvar_get(sst_nc3, "sst", verbose = F) # keep as Kelvin
    dim(sst3) #101 lon, 33 lat, 60 months
    
    #Process
    h <- (ncvar_get(sst_nc3, "time")/24) # same steps as above
    d3 <- dates(h, origin = c(1, 1, 1900))  
    m3 <- months(d3)
    yr3 <- chron::years(d3)
    
    x3 <- ncvar_get(sst_nc3, "longitude")
    y3 <- ncvar_get(sst_nc3, "latitude")
    
    # Keep track of corresponding latitudes and longitudes of each column:
    lat <- rep(y3, length(x3))
    lon <- rep(x3, each = length(y3))
    
    #Create ice data matrix??
    sst3 <- aperm(sst3, 3:1)
    mat_sst3 <- t(matrix(sst3, nrow = dim(sst3)[1], ncol = prod(dim(sst3)[2:3])))
    
    data.frame(lon = lon, lat = lat,  mat_sst3) %>%
      pivot_longer(cols = c(3:62), names_to = "month", values_to = "sst") %>%
      mutate(month = rep(m3, 3333), year = rep(yr3, 3333), sst = sst) %>%
      na.omit()-> sst_latlon_70.74

  # 1974-2021 data
    sst_nc4 <- nc_open("./data/ERA5_sst_1975-2021.nc")
    
    sst4 <- ncvar_get(sst_nc4, "sst", verbose = F) # same steps as above
    dim(sst4) #101 lon, 33 lat, 564 months
    
    #Process
    h <- (ncvar_get(sst_nc4, "time")/24)
    d4 <- dates(h, origin = c(1, 1, 1900))  
    m4 <- months(d4)
    yr4 <- chron::years(d4)
    
    x4 <- ncvar_get(sst_nc4, "longitude")
    y4 <- ncvar_get(sst_nc4, "latitude")
    
    # Keep track of corresponding latitudes and longitudes of each column:
    lat <- rep(y4, length(x4))
    lon <- rep(x4, each = length(y4))
    
    #Create sst data matrix
    sst4 <- aperm(sst4, 3:1)
    mat_sst4 <- t(matrix(sst4, nrow = dim(sst4)[1], ncol = prod(dim(sst4)[2:3])))
    
    data.frame(lon = lon, lat = lat,  mat_sst4) %>%
      pivot_longer(cols = c(3:566), names_to = "month", values_to = "sst") %>%
      mutate(month = rep(m4, 3333), year = rep(yr4, 3333), sst = sst) %>%
      na.omit()-> sst_latlon

  #2022 data
    sst22 <- ncvar_get(nc4_22, "sst", verbose = F) 
    
    dim(sst22) #101 lon, 33 lat, 2 months
    
    #Process
    h <- (ncvar_get(nc4_22, "time")/24)
    d4 <- dates(h, origin = c(1, 1, 1900))  
    m4 <- months(d4)
    yr4 <- chron::years(d4)
    
    x4 <- ncvar_get(nc4_22, "longitude")
    y4 <- ncvar_get(nc4_22, "latitude")
    
    # Keep track of corresponding latitudes and longitudes of each column:
    lat <- rep(y4, length(x4))
    lon <- rep(x4, each = length(y4))
    
    #Create sst data matrix
    sst22 <- aperm(sst22, 3:1)
    mat_sst22 <- t(matrix(sst22, nrow = dim(sst22)[1], ncol = prod(dim(sst22)[2:3])))
    
    data.frame(lon = lon, lat = lat,  mat_sst22) %>%
      pivot_longer(cols = c(3:14), names_to = "month", values_to = "sst") %>%
      mutate(month = rep(m4, 3333), year = rep(yr4, 3333), sst = sst) %>%
      na.omit()-> sst_latlon22
    
  #bind with earlier sst data, subset to BB extent, calculate mean by month
    rbind(sst_latlon_70.74, sst_latlon, sst_latlon22) %>%
      filter(lon >= -168 & lon <= -158, lat >= 54 & lat <= 60) %>%
      group_by(year, month) %>%
      reframe(mean.sst = mean(sst)) -> sst_1970.2022_BristolBay
    
    write.csv(sst_1970.2022_BristolBay, "./data/sst_1970.2022_BristolBay")
  
  
  