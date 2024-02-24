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

## set up plots
theme_set(theme_bw())
cb <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

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
  
  ice <- ice_JanApr_1970.2022 %>%
    mutate(year = as.numeric(as.character(year)),
           month = as.factor(as.character(month)),
           mean.ice = as.numeric(mean.ice)) %>%
    pivot_wider(names_from = month,
                values_from = mean.ice) %>%
    filter(year %in% 1975:2022)

  # save raw data for plot
  plot_ice <- ice
  
  plot_ice$Jan_Feb_plot_ice  <- apply(plot_ice[,c(2,3)], 1, mean)
  plot_ice$Mar_Apr_plot_ice  <- apply(plot_ice[,c(4,5)], 1, mean)

  plot_ice <- plot_ice %>%
    select(-Jan, -Feb, -Mar, -Apr)  
  
  # scale each month
ice[,2:5] <- apply(ice[,2:5], 2, scale)
  
  ice$Jan_Feb_ice  <- apply(ice[,c(2,3)], 1, mean)
  ice$Mar_Apr_ice  <- apply(ice[,c(4,5)], 1, mean)
  
ice <- ice %>%
  select(-Jan, -Feb, -Mar, -Apr)

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
    
  #bind with earlier sst data, subset to BB extent, calculate mean by month, change K to C
    rbind(sst_latlon_70.74, sst_latlon, sst_latlon22) %>%
      filter(lon >= -168 & lon <= -158, lat >= 54 & lat <= 60) %>%
      group_by(year, month) %>%
      reframe(mean.sst = mean(sst)-273.15) -> sst_1970.2022_BristolBay
    
    write.csv(sst_1970.2022_BristolBay, "./data/sst_1970.2022_BristolBay")
    
# now process
    sst <- sst_1970.2022_BristolBay %>%
      mutate(year = as.numeric(as.character(year)),
             month = as.factor(as.character(month))) %>%
      filter(month %in% c("Jan", "Feb", "Mar", "Apr", "May", "Jun"),
             year %in% 1975:2022) %>%
      pivot_wider(names_from = month,
                  values_from = mean.sst)


    # save raw data for plot
    plot_sst <- sst 
    plot_sst$Jan_Jun_SST <- rowMeans(plot_sst[,2:7])
    
    plot_sst <- plot_sst %>%
      select(-c(2:7))
        
    
    # and scale
    scaled_sst <- sst
    scaled_sst[,2:7] <- apply(scaled_sst[,2:7], 2, scale)
    
    # and get Jan-Jun mean
mean_scaled_sst <- data.frame(year = 1975:2022,
                              Jan_Jun_SST = rowMeans(scaled_sst[,2:7]))
  
## process bottom temperature data -------------------------------

    
    # estimate annual bottom temp while 
    # correcting for missing stations and variable sampling dates
    
    library(tidyverse)
    library(lubridate)
    library(mice)
    # library(Hmisc)
    library(mgcv)

    current_year = 2023
    
    haulinfo_path <- paste0("Y:/KOD_Survey/EBS Shelf/", current_year, "/Tech Memo/Data/Haul Data/")
    list.files(haulinfo_path, ignore.case = TRUE) %>%
      map_df(~read.csv(paste0(haulinfo_path, .x))) -> dat
    
    
    as.integer(count(data.frame(unique(dat$SURVEY_YEAR))))-5-> min.years # calculates the minimum number of years that must be present to use a station
    
    
    theme_set(theme_bw())
    
    head(dat)
    
    # restrict to BB district
    strata <- read.csv("./data/STRATA_RKC_NEWTIMESERIES.csv") %>%
      filter(STRATUM == 10)
    
    dat <- dat %>%
      filter(GIS_STATION %in% unique(strata$STATION_ID),
             HAUL_TYPE != 17) # dropping retows
    
    # plot what we have
    plot.dat <- dat %>%
      dplyr::select(SURVEY_YEAR, GIS_STATION, MID_LATITUDE, MID_LONGITUDE, HAUL_TYPE, GEAR_TEMPERATURE) %>%
      rename(LATITUDE = MID_LATITUDE,
             LONGITUDE = MID_LONGITUDE)
    
    
    ggplot(filter(plot.dat, SURVEY_YEAR <= 1985), aes(LONGITUDE, LATITUDE)) +
      geom_point() +
      facet_wrap(~SURVEY_YEAR)
    
    # aside - make a Julian day file 
    
    julian.dat <- dat %>%
      mutate(julian=yday(parse_date_time(START_TIME, "d-m-y", "US/Alaska"))) %>%
      rename(year = SURVEY_YEAR, station = GIS_STATION) %>%
      filter(HAUL_TYPE != 17) %>%
      dplyr::select(year, station, julian) 
    
    # write.csv(julian.dat, "./data/year_station_julian_day.csv", row.names = F)
    
    # 
    # # limit to stations sampled in 1975
    # stations <- plot.dat %>%
    #   filter(SURVEY_YEAR == 1975)
    # 
    # use.stations <- stations$GIS_STATION
    # 
    # length(use.stations) # 136
    
    # alternate approach - stations with no more than 5 missing years
    # count.stations <- plot.dat %>%
    #   filter(HAUL_TYPE != 17,         ###change to HAUL_TYPE != 17 ????
    #          !is.na(GEAR_TEMPERATURE)) %>%
    #   group_by(GIS_STATION) %>%
    #   summarise(count = n()) %>%
    #   arrange(desc(count))
    # 
    # View(count.stations)
    # 
    # # limit to 40+
    # 
    # use.stations <- count.stations %>%
    #   filter(count >= min.years)
    # 
    # length(use.stations$GIS_STATION) # 206
    # stations <- plot.dat %>%
    #   filter(SURVEY_YEAR == 1975)
    # 
    # use.stations <- stations$GIS_STATION
    # 
    # length(use.stations) # 136
    # 
    # plot.dat <- plot.dat %>%
    #   filter(GIS_STATION %in% use.stations$GIS_STATION, 
    #          HAUL_TYPE != 17,
    #          !is.na(GEAR_TEMPERATURE))
    # 
    # ggplot(plot.dat, aes(LONGITUDE, LATITUDE)) +
    #   geom_point() +
    #   facet_wrap(~SURVEY_YEAR)
    # 
    # check <- plot.dat %>%
    #   group_by(SURVEY_YEAR) %>%
    #   summarise(missing = length(use.stations$GIS_STATION)-n())
    # 
    # check
    # sum(check$missing)
    # 
    # 
    # check <- plot.dat %>%
    #   group_by(GIS_STATION) %>%
    #   summarise(count = n()) %>%
    #   arrange(count)
    
    # check
    
    # select only haul type = 3
    
    use.dat <- dat %>%
      mutate(julian=yday(parse_date_time(START_TIME, "d-m-y", "US/Alaska"))) %>%
      rename(bottom.temp = GEAR_TEMPERATURE, year = SURVEY_YEAR, station = GIS_STATION,
             latitude = MID_LATITUDE, longitude = MID_LONGITUDE) %>%
      filter(station %in% unique(strata$STATION_ID), 
             HAUL_TYPE != 17) %>%
      dplyr::select(julian, bottom.temp, year, station, latitude, longitude) 
    
    
    
    # look at these stations sampled more than once a year 
    
    check <- use.dat %>%
      group_by(year, station) %>%
      dplyr::summarize(count = n()) %>%
      filter(count > 1)
    
    check # none!
    
    
    ## impute missing values --------------------------
    
    # set up for multiple imputation
    
    # first, date
    dat.julian <- use.dat %>%
      dplyr::select(julian, year, station) %>%
      pivot_wider(names_from = station, values_from = julian) %>%
      arrange(year) %>%
      dplyr::select(-year)
    
    # check number missing
    f <- function(x) sum(is.na(x))
    
    check <- apply(dat.julian, 1, f)
    check
    
    sum(check) # 131 missing
    
    # % missing
    sum(check) / (48*length(unique(use.dat$station)))
    
    
    # examine correlations
    # r <- rcorr(as.matrix(dat.julian))$r 
    
    r <- cor(dat.julian, use = "pairwise.complete.obs") 
    r #Cross-year correlations between each station combination
    
    
    # choose 25 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
    pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
    diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE
    
    colnames(pred) <- rownames(pred) <- colnames(dat.julian) <- str_remove_all(colnames(pred), "-")
    
    blocks <- make.blocks(dat.julian)
    
    imp <- mice(data = dat.julian, method = "norm", m=100, predictorMatrix = pred, blocks = blocks) #Using Bayesian linear regression method
    saveRDS(imp, "./output/station_julian_day_imputations.RDS")
    imp_date <- readRDS("./output/station_julian_day_imputations.RDS")
    
    # are there NAs in complete(imp)?
    
    check <- is.na(complete(imp_date))
    
    sum(check) # 0

    # also create df to save mean annual temp and sampling day for each imputed temperature data set
    
    imputed.day <- data.frame()
    
    
    for(i in 1:100){
      
      # i <- 1
      temp <- complete(imp_date, action = i) %>%
        mutate(year = c(seq(1975, 2019, 1), seq(2021, current_year, 1))) %>%
        pivot_longer(cols = -year,
                     names_to = "station",
                     values_to = "day")
      
      imputed.day <- rbind(imputed.day,
                           temp)
    }
    
    mean.imputed.day <- as.data.frame(imputed.day %>%
      dplyr::group_by(year, station) %>%
      dplyr::summarize(day = round(mean(day))))
    
    str(mean.imputed.day)

    ## now impute temperature ----------------------
    # first, clean up dat
    dat.temp <- use.dat %>%
      select(bottom.temp, year, station) %>%
      pivot_wider(names_from = station, values_from = bottom.temp) %>%
      arrange(year) %>%
      select(-year)
    
    # check number missing
    f <- function(x) sum(is.na(x))
    
    check <- apply(dat.temp, 1, f)
    check
    sum(check) # 453 missing!
    
    # examine correlations
    r <- cor(dat.temp, use = "p") 
    r #Cross-year correlations between each station combination
    
    
    # choose 25 variables with highest absolute correlation (van Buuren & Oudshoorn recommend 15-25)
    pred <- t(apply(r,1, function(x) rank(1-abs(x))<=25))# T for the 30 strongest correlations for each time series
    diag(pred) <- FALSE # and of course, drop self-correlations - make the diagonal FALSE
    
    colnames(pred) <- rownames(pred) <- colnames(dat.temp) <- str_remove_all(colnames(pred), "-")
    
    blocks <- make.blocks(dat.temp)
    
    imp <- mice(data = dat.temp, method = "norm", m=100, predictorMatrix = pred, blocks = blocks) #Using Bayesian linear regression method
    
    saveRDS(imp, "./output/station_bottom_temp_imputations.RDS")
    imp_temp <- readRDS("./output/station_bottom_temp_imputations.RDS")
    
    
    # are there NAs in complete(imp)?
    
    check <- is.na(complete(imp_temp))
    
    sum(check) # 0 !
    
    
    imputed.temp <- data.frame()
    
    
    for(i in 1:100){
      
      # i <- 1
      temp <- complete(imp_temp, action = i) %>%
        mutate(year = c(seq(1975, 2019, 1), seq(2021, current_year, 1))) %>%
        pivot_longer(cols = -year,
                     names_to = "station",
                     values_to = "temp")
      
      imputed.temp <- rbind(imputed.temp,
                           temp)
    }
    
    mean.imputed.temp <- as.data.frame(imputed.temp %>%
                                        dplyr::group_by(year, station) %>%
                                        dplyr::summarize(temp = mean(temp)))
    
    str(mean.imputed.temp)
    
  # combine bottom temp day and year for analysis  
combined.dat <- left_join(mean.imputed.temp, mean.imputed.day) %>%
  mutate(station = as.factor(station),
         year = as.factor(year))
  

  # fit model
mod <- gam(temp ~ year + s(day, by = station, k = 4), data = combined.dat)

summary(mod)    

# and predict the annual values with date kept constant among years
mean.station.day <- mean.imputed.day %>%
  group_by(station) %>%
  summarize(day = round(mean(day)))

new.dat <- data.frame(year = rep(unique(combined.dat$year), each = 136),
                      station = mean.station.day$station,
                      day = mean.station.day$day)

pred.temp <- predict(mod, newdata = new.dat, se = F)

bottom.temp <- new.dat %>%
  mutate(temp = pred.temp,
         year = as.numeric(as.character(year))) %>%
  group_by(year) %>%
  summarise(bottom_temp = mean(temp))

ggplot(bottom.temp, aes(year, bottom_temp)) +
  geom_point() +
  geom_line()

# and scale
scaled_bottom_temp <- bottom.temp %>%
  mutate(bottom_temp = scale(bottom_temp))

#### save the ice and temp time series in original units for plotting----------------------

raw_dat <- left_join(plot_ice, plot_sst) %>%
  left_join(., bottom.temp) %>%
  rename(June_bottom_temp = bottom_temp)

write.csv(raw_dat, "./data/sst_temp_original_units.csv", row.names = F)

#### combine the different time series----------------------

combined_data <- data.frame(year = 1975:2022) %>%
  left_join(., scaled_bottom_temp) %>%
  left_join(., ice) %>%
  left_join(., mean_scaled_sst)

plot_combined <- combined_data %>%
  pivot_longer(cols = -year)

ice_plot <- plot_combined %>%
  filter(name %in% c("Jan_Feb_ice", "Mar_Apr_ice"))

ggplot(ice_plot, aes(year, value, color = name)) +
  geom_hline(yintercept = 0, color = "dark grey") +
  geom_point() + 
  geom_line() +
  scale_color_manual(values = cb[c(2,4)])

temp_plot <- plot_combined %>%
  filter(name %in% c("bottom_temp", "Jan_Jun_SST"))

ggplot(temp_plot, aes(year, value, color = name)) +
  geom_hline(yintercept = 0, color = "dark grey") +
  geom_point() + 
  geom_line() +
  scale_color_manual(values = cb[c(6,7)])

#### fit DFA--------------------------------------
library(MARSS)

dfa.dat <- plot_combined %>%
  pivot_wider(names_from = name, values_from = value) %>% 
  arrange(year) %>%
  dplyr::select(-year) %>%
  t()

colnames(dfa.dat) <- 1975:2022

# and plot correlations
cors <- cor(t(dfa.dat), use = "p")
diag(cors) <- 0

max(cors)
min(cors) 

# set up forms of R matrices
levels.R = c("diagonal and equal",
             "diagonal and unequal",
             "equalvarcov",
             "unconstrained")

model.data = data.frame()

# changing convergence criterion to ensure convergence
cntl.list = list(minit=200, maxit=20000, allow.degen=FALSE, conv.test.slope.tol=0.1, abstol=0.0001)

# fit models & store results
for(R in levels.R) {
  for(m in 1) {  # find best single-trend model
    
    dfa.model = list(A="zero", R=R, m=m)
    
    kemz = MARSS(dfa.dat, model=dfa.model,
                 form="dfa", z.score=TRUE, control=cntl.list)
    
    model.data = rbind(model.data,
                       data.frame(R=R,
                                  m=m,
                                  logLik=kemz$logLik,
                                  K=kemz$num.params,
                                  AICc=kemz$AICc,
                                  stringsAsFactors=FALSE))
    
    assign(paste("kemz", m, R, sep="."), kemz)
  } # end m loop
} # end R loop

# calculate delta-AICc scores, sort in descending order, and compare
model.data$dAICc <- model.data$AICc-min(model.data$AICc)
model.data <- model.data %>%
  arrange(dAICc)
model.data # unconstrained is the best model 

# fit the best model
model.list = list(A="zero", m=1, R="unconstrained") # second-best model 

mod = MARSS(dfa.dat, model=model.list, z.score=TRUE, form="dfa", control=cntl.list)

# plot loadings and trend

CI <- MARSSparamCIs(mod)

plot.CI <- data.frame(names=rownames(dfa.dat),
                      mean=CI$par$Z[1:4],
                      upCI=CI$par.upCI$Z[1:4],
                      lowCI=CI$par.lowCI$Z[1:4])

dodge <- position_dodge(width=0.9)


plot.CI$names <- reorder(plot.CI$names, CI$par$Z[1:4])

loadings.plot <- ggplot(plot.CI, aes(x=names, y=mean)) +
  geom_bar(position=dodge, stat="identity", fill=cb[6]) +
  geom_errorbar(aes(ymax=upCI, ymin=lowCI), position=dodge, width=0.5) +
  ylab("Loading") +
  xlab("") +
  theme(axis.text.x  = element_text(angle=60, hjust=1,  size=9), legend.title = element_blank(), legend.position = 'top') +
  geom_hline(yintercept = 0)

# plot trend
trend <- data.frame(t=1975:2022,
                    estimate=as.vector(mod$states),
                    conf.low=as.vector(mod$states)-1.96*as.vector(mod$states.se),
                    conf.high=as.vector(mod$states)+1.96*as.vector(mod$states.se))


trend.plot <- ggplot(trend, aes(t, estimate)) +
  theme_bw() +
  geom_line(color=cb[6]) +
  geom_hline(yintercept = 0) +
  geom_point(color=cb[6]) +
  geom_ribbon(aes(x=t, ymin=conf.low, ymax=conf.high), linetype=2, alpha=0.1, fill=cb[6]) + xlab("") + ylab("Temperature index")


# save
png("./Figs/temperature_DFA_loadings_trend.png", width = 9, height = 3.5, units = 'in', res = 300)
ggpubr::ggarrange(loadings.plot,
                  trend.plot,
                  ncol = 2,
                  widths = c(0.35, 0.65),
                  labels = "auto")

dev.off()

# and save loadings and trend
write.csv(plot.CI, "./Output/dfa_loadings.csv", row.names = F)
write.csv(trend, "./Output/dfa_trend.csv", row.names = F)

## plot predicted / observed values

DFA_pred <- print(predict(mod)) %>%
  mutate(year = rep(1975:2022, 4))

ggplot(DFA_pred, aes(estimate, y)) +
  geom_text(aes(label = year)) +
  geom_smooth(method = "lm", se = F) +
  facet_wrap(~.rownames, ncol = 2, scale = "free") +
  labs(x = "Estimated", y = "Observed")


## plot regressions of bottom temp vs other three 

raw_dat <- read.csv("./data/sst_temp_original_units.csv") %>%
  pivot_longer(cols = c(-year, -June_bottom_temp))

ggplot(raw_dat, aes(value, June_bottom_temp)) +
  facet_wrap(~name, scales = "free_x") +
  geom_text(aes(label = year)) +
  geom_smooth(method = "l,", se = F)
 # that doesn't look bad - it looks like the March_Apr ice is the best predictor, and explains the poor
 # fit for Jan_Jun_SST

## fit as a GAM and examine residuals

raw_dat <- read.csv("./data/sst_temp_original_units.csv")

mod <- gam(June_bottom_temp ~ s(Jan_Feb_plot_ice, k = 4) + s(Mar_Apr_plot_ice, k = 4) + s(Jan_Jun_SST, k = 4),
           data = raw_dat)

summary(mod)

raw_dat <- na.omit(raw_dat)

raw_dat$resid <- resid(mod) 

ggplot(raw_dat, aes(year, resid)) +
  geom_point() +
  geom_line()

# some pattern there!