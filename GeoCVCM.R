## Main code to call functions to run GeoCVCM model
## Code author: Kalai Ramea
## Date: March 22, 2016

## set working directory

 cities <- c( "Oakland", "Berkeley", "Alameda", "Belmont",
             "Burlingame", "Daly City", "East Palo Alto", "Emeryville",
             "Foster City", "Fremont", "Hayward", "Los Altos", "Menlo Park",
             "Millbrae", "Milpitas", "Newark", "Palo Alto",
             "Redwood City", "San Bruno","San Leandro", "San Mateo",
             "Saratoga", "South San Francisco", "Sunnyvale", "Union City",
             "Los Angeles", "Beverly Hills", "Santa Monica", "Torrance", "Glendale",
             "Richmond", "Albany", "El Cerrito", "Sausalito", "Marin",
             "Brisbane", "San Leandro", "San Francisco")
 
## Set the working directory
setwd("/Users/kalaivanikubendran/Documents/Kalai-PhD-Research-Documents/Spatial Consumer Choice Model/Data/")

## This function calls the nested logit core module to estimate the purchase probability of vehicles
source("nmnl.R")

## Libraries
library(plyr)
library(Amelia)
library(ggplot2)
library(stringr)
library(maptools)  
library(ggmap)
library(geosphere)
library(GISTools)
library(rgeos)
library(proj4)
library(reshape2)
library(rgdal)

## Read relevant data
vehprices <- read.csv('vehprices.csv')
eff <- read.csv('dA.csv') ## Non-electric vehicle technology efficiency in ggepm
elec_cd_eff <- read.csv('dBe.csv') ## Electricity CD (charge depletion) efficiency
elec_cs_eff <- read.csv('dB.csv') ## Electricity CS (charge sustaining) efficiency, this is used in PHEVs to measure blended range
constant <- read.csv('constant.csv') ## Calibration constant
storage <- read.csv('storage.csv') ## Fuel storage
veh_range <- read.csv('xd.csv') ## Vehicle Range
fuelprices <- read.csv('fuelprices.csv') ## Fuel prices
modelavail <- read.csv('modelAvailability.csv') ## Make and model diversity cost
risk <- read.csv('risk.csv') ## Risk premium
trips <- read.csv("survey_place.csv") ## CHTS data on trips
hh <- read.csv("survey_households.csv") ## CHTS data on households

acs_hh <- read.csv('ACS_14_5YR_DP04.csv') ## ACS data on households
zcta_zip <- read.csv('zip_to_zcta_2015.csv') ## Zip code vs ZCTA comparison data



##******************************************************************##
## The following set cleans and reduces the ACS data to home zipcodes

acs_hh <- acs_hh[,c(2,30,34)]
acs_hh <- data.frame(lapply(acs_hh, as.character), stringsAsFactors=FALSE)
acs_hh <- acs_hh[-1,]
acs_hh$home <- as.numeric(acs_hh[,2]) + as.numeric(acs_hh[,3])
acs_hh <- acs_hh[,-c(2,3)]
colnames(acs_hh) <- c("ZCTA", "home")
acs <- merge(acs_hh, zcta_zip, by="ZCTA")

##******************************************************************##
## Extracting income data

hh_income <- hh[,c("sampno", "income", "home_zipcode")] ## we are only interested in income data from CHTS households for now

## Categorize the income groups to match our categories
## <$20K, $20K to $50K, $50K to $75K, $75K to $150K, > $150K (5 categories)
## We need to map the income groups in the household data to these 5 categories

hh_income$income[hh_income$income == 1 | hh_income$income == 2] = 1
hh_income$income[hh_income$income == 3 | hh_income$income == 4] = 2
hh_income$income[hh_income$income == 5] = 3
hh_income$income[hh_income$income == 6 | hh_income$income == 7] = 4
hh_income$income[hh_income$income == 8 | hh_income$income == 9 | hh_income$income == 10] = 5

##******************************************************************##

## Points of interest data in California

poi_ca <- read.csv('poi_ca.csv')
poi_ca <- poi_ca[,-1]
poi <- poi_ca[,c(2,1)]
colnames(poi) <- c("Longitude", "Latitude")
poi <- data.frame(lapply(poi, as.character), stringsAsFactors=FALSE)
poi <- data.frame(lapply(poi, as.numeric))
poi <- na.omit(poi)

##******************************************************************##

## This set of code reads the shapefile data 

proj <- CRS("+proj=utm +zone=10 +datum=WGS84")
## Read the ZCTA shapefile
zcta <- readShapePoly("ca_zcta.shp", proj4string=proj)

zcta <- readOGR(dsn = ".", "ca_zcta")
print(proj4string(zcta))

## Convert units from decimal degrees to sq.ft.
zcta_mile <- spTransform(zcta, CRS("+init=epsg:2225"))


zcta.f <- fortify(zcta, region = "ZCTA5CE10")
colnames(zcta.f) <- c("long", "lat", "order", "hole", "piece", "ZCTA", "group")

##******************************************************************##

## This set of code extracts the Latitude, Longitude values of 
## electric charging and hydrogen stations in California.

ev_stns <- read.csv("alt_fuel_stations_ev.csv") ## Electric charging station data
h2_stns <- read.csv("alt_fuel_stations_h2.csv") ## Hydrogen station data

ev <- ev_stns[,c(26,25)]
h2 <- h2_stns[,c(26,25)]


##******************************************************************##

## Functions used in the model

## This function calculates the Vehicle-miles traveled share of each 
## zip code and income group
calc_vmt_share <- function(zip_value, income_grp){
  zip_trips <- subset(inc_trips, home_zipcode == zip_value)
  
  ## Modes 5, 6, 7 are coded as road trips in CHTS data
  zip_car_trips <- subset(zip_trips, mode == 5 | mode == 6 | mode == 7)
  zip_car_trips <- subset(zip_car_trips, income == income_grp)
  zip_car_trips_final <- zip_car_trips[,c("sampno", "trip_distance_miles")]
  
  ## Check missing values in the data
  sapply(zip_car_trips_final, function(x) sum(is.na(x)))
  #missmap(zip_car_trips_final, main = "Missing values vs observed")
  
  ## Remove missing values
  ## (This is case-by-case. If a lot of values are missing, then other course of action is suggested)
  zip_car_trips_final <- na.omit(zip_car_trips_final)
  
  trip_distance_count <- ddply(zip_car_trips_final, .(sampno), summarise, sum = sum(trip_distance_miles))
  
  #plot_title <- paste0("Daily VMT Distribution of ", zip_value)
  #hist(trip_distance_count[,2],breaks=25,col="blue",xlab="Miles", cex.main=2,cex.lab=1.5, cex.axis=2,main=plot_title, xlim=c(0,100))
  
  
  #### Annual VMT Calculation
  ### Multiply daily trip distance by 365 to get annual VMT of the sample
  
  trip_distance_count$annual_vmt <- trip_distance_count$sum * 365
  
  
  trip_distance_count$vmt_cat <- ifelse(trip_distance_count$annual_vmt <=5000, 1,
                                        ifelse(trip_distance_count$annual_vmt > 5000 & trip_distance_count$annual_vmt <= 10000, 2, 
                                               ifelse(trip_distance_count$annual_vmt > 10000 & trip_distance_count$annual_vmt <= 15000, 3,
                                                      ifelse(trip_distance_count$annual_vmt > 15000 & trip_distance_count$annual_vmt <= 20000, 4,
                                                             ifelse(trip_distance_count$annual_vmt > 20000 & trip_distance_count$annual_vmt <= 25000, 5,
                                                                    ifelse(trip_distance_count$annual_vmt > 25000 & trip_distance_count$annual_vmt <= 30000, 6, 7))))))
  
  vmt_tot <- count(trip_distance_count, "vmt_cat")
  
  categories <- as.data.frame(c(1:7))
  colnames(categories) <- "vmt_cat"
  
  vmt_share <- prop.table(as.matrix(vmt_tot[,2]),2)
  
  ## Group the VMT data per each category
  avg_vmt <- ddply(trip_distance_count, .(vmt_cat), summarise, mean = mean(annual_vmt))
  
  tot_vmt_shares <- cbind(vmt_share, as.matrix(avg_vmt[,2]))
  
  return(tot_vmt_shares)
}

## This function calculates the infrastructure availability 
## of charging station and hydrogen station at zip code level

inf_availability <- function(zip_value){
  
  inf_total <- vector(mode="numeric", 0) ## initialize vector
  
  ##Get the spatial dataframe of the zip code area
  polygon_zcta <- zcta_mile[zcta_mile$ZCTA5CE10 == zip_value,]
  
  ## This converts the electric charging stations to spatial object
  evdf <- SpatialPointsDataFrame(coords = ev, data = ev,
                                 proj4string = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))
  evdf_mile <- spTransform(evdf, CRS("+init=epsg:2225"))
  
  ## This converts the point of interest data to spatial object
  poidf <- SpatialPointsDataFrame(coords = poi, data = poi,
                                  proj4string = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))
  poidf_mile <- spTransform(poidf, CRS("+init=epsg:2225"))
  
  ## This converts the hydrogen station data to spatial object
  h2df <- SpatialPointsDataFrame(coords = h2, data = h2,
                                 proj4string = CRS("+proj=longlat +datum=NAD83 +no_defs +ellps=GRS80 +towgs84=0,0,0"))
  h2df_mile <- spTransform(h2df, CRS("+init=epsg:2225"))
  
  ## Count the number of stations, chargers, POI inside the zip code area
  evcounts = poly.counts(evdf,zcta)
  h2counts = poly.counts(h2df,zcta)
  poicounts = poly.counts(poidf, zcta)
  
  ## 5 mile buffer of zip code area
  buffer_area <- gBuffer(polygon_zcta, width=26400, byid=TRUE )
  
  #print(poly.counts(evdf_mile, polygon_zcta))
  
  z_area <- gArea(polygon_zcta)*3.587e-8 ## convert from sqft to sq.mile
  buf_area <- gArea(buffer_area)*3.587e-8
  
  #ev_per_mile_stn <- poly.counts(evdf_mile, polygon_zcta)/max(1, poly.counts(poidf_mile, polygon_zcta))
  ev_per_mile_stn <- poly.counts(evdf_mile, buffer_area)/max(1, poly.counts(poidf_mile, buffer_area))
  
  
  #print(poly.counts(evdf_mile, polygon_zcta))
  #print(ev_per_mile_stn/4)
  
  #h2_per_mile_stn <- poly.counts(h2df_mile, polygon_zcta)/z_area
  h2_per_mile_stn <- poly.counts(h2df_mile, buffer_area)/buf_area
  
  #print(poly.counts(h2df_mile, polygon_zcta))
  #print(h2_per_mile_stn/4)
  
  inf_ev <- max(min(ev_per_mile_stn, 1), 0.000001)
  inf_h2 <- max(min(h2_per_mile_stn/4, 1), 0.000001)
  
  ## availability values of electric and hydrogen stations
  inf_total <- cbind(inf_ev, inf_h2)
  
  return(inf_total)
  
}

## This function calculates the refueling inconvenience cost of vehicles
refueling_cost_calc <- function(vmt, h2avail){
  
  eff[is.na(eff)] <- 0
  eff_ldc <- eff[c(1:20),3]
  eff_ldt <- eff[c(21:40),3]
  
  storage[is.na(storage)] <- 0
  storage_ldc <- storage[c(1:20),3]
  storage_ldt <- storage[c(21:40),3]
  
  
  h2_full_refill <- (65/26.7537092 * 0.823389532 * (h2avail ^ (-0.4886)) / 60 * 22.85367791) * 3.56
  gasoline_full_refill <- (65/26.7537092 * 0.823389532 * (1 ^ (-0.4886)) / 60 * 22.85367791) * 3.56
  diesel_full_refill <- (65/26.7537092 * 0.823389532 * (0.3 ^ (-0.4886)) / 60 * 22.85367791) * 3.56
  
  
  annual_fuel_ldc = vmt * eff_ldc
  annual_fuel_ldt = vmt * eff_ldt
  
  refueling_cost_ldc <- c(rep(0,20))
  refueling_cost_ldt <- c(rep(0,20))
  
  for(t in 1:20){
    if(t == 1 | t == 4 | t == 7 | t == 8| t == 9){
      refueling_cost_ldc[t] <- (annual_fuel_ldc[t]/storage_ldc[t]/0.75)*gasoline_full_refill * 4.152978
      refueling_cost_ldt[t] <- (annual_fuel_ldt[t]/storage_ldt[t]/0.75)*gasoline_full_refill * 4.152978
    }
    else if(t == 2 | t == 5){
      refueling_cost_ldc[t] <- (annual_fuel_ldc[t]/storage_ldc[t]/0.75)*diesel_full_refill * 4.152978
      refueling_cost_ldt[t] <- (annual_fuel_ldt[t]/storage_ldt[t]/0.75)*diesel_full_refill * 4.152978
    }
    else if(t == 11){
      refueling_cost_ldc[t] <- (annual_fuel_ldc[t]/storage_ldc[t]/0.75)*h2_full_refill * 4.152978
      refueling_cost_ldt[t] <- (annual_fuel_ldt[t]/storage_ldt[t]/0.75)*h2_full_refill * 4.152978
    }
  }
  
  refueling_cost_total <- cbind(refueling_cost_ldc, refueling_cost_ldt, annual_fuel_ldc, annual_fuel_ldt)
  return(refueling_cost_total)
}

## This function calculates the range limitation cost of vehicles
range_anxiety_calc <- function(vmt, risk, evavail, home, work){
  
  veh_range[is.na(veh_range)] <- 0
  veh_range_ldc <- veh_range[c(1:20),4]
  veh_range_ldt <- veh_range[c(21:40),4]
  
  elec_cd_eff[is.na(elec_cd_eff)] <- 0
  elec_cd_eff_ldc <- elec_cd_eff[c(1:20),4]
  elec_cd_eff_ldt <- elec_cd_eff[c(21:40),4]
  
  elec_cs_eff[is.na(elec_cs_eff)] <- 0
  elec_cs_eff_ldc <- elec_cs_eff[c(1:20),4]
  elec_cs_eff_ldt <- elec_cs_eff[c(21:40),4]
  
  gamma_shape = 1.89756
  gamma_scale = (vmt/(365*gamma_shape))
  
  ## Penalty costs are given based on risk groups
  ## $10/day for EA, $20/day for EM, $50/day for LM
  if(risk=="EA"){
    anx_cost <- 10 
  }else if(risk=="EM"){
    anx_cost <- 20 
  }else{
    anx_cost <- 50 
  }
  
  
  Rrh_ldc <- vector(mode='numeric', 20)
  Rrp_ldc <- vector(mode='numeric', 20)
  Rrw_temp_ldc <- vector(mode='numeric', 20)
  Rrw_ldc <- vector(mode='numeric', 20)
  xrcd_ldc <- vector(mode='numeric', 20)
  dCDK_ldc <- vector(mode='numeric', 20)
  dCDKK_ldc <- vector(mode='numeric', 20)
  rental_days_ldc <-vector(mode='numeric', 20)
  rental_cost_ldc <- vector(mode='numeric', 20)
  kwh_ldc <- vector(mode='numeric', 20)
  term1_ldc <- vector(mode='numeric', 20)
  phev_annual_fuel_ldc <- vector(mode='numeric', 20)
  
  for(t in 1:20){
    Rrh_ldc[t] = pmin(veh_range_ldc[t], home * home_speed*home_hours/elec_cd_eff_ldc[t])
    Rrp_ldc[t] = pmin(veh_range_ldc[t], 6*2*evavail/elec_cd_eff_ldc[t])
    if(home == 1){
      Rrw_temp_ldc[t] = veh_range_ldc[t]-(Rrh_ldc[t]-gamma_scale/4)
    }else{
      Rrw_temp_ldc[t] = veh_range_ldc[t]-(0-gamma_scale/4)
    }
    Rrw_ldc[t] <- pmin(Rrw_temp_ldc[t], work*work_speed*work_hours/elec_cd_eff_ldc[t])
    xrcd_ldc[t] <- Rrh_ldc[t] + Rrw_ldc[t] + Rrp_ldc[t]
    dCDK_ldc[t] <- pgamma(xrcd_ldc[t], shape=gamma_shape, scale=gamma_scale) 
    dCDKK_ldc[t] <- pgamma(xrcd_ldc[t], shape=gamma_shape+1, scale=gamma_scale) 
    if(t == 15 | t == 16 | t == 17) {
      rental_days_ldc[t] <- (1-dCDK_ldc[t])*365
      rental_cost_ldc[t] <- rental_days_ldc[t] * anx_cost * 4.15297829332967
    }else{
      rental_cost_ldc[t] <- 0
    }
    kwh_ldc[t] <- ((elec_cd_eff_ldc[t] * (1-dCDK_ldc[t]) * xrcd_ldc[t]) + (elec_cd_eff_ldc[t] * gamma_shape * gamma_scale * dCDKK_ldc[t]))*365
    if(t == 7 | t == 8 | t == 9 ){
      term1_ldc[t] <- (gamma_shape * gamma_scale * (eff_ldc[t] + (elec_cs_eff_ldc[t] - eff_ldc[t])*dCDKK_ldc[t]))
      phev_annual_fuel_ldc[t] <- (((elec_cs_eff_ldc[t] - eff_ldc[t])*xrcd_ldc[t]*(1-dCDK_ldc[t])) + term1_ldc[t])*365
    }
  }
  
  
  
  Rrh_ldt <- vector(mode='numeric', 20)
  Rrp_ldt <- vector(mode='numeric', 20)
  Rrw_temp_ldt <- vector(mode='numeric', 20)
  Rrw_ldt <- vector(mode='numeric', 20)
  xrcd_ldt <- vector(mode='numeric', 20)
  dCDK_ldt <- vector(mode='numeric', 20)
  dCDKK_ldt <- vector(mode='numeric', 20)
  rental_days_ldt <-vector(mode='numeric', 20)
  rental_cost_ldt <- vector(mode='numeric', 20)
  kwh_ldt <- vector(mode='numeric', 20)
  term1_ldt <- vector(mode='numeric', 20)
  phev_annual_fuel_ldt <- vector(mode='numeric', 20)
  
  for(t in 1:20){
    Rrh_ldt[t] = pmin(veh_range_ldt[t], home * home_speed*home_hours/elec_cd_eff_ldt[t])
    Rrp_ldt[t] = pmin(veh_range_ldt[t], 6*2*evavail/elec_cd_eff_ldt[t])
    if(home == 1){
      Rrw_temp_ldt[t] = veh_range_ldt[t]-(Rrh_ldt[t]-gamma_scale/4)
    }else{
      Rrw_temp_ldt[t] = veh_range_ldt[t]-(0-gamma_scale/4)
    }
    Rrw_ldt[t] <- pmin(Rrw_temp_ldt[t], work*work_speed*work_hours/elec_cd_eff_ldt[t])
    xrcd_ldt[t] <- Rrh_ldt[t] + Rrw_ldt[t] + Rrp_ldt[t]
    dCDK_ldt[t] <- pgamma(xrcd_ldt[t], shape=gamma_shape, scale=gamma_scale) 
    dCDKK_ldt[t] <- pgamma(xrcd_ldt[t], shape=gamma_shape+1, scale=gamma_scale) 
    if(t == 15 | t == 16 | t == 17) {
      rental_days_ldt[t] <- (1-dCDK_ldt[t])*365
      rental_cost_ldt[t] <- rental_days_ldt[t] * anx_cost * 4.15297829332967
    }else{
      rental_cost_ldt[t] <- 0
    }
    kwh_ldt[t] <- ((elec_cd_eff_ldt[t] * (1-dCDK_ldt[t]) * xrcd_ldt[t]) + (elec_cd_eff_ldt[t] * gamma_shape * gamma_scale * dCDKK_ldt[t]))*365
    if(t == 7 | t == 8 | t == 9 ){
      term1_ldt[t] <- (gamma_shape * gamma_scale * (eff_ldt[t] + (elec_cs_eff_ldt[t] - eff_ldt[t])*dCDKK_ldt[t]))
      phev_annual_fuel_ldt[t] <- (((elec_cs_eff_ldt[t] - eff_ldt[t])*xrcd_ldt[t]*(1-dCDK_ldt[t])) + term1_ldt[t])*365
    }
  }
  
  rentalcost_kwh_total <- cbind(rental_cost_ldc, rental_cost_ldt, kwh_ldc, kwh_ldt, phev_annual_fuel_ldc, phev_annual_fuel_ldt)
  
  return(rentalcost_kwh_total)
}


##******************************************************************##



for (city in cities){
  
  ## Read the zip, zcta values for each city in the US
  zip_zcta <- read.csv("zip_to_zcta_2015.csv")
  
  ## Subset the zip, zcta values for California
  ca_zip_zcta <- subset(zip_zcta, STATE=="CA")
  
  ## Subset the zip, zcta values for the given city
  city_zip_zcta <- subset(ca_zip_zcta, PO_NAME == city & ZIP_TYPE ==  "ZIP Code area")
  
  ## Separate zip and ZCTA values into two different files
  city_zip <- city_zip_zcta[,1]
  city_zcta <- city_zip_zcta[,5]
  
  ## Add income data to trips dataset
  inc_trips <- merge(trips, hh_income, by ="sampno")
  
  risk_grps <- c("EA", "EM", "LM")
  risk_share <- c(0.08, 0.38, 0.54)
  
  year <- 2020
  home = 0
  work = 0
  home_speed = 6
  home_hours = 8
  work_speed = 6
  work_hours = 7
  fuel_price_data <- fuelprices[,"X2020"]
  eff[is.na(eff)] <- 0
  eff_ldc <- eff[c(1:20),3]
  eff_ldt <- eff[c(21:40),3]
  
  vehprice_current <- vehprices[,"X2020"]
  vehprice_current[is.na(vehprice_current)] <- 0
  
  colnames(modelavail) <- c("Tech", 2005:2050)
  modelavail_current <- modelavail["2020"]
  
  constant_current <- constant[,4]
  full_zip_pp <- vector(mode="numeric", 0)

  for(zip_value in city_zip){
    home_charge <- subset(acs, ZIP == zip_value)
    home_charge_share <- c(1- unlist(home_charge[,2]/100), unlist(home_charge[,2]/100))
    inf_total_stn = inf_availability(zip_value)
    final_pp <- vector(mode="numeric", 0)
    print(zip_value)
    zip_income <- subset(inc_trips, home_zipcode == zip_value)
    zip_income <- subset(zip_income, mode == 5 | mode == 6 | mode == 7)
    
    if(nrow(zip_income) == 0){
      final_tot_pp <-  as.data.frame(t(c(rep(NA, 40))))
      colnames(final_tot_pp) <- as.character(vehprices[,1])
    }else{
      inc_shares <- as.data.frame(ddply(zip_income, .(income), summarize, count = length(sampno)))
      inc_shares$income <- as.character(inc_shares$income)
      
      inc_grp <- as.data.frame(c(1:5))
      colnames(inc_grp) <- "income"
      inc_grp$income <- as.character(inc_grp$income)
      
      inc_shares_final <- join(inc_grp, inc_shares, by = "income")
      inc_shares_final[is.na(inc_shares_final)] <- 0
      if(sum(inc_shares_final$count) == 0){
        inc_shares_final$share <- 0
      }else{
        inc_shares_final$share <- inc_shares_final$count / sum(inc_shares_final$count)
      }
      
      inc_shares_final$avginc <- c(15000, 30000, 62500, 112500, 150000)
      
      for(ii in 1:5){
        if(inc_shares_final$share[ii] > 0){
          v <- calc_vmt_share(zip_value, ii)
          inf_total_stn = inf_availability(zip_value)
          h2avail <- inf_total_stn[2]
          evavail <- inf_total_stn[1]
          vmt_data <- v[,2]
          vmt_share <- v[,1]
          for(vmt in 1:length(vmt_data)){
            ref_cost <- refueling_cost_calc(vmt_data[vmt], h2avail)
            annual_fuel_ldc <- ref_cost[,3]
            annual_fuel_ldt <- ref_cost[,4]
            for(risk in 1:length(risk_grps)){
              for(h in 0:1){
                rental_cost <- range_anxiety_calc(vmt_data[vmt], risk_grps[risk], evavail, h, work)
                
                kwh_ldc <- rental_cost[,3]
                kwh_ldc[is.na(kwh_ldc)] <- 0
                
                kwh_ldt <- rental_cost[,4]
                kwh_ldt[is.na(kwh_ldt)] <- 0
                
                phev_annual_fuel_ldc <- rental_cost[,5]
                phev_annual_fuel_ldt <- rental_cost[,6]
                
                fuel_cost_ldc <- vector(mode="numeric", 20)
                fuel_cost_ldt <- vector(mode="numeric", 20)
                
                for (i in 1:20){
                  if(i == 2 | i == 5){
                    fuel_cost_ldc[i] <- fuel_price_data[2]*annual_fuel_ldc[i]* 4.152978
                    fuel_cost_ldt[i] <- fuel_price_data[2]*annual_fuel_ldt[i]* 4.152978
                  }else if(i == 3 | i == 6 ){
                    fuel_cost_ldc[i] <- 0.001 ## Placeholder for natural gas
                    fuel_cost_ldt[i] <- 0.001 ## Placeholder for natural gas
                  }else if(i == 10 | i == 11 | i == 12 | i == 13 | i == 14 ){
                    fuel_cost_ldc[i] <- fuel_price_data[4]*annual_fuel_ldc[i]* 4.152978
                    fuel_cost_ldt[i] <- fuel_price_data[4]*annual_fuel_ldt[i]* 4.152978
                  }else if(i == 15 | i == 16 | i == 17 ){
                    fuel_cost_ldc[i] <- 0 ## Electricity cost will be given for BEVs
                    fuel_cost_ldt[i] <- 0 ## Electricity cost will be given for BEVs
                  }else if(i == 7 | i == 8 | i == 9){
                    fuel_cost_ldc[i] <- fuel_price_data[1]*phev_annual_fuel_ldc[i]* 4.152978 ## Gasoline vehicles
                    fuel_cost_ldt[i] <- fuel_price_data[1]*phev_annual_fuel_ldt[i]* 4.152978 ## Gasoline vehicles
                  }
                  else {
                    fuel_cost_ldc[i] <- fuel_price_data[1]*annual_fuel_ldc[i]* 4.152978 ## Gasoline vehicles
                    fuel_cost_ldt[i] <- fuel_price_data[1]*annual_fuel_ldt[i]* 4.152978 ## Gasoline vehicles
                  }
                }
                
                electricity_cost_ldc <- c(rep(0,20))
                electricity_cost_ldt <- c(rep(0,20))
                
                for (i in 1:20){
                  if(i == 7 | i ==8 | i == 9| i == 12| i == 13 | i == 14 |
                     i == 15 | i == 16 | i == 17){
                    electricity_cost_ldc[i] <- fuel_price_data[3]*kwh_ldc[i]* 4.152978 ## Electricity cost will be given for BEVs
                    electricity_cost_ldt[i] <- fuel_price_data[3]*kwh_ldt[i]* 4.152978 ## Electricity cost will be given for BEVs
                  }
                }
                
                fuel_cost_current <- c(fuel_cost_ldc, fuel_cost_ldt)
                electricity_cost_current <- c(electricity_cost_ldc, electricity_cost_ldt)
                
                refueling_cost_current <- c(ref_cost[,1], ref_cost[,2]) 
                refueling_cost_current[is.na(refueling_cost_current)] <- 0
                
                rental_cost_current <- c(rental_cost[,1], rental_cost[,2]) 
                rental_cost_current[is.na(rental_cost_current)] <- 0
                
                income_pref_current <- (vehprice_current / inc_shares_final$avginc[ii])*5000
                
                tot_cost <- vehprice_current + constant_current + rental_cost_current + 
                  refueling_cost_current + fuel_cost_current + electricity_cost_current + modelavail_current + 
                  income_pref_current
                
                purchase_probabilities <- nmnl(as.data.frame(tot_cost))
                
                dem_share <- vmt_share[vmt]*risk_share[risk] * inc_shares_final$share[ii] * home_charge_share[h+1]
                
                purchase_probabilities <- dem_share * purchase_probabilities
                
                final_pp <- rbind(final_pp, purchase_probabilities) 
                #print(purchase_probabilities)
                
              }
            }
          }
        }
      }
      
      final_tot_pp <- colSums(final_pp)
      final_tot_pp <- as.data.frame(t(final_tot_pp))
    }
    
    
    
    #print(final_tot_pp)
    final_tot_pp <- cbind(zip_value, final_tot_pp)
    full_zip_pp <- rbind(final_tot_pp, full_zip_pp)
  }
  
  write.csv(full_zip_pp, paste0("pp_",city,".csv"))
  
}

