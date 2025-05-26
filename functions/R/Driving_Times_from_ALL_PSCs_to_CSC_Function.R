# This script contains functions to help with the driving time to the closest CSC hospital 
# using https://openrouteservice.org/

#' Creating function to generate driving time to the closest CSC hospital.
#' @param facilities  Data frame with hospitals.
#' @return final_output 
#' @import openrouteservice 
#' @import units
#' @import dplyr

options(openrouteservice.url = "http://ph-gb.iowa.uiowa.edu:8080/ors")


GetDrivingTimeAllPSCtoCSC<-function(facilities,warn = TRUE){
  
  
  # Checking warn
  if(!is.logical(warn) || length(warn) != 1){
    stop("warn must be either TRUE or FALSE")
  }
  
  
  # Checking for valid values of facilities
  if(!("data.frame" %in% class(facilities))){
    stop("facilities must be a data frame")
  }
  
  names(facilities) <- tolower(names(facilities))
  if(!all(c("longitude","latitude","type",'id')%in% names(facilities) )){
    stop("The facilities data frame must have 'longitude', 'latitude', 'type', and 'id'")
  }
  
  
  # Facilities with longitude and latitude
  facilities_CSC<-facilities%>%dplyr::filter(type=="CSC")
  facilities_PSC<-facilities%>%dplyr::filter(type=="PSC")
  
  facilities_by_type_long_lat_PSC<-facilities_PSC%>%dplyr::select(c("longitude","latitude"))
  facilities_by_type_long_lat_CSC<-facilities_CSC%>%dplyr::select(c("longitude","latitude"))
  
  
  
  # Great-circle function
  gcd.hf <- function(long1, lat1, long2, lat2) {
    long1 <- long1 * pi / 180
    lat1 <- lat1 * pi / 180
    long2 <- long2 * pi / 180
    lat2 <- lat2 * pi / 180
    R <- 6371 # Earth mean radius [km]
    delta.long <- (long2 - long1)
    delta.lat <- (lat2 - lat1)
    a <- sin(delta.lat/2)^2 + cos(lat1) * cos(lat2) * sin(delta.long/2)^2
    c <- 2 * asin(pmin(1,sqrt(a)))
    d <- R * c
    return(d) # Distance in km
  }
  
  
  # If a matrix or data.frame is supplied the first column should be the longitude and 
  # the second column should be the latitude. If a numeric vector of length two is 
  # supplied, the first element should be the longitude and the second element should 
  # be the latitude.
  
  
  GCD <- function(from, to, ...){
    require(units)
    UseMethod("GCD")
  }
  
  GCD.numeric <- function(from, to){
    dist <- gcd.hf(from[1], from[2], to[1], to[2]) * 1000
    units(dist) <- "m"
    dist
  }
  
  GCD.data.frame <- function(from, to){
    dist <- gcd.hf(from[, 1], from[, 2], to[, 1], to[, 2]) * 1000
    units(dist) <- "m"
    dist
  }
  
  GCD.matrix <- function(from, to){
    dist <- gcd.hf(from[, 1], from[, 2], to[, 1], to[, 2]) * 1000
    units(dist) <- "m"
    dist
  }
  
  GCD.sf <- function(from, to){
    # Checking for the correct coordinates 
    if(st_crs(from)[[1]] != "EPSG:4326" || st_crs(to)[[1]] != "EPSG:4326"){
      stop("from and to must both have a crs of 4326")
    }
    # Getting great circle distance
    Dist <- st_distance(from, to, by_element = TRUE, which = "Great Circle")
    return(Dist)
  }
  
  final_output<-list()
  for (i in 1:nrow(facilities_by_type_long_lat_PSC)) {
  facility_by_type_long_lat_PSC<-facilities_by_type_long_lat_PSC[i,]
  facility_PSC<-facilities_PSC[i,]
  
  # Reducing the number hospitals
  # Be careful with GCD, we are finding differences between one number and one column
  meters<-GCD(as.numeric(facility_by_type_long_lat_PSC),facilities_by_type_long_lat_CSC)
  # Using the new set of hospitals
  facility_by_type_long_lat_CSC<-facilities_by_type_long_lat_CSC[which(meters<=sort(meters)[5]),]
  facility_CSC<-facilities_CSC[which(meters<=sort(meters)[5]),]
  
  
  # Coordinates for openrouteservice
  coordinates<-rbind(facility_by_type_long_lat_PSC,facility_by_type_long_lat_CSC)
    
    
  # Sources and destinations for openrouteservice
  # Indices start on zero
  sources<- c(0:(nrow(facility_by_type_long_lat_PSC)-1))
  destinations<- c(nrow(facility_by_type_long_lat_PSC):(nrow(coordinates)-1))
    
  # Results from openrouteservice
  openrouteservice::ors_profile("car") 
  res <- openrouteservice::ors_matrix(coordinates, sources = sources, destinations = destinations, metrics = c("duration", "distance"),units="km")
  
  # Time matrix 
  time <- units::set_units(res$durations,"secs")
  units(time)<-"min"
  
  # Distance matrix
  distance<-units::set_units(res$distances,"km")   
  
  if(all(is.na(time))){next}else{
  # Index for the minimum value  
  index<-apply(time,1,function(x){which.min(x)})
  
  # Saving time and distance
  # Recovering information of the hospitals
  CSC<-data.frame(time=time[1,index],
                     distance=distance[1,index],
                     facility_CSC[index,])
    
  
  names(CSC)<-paste0(names(CSC),"_CSC")
  names(facility_PSC)<-paste0(names(facility_PSC),"_PSC")

  output<-cbind(CSC,facility_PSC)
  final_output<-rbind(final_output,output)}
  }
 
  return(final_output)
}

