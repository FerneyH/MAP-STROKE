# This script contains functions to help with the driving time to hospital 
# using https://openrouteservice.org/

#' Creating function to generate driving time to hospital.
#' @param random_stroke  Data frame with random stroke events.
#' @param facilities  Data frame with hospitals by state.
#' @param number_hostipals Number of closest hospitals to include. 
#' @param hospital_type The destination is PSC or CSC.
#' @return final_output 
#' @import openrouteservice 
#' @import units
#' @import dplyr

options(openrouteservice.url = "http://#####/ors")


GetDrivingTime<-function(random_stroke,facilities,number_hostipals=5, hospital_type="PSC", warn = TRUE){
  
  
  # Checking warn
  if(!is.logical(warn) || length(warn) != 1){
    stop("warn must be either TRUE or FALSE")
  }
  
  # Checking for valid values of random_stroke
  if(!("data.frame" %in%class(random_stroke))){
    stop("random_stroke must be a data frame")
  }
  
  if(!(nrow(random_stroke)==1)){
    stop("random_stroke must have one row")
  }
  
  names(random_stroke) <- tolower(names(random_stroke))
  if(!all(names(random_stroke) %in% c("cnty_fips", "longitude","latitude", "stroke_type","category", "display_name"))){
    stop("The random_stroke data frame must have 'cnty_fips', 'longitude','latitude', 'Stroke_type','category', 'display_name'")
  }
  
  if(!all(nchar(random_stroke$cnty_fips)==5)){
    stop("Invalid 'cnty_fips', this is a example of a valid 'cnty_fips': '01003'")
  }
  
  # Checking for valid values of facilities
  if(!("data.frame" %in% class(facilities))){
    stop("facilities must be a data frame")
  }
  
  names(facilities) <- tolower(names(facilities))
  if(!all(c("longitude","latitude","type",'id')%in% names(facilities) )){
    stop("The facilities data frame must have 'longitude', 'latitude', 'type', and 'id'")
  }
  
  # Checking for valid values of hospital_type
  if(!is.character(hospital_type)){
    stop("hospital_type must be either PSC or CSC")
  }
  
  if(!(hospital_type=="PSC"|| hospital_type=="CSC")){
    stop("hospital_type must be either PSC or CSC")
  }
  
  if((hospital_type=="CSC" && number_hostipals>1)){
    stop("number_hostipals must be 1")
  }
  
  
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
  
  
  # Facilities with longitude and latitude
  facilities<-facilities%>%dplyr::filter(type==hospital_type)
  
  
  facilities_by_type_long_lat<-facilities%>%dplyr::select(c("longitude","latitude"))
  random_stroke_event_long_lat<-random_stroke%>%dplyr::select(c("longitude","latitude"))%>%as.data.frame()
  
  # Reducing the number hospitals
  meters<-GCD(random_stroke_event_long_lat,facilities_by_type_long_lat)
  # Using the new set of hospitals
  facilities_by_type_long_lat<-facilities_by_type_long_lat[which(meters<sort(meters)[15]),]
  facilities<-facilities[which(meters<sort(meters)[15]),]

 
  # Coordinates for openrouteservice
  coordinates<-rbind(random_stroke_event_long_lat,facilities_by_type_long_lat)
    
    
  # Sources and destinations for openrouteservice
  # Indices start on zero
  sources<- c(0)
  destinations<- c(1:(nrow(coordinates)-1))
    
  final_output<-list()
    
  # Results from openrouteservice
  openrouteservice::ors_profile("car") 
  res <- openrouteservice::ors_matrix(coordinates, sources = sources, destinations = destinations, metrics = c("duration", "distance"),units="km")
  
  # Time matrix where columns represent the hospitals and rows represent the random stroke event
  time <- units::set_units(res$durations,"secs")
  units(time)<-"min"
  
  
  # Distance matrix where columns represent the hospitals and rows represent the random stroke event
  distance<-units::set_units(res$distances,"km")   
    
  if(all(is.na(time))){return(
    data.frame(id="Openrouteservice not available",
               name="Openrouteservice not available",
               city="Openrouteservice not available",
               state="Openrouteservice not available",
               zip="Openrouteservice not available",
               address="Openrouteservice not available",
               latitude="Openrouteservice not available",
               longitude="Openrouteservice not available",
               type="Openrouteservice not available",
               time=0,
               distance=0,
               latitude_origen=random_stroke_event_long_lat$latitude[1],
               longitude_origen=random_stroke_event_long_lat$longitude[1])
  )}else{
  
  # Index for the closest hospitals 
  index<-as.matrix(apply(time,1, function(x)(sort(as.vector(x), index.return = TRUE)$ix[1:number_hostipals])))
    
  for (i in 1:number_hostipals) {
    indexh<-index[i,1]
    # Saving time and distance
    # Recovering information of the hospitals
    output<-data.frame(
              facilities[indexh,],
              time=time[1,indexh],
              distance=distance[1,indexh],
              latitude_origen=random_stroke_event_long_lat$latitude[1],
              longitude_origen=random_stroke_event_long_lat$longitude[1])
    
    final_output<-rbind(final_output,output)
  }
  
 return(final_output)
}}

