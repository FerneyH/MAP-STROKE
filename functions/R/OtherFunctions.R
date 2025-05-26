########################################
# Get RUCA Codes from Stroke Locations #
########################################

## Taken from (https://stackoverflow.com/questions/8751497/latitude-longitude-coordinates-to-state-code-in-r/8751965#8751965)
#' @param pointsDF A data.frame whose first column contains longitudes and whose 
#' second column contains latitudes.
#' @param state The two-digit FIPS code (string) of the state you want. Can also be state name or state abbreviation.
#' @param counties The three-digit FIPS code (string) of the county you'd like to 
#' subset for, or a vector of FIPS codes if you desire multiple counties. Can also be a county name or vector of names.
#' @param tracts An sf MULTIPOLYGON object with tracts.
#' @param name_col Name of a column in `tracts` that supplies the tracts' names.
#' @details 2010 census tracts are used to be compatible with the RUCA dataset.
lonlat_to_tract <- function(pointsDF, state, counties,
                            tracts = NULL,
                            name_col = "GEOID10") {
  if(is.null(tracts)){
    # Using 2010 census tracts to be compatible with RUCA dataset
    #tracts <- tigris::tracts(state = state, county = counties, progress_bar = FALSE, 
    #                         year = 2010)
    load(paste0("datasets/", state, "_tracts.rda"))
  }
  
  ## Convert points data.frame to an sf POINTS object
  pts <- st_as_sf(pointsDF, coords = 2:1, crs = st_crs(tracts))
  
  ## Performing spatial join
  tracts <- st_join(pts, tracts, join = st_nearest_feature)
  tracts[[name_col]]
}

## Taken from (https://stackoverflow.com/questions/8751497/latitude-longitude-coordinates-to-state-code-in-r/8751965#8751965)
#' @param pointsDF A data.frame whose first column contains longitudes and
#'                whose second column contains latitudes.
#' @param tracts An sf MULTIPOLYGON object with tracts.
#' @param name_col Name of a column in `tracts` that supplies the tracts' names.
all_lonlat_to_tract <- function(pointsDF, states, counties,
                                tracts = NULL,
                                name_col = "GEOID10") {
  res <- c()
  for(state in unique(states)){
    ind <- which(states == state)
    res[ind] <- lonlat_to_tract(pointsDF[ind, ], state, counties[ind], 
                                tracts = tracts, name_col = name_col)
  }
  return(res)
}

# Gets RUCA codes from long-latitude pairs
#' @param RUCA a data.frame with the RUCA code for each census tract.
#' @param ... arguments passed to [all_lonlat_to_tract].
lonlat_to_RUCA <- function(RUCA, ...){
  ## Getting census tracts
  tracts <- all_lonlat_to_tract(...)
  
  ## Getting RUCA codes                              
  RUCA$RUCA[match(tracts, RUCA$TractFips)]
}

###################################################
# Get Probability of Receiving EVT based on Times #
###################################################

# Could translate this to C++ for speed
get_mindices <- function(a,b){
  for (i in 1:length(a)){
    a[i] <- min(which(b >= a[i]))
  }
  a
}

# Gets Probability of Receiving EVT based on time to hospital
#' @param TimeToHosp A numeric object with the time it takes to get to a CSC.
EVTProb <- function(TimeToHosp){
  
  # Global variable: EVT_logit must exist
  idx <- get_mindices(TimeToHosp, EVT_logit$Time)
  idx[!is.finite(idx)] <- 1
  
  # Faster than plogis for some reason
  probs <- 1/(1+exp(-(EVT_logit[idx,]$LogitFit + 
                        rnorm(length(idx), 
                        mean = 0, 
                        sd = EVT_logit[idx,]$LogitSE))))
  
  ## truncating at 24 hours
  probs[TimeToHosp > 24] <- 0
  
  ## Returning results
  return(probs)
}
