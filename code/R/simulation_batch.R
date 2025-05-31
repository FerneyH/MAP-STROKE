# Loading stuff here
## R packages
library(sf)       # Maps 
library(dplyr)
library(readxl)
library(units)

## Simulation and prediction functions
source("functions/R/SimulateTimeFunctions.R")
source("functions/R/PredictionFunctions.R")
source("functions/R/OtherFunctions.R")
source("functions/R/SimulatePatientProfileFunctions.R")
Rcpp::sourceCpp("functions/Cpp/Multinomial.cpp")
Rcpp::sourceCpp("functions/Cpp/TruncatedPoisson.cpp")
Rcpp::sourceCpp("functions/Cpp/OtherFunctions.cpp")

## Datasets for simulating patient profiles and treatment times
BPQuants <- read.csv("datasets/BPPatientDistributions.csv")
OtherQuants <- read.csv("datasets/OtherPatientDistributions.csv")
timesState <- read.csv("datasets/SimulationTimesState.csv")
timesRUCA <- read.csv("datasets/SimulationTimesRUCA.csv")
times <- list("state" = timesState, "RUCA" = timesRUCA)
EVT_logit <- read.csv("datasets/EVT_Receipt_Model_GB_2024_11_22.csv")

## Model object
model <- readRDS("datasets/BothModels.RDS")

# Getting RUCA dataset
RUCA <- read_xlsx("datasets/ruca2010revised.xlsx", sheet = 1, skip = 1)
colnames(RUCA) <- c("CountyFips", "State", "County", "TractFips", 
                    "RUCA1", "RUCA2", "TractPop", "Area", "PopDens")
RUCA$RUCA <- RUCA$RUCA1

# Setting chunk size
chunksize <- 1000

# Setting number of mRS replicates to generate
nmRS <- 100

DoEverything <- function(fips, chunksize = 1000, nmRS = 100, 
                         results_dir = "results/", generated_dir = "generated/", ...){
  # Getting name for state from fips code
  state <- unique(times$state$state_name[times$state$state_code == fips])
  state <- toupper(state)
  
  # Getting datasets for strokes
  ## Checking if this state has simulations
  if(!file.exists(paste0(generated_dir, state, "_Random_Strokes.rda"))){
    # stop(paste0("fips code ", state, " not supported"))
    # exit function since we don't have data for this state
    return()
  }
  
  ## Loading simulated times
  load(paste0(generated_dir, state, "_Random_Strokes.rda"))
  load(paste0(generated_dir, state, "_closest_CSC.rda"))
  load(paste0(generated_dir, state, "_closest_five_PSCs.rda"))
  load(paste0(generated_dir, state, "_from_PSC_to_CSC.rda"))
  
  ## Replacing NA states
  ### closest CSC
  if(any(is.na(closest_CSC$state))){
    closest_CSC$state <- ifelse(is.na(closest_CSC$state) & (closest_CSC$city == "WASHINGTON"), 
                                "DISTRICT OF COLUMBIA", 
                                closest_CSC$state)
  }
  
  ### closest 5 PSCs
  if(any(is.na(closest_five_PSCs$state))){
    closest_five_PSCs$state <- ifelse(is.na(closest_five_PSCs$state) & 
                                      (closest_five_PSCs$city == "WASHINGTON"), 
                                      "DISTRICT OF COLUMBIA", 
                                      closest_five_PSCs$state)
  }
  
  ### from PSC to CSC
  #### PSC
  if(any(is.na(from_PSC_to_CSC$state_PSC))){
    from_PSC_to_CSC$state_PSC <- ifelse(is.na(from_PSC_to_CSC$state_PSC) & 
                                        (from_PSC_to_CSC$city_PSC == "WASHINGTON"), 
                                        "DISTRICT OF COLUMBIA", 
                                        from_PSC_to_CSC$state_PSC)
  }
  #### CSC
  if(any(is.na(from_PSC_to_CSC$state_CSC))){
    from_PSC_to_CSC$state_CSC <- ifelse(is.na(from_PSC_to_CSC$state_CSC) & 
                                          (from_PSC_to_CSC$city_CSC == "WASHINGTON"), 
                                        "DISTRICT OF COLUMBIA", 
                                        from_PSC_to_CSC$state_CSC)
  }
  
  # Mothership
  ## Getting RUCA codes for stroke locations
  msPointsDF <- data.frame("latitude" = Strokes$latitude, 
                           "longitude" = Strokes$longitude)
  msStates <- substr(Strokes$cnty_fips, 1, 2)
  msCounties <- substr(Strokes$cnty_fips, 3, 5)
  msRUCA <- lonlat_to_RUCA(RUCA, pointsDF = msPointsDF, states = msStates, 
                           counties = msCounties)
  
  ## Getting rurality for each stroke location
  msRurality <- ifelse(msRUCA < 2, "Urban", 
                       ifelse(msRUCA < 4, "Suburban", "Rural"))
  
  ## Geting 4-category RUCA codes
  msRUCA <- 1 + (msRUCA > 3) + (msRUCA > 6) + (msRUCA >= 10) + (msRUCA == 99)
  msRUCA[msRUCA == 5] <- NA
  
  ## Making dataset to perform simulations and get predictions
  ### setting hospRUCA to -1 since it is only used in drip-and-ship
  msdata <- inner_join(Strokes, closest_CSC, 
                       by = c("longitude" = "longitude_origen", 
                              "latitude" = "latitude_origen")) |> 
    mutate(state_CSC = state, 
           state = substr(cnty_fips, 1, 2),
           Type = case_when(category == "LVO" ~ "lvo", 
                            category == "Non-LVO" ~ "not lvo", 
                            category == "Hemorrhagic" ~ "hemorrhagic", 
                            category == "Mimic" ~ "mimic", 
                            TRUE ~ "other"), 
           TimeToHosp = set_units(time / 60, NULL), regime = "ms", 
           TransferTime = 0, RUCA = msRUCA, hospRUCA = -1, 
           rurality = msRurality)
  
  ### Converting latitude.y and longitude to numeric values since this is sometimes an issue
  msdata$longitude.y <- as.numeric(msdata$longitude.y)
  msdata$latitude.y <- as.numeric(msdata$latitude.y)
  
  # Drip-and-Ship
  ## Making dataset to perform simulations and get predictions
  dsdata <- inner_join(Strokes, closest_five_PSCs, 
                       by = c("longitude" = "longitude_origen", 
                              "latitude" = "latitude_origen")) |> 
    mutate(state_PSC = state,
           state_CSC = from_PSC_to_CSC$state_CSC,
           state = substr(cnty_fips, 1, 2), 
           Type = case_when(category == "LVO" ~ "lvo", 
                            category == "Non-LVO" ~ "not lvo", 
                            category == "Hemorrhagic" ~ "hemorrhagic", 
                            category == "Mimic" ~ "mimic", 
                            TRUE ~ "other"), 
           TimeToHosp = set_units(time / 60, NULL), regime = "ds", 
           TransferTime = set_units(from_PSC_to_CSC$time_CSC / 60, NULL))
  
  ### Converting latitude.y and longitude to numeric values since this is sometimes an issue
  dsdata$longitude.y <- as.numeric(dsdata$longitude.y)
  dsdata$latitude.y <- as.numeric(dsdata$latitude.y)
  
  ## getting RUCA codes
  ### for stroke locations
  dsStrokePointsDF <- data.frame("latitude" = dsdata$latitude, 
                                 "longitude" = dsdata$longitude)
  dsStrokeStates <- substr(dsdata$cnty_fips, 1, 2)
  dsStrokeCounties <- substr(dsdata$cnty_fips, 3, 5)
  dsStrokeRUCA <- lonlat_to_RUCA(RUCA, pointsDF = dsStrokePointsDF, 
                                 states = dsStrokeStates, 
                                 counties = dsStrokeCounties)
  
  ## Getting rurality for each stroke location
  dsRurality <- ifelse(dsStrokeRUCA < 2, "Urban", 
                       ifelse(dsStrokeRUCA < 4, "Suburban", "Rural"))
  
  ## Getting 4-category RUCA codes
  dsStrokeRUCA <- 1 + (dsStrokeRUCA > 3) + (dsStrokeRUCA > 6) + 
    (dsStrokeRUCA >= 10) + (dsStrokeRUCA == 99)
  dsStrokeRUCA[dsStrokeRUCA == 5] <- NA
  
  dsdata$RUCA <- dsStrokeRUCA
  dsdata$rurality <- dsRurality
  
  ### for PSCs
  dsPSCPointsDF <- data.frame("latitude" = dsdata$latitude.y, 
                              "longitude" = dsdata$longitude.y)
  dsPSCStates <- dsdata$state_PSC
  dsPSCCounties <- NULL
  dsPSCRUCA <- lonlat_to_RUCA(RUCA, pointsDF = dsPSCPointsDF, 
                              states = dsPSCStates, 
                              counties = dsPSCCounties)
  
  ## Getting 4-category RUCA codes
  dsPSCRUCA <- 1 + (dsPSCRUCA > 3) + (dsPSCRUCA > 6) + 
    (dsPSCRUCA >= 10) + (dsPSCRUCA == 99)
  dsPSCRUCA[dsPSCRUCA == 5] <- NA
  
  dsdata$hospRUCA <- dsPSCRUCA
  
  # Merging datasets
  Data <- bind_rows(msdata, dsdata) |> 
    group_by(latitude, longitude) |> 
    mutate(stroke_ID = cur_group_id()) |> 
    ungroup() |> 
    arrange(stroke_ID)
  
  ## Setting MT and TPA to TRUE to signify that all patients 
  ## should receive these
  MT <- rep(TRUE, nrow(Data))
  TPA <- rep(TRUE, nrow(Data))
  
  ## Changing state_PSC and state_CSC to state when missing for now
  Data$state_PSC <- ifelse(is.na(Data$state_PSC), state, Data$state_PSC)
  Data$state_CSC <- ifelse(is.na(Data$state_CSC), state, Data$state_CSC)
  
  ## Performing simulations and getting predictions
  ### Need to make id variable that identifies each stroke patient
  ### Might need to change in the future if we change how many replicates
  id <- Data$stroke_ID
  Preds <- ChunkPerformSimulations(chunksize, data = Data, model = model, 
                          BPQuants = BPQuants, OtherQuants = OtherQuants, 
                          TimeQuants = times, TPA = TPA, MT = MT, id = id, 
                          warn = FALSE, nmRS = nmRS, ...)
  data <- Preds$data
  mRS_simul <- Preds$mRS_simul
  
  ### Caching results
  #### Creating results_dir if it doesn't exist
  if(!dir.exists(results_dir)){
    dir.create(results_dir)
  }
  
  #### Creating path to store results for this state
  path <- paste0(results_dir, state)
  if(!dir.exists(path)){
    dir.create(path)
  }
  
  ### Saving results
  saveRDS(data, file = paste0(path, '/data.rds'))
  saveRDS(mRS_simul, file = paste0(path, '/mRS_simul.rds'))
  #fwrite(data, file = paste0(path, '/data.csv'), showProgress = FALSE)
  #fwrite(mRS_simul, file = paste0(path, '/mRS_simul.csv'), showProgress = FALSE)
}

# Get fips code
args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
  stop("no fips code found")
} else {
  fips <- as.numeric(args[[1]])
}

# Do simulations
set.seed(123123 + fips)
generated_dir <- paste0("generated", "/")
results_dir <- paste0("results", "/")
DoEverything(fips, chunksize, nmRS, generated_dir = generated_dir, 
             results_dir = results_dir, showprogress = FALSE)
