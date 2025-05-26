##############################################
# Simulate Patient Characteristics from GWTG #
##############################################

#' Generate n sets of simulated of patient profiles
#' @param n the number of simulated patient profiles to generate.
#' @param type one of 'lvo', 'not lvo', 'hemorrhagic', or 'mimic'.
#' @param BPQuants a data.frame with the necessary BP variable quantiles from GWTG.
#' @param OtherQuants a data.frame with the necessary RACE and age quantiles from GWTG.
#' @return A data.frame with the following columns 
#' \item{`Type`}{ the specified stroke type}
#' \item{`sysBP`}{ the simulated systolic blood pressure used for predictions ((sysBP - 120) / 80)}
#' \item{`diasBP`}{ the simulated diasolic blood pressure used for predictions ((diasBP - 80) / 60)}
#' \item{`RACE_scale`}{ the simulated RACE score used for predictions (RACE - 6.347258)}
#' \item{`age`}{ the simulated age used for predictions ((years - 68.4208) / 10)}
#' \item{`sex`}{ the simulated sex}
GetNProfiles <- function(n, type, BPQuants, OtherQuants){
  # Checking n
  if(!is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n)){
    stop("n must be a positive integer")
  }
  
  # Checking for valid values of type
  type <- tolower(type)
  if(!(type %in% c("lvo", "not lvo", "hemorrhagic", "mimic"))){
    stop("type must be one of 'lvo', 'not lvo', 'hemorrhagic', or 'mimic'")
  }
  if(type == "lvo"){
    type <- "LVO"
  }else if(type == "not lvo"){
    type <- "Not LVO"
  }else if(type == "hemorrhagic"){
    type <- "Hemorrhagic"
  }else{
    type <- "Stroke Mimic"
  }
  
  # Getting profiles
  ## Sex
  sex <- sample(c("Male", "Female"), n, replace = TRUE)
  
  ## Blood pressure
  BPquantileind <- grep("^X.+", colnames(BPQuants))
  BPquantiles <- as.matrix(BPQuants[, BPquantileind])
  
  ### Systolic
  #### Removing min and max quantiles
  sysBPQuantiles <- BPquantiles[BPQuants$variable == "sysBP" & BPQuants$type == type, ]
  sysBPQuantiles <- sysBPQuantiles[-c(1, length(sysBPQuantiles))]
  sysBP <- sample(sysBPQuantiles, n, replace = TRUE)
  
  ### Diastolic
  #### Removing min and max quantiles
  diasBPInd <- BPQuants[BPQuants$variable == "diasBP" & BPQuants$type == type, ]
  diasBPInd <- cbind("sysBPLow" = diasBPInd$sysBPLow, "sysBPHi" = diasBPInd$sysBPHi)
  diasBPQuantiles <- BPquantiles[BPQuants$variable == "diasBP" & BPQuants$type == type, ]
  diasBPQuantiles <- diasBPQuantiles[, -c(1, ncol(diasBPQuantiles))]
  
  #### Doing conditional sampling based on sysBP
  diasBP <- sapply(sysBP, function(x){
    ind <- min(which(x <= diasBPInd[, 2]))
    sample(diasBPQuantiles[ind, ], 1)
  })
  
  ## Others
  Otherquantileind <- grep("^X.+", colnames(OtherQuants))
  Otherquantiles <- as.matrix(OtherQuants[, Otherquantileind])
  Otherquantiles <- Otherquantiles[, -c(1, ncol(Otherquantiles))]
  
  ### RACE
  RACE <- sample(Otherquantiles[OtherQuants$variable == "RACE" & 
                                OtherQuants$type == type, ], n, replace = TRUE)
  
  ### age
  age <- sample(Otherquantiles[OtherQuants$variable == "age" & 
                               OtherQuants$type == type, ], n, replace = TRUE)
  
  ## Returning results
  results <- data.frame("Type" = type, 
                        "sysBP" = sysBP, 
                        "diasBP" = diasBP, 
                        "RACE_scale" = RACE, 
                        "age" = age, 
                        "sex" = sex)
  rownames(results) <- NULL
  return(results)
}

#' Generate simulated patient profiles for each subject
#' @param data a data.frame with at least a column named 'Type'.
#' @param ... further arguments passed to [GetNProfiles]
#' @return A data.frame with the simulated patient profiles.
GetMultipleProfiles <- function(data, ...){
  # Checking data
  if(!is.data.frame(data)){
    stop("data must be a data.frame")
  }else if(!("Type" %in% colnames(data))){
    stop("data must have a column named 'Type'")
  }
  
  # Getting results
  Type <- data$Type
  res <- lapply(unique(Type), function(x){
    ind <- which(Type == x)
    cbind.data.frame("index" = ind, GetNProfiles(length(ind), x, ...))
  }) |> do.call(what = rbind)
  
  # Combining results with data
  datarownames <- rownames(data)
  data <- res[order(res$index), -1, drop = FALSE]
    
  ## Returning results
  return(data)
}

################################################################
# Simulate Patient Characteristics, Times, and Get Predictions #
################################################################

#' Generate simulated patient profiles, treatment times, and gets predictions
#' @param data a data.frame with at least columns named 'Type', 'TimeToHosp', 
#' 'TransferTime', 'regime', 'state', and 'RUCA'. Can also have a column labelled 
#' 'hospRUCA' for the RUCA of the first hospital with drip-and-ship.
#' @param model an object of class `TwoModels` or `Model`.
#' @param BPQuants see BPQuants argument of [GetNProfiles].
#' @param OtherQuants see OtherQuants argument of [GetNProfiles].
#' @param TimeQuants see times argument of [GetNTimes].
#' @param TPA a logical vector indicating whether subjects should receive TPA or not.
#' @param MT a logical vector indicating whether subjects should receive MT or not.
#' @param id a numeric vector denoting which stroke events should have the same 
#' patient characteristics.
#' @param nmRS the number of simulated mRS values to generate for each stroke.
#' @param ... further arguments passed to [GetNProfiles]
#' @return A data.frame with the simulated patient profiles appended to `data`.
PerformSimulations <- function(data, model, BPQuants, OtherQuants, TimeQuants, 
                               TPA, MT, id, nMRS, ...){
  # Checking data
  if(!is.data.frame(data)){
    stop("data must be a data.frame")
  }else if(!all(c("Type", "TimeToHosp", "TransferTime", "state", "regime", "RUCA") %in% colnames(data))){
    stop("data must have columns named 'Type','TimeToHosp', 'TransferTime', 'state', 'RUCA', and 'regime'")
  }
  
  # Checking model
  if(!is(model, "TwoModels")){
    stop("model must be an object of class TwoModels")
  }
  
  # Checking TPA and MT
  if(missing(TPA)){
    TPA <- rep("TRUE", nrow(data))
  }else if(length(TPA) != nrow(data)){
    stop("TPA must be missing or have length equal to the number of rows in data")
  }else if(!is.logical(TPA)){
    stop("TPA must be a logical vector")
  }
  if(missing(MT)){
    MT <- rep("TRUE", nrow(data))
  }else if(length(MT) != nrow(data)){
    stop("MT must be missing or have length equal to the number of rows in data")
  }else if(!is.logical(MT)){
    stop("MT must be a logical vector")
  }
  
  # Checking id
  if(missing(id)){
    id <- 1:nrow(data)
  }else if(length(id) != nrow(data)){
    stop("id must be missing or have length equal to the number of rows in data")
  }
  firstid <- !(duplicated(id))
  firstidvals <- id[firstid]
  idindices <- match(id, firstidvals)
  
  # Getting hospRUCA
  if("hospRUCA" %in% colnames(data)){
    hospRUCA <- data$hospRUCA
  }else{
    hospRUCA <- data$RUCA
  }
  
  # Adding patient profiles
  TypeInd <- which(colnames(data) == "Type")
  profileData <- GetMultipleProfiles(data[firstid, ], BPQuants = BPQuants, OtherQuants = OtherQuants)
  data <- cbind(profileData[idindices, ], data[, -TypeInd])
  
  # Adding times
  ## Adding time to hospital and transfer times
  data$TravelTime <- data$TimeToHosp
  regime <- data$regime
  TimeToHospInd <- which(colnames(data) == "TimeToHosp")
  timesData <- GetTreatmentTimes(data$TimeToHosp, 
                                 regime = regime, 
                                 state_PSC = data$state_PSC,
                                 state_CSC = data$state_CSC,
                                 RUCA = data$RUCA, hospRUCA = hospRUCA, 
                                 times = TimeQuants, id = id, TPA = TPA, MT = MT, 
                                 TransferTime = data$TransferTime)
  data <- cbind(timesData[, -(10:14)], data[, -TimeToHospInd])
  
  # Getting probability of EVT
  EVTProbVec <- EVTProb(data$TimeToCSC)
  
  # Modifying TPA if their time is 0
  data$TPA <- (data$TPATime > 0)
  rownames(data) <- NULL
  
  # Temporarily setting MT to TRUE
  data$MT <- (data$MTTime > 0)
  
  # Getting predictions
  ## Getting predicted stroke types
  PredictedTypes <- predict(model$diagnosticModel, newdata = data, type = "draws")
  
  
  ## Getting predicted stroke types
  mRSPredsUnknown <- colMeans(predict(model, data, 
                              predsDiagnostic = PredictedTypes, type = "draws", 
                              EVTProb = EVTProbVec))
  
  ## Getting posterior mean of PredictedTypes
  TypeNames <- dimnames(PredictedTypes)[[3]]
  PredictedTypes <- cubeMeanCpp(PredictedTypes)
  colnames(PredictedTypes) <- TypeNames
  
  # Deciding if subjects receive EVT based on their probabilities
  data$MT <- as.logical(rbinom(length(EVTProbVec), 1, prob = EVTProbVec)) * (data$Type == "LVO")
  data$MTTime <- data$MTTime * data$MT
  
  ## Adding predictions and EVTProb
  data <- cbind(data, PredictedTypes, "PredsUnknown" = mRSPredsUnknown, 
                "EVTProb" = EVTProbVec)
  
  ## Getting simulated mRS values
  mRS_simul <- simulatemRS(model$outcomeModel, data, nmRS)
  
  ## Returning results
  return(list("data" = data, "mRS_simul" = mRS_simul))
}

#' Generate simulated patient profiles, treatment times, and gets predictions in chunks
#' @param chunksize a positive integer denoting the chunk size.
#' @param data see data argument of [PerformSimulations].
#' @param TPA see TPA argument of [PerformSimulations].
#' @param MT see MT argument of [PerformSimulations].
#' @param id see id argument of [PerformSimulations].
#' @param showprogress a logical value indicating if a progress bar should be displayed.
#' @param ... further arguments passed to [PerformSimulations]
#' @return A data.frame with the simulated patient profiles appended to `data`.
ChunkPerformSimulations <- function(chunksize, data, TPA, MT, id, showprogress = TRUE, 
                                    ...){
  # Making id if it is missing
  if(missing(id)){
    id <- 1:nrow(data)
  }
  maxrep <- max(table(id))
  
  # Getting chunks
  chunks <- ceiling((id * maxrep) / (chunksize))
  nchunks <- length(unique(chunks))
  
  # Starting progress bar
  if(showprogress){
    pb <- txtProgressBar(0, nrow(data), style = 3)
    on.exit(close(pb))
  }
  
  # Doing stuff in chunks
  res <- lapply(1:nchunks, function(i){
    ## Getting indices for current chunk
    ind <- which(chunks == i)
    
    ## Getting results
    results <- PerformSimulations(data[ind, , drop = FALSE], TPA[ind], MT[ind], 
                                  id[ind], ...)
    
    ## Updating progress bar
    if(showprogress){
      setTxtProgressBar(pb, sum(chunks <= i))
    }
    # Returning results
    return(results)
  })
  ## Separating components
  finalres <- list("data" = bind_rows(lapply(res, function(x) x$data)), 
                   "mRS_simul" = bind_rows(lapply(res, function(x) x$mRS_simul)))
  
  return(finalres)
}
