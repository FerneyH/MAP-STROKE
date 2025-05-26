# This script contains functions to help with the simulation of treatment times 
# based on results from GWTG

#' Checks the supplied state
#' @param state a string with the name of the state or abbreviation. Can also supply a 
#' numeric value indicating the fips code.
#' @param times see times argument of [GetNTimes].
#' @param warn see warn argument of [GetNTimes].
#' @return The rows from times that match with the supplied state.
checkState <- function(state, times, warn){
  # Checking state
  if(is.numeric(state) || grepl('^(?=.)([+-]?([0-9]*)(\\.([0-9]+))?)$', state, perl = TRUE)){
    ## Converting character to numeric
    if(is.character(state)){
      state <- as.numeric(state)
    }
    
    ## converting fips to state
    if(nchar(state) > 2){
      state <- as.numeric(substr(state, 1, 2))
      if(warn){
        warning("only using first 2 digits from state for fips code")
      }
    }else if(nchar(state) == 0){
      stop("empty numeric value provided for state")
    }
    StateInd <- which(times$state$state_code == state)
    if(length(StateInd) == 0){
      stop(paste0("supplied state(" , state, ") is not supported"))
    }
    StateTimes <- times$state[StateInd, , drop = FALSE]
  }else if(!is.character(state) || length(state) != 1){
    stop("state must be a character vector of length one")
  }else if(nchar(state) == 2){
    StateInd <- which(times$state$state == state)
    if(length(StateInd) == 0){
      stop(paste0("supplied state(" , state, ") is not supported"))
    }
    StateTimes <- times$state[StateInd, , drop = FALSE]
  }else{
    s <- strsplit(state, " ")[[1]]
    state <- paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)), sep="", collapse=" ")
    ## Checking for DC
    if(state == "District Of Columbia"){
      state <- "District of Columbia"
    }
    
    StateInd <- which(times$state$state_name == state)
    if(length(StateInd) == 0){
      stop(paste0("supplied state(" , state, ") is not supported"))
    }
    StateTimes <- times$state[StateInd, , drop = FALSE]
  }
  return(StateTimes)
}

#' Generate n sets of simulated of times
#' @param n the number of sets of simulated times to generate.
#' @param regime one of 'ds' or 'ms' to denote the desired treatment regime.
#' @param state_PSC a character vector with the name of the state or abbreviation. Can also supply a 
#' numeric value indicating the fips code. This is for the PSC.
#' @param state_CSC a character vector with the name of the state or abbreviation. Can also supply a 
#' numeric value indicating the fips code. This is for the CSC.
#' @param RUCA a string or numeric value with the RUCA code for the stroke location.
#' @param hospRUCA a string or numeric value with the RUCA code for the first hospital 
#' in drip-and-ship. Can also be a vector of these values with a length of n.
#' @param times a data.frame with the time necessary time variables from GWTG.
#' @param TPA a logical value indicating if the simulated subjects should receive TPA.
#' @param warn a logical value indicating whether warnings should be displayed.
#' @param timeToHosp NULL or a numeric vector with a time to hospital after EMS 
#' attention (hours), if this is not supplied, then a time of 0 is assumed. 
#' @param TPACutoff the specified cutoff used to decide if a patient should receive 
#' TPA. This is only used if timeToHosp is supplied.
#' @param id a numeric vector denoting which stroke events should have the same 
#' time from LKW to EMS attention.
#' @return A data.frame with the following columns 
#' \item{`LKWToEMS`}{ the simulated time from last known well to EMS attention (hours)}
#' \item{`TimeToHosp`}{ the simulated time from last known well to hospital arrival (hours)}
#' \item{`TPA`}{ the simulated time from hospital arrival to TPA (hours)}
#' \item{`DIDO`}{ the simulated door-in-door-out time (hours) if drip-and-ship and 0 otherwise}
#' \item{`MT`}{ the simulated time from hospital arrival to groin puncture (hours) if drip-and-ship or the simulated time from TPA to groin puncture (hours) if mothership}
GetNTimes <- function(n, regime, state_PSC, state_CSC, RUCA, hospRUCA, times, TPA = TRUE, warn = TRUE, 
                      timeToHosp = NULL, TPACutoff = 4.5, id){
  # Checking n
  if(!is.numeric(n) || length(n) != 1 || n <= 0 || n != as.integer(n)){
    stop("n must be a positive integer")
  }
  
  # Checking warn
  if(!is.logical(warn) || length(warn) != 1){
    stop("warn must be either TRUE or FALSE")
  }
  
  # Checking TPA
  if((!is.logical(TPA)) || !(length(TPA) %in% c(1, n))){
    stop("TPA must be a logical vector of length 1 or n")
  }
  if(length(TPA) == 1){
    TPA <- rep(TPA, n)
  }
  
  # Checking id
  if(missing(id)){
    id <- 1:n
  }else if(length(id) != n){
    stop("id must be missing or have length equal n")
  }
  firstid <- !(duplicated(id))
  firstidvals <- id[firstid]
  idindices <- match(id, firstidvals)
  
  # Checking timeToHosp
  if(!is.null(timeToHosp) && (length(timeToHosp) != n || !is.numeric(timeToHosp))){
    stop("timeToHosp must be a numeric vector of length n or NULL")
  }else if(any(timeToHosp < 0) || any(timeToHosp > 48)){
    stop("timeToHosp should be non-negative and in hours")
  }
  
  # Checking for valid values of regime
  regime <- tolower(regime)
  if(!all(regime %in% c("ds", "ms"))){
    stop("regime must only contain 'ds' or 'ms'")
  }else if(!(length(regime) %in% c(1, n))){
    stop("regime mus have a length of 1 or n")
  }
  if(length(regime) == 1){
    regime <- rep(regime, n)
  }
  dsInd <- which(regime == "ds")
  nds <- length(dsInd)
  msInd <- which(regime == "ms")
  nms <- length(msInd)
  
  # Checking times
  if(!is.list(times)){
    stop("times must be a list")
  }else if(!all(c("state", "RUCA") %in% names(times))){
    stop("times must have a state component and RUCA component")
  }
  
  # Checking RUCA
  if(length(RUCA) != 1){
    stop("RUCA must have a length of 1")
  }else if(!(RUCA %in% (1:4)) && !is.na(RUCA)){
    stop("RUCA must be a numeric value or character equal to 1, 2, 3, or 4")
  }else{
    RUCATimes <- times$RUCA[which(times$RUCA$RUCA == RUCA), , drop = FALSE]
  }
  
  # Checking TPACutoff
  if(length(TPACutoff) != 1 || !is.numeric(TPACutoff)){
    stop("TPACutoff must be a numeric value of length one")
  }else if(TPACutoff < 0){
    stop("TPACutoff must be a non-negative value")
  }
  
  # Checking hospRUCA
  if(!(length(hospRUCA) %in% c(1, n))){
    stop("hospRUCA must have a length of 1 or n")
  }else if(any(!(hospRUCA %in% (1:4)) & (!is.na(hospRUCA)) & (regime == "ds"))){
    stop("hospRUCA must be a numeric value or character equal to 1, 2, 3, or 4 for ds")
  }else if (length(hospRUCA) == 1){
    hospRUCA <- rep(hospRUCA, n)
  }
  
  # Checking states
  ## PSC
  ### This can be a vector
  if(!(length(state_PSC)) %in% c(1, n)){
    stop("state_PSC must be a vector of length 1 or n")
  }
  if(length(state_PSC) == 1){
    state_PSC <- rep(state_PSC, n)
  }
  
  ### Getting time quantiles
  PSCTimes <- lapply(unique(state_PSC), function(x){
    statequantiles <- checkState(x, times, warn)
    quantileind <- grep("^Quantile.+", colnames(statequantiles))
    quantiles <- as.matrix(statequantiles[, quantileind])
    quantiles <- quantiles[, -c(1, ncol(quantiles))]
    list(ind = (state_PSC == x),
         others = statequantiles[, -quantileind],
         quantiles = quantiles)
  })
  
  ## CSC
  ### This can be a vector
  if(!(length(state_CSC)) %in% c(1, n)){
    stop("state_CSC must be a vector of length 1 or n")
  }
  if(length(state_CSC) == 1){
    state_CSC <- rep(state_CSC, n)
  }
  
  ### Getting time quantiles
  CSCTimes <- lapply(unique(state_CSC), function(x){
    statequantiles <- checkState(x, times, warn)
    quantileind <- grep("^Quantile.+", colnames(statequantiles))
    quantiles <- as.matrix(statequantiles[, quantileind])
    quantiles <- quantiles[, -c(1, ncol(quantiles))]
    list(ind = (state_CSC == x),
         others = statequantiles[, -quantileind],
         quantiles = quantiles)
  })
  
  # Getting times
  ## DIDO and LKW
  RUCAquantileind <- grep("^Quantile.+", colnames(times$RUCA))
  quantiles <- as.matrix(RUCATimes[, RUCAquantileind])
  quantiles <- quantiles[, -c(1, ncol(quantiles))]
  
  ### Making sure LKW to EMS time is the same for identical ids
  LKWToEMSTime <- sample(quantiles[RUCATimes$Type == "LKW to EMS" & 
                                   RUCATimes$TimeType == "LKW to EMS", ], sum(firstid), replace = TRUE)
  LKWToEMSTime <- LKWToEMSTime[idindices]
  LKWToEMSTime2 <- LKWToEMSTime
  if(!is.null(timeToHosp)){
    LKWToEMSTime <- LKWToEMSTime + timeToHosp
  }
  
  ## Using TPA cutoff
  TPA <- TPA & (LKWToEMSTime <= TPACutoff)
  
  ### Making data.frame to store results
  results <- data.frame("LKWToEMS" = numeric(n), 
                        "TimeToHosp" = numeric(n), 
                        "TPA" = numeric(n), 
                        "DIDO" = numeric(n), 
                        "MT" = numeric(n))
  if(nds > 0){
    ## drip-and-ship
    ### Getting DIDO time for each hosp RUCA
    DIDOTime <- rep(0, n)
    for(i in unique(hospRUCA)){
      ### Not getting DIDO for hospRUCA encoded as less than one
      if(!is.na(i) && i >= 1){
        hospRUCAInd <- which(hospRUCA == i)
        hospRUCATimes <- times$RUCA[which(times$RUCA$RUCA == i), , drop = FALSE]
        DIDOTime[hospRUCAInd] <- sample(as.matrix(hospRUCATimes[hospRUCATimes$Type == "DIDO" & 
                                                      hospRUCATimes$TimeType == "DIDO", 
                                                      RUCAquantileind]), 
                                        length(hospRUCAInd), replace = TRUE)
      }
    }
    
    #### Making useful variables
    TPAds <- (TPA & (regime == "ds"))
    nTPAds <- sum(TPAds)
    noTPAds <- (!TPA & (regime == "ds"))
    nNoTPAds <- sum(noTPAds)
    if(nTPAds > 0){
      ### If they receive TPA
      TPATime <- rep(0, nTPAds)
      for(i in 1:length(PSCTimes)){
        tempPSCTimes <- PSCTimes[[i]]
        tempInd <- tempPSCTimes$ind[TPAds]
        ntemp <- sum(tempInd)
        if(ntemp > 0){
          TPATime[tempInd] <- sample(tempPSCTimes$quantiles[tempPSCTimes$others$Type == "PSC" & 
                                                            tempPSCTimes$others$TimeType == "TPA", ], 
                                     ntemp, replace = TRUE)
        }
      }
      
      MTTime <- rep(0, nTPAds)
      for(i in 1:length(CSCTimes)){
        tempCSCTimes <- CSCTimes[[i]]
        tempInd <- tempCSCTimes$ind[TPAds]
        ntemp <- sum(tempInd)
        if(ntemp > 0){
          MTTime[tempInd] <- sample(tempCSCTimes$quantiles[tempCSCTimes$others$Type == "Transfer" & 
                                                           tempCSCTimes$others$TimeType == "MT", ], 
                                    ntemp, replace = TRUE)
        }
      }
      results[TPAds, ] <- data.frame("LKWToEMS" = drop(LKWToEMSTime2[TPAds]), 
                                     "TimeToHosp" = drop(LKWToEMSTime[TPAds]), 
                                     "TPA" = drop(TPATime), 
                                     "DIDO" = drop(DIDOTime[TPAds]), 
                                     "MT" = drop(MTTime))
    }
    if(nNoTPAds > 0){
      ### If they don't receive TPA
      MTTime <- rep(0, nNoTPAds)
      for(i in 1:length(CSCTimes)){
        tempCSCTimes <- CSCTimes[[i]]
        tempInd <- tempCSCTimes$ind[noTPAds]
        ntemp <- sum(tempInd)
        if(ntemp > 0){
          MTTime[tempInd] <- sample(tempCSCTimes$quantiles[tempCSCTimes$others$Type == "Transfer" & 
                                                           tempCSCTimes$others$TimeType == "MT", ], 
                                    ntemp, replace = TRUE)
        }
      }
      results[noTPAds, ] <- data.frame("LKWToEMS" = drop(LKWToEMSTime2[noTPAds]), 
                                       "TimeToHosp" = drop(LKWToEMSTime[noTPAds]), 
                                       "TPA" = 0, 
                                       "DIDO" = drop(DIDOTime[noTPAds]), 
                                       "MT" = drop(MTTime))
    }
  }
  if(nms > 0){
    ## mothership
    #### Making useful variables
    TPAms <- (TPA & (regime == "ms"))
    nTPAms <- sum(TPAms)
    noTPAms <- (!TPA & (regime == "ms"))
    nNoTPAms <- sum(noTPAms)
    if(nTPAms > 0){
      TPATime <- rep(0, nTPAms)
      MTTime <- rep(0, nTPAms)
      for(i in 1:length(CSCTimes)){
        tempCSCTimes <- CSCTimes[[i]]
        tempInd <- tempCSCTimes$ind[TPAms]
        ntemp <- sum(tempInd)
        if(ntemp > 0){
          TPATime[tempInd] <- sample(tempCSCTimes$quantiles[tempCSCTimes$others$Type == "Mothership" & 
                                                            tempCSCTimes$others$TimeType == "TPA", ],
                                     ntemp, replace = TRUE)
          MTTime[tempInd] <- sample(tempCSCTimes$quantiles[tempCSCTimes$others$Type == "Mothership" & 
                                                           tempCSCTimes$others$TimeType == "MT(TPA)", ], 
                                    ntemp, replace = TRUE)
        }
      }
      
      results[TPAms, ] <- data.frame("LKWToEMS" = unname(LKWToEMSTime2[TPAms]), 
                                     "TimeToHosp" = unname(LKWToEMSTime[TPAms]), 
                                     "TPA" = unname(TPATime), 
                                     "DIDO" = 0, 
                                     "MT" = unname(MTTime))
    }
    if(nNoTPAms > 0){
      MTTime <- rep(0, nNoTPAms)
      for(i in 1:length(CSCTimes)){
        tempCSCTimes <- CSCTimes[[i]]
        tempInd <- tempCSCTimes$ind[noTPAms]
        ntemp <- sum(tempInd)
        if(ntemp > 0){
          MTTime[tempInd] <- sample(tempCSCTimes$quantiles[tempCSCTimes$others$Type == "Mothership" & 
                                                           tempCSCTimes$others$TimeType == "MT(No TPA)", ], 
                                    ntemp, replace = TRUE)
        }
      }
      
      results[noTPAms, ] <- data.frame("LKWToEMS" = unname(LKWToEMSTime2[noTPAms]), 
                                       "TimeToHosp" = unname(LKWToEMSTime[noTPAms]), 
                                       "TPA" = 0, 
                                       "DIDO" = 0, 
                                       "MT" = unname(MTTime))
    }
  } 
  
  ## Returning results
  return(results)
}

#' Generate simulated treatment times for each subject
#' @param TimeToHosp a numeric vector with the time from EMS attendance to hospital 
#' arrival for each subject in hours.
#' @param regime one of 'ds' or 'ms' to denote the desired treatment regime. Can also 
#' supply a vector with different values of 'ds' or 'ms' for each subject.
#' @param state a string with the name of the state or abbreviation. Can also supply a 
#' numeric value indicating the fips code. Can also supply a vector with a state 
#' for each subject.
#' @param RUCA a string or numeric value with the RUCA code for the stroke location.
#' @param hospRUCA a string or numeric value with the RUCA code for the first hospital 
#' in drip-and-ship.  
#' @param TransferTime a numeric vector with the transfer time for each subject in hours. 
#' This is only used for subjects with regime = "ds".
#' @param TPA a logical vector indicating whether subjects should receive TPA or not.
#' @param MT a logical vector indicating whether subjects should receive MT or not.
#' @param id a numeric vector denoting which stroke events should have the same 
#' time from LKW to EMS attention.
#' @param ... further arguments passed to [GetNTimes]
#' @return A data.frame with the following columns 
#' \item{`LKWToEMSTime`}{ the simulated time from last known well to EMS attention (hours)}
#' \item{`DIDOTime`}{ the simulated door-in-door-out (dido) time (hours)}
#' \item{`ArrivalToTPA`}{ the simulated time from hospital arrival to TPA (hours)}
#' \item{`ArrivalToMT`}{ the simulated time from hospital arrival to groin puncture (hours) if TPA is not administered at the same hospital}
#' \item{`TPAToMT`}{ the simulated time from TPA administration to groin puncture (hours) if TPA is administered at the same hospital}
#' \item{`TimeToHosp`}{ the simulated time from last known well to hospital arrival (hours)}
#' \item{`TPATime`}{ the simulated time from last known well to TPA (hours)}
#' \item{`MTTime`}{ the simulated time from last known well to groin puncture (hours)}
#' \item{`regime`}{ the regime used to create the simulated times}
#' \item{`state`}{ the state used to create the simulated times}
#' \item{`RUCA`}{ the RUCA code for stroke location used to create the simulated times}
#' \item{`hospRUCA`}{ the RUCA code for the first hospital used to create the simulated times}
GetTreatmentTimes <- function(TimeToHosp, TransferTime, regime, state_PSC, state_CSC, RUCA, hospRUCA, TPA, MT, id, ...){
  # Checking TimeToHosp
  if(!is.numeric(TimeToHosp) || length(TimeToHosp) == 0){
    stop("TimeToHosp must be a numeric vector")
  }else if(any(TimeToHosp < 0)){
    stop("TimeToHosp should be positive")
  }else if(any(TimeToHosp > 48)){
    stop("TimeToHosp should be less than 48 since it should be in hours")
  }
  
  # Checking TransferTime
  if(!is.numeric(TransferTime) || length(TransferTime) == 0){
    stop("TransferTime must be a numeric vector")
  }else if(length(TransferTime) != length(TimeToHosp)){
    stop("TransferTime must have the same length as TimeToHosp")
  }else if(any(TransferTime < 0)){
    stop("TransferTime should be positive")
  }else if(any(TransferTime > 48)){
    stop("TransferTime should be less than 48 since it should be in hours")
  }
  
  # Checking TPA
  if((!is.logical(TPA)) || !(length(TPA) %in% c(1, length(TimeToHosp)))){
    stop("TPA must be a logical vector of length 1 or the length of TimeToHosp")
  }
  if(length(TPA) == 1){
    TPA <- rep(TPA, length(TimeToHosp))
  }
  
  # Checking for valid values of regime
  regime <- tolower(regime)
  if(length(regime) == 1){
    regime <- rep(regime, length(TimeToHosp))
  }else if(length(regime) != length(TimeToHosp)){
    stop("length of regime must be either 1 or the length of TimeToHosp")
  }
  if(!is.character(regime) | !all(regime %in% c("ds", "ms"))){
    stop("regime must only conain 'ds' or 'ms'")
  }
  
  # Checking state_PSC and state_CSC
  ## state_PSC
  if(length(state_PSC) == 1){
    state_PSC <- rep(state_PSC, length(TimeToHosp))
  }else if(length(state_PSC) != length(TimeToHosp)){
    stop("length of state_PSC must be either 1 or the length of TimeToHosp")
  }
  
  ## state_CSC
  if(length(state_CSC) == 1){
    state_CSC <- rep(state_CSC, length(TimeToHosp))
  }else if(length(state_CSC) != length(TimeToHosp)){
    stop("length of state_CSC must be either 1 or the length of TimeToHosp")
  }
  
  # Checking hospRUCA
  if(length(hospRUCA) == 1){
    hospRUCA <- rep(hospRUCA, length(TimeToHosp))
  }else if(length(hospRUCA) != length(TimeToHosp)){
    stop("length of hospRUCA must be either 1 or the length of TimeToHosp")
  }
  
  # Checking RUCA
  if(length(RUCA) == 1){
    RUCA <- rep(RUCA, length(TimeToHosp))
  }else if(length(RUCA) != length(TimeToHosp)){
    stop("length of RUCA must be either 1 or the length of TimeToHosp")
  }
  
  # Checking id
  if(missing(id)){
    id <- 1:length(TimeToHosp)
  }else if(length(id) != length(TimeToHosp)){
    stop("id must be missing or have length equal to the length of TimeToHosp")
  }
  
  # Getting times
  ## Creating useful objects
  res <- data.frame("LKWToEMSTime" = numeric(), 
                    "DIDOTime" = numeric(), 
                    "ArrivalToTPA" = numeric(),
                    "ArrivalToMT" = numeric(),
                    "TPAToMT" = numeric(),
                    "TimeToHosp" = numeric(),
                    "TimeToCSC" = numeric(),
                    "TPATime" = numeric(), 
                    "MTTime" = numeric(), 
                    "regime" = numeric(), 
                    "state_PSC" = character(),
                    "state_CSC" = character(),
                    "RUCA" = do.call(class(RUCA), list()), 
                    "hospRUCA" = do.call(class(hospRUCA), list()))
  res[1:length(TimeToHosp), ] <- NA
  msbool <- (regime == "ms")
  dsbool <- (regime == "ds")
  
  ## Getting the rest of the times
  for(y in unique(RUCA)){
    ## Getting indices for this state and RUCA
    RUCAInd <- (!is.na(RUCA) & (RUCA == y))
    
    ### skipping this iteration if no one exists with this combination
    if(sum(RUCAInd) == 0 || is.na(y)){
      next
    }
      
    ## Getting times
    RUCATimes <- GetNTimes(sum(RUCAInd), 
                            regime = regime[RUCAInd], 
                            state_PSC = state_PSC[RUCAInd], state_CSC = state_CSC[RUCAInd], RUCA = y, 
                            hospRUCA = hospRUCA[RUCAInd], 
                            TPA = TPA[RUCAInd], 
                            id = id[RUCAInd], 
                            timeToHosp = TimeToHosp[RUCAInd], ...)
    
    ## Getting TPA from cutoff
    stateTPA <- (RUCATimes$TPA > 0)
    stateRegime <- regime[RUCAInd]
    
    ## Getting time to CSC
    TimeToCSC <- rep(NA, sum(RUCAInd))
    TimeToCSC[stateRegime == "ms"] <- RUCATimes$TimeToHosp[stateRegime == "ms"]
    TimeToCSC[stateRegime == "ds"] <- (RUCATimes$TimeToHosp + RUCATimes$DIDO + 
      RUCATimes$MT + TransferTime[RUCAInd])[stateRegime == "ds"]
    
    ## Adding results to res
    res[RUCAInd, ] <-  data.frame("LKWToEMSTime" = RUCATimes$LKWToEMS, 
                                   "DIDOTime" = RUCATimes$DIDO, 
                                   "ArrivalToTPA" = RUCATimes$TPA,
                                   "ArrivalToMT" = RUCATimes$MT * (stateRegime == "ds" | (!stateTPA)),
                                   "TPAToMT" = RUCATimes$MT * (stateRegime == "ms" & stateTPA),
                                   "TimeToHosp" = RUCATimes$TimeToHosp,
                                   "TimeToCSC" = TimeToCSC,
                                   "TPATime" = (RUCATimes$TimeToHosp + RUCATimes$TPA) * stateTPA, 
                                   "MTTime" = RUCATimes$TimeToHosp + RUCATimes$DIDO + 
                                     RUCATimes$MT + TransferTime[RUCAInd] + 
                                     (RUCATimes$TPA * (stateRegime == "ms")), 
                                   "regime" = stateRegime, "state_PSC" = state_PSC[RUCAInd],
                                   "state_CSC" = state_CSC[RUCAInd],
                                    "RUCA" = y, "hospRUCA" = hospRUCA[RUCAInd])
                                  
  }
  
  # Changing Time to groin puncture based on MT
  res$MTTime[!MT] <- 0
  
  ## Returning results
  return(res)
}
