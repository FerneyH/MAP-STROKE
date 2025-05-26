# Defining predict function for a Model object
## When using multilogitlink and type = "response", the 
## rows correspond to an individual.
## When using multilogitlink and type = "draws", the 
## 3rd dimension corresponds to individuals.
predict.Model <- function(object, newdata, type = "draws", ...){
  if(type != "draws"){
    stop("type must be 'draws'")
  }
  # Getting data
  mf <- model.frame(object$formula, data = newdata, xlev = object$xlev, 
                    na.action = "na.pass")
  x <- model.matrix(terms(mf), data = mf)
  
  # Getting predictions for each row of betas
  if(!is.matrix(object$draws) && is.array(object$draws)){
    ## Getting all linear predictors
    #res <- lapply(object$draws, function(temp){
    #  return(tcrossprod(temp, x))
    #})
    #res <- abind(res, along = 3)
    #res <- aperm(res, c(2, 3, 1))
    res <- arrayMultCpp(object$draws, x)
    res <- multinomialCubeCpp(res)
    dimnames(res)[[3]] <- object$names
    return(res)
  }else{
    ## Getting indices of columns that are not all zeros
    ### This isn't needed for multinomial case because it shouldn't have any 
    ### all zero columns
    ind <- which(colSums(abs(x), na.rm = TRUE) != 0)
    
    ## Getting all linear predictors and applying link function
    ### Only using variables that have any non-zero elements in x
    res <- tcrossprod(object$draws[, ind], x[, ind])
    if(identical(object$link, exp)){
      #res <- armaExpCpp(res)
      ## Testing in place exponentiation
      armaExpCpp3(res)
    }else{
      res <- object$link(res)
    }
    return(res)
  }
}

# Defining function to get predictions from TwoModels object
## This takes into account that some LVO patients will not receive EVT
predict.TwoModels <- function(models, newdata, EVTProb, type = "draws", 
                              predsDiagnostic = NULL, ...){
  
  if(type != "draws"){
    stop("type must be 'draws'")
  }
    
  # Getting predictions with each draw of coefficients
  # Each column contains the predictions with each set of coefficients for a single individual
  ## Getting predictions from first model
  if(is.null(predsDiagnostic)){
    predsDiagnostic <- predict(models$diagnosticModel, newdata = newdata, type = "draws")
  }
  
  ## Checking that predictions are between 0 and 1
  ### Not checking this for now because it is quite slow
  #if(any(!is.na(predsDiagnostic) & predsDiagnostic > 1) || 
  #   any(!is.na(predsDiagnostic) & predsDiagnostic < 0)){
  #  stop("Predictions from model 1 must be between 0 and 1")
  #}
  
  ## Getting predictions from second model
  preds <- lapply(1:dim(predsDiagnostic)[3], function(i){
    StrokeTypeName <- models$diagnosticModel$names[i]
    if(StrokeTypeName != "LVO"){
      tempData <- mutate(newdata, Type = StrokeTypeName) |>
        mutate(TPA = TPA * (Type == "LVO" | Type == "Not LVO"),
               MT = MT * (Type == "LVO"))
      predsOutcome <- predict(models$outcomeModel, newdata = tempData, type = "draws")
    }else{
      if(any(EVTProb > 0)){
        ### With EVT
        tempData <- mutate(newdata, Type = StrokeTypeName) |>
          mutate(TPA = TPA * (Type == "LVO" | Type == "Not LVO"),
                 MT = MT * (Type == "LVO"))
        predsOutcome <- predict(models$outcomeModel, newdata = tempData, type = "draws")
        #predsOutcome <- t(t(predsOutcome) * EVTProb)
        predsOutcome <- colwiseMultCpp(predsOutcome, EVTProb)
      }else{
        predsOutcome <- 0
      }
      
      if(any(EVTProb < 1)){
        ### Without EVT
        tempData <- mutate(newdata, Type = StrokeTypeName) |>
          mutate(TPA = TPA * (Type == "LVO" | Type == "Not LVO"),
                 MT = FALSE, MTTime = 0)
        predsOutcome <- colwiseMultCpp(predict(models$outcomeModel, newdata = tempData, 
                                    type = "draws"), (1 - EVTProb)) + predsOutcome
        #predsOutcome <- t(t(predict(models$outcomeModel, newdata = tempData, 
        #                            type = "draws")) * (1 - EVTProb)) + predsOutcome
      }
    }
    # Getting mean for truncated Poisson
    truncPoisMeanCpp3(predsOutcome, 6)
    predsDiagnostic[, , i] * predsOutcome
  })
  
  ## Combining predictions to get final predictions
  Reduce("+", preds)
}

# This function is used to generate replicates of mRS
## Currently uses the same nmRS randomly selected sets of draws for each stroke
simulatemRS <- function(model, data, nmRS){
  # simulate nmRS numbers
  sampleID <- sample(nrow(model$draws), nmRS, replace = TRUE)
  draws <- model$draws[sampleID, ]
  
  # Getting model matrix
  mf <- model.frame(model$formula, data = data, xlev = model$xlev, 
                    na.action = "na.pass")
  x <- model.matrix(terms(mf), data = mf)
  
  # Getting mean parameters
  lambdas <- tcrossprod(draws, x)
  armaExpCpp3(lambdas)
  
  # Generating samples
  mRS <- rtruncPoisCpp(lambdas, 6)
  
  # Combining samples with data
  mRS <- as.data.frame(t(mRS))
  colnames(mRS) <- paste0("SimulatedmRS_", 1:nmRS)
  return(mRS)
}
