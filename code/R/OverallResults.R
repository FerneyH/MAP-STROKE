# Loading parallel package
library(parallel)
library(dplyr)
library(tidyr)
library(units)
library(rlang)



# Reading in times dataset to use for fips conversion to state names
timesState <- read.csv("datasets/SimulationTimesState.csv")
timesRUCA <- read.csv("datasets/SimulationTimesRUCA.csv")
times <- list("state" = timesState, "RUCA" = timesRUCA)


# Creating functions to be used in for loop
naMean <- function(x) mean(x, na.rm = TRUE)
Typefns <- list("LVO" = function(x) sum(x == "LVO"), 
                "Not LVO" = function(x) sum(x == "No LVO"),
                "Hemorrhagic" = function(x) sum(x == "Hemorrhagic"),
                "Stroke Mimic" = function(x) sum(x == "Stroke Mimic"))

clFunc <- function(ncores = 8,sim_idx){
  
  # Check if we have a full set of replicates
  have_results <- sapply(1:56, function(i){
    ## Converting FIPS to state name
    state <- unique(times$state$state_name[times$state$state_code == i])
    state <- toupper(state)
    
    ## Getting j as a string with padded 0s
    # this is just sim_idx
    j <- sim_idx
    j_str <- paste0(paste0(rep(0, 3 - nchar(j)), collapse = ""), j)
    
    ## Checking if this state has results
    cur_dir <- paste0("results/", state, "_", j_str)
    return(file.exists(paste0(cur_dir, "/data.rds")) && 
             file.exists(paste0(cur_dir, "/mRS_simul.rds")))
  })
  if (sum(have_results)!=49){
    return(NULL)
  }
  
  ## check if the results file is already there and we can skip the work
  j <- sim_idx
  j_str <- paste0(paste0(rep(0, 3 - nchar(j)), collapse = ""), j)
  
  if (file.exists(paste0("./mRS_percentages/table_rep_", j_str, ".rda"))){
    return(NULL)
  }
  
  
  ## Ok, we have data and no file yet - run the job:
  
  # Creating cluster
  cl <- makeCluster(ncores)
  clusterExport(cl, varlist = list("times", "naMean", "Typefns", "sim_idx"))
  parallel::clusterSetRNGStream(cl, 123123 + sim_idx )
  clusterEvalQ(cl, {
    ## Loading in packages
    library(dplyr)
    library(tidyr)
    library(rlang)
    library(units)}
  )
  
  # closing cl 
  on.exit(stopCluster(cl))
  
  # allres <- clusterApplyLB(cl, 1:200, function(j){
  res <- parLapply(cl, 1:56, function(i){
    ## Converting FIPS to state name
    state <- unique(times$state$state_name[times$state$state_code == i])
    state <- toupper(state)
    
    ## Getting j as a string with padded 0s
    # this is just sim_idx
    j <- sim_idx
    j_str <- paste0(paste0(rep(0, 3 - nchar(j)), collapse = ""), j)
    
    ## Checking if this state has results
    cur_dir <- paste0("results_big/", state, "_", j_str)
    
    if(!dir.exists(cur_dir)){
      return(NULL)
    }else{
      if(file.exists(paste0(cur_dir, "/data.rds")) && 
         file.exists(paste0(cur_dir, "/mRS_simul.rds"))){
        Preds <- readRDS(paste0(cur_dir, "/data.rds"))
        mRS_simul <- readRDS(paste0(cur_dir, "/mRS_simul.rds"))
      }
      else{
        return(NULL)
      }
    }
    
    # Making Preds
    samples <- sample(ncol(mRS_simul), nrow(mRS_simul), replace = TRUE)
    Preds <- cbind(Preds, "PredsKnown" = mRS_simul[cbind(1:nrow(mRS_simul), samples)])
    Preds$PredsKnownBin <- (Preds$PredsKnown <= 2)
    Preds$Cat0 <- (Preds$PredsKnown==0)
    Preds$Cat1 <- (Preds$PredsKnown==1)
    Preds$Cat2 <- (Preds$PredsKnown==2)
    Preds$Cat3 <- (Preds$PredsKnown==3)
    Preds$Cat4 <- (Preds$PredsKnown==4)
    Preds$Cat5 <- (Preds$PredsKnown==5)
    Preds$Cat6 <- (Preds$PredsKnown==6)
    
    
    
    ## Adding unique ID for each triage situation
    Preds <- group_by(Preds, longitude, latitude, regime) |> 
      mutate(triage_ID = 1:n()) |> 
      ungroup()
    
    # Getting summaries
    Preds <- drop_na(Preds, PredsKnown, PredsUnknown, TravelTime)
    Preds$TPA <- Preds$TPA * (Preds$Type %in% c("LVO", "Not LVO"))
    Preds$MT <- Preds$MT * (Preds$Type == "LVO")
    Preds$NoTreatment <- !(Preds$TPA | Preds$MT)
    Preds$TPATime <- Preds$TPATime * (Preds$Type %in% c("LVO", "Not LVO"))
    Preds$MTTime <- Preds$MTTime * (Preds$Type == "LVO")
    
    ## Generate random outcome for proportion <= 2
    ## Summaries to compare different triage algorithms
    ### Naive
    naivePreds <- group_by(Preds, latitude, longitude) |> 
      slice_min(TravelTime, n = 1, with_ties = FALSE) |> 
      ungroup() |> 
      mutate(triage = "Naive")
    
    ### Our decisions
    ourPreds <- group_by(Preds, latitude, longitude) |> 
      slice_min(PredsUnknown, n = 1, with_ties = FALSE) |> 
      ungroup() |> 
      mutate(triage = "MapStroke")
    
    ### Our decisions with ROPE
    ourROPEPreds <- mutate(Preds, triage = "MapStroke(ROPE)", 
                           selector = ifelse(regime == "ds", 
                                             PredsUnknown, 
                                             PredsUnknown - log(0.95))) |> 
      group_by(latitude, longitude) |> 
      slice_min(selector, n = 1, with_ties = FALSE) |> 
      ungroup()
    
    
    
    ### AHA decisions
    AHAPreds <- group_by(Preds, latitude, longitude, regime) |> 
      slice_min(TravelTime, n = 1, with_ties = FALSE) |> 
      ungroup() |> 
      select(latitude, longitude, regime, LVO, `Not LVO`, Hemorrhagic, `Stroke Mimic`, 
             state, Type, RUCA, rurality, age, RACE_scale, sysBP, diasBP, 
             TravelTime, TimeToHosp, TPATime, MTTime, PredsKnown, 
             PredsUnknown, PredsKnownBin, Cat0,Cat1,Cat2,Cat3,Cat4,Cat5,Cat6, TPA, MT, NoTreatment, triage_ID) |> 
      pivot_wider(names_from = regime, 
                  values_from = -c(LVO, `Not LVO`, Hemorrhagic, `Stroke Mimic`, 
                                   regime, state, Type,
                                   latitude, longitude, 
                                   RUCA, rurality, age, 
                                   RACE_scale, sysBP, diasBP)) |> 
      mutate(triage = "AHA", 
             timeDiff = TravelTime_ms - TravelTime_ds, 
             regime = case_when(
               is.na(TravelTime_ms) & !is.na(TravelTime_ds) ~ "ds", 
               is.na(TravelTime_ds) & !is.na(TravelTime_ms) ~ "ms",
               RACE_scale <= (4 - 6.347258) ~ "ds", 
               rurality == "Urban" & TravelTime_ms < 0.5 ~ "ms", 
               rurality == "Urban" & TravelTime_ms >= 0.5 ~ "ds", 
               rurality == "Suburban" & timeDiff < 0.5 & TravelTime_ms < 0.75 ~ "ms",
               rurality == "Suburban" & (timeDiff >= 0.5 | TravelTime_ms >= 0.75) ~ "ds", 
               rurality == "Rural" & timeDiff < 0.5 & TravelTime_ms < 1 ~ "ms",
               rurality == "Rural" & (timeDiff >= 0.5 | TravelTime_ms >= 1) ~ "ds", 
               TRUE ~ "other")) 
    
    AHAPreds$PredsKnown <- ifelse(AHAPreds$regime == "ms",
                                  AHAPreds$PredsKnown_ms, 
                                  AHAPreds$PredsKnown_ds)
    
    AHAPreds$PredsKnownBin <- ifelse(AHAPreds$regime == "ms",
                                     AHAPreds$PredsKnownBin_ms, 
                                     AHAPreds$PredsKnownBin_ds)
    
    AHAPreds$Cat0 <- ifelse(AHAPreds$regime == "ms",
                                     AHAPreds$Cat0_ms, 
                                     AHAPreds$Cat0_ds)
    
    AHAPreds$Cat1 <- ifelse(AHAPreds$regime == "ms",
                              AHAPreds$Cat1_ms, 
                              AHAPreds$Cat1_ds)
    
    AHAPreds$Cat2 <- ifelse(AHAPreds$regime == "ms",
                              AHAPreds$Cat2_ms, 
                              AHAPreds$Cat2_ds)
    
    AHAPreds$Cat3 <- ifelse(AHAPreds$regime == "ms",
                              AHAPreds$Cat3_ms, 
                              AHAPreds$Cat3_ds)
    
    AHAPreds$Cat4 <- ifelse(AHAPreds$regime == "ms",
                              AHAPreds$Cat4_ms, 
                              AHAPreds$Cat4_ds)
    
    AHAPreds$Cat5 <- ifelse(AHAPreds$regime == "ms",
                              AHAPreds$Cat5_ms, 
                              AHAPreds$Cat5_ds)
    
    AHAPreds$Cat6 <- ifelse(AHAPreds$regime == "ms",
                            AHAPreds$Cat6_ms, 
                            AHAPreds$Cat6_ds)
    
    AHAPreds$PredsUnknown <- ifelse(AHAPreds$regime == "ms",
                                    AHAPreds$PredsUnknown_ms, 
                                    AHAPreds$PredsUnknown_ds)
    AHAPreds$TimeToHosp <- ifelse(AHAPreds$regime == "ms",
                                  AHAPreds$TimeToHosp_ms, 
                                  AHAPreds$TimeToHosp_ds)
    AHAPreds$TPATime <- ifelse(AHAPreds$regime == "ms",
                               AHAPreds$TPATime_ms, 
                               AHAPreds$TPATime_ds)
    AHAPreds$MTTime <- ifelse(AHAPreds$regime == "ms",
                              AHAPreds$MTTime_ms, 
                              AHAPreds$MTTime_ds)
    AHAPreds$TPA <- ifelse(AHAPreds$regime == "ms",
                           AHAPreds$TPA_ms, 
                           AHAPreds$TPA_ds)
    AHAPreds$MT <- ifelse(AHAPreds$regime == "ms",
                          AHAPreds$MT_ms, 
                          AHAPreds$MT_ds)
    AHAPreds$NoTreatment <- ifelse(AHAPreds$regime == "ms",
                                   AHAPreds$NoTreatment_ms, 
                                   AHAPreds$NoTreatment_ds)
    AHAPreds$triage_ID <- ifelse(AHAPreds$regime == "ms",
                                 AHAPreds$triage_ID_ms, 
                                 AHAPreds$triage_ID_ds)
    
    ### Combining results
    combinedPreds <- bind_rows(ourPreds, ourROPEPreds, AHAPreds, naivePreds) |> 
      mutate(MS = (regime == "ms"))
    
    #### Getting concordance between each triage algorithm and AHA 
    concordanceAHA <- left_join(combinedPreds, y = AHAPreds, suffix = c("_our", "_AHA"),
                             by = c("longitude", "latitude")) |>
      mutate(concordance_regime = (regime_our == regime_AHA),
             concordance_hospital = (concordance_regime & 
                                       (triage_ID_our == triage_ID_AHA))) |> 
      select(longitude, latitude, triage_our, concordance_regime, concordance_hospital)
    
    #### Getting concordance between each triage algorithm and Naive routing
    concordanceNaive <- left_join(combinedPreds, y = naivePreds, suffix = c("_our", "_Naive"),
                             by = c("longitude", "latitude")) |>
      mutate(concordance_regime_Naive = (regime_our == regime_Naive),
             concordance_hospital_Naive = (concordance_regime_Naive & 
                                       (triage_ID_our == triage_ID_Naive))) |> 
      select(longitude, latitude, triage_our, concordance_regime_Naive, concordance_hospital_Naive)
    
    #### Combining concordance
    concordance<-left_join(concordanceAHA, concordanceNaive, 
                           by = c("latitude" = "latitude", "longitude" = "longitude", 
                                  "triage_our" = "triage_our"))
    
    #### Combining results with concordance and summarizing
    combinedPreds <- left_join(combinedPreds, concordance, 
                               by = c("latitude" = "latitude", "longitude" = "longitude", 
                                      "triage" = "triage_our")) |>    
      group_by(state, triage, rurality, Type) |> 
      summarize(Count = n(),
                across(c(TimeToHosp, TPATime, MTTime), 
                       .fns = function(x) mean(x[x > 0], na.rm = TRUE)), 
                across(c(PredsKnown, PredsUnknown, PredsKnownBin, Cat0,Cat1,Cat2,Cat3,Cat4,Cat5,Cat6, LVO, `Not LVO`, 
                         Hemorrhagic, `Stroke Mimic`, age, 
                         sysBP, diasBP, TPA, MT, NoTreatment, RACE_scale, 
                         concordance_regime, concordance_hospital, concordance_regime_Naive, concordance_hospital_Naive, 
                         MS), .fns = naMean), 
                #sum(TimeToHosp > 0), 
                TimeToHospCount = sum(Type %in% c("Hemorrhagic",
                                                  "LVO",
                                                  "Not LVO",
                                                  "Stroke Mimic")), # Reporting for all types
                TPACount = sum(Type %in% c("LVO", "Not LVO")),  # Changing to be TPA Eligible?
                TPACount2 = sum(Type %in% c("Not LVO")),  # Just among non-LVO
                MTCount = sum(Type %in% c("LVO")), # Same change - eligible to receive MT
                .groups = "keep") |> 
      ungroup()
    
    ### Saving combinedPreds and concordance in a list
    combinedPreds
  })
  ### Combine states here
  res <- bind_rows(res)
  if(nrow(res) > 0){
    res_all <- res |> group_by(triage, rurality, Type) |> 
      summarize(TimeToHosp = sum(TimeToHosp * TimeToHospCount, na.rm = TRUE) 
                / sum(TimeToHospCount, na.rm = TRUE), 
                MTTime = sum(MTTime * MTCount, na.rm = TRUE) / 
                  sum(MTCount, na.rm = TRUE),
                TPATime = sum(TPATime * TPACount, na.rm = TRUE) / 
                  sum(TPACount, na.rm = TRUE), 
                across(c(PredsKnown, PredsUnknown, PredsKnownBin, Cat0,Cat1,Cat2,Cat3,Cat4,Cat5,Cat6, LVO, `Not LVO`, 
                         Hemorrhagic, `Stroke Mimic`, age, 
                         sysBP, diasBP, TPA, MT, NoTreatment, RACE_scale, 
                         concordance_regime, concordance_hospital, concordance_regime_Naive, concordance_hospital_Naive, 
                         MS), .fns = function(x){sum(x * Count, na.rm = TRUE) / 
                             sum(Count, na.rm = TRUE)}), 
                Count = sum(Count),
                .groups = "keep")
    
    res_all_not_categories <- res |> group_by(triage) |> 
      summarize(TimeToHosp = sum(TimeToHosp * TimeToHospCount, na.rm = TRUE) 
                / sum(TimeToHospCount, na.rm = TRUE), 
                MTTime = sum(MTTime * MTCount, na.rm = TRUE) / 
                  sum(MTCount, na.rm = TRUE),
                TPATime = sum(TPATime * TPACount, na.rm = TRUE) / 
                  sum(TPACount, na.rm = TRUE),
                MT = sum(MT * MTCount, na.rm = TRUE) / 
                  sum(MTCount, na.rm = TRUE),
                TPA = sum(TPA * TPACount, na.rm = TRUE) / 
                  sum(TPACount, na.rm = TRUE),
                across(c(PredsKnown, PredsUnknown, PredsKnownBin,Cat0,Cat1,Cat2,Cat3,Cat4,Cat5,Cat6, LVO, `Not LVO`, 
                         Hemorrhagic, `Stroke Mimic`, age, 
                         sysBP, diasBP, NoTreatment, RACE_scale, 
                         concordance_regime, concordance_hospital, concordance_regime_Naive, concordance_hospital_Naive,
                         MS), .fns = function(x){sum(x * Count, na.rm = TRUE) / 
                             sum(Count, na.rm = TRUE)}), 
                Count = sum(Count),
                .groups = "keep")
    
    res_all_rurality <- res |> group_by(triage,rurality) |> 
      summarize(TimeToHosp = sum(TimeToHosp * TimeToHospCount, na.rm = TRUE) 
                / sum(TimeToHospCount, na.rm = TRUE), 
                MTTime = sum(MTTime * MTCount, na.rm = TRUE) / 
                  sum(MTCount, na.rm = TRUE),
                TPATime = sum(TPATime * TPACount, na.rm = TRUE) / 
                  sum(TPACount, na.rm = TRUE),
                MT = sum(MT * MTCount, na.rm = TRUE) / 
                  sum(MTCount, na.rm = TRUE),
                TPA = sum(TPA * TPACount, na.rm = TRUE) / 
                  sum(TPACount, na.rm = TRUE),
                across(c(PredsKnown, PredsUnknown, PredsKnownBin,Cat0,Cat1,Cat2,Cat3,Cat4,Cat5,Cat6, LVO, `Not LVO`, 
                         Hemorrhagic, `Stroke Mimic`, age, 
                         sysBP, diasBP, NoTreatment, RACE_scale, 
                         concordance_regime, concordance_hospital, concordance_regime_Naive, concordance_hospital_Naive,
                         MS), .fns = function(x){sum(x * Count, na.rm = TRUE) / 
                             sum(Count, na.rm = TRUE)}), 
                Count = sum(Count),
                .groups = "keep")
    
    res_all_type <- res |> group_by(triage,Type) |> 
      summarize(TimeToHosp = sum(TimeToHosp * TimeToHospCount, na.rm = TRUE) 
                / sum(TimeToHospCount, na.rm = TRUE), 
                MTTime = sum(MTTime * MTCount, na.rm = TRUE) / 
                  sum(MTCount, na.rm = TRUE),
                TPATime = sum(TPATime * TPACount, na.rm = TRUE) / 
                  sum(TPACount, na.rm = TRUE),
                MT = sum(MT * MTCount, na.rm = TRUE) / 
                  sum(MTCount, na.rm = TRUE),
                TPA = sum(TPA * TPACount, na.rm = TRUE) / 
                  sum(TPACount, na.rm = TRUE),
                across(c(PredsKnown, PredsUnknown, PredsKnownBin,Cat0,Cat1,Cat2,Cat3,Cat4,Cat5,Cat6, LVO, `Not LVO`, 
                         Hemorrhagic, `Stroke Mimic`, age, 
                         sysBP, diasBP, NoTreatment, RACE_scale, 
                         concordance_regime, concordance_hospital, concordance_regime_Naive, concordance_hospital_Naive,
                         MS), .fns = function(x){sum(x * Count, na.rm = TRUE) / 
                             sum(Count, na.rm = TRUE)}), 
                Count = sum(Count),
                .groups = "keep")
    
    # Save result to a file:
    if (!dir.exists("./big_tables_type/")){
      dir.create("./big_tables_type/")
    }
    # Save result to a file:
    if (!dir.exists("./big_tables_rurality/")){
      dir.create("./big_tables_rurality/")
    }
    # Save result to a file:
    if (!dir.exists("./big_tables_not_categories/")){
      dir.create("./big_tables_not_categories/")
    }
    # Save result to a file:
    if (!dir.exists("./big_tables_type_rurality/")){
      dir.create("./big_tables_type_rurality/")
    }
    
    ### Removing this for now 
    # if (!dir.exists("./big_tables_state_type_rurality/")){
    #  dir.create("./big_tables_state_type_rurality/")
    # }
    
    j <- sim_idx
    j_str <- paste0(paste0(rep(0, 3 - nchar(j)), collapse = ""), j)
    save(res_all, file = paste0("./big_tables_type_rurality/table_rep_", j_str, ".rda"))
    save(res_all_not_categories, file = paste0("./big_tables_not_categories/table_rep_", j_str, ".rda"))
    save(res_all_rurality, file = paste0("./big_tables_rurality/table_rep_", j_str, ".rda"))
    save(res_all_type, file = paste0("./big_tables_type/table_rep_", j_str, ".rda"))
    return(res)
  }else{
    stop("Should have results here!")
  }
}

args = commandArgs(trailingOnly=TRUE)
if(length(args) == 0){
  stop("no fips code found")
} else {
  sim_idx <- as.numeric(args[[1]])
}

# Running function
clFunc(8,sim_idx)
