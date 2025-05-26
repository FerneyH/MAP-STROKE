# Loading time datasets
timesState <- read.csv("datasets/SimulationTimesState.csv")
timesRUCA <- read.csv("datasets/SimulationTimesRUCA.csv")
times <- list("state" = timesState, "RUCA" = timesRUCA)

# Creating and saving datasets
for(i in 1:56){
  fips <- i
  if(nchar(i) == 1){
    fips <- paste0(0, fips)
  }
  state <- unique(times$state$state_name[times$state$state_code == i])
  state <- toupper(state)
  if(!file.exists(paste0("generated/", state, "_Random_Strokes.rda"))){
    #stop(paste0("fips code ", state, " not supported"))
    next
  }
  tracts <- tigris::tracts(state = state, progress_bar = FALSE, 
                           year = 2010)
  save(tracts, file = paste0("datasets/", fips, "_tracts.rda"))
  save(tracts, file = paste0("datasets/", state, "_tracts.rda"))
}
