Here is a description of all of the files in this directory
------------------------------------------------------------

datasets - Contains census tract datasets for each state, summary data from GWTG, and census tract dataset.

functions
----------
Cpp
----------
Multinomial.Cpp - Contains Rcpp functions that much more efficiently perform tasks for multinomial regression.

OtherFunctions.Cpp - Contains Rcpp functions that much more efficiently perform tasks that don't fit in the other files.

TruncatedPoisson.Cpp - Contains Rcpp functions that much more efficiently perform tasks for truncated Poisson distributions.

R
----------

OtherFunctions.R - Contains code for getting RUCA codes and generating probability of receiving EVT.

PredictionFunctions.R - Contains code for getting predictions from our model.

SimulatePatientProfiles.R - Contains code for simulating patient characteristics from GWTG. Also has a 
function that is used to generate patient profiles and get predictions.

SimulateTimeFunctions.R - Contains code for generating treatment times for patients.

----------
generated - Contains simulated stroke location and travel time data from Ferney.

results - Contains results from the simulations.

CreatingTractDatasets.R - Creates census tract datasets to help get RUCA codes.

simulation_batch.job - job file for HPC

simulation_batch.R - R file used to do simulations on the HPC

submit_job.sh - Used to submit array job to the HPC






