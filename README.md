<p align="center" >
  <img src="img/Map-Stroke_Logo.png" alt="logo" width="200" />
</p>

# MAP-STROKE

We developed a novel personalized prehospital destination selection algorithm using jointly fit Bayesian models, trained on over 10,000 patient records from clinical trials and real-world data. The algorithm incorporates transport times, locations, and hospital capabilities to optimize prehospital triage decisions.

The repository contains the main functions and simulated datasets used in the MAP-STROKE project.

1. _**Data**_: Contains census tract datasets for each state and summary data from GWTG.
2. _**Generated**_: Contains simulated stroke locations and driving times to the EVT-capable stroke centers and Local Hospitals.
3. _**Code**_: Contains the functions used:
   
   -`Multinomial.Cpp`: Contains Rcpp functions for performing tasks related to multinomial regression;

   -`OtherFunctions.Cpp`: Contains Rcpp functions that much more efficiently perform tasks that do not fit in the other files;

   -`TruncatedPoisson.Cpp`: Contains Rcpp functions for performing tasks related to  truncated Poisson distributions;

   -`OtherFunctions.R`: Contains code for getting RUCA codes and generating probability of receiving EVT;

   -`PredictionFunctions.R`: Contains code for getting predictions from our model;

   -`SimulatePatientProfiles.R`: Contains code for simulating patient characteristics from GWTG, generating patient profiles, and getting predictions;
  
   -`SimulateTimeFunctions.R`: Contains code for generating treatment times for patients;

   -`DrivingTimestoCSC.R`: Contains code for generating driving times to EVT-capable stroke centers;

   -`DrivingTimestoHospital.R`: Contains code for generating driving times to local hospitals;

   -`RandomWeightedStrokeEvents.R`: Contains code for generating random stroke event locations based on population density and county-level rates from the CDC.
   
    

