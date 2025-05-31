<p align="center" >
  <img src="img/Map-Stroke_Logo.png" alt="logo" width="200" />
</p>

# MAP-STROKE

We developed a novel personalized prehospital destination selection algorithm using jointly fit Bayesian models, trained on over 10,000 patient records from clinical trials and real-world data. The algorithm incorporates transport times, locations, and hospital capabilities to optimize prehospital triage decisions.

The repository has the main functions and simulated data sets used in the MAP-STROKE project.

1. _**Data**_: Contains census tract datasets for each state and summary data from GWTG.
2. _**Generated**_: Contains simulated stroke locations and driving times to the EVT-capable stroke centers and Local Hospitals.
3. _**Code**_: Contains the functions used:
-`Multinomial.Cpp`: Contains Rcpp functions that much more efficiently perform tasks for multinomial regression;
-`OtherFunctions.Cpp`: Contains Rcpp functions that much more efficiently perform tasks that don't fit in the other files;
-`TruncatedPoisson.Cpp`: Contains Rcpp functions that much more efficiently perform tasks for truncated Poisson distributions;
-`therFunctions.R` - Contains code for getting RUCA codes and generating probability of receiving EVT;
-`PredictionFunctions.R` - Contains code for getting predictions from our model;
-`SimulatePatientProfiles.R` - Contains code for simulating patient characteristics from GWTG, generate patient profiles, and get predictions;
-`SimulateTimeFunctions.R` - Contains code for generating treatment times for patients.
    

