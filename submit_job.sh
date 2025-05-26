#!/bin/bash

# Loading modules
module load stack/2022.2
module load r/4.2.2_gcc-9.5.0
module load r-abind/1.4-5_gcc-9.5.0
module load r-dplyr/1.0.10_gcc-9.5.0
module load r-rcpp/1.0.9_gcc-9.5.0
module load r-rcpparmadillo/0.11.4.0.1_gcc-9.5.0
module load r-readxl/1.4.1_gcc-9.5.0
module load r-sf/1.0-9_gcc-9.5.0
module load r-units/0.8-0_gcc-9.5.0

# "-pe smp 2": Single threaded program, just request 2 slot to get one full core
# "-cwd": Run this job in the current working directory
qsub -pe smp 8 -cwd -e /dev/null -o /dev/null -q UI,AIS -t 1-56 simulation_batch.job




