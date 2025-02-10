The conformal prediction R code needs the following input:

1. MSI chl-a estimates (see .csv file for formatting (CV(not essential), OWT (not essential), Algorithm (not essential), Date, MSI_chl, AC (not essential, Lake (not essential))
2. Observed chl-a values (see .csv file for formatting (Date, Chlorophyll))
     - Used for (3-day) Matchup with MSI chl-a estimates
     - Used to determine estimates for MSI chl-a estimates without available matchups
  
The output is a csv file with iteration year, Date, match up, MSI chl-a, conformal (yes/no), prediction intervals (Upper and Lower), 
          non-conformity scores (uncertainty) and alpha_test

More detail about the code can be found in ... paper 
