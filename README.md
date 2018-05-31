# calanda_range_shifts
Adaptations of the coexistence_range_shifts code for the data from the Calanda experiment

R code for simulating species range shifts and calculating the strength of coexistence between competitors. 
This code is based on the coexistence_range_shifts code, but includes several tweaks including: 

1. Statistical explorations of probability of survival and flowering, as a function of elevation and of plant size. 
2. Kriging based on the statistical fits. 
3. Dispersal kernels that are based on the WALD approach (see code for references) 
4. Analysis of persistence with dispersal kernels that vary in space (i.e. are per-site dispersal kernels). 

NOTE: calanda_pop_krigALL_prob_wald2.R is the cleaner version of the code, where all of the statistical exploration has been removed. The original file, calanda_pop_krigALL_prob_wald.R, is most useful for looking at the statistical analyses of demographic parameters.  


