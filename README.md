# calanda_range_shifts
Adaptations of the coexistence_range_shifts code for the data from the Calanda experiment

R code for simulating species range shifts and calculating the strength of coexistence between competitors. 
This code is based on the coexistence_range_shifts code from an earlier repository, but includes several changes/additions including: 

1. Statistical explorations of probability of survival and flowering, as a function of elevation, abiotic variables, and of plant size. 
2. Kriging based on the statistical fits. 
3. Dispersal kernels that are based on the WALD approach (see code for references) 
4. Analysis of persistence with dispersal kernels that vary in space (i.e. are per-site dispersal kernels). 

calanda_pop_krigALL_prob_wald.R Is the original file that includes the statistical explorations of models (including tests for effects of plant size -- as in an Integral Projection Modeling (IPM) approach), simulations of population dynamics, and analysis of spatial coexistence mechanisms. 

calanda_pop_krigALL_prob_env6_wald2.R Builds from the original analyses to try and uncover the main drivers of species intrinsic ranges. It includes a number of additional environmental variables (e.g. mean temp, max/min temp, soil moisture, moisture deficit, GDD, etc.) that have been built into the probabilty and number of flowers which determine the overall intrinsic range.

calanda_pop_krigALL_prob_env7_wald2S.R Is the cleaner version of the code, where all of the statistical exploration has been removed. 

range_coexistence_functionsWALD.R This is the heart of the population simulations and the coexistence analyses and includes the main functions for these portions of the analysis. This code corresponds to that in coexistence_range_shifts/range_coexistence_functions.R

wald_functions1.R Functions for calculating the WALD dispersal kernels. See Katul et al. 2005 for the background and analytical description of this approach. These are mechanistic dispersal kernels which account for turbulent updrafts and require detailed knowledge of the windspeed profiles at a site for calculation. 
