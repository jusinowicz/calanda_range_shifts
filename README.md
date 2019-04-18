# calanda_range_shifts
Adaptations of the coexistence_range_shifts code for the data from the Calanda experiment

R code for simulating species range shifts and calculating the strength of coexistence between competitors. 
This code is based on the coexistence_range_shifts code from an earlier repository, but includes several changes/additions including: 

1. Statistical explorations of probability of survival and flowering, as a function of elevation, abiotic variables, and of plant size. 
2. Reproduction is a function of two processes that are both a function of spatial location: probability of flowering, total flowers 
3. GAMMs fit to abiotic covariates of elevation using R package gamm4, i.e. mean temp, min temp,GDD, and soil moisture. 
4. Kriging of intrinsic ranges based on the GAMM statistical fits.
5. Dispersal kernels that are based on the WALD approach (see code for references)
6. Environmental change modeled on an annual time-step with IPCC climate change scenarios 8.5 and 2.6. This code calculates population spread rates and checks whether species fail to track their intrinsic ranges. 

#Main code files: 
calanda_pop_ccs4_ciwald2.R Gives results for incrememntal environmental change based on IPCC scenarios for the average mean temperature increae per year, calculates population spread rates, and modifies species intrinsic ranges if they fail to track their ideal habitat conditions. 

calanda_pop_ccs4B_ciwald2.R Is similar to the previous file, except that it loads IPCC scenarios for the average mean temperature increae per year, calculates population spread rates, and modifies species intrinsic ranges if they fail to track their ideal habitat conditions. 

#The older main code files for statistical exploration of species intrinsic fitness and competition: 
calanda_pop_krigALL_prob_wald.R Is the original file that includes the statistical explorations of models (including tests for effects of plant size -- as in an Integral Projection Modeling (IPM) approach), simulations of population dynamics, and analysis of spatial coexistence mechanisms. 

calanda_pop_krigALL_prob_env7_wald2S.R Builds from the original analyses to try and uncover the main drivers of species intrinsic ranges. It includes a number of additional environmental variables (e.g. mean temp, max/min temp, soil moisture, moisture deficit, GDD, etc.) that have been built into the probabilty and number of flowers which determine the overall intrinsic range.

calanda_pop_krigALL_prob_env10_wald2.R Is the cleaner version of env7. Each species has a specific model of its intrinsic range that has been selected according to the following criteria: 1) Models are first selected that give a species a closed range, i.e. that reach 0 between the elevations of 0 and 3500m 2) Lowest AIC model is chosen from amongst these. 

#Code for figures:
calanda_figures5.R Code to produce plots of intrinsic ranges, environmental novelty, and the spatial LGR and its components (using data generated from . 

calanda_gcs_figures5.R Modified figure code to produce figures based on the output of the IPCC scenarios.

#Important functions called by the main code: 
range_coexistence_functionsWALD.R This is the heart of the population simulations and the coexistence analyses and includes the main functions for these portions of the analysis. This code corresponds to that in coexistence_range_shifts/range_coexistence_functions.R

wald_functions1.R Functions for calculating the WALD dispersal kernels. See Katul et al. 2005 for the background and analytical description of this approach. These are mechanistic dispersal kernels which account for turbulent updrafts and require detailed knowledge of the windspeed profiles at a site for calculation. 
