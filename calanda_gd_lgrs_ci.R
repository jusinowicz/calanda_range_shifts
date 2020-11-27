#=============================================================================
# R code for simulating climate-driven range shifts and calculating the 
# strength of coexistence between competitors. This code is built around a 
# specific spatially excplicit model of population dynamics, the Leslie-Gower 
# model (with some slight modifications). It includes a number of statistical 
# models for inferring Leslie-Gower parameters from experimental data. 
#
# The main loop of this code increments over yearly changes in environmental
# conditions taken from one of two IPCC scenarios that bookend expectations 
# (2.6 and 8.5). 
# 
# 1. Underlying environmental data come from a combination of sources: HOBO
#  data loggers placed in the field from 2015 -2017, and TerraClimate maps.
#  See calanda_env2017B_2.R for the pre-processing approach to these data. 
#  Here, the data are assumed to be available already in the main data
#  variable (allB2). These data are loaded and fit by simple GAMMs (with 
#  year included as a random effect) with elevation as the only covariate,
#  which are used to interpolate the environmental variables at a finer 
#  spatial scale.
# 
# 2. Survival is fit by a GAMM (year as a random effect) to look for significant
#  variation over elevation. Since it is fairly minimal over the sampled 
#  elevation range (1000-2000m) an average value is used in the Leslie-Gower
#    model. 
# 
# 3. Competition coefficients are fit through a combination of bootstrapping 
#  and NLS. See calanda_stats_tests1.R for an exploration of different 
#  approaches for finding the competition coefficients. 
#
# 4. Dispersal kernels are based on the WALD approach (Katul, G. G., et al. 2005. 
#    Mechanistic Analytical Models for Long‐Distance Seed Dispersal by Wind. 
#  The American Naturalist 166:368–381). This requires calculating the height
#  of plants from field measurements. Average height is used across sampled
#  elevations. A number of other parameters are determined either from the 
#  primary literature or field measurements of wind. See wald_functions1.R for 
#    more information and references. 
#
# 5. Intrinsic ranges are calculated as a product of two separate underlying 
#    processes: probability of flowering, and number of flowers. Each process
#  is fit using a GAMM (year as a random effect). Intrinsic ranges can be 
#  fit with a number of environmental covariates, and this code allows a number
#  of combinations. However, the covariates of primary interest are the 
#    mean growing season temp, min temp, number of growing degree days (GDD), 
#  and monthly average soil moisture. In order to produce nice ranges that
#  taper towards zero at range margins, the underlying smooth fits of covariates
#  have been constrained to go to 0 at upper and lower intervals. This has 
#  necessitated writing a bit of additional code by hand. 
#    
#  Range shifts are modeled by using the GAMMs (which have been fit to baseline 
#    conditions) to project a new range given the underlying change in covariates 
#    with each IPCC scenario. 
#
# 6. Coexistence is calculated and the coexistence mehanisms are parsed out in 
#  two phases. Resident stationary distributions are calculate, and the invader
#    low-density stationary distribution is calculated. A combination of analytical
#  approximations and simulations are used to calculate the LDG and the 
#  fitness-density covariance. See range_coexistence_functionsWALD.R. 
#
# 7. The full community stationary distribution is calculted numerically. 
#=============================================================================
#=============================================================================
# Load these libraries
#=============================================================================
library(mgcv)
library(gamm4)
library(lme4)
library(plyr)
library(MASS)
library(fields)
library(viridis)
source("./wald_functions2.R")
source("./range_coexistence_functionsWALD.R")
#source("./range_coexistence_functionsD.R") #The old library

#=============================================================================
# Do this for lme4
#=============================================================================
options(lmerControl=list(check.nobs.vs.rankZ = "warning",
  check.nobs.vs.nlev = "warning",
  check.nobs.vs.nRE = "warning",
  check.nlev.gtreq.5 = "warning",
  check.nlev.gtr.1 = "warning"))
options(glmerControl=list(check.nobs.vs.rankZ = "warning",
  check.nobs.vs.nlev = "warning",
  check.nobs.vs.nRE = "warning",
  check.nlev.gtreq.5 = "warning",
  check.nlev.gtr.1 = "warning"))
#=============================================================================
#For naming files
#=============================================================================
#f.name1=c("calanda_ccs26gd_temp2_2017")
f.name1=c("calanda_ccs26gd_tsmfull4_2017")
#f.name1=c("calanda_ccs85gd_temp2_2017")
#f.name1=c("calanda_ccs85gd_tsmfull3_2017")
#=============================================================================
#Data: (see calanda2017.R)
#=============================================================================
load(file="calanda_allB2017_2.var")
rm(allB)

#With climate-change scenarios:
#ccs_temp = read.csv("rcp85_2019_2100.csv")[2:3]
ccs_temp = read.csv("rcp26_2019_2100.csv")[2:3]

#ccs_soil = read.csv("rcp85_mrso_2019_2100.csv")[1:2]
ccs_soil = read.csv("rcp26_mrso_2019_2100.csv")[1:2]

#=============================================================================
#Some necessary tweaks to the data set to fix various things 
#=============================================================================
#Fix gs_mean_temp for allB2 in year 2015 by fitting a model to the other years
#by elevation, then projecting based on the mean temp at 1400m (the only reliable
#data point from 2015):

elevations = unique(allB2$elevation)[2:6]
gs_tmp = subset(allB2, year != 2015)
gs_tmp2015 = subset(allB2, year == 2015)
gs_tmp2 = subset(allB2, year == 2015 & elevation == 1400 )

#Fit a model and predict the new temperatures. 
lm_gs = gamm4(gs_mean_temp~s(elevation,k=4), random = ~(1|year),data=gs_tmp)
pgs =  predict.gam(lm_gs$gam, newdata=data.frame(elevation=elevations),type="response")

#This will serve as a lookup table for the new values of gs_mean_temp
pgs2015 = data.frame(elevation = elevations, gs_mean_temp =pgs - (pgs[2]- unique(gs_tmp2$gs_mean_temp)))

#Make a new object with a column of gs_mean_temp to replace the one in allB2
temp.all=pgs2015[match(interaction(gs_tmp2015$elevation,"2015"), interaction(pgs2015$elevation,"2015") ),]
allB2[1:1200,]$gs_mean_temp = temp.all$gs_mean_temp

#Fix max_gdd for allB2 in year 2015 by fitting a model to the other years
#by elevation, then projecting based on the mean temp at 1400m (the only reliable
#data point from 2015):

elevations = unique(allB2$elevation)[2:6]
gdd_tmp = subset(allB2, year != 2015)
gdd_tmp2015 = subset(allB2, year == 2015)
gdd_tmp2 = subset(allB2, year == 2015 & elevation == 1400 )

#Fit a model and predict the new temperatures. 
lm_gdd = gamm4(max_gdd~s(elevation,k=4), random = ~(1|year),data=gdd_tmp)
pgdd =  predict.gam(lm_gdd$gam, newdata=data.frame(elevation=elevations),type="response")

#This will serve as a lookup table for the new values of gs_mean_temp
pgdd2015 = data.frame(elevation = elevations, max_gdd =pgdd - (pgdd[2]- unique(gdd_tmp2$max_gdd)))

#Make a new object with a column of gs_mean_temp to replace the one in allB2
temp.all=pgdd2015[match(interaction(gdd_tmp2015$elevation,"2015"), interaction(pgdd2015$elevation,"2015") ),]
allB2[1:1200,]$max_gdd = temp.all$max_gdd

#Modify allB2 so that NA max_gdd counts match extremes: 
allB2[allB2$elevation==600,]$max_gdd = subset(allB2,allB2$elevation==1000)$max_gdd
allB2[allB2$elevation==2100,]$max_gdd = subset(allB2,allB2$elevation==2000)$max_gdd
allB2[allB2$elevation==2700,]$max_gdd = subset(allB2,allB2$elevation==2000)$max_gdd

#Modify allB2 so that NA gs_min_temp counts match extremes: 
gmt_max1400 = max(subset(allB2, elevation == 1400 )$gs_min_temp, na.rm=T)

allB2[allB2$elevation==1400,]$gs_min_temp = gmt_max1400 
allB2[allB2$elevation==600,]$gs_min_temp = subset(allB2,allB2$elevation==1000)$gs_min_temp
allB2[allB2$elevation==2100,]$gs_min_temp = subset(allB2,allB2$elevation==2000)$gs_min_temp
allB2[allB2$elevation==2700,]$gs_min_temp = subset(allB2,allB2$elevation==2000)$gs_min_temp

#Modify allB2 so that NA inflor counts are 0: 
allB2[allB2$elevation==600,][2:4] = subset(allB2,allB2$elevation==1000)[,2:4] 
allB2[allB2$elevation==600,]$nr.inflor[is.na(allB2[allB2$elevation==600,]$nr.inflor)] = 0

allB2[allB2$elevation==2100,][2:4] = subset(allB2,allB2$elevation==1000)[,2:4] 
allB2[allB2$elevation==2100,]$nr.inflor[is.na(allB2[allB2$elevation==2100,]$nr.inflor)] = 0

allB2[allB2$elevation==2700,][2:4] = subset(allB2,allB2$elevation==1000)[,2:4] 
allB2[allB2$elevation==2700,]$nr.inflor[is.na(allB2[allB2$elevation==2700,]$nr.inflor)] = 0

allB2$flyes[is.na(allB2$flyes)] = 0
head(allB2)

#Dactylis glomerata conversion factor for flowers to seeds: 
lfbft= c(9.5,11.4,17.5,16.6,6,11.8,14.5,11.5,10.5,6.6)
lsbft=c(6.5,11.5,10.5,9.5,4.9,7.2,8.6,6.7,7.7,5.3)
wls = c(0.8,1,1.2,1.8,0.7,1.6,1.3,0.9,0.9,0.7)
lfb = c(3.5,4,5.5,6,1,4.1,5.5,4.5,2.3,1.2)

lfbft_lg = log(lfbft)
lsbft_lg = log(lsbft)
wls_lg = log(wls)
lfb_lg = log(lfb)

#In order of best explanatory fit to worst. All are significant 
lm3 = lm(lfbft_lg~lfb_lg) #R2 of 0.8792
lm2 = lm(lfbft_lg~wls_lg) #R2 of 0.6367
lm4 = lm(lsbft_lg~lfb_lg) #R2 of 0.5743
lm1 = lm(lsbft_lg~wls_lg) #R2 of 0.3307

#An average of 58.3 per cm. See Tormo-Molina et al. 2015
# dg_conv = 58
# allB2[ (rownames(subset(allB2,Sp =="DG") )), ]$nr.inflor = dg_conv*
#   	exp(log(allB2[(row.names(subset(allB2,Sp =="DG") )), ]$nr.inflor) * lm2$coefficients[2] +
#   	lm2$coefficients[1]) 

allB2$elevation =c(allB2$elevation)
allB2$e2 = c(allB2$e2)

allB2=data.frame(allB2, year_numeric = as.numeric(allB2$year))

nspp = length(unique(allB2$Sp))

#=============================================================================
#Which version of the intrinsic ranges? Which abiotic factors to use? 
#=============================================================================
# All of the possible variables for the model. Right now, this is 
# elevation, year, max_gdd, gs_min_temp,gs_max_temp, gs_mean_temp, soil_moist
variable_mat = colnames(allB2)[c(14,17,30,32,33,34,37)]
v_use_prob = vector("list", nspp)
v_use_flow = vector("list", nspp)

#Define species-specific variable combinations for both probability
#of flowering and number of flowers
v_use_prob [[1]] = variable_mat[c(2,4,6,7)]
v_use_prob [[2]] = variable_mat[c(2,3,6,7)]
v_use_prob [[3]] = variable_mat[c(2,3,4,6,7)]

v_use_flow [[1]] = variable_mat[c(2,3,6,7)]
v_use_flow [[2]] = variable_mat[c(2,3,6,7)]
v_use_flow [[3]] = variable_mat[c(2,3,4,6,7)]

nvar_probs = matrix(0,nspp,1)
for (n in 1:nspp) {nvar_probs[n] = length(v_use_prob [[n]])-1 }
nvar_flows = matrix(0,nspp,1)
for (n in 1:nspp) {nvar_flows[n] = length(v_use_flow [[n]])-1 }

#Speces-specific, variable-specific knots for smooth fits.
#These are based on analysis of underlying factor fits (see file: )
#Each species has 2 columns, one for flower probability model and one for total
#flowers model:

dgk = matrix(3,5,2); axk = dgk; hnk = dgk
dgk[2,1] = 2 #Soil moisture for DG prob
#axk[2,1] = 5 #Soil moisture for AX prob
#hnk[2:3,1] = 5 #Soil moisture and min temp for HN prob
var_spp_knots = list(dgk,axk,hnk)

#=============================================================================
#Tunable lattice parameters (space and time)
#=============================================================================
#If fast.bys = TRUE, then populations will be allowed to spread and track their
#intrinsic ranges. 
fast.bys=FALSE
c.tol = 1e-2 #Tolerance for finite inidividuals, sets left/right range edges

#For the climate change scenarios, tune the change in temp to the current 
#average from field data: 
temp_field = mean(allB2$gs_mean_temp,na.rm=T)
ccs_temp = ccs_temp[,2]- temp_field-(ccs_temp[1,2]- temp_field)

ccs_soil=ccs_soil/100
soil_field = mean(allB2$soil_moist,na.rm=T)
ccs_soil = ccs_soil[,2]- soil_field - (ccs_soil[1,2]-soil_field)


#For standard increments of time through the ccs: 
tstart = 1
tstop = 61
tinc = 10
ttot= ceiling((tstop-tstart)/tinc)

ngens=ttot #Generations for environmental change
iconfig=1 #Interattions of initial configuration of species
n_draws = 10000 #Draws from posterior distribution of flowering
n_runs = 1 #Draws from PD of flowering for LGR prediction intervals

#Spatial scale: Assume that the gradient is from 500m to 2500m in units of 10m
elevations=c(600,1000,1400,1600,1800,2000,2100,2700)
el_real=c(1000,1400,1600,1800,2000)
el1=0
el2=3500

xx1= (seq(el1,el2,1))
xx2= (xx1^2)

#How to handle years?
years= unique((allB2$year))
#years= unique(as.character(allB$year))
#levels(xx$year) = years

yx= matrix(2017,length(xx1),1)
#yx= c(matrix(3,length(xx1),1))
xx=data.frame( elevation=xx1, year=years[1], e2=xx2)

#=============================================================================
#Environmental variation -- defining the abiotic influence on intrinsic growth
#rates and projecting it into the future. 
#=============================================================================

#=============================================================================
#Fit a GAM with each of the environmental variables against elevation
#Note: B2 is the subset of B without the first year of data, which
#is not trustworthy for the hobos. 

#Temp
temp_gam = gamm4(gs_mean_temp~s(elevation,k=4), random = ~(1|year),data=subset(allB2, !is.na(allB2$gs_mean_temp)))

#soil moisture 
sm_gam = gamm4(soil_moist~s(elevation,k=4), random = ~(1|year),data=subset(allB2, !is.na(allB2$gs_mean_temp)))


#PET
sp_gam = gamm4( soil_pet~ elevation+I(elevation^2),random = ~(1|year), data=allB2)


#Moisture deficit
sd_gam = gamm4( soil_def~ elevation+I(elevation^2),random = ~(1|year), data=allB2)

#annual mean light
ml_gam = gamm4( an_mean_light~ I(elevation^2), random = ~(1|year), data=subset(allB2, !is.na(allB2$gs_mean_temp)) )

#GDD
gdd_gam = gamm4(max_gdd~s(elevation,k=5), random = ~(1|year),data=allB2)

#Minimum growing season temp
gmt_gam = gamm4( gs_min_temp~ elevation+ I(elevation^2), random = ~(1|year), data=subset(allB2, !is.na(allB2$gs_mean_temp)) )

#=============================================================================
# Krig:
xx_mt = as.vector(predict.gam(temp_gam$gam, newdata=xx,type= "response"))

#The soil moisture vector: 
xx_sm = as.vector(predict.gam(sm_gam$gam, newdata=xx,type= "response"))

#The soil moisture vector: 
xx_sp = as.vector(predict.gam(sp_gam$gam, newdata=xx,type= "response"))

#The soil moisture vector: 
xx_sd = as.vector(predict.gam(sd_gam$gam, newdata=xx,type= "response"))

#Mean light
xx_ml = as.vector(predict.gam(ml_gam$gam, newdata=xx,type= "response"))

#GDD
xx_gdd = as.vector(predict.gam(gdd_gam$gam, newdata=xx,type= "response"))

#Min temp 
xx_gmt = as.vector(predict.gam(gmt_gam$gam, newdata=xx,type= "response"))
xx_gmt = predict( lm(xx_gmt~xx$elevation)) #Uncomment to linearize xx_gmt

#Make the internal spatial variables
xx=cbind(xx,xx_mt,xx_sm,xx_sp,xx_sd,xx_ml,xx_gdd, xx_gmt)
colnames(xx)[4:10] = c("gs_mean_temp","soil_moist","soil_pet",
	"soil_def","an_mean_light","max_gdd","gs_min_temp")


xx$gs2 = as.matrix(xx_mt^2)
xx$sm2 = as.matrix(xx_sm^2)
xx$sp2 = as.matrix(xx_sp^2)
xx$sd2 = as.matrix(xx_sd^2)

#Spatial scale
x1=  as.matrix(seq(el1,el2,1))
colnames(x1)=c('elevation')
x1=data.frame(x1)

spp = unique(as.character(allB2$Sp[!is.na(allB2$Sp)]))
nspp = length(spp)
s.index1= 1:nspp

bgs = unique(as.character(allB2$bg))
pop_labels = paste(expand.grid(spp,bgs,elevations)[,1],expand.grid(spp,bgs,elevations)[,
		2],expand.grid(spp,bgs,elevations)[,3],sep="")

ns=dim(xx)[1]-1 #Lattice width/height -- This should be even

#Make these automatically:
xx0=matrix(seq((-ns/2),(ns/2)))
np=length(xx0)
ngenst=iconfig+ngens


#=============================================================================
#survival 
#=============================================================================
#-- this section includes data exploration across the backgrounds, looking
#	for varios inter/intraspecific effects on survival. 
#=============================================================================
#See calanda_pop_krigALL_wald.R for the full statistical exploration

#variable kriged survival on only bare
sr_krigB = matrix(0,dim(xx)[1], length(spp))
sr_probB=NULL
kn = c(4,4,4)
sr=matrix(0,nspp,1)
raw_sr_all=vector("list", nspp)
for(sp in 1:length(spp)){
	allsp = subset(allB2, Sp == spp[sp])

	##1) Survival with background as fixed effect
	# sr_dat = allsp
	# sr_dat = subset(sr_dat, elevation>=el_real[1] & elevation<=el_real[5])
	# sr_gam = gam(survival~bg+s(elevation,k=kn[sp])+s(year,bs="re"),family=binomial(link='logit'), data=sr_dat)
	# sr_tmp=as.vector(predict.gam(sr_gam, newdata=data.frame(xx,bg = "B"),type= "response"))

	###2) Survival conditione only on non-competition sites
	sr_dat=subset(allsp, bg =='B')	
	sr_dat = subset(sr_dat, elevation>=el_real[1] & elevation<=el_real[5])
	sr_gam = gam(survival~s(elevation,k=kn[sp])+s(year,bs="re"),family=binomial(link='logit'), data=sr_dat)
	sr_tmp=as.vector(predict.gam(sr_gam, newdata=xx,type= "response"))

	#Remove 0s
	sr_tmp[sr_tmp<0] = 0
	sr_krigB[,sp] = sr_tmp
	sr_probB[[sp]] = sr_gam

	### Raw survival probabilities
 	raw_sr= table(sr_dat$survival,sr_dat$elevation)
 	raw_sr = as.data.frame(raw_sr/matrix(colSums(raw_sr),2,dim(raw_sr)[2],byrow=T ))
	raw_sr = subset(raw_sr, Var1 == 1 )
	raw_sr[,1] = as.numeric(levels(raw_sr[,1]))[raw_sr[,1]]
	raw_sr[,2] = as.numeric(levels(raw_sr[,2]))[raw_sr[,2]]
	raw_sr_all[[sp]] = raw_sr

}

#Plot kriged values
# par(mfrow=c(3,1))
# for(sp in 1:3) {
# 	plot(sr_krigB[,sp],ylim=c(0,1))
# 	#points(sr_krig[,sp],col="red")
# }


#=============================================================================
#Intraspecific competition
#=============================================================================
#See calanda_pop_krigALL_wald.R for the full statistical exploration
#=============================================================================
#Use NLS to fit the functional relationship (leslie-gower?). Try it per elevation
#too, to see if it is significant. If it is, Krige it.
#This version allows for the matching of a the per-plot background individual
#to the intraspecific treatment in its plot (i.e. matching the bg 'B' to the 
#IGR of the competitve treatment with the same 7 character prefix)

Crr_nls = matrix(0,length(spp),1) #Elevations grouped
Crr_nls_krig = matrix(0,dim(xx)[1], length(spp)) #Krig across elevations
Crr_nls_el = matrix(0,length(el_real), length(spp)) #Elevations individually

kn = c(4,4,4)
bsamp = 10000

for(sp in 1:length(spp)){

	allsp = subset(allB2, Sp == spp[sp] )	
	Cr_dat =subset(allsp, bg == substr(as.character(spp[sp]),1,1))
	Cr_dat= subset(Cr_dat, elevation>=el_real[1] & elevation<=el_real[5])
	
	rr_sub =subset(allsp, bg =='B')
	rr_sub = subset(rr_sub, elevation>=el_real[1] & elevation<=el_real[5])
	
	if(sp ==1 ) { 
	
		Cr_dat$nr.inflor = Cr_dat$nr.inflor*Cr_dat$stalks	
		rr_sub$nr.inflor = rr_sub$nr.inflor*rr_sub$stalks


	}
	nsamp = dim(rr_sub)[1]

	#Bootstrap this by randomly pairing bare and competitive backgrounds
	#This would ignore elevation as a treatment
	crdat= data.frame(cbind(Cr_dat$nr.inflor[ceiling(runif(bsamp)*nsamp)], 
		rr_sub$nr.inflor[ceiling(runif(bsamp)*nsamp)]))
	
	colnames(crdat) = c("lrr","Rr")
	
	#Fit with nls
	Cr_nls = nls(lrr~ Rr/(1+cr),data=crdat, start=list(cr=0.5))
	print(summary(Cr_nls))
	Crr_nls[sp] = Cr_nls$m$getPars()
	
	#Now, per elevation: 
	#plot(NA, NA, ylim=c(0, 30),xlim = c(0, 20)) 
	for(ee in 1:length(el_real)){
		Cr_dat_e = subset(Cr_dat, elevation == el_real[ee])
		rr_sub_e = subset(rr_sub, elevation == el_real[ee])
		nsamp = dim(rr_sub_e)[1]

		#Bootstrap this by randomly pairing bare and competitive backgrounds
		#This would maintain within-site pairing
		crdat_e= data.frame(cbind(Cr_dat_e$nr.inflor[ceiling(runif(bsamp)*nsamp)], 
			rr_sub_e$nr.inflor[ceiling(runif(bsamp)*nsamp)]))
		colnames(crdat_e) = c("lrr","Rr")
		#points(crdat_e$Rr ,crdat_e$lrr )

		#Fit with nls
		tryCatch( {Cr_nls_e = nls(lrr~ Rr/(1+cr),data=crdat_e, start=list(cr=0.5),control = list(maxiter = 500))
			Crr_nls_el[ee, sp] = Cr_nls_e$m$getPars()}, error = function(e){})
		
	}
	
	#Now use the per elevation data to krige
	Cr_k = data.frame(cbind(Crr_nls_el[, sp], el_real))
	colnames(Cr_k) = c("cr","elevation")
	Cr_gam = gam(cr~s(elevation,k=kn[sp]),data=Cr_k )
	Cr_tmp=as.vector(predict.gam(Cr_gam, newdata=x1))

	#Remove 0s
	Cr_tmp[Cr_tmp<0] = 0
	Crr_nls_krig[,sp] = Cr_tmp
}

#=============================================================================
#interspecific growth
#=============================================================================
#See calanda_pop_krigALL_wald.R for the full statistical exploration
#=============================================================================
#Use NLS to fit the functional relationship (leslie-gower?). Try it per elevation
#too, to see if it is significant. If it is, Krige it. If not...? 
#Allows for the matching of a the per-plot background individual
#to each interspecific treatment in its plot in turn (i.e. matching the bg 'B' to the 
#IGR of the competitve treatment with the same 7 character prefix)


Cir_nls = matrix(0,length(spp),length(spp)) #Elevations grouped
Cir_nls_all = NULL
Cir_el_all = NULL

kn = c(4,4,4)

for(sp in 1:length(spp)){

	Cir_nls_krig = matrix(0,dim(xx)[1], length(spp)) #Krig across elevations
	Cir_nls_el = matrix(0,length(el_real), length(spp)) #Elevations individually

	sp_use = 1:(length(spp))
	sp_use= sp_use[-sp]

	for(sa in 1: (length(spp)-1)){
		allsp = subset(allB2, Sp == spp[sp])
		rr_sub =subset(allsp, bg =='B')
		rr_sub = subset(rr_sub, elevation>=el_real[1] & elevation<=el_real[5])

		nsamp = dim(rr_sub)[1]

		sb= sp_use[sa]
		Ci_dat =subset(allsp, bg == substr(as.character(spp[sb]),1,1))
		Ci_dat= subset(Ci_dat, elevation>=el_real[1] & elevation<=el_real[5])
		
		if(sp ==1 ) { 
	
			Ci_dat$nr.inflor = Ci_dat$nr.inflor*Ci_dat$stalks	
			rr_sub$nr.inflor = rr_sub$nr.inflor*rr_sub$stalks	
		}

		#Bootstrap this by randomly pairing bare and competitive backgrounds
		#This would ignore elevation as a treatment
		cidat= data.frame(cbind(Ci_dat$nr.inflor[ceiling(runif(bsamp)*nsamp)],
		 rr_sub$nr.inflor[ceiling(runif(bsamp)*nsamp)]))
		
		colnames(cidat) = c("lir","Ri")
		
		#Fit with nls
		tryCatch( {		Ci_nls = nls(lir~ Ri/(1+ci),data=cidat, start=list(ci=Crr_nls[sp]),control = list(maxiter = 500))
				Cir_nls[sp, sb] = Ci_nls$m$getPars()
				}, error = function(e){})
		#Ci_nls = nls(lir~ Ri/(1+ci),data=cidat, start=list(ci=Crr_nls[sp]),control = list(maxiter = 500))
		print(summary(Ci_nls))
		Cir_nls[sp,sb] = Ci_nls$m$getPars()
		
		#Now, per elevation: 
		#plot(NA, NA, ylim=c(0, 30),xlim = c(0, 20)) 
		for(ee in 1:length(el_real)){
			Ci_dat_e = subset(Ci_dat, elevation == el_real[ee])
			rr_sub_e = subset(rr_sub, elevation == el_real[ee])
			nsamp = dim(rr_sub_e)[1]

			#Bootstrap this by randomly pairing bare and competitive backgrounds
			#This would maintain within-site pairing
			cidat_e= data.frame(cbind(Ci_dat_e$nr.inflor[ceiling(runif(bsamp)*nsamp)],
			 rr_sub_e$nr.inflor[ceiling(runif(bsamp)*nsamp)]))
			colnames(cidat_e) = c("lir","Ri")
			#points(cidat_e$Ri ,cidat_e$lir )
			
			#Fit with nls
			tryCatch( {Ci_nls_e = nls(lir~ Ri/(1+ci),data=cidat_e, start=list(ci=Crr_nls[sp]),control = list(maxiter = 500))
				Cir_nls_el[ee, sb] = Ci_nls_e$m$getPars()
				}, error = function(e){})
			
		}
		mean(Ci_dat$nr.inflor,na.rm=T)
		#Now use the per elevation data to krige
		Ci_k = data.frame(cbind(Cir_nls_el[, sb], el_real))
		colnames(Ci_k) = c("ci","elevation")
		Ci_gam = gam(ci~s(elevation,k=kn[sp]),data=Ci_k )
		Ci_tmp=as.vector(predict.gam(Ci_gam, newdata=x1))
		#Remove 0s
		Ci_tmp[Ci_tmp<0] = 0
		Cir_nls_krig[,sb] = Ci_tmp
	}
	Cir_nls_all[[sp]] = list(Cir_nls_krig)
	Cir_el_all[[sp]] = list(Cir_nls_el)
}



#=============================================================================
#air/arr
#=============================================================================
#See calanda_pop_krigALL_wald.R for the full statistical exploration
#=============================================================================
#=============================================================================
#Get competition (flowers) per site. This is based on solving each IGR for arr*Nr 
#and air*Nr, then taking the ratio air/arr, setting arr=1
#For flowers, this doesn't include the survival term so: 
#air/arr = (m*(n+li))/((m+lr)*n)
#Version 3: use the constant alphas fit with NLS

arat3 = matrix(0,length(spp),length(spp))

for(sp in 1:length(spp)){
	sp_use = 1:(length(spp))
	sp_use= sp_use[-sp]

	for(sa in 1: (length(spp)-1)){
		sb= sp_use[sa]
		#Use medians of IGRs
		print(Crr_nls[sb])
		arat3[sp,sb]=(Cir_nls[sp,sb])/(Crr_nls[sb])
		
	}
}		
	

#=============================================================================
#Get competition (flowers) per site. This is based on solving each IGR for arr*Nr 
#and air*Nr, then taking the ration air/arr, setting arr=1
#For flowers, this doesn't include the survival term so: 
#air/arr = (m*(n+li))/((m+lr)*n)
#Version 3B: use an average over the elvational alphas fit with NLS

arat3B = matrix(0,length(spp),length(spp))

for(sp in 1:length(spp)){
	sp_use = 1:(length(spp))
	sp_use= sp_use[-sp]
	
	
	for(sa in 1: (length(spp)-1)){
		sb= sp_use[sa]
		crr = Crr_nls_el[,sb]
		cir_tmp = Cir_el_all[[sp]][[1]][,sb]
		cir_use = cbind(crr,cir_tmp)
		cir_use =cir_use [cir_use[,1]>0,]
		cir_e = cir_use[,2]/cir_use[,1]
		#cir = mean(cir_tmp[cir_tmp>0])
		arat3B[sp,sb]=mean(cir_e,na.rm=T)
		
	}
}		
	

#=============================================================================
#Invasion growth rates
#=============================================================================
#=============================================================================
#Tunable Species parameters
#=============================================================================
#Calculate spatiotemporal Fr (fundamental niche)
###Maximum reproduction rates, corresponding to spatio-temporal ideal conditions

###
###Survival
###
sr=colMeans(sr_krigB[800:2000,])
#sr[3]=sr[3]+0.05
#From raw probability data: 
#sr =c(mean(raw_sr_all[[1]][,3]), mean(raw_sr_all[[2]][,3]), mean(raw_sr_all[[3]][,3]))

###
###Competition coeffcients
###
#3B is calculating a per-site alpha, then averaging 
# diag(arat3B)=c(1,1,1)
# alphas=arat3B
# alphas[2,1] = arat3[2,1]

#3 is calculating an alpha that has been calculated regardless
#of elevation
diag(arat3)=c(1,1,1)
alphas=arat3

###
###Competition distance -- These are unknown, and chosen just to be limited or not
###
#b_rr=1/(100*np) #Essentially global
b_rr=c(.1,.1,.1) #10 cm plots

###
###Dispersal distance. Used for global dispersal in this version. 
a_rr= c(1/(100*np),1/(100*np),1/(100*np) ) #essentially global

#=============================================================================
# Internal variables: 1D 
#=============================================================================

#Make the full array of space-time coordinates (for 1D space)
#meshgrid will produce an R data frame
stc.temp=meshgrid(1:ngenst,xx0)
#Convert the coordinates into an array for easier access
stc=array(c(matrix(stc.temp$x,ngenst,np,byrow=T),matrix(stc.temp$y,ngenst,np,byrow=T)),dim=c(ngenst,np,2)) 
#Note: in 1d, things are a bit odd. stc.temp$x ends up as time, and $y as space. So in stc, stc[,,1] is time
# and stc[,,2] is space. stc[,,1] is the same as the column vector 1:ngenst repeated over np columns. Then
# stc[,,2] is space as a row vector (xx) repeated over ngenst rows. 


#=============================================================================
#dispersal
#=============================================================================

#Exponential kernel
kd=matrix(0,np,nspp)
fkd=matrix(0,np,nspp)

for( s in 1:nspp){ 
	kd[,s] = a_rr[s]/2*exp(-a_rr[s]*abs(xx0))
	kd[,s]=kd[,s]/(sum(kd[,s]))
	fkd[,s]=fft(kd[,s])#/(np+1)
	fkd.yes = TRUE #Pass the transformed kd to functions later for speed

}


# kd=matrix(0,np,nspp)
# fkd=matrix(0,np,nspp)

# #Delta function kernel
# for( s in 1:nspp){ 
# 	kd[ceiling(np/2) ,s] = 1
# 	fkd[,s]=fft(kd[,s])#/(np+1)
# 	fkd.yes = TRUE #Pass the transformed kd to functions later for speed
# }


####Competition kernels and their Fourier transforms 

#Exponential
kc=matrix(0,np,nspp)
fkc=matrix(0,np,nspp)
for( s in 1:nspp){ 
	kc[,s] = b_rr[s]/2*exp(-b_rr[s]*abs(xx0))
	kc[,s]=kc[,s]/(sum(kc[,s]))
	fkc[,s]=fft(kc[,s])#/(np+1)
}

#Delta function
# kc=matrix(0,np,nspp)
# fkc=matrix(0,np,nspp)
# #Delta function kernel
# for( s in 1:nspp){ 
# 	kc[ceiling(np/2) ,s] = 1
# 	fkc[,s]=fft(kc[,s])#/(np+1)
	
# }

#=============================================================================
# Key variables for output: the invasion growth rates, their components
#=============================================================================
#Intrinsic ranges
Frs = array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp)) 
Frs_use = array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp)) 
Frs_mlp = array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp)) 
Frs_post = vector("list", ngenst) 
Frs_ci = vector("list", ngenst) 
#The maximum of species' fitted reproduction, for scaling
frs_max = matrix(0,ngenst,nspp)

#Realized ranges in full commmunity
nf = array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp))  

#Equilibrium populations of residents. Columns represent space, rows are time. 
nr = array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(n_runs,ngenst,np,nspp)) 

#Simulation IGR
ldg.sim = array( c( matrix(0,ngenst, nspp),matrix(0,ngenst, nspp)),dim=c(n_runs,ngenst,nspp) )
ldgIMP.sim = array( c( matrix(0,ngenst, nspp),matrix(0,ngenst, nspp)),dim=c(n_runs,ngenst,nspp) )
covIMP.sim =array( c( matrix(0,ngenst, nspp),matrix(0,ngenst, nspp)),dim=c(n_runs,ngenst,nspp) )

#Components of IGR

#Keep track of total interspecific competition of species when it is invader
cc=array(c(matrix(0,ngenst,np),matrix(0,ngenst,np)),dim=c(n_runs,ngenst,np,nspp)) 

l1=array( c( matrix(0,ngenst, nspp),matrix(0,ngenst, nspp)),dim=c(n_runs,ngenst,nspp) )
D = array( c( matrix(0,ngenst, nspp),matrix(0,ngenst, nspp)),dim=c(n_runs,ngenst,nspp) )
var_mu_Us=array( c( matrix(0,ngenst, nspp),matrix(0,ngenst, nspp)),dim=c(n_runs,ngenst,nspp) )
cov_e_mu_Us=array( c( matrix(0,ngenst, nspp),matrix(0,ngenst, nspp)),dim=c(n_runs,ngenst,nspp) )
cov_lam_vc=array( c( matrix(0,ngenst, nspp),matrix(0,ngenst, nspp)),dim=c(n_runs,ngenst,nspp) )
cov_lam_vc2=array( c( matrix(0,ngenst, nspp),matrix(0,ngenst, nspp)),dim=c(n_runs,ngenst,nspp) )
Elam1=array( c( matrix(0,ngenst, nspp),matrix(0,ngenst, nspp)),dim=c(n_runs,ngenst,nspp) )
Elam2=array( c( matrix(0,ngenst, nspp),matrix(0,ngenst, nspp)),dim=c(n_runs,ngenst,nspp) )
gr1.n=array( c( matrix(0,ngenst, nspp),matrix(0,ngenst, nspp)),dim=c(n_runs,ngenst,nspp) )
gr1=array( c( matrix(0,ngenst, nspp),matrix(0,ngenst, nspp)),dim=c(n_runs,ngenst,nspp) )

#Components kept as lists
y.full=vector(ttot, mode="list")
w.eq.full=vector(ttot, mode="list")

#=============================================================================
#Outermost loop: treat each species as invader
#=============================================================================

for( t in 1: ngenst){ 

	#set.seed(2)
	ti = seq(tstart,tstop,tinc)[t]
	temp_i = ccs_temp[ti]
	soil_i = ccs_soil[ti]

	xx_new = xx

	xx_new$gs_mean_temp = as.numeric(xx_new$gs_mean_temp+temp_i)
	xx_new$gs2 = as.numeric(xx_new$gs_mean_temp^2)

	xx_new$soil_moist = as.numeric(xx_new$soil_moist+soil_i)

	#This code takes the temp at sites in the future and matches it to
	#the equivalent current site based on that temp.
	new_ind=sapply(xx_new$gs_mean_temp,function(x) which.min(abs(x-xx$gs_mean_temp)))
	xx_new$max_gdd = as.numeric(xx$max_gdd[new_ind])
	xx_new$gs_min_temp = as.numeric(xx$gs_min_temp[new_ind])

	#env.ind = c(1,4,7,9,10)
	env.ind = c(4,5)

	#Environmental distance, based on novelty metric in 
	#env_distance = get_env_distance(xx,xx_new,env.ind)
	#env_analogue = get_analogue(xx,xx_new,env.ind)

	#=============================================================================
	#FLOWERS
	#=============================================================================

	#=============================================================================
	#intrinsic growth
	#=============================================================================
	###
	###Intrinsic probability of flowering
	###
	flower_probK=matrix(0,dim(xx)[1], length(spp))
	flower_prob_ci=NULL
	flower_prob_bu=NULL
	flower_prob=NULL
	flower_prob_post=NULL
	flower_prob_act = matrix(0,length(el_real),3)

	for(sp in 1:length(spp)){
		allsp = subset(allB2, Sp == spp[sp])
		rr_dat =subset(allsp, bg =='B')
		rr_dat = subset (rr_dat, elevation <= el_real[5] & elevation >= el_real[1])
		
		nvar = nvar_probs[sp]
		v_use =  v_use_prob [[sp]]

		if(sp ==1 ) { 
		
			allsp$nr.inflor = allsp$nr.inflor*allsp$stalks	
			rr_dat$nr.inflor = rr_dat$nr.inflor*rr_dat$stalks	
		}

		### Make three important lists: spline fits, knots, penalty matrixes
		sm=list()
		knots_sm = list()
		para_pen = list()

		#Build the main data set with all of the variables
		y = c(matrix(0, 1,1),rr_dat$flyes,matrix(0, 1,1))
		yearx=unlist(list( factor(matrix(years,1,1)),rr_dat$year,factor(matrix(years,1,1))))
		dat = data.frame (flyes = rr_dat$flyes, rr_dat[paste(v_use[(1:length(v_use))])])
		dzeros = data.frame ( matrix(0, 1, ncol(dat)))
		colnames(dzeros) = colnames(dat)
		dat = rbind(dzeros,dat,dzeros)
		dat$year = yearx

		#This section corresponds to fitting the constrained GAMMs. First, 
		#a series of knots are specified for each smooth fit in each 
		#species' full model in order to identify the appropriate offsets 
		#thatconstrain the intercepts.This loop goes through each variable 
		#in sequence:  

		for( v in 1:nvar){ 
			#Pick out the variable of interest:
			iv_name = paste("rr_dat[[\"",v_use[(v+1)],"\"]]", sep="")
			iv_name_pred = paste("xx_new[[\"",v_use[(v+1)],"\"]]", sep="")
			#iv1 =  paste(v_use[(v+1)],sep="")
			iv1=eval(parse(text=iv_name))
			#assign(iv1,eval(as.name(iv_name)))

			#### Ensuring that smooth fits taper to 0 at upper and lower values
			#By default, mgcv::gam places a knot at the extremes of the data and then the 
			#remaining "knots" are spread evenly over the interval.
			kuse= var_spp_knots[[sp]][v,1] #Effectively the number of knots
			nk = 2 #How many knots to add at ends to contstrain smooth? 
			k1 = unique(iv1)
			knots = seq(min(k1),max(k1),length=kuse)
			kby = diff(knots)[1]

			###Determine the ideal spacing for knots, and especially the constraints at
			###the endpoints: 
			iv1_pred = eval(parse(text=iv_name_pred))

			### Find the largest and smallest values 
			# iv1_lims = c(min( min(iv1_pred),min(iv1)),max(max(iv1_pred),max(iv1)))
			
			# #If the smallest/largest values come from the actual dataset, then use 
			# #the kby value to set limits:
			# # if( min (iv1_lims) == min(iv1)){iv1_lims[1] = min(iv1) - kby }
			# # if( max (iv1_lims) == max(iv1)){iv1_lims[2] = max(iv1) + kby }

			# iv1_lims[1] = min(iv1) - kby
			# iv1_lims[2] = max(iv1) + kby 
			
			iv1_lims = c(min( min(iv1)), max(iv1))
			iv1_lims[1] = min(iv1) - 1
			iv1_lims[2] = max(iv1) + 1 

			knots= data.frame(x=c( iv1_lims[1],knots,iv1_lims[2] ))
			#knots= data.frame(x=c( knots[1]-1,knots,knots[kuse]+1 ))
			xk = dim(knots)[1]
			knots_sm[[paste(v_use[v+1])]] = list(knots$x)


			### Build the data set
			x = c( matrix(knots$x[1],1,1), iv1, matrix(knots$x[xk],1,1))
			dat.tmp = data.frame(x=x)
			dat[[paste(v_use[v+1])]] = x #Add the edited column to the data frame
			
			cp_x = c(1, xk ) #x Values to drop, i.e. the extremes

			## set up smoother... 
			sm[[v]] = smoothCon(s(x,k=dim(knots)[1],bs="cr"),dat.tmp,knots=knots)[[1]] 
			sm[[v]]$term = 	v_use[v+1]
			sm[[v]]$label = paste("s(",v_use[v+1], ")" )
			## set points to control the to approach 0 by dropping... 
			X.name = paste("X",v,sep="")
			S.name = paste("S",v,sep="")
			assign(X.name,sm[[v]]$X[,-cp_x])        ## spline basis 
			assign(S.name,sm[[v]]$S[[1]][-(cp_x),-(cp_x)]) ## spline penalty 
			S = sm[[v]]$S[[1]][-(cp_x),-(cp_x)]
			para_pen[[paste(X.name)]] = list(S)

		}

		factors = matrix("f",nvar)
		factors2 = matrix("f",nvar)

		for (f in 1:nvar){
			factors[f] = paste("X", f, sep="")
			factors2[f] = paste("s(",v_use[(f+1)], ", k = ", var_spp_knots[[sp]][f,1]+nk, ",  bs=\"cr\")", sep="")

		} 
		## In order to find the minimum intercept (which is not zero with the logit link function): 
		## First fit a model to the data: 
		b.u = gam(as.formula(paste("flyes~", paste( paste(factors2, collapse="+"),"+s(year, bs=\"re\")"))), data=dat, knots = knots_sm, 
			family=binomial(link='logit') )
		bu_tmp = predict(b.u, exclude = "s(year)") #We need the minimum of this 
		cp_y = c(matrix(min(bu_tmp),length(cp_x),1)) #y Values to force the intercept to
		off = dat$flyes*0 + cp_y[1] ## offset term to force curve through a point

		rr_gam = gam(as.formula(paste("flyes~", paste( paste(factors, collapse="+"),"+offset(off)+s(year, bs=\"re\")-1"))), 
			paraPen=para_pen, data=dat, family=binomial(link='logit'),method = "REML" )
		# rr_gam = gam(as.formula(paste("flyes~", paste( paste(factors, collapse="+"),"+offset(off)+s(year, bs=\"re\")-1"))), 
		# 	paraPen=para_pen, data=dat, family=binomial(link='logit'))
	
		rr_tmp = predict(rr_gam, exclude = "s(year)")
		#rr_tmp = predict(rr_gam, type="response", exclude = "s(year)")
		
		### coefficients from the penalized regression	
		# ncs = c(0,0)
		# if (coef(rr_gam)[2]>=0){ ncs[1] = 0} else {ncs[1] = min(coef(rr_gam)[2:(kuse+1)])}
		# if (coef(rr_gam)[(kuse+1)]>=0){ ncs[2] = 0} else {ncs[2] = min(coef(rr_gam)[2:(kuse+1)])}
		# #beta = c(ncs[1],coef(rr_gam)[2:(kuse+1)],ncs[2])
		#beta = c(0,coef(rr_gam)[2:(kuse+1)],0) #With intercept
		kn = var_spp_knots[[sp]][ ,1]
		kuse = sum(kn[1:nvar])
		beta = matrix(0, 1, (kuse+nvar*nk) )
		#Covariance matrix for getting CIs 
		Vp = matrix(0,(kuse+nvar*nk),(kuse+nvar*nk) )

		pnk = nk/2
		knc=0
		for( b in 1:nvar){
			if(b>1) {knc = sum(kn[1:(b-1)])}
			bpos = (b-1)*(nk)+knc+pnk+1
			kpos = knc+1 #Intercept removed from rr_gam
			beta[ bpos:(bpos+kn[b]-1)] = coef(rr_gam) [kpos:(kpos+kn[b]-1)]
			Vp [bpos:(bpos+kn[b]-1),bpos:(bpos+kn[b]-1) ] = 
				rr_gam$Vp[kpos:(kpos+kn[b]-1),kpos:(kpos+kn[b]-1)]
	
		} 

		### prediction matrix --This step is equivalent to using gam.predict(...type = "lpmatrix")
		### It is necessary because our constrained model is not written in terms of the abiotic
		### variables, but in terms of the knots (e.g. do summary(b.u) vs. summary(rr_gam) )
		Xp = matrix(0, dim(xx_new)[1], (kuse+nvar*nk))
		#Xp = matrix(0, dim(dat)[1], nvar*(kuse+nk))
		knc=0
		for( b in 1:nvar){
			if(b>1) {knc = sum(kn[1:(b-1)])}
			bpos = (b-1)*(nk)+knc+pnk
			Xp_tmp = PredictMat(sm[[b]], xx_new)
			#Xp_tmp = PredictMat(sm[[b]], dat)
			Xp[,bpos:(bpos+kn[b]+nk-1)] = Xp_tmp
		}
		

		### the predicted smooths, i.e. species' distribution
		#flower_tmp = Xp%*%beta+coef(rr_gam)[1] #These are not on the response scale
		flower_tmp = Xp%*%t(beta)+ cp_y[1] #These are not on the response scale
		flower_tmp = exp(flower_tmp )/(1+exp(flower_tmp ) ) #Response scale for logistic regression
		
		#The unconstrained model, to compare the method for extracting SDs, etc.: 
		# bu_tmp=as.vector(predict(b.u,newdata=xx_new,type= "response",exclude = "s(year)", se.fit=TRUE ))
		# bu_tmp2=predict(b.u,newdata=xx_new,exclude = "s(year)", type="lpmatrix" )
		# vb = vcov(b.u)
		# b1 = coef(b.u)
		# br2=mvrnorm(n=n_draws,b1, vb)
		# ##Model SE: 
		# bu_sims = exp(bu_tmp2%*%t(br2) )/(1+exp(bu_tmp2%*%t(br2) ) )
		# bu_sd = apply(bu_sims, 1, sd)
		# ##Mean model: 
		# #bu_mean = rowMeans(bu_sims)
		# bu_m1 = rowMeans(bu_tmp2%*%t(br2))
		# bu_mean = exp(bu_m1 )/(1+exp(bu_m1  ) )


		#### To get the mean of the constrained model, with confidence intervals: 
		br1=mvrnorm(n=n_draws, beta, Vp)
		##This produces the mean model: 
		mean.model= rowMeans(Xp%*%t(br1)+cp_y[1])
		mean.model = exp(mean.model)/(1+exp(mean.model))
		flower_tmp = mean.model

		##To keep each posterior simulation and get CIs: 
		model_post= Xp%*%t(br1)+cp_y[1]
		model_post = exp(model_post)/(1+exp(model_post))
		#q1=quantile(mean.model,c(.025,.975))
		flower_ci_tmp = t(apply(model_post, 1, quantile, c(0.025,0.975)))
		flower_se_tmp = apply(model_post, 1, sd ) 

		#Mean model that matches CIs
		#flower_tmp = rowMeans(mean.model)

		#Save variables: 
		flower_probK[,sp] = flower_tmp
		#flower_probK[,sp] = bu_tmp
		flower_prob[[sp]] = rr_gam
		flower_prob_bu[[sp]] = b.u
		flower_prob_post[[sp]] = model_post
		flower_prob_ci[[sp]] = flower_ci_tmp
	}


	###
	### Total Flowers
	###

	#Krig a distribution of Rs (intrinsic growth) from flower counts of plants 
	#grown alone, per site 
	#Use only those that have flowered
	#The best approach here seems to be to use the probability by elevation, 
	#then each R is actually a constant.

	rr_krig = matrix(0,dim(xx)[1], length(spp))
	rr_ci=NULL
	rr_gams=NULL
	rr_gams_bu=NULL
	rr_post = NULL

	flower_act =matrix(0,length(el_real),3)

	for(sp in 1:length(spp)){
		allsp = subset(allB2, Sp == spp[sp])
		rr_dat_bg =subset(allsp, bg =='B')
		rr_dat =subset(rr_dat_bg, flyes ==1)
		nvar = nvar_flows[sp]
		v_use =  v_use_flow [[sp]]

		if(sp ==1 ) { 
			#Ci_dat$nr.inflor = Ci_dat$nr.inflor*Ci_dat$stalks+Ci_dat$nr.stms	
			#rr_sub$nr.inflor = rr_sub$nr.inflor*rr_sub$stalks+rr_sub$nr.stms	
			allsp$nr.inflor = allsp$nr.inflor*allsp$stalks	
			rr_dat$nr.inflor = rr_dat$nr.inflor*rr_dat$stalks	
		}


		for (n in 1:length(el_real)) { 
			tfp = subset(rr_dat, elevation == elevations[n])
			flower_act[n,sp]= mean(tfp$nr.inflor,na.rm=T)
		}
		flower_act[,sp][is.na(flower_act[,sp])] = 0


		### Make three important lists: spline fits, knots, penalty matrixes
		sm=list()
		knots_sm = list()
		para_pen = list()

		#Build the main data set with all of the variables
		y = c(matrix(0, 1,1),rr_dat$nr.inflor,matrix(0, 1,1))
		yearx=unlist(list( factor(matrix(years,1,1)),rr_dat$year,factor(matrix(years,1,1))))
		dat = data.frame (nr.inflor = rr_dat$nr.inflor, rr_dat[paste(v_use[(1:length(v_use))])])
		dzeros = data.frame ( matrix(0, 1, ncol(dat)))
		colnames(dzeros) = colnames(dat)
		dat = rbind(dzeros,dat,dzeros)
		dat$year = yearx

		#Are there any influential outliers? 
		cd_cut = 0.4 #Set a threhshold 
		cdlm = cooks.distance(lm(rr_dat$nr.inflor~rr_dat$gs_mean_temp))
		cd_remove = cdlm >=cd_cut
		names(cd_remove) = NULL
		#cd_remove =  c(matrix(FALSE, 3,1),cd_remove,matrix(FALSE, 3,1))
		cd_remove =  c(matrix(FALSE, 1,1),cd_remove,matrix(FALSE, 1,1))

		#Replace with the mean of remaining points to remove its influence
		if(sum(cd_remove)>0){dat[cd_remove,]$nr.inflor = mean(dat$nr.inflor[!cd_remove])}		


		#This section corresponds to fitting the constrained GAMMs. First, 
		#a series of knots are specified for each smooth fit in each 
		#species' full model in order to identify the appropriate offsets 
		#thatconstrain the intercepts.This loop goes through each variable 
		#in sequence:  

		for( v in 1:nvar){ 
			#Pick out the variable of interest:
			iv_name = paste("rr_dat[[\"",v_use[(v+1)],"\"]]", sep="")
			iv_name_pred = paste("xx_new[[\"",v_use[(v+1)],"\"]]", sep="")
			#iv1 =  paste(v_use[(v+1)],sep="")
			iv1=eval(parse(text=iv_name))
			#assign(iv1,eval(as.name(iv_name)))

			#### Ensuring that smooth fits taper to 0 at upper and lower values
			#By default, mgcv::gam places a knot at the extremes of the data and then the 
			#remaining "knots" are spread evenly over the interval.
			kuse= var_spp_knots[[sp]][v,2] #Effectively the number of knots
			nk = 2 #How many knots to add at ends to contstrain smooth? 
			k1 = unique(iv1)
			knots = seq(min(k1),max(k1),length=kuse)
			kby = diff(knots)[1]

			###Determine the ideal spacing for knots, and especially the constraints at
			###the endpoints: 
			iv1_pred = eval(parse(text=iv_name_pred))
			
			### Find the largest and smallest values 
			# iv1_lims = c(min( min(iv1_pred),min(iv1)),max(max(iv1_pred),max(iv1)))
			
			# #If the smallest/largest values come from the actual dataset, then use 
			# #the kby value to set limits:
			# # if( min (iv1_lims) == min(iv1)){iv1_lims[1] = min(iv1) - kby }
			# # if( max (iv1_lims) == max(iv1)){iv1_lims[2] = max(iv1) + kby }

			# iv1_lims[1] = min(iv1) - kby
			# iv1_lims[2] = max(iv1) + kby 
			
			iv1_lims = c(min( min(iv1)), max(iv1))
			iv1_lims[1] = min(iv1) - 1
			iv1_lims[2] = max(iv1) + 1 

			knots= data.frame(x=c( iv1_lims[1],knots,iv1_lims[2] ))
			#knots= data.frame(x=c( knots[1]-1,knots,knots[kuse]+1 ))
			xk = dim(knots)[1]
			knots_sm[[paste(v_use[v+1])]] = list(knots$x)

			### Build the data set
			x = c( matrix(knots$x[1],1,1), iv1, matrix(knots$x[xk],1,1))
			dat.tmp = data.frame(x=x)
			dat[[paste(v_use[v+1])]] = x #Add the edited column to the data frame
			
			cp_x = c(1, xk ) #x Values to drop, i.e. the extremes

			## set up smoother... 
			sm[[v]] = smoothCon(s(x,k=dim(knots)[1],bs="cr"),dat.tmp,knots=knots)[[1]] 
			sm[[v]]$term = 	v_use[v+1]
			sm[[v]]$label = paste("s(",v_use[v+1], ")" )
			## set points to control the to approach 0 by dropping... 
			X.name = paste("X",v,sep="")
			S.name = paste("S",v,sep="")
			assign(X.name,sm[[v]]$X[,-cp_x])        ## spline basis 
			assign(S.name,sm[[v]]$S[[1]][-(cp_x),-(cp_x)]) ## spline penalty 
			S = sm[[v]]$S[[1]][-(cp_x),-(cp_x)]
			para_pen[[paste(X.name)]] = list(S)

		}

		factors = matrix("f",nvar)
		factors2 = matrix("f",nvar)

		for (f in 1:nvar){
			factors[f] = paste("X", f, sep="")
			factors2[f] = paste("s(",v_use[(f+1)], ", k = ",var_spp_knots[[sp]][f,2]+nk, ",  bs=\"cr\")", sep="")
		} 

		## In order to find the minimum intercept (which is not zero with the logit link function): 
		## First fit a model to the data: 
		b.u = gam(as.formula(paste("nr.inflor~", paste( paste(factors2, collapse="+"),"+s(year, bs=\"re\")"))), data=dat, knots = knots_sm,
			family=nb() )
		# b.u = gam(as.formula(paste("nr.inflor~", paste( paste(factors2, collapse="+"),"+s(year, bs=\"re\")"))), data=dat, knots = knots_sm)
		bu_tmp = predict(b.u, exclude = "s(year)") #We need the minimum of this 
		cp_y = c(matrix(0,length(cp_x),1)) #y Values to force the intercept to
		off = dat$nr.inflor*0 + cp_y[1] ## offset term to force curve through a point

		rr_gam= gam(as.formula(paste("nr.inflor~", paste( paste(factors, collapse="+"),"+offset(off)+s(year, bs=\"re\")-1"))), 
			 paraPen=para_pen, data=dat,method = "REML",family=nb())
		# rr_gam= gam(as.formula(paste("nr.inflor~", paste( paste(factors, collapse="+"),"+offset(off)+s(year, bs=\"re\")-1"))), 
		# 	 paraPen=para_pen, data=dat,method = "REML")
	
		rr_tmp = predict(rr_gam, exclude = "s(year)")
	
		### coefficients from the penalized regression	
		# ncs = c(0,0)
		# if (coef(rr_gam)[2]>=0){ ncs[1] = 0} else {ncs[1] = min(coef(rr_gam)[2:(kuse+1)])}
		# if (coef(rr_gam)[(kuse+1)]>=0){ ncs[2] = 0} else {ncs[2] = min(coef(rr_gam)[2:(kuse+1)])}
		# #beta = c(ncs[1],coef(rr_gam)[2:(kuse+1)],ncs[2])
		#beta = c(0,coef(rr_gam)[2:(kuse+1)],0) #With intercept
		kn = var_spp_knots[[sp]][ ,2]
		kuse = sum(kn[1:nvar])
		beta = matrix(0, 1, (kuse+nvar*nk) )
		#Covariance matrix for getting CIs 
		Vp = matrix(0,(kuse+nvar*nk),(kuse+nvar*nk) )

		pnk = nk/2
		knc=0
		for( b in 1:nvar){
			if(b>1) {knc = sum(kn[1:(b-1)])}
			bpos = (b-1)*(nk)+knc+pnk+1
			kpos = knc+1 #Intercept removed from rr_gam
			beta[ bpos:(bpos+kn[b]-1)] = coef(rr_gam) [kpos:(kpos+kn[b]-1)]
			Vp [bpos:(bpos+kn[b]-1),bpos:(bpos+kn[b]-1) ] = 
				rr_gam$Vp[kpos:(kpos+kn[b]-1),kpos:(kpos+kn[b]-1)]
	
		} 

		### prediction matrix --This step is equivalent to using gam.predict(...type = "lpmatrix")
		### It is necessary because our constrained model is not written in terms of the abiotic
		### variables, but in terms of the knots (e.g. do summary(b.u) vs. summary(rr_gam) )
		Xp = matrix(0, dim(xx_new)[1], (kuse+nvar*nk))
		#Xp = matrix(0, dim(dat)[1], nvar*(kuse+nk))
		knc=0
		for( b in 1:nvar){
			if(b>1) {knc = sum(kn[1:(b-1)])}
			bpos = (b-1)*(nk)+knc+pnk
			Xp_tmp = PredictMat(sm[[b]], xx_new)
			#Xp_tmp = PredictMat(sm[[b]], dat)
			Xp[,bpos:(bpos+kn[b]+nk-1)] = Xp_tmp
		}
		
		### the predicted smooth, i.e. species' distribution
		#rr_tmp = Xp%*%beta+coef(rr_gam)[1] #These are not on the response scale
		#rr_tmp = Xp%*%t(beta) + cp_y[1]

		#The unconstrained model, to compare the method for extracting SDs, etc.: 
		bu_tmp=as.vector(predict(b.u,newdata=xx_new,type= "response",exclude = "s(year)", se.fit=TRUE ))
		# bu_tmp2=predict(b.u,newdata=xx_new,exclude = "s(year)", type="lpmatrix" )
		# vb = vcov(b.u)
		# b1 = coef(b.u)
		# br2=mvrnorm(n=n_draws,b1, vb)
		# ##Model SE: 
		# bu_sims = exp(bu_tmp2%*%t(br2) )
		# bu_sd = apply(bu_sims, 1, sd)
		# ##Mean model: 
		# #bu_mean = rowMeans(bu_sims)
		# bu_m1 = rowMeans(bu_tmp2%*%t(br2))
		# bu_mean = exp(bu_m1 )

		#### To get the mean of the constrained model, with confidence intervals: 
		br1=mvrnorm(n=n_draws, beta, Vp)
		##This produces the mean model: 
		mean.model= rowMeans(Xp%*%t(br1)+cp_y[1])
		mean.model = exp(mean.model)
		rr_tmp = mean.model

		##To keep each posterior simulation and get CIs: 
		model_post= Xp%*%t(br1)+cp_y[1]
		model_post = exp(model_post)
		#q1=quantile(mean.model,c(.025,.975))
		rr_ci_tmp = t(apply(model_post, 1, quantile, c(0.025,0.975)))
		rr_se_tmp = apply(model_post, 1, sd ) 

		#Save variables
		rr_tmp[rr_tmp<0] = 0
		rr_krig[,sp] = rr_tmp
		#rr_krig[,sp] = bu_tmp

		rr_gams[[sp]]=rr_gam
		rr_gams_bu[[sp]]=b.u
		rr_post[[sp]] = model_post
		rr_ci[[sp]] = rr_ci_tmp
	}


	###
	### Intrinsic ranges
	###
	# Frs=array(c(rr_krig*flower_probK),dim=c(ngenst,np,nspp)) 
	# #Frs[,,1] = Frs[,,1]*12
	# Frs[,,1] = Frs[,,1]*8
	#Frs[t,,] =matrix(rr_krig*flower_probK,np,nspp)
	Frs_temp =matrix(rr_krig*flower_probK,np,nspp)
	Frs[t,,] = matrix(rr_krig*flower_probK,np,nspp)
	for( sa in 1:nspp){

		print(max(flower_act[,sa],na.rm=T))
		frs_max[t,sa] = max(Frs_temp[,sa],na.rm=T)
		#Make the scaling factor a percentage of the non-shifted scaling
		max_scale =  frs_max[t,sa]/frs_max[1,sa]*max(flower_act[,sa],na.rm=T)
		#Scale all of the posterior distributions.
		Frs_post[[t]][[sa]] = flower_prob_post[[sa]]*rr_post[[sa]]*max_scale/frs_max[t,sa]
		#Mean value 
		Frs[t,,sa] = Frs[t,,sa]*max_scale/frs_max[t,sa]
		#Prediction intervals
		Frs_q1 = apply(Frs_post[[t]][[sa]], 1, quantile, c(0.025,0.975))
		Frs_ci[[t]][[sa]] = Frs_q1
	
	}

	Frs[t,,][Frs[t,,]<0] = 0
	Frs[t,,1] = Frs[t,,1]*3
	Frs_post[[t]][[1]] = Frs_post[[t]][[1]]*3

	# fig.name = paste("calanda_ranges",paste(v_use[-1],collapse=""),".pdf",sep="")
	# pdf(file=fig.name, height=11, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)
	col_use = c("black","red","blue")
	par(mfrow=c(3,1))
	plot(Frs[t,,2],t="l")
	#plot(mean_post[[1]],t="l",ylim=c(-20,60))

	for (tt in 1:3){
		lines(Frs[t,,tt],col=col_use[tt])
		#lines(mean_post[[tt]],col=col_use[tt])
		# lines(Frs_ci[[t]][[tt]][1,], col=col_use[tt])
		# lines(Frs_ci[[t]][[tt]][2,], col=col_use[tt])	
	
	}

	plot(flower_probK[,2],t="l",ylim=c(0,1))
	for (tt in 1:3){
		lines(flower_probK[,tt],col=col_use[tt])
		# lines(flower_prob_ci[[tt]][,1],col=col_use[tt] )	
		# lines(flower_prob_ci[[tt]][,2],col=col_use[tt] )	

	}


	plot(rr_krig[,2],t="l")
	for (tt in 1:3){
		lines(rr_krig[,tt],col=col_use[tt])	
		# lines(rr_ci[[tt]][,1],col=col_use[tt] )	
		# lines(rr_ci[[tt]][,2],col=col_use[tt] )	

	}


	#dev.off()

	# env_var[[t]] = xx_new[,env.ind] }

	#=============================================================================
	#=============================================================================

	#=============================================================================
	#Invader and resident stationary distributions
	#=============================================================================
	
	###
	### In this version, draw from the posterior simulations of 
	### intrinsic reproduction and repeat to generate mean and PIs
	###
	ptm = proc.time()
	draw_from = floor(runif(n_runs,0, n_draws) )
	for (d in 1:n_runs){

	
		for( sa in (1:nspp)){ 
			#Make the first run the mean model: 
			if( d == 1){
				#Use the mean: 
				Frs_use[t,,sa] = Frs[t,,sa] 

			} else { 
				#Use a random draw: 
				Frs_use[t,,sa]  = Frs_post[[t]][[sa]][,draw_from[d]]
				#Use the mean: 
				#Frs[t,,sa] = rowMeans(Frs_post[[t]][[sa]])
				Frs_use[t,,sa] [Frs_use[t,,sa] <0] = 0
			}

		}

		y = array( c(matrix(0,nspp, nspp-1),matrix(0,nspp, nspp-1) ), dim = c(n_runs,nspp,nspp-1) )
		w.eq = array( c( matrix(0,nspp, nspp-1),matrix(0,nspp, nspp-1)),dim = c(n_runs,nspp,nspp-1) )

		for ( s in 1:nspp) { 

			if( fast.bys==TRUE & (t > iconfig) ){   
				
					#Calculate the IGR for species 1 for only this time step
					gr1.fast=get.fast.igr(Frs[(t-1),,],nr[(t-1),,], sr, alphas, fkd,fkc,fkd.yes=TRUE, n.inv=s)
					#Now use the IGR to calculate the spread rate
					cs_all = get.spread.rate(gr1.fast,a_rr,sr)	

					# Make the new intrinsic fitness distribution for the next timestep
					# First, identify a leading edge of the intrinsic range based on  
					# a "finite individial" limit 
					Fr.inv=Frs[(t-1),,s]
					Fr.finite=as.numeric(Fr.inv>c.tol)
					fforw.up = max(which(Fr.finite ==1) )
					fforw.down = min(which(Fr.finite ==1) )
					
					#Shift it by the spread rate cs_all
					if ((fforw.up+ceiling(cs_all))<np){ 
						Fr.finite[fforw.up:(fforw.up+ceiling(cs_all))] = 1}
					if ((fforw.down-ceiling(cs_all))>0){ 
						Fr.finite[(fforw.down-ceiling(cs_all)):fforw.down] = 1}

					# Determine how much of the intrinsic range is being 
					# chopped off because low population spread rates prevent species from
					# tracking geographic shifts in their intrinsic ranges.

					Frs[t,,s] = Frs[t,,s]*Fr.finite

					print("Fast.bys")
						par(mfrow=c(1,1))
						plot(Frs[t,,1],t="l")
						for (tt in 2:3){
							lines(Frs[t,,tt])	
						}		

			}

			s.index = s.index1[-s]

			# #####Step 1
			# Frs.spp=Frs[t,,]
			# sr.spp =sr
			# alpha.spp=alphas
			# kd.spp = fkd 
			# kc.spp = fkc
			# fkd.yes
			# fast=TRUE
			# burns = 5

			# ###Step 2
			# np = dim(as.matrix(Frs.spp))[1] #Space
			# nspp = length(s.index) #Number of residents
			# nr.spp = s.index
			# Frs.spp = Frs.spp[,s.index]
			# sr.spp = sr.spp[s.index]
			# nr.burns=array(c(matrix(0.0,burns,np),matrix(0.0,burns,np)),dim=c(burns,np,nspp)) 
			# #Initialize residents
			# for ( sa in 1:nspp) { nr.burns[1,,sa] = matrix(0.01,1,np)}
			# b=1


			#Get the equilibrium of the resident community when invader is absent
			nr[d,t,,(s.index)] = get.res.eq(Frs_use[t,,],s.index,sr,alphas, fkd,fkc, fkd.yes,fast=TRUE, burns = 5 )
		
			#Get the low-density equilibrium density of the invader against the resident community
			nr[d,t,,s] = get.inv.ldeq(Frs_use[t,,], nr[d,t,,], s, sr,alphas, fkd, fkc,fkd.yes,fast=TRUE, burns = 5 )
		
			#=============================================================================
			# Low-density growth rates -- using simulation data, all 3 spp
			#=============================================================================
			#Simple numerical LGR: use the invader low-density equilibrium density, 
			#then let it invade against the community: 

			sp.ids =c(s, s.index) #IDs of resident, and invaders (this way uses all species)

			inv.one =pop_lg(Frs_use[t,,],nr[d,t,,], sr, alphas, fkd,fkc, fkd.yes )[,s]
			ldg.sim[d,t,s]=mean(inv.one)/mean(nr[d,t,,s])

			#Calculate the spatially implicit portion and the fitness-density covariance
			#separately, still using simulations. This just uses a "global" dispersal kernel
			####Dispersal kernels and their Fourier transforms 
			inv.id =1e-5

			kdg=matrix(0,np,nspp)
			fkdg=matrix(0,np,nspp)
			for( sa in 1:nspp){ 
				kdg[,sa] = 1/np
				fkdg[,sa]=fft(kdg[,sa])#/(np+1)
			}
			
			nr2=nr
			#Make the invader spatially homogenous, at low density
			nr2[d,,,s] = matrix(inv.id,ngenst,np)

			inv.oneIMP =pop_lg(Frs_use[t,,], nr2[d,t,,], sr, alphas, kdg,kc )[,s]
			ldgIMP.sim[d,t,s]=mean(inv.oneIMP)/mean(nr2[d,t,,s])

			covIMP.sim[d,t,s] = ldg.sim[d,t,s]-ldgIMP.sim[d,t,s]

			#=============================================================================
			#Calculate standard spatial persistence/coexistence terms
			#=============================================================================
			nrs= nr[d,t,,s.index] #Residents 
			
			#isd.tmp = get.isd(Frs[t,,s], nrs,s, alphas, sr[s], kc, kd[,s])

			y[d,s,] = colMeans(nrs,na.rm=T) #Means of resident densities 
			w.eq[d,s,] = alphas[s,s.index]*y[d,s,] #Weights for the LDG

			l1[d,t,s]=mean(Frs_use[t,,s])
			D[d,t,s]= 1+sum(w.eq[d,s,])


			#Calculate the standardized competition from residents
			muj=nrs/(matrix(y[d,s,],np,( nspp-1),byrow=T))-1
			uijmuj = matrix(0,np,nspp-1)
			for(a in 1:(nspp-1)) {
					ns = s.index[a]
					uijmuj[,a] = convolve(muj[,a],kc[,ns])
					uijmuj[,a]= c(uijmuj[ceiling(np/2):(np),a], uijmuj[1:floor(np/2),a] )
					uijmuj[,a] = w.eq[d,s,a]*(uijmuj[,a]) #+abs(min(uijmuj[,a])))

			}

			uijmuj[uijmuj<0]=0

			#Total (standardized) competition experienced by invader used to calculate invasion growth rates
			cc[d,t,,s] = apply(uijmuj,1,sum)/D[d,t,s]
			
			#Non-linear competitive variance and the spatial storage effect
			#These terms represent the perturbation terms from Snyder 2008 
			#NOTE: the /D[t,s] and /D[t,s]^2 terms are missing from the 
			#covariance and variance below because they are applied directly
			#to the standardization of cc above. 

			var_mu_Us[d,t,s] = var(cc[d,t,,s])

			cov_e_mu_Us[d,t,s] = cov(Frs_use[t,,s]/mean(Frs_use[t,,s])-1, cc[d,t,,s])

			Elam1[d,t,s]=l1[d,t,s]/D[d,t,s]*(1+(var_mu_Us[d,t,s])-
					(cov_e_mu_Us[d,t,s]))+sr[s]-1

			Elam2[d,t,s]=(l1[d,t,s]*(1/D[d,t,s]) +sr[s]-1)^2+2*(l1[d,t,s]*(1/D[d,t,s]) +sr[s]-1)*
				(var_mu_Us[d,t,s]-cov_e_mu_Us[d,t,s])/D[d,t,s]^4

			#The spatially implicit portion
			gr1.n[d,t,s] = exp(Elam1[d,t,s]-0.5*Elam2[d,t,s])
			
			#The fitness-density covariance 
			# tryCatch( {cov_lam_vc[t,s]=get.fd.cov(Frs_use[,s],sp.ids, 
			# 	nr[t,,],sr[s], alphas[s,], a_rr[s], kc)} , error=function(e){} )
									
			# cov_lam_vc2[t,s]=get.fd.cov2(Frs_use[,s],sp.ids, 
			# 	nr[t,,],sr[s], alphas, kd, kc)		
			#The full LDG
			#gr1[t,s] = gr1.n[d,t,s]+cov_lam_vc[d,t,s]
			gr1[d,t,s] = gr1.n[d,t,s]+covIMP.sim[d,t,s]


			print(s) #End loop through species
		}
		print(d)
		proc.time() - ptm
	}

		#=============================================================================
		#Stationary multispecies distribution through numerical integration
		#This part is very time-consuming. Implemented here so that only species with
		#positive LDG are allowed to equilibrate
		#=============================================================================
		for( sa in (1:nspp)){ 
			#Use the mean of reproduction: 
			Frs_use[t,,sa] = Frs[t,,sa]
		}
		#Frs[t,,1] = Frs[t,,1]*3


		print(ldg.sim[d,t,])
		s.index1= 1:nspp
		burns=500
		#Equilibrium populations of all spp. Columns represent space, rows are time. 
		nf.tmp=array(c(matrix(0.0,burns,np),matrix(0.0,burns,np)),dim=c(burns,np,nspp)) 
		for ( sa in 1:nspp) { 
			if (mean(ldg.sim [,t,sa]) >= 1){
			nf.tmp[1,,sa] = matrix(0.01,1,np)
			}
		}

		#Get the stationary multispecies distribution through numerical integration
		tcomp=0 #Condition for stationary joint distribution of community
		ts = 1
		comp_thresh = 2 #Sensitivity threshold. Theoretically, should be 0, but allow for numerical discrepency

		while (tcomp == 0 ){

			nf.tmp[ts+1,,] = pop_lg(Frs_use[t,,], nf.tmp[ts,,], sr,alphas, fkd, fkc,fkd.yes  )
			
			#Test for stationarity of community
			#if ( (mean(mean(nf[(ts+1),,]/nf[(ts),,]))-1) <= comp_thresh){ tcomp=1; is.final=T}
			if ( round(mean(mean(nf.tmp[(ts+1),,]/nf.tmp[(ts),,],na.rm=T),na.rm=T)-1,comp_thresh) == 0){ tcomp=1; is.final=T}

			#Stop if burn limit is reached
			if(ts >= burns) {tcomp=1; is.final=F}
			ts = ts+1

		}

		nf[t,,] = nf.tmp[(ts-1),,]

		y.full[[t]]= y 
		w.eq.full[[t]]=w.eq
		
	print(t) #End loop through time

	#=============================================================================
	#Data saving
	#=============================================================================
	#Save all of the meaningful invasion growth rate variables
	#in a file with an informative but awkwardly long file name. 
	#This includes information about the temporal and spatial extent,
	#and the details of the scenario type. 

	file.name = (paste(f.name1, "_waldmean.var",sep=""))
	save(file=file.name, "l1", "D", "var_mu_Us","cov_e_mu_Us",
	"cov_lam_vc", "cov_lam_vc2", "Elam1", "Elam2", "gr1.n", "gr1", "y.full", 
	"w.eq.full","ldg.sim", "ldgIMP.sim", "covIMP.sim", "Frs", "Frs_post", "nf","is.final","Frs_ci")

	
}

