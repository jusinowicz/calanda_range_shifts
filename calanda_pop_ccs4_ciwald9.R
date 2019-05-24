#=============================================================================
# R code for simulating species range shifts and calculating the strength of 
# coexistence between competitors. This code is based on the 
# coexistence_range_shifts code from an earlier repository, which means that 
# the underlying population dynamics are still modeled with a spatially explicit
# Leslie-Gower function. But the new code includes severa changes/additions 
# including:
#
# 1. Reproduction is a function of two processes that are both a function of 
# 	spatial location: probability of flowering, total flowers 
# 2. GAMs fit to abiotic covariates of elevation, i.e. mean temp, min temp, 
#    GDD, and soil moisture. 
# 3. Kriging of intrinsic ranges based on the GAM statistical fits.
# 4. Dispersal kernels that are based on the WALD approach (see code for references)
# 5. Analysis of persistence with dispersal kernels that vary in space 
# 	(i.e. are per-site dispersal kernels).
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
source("./wald_functions1.R")
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
#f.name1=c("calanda_ccs26_temp2_2017")
#f.name1=c("calanda_ccs26_tsmfull2_2017")
#f.name1=c("calanda_ccs85_temp2_2017")
f.name1=c("calanda_ccs85_tsmfull2_2017")
#=============================================================================
#Data: (see calanda2017.R)
#=============================================================================
load(file="calanda_allB2017_2.var")
rm(allB)

#With climate-change scenarios:
ccs_temp = read.csv("rcp85_2019_2100.csv")[2:3]
#ccs_temp = read.csv("rcp26_2019_2100.csv")[2:3]

ccs_soil = read.csv("rcp85_mrso_2019_2100.csv")[1:2]
#ccs_soil = read.csv("rcp26_mrso_2019_2100.csv")[1:2]

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

#Fix gs_min_tempfor allB2 in year 2015 by fitting a model to the other years
#by elevation, then projecting based on the mean temp at 1400m (the only reliable
#data point from 2015):

# elevations = unique(allB2$elevation)[2:6]
# gmt_tmp = subset(allB2, year != 2015)
# gmt_tmp2015 = subset(allB2, year == 2015)
# gmt_tmp2 = subset(allB2, year == 2015 & elevation == 1400 )

# #Fit a model and predict the new temperatures. 
# lm_gmt = gamm4(gs_min_temp~s(elevation,k=4), random = ~(1|year),data=gmt_tmp)
# pgmt =  predict.gam(lm_gmt$gam, newdata=data.frame(elevation=elevations),type="response")

# #This will serve as a lookup table for the new values of gs_mean_temp
# pgmt2015 = data.frame(elevation = elevations, gs_min_temp =pgmt - (pgmt[2]- unique(gmt_tmp2$gs_min_temp)))

# #Make a new object with a column of gs_mean_temp to replace the one in allB2
# temp.all=pgmt2015[match(interaction(gmt_tmp2015$elevation,"2015"), interaction(pgmt2015$elevation,"2015") ),]
# allB2[1:1200,]$gs_min_temp= temp.all$gs_min_temp

#Modify allB2 so that NA gs_min_temp counts match extremes: 
gmt_max1400 = max(subset(allB2, elevation == 1400 )$gs_min_temp, na.rm=T)

allB2[allB2$elevation==1400,]$gs_min_temp = gmt_max1400 
allB2[allB2$elevation==600,]$gs_min_temp = subset(allB2,allB2$elevation==1000)$gs_min_temp
allB2[allB2$elevation==2100,]$gs_min_temp = subset(allB2,allB2$elevation==2000)$gs_min_temp
allB2[allB2$elevation==2700,]$gs_min_temp = subset(allB2,allB2$elevation==2000)$gs_min_temp



# allB2[allB2$elevation>600 & allB2$elevation < 2100,]$gs_mean_temp = allB$gs_mean_temp

# allB2[allB2$elevation==1000 & allB2$year == 2015,]$gs_mean_temp = 
# 	allB2[allB2$elevation==1000 & allB2$year == 2015,]$gs_mean_temp - 3

#Modify allB2 so that NA inflor counts are 0: 
allB2[allB2$elevation==600,][2:4] = subset(allB2,allB2$elevation==1000)[,2:4] 
allB2[allB2$elevation==600,]$nr.inflor[is.na(allB2[allB2$elevation==600,]$nr.inflor)] = 0

allB2[allB2$elevation==2100,][2:4] = subset(allB2,allB2$elevation==1000)[,2:4] 
allB2[allB2$elevation==2100,]$nr.inflor[is.na(allB2[allB2$elevation==2100,]$nr.inflor)] = 0

allB2[allB2$elevation==2700,][2:4] = subset(allB2,allB2$elevation==1000)[,2:4] 
allB2[allB2$elevation==2700,]$nr.inflor[is.na(allB2[allB2$elevation==2700,]$nr.inflor)] = 0

#allB2 = allB2[allB2$elevation>600 & allB2$elevation < 2100,]

allB2$flyes[is.na(allB2$flyes)] = 0
head(allB2)

#allB2$year=as.factor(as.numeric(allB2$year))

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


#=============================================================================
#Which version of the intrinsic ranges? Which abiotic factors to use? 
#=============================================================================
# All of the possible variables for the model. Right now, this is 
# elevation, year, max_gdd, gs_min_temp,gs_max_temp, gs_mean_temp, soil_moist
variable_mat = colnames(allB2)[c(14,17,30,32,33,34,37)]

# Control which variables to use. Should include at least year and mean temp
#v_use= variable_mat[c(2,6)] #Temp only 
#v_use= variable_mat[c(2,7)] #Soil moisture only 	
#v_use= variable_mat[c(2,6,7)] #Temp + soil moisture
#v_use= variable_mat[c(2,4,6,7)] #Temp + soil moisture+min_temp
#v_use= variable_mat[c(2,3,6,7)] #Temp + soil moisture+max_gdd
v_use= variable_mat[c(2,3,4,6,7)] #Temp + soil moisture + max_gdd+ gs_mean_temp 
#v_use= variable_mat[c(2,1,3,4,6,7)] #elevation+Temp + soil moisture + max_gdd+ gs_mean_temp 

nvar = length(v_use)-1 # Exclude year from this, as it will always be a r.e.

#Speces-specific, variable-specific knots for smooth fits.
#These are based on analysis of underlying factor fits (see file: )
#Each species has 2 columns, one for flower probability model and one for total
#flowers model:

dgk = matrix(3,5,2); axk = dgk; hnk = dgk
#axk[2,1] = 5 #Soil moisture for AX
#hnk[2:3,1] = 5 #Soil moisture and min temp for HN
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

ngens=ttot #Number of generations for environmental change
iconfig=1 #NUmber of interattions of initial configuration of species

#Spatial scale: Assume that the gradient is from 500m to 2500m in units of 10m
elevations=c(600,1000,1400,1600,1800,2000,2100,2700)
el_real=c(1000,1400,1600,1800,2000)
el1=0
el2=3500

#Make some internal spatial variables:
# xx1= matrix(seq(el1,el2,1))
# xx2= matrix(xx1^2)
# #yx= matrix(2017,length(xx1),1)
# yx= matrix(3,length(xx1),1)

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
#temp_gam = gam( gs_mean_temp~ s(year,bs="re")+s(elevation,k=3),data=allB2)
#temp_gam1 = gam( gs_mean_temp~ s(year,bs="re")+I(elevation^2),data=allB2)
#temp_gam1 = gam( gs_mean_temp~ 1+s(year,bs="re")+I(elevation^2),data=subset(allB2, !is.na(allB2$gs_mean_temp)))
#temp_gam = gam( gs_mean_temp~ elevation+ I(elevation^2),data=allB2)
#temp_gam = gamm4( gs_mean_temp~ elevation, random = ~(1|year), data=subset(allB2, !is.na(allB2$gs_mean_temp)) )
#temp_gam3 = lmer( gs_mean_temp~ I(elevation^2)+(1|year), data=subset(allB2, !is.na(allB2$gs_mean_temp)) )
temp_gam = gamm4(gs_mean_temp~s(elevation,k=4), random = ~(1|year),data=subset(allB2, !is.na(allB2$gs_mean_temp)))

#temp_gam2 = gamm4( gs_mean_temp~ elevation, random = ~(1|year), data=subset(allB2, !is.na(allB2$gs_mean_temp)) )

#soil moisture 
#sm_gam = gam( soil_moist~ s(year,bs="re")+s(elevation,k=3),data=allB2)
#sm_gam = gam( soil_moist~ s(year,bs="re")+elevation+I(elevation^2),data=allB2)
#sm_gam = gamm4( soil_moist~ elevation+I(elevation^2),random = ~(1|year), data=allB2)
sm_gam = gamm4(soil_moist~s(elevation,k=4), random = ~(1|year),data=subset(allB2, !is.na(allB2$gs_mean_temp)))


#PET
#sp_gam = gam( soil_pet~ s(year,bs="re")+s(elevation,k=3),data=allB2)
#sp_gam = gam( soil_pet~ s(year,bs="re")+I(elevation^2),data=allB2)
sp_gam = gamm4( soil_pet~ elevation+I(elevation^2),random = ~(1|year), data=allB2)


#Moisture deficit
#sd_gam = gam( soil_def~ s(year,bs="re")+s(elevation,k=3),data=allB2)
#sd_gam = gam( soil_def~ s(year,bs="re")+I(elevation^2),data=allB2)
sd_gam = gamm4( soil_def~ elevation+I(elevation^2),random = ~(1|year), data=allB2)

#annual mean light
#ml_gam = gam( an_mean_light~ s(year,bs="re")+s(elevation,k=3),data=allB2)
#ml_gam = gam( an_mean_light~ s(year,bs="re")+I(elevation^2),data=allB2)
ml_gam = gamm4( an_mean_light~ I(elevation^2), random = ~(1|year), data=subset(allB2, !is.na(allB2$gs_mean_temp)) )

#GDD
#gdd_gam = gam( max_gdd~ s(year,bs="re")+s(elevation,k=3),data=allB2)
#gdd_gam = gam( max_gdd~ s(year,bs="re")+I(elevation^2),data=allB2)
#gdd_gam = gamm4( max_gdd~ I(elevation^2), random = ~(1|year), data=subset(allB2, !is.na(allB2$gs_mean_temp)) )
#gdd_gam = gamm4(max_gdd~s(elevation,k=4), random = ~(1|year),data=subset(allB2, !is.na(allB2$gs_mean_temp)))
gdd_gam = gamm4(max_gdd~s(elevation,k=5), random = ~(1|year),data=allB2)

#Minimum growing season temp
#gmt_gam = gam( gs_min_temp~ s(year,bs="re")+s(elevation,k=3),data=allB2)
#gmt_gam = gam( an_min_temp~ s(year,bs="re")+I(elevation^2),data=allB2)
#gmt_gam = gamm4( gs_min_temp~ I(elevation^2), random = ~(1|year), data=subset(allB2, !is.na(allB2$gs_mean_temp)) )
gmt_gam = gamm4(gs_min_temp~s(elevation,k=5), random = ~(1|year),data=allB2)

#=============================================================================
# #Equivalent but wrong: 
# #The mean temp vector: 
# xx_mt1 = as.vector(predict.gam(temp_gam1, newdata=xx,type= "response"))
# #xx_mt3 without re.form gives the same fit as xx_mt1 with a specific year specified
# xx_mt3 = as.vector(predict(temp_gam3, newdata=xx,type= "response"))

# #Different and wrong:
# #xx_mt1 with the random effect fit but then excluded is not correct (it gives the lowest prediction)
#xx_mt1 = as.vector(predict.gam(temp_gam1, exclude ="s(year)", newdata=xx,type= "response"))

#Equivalent and correct:
xx_mt = as.vector(predict.gam(temp_gam$gam, newdata=xx,type= "response"))
#xx_mt3 with re.form=0 gives the same fit as just the gam part of xx_mt2 (temp_gam2$gam)
#xx_mt3 = as.vector(predict(temp_gam3,re.form= NA, newdata=xx,type= "response"))

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
#xx_gmt = predict( lm(xx_gmt~xx$elevation)) #Uncomment to linearize xx_gmt

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
for(sp in 1:length(spp)){
	allsp = subset(allB2, Sp == spp[sp])
	sr_dat=subset(allsp, bg =='B')
	sr_dat = subset(sr_dat, elevation>=el_real[1] & elevation<=el_real[5])
	#sr[sp] = mean(allB2$survival)
	sr_gam = gam(survival~s(elevation,k=kn[sp])+s(year,bs="re"),family=binomial(link='logit'), data=sr_dat)
	sr_tmp=as.vector(predict.gam(sr_gam, newdata=xx,type= "response"))
	#Remove 0s
	sr_tmp[sr_tmp<0] = 0
	sr_krigB[,sp] = sr_tmp
	sr_probB[[sp]] = sr_gam
}

#Plot kriged values
par(mfrow=c(3,1))
for(sp in 1:3) {
	plot(sr_krigB[,sp],ylim=c(0,1))
	#points(sr_krig[,sp],col="red")
}


#=============================================================================
#Intraspecific competition
#=============================================================================
#See calanda_pop_krigALL_wald.R for the full statistical exploration
#=============================================================================
#Use NLS to fit the functional relationship (leslie-gower?). Try it per elevation
#too, to see if it is significant. If it is, Krige it. If not...? 
#This version follows the matching of a the per-plot background individual
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
		# Cr_dat$nr.inflor = Cr_dat$nr.inflor*Cr_dat$stalks+Cr_dat$nr.stms	
		# rr_sub$nr.inflor = rr_sub$nr.inflor*rr_sub$stalks+rr_sub$nr.stms	
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
	plot(NA, NA, ylim=c(0, 30),xlim = c(0, 20)) 
	for(ee in 1:length(el_real)){
		Cr_dat_e = subset(Cr_dat, elevation == el_real[ee])
		rr_sub_e = subset(rr_sub, elevation == el_real[ee])
		nsamp = dim(rr_sub_e)[1]

		#Bootstrap this by randomly pairing bare and competitive backgrounds
		#This would maintain within-site pairing
		crdat_e= data.frame(cbind(Cr_dat_e$nr.inflor[ceiling(runif(bsamp)*nsamp)], 
			rr_sub_e$nr.inflor[ceiling(runif(bsamp)*nsamp)]))
		colnames(crdat_e) = c("lrr","Rr")
		points(crdat_e$Rr ,crdat_e$lrr )

		#Fit with nls
		tryCatch( {Cr_nls_e = nls(lrr~ Rr/(1+cr),data=crdat_e, start=list(cr=0.5),control = list(maxiter = 500))
			Crr_nls_el[ee, sp] = Cr_nls_e$m$getPars()}, error = function(e){})
		
	}
	
	#Now use the per elevation data to Krige other values
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
#This version follows the matching of a the per-plot background individual
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
			#Ci_dat$nr.inflor = Ci_dat$nr.inflor*Ci_dat$stalks+Ci_dat$nr.stms	
			#rr_sub$nr.inflor = rr_sub$nr.inflor*rr_sub$stalks+rr_sub$nr.stms	
			Ci_dat$nr.inflor = Ci_dat$nr.inflor*Ci_dat$stalks	
			rr_sub$nr.inflor = rr_sub$nr.inflor*rr_sub$stalks	
		}

		#Bootstrap this by randomly pairing bare and competitive backgrounds
		#This would ignore elevation as a treatment
		cidat= data.frame(cbind(Ci_dat$nr.inflor[ceiling(runif(bsamp)*nsamp)],
		 rr_sub$nr.inflor[ceiling(runif(bsamp)*nsamp)]))
		
		colnames(cidat) = c("lir","Ri")
		
		#Fit with nls
		#Fit with nls
		tryCatch( {		Ci_nls = nls(lir~ Ri/(1+ci),data=cidat, start=list(ci=Crr_nls[sp]),control = list(maxiter = 500))
				Cir_nls[sp, sb] = Ci_nls$m$getPars()
				}, error = function(e){})
		#Ci_nls = nls(lir~ Ri/(1+ci),data=cidat, start=list(ci=Crr_nls[sp]),control = list(maxiter = 500))
		#print(summary(Ci_nls))
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
			points(cidat_e$Ri ,cidat_e$lir )
			
			#Fit with nls
			tryCatch( {Ci_nls_e = nls(lir~ Ri/(1+ci),data=cidat_e, start=list(ci=Crr_nls[sp]),control = list(maxiter = 500))
				Cir_nls_el[ee, sb] = Ci_nls_e$m$getPars()
				}, error = function(e){})
			
		}
		mean(Ci_dat$nr.inflor,na.rm=T)
		#Now use the per elevation data to Krige other values
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
#and air*Nr, then taking the ration air/arr, setting arr=1
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
#dispersal
#=============================================================================
#Since all species' dispersal kernels are so local, no dispersal happens 
#across sites. 
#=============================================================================
#Krig a distribution of heights as a function of background type
#Use only those that have flowered

#=============================================================================
#Invasion growth rates
#=============================================================================
#=============================================================================
#Tunable Species parameters
#=============================================================================
#Calculate spatiotemporal Fr (fundamental niche)
###Maximum reproduction rates, corresponding to spatio-temporal ideal conditions

###Survival
#sr=surv_spp
#sr=sr_krigB
#Or take the average (since there's not that much variation), but load the
#other function file:
sr=colMeans(sr_krigB[800:2000,])
sr[3]=sr[3]+0.05

###Competition coeffcients
#3B is calculating a per-site alpha, then averaging 
# diag(arat3B)=c(1,1,1)
# alphas=arat3B
# alphas[2,1] = arat3[2,1]

#3 is calculating an alpha that has been calculated regardless
#of elevation
diag(arat3)=c(1,1,1)
alphas=arat3


###Competition distance -- These are unknown, and chosen just to be limited or not
#b_rr=1/(100*np) #Essentially global
b_rr=c(.1,.1,.1) #10 cm plots

###Dispersal distance -- These are calculated using the WALD model --see 
#wald_model1.R

#a_rr= c(1/(100*np),1/(100*np),1/(100*np) #essentially global
#a_rr=c(1/50,1/50,1/50)
#Set the means (for an exponential) to what they are in the WALD kernel below 

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
#-- this section includes data exploration across elevations and an effect of 
# elevation on dispersal via the height
#=============================================================================
#Krig a distribution of heights as a function of background type
#Use only those that have flowered

#Height in the bare soil (intrinsic)
height_krig = matrix(0,dim(xx)[1], length(spp))
height_gams=NULL
kn = c(5,5,5)
for(sp in 1:length(spp)){
	allsp = subset(allB2, Sp == spp[sp])
	ht_dat_bg =subset(allsp, bg =='B')
	ht_dat_bg =subset(ht_dat_bg, elevation>=el_real[1] & elevation<=el_real[5])
	ht_dat =subset(ht_dat_bg, flyes ==1)
	ht_gam = gam(stalk.height~s(year,bs="re")+s(elevation,k=kn[sp]), data=ht_dat)
	
	#ht_dat =subset(allsp, flyes ==1)	
	#ht_gam = gam(stalk.height~bg+s(year,bs="re")+s(elevation,k=kn[sp]), data=ht_dat)
	
	ht_tmp=as.vector(predict.gam(ht_gam, newdata=xx,type= "response"))
	#Remove 0s
	ht_tmp[ht_tmp<0] = 0
	height_krig[,sp] = ht_tmp
	height_gams[[sp]]=ht_gam
}



#Height as a function of background type -- The statistical models don't 
#seem to support this as a signficant effect (whether because of sufficient data
#is unclear)

####Dispersal kernels and their Fourier transforms 

#Get the WALD kernel based on field measurements
#Find the mean windspeed: WARNING, do not do this every time if these files are
#very large (i.e. >500 MB)!
#write.csv(stmp,file="nswind2017s.csv" ) #File with only first 14 columns
# site_files=c("/home/jacob/labshare/Jacob/wind_data/2017/Aurella/arwind2017.csv",
# 			"/home/jacob/labshare/Jacob/wind_data/2017/Neselboden/nswind2017.csv",
# 			"/home/jacob/labshare/Jacob/wind_data/2017/Barenmos/bmwind2017.csv",
# 			"/home/jacob/labshare/Jacob/wind_data/2017/Neusass/newind2017.csv",
# 			"/home/jacob/labshare/Jacob/wind_data/2017/Calanda/cawind2017.csv")

#u_mean = get.u_mean(site_files, headers=FALSE, col_name = )

#Site-specific means from get.u_mean
u_mean_sites=c(1.157264,1.253033,1.070176, 1.045009, 2.845403)
u_var_sites = c(0.5265219,0.6048379,0.4641298,0.5813984,5.014292 )

u_mean=mean(u_mean_sites)

#u_mean = 2.85 #Mean windspeed above canopy NOTE: Old value
Vt=c(2.6,3,3.3) #terminal velocity DG maybe 2.6, AX maybe 3, HN maybe 3.3

#Model parameters for the Beta function: Based on Su et al. 2001
a1=1.05
a2=2
a3=0.1

kd=array(c(matrix(0,np,np),matrix(0,np,np)),dim=c(np,np,nspp)) 
fkd=kd

wald.list = NULL
for( s in 1:nspp){ 

	#This version just uses the average height
	height_k2 = matrix(c(colMeans(height_krig)),np,nspp, byrow=T)


	#The variable wald.list[[s]] will contain a list with: 
	#site.kernels[[1]][[1]]		the dispersal kernel
	#site.kernels[[1]][[2]]		the mean dispersal distance
	#site.kernels[[1]][[3]]		the 90th percentile distance 

	#Pick version of height for kernel:
	#Variable height
	#wald.list[[s]] = get.WALD.kernel(u_mean, xx0, height_krig[,s]/100, Vt[s],a1,a2,a3 )[[1]][[1]]

	#Mean height
	wald.list[[s]] = get.WALD.kernel(u_mean, xx0, height_k2[,s]/100, Vt[s],a1,a2,a3 )[[1]][[1]]

	kd[2:ceiling(np/2),,s] = wald.list[[s]][floor(np/2):1,]
	kd[ceiling(np/2):(np-1),,s] = wald.list[[s]][1:floor(np/2),]


	#Normalize again. 
	#Because get.WALD.kernel normalizes the right half of the kernel to 1
	kd[,,s]=kd[,,s]/matrix(colSums(kd[,,s],na.rm=T),np,np,byrow=T)
	kd[,,s][is.na(kd[,,s])] = 0

	#FFT
	fkd[,,s]=mvfft(kd[,,s])#/(np+1)
	fkd.yes=TRUE #Pass the transformed kd to functions later for speed

}

#Need these for testing the effect of each site on LGR later on
kd.n = kd[-np,,]
fkd.n = kd.n
for (sa in 1:nspp) {fkd.n[,,sa]=mvfft(kd.n[,,sa]) }

#Set the means (for an exponential) to what they are in the WALD kernel: 
a_rr = c(1/wald.list[[1]][[2]],1/wald.list[[2]][[2]],1/wald.list[[3]][[2]])

#Exponential kernel
# kd=matrix(0,np,nspp)
# fkd=matrix(0,np,nspp)

# for( s in 1:nspp){ 
# 	kd[,s] = a_rr[s]/2*exp(-a_rr[s]*abs(xx0))
# 	kd[,s]=kd[,s]/(sum(kd[,s]))
# 	fkd[,s]=fft(kd[,s])#/(np+1)
# 	fkd.yes = TRUE #Pass the transformed kd to functions later for speed

# }


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
Frs_mlp = array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp)) 
Frs_ci = vector("list", ngenst) 
#The maximum of species' fitted reproduction, for scaling
frs_max = matrix(0,ngenst,nspp)

#Realized ranges in full commmunity
nf = array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp)) 

#Equilibrium populations of residents. Columns represent space, rows are time. 
nr = array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp)) 


#Simulation IGR
ldg.sim = matrix(0,ngenst, nspp)
ldgIMP.sim = matrix(0,ngenst, nspp)
covIMP.sim = matrix(0,ngenst, nspp)

#Components of IGR
l1=matrix(0,ngenst,nspp)
D = matrix(0,ngenst,nspp)
var_mu_Us=matrix(0,ngenst,nspp) 
cov_e_mu_Us=matrix(0,ngenst,nspp)
cov_lam_vc=matrix(0,ngenst,nspp)
cov_lam_vc2=matrix(0,ngenst,nspp)
Elam1=matrix(0,ngenst,nspp) 
Elam2=matrix(0,ngenst,nspp)
gr1.n=matrix(0,ngenst,nspp)
gr1=matrix(0,ngenst,nspp)

#Components kept as lists
y.full=vector(ttot, mode="list")
w.eq.full=vector(ttot, mode="list")

#Site impacts on the LGR
sim.impacts=array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp))
sim.impactsIMP=array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp))
covIMP.impacts=array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp))
lded=array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp)) 

env_analogue_all = vector(ttot, mode="list")
env_var = vector(ttot, mode="list")
#=============================================================================
#Outermost loop: treat each species as invader
#=============================================================================

#Keep track of total interspecific competition of species when it is invader
cc=array(c(matrix(0,ngenst,np),matrix(0,ngenst,np)),dim=c(ngenst,np,nspp)) 

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

	#env_distance = get_env_distance(xx,xx_new,env.ind)
	env_analogue = get_analogue(xx,xx_new,env.ind)
	# var_dist = env_analogue[[4]]
	# colnames(var_dist)=colnames(xx)[env.ind]
	# an_dist = sqrt( (env_analogue[[3]] - xx_new$elevation)^2)

	#=============================================================================
	#FLOWERS
	#=============================================================================

	#=============================================================================
	#intrinsic growth
	#=============================================================================
	#Intrinsic probability of flowering

	flower_probK=matrix(0,dim(xx)[1], length(spp))
	flower_prob_ci=matrix(0,length(spp),2)
	flower_prob=NULL
	flower_prob_post=NULL
	flower_prob_act = matrix(0,length(el_real),3)

	for(sp in 1:length(spp)){
		allsp = subset(allB2, Sp == spp[sp])
		rr_dat =subset(allsp, bg =='B')
		rr_dat = subset (rr_dat, elevation <= el_real[5] & elevation >= el_real[1])
		
		if(sp ==1 ) { 
			#Ci_dat$nr.inflor = Ci_dat$nr.inflor*Ci_dat$stalks+Ci_dat$nr.stms	
			#rr_sub$nr.inflor = rr_sub$nr.inflor*rr_sub$stalks+rr_sub$nr.stms	
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
			#factors2[f] = paste("s(",v_use[(f+1)], ", k = ", kn[sp], ",  bs=\"cr\")", sep="")
			factors2[f] = paste("s(",v_use[(f+1)], ", k = ", var_spp_knots[[sp]][f,1]+nk, ",  bs=\"cr\")", sep="")
			#factors2[f] = paste("s(",v_use[(f+1)], ", k = ", kn[sp]+nk, ")", sep="")

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
		pnk = nk/2
		kc=0
		for( b in 1:nvar){
			if(b>1) {kc = sum(kn[1:(b-1)])}
			bpos = (b-1)*(nk)+kc+pnk+1
			kpos = kc+1 #Intercept removed from rr_gam
			beta[ bpos:(bpos+kn[b]-1)] = coef(rr_gam) [kpos:(kpos+kn[b]-1)]
		} 

		### prediction matrix
		Xp = matrix(0, dim(xx_new)[1], (kuse+nvar*nk))
		#Xp = matrix(0, dim(dat)[1], nvar*(kuse+nk))
		kc=0
		for( b in 1:nvar){
			if(b>1) {kc = sum(kn[1:(b-1)])}
			bpos = (b-1)*(nk)+kc+pnk
			Xp_tmp = PredictMat(sm[[b]], xx_new)
			#Xp_tmp = PredictMat(sm[[b]], dat)
			Xp[,bpos:(bpos+kn[b]+nk-1)] = Xp_tmp
		}
		

		### the predicted smooth, i.e. species' distribution
		#flower_tmp = Xp%*%beta+coef(rr_gam)[1] #These are not on the response scale
		flower_tmp = Xp%*%t(beta)+ cp_y[1] #These are not on the response scale
		flower_tmp = exp(flower_tmp )/(1+exp(flower_tmp ) ) #Response scale for logistic regression
		
		#The unconstrained model: 
		bu_tmp=as.vector(predict(b.u,newdata=xx_new,type= "response",exclude = "s(year)") )
		#flower_tmp=as.vector(predict.gam(flower_ptmp$gam, newdata=xx_new,type= "response"))
		#flower_tmp=as.vector(predict.gam(flower_ptmp, newdata=xx_new,type= "response"))

		#Get confidence intervals: 
		#Xp =predict.gam(flower_ptmp$gam, newdata=xx_new,type= "lpmatrix")
		#br1=mvrnorm(n=10000, coef(flower_ptmp$gam), flower_ptmp$gam$Vp)
		#This produces the mean model: 
		#mean.model= rowMeans(exp(Xp%*%t(br1)))
		#But to get the CIs:
		#mean.model= colMeans(exp(Xp%*%t(br1)))
		#q1=quantile(mean.model,c(.025,.975))
		
		#Remove 0s
		#flower_tmp[flower_tmp<0] = 0
		flower_probK[,sp] = flower_tmp
		#flower_probK[,sp] = bu_tmp
		flower_prob[[sp]] = rr_gam
		#flower_prob_post[[sp]] = exp(Xp%*%t(br1))
	
		#flower_prob_ci [sp,] = q1 
		#plot(rr_dat$gs_mean_temp, residuals(flower_prob[[sp]]$mer))

	}


	#Krig a distribution of Rs (intrinsic growth) from flower counts of plants 
	#grown alone, per site 
	#Use only those that have flowered
	#The best approach here seems to be to use the probability by elevation, 
	#then each R is actually a constant.

	rr_krig = matrix(0,dim(xx)[1], length(spp))
	rr_ci=matrix(0,length(spp),2)
	rr_gams=NULL
	rr_post = NULL

	flower_act =matrix(0,length(el_real),3)

	for(sp in 1:length(spp)){
		allsp = subset(allB2, Sp == spp[sp])
		rr_dat_bg =subset(allsp, bg =='B')
		rr_dat =subset(rr_dat_bg, flyes ==1)

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
			#factors2[f] = paste("s(",v_use[(f+1)], ", k = ", kn[sp], ",  bs=\"cr\")", sep="")
			factors2[f] = paste("s(",v_use[(f+1)], ", k = ",var_spp_knots[[sp]][f,2]+nk, ",  bs=\"cr\")", sep="")
			#factors2[f] = paste("s(",v_use[(f+1)], ", k = ", kn[sp]+nk, ")", sep="")

		} 

		## In order to find the minimum intercept (which is not zero with the logit link function): 
		## First fit a model to the data: 
		b.u = gam(as.formula(paste("nr.inflor~", paste( paste(factors2, collapse="+"),"+s(year, bs=\"re\")"))), data=dat, knots = knots_sm )
		bu_tmp = predict(b.u, exclude = "s(year)") #We need the minimum of this 
		cp_y = c(matrix(0,length(cp_x),1)) #y Values to force the intercept to
		off = dat$nr.inflor*0 + cp_y[1] ## offset term to force curve through a point

		rr_gam= gam(as.formula(paste("nr.inflor~", paste( paste(factors, collapse="+"),"+offset(off)+s(year, bs=\"re\")-1"))), 
			 paraPen=para_pen, data=dat,method = "REML")
		# rr_gam= gam(as.formula(paste("nr.inflor~", paste( paste(factors, collapse="+"),"+offset(off)+s(year, bs=\"re\")-1"))), 
		# 	paraPen=para_pen, data=dat)
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
		pnk = nk/2
		kc=0
		for( b in 1:nvar){
			if(b>1) {kc = sum(kn[1:(b-1)])}
			bpos = (b-1)*(nk)+kc+pnk+1
			kpos = kc+1 #Intercept removed from rr_gam
			beta[ bpos:(bpos+kn[b]-1)] = coef(rr_gam) [kpos:(kpos+kn[b]-1)]
		} 

		### prediction matrix
		Xp = matrix(0, dim(xx_new)[1], (kuse+nvar*nk))
		#Xp = matrix(0, dim(dat)[1], nvar*(kuse+nk))
		kc=0
		for( b in 1:nvar){
			if(b>1) {kc = sum(kn[1:(b-1)])}
			bpos = (b-1)*(nk)+kc+pnk
			Xp_tmp = PredictMat(sm[[b]], xx_new)
			#Xp_tmp = PredictMat(sm[[b]], dat)
			Xp[,bpos:(bpos+kn[b]+nk-1)] = Xp_tmp
		}
		
		### the predicted smooth, i.e. species' distribution
		#rr_tmp = Xp%*%beta+coef(rr_gam)[1] #These are not on the response scale
		rr_tmp = Xp%*%t(beta) + cp_y[1]

		#The unconstrained model: 
		bu_tmp=as.vector(predict(b.u,newdata=xx_new,type= "response",exclude = "s(year)") )
	
		#Get confidence intervals: 
		# Xp =predict.gam(rr_gam$gam, newdata=xx_new,type= "lpmatrix")
		# br1=mvrnorm(n=10000, coef(rr_gam$gam), rr_gam$gam$Vp)
		#This produces the mean model: 
		#mean.model= rowMeans(exp(Xp%*%t(br1)))
		#But to get the CIs:
		#mean.model= colMeans((Xp%*%t(br1)))
		#q1=quantile(mean.model,c(.025,.975))
		

		#Remove 0s
		rr_tmp[rr_tmp<0] = 0
		rr_krig[,sp] = rr_tmp
		#rr_krig[,sp] = bu_tmp

		rr_gams[[sp]]=rr_gam
		#rr_post[[sp]] = (Xp%*%t(br1))
		#rr_ci [sp,] = q1 

		#plot(rr_dat[[sp]]gs_mean_temp, residuals(rr_gams[[sp]]))

	}

	####Intrinsic ranges
	# Frs=array(c(rr_krig*flower_probK),dim=c(ngenst,np,nspp)) 
	# #Frs[,,1] = Frs[,,1]*12
	# Frs[,,1] = Frs[,,1]*8
	####Intrinsic ranges
	Frs[t,,] =matrix(rr_krig*flower_probK,np,nspp)

	#Rescaled from data:
	max(flower_act[,sa],na.rm=T)
	for( sa in 1:nspp){

		frs_max[t,sa] = max(Frs[t,,sa],na.rm=T)
		#Make the scaling factor a percentage of the non-shifted scaling
		max_scale =  frs_max[t,sa]/frs_max[1,sa]*max(flower_act[,sa],na.rm=T)
		Frs[t,,sa]=Frs[t,,sa]*max_scale/max(Frs[t,,sa],na.rm=T)
	}

	Frs[t,,1] = Frs[t,,1]*3


	fig.name = paste("calanda_ranges",paste(v_use[-1],collapse=""),".pdf",sep="")
	pdf(file=fig.name, height=11, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)
	col_use = c("black","red","blue")
	par(mfrow=c(3,1))
	plot(Frs[t,,2],t="l")
	for (tt in 1:3){
		lines(Frs[t,,tt],col=col_use[tt])	
	}

	plot(flower_probK[,2],t="l",ylim=c(0,1))
	for (tt in 1:3){
		lines(flower_probK[,tt],col=col_use[tt])	
	}


	plot(rr_krig[,2],t="l")
	for (tt in 1:3){
		lines(rr_krig[,tt],col=col_use[tt])	
	}


	dev.off()
	
	# #Frs[t,,1] = Frs[t,,1]*8
	# for (sp in 1:length(spp)) {
	# 	Frs_q1 = apply(flower_prob_post[[sp]]*rr_post[[sp]], 1, quantile, c(0.025,0.975))
	# 	Frs_ci[[t]][[sp]] = Frs_q1
	# 	#Frs_mlp[t,,sp] = rowMeans(flower_prob_post[[sp]]*rr_post[[sp]])
	# }
	
	# #Frs[,,1] =Frs[,,1] - matrix(apply(matrix(Frs[,,1]),2,min),np, ngenst, byrow=T)


	# env_var[[t]] = xx_new[,env.ind] }

	#=============================================================================
	#=============================================================================

	#=============================================================================
	#Invader and resident stationary distributions
	#=============================================================================
	y = matrix(0,nspp, nspp-1)
	w.eq = matrix(0,nspp, nspp-1)

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

		#Get the equilibrium of the resident community when invader is absent
		nr[t,,(s.index)] = get.res.eq(Frs[t,,],s.index,sr,alphas, fkd,fkc, fkd.yes,fast=TRUE, burns = 5 )
	
		#Get the low-density equilibrium density of the invader against the resident community
		nr[t,,s] = get.inv.ldeq(Frs[t,,], nr[t,,], s, sr,alphas, fkd, fkc,fkd.yes,fast=TRUE, burns = 5 )
		lded [t,,s] = nr[t,,s]
	
		#=============================================================================
		# Low-density growth rates -- using simulation data, all 3 spp
		#=============================================================================
		#Simple numerical LGR: use the invader low-density equilibrium density, 
		#then let it invade against the community: 

		sp.ids =c(s, s.index) #IDs of resident, and invaders (this way uses all species)

		inv.one =pop_lg(Frs[t,,],nr[t,,], sr, alphas, fkd,fkc, fkd.yes )[,s]
		ldg.sim[t,s]=mean(inv.one)/mean(nr[t,,s])

		# #Numerically calculate the contribution of each site to the LGR
		# # First, resize competition kernel
		# kc.n=matrix(0,np-1,nspp)
		# fkc.n=matrix(0,np-1,nspp)
		# xx0B=matrix(seq((-(ns-1)/2),((ns-1)/2)))
		# for( sa in 1:nspp){ 
		# 	kc.n[,sa] = b_rr[sa]/2*exp(-b_rr[sa]*abs(xx0B))
		# 	kc.n[,sa]=kc.n[,sa]/(sum(kc.n[,sa]))
		# 	fkc.n[,sa]=fft(kc.n[,sa])#/(np+1)

		# }
			
		# #Do the calculations: 
		# for (t in 1:ngenst) { 
		# 	sim.impacts[t,,s] = site_impact_sim(Frs[t,,],nr[t,,], sr, alphas, fkd.n,fkc.n, fkd.yes,ngenst,sp.ids )
		# }

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
		nr2[,,s] = matrix(inv.id,ngenst,np)

		inv.oneIMP =pop_lg(Frs[t,,], nr2[t,,], sr, alphas, kdg,kc )[,s]
		ldgIMP.sim[t,s]=mean(inv.oneIMP)/mean(nr2[t,,s])

		covIMP.sim[t,s] = ldg.sim[t,s]-ldgIMP.sim[t,s]

		#Numerically calculate the contribution of each site to the LGR
		# First, resize competition kernel
		# kdg.n=matrix(0,np-1,nspp)
		# fkdg.n=matrix(0,np-1,nspp)
		# for( sa in 1:nspp){ 
		# 	kdg.n[,sa] = 1/(np-1)
		# 	fkdg.n[,sa]=fft(kdg.n[,sa])#/(np+1)
		# }
		
		# #Do the calculations: 
		# for (t in 1:ngenst) { 
		
		# 	sim.impactsIMP[t,,s] = site_impact_sim(Frs[t,,],nr2[t,,], sr, alphas, kdg,kc.n,fkd.yes=FALSE,ngenst,sp.ids )
		
		# }


		# covIMP.impacts[,,s] = sim.impacts[,,s]-sim.impactsIMP[,,s]

		#Calculate the LDG for species vs. community by component

		#=============================================================================
		#Calculate standard spatial terms
		#=============================================================================
		nrs= nr[t,,s.index] #Residents 
		
		#isd.tmp = get.isd(Frs[t,,s], nrs,s, alphas, sr[s], kc, kd[,s])

		y[s,] = colMeans(nrs,na.rm=T) #Means of resident densities 
		w.eq[s,] = alphas[s,s.index]*y[s,] #Weights for the LDG

		l1[t,s]=mean(Frs[t,,s])
		D[t,s]= 1+sum(w.eq[s,])


		#Calculate the standardized competition from residents
		muj=nrs/(matrix(y[s,],np,( nspp-1),byrow=T))-1
		uijmuj = matrix(0,np,nspp-1)
		for(a in 1:(nspp-1)) {
				ns = s.index[a]
				uijmuj[,a] = convolve(muj[,a],kc[,ns])
				uijmuj[,a]= c(uijmuj[ceiling(np/2):(np),a], uijmuj[1:floor(np/2),a] )
				uijmuj[,a] = w.eq[s,a]*(uijmuj[,a]) #+abs(min(uijmuj[,a])))

		}

		uijmuj[uijmuj<0]=0

		#Total (standardized) competition experienced by invader used to calculate invasion growth rates
		cc[t,,s] = apply(uijmuj,1,sum)/D[t,s]
		
		#Non-linear competitive variance and the spatial storage effect
		#These terms represent the perturbation terms from Snyder 2008 
		#NOTE: the /D[t,s] and /D[t,s]^2 terms are missing from the 
		#covariance and variance below because they are applied directly
		#to the standardization of cc above. 

		var_mu_Us[t,s] = var(cc[t,,s])

		cov_e_mu_Us[t,s] = cov(Frs[t,,s]/mean(Frs[t,,s])-1, cc[t,,s])

		Elam1[t,s]=l1[t,s]/D[t,s]*(1+(var_mu_Us[t,s])-
				(cov_e_mu_Us[t,s]))+sr[s]-1

		Elam2[t,s]=(l1[t,s]*(1/D[t,s]) +sr[s]-1)^2+2*(l1[t,s]*(1/D[t,s]) +sr[s]-1)*
			(var_mu_Us[t,s]-cov_e_mu_Us[t,s])/D[t,s]^4

		#The spatially implicit portion
		gr1.n[t,s] = exp(Elam1[t,s]-0.5*Elam2[t,s])
		
		#The fitness-density covariance 
		# tryCatch( {cov_lam_vc[t,s]=get.fd.cov(Frs[t,,s],sp.ids, 
		# 	nr[t,,],sr[s], alphas[s,], a_rr[s], kc)} , error=function(e){} )
								
		# cov_lam_vc2[t,s]=get.fd.cov2(Frs[t,,s],sp.ids, 
		# 	nr[t,,],sr[s], alphas, kd, kc)		
		#The full LDG
		#gr1[t,s] = gr1.n[t,s]+cov_lam_vc[t,s]
		gr1[t,s] = gr1.n[t,s]+covIMP.sim[t,s]


	print(s) #End loop through species
	}


	#=============================================================================
	#Stationary multispecies distribution through numerical integration
	#This part is very time-consuming. Implemented here so that only species with
	#positive LDG are allowed to equilibrate
	#=============================================================================


	s.index1= 1:nspp
	burns=500
	#Equilibrium populations of all spp. Columns represent space, rows are time. 
	nf.tmp=array(c(matrix(0.0,burns,np),matrix(0.0,burns,np)),dim=c(burns,np,nspp)) 
	for ( sa in 1:nspp) { 
		if (ldg.sim [t,sa] >= 1){
		nf.tmp[1,,sa] = matrix(0.01,1,np)
		}
	}

	#Get the stationary multispecies distribution through numerical integration
	tcomp=0 #Condition for stationary joint distribution of community
	ts = 1
	comp_thresh = 2 #Sensitivity threshold. Theoretically, should be 0, but allow for numerical discrepency

	while (tcomp == 0 ){

		nf.tmp[ts+1,,] = pop_lg(Frs[t,,], nf.tmp[ts,,], sr,alphas, fkd, fkc,fkd.yes  )
		
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
	#env_analogue_all[[t]] = env_analogue 
	env_var[[t]] = xx_new[,env.ind]

print(t) #End loop through time
}

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
"w.eq.full","ldg.sim", "ldgIMP.sim", "covIMP.sim", "Frs", "env_analogue_all",
"sim.impacts", "sim.impactsIMP","covIMP.impacts","lded","nf","is.final","Frs_ci")


#=============================================================================
#BIOMASS
#=============================================================================

#Calculate rr (intrinsic growth) from biomass measures of plants grown alone, 
#per site 
#Calculate average intraspecific invasion growth rate (biomass) per site
#Calculate resident stationary distribution from intraspecific competition
#Get competition (biomass) per site
#Calculate full invasion growth rate. 


#=============================================================================
#PLOTS
#=============================================================================
epoints = match( elevations,xx$elevation)

fig.name = paste("seeds_kriged.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

#par(mfrow=c(1,1),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
ylim=c(-1,30)
plot(xx$elevation, rr_krig[,1], t="l", ylab="Seeds (per-capita, Kriged) ", xlab="Elevation", xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim)
lines(xx$elevation, rr_krig[,1])
lines(xx$elevation, rr_krig[,2],col="red")
lines(xx$elevation, rr_krig[,3],col="blue")
abline(v=c(elevations),lwd=10,col="grey80")
points(xx$elevation[epoints], rr_krig[epoints,1])
points(xx$elevation[epoints], rr_krig[epoints,2],col="red")
points(xx$elevation[epoints], rr_krig[epoints,3],col="blue")

dev.off()
#=============================================================================

epoints = match( elevations,xx$elevation)

fig.name = paste("seeds_kriged2040.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

#par(mfrow=c(1,1),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
ylim=c(-1,70)
plot(xx$elevation, Frs[1,,1], t="l", ylab="Seeds (per-capita, Kriged) ", xlab="Elevation", xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim)
lines(xx$elevation, Frs[1,,1])
lines(xx$elevation, Frs[1,,2],col="red")
lines(xx$elevation, Frs[1,,3],col="blue")
abline(v=c(elevations),lwd=10,col="grey80")
points(xx$elevation[epoints], Frs[1,epoints,1])
points(xx$elevation[epoints], Frs[1,epoints,2],col="red")
points(xx$elevation[epoints], Frs[1,epoints,3],col="blue")

f1 = flower_prob_act*flower_act
points(xx$elevation[epoints], f1[,1],pch=3)
points(xx$elevation[epoints], f1[,2],col="red",pch=3)
points(xx$elevation[epoints], f1[,3],col="blue",pch=3)


dev.off()



#plot(xx$elevation, rr_krig[,2])
#plot(xx$elevation, rr_krig[,3])
#=============================================================================
fig.name = paste("intraLR_kriged.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

ylim=c(-1,6)
plot(xx$elevation, lrr_krig[,1], t="l", ylab="Seeds (per-capita, Kriged) ", xlab="Elevation", xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim)
lines(xx$elevation, lrr_krig[,1])
lines(xx$elevation, lrr_krig[,2],col="red")
lines(xx$elevation, lrr_krig[,3],col="blue")
abline(v=c(elevations),lwd=10,col="grey80")
points(xx$elevation[epoints], lrr_krig[epoints,1])
points(xx$elevation[epoints], lrr_krig[epoints,2],col="red")
points(xx$elevation[epoints], lrr_krig[epoints,3],col="blue")

dev.off()


#=============================================================================
#fig.name = paste("all4pairs_linear_widthOverlap_2Dgr4.pdf",sep="")
#pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow=c(3,1),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75)))
plot(xx$elevation, sr_krig[,1])
plot(xx$elevation, sr_krig[,2])
plot(xx$elevation, sr_krig[,3])
#=============================================================================
#fig.name = paste("all4pairs_linear_widthOverlap_2Dgr4.pdf",sep="")
#pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow=c(3,1),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75)))
plot(xx$elevation, surv_spp[1]+lrr_krig[,1])
plot(xx$elevation, surv_spp[2]+lrr_krig[,2])
plot(xx$elevation, surv_spp[3]+lrr_krig[,3])
#=============================================================================
#fig.name = paste("all4pairs_linear_widthOverlap_2Dgr4.pdf",sep="")
#pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow=c(3,1),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75)))
plot(xx$elevation, Crr_krig[,1])
plot(xx$elevation, Crr_krig[,2])
plot(xx$elevation, Crr_krig[,3])
#=============================================================================
fig.name = paste("interLR_kriged.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow=c(3,1))
col_use = c("black","red", "blue")
ylimits = matrix(c(0,15,0,8,0,1.5),3,2, byrow=T)
for(s in 1:length(spp)){
		sp_use = 1:(length(spp))
		sp_use= sp_use[-s]
		col1 = col_use[-s]	

		plot(xx$elevation, lii_all[[s]][[1]][,1], col= col1[1], t="l", ylab="Seeds (per-capita, Kriged) ", xlab="Elevation", xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylimits[s,])
		m1 = matrix(0,5,1)
		lines(xx$elevation, lii_all[[s]][[1]][,2],col=col1[2])
		abline(v=c(elevations),lwd=10,col="grey80")
		points(xx$elevation[epoints], lii_all[[s]][[1]][epoints,1],col=col1[1])
		points(xx$elevation[epoints], lii_all[[s]][[1]][epoints,2],col=col1[2])
		#for(sp in 1:2){
		# for(ee in 1:5){
		# 	m1[ee] = mean(subset(allB2, Sp == spp[s] & bg == substr(as.character(spp[sp_use[sp]]),1,1) & elevation == elevations[ee])$nr.inflor,na.rm=T)
		# 	}
		# points(elevations,m1)
}	

dev.off()
#=============================================================================
#fig.name = paste("all4pairs_linear_widthOverlap_2Dgr4.pdf",sep="")
#pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

fig.name = paste("arat1_kriged.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow=c(3,1))
col_use = c("black","red", "blue")
ylimits = matrix(c(0,1,0,5,0,5),3,2, byrow=T)
for(s in 1:length(spp)){
		sp_use = 1:(length(spp))
		sp_use= sp_use[-s]
		col1 = col_use[-s]	

		plot(xx$elevation, arat_all[[s]][[1]][,1], col= col1[1], t="l", ylab="Seeds (per-capita, Kriged) ", xlab="Elevation", xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylimits[s,])
		m1 = matrix(0,5,1)
		lines(xx$elevation, arat_all[[s]][[1]][,2],col=col1[2])
		abline(v=c(elevations),lwd=10,col="grey80")
		points(xx$elevation[epoints],arat_all[[s]][[1]][epoints,1],col=col1[1])
		points(xx$elevation[epoints], arat_all[[s]][[1]][epoints,2],col=col1[2])
		#for(sp in 1:2){
		# for(ee in 1:5){
		# 	m1[ee] = mean(subset(allB2, Sp == spp[s] & bg == substr(as.character(spp[sp_use[sp]]),1,1) & elevation == elevations[ee])$nr.inflor,na.rm=T)
		# 	}
		# points(elevations,m1)
	}	
dev.off()


#=============================================================================
#fig.name = paste("all4pairs_linear_widthOverlap_2Dgr4.pdf",sep="")
#pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow=c(3,2),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75)))
for(s in 1:length(spp)){
		sp_use = 1:(length(spp))
		sp_use= sp_use[-sp]

	for(sp in 1:2){ 
		plot(xx$elevation, arat_all2[[s]][[1]][,sp])
		print(median( arat_all2[[s]][[1]][,sp],na.rm=T ))
				
}}	
#=============================================================================
#fig.name = paste("all4pairs_linear_widthOverlap_2Dgr4.pdf",sep="")
#pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow=c(3,2),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75)))
for(s in 1:length(spp)){
		sp_use = 1:(length(spp))
		sp_use= sp_use[-sp]

	for(sp in 1:2){ 
		plot(xx$elevation, arat_all[[s]][[1]][,sp])
		points(xx$elevation, arat_all2[[s]][[1]][,sp],col="red")
		#print(median( arat_all[[s]][[1]][,sp],na.rm=T ))
				
}}	
	
#=============================================================================
#fig.name = paste("all4pairs_linear_widthOverlap_2Dgr4.pdf",sep="")
#pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow=c(3,2),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75)))
for(s in 1:length(spp)){
		sp_use = 1:(length(spp))
		sp_use= sp_use[-sp]

	for(sp in 1:2){ 
		plot(xx$elevation, Cir_all[[s]][[1]][,sp])
		m1 = matrix(0,5,1)

		for(ee in 1:5){
			m1[ee] = mean(subset(allB2, Sp == spp[s] & bg == substr(as.character(spp[sp_use[sp]]),1,1) & elevation == elevations[ee])$inter_flor,na.rm=T)
			}
		points(elevations,m1)
}}	
#=============================================================================
#fig.name = paste("all4pairs_linear_widthOverlap_2Dgr4.pdf",sep="")
#pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow=c(3,2),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75)))
for(s in 1:length(spp)){
		sp_use = 1:(length(spp))
		sp_use= sp_use[-sp]

	for(sp in 1:2){ 
		plot(xx$elevation, Cir_all[[s]][[1]][,sp])
		m1 = matrix(0,5,1)

		for(ee in 1:5){
			m1[ee] = mean(subset(allB2, Sp == spp[s] & bg == substr(as.character(spp[sp_use[sp]]),1,1) & elevation == elevations[ee])$inter_flor,na.rm=T)
			}
		points(elevations,m1)
}}	

#=============================================================================


sub1=((lrr_krig[,sb]) - surv_spp[sb]) / 
(surv_spp[sp] - ((lii_all[[sp]][[1]][,sb])))
sub2 =(-surv_spp[sp]+(lii_all[[sp]][[1]][,sb]) - rr_krig[,sp])/		
(surv_spp[sb]-(lrr_krig[,sb]) + rr_krig[,sb])

