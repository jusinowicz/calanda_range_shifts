#=============================================================================
# Determine species intrinsic ranges, do some analyses, and create figures
# of individual smooth fit terms.
#
# This code is mostly for diagnostics and producing useful plots that could
# end up in an appendix. 
#
# This code is for considering each abiotic variable in turn, as opposed to 
# fiting full models. It is most useful for visualizing the environmental space
# covered by sampling, and the response of ranges to each variable.
#
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
f.name1 = c("calanda_ranges_ccs85_tsmfull")
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
# In this version of the code, lists are made of the variables to incorporate,
# creating a hierarchy of models: 

v_use_list = list(
	variable_mat[c(2,6)], # gs_mean_temp   
	variable_mat[c(2,7)], # soil moisture
	variable_mat[c(2,4)], # min_temp
	variable_mat[c(2,3)], #  max_gdd
	variable_mat[c(2,1)] # elevation 
)

v_num = length(v_use_list)

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
# tstart = 1
# tstop = 61
# tinc = 10
# ttot= ceiling((tstop-tstart)/tinc)

# ngens=ttot #Number of generations for environmental change
ngens = 0 #Just initial ranges
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
# Key variables for output
# Note: These variables could be declared as 4D arrays, or as lists of arrays. 
# 4D arrays seem to use less memory:
# flower_probK_all = rep( list( array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),
# 		dim=c(ngenst,np,nspp)) ), v_num ) 
# flower_probK_all = array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),
#		dim=c(v_num,ngenst,np,nspp)) 
# Is this system dependent? Use object.size() to find out
# Using list approach, since length of data set to fit underlying environmental
# variables is unknown ahead of time. Could make space a 4D array still, but it
# seems easier to make all of the underlying code based on same indexing (i.e.
# list) 
#=============================================================================
### Intrinsic ranges
# Frs = rep( list( matrix(np,nspp) ), v_num ) # Mean model from posteriors: 
# Frs_mlp = array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp)) 
# Confidence intervals:
# Frs_ci = vector("list", ngenst)

### The fitted variables against environmental space
#Models
flower_prob_all = vector("list",v_num)
rr_gams_all = vector("list",v_num)
#Fits
flower_fit = vector("list",v_num)
rr_fit = vector("list",v_num)
flower_unc=vector("list",v_num) #Unconstrained fit
rr_unc = vector("list",v_num) #Unconstrained fit
#CIs
flower_var_ci=vector("list",v_num)
rr_var_ci=vector("list",v_num)

### Raw probability tables
raw_prob_all = vector("list", v_num)
raw_prob_dat = vector("list", v_num)
### Raw flower number data: 
raw_rr_all = vector("list", v_num)

### The fitted variables against elevation. This can be declared ahead of time 
### to speed up code since the spatial extent is known and fixed.  
#Fits
flower_probK_all = rep( list( matrix(0,np,nspp) ), v_num )
rr_krig_all = rep( list( matrix(0,np,nspp) ), v_num )
#CIs
flower_space_ci = rep( list( matrix(0,np,2*nspp) ), v_num )
rr_space_ci = rep( list( matrix(0,np,2*nspp) ), v_num )

#The maximum of species' fitted reproduction, for scaling
frs_max = matrix(0,ngenst,nspp)

#=============================================================================
# Outermost loop: through number of abiotic variables
#=============================================================================
par(mfrow = c(3,2))
for( nv in 1:v_num){ 

	v_use = v_use_list[[nv]]
	nvar = length(v_use) - 1
#=============================================================================
# Second loop: time
# This is just a single time step now, to get intrinsic ranges. 
# Is it worth implementing the shift with time? 
# 	for( t in 1:1){ 
# Ignoring time for now 	
#=============================================================================
	xx_new = xx
#=============================================================================
#FLOWERS
#=============================================================================

#=============================================================================
#intrinsic growth
#=============================================================================
#Intrinsic probability of flowering
	#Knots to use by default: 
	kn = c(3,3,3)
	
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

		### Find the elevation at which flowering peaks, then get the 
		### average value of variables corresponding to the peak elevation. 
				
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
			kuse= kn[sp] #Effectively the number of knots
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

			### Raw probabilities: 
			raw_prob = table(rr_dat$flyes,iv1)
			raw_prob = as.data.frame(raw_prob/colSums(raw_prob) )
			raw_prob = subset(raw_prob, Var1 == 1 )
			raw_prob[,1] = as.numeric(levels(raw_prob[,1]))[raw_prob[,1]]
			raw_prob[,2] = as.numeric(levels(raw_prob[,2]))[raw_prob[,2]]
			raw_prob_all[[nv]][[sp]] = raw_prob

		}

		factors = matrix("f",nvar)
		factors2 = matrix("f",nvar)

		for (f in 1:nvar){
			factors[f] = paste("X", f, sep="")
			#factors2[f] = paste("s(",v_use[(f+1)], ", k = ", kn[sp], ",  bs=\"cr\")", sep="")
			factors2[f] = paste("s(",v_use[(f+1)], ", k = ", kn[sp]+nk, ",  bs=\"cr\")", sep="")
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
		rr_tmp = predict(rr_gam, type="response", exclude = "s(year)")
		
		### Get confidence intervals for plot vs. variables
		### Coefficients from the penalized regression:	
		# Get the Xp matrix. 
		# Note: exclude = "s(year)" does not work with lpmatrix
		Xp =predict(rr_gam, type= "lpmatrix", exclude = "s(year)")[,1:kuse]
		# Posterior draw of model coefficients
		beta = coef(rr_gam)[1:kuse]
		Vp = rr_gam$Vp[1:kuse,1:kuse]
		br1=mvrnorm(n=10000, beta, Vp)
		#Get the CIs:
		mean.model= Xp%*%t(br1)+cp_y[1] #Need to add the offset
		mean.model = exp(mean.model)/(1+exp(mean.model))
		#Save the results
		flower_fit[[nv]][[sp]] = rr_tmp
		flower_var_ci[[nv]][[sp]] = apply(mean.model, 1, quantile, c(0.025,0.975)) 
		raw_prob_dat[[nv]][[sp]] = dat
		### Could turn this on for visuals: 
		plot(raw_prob$iv1,raw_prob$Freq,ylim=c(0,1))
		#amm = matrix(rowMeans(mean.model),dim(dat)[1], 1)
		lines(dat[,3][order(dat[,3])],rr_tmp[order(dat[,3])],col="red",t="l")
		bu_tmp = predict(b.u, type="response", exclude = "s(year)") #We need the minimum of this 
		lines(dat[,3][order(dat[,3])],bu_tmp[order(dat[,3])])
    fci=flower_var_ci[[nv]][[sp]]
		lines(dat[,3][order(dat[,3])],fci[,order(dat[,3])][1,] ,col="red",t="l")
		lines(dat[,3][order(dat[,3])],fci[,order(dat[,3])][2,] ,col="red",t="l")
		
    #Save the unconstrained fit
    flower_unc[[nv]][[sp]] = bu_tmp

		### coefficients from the penalized regression	
		#There are two approaches to getting beta and Xp.
		#In this version, the coefficients are used directly and 
		#the columns that would correspond to 0 coefficients in beta
		#have been dropped from Xp. The rationale for this is 
		beta = matrix(0, 1, (kuse*nvar+nvar*nk) )
		#Covariance matrix for getting CIs 
		Vp = matrix(0,(kuse*nvar+nvar*nk),(kuse*nvar+nvar*nk) )
		pnk = nk/2
		for( b in 1:nvar){
			bpos = (b-1)*(kuse+nk)+pnk+1
			kpos = (b-1)*kuse+1 #Intercept removed from rr_gam
			beta[ bpos:(bpos+kuse-1)] = coef(rr_gam) [kpos:(kpos+kuse-1)]
			Vp [bpos:(bpos+kuse-1),bpos:(bpos+kuse-1) ] = 
				rr_gam$Vp[kpos:(kpos+kuse-1),kpos:(kpos+kuse-1)]
		} 
		#beta1 = coef(rr_gam)[1:kuse]
		#Vp1 = rr_gam$Vp[1:kuse,1:kuse]
		
		### prediction matrix
		Xp = matrix(0, dim(xx_new)[1], nvar*(kuse+nk))
		#Xp1 = matrix(0, dim(xx_new)[1], nvar*kuse)
		#Xp = matrix(0, dim(dat)[1], nvar*(kuse+nk))

		for( b in 1:nvar){
			kpos = (b-1)*(kuse+nk)+1
			kpos1 = (b-1)*(kuse)+1
			Xp_tmp = PredictMat(sm[[b]], xx_new)
			#Xp_tmp = PredictMat(sm[[b]], dat)
			Xp[,kpos:(kpos+nk+kuse-1)] = Xp_tmp
			#Xp1[,kpos1:(kpos+kuse-1)] = PredictMat(sm[[b]], xx_new)[,(kpos+pnk):(kpos+nk+kuse-1-pnk) ]
		}
		
	

		### The predicted smooth, i.e. species' distribution
		#flower_tmp = Xp%*%beta+coef(rr_gam)[1] #These are not on the response scale
		flower_tmp = Xp%*%t(beta)+cp_y[1] #These are not on the response scale
		flower_tmp = exp(flower_tmp )/(1+exp(flower_tmp ) ) #Response scale for logistic regression
		
		### The unconstrained model: 
		#bu_tmp=as.vector(predict(b.u,newdata=xx_new,type= "response",exclude = "s(year)") )
		#flower_tmp=as.vector(predict.gam(flower_ptmp$gam, newdata=xx_new,type= "response"))
		#flower_tmp=as.vector(predict.gam(flower_ptmp, newdata=xx_new,type= "response"))

		### Get confidence intervals: 
		br1=mvrnorm(n=10000, beta, Vp)
		##This produces the mean model: 
		# mean.model= Xp%*%t(br1)+cp_y[1]
		# mean.model = rowMeans(exp(mean.model)/(1+exp(mean.model)))

		##But to get the CIs:
		mean.model= Xp%*%t(br1)+cp_y[1]
		mean.model = exp(mean.model)/(1+exp(mean.model))
		#q1=quantile(mean.model,c(.025,.975))
		flower_space_ci[[nv]][,((sp-1)*2+1):((sp-1)*2+2)] = t(apply(mean.model, 1, quantile, c(0.025,0.975)))

		### Keep the data
		flower_probK_all[[nv]][,sp] = flower_tmp
		#flower_probK[,sp] = bu_tmp
		flower_prob_all[[nv]][[sp]] = rr_gam

	}


	#Krig a distribution of Rs (intrinsic growth) from flower counts of plants 
	#grown alone, per site 
	#Use only those that have flowered
	#The best approach here seems to be to use the probability by elevation, 
	#then each R is actually a constant.

	#Default number of knots
	kn = c(3,3,3)

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
			kuse= kn[sp] #Effectively the number of knots
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
			factors2[f] = paste("s(",v_use[(f+1)], ", k = ", kn[sp]+nk, ",  bs=\"cr\")", sep="")
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
		rr_tmp = predict(rr_gam, type ="response", exclude = "s(year)")

		### Get confidence intervals for plot vs. variables
		### Coefficients from the penalized regression:	
		# Get the Xp matrix. 
		# Note: exclude = "s(year)" does not work with lpmatrix
		Xp =predict(rr_gam, type= "lpmatrix", exclude = "s(year)")[,1:kuse]
		# Posterior draw of model coefficients
		beta = coef(rr_gam)[1:kuse]
		Vp = rr_gam$Vp[1:kuse,1:kuse]
		br1=mvrnorm(n=10000, beta, Vp)

		#Get the CIs:
		mean.model= Xp%*%t(br1)+cp_y[1] #Need to add the offset

		#Save the results
		rr_fit[[nv]][[sp]] = rr_tmp
		rr_var_ci[[nv]][[sp]] = apply(mean.model, 1, quantile, c(0.025,0.975)) 
		raw_rr_all[[nv]][[sp]] = dat

		### Could turn this on for visuals: 
		plot(dat[,3],dat[,1])
		amm = matrix(rowMeans(mean.model),dim(dat)[1], 1)
		lines(dat[,3][order(dat[,3])],rr_tmp[order(dat[,3])],col="red",t="l")
		bu_tmp = predict(b.u, type="response", exclude = "s(year)") #We need the minimum of this 
		lines(dat[,3][order(dat[,3])],bu_tmp[order(dat[,3])])
		fci = rr_var_ci[[nv]][[sp]] 
		lines(dat[,3][order(dat[,3])],fci[,order(dat[,3])][1,] ,col="red",t="l")
		lines(dat[,3][order(dat[,3])],fci[,order(dat[,3])][2,] ,col="red",t="l")

   #Save the unconstrained fit
    rr_unc[[nv]][[sp]] = bu_tmp


		### coefficients from the penalized regression	
		beta = matrix(0, 1, (kuse*nvar+nvar*nk) )
		Vp = matrix(0,(kuse*nvar+nvar*nk),(kuse*nvar+nvar*nk) )
		pnk = nk/2
		for( b in 1:nvar){
			bpos = (b-1)*(kuse+nk)+pnk+1
			kpos = (b-1)*kuse+1 #Intercept removed from rr_gam
			beta[ bpos:(bpos+kuse-1)] = coef(rr_gam) [kpos:(kpos+kuse-1)]
			Vp [bpos:(bpos+kuse-1),bpos:(bpos+kuse-1) ] = 
				rr_gam$Vp[kpos:(kpos+kuse-1),kpos:(kpos+kuse-1)]
		} 
		

		### prediction matrix
		Xp = matrix(0, dim(xx_new)[1], nvar*(kuse+nk))
		#Xp = matrix(0, dim(dat)[1], nvar*(kuse+nk))

		for( b in 1:nvar){
			kpos = (b-1)*(kuse+nk)+1
			Xp_tmp = PredictMat(sm[[b]], xx_new)
			#Xp_tmp = PredictMat(sm[[b]], dat)
			Xp[,kpos:(kpos+nk+kuse-1)] = Xp_tmp
		}
		
		### the predicted smooth, i.e. species' distribution
		#rr_tmp = Xp%*%beta+coef(rr_gam)[1] #These are not on the response scale
		rr_tmp = Xp%*%t(beta) 

		#The unconstrained model: 
		bu_tmp=as.vector(predict(b.u,newdata=xx_new,type= "response",exclude = "s(year)") )
	
		### Get confidence intervals: 
		br1=mvrnorm(n=10000, beta, Vp)
		##This produces the mean model: 
		# mean.model= Xp%*%t(br1)+cp_y[1]
		# mean.model = rowMeans(exp(mean.model)/(1+exp(mean.model)))

		##But to get the CIs:
		mean.model= Xp%*%t(br1)+cp_y[1]
		#q1=quantile(mean.model,c(.025,.975))
		rr_space_ci[[nv]][,((sp-1)*2+1):((sp-1)*2+2)] = t(apply(mean.model, 1, quantile, c(0.025,0.975)))


		#Remove 0s
		rr_tmp[rr_tmp<0] = 0
		rr_krig_all[[nv]][,sp] = rr_tmp
		#rr_krig[,sp] = bu_tmp
		rr_gams_all[[nv]][[sp]]=rr_gam
		
	}


	# ####Intrinsic ranges
	# # Frs=array(c(rr_krig*flower_probK),dim=c(ngenst,np,nspp)) 
	# # #Frs[,,1] = Frs[,,1]*12
	# # Frs[,,1] = Frs[,,1]*8
	# ####Intrinsic ranges
	# Frs[]=matrix(rr_krig*flower_probK,np,nspp)

	# #Rescaled from data:
	# max(flower_act[,sa],na.rm=T)
	# for( sa in 1:nspp){

	# 	frs_max[t,sa] = max(Frs[t,,sa],na.rm=T)
	# 	#Make the scaling factor a percentage of the non-shifted scaling
	# 	max_scale =  frs_max[t,sa]/frs_max[1,sa]*max(flower_act[,sa],na.rm=T)
	# 	Frs[t,,sa]=Frs[t,,sa]*max_scale/max(Frs[t,,sa],na.rm=T)
	# }

	# Frs[t,,1] = Frs[t,,1]*3

	# par(mfrow=c(3,1))
	# plot(Frs[t,,2],t="l")
	# for (tt in 1:3){
	# 	lines(Frs[t,,tt])	
	# }

	# plot(rr_krig[,2],t="l")
	# for (tt in 1:3){
	# 	lines(rr_krig[,tt])	
	# }

	# plot(flower_probK[,2],t="l")
	# for (tt in 1:3){
	# 	lines(flower_probK[,tt])	
	# }

#} #Second loop, time 


} #Outermost loop, number of variables



spp_codes = c("DG","AX","HN")

for (sp in 1:nspp){

fig.name = paste("calanda_flower",spp_codes[sp],"_vars.pdf",sep="")
pdf(file=fig.name, height=11, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow=c(5,2), mar= c( 3, 0, 1, 0), oma=c(4,5,4,5) )

v_names = list( "Mean temperature", "Soil moisture", 
      "Minimum temperature", "Growing-degree days", "Elevation")

  for(nv in 1:v_num ) { 

    ### for prob variables, the probability of flowering
    raw_prob = raw_prob_all[[nv]][[sp]]
    raw_dat = raw_prob_dat[[nv]][[sp]]
    rr_tmp = flower_fit[[nv]][[sp]]
    bu_tmp = flower_unc[[nv]][[sp]]
    fci=flower_var_ci[[nv]][[sp]]

    plot(raw_prob$iv1,raw_prob$Freq,ylim=c(0,1), yaxt="n", xlab="", ylab = "")
    axis(2,cex.lab = 1.5)
    mtext(paste(v_names[[nv]]),side=1, line=2.5, at=par("usr")[1]+1*diff(par("usr")[1:2]), cex.lab=1)
    axis(1, cex.lab = 1.5)
    if(nv ==3 ){ mtext("Probability of flowering",side=2,line=2.5)}
    lines(raw_dat[,3][order(raw_dat[,3])],rr_tmp[order(raw_dat[,3])],col="red",t="l")
    lines(raw_dat[,3][order(raw_dat[,3])],bu_tmp[order(raw_dat[,3])])
    lines(raw_dat[,3][order(raw_dat[,3])],fci[,order(raw_dat[,3])][1,] ,col="red",t="l",lty =3)
    lines(raw_dat[,3][order(raw_dat[,3])],fci[,order(raw_dat[,3])][2,] ,col="red",t="l",lty =3)



  ### for rr variables, the number of flowers
    
    raw_rr = raw_rr_all[[nv]][[sp]]
    rr_tmp = rr_fit[[nv]][[sp]]
    bu_tmp = rr_unc[[nv]][[sp]]
    fci=rr_var_ci[[nv]][[sp]]
    plot(raw_rr[,3],raw_rr[,1], yaxt="n", ylab = "", xlab="")
    axis(4,cex.lab = 1.5)
    if(nv ==3 ){ mtext("Number of flowers", side = 4, line=2.5, cex.lab=1)}
    lines(raw_rr[,3][order(raw_rr[,3])],rr_tmp[order(raw_rr[,3])],col="red",t="l")
    lines(raw_rr[,3][order(raw_rr[,3])],bu_tmp[order(raw_rr[,3])])
    lines(raw_rr[,3][order(raw_rr[,3])],fci[,order(raw_rr[,3])][1,] ,col="red",t="l",lty =3)
    lines(raw_rr[,3][order(raw_rr[,3])],fci[,order(raw_rr[,3])][2,] ,col="red",t="l",lty =3)



  }
dev.off()
}