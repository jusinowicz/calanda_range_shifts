#=============================================================================
# In this version, fit probability of flowering and survival by size class. 
# Fit the parameters only to plants that have survived. 
# Fit reproduction only to those that have flowered. 
#=============================================================================
#=============================================================================
# Load these libraries
#=============================================================================
library(mgcv)
library(MASS)
library(fields)
library(viridis)
source("./wald_functions1.R")
source("./range_coexistence_functionsWALD.R")
#source("./range_coexistence_functionsD.R") #The old library

#=============================================================================
#For naming files
#=============================================================================
f.name1=c("calanda_igrsites12it_rcp85_multi_2017")
#f.name1=c("calanda_igrsites12it_rcp26_multi_2017")
#=============================================================================
#Data: (see calanda2017.R)
#=============================================================================
load(file="calanda_allB2017.var")
head(allB)
allB$year=as.numeric(allB$year)
#Dactylis glomerata conversion factor for flowers to seeds: 
#dg_conv = 1/10
# allB[ as.numeric(rownames(subset(allB,Sp =="DG") )), ]$nr.inflor = 
# 	allB[ as.numeric(row.names(subset(allB,Sp =="DG") )), ]$nr.inflor *dg_conv 

#With climate-change scenarios:
ccs_temp = read.csv("rcp85_2019_2100.csv")[2:3]
#ccs_temp = read.csv("rcp26_2019_2100.csv")[2:3]

#=============================================================================
#Tunable lattice parameters (space and time)
#=============================================================================
#If fast.bys = TRUE, then populations will be allowed to spread and track their
#intrinsic ranges. 
fast.bys=TRUE
c.tol = 1e-2 #Tolerance for finite inidividuals, sets left/right range edges

#For the climate change scenarios, tune the change in temp to the current 
#average from field data: 
temp_field = mean(allB$gs_mean_temp)
ccs_temp = ccs_temp[,2]- temp_field
ccs_temp=c(0,ccs_temp)
mstart = 1
mstop = length(ccs_temp)
minc = 1
mtot= length(ccs_temp)

#For standard increments: 
# mstart = 0
# mstop = 6
# minc = 0.5
# mtot= ceiling((mstop-mstart)/minc)

ngens=mtot-1 #Number of generations for environmental change
iconfig=1 #NUmber of interattions of initial configuration of species

#Spatial scale: Assume that the gradient is from 500m to 2500m in units of 10m
elevations=c(1000,1400,1600,1800,2000)
el1=0
el2=3500

#Make some internal spatial variables:
xx1= matrix(seq(el1,el2,1))
xx2= matrix(xx1^2)
#yx= matrix(2017,length(xx1),1)
yx= matrix(3,length(xx1),1)


xx=cbind(xx1,yx,xx2)
colnames(xx)=c('elevation','year','e2')
xx=data.frame(xx)
xx$elevation = as.matrix(xx$elevation)
xx$year = as.numeric(xx$year)
xx$e2 = as.matrix(xx$e2)

#=============================================================================
#Environmental variation -- defining the abiotic influence on intrinsic growth
#rates and projecting it into the future. 
#=============================================================================

#=============================================================================
#Fit a GAM with each of the environmental variables against elevation
#Note: B2 is the subset of B without the first year of data, which
#is not trustworthy for the hobos. 

#Temp
temp_gam = gam( gs_mean_temp~ s(year,bs="re")+s(elevation,k=3),data=allB)
#temp_gam = gam( gs_mean_temp~ s(year,bs="re")+elevation,data=allB)

#soil moisture 
sm_gam = gam( soil_moist~ s(year,bs="re")+s(elevation,k=3),data=allB)
#sm_gam = gam( soil_moist~ s(year,bs="re")+elevation,data=allB)

#PET
sp_gam = gam( soil_pet~ s(year,bs="re")+s(elevation,k=3),data=allB)
#sp_gam = gam( soil_pet~ s(year,bs="re")+elevation,data=allB)

#Moisture deficit
sd_gam = gam( soil_def~ s(year,bs="re")+s(elevation,k=3),data=allB)
#sd_gam = gam( soil_def~ s(year,bs="re")+elevation,data=allB)

#annual mean light
ml_gam = gam( an_mean_light~ s(year,bs="re")+s(elevation,k=3),data=allB)
#ml_gam = gam( an_mean_light~ s(year,bs="re")+elevation,data=allB)

#GDD
gdd_gam = gam( max_gdd~ s(year,bs="re")+s(elevation,k=3),data=allB)
#gdd_gam = gam( max_gdd~ s(year,bs="re")+elevation,data=allB)

#Minimum growing season temp
gmt_gam = gam( gs_min_temp~ s(year,bs="re")+s(elevation,k=3),data=allB)
#gmt_gam = gam( gs_min_temp~ s(year,bs="re")+elevation,data=allB)

#=============================================================================


#The mean temp vector: 
xx_mt = as.vector(predict.gam(temp_gam, newdata=xx,type= "response"))

#The soil moisture vector: 
xx_sm = as.vector(predict.gam(sm_gam, newdata=xx,type= "response"))

#The soil moisture vector: 
xx_sp = as.vector(predict.gam(sp_gam, newdata=xx,type= "response"))

#The soil moisture vector: 
xx_sd = as.vector(predict.gam(sd_gam, newdata=xx,type= "response"))

#Mean light
xx_ml = as.vector(predict.gam(ml_gam, newdata=xx,type= "response"))

#GDD
xx_gdd = as.vector(predict.gam(gdd_gam, newdata=xx,type= "response"))

#Min temp 
xx_gmt = as.vector(predict.gam(gmt_gam, newdata=xx,type= "response"))

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

spp = unique(as.character(allB$Sp))
nspp = length(spp)

bgs = unique(as.character(allB$bg))
pop_labels = paste(expand.grid(spp,bgs,elevations)[,1],expand.grid(spp,bgs,elevations)[,
		2],expand.grid(spp,bgs,elevations)[,3],sep="")
years= unique(as.character(allB$year))

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
for(sp in 1:length(spp)){
	allsp = subset(allB, Sp == spp[sp])
	sr_dat=subset(allsp, bg =='B')
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
Crr_nls_el = matrix(0,length(elevations), length(spp)) #Elevations individually

kn = c(4,4,4)
bsamp = 10000

for(sp in 1:length(spp)){

	allsp = subset(allB, Sp == spp[sp] )	
	Cr_dat =subset(allsp, bg == substr(as.character(spp[sp]),1,1))
	rr_sub =subset(allsp, bg =='B')
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
	for(ee in 1:length(elevations)){
		Cr_dat_e = subset(Cr_dat, elevation == elevations[ee])
		rr_sub_e = subset(rr_sub, elevation == elevations[ee])
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
	Cr_k = data.frame(cbind(Crr_nls_el[, sp], elevations))
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
	Cir_nls_el = matrix(0,length(elevations), length(spp)) #Elevations individually

	allsp = subset(allB, Sp == spp[sp])
	rr_sub =subset(allsp, bg =='B')
	nsamp = dim(rr_sub)[1]

	sp_use = 1:(length(spp))
	sp_use= sp_use[-sp]

	for(sa in 1: (length(spp)-1)){
		
		sb= sp_use[sa]
		Ci_dat =subset(allsp, bg == substr(as.character(spp[sb]),1,1))
		
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
		for(ee in 1:length(elevations)){
			Ci_dat_e = subset(Ci_dat, elevation == elevations[ee])
			rr_sub_e = subset(rr_sub, elevation == elevations[ee])
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
	
		#Now use the per elevation data to Krige other values
		Ci_k = data.frame(cbind(Cir_nls_el[, sb], elevations))
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
		arat3[sp,sb]=(Cir_nls[sp,sb])/(Crr_nls[sp])
		
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
	crr = Crr_nls_el[,sp]
	
	for(sa in 1: (length(spp)-1)){
		sb= sp_use[sa]
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
sr=colMeans(sr_krigB)

###Competition coeffcients
#3B is calculating a per-site alpha, then averaging 
diag(arat3B)=c(1,1,1)
alphas=arat3B
alphas[2,1] = arat3[2,1]

#3 is calculating an alpha that has been calculated regardless
#of elevation
#diag(arat3)=c(1,1,1)
#alphas=arat3


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
	allsp = subset(allB, Sp == spp[sp])
	ht_dat_bg =subset(allsp, bg =='B')
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
y.full=vector(mtot, mode="list")
w.eq.full=vector(mtot, mode="list")

#Site impacts on the LGR
sim.impacts=array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp))
sim.impactsIMP=array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp))
covIMP.impacts=array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp))
lded=array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp)) 

env_analogue_all = vector(mtot, mode="list")
	
#=============================================================================
#Outermost loop: treat each species as invader
#=============================================================================

#Keep track of total interspecific competition of species when it is invader
cc=array(c(matrix(0,ngenst,np),matrix(0,ngenst,np)),dim=c(ngenst,np,nspp)) 

for( t in 1: ngenst){ 

	#ti = seq(mstart,mstop,minc)[t]
	ti = ccs_temp[t]

	xx_new = xx

	xx_new$gs_mean_temp = as.matrix(xx_new$gs_mean_temp+ti)
	xx_new$gs2 = as.matrix(xx_new$gs_mean_temp^2)

	#This code takes the temp at sites in the future and matches it to
	#the equivalent current site based on that temp.
	new_ind=sapply(xx_new$gs_mean_temp,function(x) which.min(abs(x-xx$gs_mean_temp)))
	xx_new$max_gdd = as.matrix(xx$max_gdd[new_ind])
	xx_new$gs_min_temp = as.matrix(xx$gs_min_temp[new_ind])

	env.ind = c(1,4,7,9,10)
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
	flower_prob=NULL
	flower_prob_act = matrix(0,5,3)

	kn = c(5,5,3)

	env_distance_spp=list()
	env_analogue_spp=list()

	for(sp in 1:length(spp)){
		allsp = subset(allB, Sp == spp[sp])
		rr_dat =subset(allsp, bg =='B')
		
		#What are actual per-site probs? 
		for (n in 1:5) { 
			tfp = subset(rr_dat, elevation == elevations[n])
			flower_prob_act[n,sp]= sum(tfp$flyes)/nrow(tfp)}

		#Run different models for species 1, vs 2 and 3
		if( sp == 1) { 
		
			#flower_ptmp= gam(flyes ~ year+s(elevation,k=kn[sp])+te(gs_mean_temp,soil_def,gs_min_temp),family=binomial(link='logit'),data=rr_dat)
			flower_ptmp= gam(flyes ~ year+s(elevation,k=kn[sp])+te(gs_mean_temp,soil_def,gs_min_temp),family=binomial(link='logit'),data=rr_dat)

			env.ind.spp = c(1, 4, 7,10) 
			env_distance_spp[[sp]] = get_env_distance(xx,xx_new,env.ind.spp)
			env_analogue_spp[[sp]] = get_analogue(xx,xx_new,env.ind.spp)

		} else if (sp == 2) {

			flower_ptmp= gam(flyes ~ year+s(elevation,k=kn[sp])+s(gs_mean_temp,k=kn[sp])+s(soil_def,k=kn[sp])+s(max_gdd,k=kn[sp])+s(gs_min_temp,k=kn[sp]),family=binomial(link='logit'),data=rr_dat)
			
			env.ind.spp = c(1, 4, 7,9,10) 
			env_distance_spp[[sp]] = get_env_distance(xx,xx_new,env.ind.spp)
			env_analogue_spp[[sp]] = get_analogue(xx,xx_new,env.ind.spp)


		} else if (sp == 3) {

			flower_ptmp= gam(flyes ~ year+s(gs_mean_temp,k=kn[sp])+s(soil_def,k=kn[sp])+s(gs_min_temp,k=kn[sp]),family=binomial(link='logit'),data=rr_dat)
			
			env.ind.spp = c(4,7,10) 
			env_distance_spp[[sp]] = get_env_distance(xx,xx_new,env.ind.spp)
			env_analogue_spp[[sp]] = get_analogue(xx,xx_new,env.ind.spp)

		}


		#Use the model to predict values: 
		#Under warmer conditions, uncomment:
		flower_tmp=as.vector(predict.gam(flower_ptmp, newdata=xx_new,type= "response"))

		#Under current conditions, uncomment:
		#flower_tmp=as.vector(predict.gam(flower_ptmp, newdata=xx,type= "response"))

		#Remove 0s
		#flower_tmp[flower_tmp<0] = 0
		flower_probK[,sp] = flower_tmp
		flower_prob[[sp]] = flower_ptmp
	}


	#Krig a distribution of Rs (intrinsic growth) from flower counts of plants 
	#grown alone, per site 
	#Use only those that have flowered
	#The best approach here seems to be to use the probability by elevation, 
	#then each R is actually a constant.

	rr_krig = matrix(0,dim(xx)[1], length(spp))
	rr_gams=NULL
	kn = c(3,3,3)
	flower_act =matrix(0,5,3)

	for(sp in 1:length(spp)){
		allsp = subset(allB, Sp == spp[sp])
		rr_dat_bg =subset(allsp, bg =='B')
		rr_dat =subset(rr_dat_bg, flyes ==1)
		
		for (n in 1:5) { 
			tfp = subset(rr_dat, elevation == elevations[n])
			flower_act[n,sp]= mean(tfp$nr.inflor)}


		#1)
		rr_gam= gam(nr.inflor ~ year+s(elevation,k=kn[sp])+s(gs_mean_temp,k=kn[sp])+s(soil_def,k=kn[sp])+s(gs_min_temp,k=kn[sp]),data=rr_dat)
		
		#Use the model to predict values: 
		#Under warmer conditions, uncomment:
		rr_tmp=as.vector(predict.gam(rr_gam, newdata=xx_new,type= "response"))

		#Under current conditions, uncomment:
		#rr_tmp=as.vector(predict.gam(rr_gam, newdata=xx,type= "response"))
		

		#Remove 0s
		rr_tmp[rr_tmp<0] = 0
		rr_krig[,sp] = rr_tmp
		rr_gams[[sp]]=rr_gam
	}

	####Intrinsic ranges
	# Frs=array(c(rr_krig*flower_probK),dim=c(ngenst,np,nspp)) 
	# #Frs[,,1] = Frs[,,1]*12
	# Frs[,,1] = Frs[,,1]*8
	####Intrinsic ranges
	Frs[t,,] =matrix(rr_krig*flower_probK,np,nspp) 
	#Frs[,,1] = Frs[,,1]*8*1.5
	Frs[t,,1] = Frs[t,,1]*8

	#Frs[,,1] =Frs[,,1] - matrix(apply(matrix(Frs[,,1]),2,min),np, ngenst, byrow=T)
	par(mfrow=c(1,1))
	plot(Frs[t,,1],t="l")
	for (tt in 2:3){
		lines(Frs[t,,tt])	
	}

	#=============================================================================
	#Stationary multispecies distribution through numerical integration
	#=============================================================================

	s.index1= 1:nspp
	burns=1000
	#Equilibrium populations of all spp. Columns represent space, rows are time. 
	nf.tmp=array(c(matrix(0.0,burns,np),matrix(0.0,burns,np)),dim=c(burns,np,nspp)) 
	for ( sa in 1:nspp) { nf.tmp[1,,sa] = matrix(0.01,1,np)}
	#Get the stationary multispecies distribution through numerical integration
	tcomp=0 #Condition for stationary (joint) distribution of community
	ts = 1
	comp_thresh = 1 #Sensitivity threshold. Theoretically, should be 0, but allow for numerical discrepency

	while (tcomp == 0 ){

		nf.tmp[ts+1,,] = pop_lg(Frs[t,,], nf.tmp[ts,,], sr,alphas, fkd, fkc,fkd.yes  )
		
		#Test for stationarity of community
		#if ( (mean(mean(nf[(ts+1),,]/nf[(ts),,]))-1) <= comp_thresh){ tcomp=1; is.final=T}
		if ( round(mean(mean(nf.tmp[(ts+1),,]/nf.tmp[(ts),,]))-1,comp_thresh) == 0){ tcomp=1; is.final=T}

		#Stop if burn limit is reached
		if(ts >= burns) {tcomp=1; is.final=F}
		ts = ts+1

	}

	nf[t,,] = nf.tmp[(ts-1),,]

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
					for (tt in 1:3){
					lines(Frs[t,,tt], lty=3,col="blue")}	

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

	y.full[[t]]= y 
	w.eq.full[[t]]=w.eq
	env_analogue_all[[t]] = env_analogue 

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
"w.eq.full","ldg.sim", "ldgIMP.sim", "covIMP.sim", "Frs", "env_analogue",
"sim.impacts", "sim.impactsIMP","covIMP.impacts","lded","nf","is.final")


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
		# 	m1[ee] = mean(subset(allB, Sp == spp[s] & bg == substr(as.character(spp[sp_use[sp]]),1,1) & elevation == elevations[ee])$nr.inflor,na.rm=T)
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
		# 	m1[ee] = mean(subset(allB, Sp == spp[s] & bg == substr(as.character(spp[sp_use[sp]]),1,1) & elevation == elevations[ee])$nr.inflor,na.rm=T)
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
			m1[ee] = mean(subset(allB, Sp == spp[s] & bg == substr(as.character(spp[sp_use[sp]]),1,1) & elevation == elevations[ee])$inter_flor,na.rm=T)
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
			m1[ee] = mean(subset(allB, Sp == spp[s] & bg == substr(as.character(spp[sp_use[sp]]),1,1) & elevation == elevations[ee])$inter_flor,na.rm=T)
			}
		points(elevations,m1)
}}	

#=============================================================================


sub1=((lrr_krig[,sb]) - surv_spp[sb]) / 
(surv_spp[sp] - ((lii_all[[sp]][[1]][,sb])))
sub2 =(-surv_spp[sp]+(lii_all[[sp]][[1]][,sb]) - rr_krig[,sp])/		
(surv_spp[sb]-(lrr_krig[,sb]) + rr_krig[,sb])

