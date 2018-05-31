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
source("./wald_functions1.R")
source("./range_coexistence_functionsWALD.R")
#source("./range_coexistence_functionsD.R") #The old library

#=============================================================================
#For naming files
#=============================================================================
f.name1=c("calanda_igr2_2017")
#=============================================================================
#Data: (see alanda2017.R)
#=============================================================================
head(allB)

#Dactylis glomerata conversion factor for flowers to seeds: 
dg_conv = 10
allB[ as.numeric(row.names(subset(allB,Sp =="DG") )), ]$nr.inflor = 
	allB[ as.numeric(row.names(subset(allB,Sp =="DG") )), ]$nr.inflor *dg_conv 

#=============================================================================
#Tunable lattice parameters (space and time)
#=============================================================================

ngens=0 #Number of generations for environmental change
iconfig=1 #NUmber of interattions of initial configuration of species

#Spatial scale: Assume that the gradient is from 500m to 2500m in units of 10m
elevations=c(1000,1400,1600,1800,2000)
xx1= matrix(seq(900,2100,1))
xx2= matrix(xx1^2)
yx= matrix(2016,length(xx1),1)

xx=cbind(xx1,yx,xx2)
colnames(xx)=c('elevation','year','e2')
xx=data.frame(xx)
xx$elevation = as.matrix(xx$elevation)
xx$year = as.matrix(xx$year)
xx$e2 = as.matrix(xx$e2)


#Spatial scale: Assume that the gradient is from 500m to 2500m in units of 10m
x1=  as.matrix(seq(900,2100,1))
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
#FLOWERS
#=============================================================================

#=============================================================================
#intrinsic growth
#=============================================================================
#Intrinsic probability of flowering

flower_probK=matrix(0,dim(xx)[1], length(spp))
flower_prob=NULL
kn = c(4,4,4)

for(sp in 1:length(spp)){
	allsp = subset(allB, Sp == spp[sp])
	rr_dat =subset(allsp, bg =='B')
	
	#With only the effect of year removed
	#flower_ptmp= gam(flyes ~ s(year,bs="re")+s(elevation,k=kn[sp]),family=binomial(link='logit'),data=rr_dat)
	
	#Use splines to account for spatial distribution
	#flower_ptmp= gam(flyes ~ year+area+s(elevation,k=5),family=binomial(link='logit'),data=rr_dat)
	
	#Use polynomial model (assumes elevation^2 exists)
	flower_ptmp= gam(flyes ~ s(year,bs="re")+elevation +e2,family=binomial(link='logit'),data=rr_dat)
	flower_tmp=as.vector(predict.gam(flower_ptmp, newdata=xx,type= "response"))
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
kn = c(5,5,5)
for(sp in 1:length(spp)){
	allsp = subset(allB, Sp == spp[sp])
	rr_dat_bg =subset(allsp, bg =='B')
	rr_dat =subset(rr_dat_bg, flyes ==1)
	
	#With only the effect of year removed
	#rr_gam = gam(nr.inflor~s(year,bs="re"), data=rr_dat)
	
	#Use splines to account for spatial distribution
	#rr_gam = gam(nr.inflor~s(year,bs="re")+s(elevation,k=kn[sp]), data=rr_dat)
	
	#Use polynomial model (assumes elevation^2 exists)
	rr_gam = gam(nr.inflor~s(year,bs="re")+elevation+e2, data=rr_dat)
	rr_tmp=as.vector(predict.gam(rr_gam, newdata=xx,type= "response"))
	#Remove 0s
	rr_tmp[rr_tmp<0] = 0
	rr_krig[,sp] = rr_tmp
	rr_gams[[sp]]=rr_gam
}

#Fix the lower-elevation tail of DG and AX (check if this is still necessary): 
# lt = c(120,100)
# for (le in 1:2) { 
# 	nl = min(rr_krig[1:(lt[le]),le])
# 	nl_min = which(rr_krig[1:(lt[le]),le] == nl )
# 	#Replace higher values to the left with nl
# 	rr_krig[1:nl_min,le][ rr_krig[1:nl_min,le] > nl] = nl

# 	#For flower probabilities:
# 	nl = min(flower_probK[1:(lt[le]),le])
# 	nl_min = which(flower_probK[1:(lt[le]),le] == nl )
# 	#Replace higher values to the left with nl
# 	flower_probK[1:nl_min,le][ flower_probK[1:nl_min,le] > nl] = nl
# }

# #Tweaking these more so that they don't stay level, but tape exponentially
# #Species 1:
# le=1
# 	infl1 = which(flower_probK[,le] ==max(flower_probK[,le]))
# 	to_flip = flower_probK[1:infl1]


# #Species 2: 
# le =2 
# #Exponential taper through inflection points on left side: 
# 	#point 1: Inflection point 
# 	infl1 = which(diff(rr_krig[,le]) ==max(diff(rr_krig[,le]))) 
# 	xe1 = xx[infl1,1]
# 	ye1 = rr_krig[infl1,le]

# 	#point 2: Minimum
# 	infl2 = min(which((diff(rr_krig[,2]))>0))
# 	xe2 = xx[infl2,1]
# 	ye2 = rr_krig[infl2,le]

# 	ae=1

# 	#Use a log-linear slope to fit the exponential between the two points
# 	m = (log(ye1)-log(ye2))/(xe1-xe2)
# 	b = -(xe2*log(ye1)-xe1*log(ye2))/(xe1-xe2)

# 	new_rr = exp(m*xx[,1]+b)
# 	rr_krig [1:infl2,le] = new_rr[1:infl2]


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
#This version retains the plot-specific lrr, but uses an average Rr across plots
#at an elevation based on the fitted models (rr_krig and flower_probK)

Crr_nls2 = matrix(0,length(spp),1) #Elevations grouped
Crr_nls_krig2 = matrix(0,dim(xx)[1], length(spp)) #Krig across elevations
Crr_nls_el2 = matrix(0,length(elevations), length(spp)) #Elevations individually

kn = c(4,4,4)

par(mfrow=c(3,1))
for(sp in 1:length(spp)){
	allsp = subset(allB, Sp == spp[sp])
	Cr_dat =subset(allsp, bg == substr(as.character(spp[sp]),1,1))
	rr_sub =cbind(xx[,1], rr_krig[,sp], flower_probK[,sp])
	rr_table = rr_sub[rr_sub[,1]%in%Cr_dat$elevation,]
	crdat= data.frame(cbind(Cr_dat$nr.inflor, 	rr_table[match(Cr_dat$elevation,rr_table),]))
	colnames(crdat) = c("lrr","elevation","Rr", "fp")
	
	Cr_nls = nls(lrr~ (Rr*fp)/(1+cr),data=crdat, start=list(cr=0.5))
	#print(paste(sp))
	#print(summary(Cr_nls))
	Crr_nls2[sp] = Cr_nls$m$getPars()
	
	#Now, per elevation: 
	plot(NA, NA, xlab="Average intrinsic R per site",ylab="LDG per plot", ylim=c(0, max(Cr_dat$nr.inflor,na.rm=T)+10 ),xlim = c(0, max(rr_krig)+10)) 
	for(ee in 1:length(elevations)){
		crdat_e = subset(crdat, elevation == elevations[ee])
		points(crdat_e$Rr ,crdat_e$lrr )
		tryCatch( {Cr_nls_e = nls(lrr~ (Rr*fp)/(1+cr),data=crdat_e, start=list(cr=0.5),control = list(maxiter = 500))
			Crr_nls_el2[ee, sp] = Cr_nls_e$m$getPars()
			#print(paste(sp))
			#print(summary(Cr_nls_e))
			}, error = function(e){})
		
	}
	
	#Now use the per elevation data to Krige other values
	Cr_k = data.frame(cbind(Crr_nls_el2[, sp], elevations))
	colnames(Cr_k) = c("cr","elevation")
	Cr_gam = gam(cr~s(elevation,k=kn[sp]),data=Cr_k )
	Cr_tmp=as.vector(predict.gam(Cr_gam, newdata=x1))
	#Remove 0s
	Cr_tmp[Cr_tmp<0] = 0
	Crr_nls_krig2[,sp] = Cr_tmp
}



#=============================================================================
#interspecific growth
#=============================================================================
#See calanda_pop_krigALL_wald.R for the full statistical exploration
#=============================================================================
#Use NLS to fit the functional relationship (leslie-gower?). Try it per elevation
#too, to see if it is significant. If it is, Krige it. If not...? 
#This version retains the plot-specific lir, but uses an average Rr across plots
#at an elevation based on the fitted models (rr_krig and flower_probK)


Cir_nls2 = matrix(0,length(spp),length(spp)) #Elevations grouped
Cir_nls_all2 = NULL
Cir_el_all2 = NULL

kn = c(4,4,4)
par(mfrow=c(3,2))
for(sp in 1:length(spp)){

	Cir_nls_krig = matrix(0,dim(xx)[1], length(spp)) #Krig across elevations
	Cir_nls_el = matrix(0,length(elevations), length(spp)) #Elevations individually

	allsp = subset(allB, Sp == spp[sp])
	rr_sub =cbind(xx[,1], rr_krig[,sp], flower_probK)
	rr_table = rr_sub[rr_sub[,1]%in%Cr_dat$elevation,]

	sp_use = 1:(length(spp))
	sp_use= sp_use[-sp]

	for(sa in 1: (length(spp)-1)){
		
		sb= sp_use[sa]
		Ci_dat =subset(allsp, bg == substr(as.character(spp[sb]),1,1))
		cidat= data.frame(cbind(Ci_dat$nr.inflor, rr_table[match(Ci_dat$elevation,rr_table),]))
		colnames(cidat) = c("lir","elevation","Ri", "fp")
		Ci_nls = nls(lir~ (Ri*fp)/(1+ci),data=cidat, start=list(ci=0.5))
		#print(paste(sp, sb))
		#print(summary(Ci_nls))
		Cir_nls2[sp,sb] = Ci_nls$m$getPars()
		
		#Now, per elevation: 
		plot(NA, NA, xlab="Average intrinsic R per site",ylab=paste("LDG per plot for",spp[sp] ,sep=" "), ylim=c(0, max(Cr_dat$nr.inflor,na.rm=T)+10 ),xlim = c(0, max(rr_krig,na.rm=T)+10)) 
		for(ee in 1:length(elevations)){
			cidat_e = subset(cidat, elevation == elevations[ee])
			points(cidat_e$Ri ,cidat_e$lir )
			tryCatch( {Ci_nls_e = nls(lir~ (Ri*fp)/(1+ci),data=cidat_e, start=list(ci=0.5),control = list(maxiter = 500))
				Cir_nls_el[ee, sb] = Ci_nls_e$m$getPars()
				#print(paste(sp,sb))
				#print(summary(Ci_nls_e))
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
	Cir_nls_all2[[sp]] = list(Cir_nls_krig)
	Cir_el_all2[[sp]] = list(Cir_nls_el)
}

#=============================================================================
#air/arr
#=============================================================================
#See calanda_pop_krigALL_wald.R for the full statistical exploration
#=============================================================================
#Get competition (flowers) per site. This is based on solving each IGR for arr*Nr 
#and air*Nr, then taking the ration air/arr, setting arr=1
#For flowers, this doesn't include the survival term so: 
#air/arr = (m*(n+li))/((m+lr)*n)
#Version 4: use the constant alphas fit with NLS

arat4 = matrix(0,length(spp),length(spp))

for(sp in 1:length(spp)){
	sp_use = 1:(length(spp))
	sp_use= sp_use[-sp]

	for(sa in 1: (length(spp)-1)){
		sb= sp_use[sa]
		#Use medians of IGRs
		arat4[sp,sb]=(Cir_nls2[sp,sb])/(Crr_nls2[sp])
		
	}
}		

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
# diag(arat3)=c(1,1,1)
# alphas=arat3
diag(arat4)=c(1,1,1)
alphas=arat4

####Intrinsic ranges
Frs=array(c(rr_krig*flower_probK),dim=c(ngenst,np,nspp)) 
Frs[,,1] = Frs[,,1]*2
#Frs[,,1] =Frs[,,1] - matrix(apply(matrix(Frs[,,1]),2,min),np, ngenst, byrow=T)

###Competition distance -- These are unknown, and chosen just to be limited or not
#b_rr=1/(100*np) #Essentially global
b_rr=c(.1,.1,.1) #10 cm plots

###Dispersal distance -- These are calculated using the WALD model --see 
#wald_model1.R

#a_rr=1/(100*np) #essentially global
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

####Dispersal kernels and their Fourier transforms 

fast.bys=FALSE


#Get the WALD kernel based on field measurements
u_mean = 2.85 #Mean windspeed above canopy
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

#Set the means (for an exponential) to what they are in the WALD kernel: 
a_rr = c(1/wald.list[[1]][[2]],1/wald.list[[2]][[2]],1/wald.list[[3]][[2]])

#Exponential kernel
kd=matrix(0,np,nspp)
fkd=matrix(0,np,nspp)

for( s in 1:nspp){ 
	kd[,s] = a_rr[s]/2*exp(-a_rr[s]*abs(xx0))
	kd[,s]=kd[,s]/(sum(kd[,s]))
	fkd[,s]=fft(kd[,s])#/(np+1)
	fkd.yes = TRUE #Pass the transformed kd to functions later for speed

}


####Competition kernels and their Fourier transforms 
kc=matrix(0,np,nspp)
fkc=matrix(0,np,nspp)
for( s in 1:nspp){ 
	kc[,s] = b_rr[s]/2*exp(-b_rr[s]*abs(xx0))
	kc[,s]=kc[,s]/(sum(kc[,s]))
	fkc[,s]=fft(kc[,s])#/(np+1)
}
	

#=============================================================================
# Key variables for output: the invasion growth rates, their components
#=============================================================================
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
y.full=NULL
w.eq.full=NULL

###################################################################
#These are for each of the pairwise scenarios. Same variables, but for 
#1 vs 2, 1 vs 3, etc. 

#Simulation IGR
ldg.sim.p = matrix(0,ngenst, nspp*(nspp-1))
ldgIMP.sim.p = matrix(0,ngenst, nspp*(nspp-1))
covIMP.sim.p= matrix(0,ngenst, nspp*(nspp-1))

#Components of IGR
l1.p=matrix(0,ngenst,nspp*(nspp-1))
y.p=matrix(0,ngenst,nspp*(nspp-1))
w.eq.p=matrix(0,ngenst,nspp*(nspp-1))
D.p = matrix(0,ngenst,nspp*(nspp-1))
var_mu_Us.p=matrix(0,ngenst,nspp*(nspp-1)) 
cov_e_mu_Us.p=matrix(0,ngenst,nspp*(nspp-1))
cov_lam_vc.p=matrix(0,ngenst,nspp*(nspp-1))
cov_lam_vc.p2=matrix(0,ngenst,nspp*(nspp-1))
Elam1.p=matrix(0,ngenst,nspp*(nspp-1)) 
Elam2.p=matrix(0,ngenst,nspp*(nspp-1))
gr1.n.p=matrix(0,ngenst,nspp*(nspp-1))
gr1.p=matrix(0,ngenst,nspp*(nspp-1))

#=============================================================================
#Outermost loop: treat each species as invader
#=============================================================================

#Keep track of total interspecific competition of species when it is invader
cc=array(c(matrix(0,ngenst,np),matrix(0,ngenst,np)),dim=c(ngenst,np,nspp)) 

for ( s in 1:nspp) { 

	s.index= 1:nspp
	s.index = s.index[-s]

	#Equilibrium populations of residents. Columns represent space, rows are time. 
	nr=array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp)) 
	nr.p=array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp)) 
	ni.p=array(c(matrix(0.0,ngenst,np),matrix(0.0,ngenst,np)),dim=c(ngenst,np,nspp)) 
	
	#Theoretical pairwise equilibrium of residents
	nr2_eq_all = array(c(matrix(0,ngenst,np),matrix(0,ngenst,np)),dim=c(ngenst,np,(nspp))) 

	#Initialize residents
	for ( sa in s.index) { nr[1,,sa] = matrix(0.01,1,np)}

	
	#=============================================================================
	#Invader and resident stationary distributions
	#=============================================================================


	for( t in 1: ngenst){ 

		#If the "FAST" scenario is chosen, need to do a special set of 
		#calculations of the invader.
		if(fast.bys==TRUE ){   
			#Declare new variables for range 
	

			if (t==(iconfig)){
				#Remake the range of fast species
				if ( sum(peaks[1]-peaks.end[1]) != 0 ) { 
				pks = c( peaks[1], peaks.end[1], peaks.by[1]) } else {
				pks = peaks[1]} 

				if ( sum(Ds[1]-Ds.end[1]) != 0) { 
				Ds.tmp = c( Ds[1], Ds.end[1], Ds.by[1]) } else {
				Ds.tmp = Ds[1]} 
				
				Frs[,,1] = make.range(Fr[1], pks, Ds.tmp, stc, ngens, iconfig, fast.bys,dnwidth)
				#Declare new variables for range 
				peak.new=which(nr[(t-1),,s] == max(nr[(t-1),,s]),arr.ind=T)-floor(np/2)
				Fr.inv=Frs[,,1]		
				print(t)
			} 

			if (t >(iconfig)){
				#Calculate the IGR for species 1 for only this time step
				gr1.fast=get.fast.igr(Frs[(t-1),,],nr[(t-1),,], sr, alphas, kd,kc,n.inv=1 )
				#Now use the IGR to calculate the spread rate
				cs_all = get.spread.rate(gr1.fast,a_rr,sr)	

				# Make the new intrinsic fitness distribution for the next timestep
				# First, make a "pretend" fitness distribution based on how far 
				# the population of the invader can spread in a single time step
				Fr.inv[t,] = get.fast.peak(peak.new,Fr,Ds,cs_all,stc[(t-1),,],np)
				# Then filter the actual intrinsic range according to the amount 
				# of overlap
				Frs[t,,1] = Frs[t,,1]*as.numeric(Fr.inv[t,]>1e-2) 
				#Frs[t,,1] = pmin(Frs[t,,1],Fr.inv[t,])
				peak.new=peak.new+cs_all
				print(t)

			}
		}
		
		
		#Get the equilibrium of the resident community when invader is absent
		nr[t,,(s.index)] = get.res.eq(Frs[t,,],s.index,sr,alphas, fkd,kc, fkd.yes,fast=TRUE, burns = 5 )
		for (sp in 1:(nspp-1)){ 
			#Get the pairwise resident equilibrium
			nr.p[t,,s.index[sp]] = get.res.eq(Frs[t,,],s.index[sp],sr,alphas, fkd,kc,fkd.yes,fast=TRUE, burns = 5)
		}

		#This is the analytical version of each resident's equilibrium density when it is alone. 
		# for( sa in 1:(length(s.index)) ) {
		# 	sb=s.index[sa]
		# 	nr2_eq_all[t,,sb] = get.rsd.pair(Frs[t,,sb],sr[sb], alphas[s,sb],kd[,sb],kc[,sb])
		# 	}


			#Get the low-density equilibrium density of the invader against the resident community
			nr[t,,s] = get.inv.ldeq(Frs[t,,], nr[t,,], s, sr,alphas, fkd, kc,fkd.yes,fast=TRUE, burns = 5 )
			for (sp in 1:(nspp-1)){ 
				if(s<s.index[sp]) {s.inv=1}else{s.inv=2}
				#Get the pairwise invader low-density equilibrium 
				ni.p[t,,s.index[sp]] =get.inv.ldeq(Frs[t,,sort(c(s,s.index[sp]))], nr.p[t,,sort(c(s,s.index[sp]))], s.inv, sr,alphas, fkd,kc,fkd.yes,fast=TRUE, burns = 5 )
			}


	print(t)
	} 

	#=============================================================================
	# Low-density growth rates -- using simulation data, all 3 spp
	#=============================================================================
	#Simple numerical LDG: use the invader low-density equilibrium density, 
	#then let it invade against the community: 
	
	for (t in 1:ngenst) { 
		inv.one =pop_lg(Frs[t,,],nr[t,,], sr, alphas, fkd,kc, fkd.yes )[,s]
		ldg.sim[t,s]=mean(inv.one)/mean(nr[t,,s])

	}

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


	for (t in 1:ngenst) { 
		inv.oneIMP =pop_lg(Frs[t,,], nr2[t,,], sr, alphas, kdg,kc )[,s]
		ldgIMP.sim[t,s]=mean(inv.oneIMP)/mean(nr2[t,,s])

	}

	covIMP.sim[,s] = ldg.sim[,s]-ldgIMP.sim[,s]

	#Calculate the LDG for species vs. community by component

	y = matrix(0,ngenst, nspp-1)
	w.eq = matrix(0,ngenst, nspp-1)

	#Calculate standard spatial terms
	for (t in 1:ngenst) { 
		
		sp.ids =c(s, s.index) #IDs of resident, and invaders (this way uses all species)
		nrs= nr[t,,s.index] #Residents 
		
		#isd.tmp = get.isd(Frs[t,,s], nrs,s, alphas, sr[s], kc, kd[,s])

		y[t,] = colMeans(nrs,na.rm=T) #Means of resident densities 
		w.eq[t,] = alphas[s,s.index]*y[t,] #Weights for the LDG

		l1[t,s]=mean(Frs[t,,s])
		D[t,s]= 1+sum(w.eq[t,])


		#Calculate the standardized competition from residents
		muj=nrs/(matrix(y[t,],np,( nspp-1),byrow=T))-1
		uijmuj = matrix(0,np,nspp-1)
		for(a in 1:(nspp-1)) {
				ns = s.index[a]
				uijmuj[,a] = convolve(muj[,a],kc[,ns])
				uijmuj[,a]= c(uijmuj[ceiling(np/2):(np),a], uijmuj[1:floor(np/2),a] )
				uijmuj[,a] = w.eq[t,a]*(uijmuj[,a]) #+abs(min(uijmuj[,a])))

		}

		uijmuj[uijmuj<0]=0
		#uijmuj = uijmuj*(matrix(c(w.eq[t,]/D[t,s]),np,nspp-1,byrow=T))

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

	} #End spatial terms

	plot(gr1.n[,s],ylim= c(0.5,2))
	points(ldgIMP.sim[,s],col="red")

	plot(gr1[,s],ylim= c(0.5,2))
	points(ldg.sim[,s],col="red")

	y.full=c(y.full, list(y) )
	w.eq.full=c(w.eq.full,list(w.eq))
	#=============================================================================
	# Low-density growth rates for pairwise cases!  -- using theoretical 2spp equilibrium
	#=============================================================================


	#Outer loop through residents 

	for (rs in 1:(nspp-1)){
		#These are all of the necessary terms to calculate the standard invasion growth rates

		enter1 = (nspp-1)*(s-1)+rs
		rez.no = s.index[rs]

		#Calculate standard spatial terms
		for (t in 1:ngenst) { 

			#Simple numerical pairwise LDG: use the invader low-density equilibrium density, 
			#then let it invade against the community: 
			nr.inv = matrix(0,np,2)
			nr.inv[,1] = ni.p[t,,rez.no] #Invader low-density distribution
			nr.inv[,2] = nr.p[t,,rez.no] #Resident stationary distribution
			inv.one.p =pop_lg(Frs[t,,c(s,rez.no)],nr.inv,sr, alphas,fkd, kc, fkd.yes)[,1]
			ldg.sim.p[t,enter1]=mean(inv.one.p)/mean(nr.inv[,1])

			#Calculate the spatially implicit portion and the fitness-density covariance
			#separately, still using simulations. This just uses a "global" dispersal kernel
			####Dispersal kernels and their Fourier transforms 
			#Make the invader spatially homogenous, at low density
			nr.inv[,1] = matrix(inv.id,np) 
			inv.oneIMP.p =pop_lg(Frs[t,,c(s,rez.no)],nr.inv,sr, alphas, kdg,kc )[,1]
			ldgIMP.sim.p[t,enter1]=mean(inv.oneIMP.p)/mean(nr.inv[,1])

			#pairwise fd-covariance from simulations
			covIMP.sim.p[t,enter1] = ldg.sim.p[t,enter1]-ldgIMP.sim.p[t,enter1]

			l1.p[t,enter1]=mean(Frs[t,,s])
			y.p[t,enter1]=mean(nr2_eq_all[t,,rez.no])#Means of resident densities 
			w.eq.p[t,enter1]=alphas[s,rez.no]*y.p[t,enter1]#Weights for the LDG
			D.p[t,enter1]= 1+sum(w.eq.p[t,enter1])

			#Calculate the standardized competition from residents
			nrs.p=Re(fft((alphas[s,rez.no]*fft(nr2_eq_all[t,,rez.no])*fkc[,rez.no]),inverse=T)/(np+1))
			nrs.p=c(nrs.p[ceiling(np/2):(np)],nrs.p[1:floor(np/2)])
		
			#Calculate the standardized competition from residents
			muj.p = nr2_eq_all[t,,rez.no]/y.p[t,enter1]-1
			uijmuj.p = convolve(muj.p,kc[,rez.no])
			uijmuj.p= c(uijmuj.p[ceiling(np/2):(np)], uijmuj.p[1:floor(np/2)] )
			uijmuj.p = w.eq.p[t,enter1]*(uijmuj.p)#+abs(min(uijmuj.p)))	
			#uijmuj.p[uijmuj.p<0]=0

			cc.p=uijmuj.p

			#Non-linear competitive variance and the spatial storage effect
			#These terms represent the perturbation terms from Snyder 2008 
			var_mu_Us.p[t,enter1] = var(cc.p)
			cov_e_mu_Us.p[t,enter1] = cov(Frs[t,,s]/mean(Frs[t,,s])-1, cc.p)
			
			#First-order approximation
			Elam1.p[t,enter1]=l1.p[t,enter1]/(D.p[t,enter1])*( 1+1/(D.p[t,enter1])^2*var_mu_Us.p[t,enter1]-
					(1/D.p[t,enter1])*cov_e_mu_Us.p[t,enter1] )+sr[s]-1
			#Second-order approximation
			Elam2.p[t,enter1]=(l1.p[t,enter1]*(1/(D.p[t,enter1])) +
				sr[s]-1)^2+2*(l1.p[t,enter1]*(1/(D.p[t,enter1])) +
				sr[s]-1)*(var_mu_Us.p[t,enter1]-cov_e_mu_Us.p[t,enter1])/(D.p[t,enter1])^4
					
			gr1.n.p[t,enter1] = exp(Elam1.p[t,enter1]-0.5*Elam2.p[t,enter1])
			
			#The fitness-density covariance 
			#Method 1: 
			# tryCatch( {cov_lam_vc.p[t,enter1]=get.fd.cov.p(Frs[t,,s], nr2_eq_all[t,,rez.no], 
			# 	sr[s], alphas[s,rez.no], a_rr[s], kc[,s])} , error=function(e){} )
			# #Method 2: 
			# #calculate the invader stationary distribution: 
			# cov_lam_vc.p2[t,enter1]=get.fd.cov.p2(Frs[t,,s],
			# 	nr2_eq_all[t,,rez.no], sr[s], alphas[s,rez.no], kd[,s], kc[,s])
			

			#The full LDG
			#gr1.p[t,enter1] = gr1.n.p[t,enter1]+cov_lam_vc.p[t,enter1]
			gr1.p[t,enter1] = gr1.n.p[t,enter1]+covIMP.sim.p[t,enter1]


		} #End spatial terms


	}#End outer species loop through residents
 print(s)

} #End main loop through species

#=============================================================================
#Data saving
#=============================================================================


#Save all of the meaningful invasion growth rate variables
#in a file with an informative but awkwardly long file name. 
#This includes information about the temporal and spatial extent,
#and the details of the scenario type. 

file.name = (paste(f.name1, "_waldmean,var",sep=""))
save(file=file.name, "l1", "D", "var_mu_Us","cov_e_mu_Us",
"cov_lam_vc", "cov_lam_vc2", "Elam1", "Elam2", "gr1.n", "gr1", "y.full", "w.eq.full", "l1.p",
"y.p", "w.eq.p", "D.p", "var_mu_Us.p", "cov_e_mu_Us.p","cov_lam_vc.p", "cov_lam_vc.p2",
"Elam1.p","Elam2.p", "gr1.n.p", "gr1.p","ldg.sim","ldgIMP.sim","covIMP.sim","ldg.sim.p","ldgIMP.sim.p","covIMP.sim.p","Frs")


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