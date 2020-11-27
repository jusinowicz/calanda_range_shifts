#=============================================================================
# This code contains basic functions to calculate species intrinsic ranges
# and determine components of the invasion growth rates for species.
#
# The model is assumed to be the spatially explicit Leslie-Gower model, where
# the intrinsic rate of reproduction is assumed to be the spatially-varying
# parameter. 
#
# Dispersal kernels are allowed to vary in space and be location-specific
#
# The functions are broken into 3 sections: 
# 1. Miscellanerous functions: Mostly helper functions for later routines.
#		Also some functions to set up theoretical intrinsic ranges for 
#		simulations. 
# 2. Invader and resident equilibrium densities: Functions to calculate 
#		components of the Low-density growth rates (LGR) in various ways
#		(including both analytical and numerical approaches).
#		Two subsections: 1) Multispecies (full community) case and 2) Pairwise
#		case. 
# 3. Environmental distance and site-LGR impact: Functions to determine 1)
#		1) the environmental distance between a site's current and future state
#		(environmental distance is the Euclidean distance in the multivariate
#		abiotic space), and its environmental/physical distance from the 
#		nearest analogue and, 
#		2) the incluence that each site has on the various LGR metrics of (2)
#
#=============================================================================


#=============================================================================
# FUNCTION DEFINITIONS
#=============================================================================

#=============================================================================
#Miscellaneous functions
#=============================================================================

#R equivalent for matlab's meshgrid
meshgrid=function(a,b) {
  list(
       x=outer(b*0,a,FUN="+"),
       y=outer(b,a*0,FUN="+")
       )
} 
#These two functions can be used to find upper and lower edges
#of the range distributions. Needs the peak position, variance
#of the distribution, and a tolerance specifying number
#of standard deviations to include. 
get.upper = function (np, burns, ngens, pks, Ds, D.tol= 3) {
	
	ngenst=burns+ngens
	#Test if the peak and spatial width parameters are a single value
	#(i.e. stationary or moving)
	if(length(pks)>1){
		pk1 = pks[1]
		pk.end = pks[2]
		pk.by = pks[3]
		peak.stc.burn=c(matrix(pk1, ngens,1))
		peak.stc.gens=seq(pk1, pk.end, by=pk.by)[1:burns]
		peak.stc = c(peak.stc.burn, peak.stc.gens)} else {
		pk1 = pks
		peak.stc=c(matrix(pk1, ngenst,1))
	}


	if(length(Ds)>1){
		Ds1 = Ds[1]
		Ds.end = Ds[2]
		Ds.by = Ds[3]
		Ds.stc.burn=c(matrix(Ds1, ngens,1))
		Ds.stc.gens=seq(Ds1, Ds.end, by=Ds.by)[1:ngens]
		Ds.stc= c(Ds.stc.burn, Ds.stc.gens)
		} else { 
				
		Ds1 = Ds
		Ds.stc=c(matrix(Ds1, ngens,1))
	} 
	
		sds.up=round(np/2)+peak.stc+sqrt(Ds.stc)*D.tol
		return(sds.up)

}

get.lower = function (np, burns, ngens, pks, Ds, D.tol= 3) {
	
	ngenst=burns+ngens
	#Test if the peak and spatial width parameters are a single value
	#(i.e. stationary or moving)
	if(length(pks)>1){
		pk1 = pks[1]
		pk.end = pks[2]
		pk.by = pks[3]
		peak.stc.burn=c(matrix(pk1, ngens,1))
		peak.stc.gens=seq(pk1, pk.end, by=pk.by)[1:burns]
		peak.stc = c(peak.stc.burn, peak.stc.gens)} else {
		pk1 = pks
		peak.stc=c(matrix(pk1, ngenst,1))
	}


	if(length(Ds)>1){
		Ds1 = Ds[1]
		Ds.end = Ds[2]
		Ds.by = Ds[3]
		Ds.stc.burn=c(matrix(Ds1, ngens,1))
		Ds.stc.gens=seq(Ds1, Ds.end, by=Ds.by)[1:ngens]
		Ds.stc= c(Ds.stc.burn, Ds.stc.gens)

		} else { 
				
		Ds1 = Ds
		Ds.stc=c(matrix(Ds1, ngens,1))
	} 
	
		sds.low=round(np/2)+peak.stc-sqrt(Ds.stc)*D.tol
		return(sds.low)

}


#=============================================================================
#Ranges and range shifts 
#=============================================================================
#Make a species' intrinsic range using a Gaussian distribution. 
#Variables correspond to: 
#	Fr	Max height of distribution (max reproduction)
#	pks	Vector with: 
#		 1. Initial spatial mean (position)
#		 2. pk.ends: Ending spatial mean
#		 3. by: Rate at which mean moves 
#		If not a vector, it is assumed stationary (no change in position) 
#
#	Ds	Vector with:
#		 1. Spatial variance (width)
#		 2. Ds.end: ending spatial variance
#		 3. Rate at which variance changes
#		If not a vector, it is assumed stationary (no change in position) 
#
#	stc	spatial array with space and time coordinates
#	burns	time period where distribution is stationary
#	ngens	time period of change for distribution
#=============================================================================

make.range=function(Fr, pks, Ds, stc, ngens, burns, fast.by = FALSE, dnwidth=NULL) { 


np=dim(stc)[2]
nt=dim(stc)[1]
Fr.tmp = matrix(0, nt, np)

#Test if the peak and spatial width parameters are a single value
if(length(pks)>1){pk1 = pks[1]} else { 
			pk1 = pks} 

if(length(Ds)>1){Ds1 = Ds[1]} else { 
			Ds1 = Ds} 

#Burn phase
X1=stc[1:burns,(1:(np)),2]
G1D=exp(-((X1-pk1)^2)/(2*Ds1))
G1D[G1D<1e-7]=0
G1D=G1D/(max(G1D))*Fr
Fr.tmp[1:burns,]=G1D

#Test if it is the fast or normal scenario
# Note: fast scenario inherently assumes changing peak position
# BUT NO change in width. 

if(fast.by==FALSE) {
	#Moving phase
	#Test if mean and variance change and create appropriate variables if
	#they do
	if(length(pks)>1){ pk.end = pks[2]
				pk.by = pks[3]
				peak.stc=seq(pk1, pk.end, by=pk.by)[1:ngens]} else {
				peak.stc=c(matrix(pk1, ngens,1))}
	 

	if(length(Ds)>1){  Ds.end = Ds[2]
				Ds.by = Ds[3]
				Ds.stc=seq(Ds1, Ds.end, by=Ds.by)[1:ngens]} else {
				Ds.stc=c(matrix(Ds1, ngens,1))}
	  

	X2=stc[(1+burns):nt,(1:(np)),2]
	G1D=exp(-((X2-peak.stc)^2)/(2*Ds.stc))
	G1D[G1D<1e-7]=0
	G1D=G1D/(max(G1D))*Fr
	Fr.tmp[(1+burns):nt,]=G1D
} else { 

	pk.end = pks[2]
	pk.by = pks[3]
	peak.stc1=seq(pk1, pk.end, by=pk.by)[1:dnwidth]
	peak.stc2=matrix(pk.end, (ngens-dnwidth),1)
	peak.stc=c(peak.stc1,peak.stc2)

	Ds.stc=c(matrix(Ds1, ngens,1))

	X2=stc[(1+burns):nt,(1:(np)),2]
	G1D=exp(-((X2-peak.stc)^2)/(2*Ds.stc))
	G1D[G1D<1e-7]=0
	G1D=G1D/(max(G1D))*Fr
	Fr.tmp[(1+burns):nt,]=G1D

}


return(Fr.tmp)

}

#=============================================================================
#Calculate the population spread rate from an IGR. This is based on
#the math from Neubert and Caswell 00, but originating with Weinberger 78, 
#and Kot 92. 
#=============================================================================

get.spread.rate = function(gr1.fast,a_rr,sr) {
	u=seq(0,.1,0.0001)
	#For Lalpacian kernel 
	r1=gr1.fast*1/(1-(1/a_rr[1])^2*u^2) 
	r2=sr[1]*1
	cs=1/u*log(r1+r2)
	cs_all=min(cs,na.rm=T)/10
	
	return(cs_all)
	}


#=============================================================================
#Calculate a new fitness distribution for a single time step. 
#=============================================================================

get.fast.peak = function(peak.new,Fr,Ds,cs_all,stc,np){
	peak.new=peak.new+cs_all 
	max.new=Fr[1]

	X2.new=stc[(1:(np)),2]

	G1D.new=exp(-((X2.new-peak.new)^2)/(2*Ds[1]));
	G1D.new[G1D.new<1e-7]=0
	G1D.new=G1D.new/(max(G1D.new))*max.new

	return(G1D.new)

	}

#=============================================================================
#Population dynamics
#=============================================================================
#One step of the spatial Leslie-Gower model for N species
#1D array of sites. 
#It takes as variables: 
#	Fr.spp 	The intrinsic ranges of all species
#	nr.spp 	Population matrix with current populations or initial conditions
#	sr.spp 	Species survival
#	alpha.spp Competition matrix
#	kd.spp 	Dispersal kernels
#	fkd.yes If this is yes, then it is the Fourier transformed kernels being
#			passed to the function
#	kc.spp 	Competition kernels
#=============================================================================

pop_lg = function (Frs.spp,nr.spp, sr.spp,alpha.spp, kd.spp, kc.spp,fkd.yes=FALSE){

	Frs.spp = as.matrix(Frs.spp)
	np = dim(as.matrix(Frs.spp))[1] #Space
	nspp = dim(as.matrix(Frs.spp))[2] #Number of species
	

	#Check the dimensionality of the dispersal kernel. If it is a spatially 
	#heteroeneous dispersal kernel, then set het_kern=1

	het_kern=0

	if (!is.na(dim(kd.spp)[3])){ het_kern=1 } 

	#When fkd.yes=FALSE: 
	#Transform dispersal and competition kernels
	if( fkd.yes == FALSE){ 
		if(het_kern==1){
			fkd = kd.spp*0

			if(dim(kd.spp)[2]<2){ 
				fkc=fft(as.matrix(kc.spp))
				
				for(ss in 1:nspp){
					fkd[,,ss]=fft(as.matrix(kd.spp[,,ss]))
				}

			}else{
				fkc=mvfft(as.matrix(kc.spp))

				for(ss in 1:nspp){
					fkd[,,ss]=mvfft(as.matrix(kd.spp[,,ss]))
				}
			}

		}else{ 

			if(dim(as.matrix(kd.spp))[2]<2){ 
				fkd=fft(as.matrix(kd.spp))
				fkc=fft(as.matrix(kc.spp))
				
			}else{

				fkd=mvfft(as.matrix(kd.spp))
				fkc=mvfft(as.matrix(kc.spp))
			}
		}
	} else {

		fkd = kd.spp
		fkc = kc.spp
	}


	nr.burns =array(c(matrix(0.0,2,np),matrix(0.0,2,np)),dim=c(2,np,nspp)) 
	nr.burns[1,,] = nr.spp
	#Seedlings germinate and produce seeds at rate Fr, weighted by competition
	if(dim(as.matrix(nr.burns[1,,]))[2]<2){
		fp = as.matrix(fft(nr.burns[1,,]))
	}else{ 
		fp = as.matrix(mvfft(nr.burns[1,,]))
	}

	#Total competition for each resident species 
	cr=matrix(0,np,nspp^2)
	cr_tot = matrix(0,np,nspp)
		for (sa in 1:nspp) { 
			for (sb in 1:nspp) {
				

					cr[,(sa*(sa-1)+sb)] = Re(fft((alpha.spp[sa,sb]*fp[,sb]*fkc[,sb]),
						inverse=T)/(np+1))

					#This if/else makes the code more robust to odd/even spatial extent
					if(ceiling(np/2) == floor (np/2) ){
						cr[,(sa*(sa-1)+sb)] = c(cr[((ceiling(np/2)+1):(np)),(sa*(sa-1)+sb)], 
							cr[(1:floor(np/2)),(sa*(sa-1)+sb)])
					}else {
						cr[,(sa*(sa-1)+sb)] = c(cr[(ceiling(np/2):(np)),(sa*(sa-1)+sb)], 
							cr[(1:floor(np/2)),(sa*(sa-1)+sb)])
					}

					cr_tot[,sa] = cr_tot[,sa] + cr[ ,(sa*(sa-1)+sb)]	
			}
		}

	#competition-weighted seed production
	lam = matrix(0,np,nspp)

	#Seed dispersal
	fd = matrix(0,np,nspp)
	fd2= array(matrix(0,np,np),dim=c(np,np,nspp) )
	nr_disp = matrix(0,np,nspp)

	for(sa in 1:nspp) { 
		#competition-weighted seed production (Leslie-Gower)
		lam[,sa]=(Frs.spp[,sa])/(1+cr_tot[,sa])
		lam[,sa]= lam[,sa]*(nr.burns[1,,sa]>1e-34)
		#If both species have zero germination, then NAs are produced: 
		lam[,sa][is.na(lam[,sa]) ] = 0

		#Seeds disperse

		if(het_kern==1){
			#The most robust way to do this is to disperse each site individually 
			#according to its numerical dispersal kernel
			lam.tmp = lam[,sa]*nr.burns[1,,sa]
			lam2 = matrix(0,np,np)
			#lam2[ceiling(np/2),] = lam2[ceiling(np/2),]+t(as.matrix(lam.tmp))
			diag(lam2) = as.matrix(lam.tmp)
			fd2[,,sa]=mvfft(lam2)
			nr_disp[,sa] = rowSums( Re(mvfft( (fd2[,,sa]*fkd[,,sa]), inverse=T)/(np+1)),na.rm=T	)

		} else { 

			fd[,sa]=fft(lam[,sa]*nr.burns[1,,sa])
			nr_disp[,sa] = Re(fft( (fd[,sa]*fkd[,sa]), inverse=T)/(np+1))
		}

		#This if/else makes the code more robust to odd/even spatial extent
		if(ceiling(np/2) == floor (np/2) ){
			nr_disp[,sa] = c(nr_disp[((ceiling(np/2)+1):(np)),sa],nr_disp[(1:floor(np/2)),sa])
		}else {
			nr_disp[,sa] = c(nr_disp[(ceiling(np/2):(np)),sa],nr_disp[(1:floor(np/2)),sa])
		}
		
		#Final population step
		nr.burns[2,,sa] = nr_disp[,sa] + nr.burns[1,,sa]*sr[sa]

	}

	return(nr.burns[2,,])
}


#=============================================================================
#Invader and resident equilibrium densities
#=============================================================================

#=======================#
#Multi-species routines
#=======================#

#=============================================================================
#Get the resident equilibrium densities using simulations. This version does
#not require any analytical calculations and works for an arbitrary number 
#of species. 
#1D array of sites. 
#It takes as variables: 
#	Fr.spp 	The intrinsic ranges of all species
#	s.index The identities of the residents
#	sr.spp 	Species survival
#	alpha.spp Competition matrix
#	kd.spp 	Dispersal kernels
#	fkd.yes If this is yes, then it is the Fourier transformed kernels being
#			passed to the function
#	kc.spp 	Competition kernels
#
# 	burns	Number of time steps at which to test equilibrium: default=500
#	burns.max Maximum number of iterations (in multiples of burns)
#			to test for equilibrium. Default = 10
#	tol 	The tolerance for allowing equilibrium, in terms of significance of
#			fit of a line through the last burns/2 points. Default alpha = 0.05
#
#	fast 	If TRUE, the check for equilibrium is switched to a different routine
#			which is much faster, but potentially less accurate. 
#	fast.tol Test for equilibrium: how close is slope to 0? Default = 1e-5
#=============================================================================

#
get.res.eq = function(Frs.spp,s.index,sr.spp,alpha.spp,kd.spp,  kc.spp,fkd.yes =FALSE,burns=500, 
	burns.max = 10,tol=1,fast=FALSE,fast.tol=1e-5){

	np = dim(as.matrix(Frs.spp))[1] #Space
	nspp = length(s.index) #Number of residents
	
	Frs.spp = Frs.spp[,s.index]
	sr.spp = sr.spp[s.index]
	#alpha.spp = alpha.spp[s.index,s.index]

	nr.burns=array(c(matrix(0.0,burns,np),matrix(0.0,burns,np)),dim=c(burns,np,nspp)) 
	#Initialize residents
	for ( sa in 1:nspp) { nr.burns[1,,sa] = matrix(0.01,1,np)}
	b=1
	#Test to see whether the maximum number of iterations has been reached
	while (b <= burns.max){ 
		#Start each burn iteration
		for (n in 1:(burns-1)) {
		
			nr.burns[(n+1),,] = pop_lg(Frs.spp ,nr.burns[n,,], sr.spp,alpha.spp,kd.spp,  kc.spp, fkd.yes)
			print(paste ("first", n, sep=" "))
		}

		#Test for equilibrium for each specie: has the equilibrium density stayed unchanged,
		#according to the significance of a linear fit through the last burns/2 points? 
		
		if(fast==FALSE) { 
			is.eq = matrix(0,nspp,1)
			tt=1:(burns/2+1)
			for(sa in 1:nspp) {
				nr.tmp = rowMeans(nr.burns[(burns/2):burns,,sa])
				sp.lm = lm(nr.tmp~tt)
				if( summary(sp.lm)$coefficients[2,4] < tol ) {is.eq[sa] = 1}
			} 
		}else { 

			is.eq = matrix(0,nspp,1)
			tt=1:(burns/2+1)
			for(sa in 1:nspp) {
				nr.tmp = rowMeans(nr.burns[(burns/2):burns,,sa])
				sp.lm = lm(nr.tmp~tt)
				if( sp.lm$coefficients[2] > fast.tol ) {is.eq[sa] = 1}
			}
		}
		
		#Test whether all species have reached equilibrium. If yes, quit and return
		#nr.burns. Otherwise, continue. 
		if( sum(is.eq[sa] == 0) ) {b = burns.max+1}else{ 
			nr.burns[1,,] =nr.burns[burns,,]
			b=b+1
			}
			print(paste ("outer", n, sep=" "))

	}

	#If this quit due to equilibrium, return equilibriums. Otherwise, return an error
	if (b > burns.max) { return(nr.burns[burns,,])} else {
			return( print("Error: Max iterations exceded without reaching equilibrium"))}

}
#============================================================================
#Get the invader's low-density equilibrium density using simulations. This version does
#not require any analytical calculations and works for an arbitrary number 
#of species. 
#1D array of sites. 
#It takes as variables: 
#	Fr.inv	The intrinsic range of the invader
# 	nr.res 	The residents' stationary distributions
#	s.inv 	The identity of the invader
#	sr.inv 	Invader survival
#	alpha.spp Competition matrix
#	kd.spp 	Dispersal kernels
#	fkd.yes If this is yes, then it is the Fourier transformed kernels being
#			passed to the function
#	kc.spp 	Competition kernels
#
# 	burns	Number of time steps at which to test equilibrium: default=500
#	burns.max Maximum number of iterations (in multiples of burns)
#			to test for equilibrium. Default = 10
#	tol 	The tolerance for allowing equilibrium, in terms of significance of
#			fit of a line through the last burns/2 points. Default alpha = 0.05
#
#	fast 	If TRUE, the check for equilibrium is switched to a different routine
#			which is much faster, but potentially less accurate. 
#	fast.tol Test for equilibrium: how close is slope to 0? Default = 1e-5
#=============================================================================

get.inv.ldeq = function(Frs.spp, nr.res,s.inv,sr.spp,alpha.spp,kd.spp,  kc.spp,fkd.yes =FALSE,burns=500, burns.max = 10,tol=1,fast=FALSE,fast.tol=1e-5){

	np = dim(as.matrix(Frs.spp))[1] #Space
	nspp = dim(as.matrix(Frs.spp))[2] #Number of species
	s.index= 1:nspp
	s.index = s.index[-s.inv]

	inv.id =1e-5

	nr.burns=array(c(matrix(0.0,burns,np),matrix(0.0,burns,np)),dim=c(burns,np,nspp)) 
	
	#Initialize invader and residents
	nr.burns[1,,s.inv] = matrix(inv.id,1,np)
	nr.burns[1,,s.index] = nr.res[,s.index]

	b=1
	#Test to see whether the maximum number of iterations has been reached
	while (b <= burns.max){ 
		#Start each burn iteration
		for (n in 1:(burns-1)) {
			nr.burns[(n+1),,s.inv] = pop_lg(Frs.spp ,nr.burns[n,,], sr.spp,alpha.spp,kd.spp, kc.spp, fkd.yes)[,s.inv]


			#Reset invader to low density
			nr.burns[(n+1),,s.inv] = (inv.id /sum(nr.burns[(n+1),,s.inv])) * nr.burns[(n+1),,s.inv]

		}

		#Test for equilibrium for invader: has the density stayed unchanged,
		#according to the significance of a linear fit through the last burns/2 points? 

		is.eq =0
		tt=1:(burns/2+1)
		if(fast==FALSE) { 
			
			nr.tmp = rowMeans(nr.burns[(burns/2):burns,,s.inv])
			sp.lm = lm(nr.tmp~tt)
			if( summary(sp.lm)$coefficients[2,4] < tol ) {is.eq = 1}
			
		}else { 
	
			nr.tmp = rowMeans(nr.burns[(burns/2):burns,,s.inv])
			sp.lm = lm(nr.tmp~tt)
			if( sp.lm$coefficients[2] > fast.tol ) {is.eq = 1}
			
		}
		
		#Test whether the invader has reached equilibrium. If yes, quit and return
		#nr.burns. Otherwise, continue. 
		if( sum(is.eq == 0) ) {b = burns.max+1}else{ 
			nr.burns[1,,s.inv] =nr.burns[burns,,s.inv]
			b=b+1
			}

	}

	#If this quit due to equilibrium, return equilibriums. Otherwise, return an error
	if (b > burns.max) { return(nr.burns[burns,,s.inv])} else {
			return( print("Error: Max iterations exceded without reaching equilibrium"))}

}
#=============================================================================
#Function to calculate the invader's stationary distribution in multispeces scenarios
#based on a combination of approximations and simulations. 
#1D array of sites. 
#It takes as variables: 
#	Fr.spp 	The intrinsic range of invader
#	nr.spp The standardized resident stationary distribution
#	sr.spp 	The survival of the invader.
#	alpha.spp Competition coefficients 
#	kd.spp 	Dispersal kernel
#	kc.spp 	Interspecific competition kernel
#=============================================================================
get.isd = function(Fr.spp,  nr.spp,s.inv, alpha.spp, sr.spp,kc.spp,kd.spp){
	np = length(Fr.spp)
	nres = dim(as.matrix(nr.spp))[2] #Number of species
	nspp=nres+1
	s.index= 1:nspp
	s.index = s.index[-s.inv]

	n.ave=colMeans(as.matrix(nr.spp))
	w.eq = alpha.spp[s.inv, s.index]*n.ave
	D.spp = 1+sum(w.eq)

	#Calculate the standardized competition from residents
	muj=nr.spp/(matrix(n.ave,np,( nres),byrow=T))-1
	uijmuj = matrix(0,np,nres)
		for(a in 1:(nres)) {
				ns = s.index[a]
				uijmuj[,a] = convolve(muj[,a],kc.spp[,ns])
				uijmuj[,a] = c(uijmuj[ceiling(np/2):(np),a], uijmuj[1:floor(np/2),a] )
				#Note that I am multiplying by the constants here
				uijmuj[,a] = w.eq[a]/D.spp*(uijmuj[,a])#+abs(min(uijmuj[,a])))

		}

	uijmuj[uijmuj<0]=0

	inv.id =1e-5
	burns=100

	fnr.its=matrix(0.0,burns,2*np+1) 
	nr.its=matrix(0.0,burns,np) 

	#Average LDG excluding survival
	lam1 = mean(Fr.spp)/(D.spp)
	#Transformed environmental response padded with 0s for fft
	eft= fft( c( (1:ceiling(np/2)*0), Fr.spp/mean(Fr.spp,na.rm=T)-1,(1:ceiling(np/2)*0) ) ) 
	fkd.spp = fft( c((1:ceiling(np/2)*0), kd.spp, (1:ceiling(np/2)*0) )) 
	#Transform resident stationary distribution
	fnr = fft(c((1:ceiling(np/2)*0), apply(uijmuj,1,sum), (1:ceiling(np/2)*0) ))
	
	# Be = lam1*fkd.spp
	# Bur = -lam1*fkd.spp
	# fnr.its [1,] = Be+Bur*fnr
	# for( b in 1:(burns-1)) {

	# 	Bu = (lam1*fkd.spp+sr[sa])
	# 	Be = lam1*fkd.spp
	# 	Bur = -lam1*fkd.spp
	# 	fnr.its[b+1,] = Bu*fnr.its[b,]+Be*eft+Bur*fnr
	# 	nr.its.tmp = Re(fft(fnr.its[b,], inverse=T))/(np+0.2*np)
	# 	nr.its[b+1,] =c(nr.its.tmp[ceiling(np+np/2+1):(np*2)], nr.its.tmp[1:ceiling(np/2)])
	# 	nr.its[b+1,] = nr.its[b+1,]+min(abs(nr.its[b+1,] ))
	# 	nr.its[b+1,][nr.its[b+1,]<0] = 0	
	# }

	#  return(nr.its)

	ir_tmp.f = 0

	#The iterative solution to the resident stationary distribution
	for( b in seq(10,1,-1)) { 
		

			Bu = (lam1*fkd.spp+sr.spp)^(b-1)
			Be = lam1*fkd.spp
			Bur = -lam1*fkd.spp
			ir_tmp.f = ir_tmp.f + Bu*(Be*eft+Bur*fnr)
		}


		ir_tmp = Re( fft(ir_tmp.f,inverse=T))/(np+0.2*np)
		
		ir_tmp= c(ir_tmp[ceiling(np+np/2+1):(np*2)], ir_tmp[1:ceiling(np/2)])
		#Some numerical adjustments -- move the minimum vale up from its current
		# value which is less than 0, and smooth out spurious edges created by fft
		ir_tmp = ir_tmp+min(abs(ir_tmp ))
		ir_tmp[ir_tmp<0] = 0

		return(ir_tmp)

}


#=============================================================================
#Pairwise routines
#=============================================================================
#Function to calculate the resident's stationary distribution in the pairwise scenarios
#based on analytical approximations. 
#1D array of sites. 
#It takes as variables: 
#	Fr.spp 	The intrinsic range of resident
#	sr.spp 	The survival of the resident.
#	alpha.rr Intraspecific competetion of resident 
#	kd.spp 	Dispersal kernel
#	kc.spp 	Intraspecific competition kernel
#=============================================================================

get.rsd.pair = function(Fr.spp,sr.spp,alpha.rr,kd.spp,kc.spp ){
	np = length(Fr.spp)

	#Theoretical average of resident density
	n.ave= mean(Fr.spp)/(1-sr.spp)-1
	#Theoretical scaling constant
	D.spp = 1+(alpha.rr*n.ave)
	#Average LDG excluding survival
	lam1 = mean(Fr.spp)/(D.spp)
	#Transformed environmental response padded with 0s for fft
	eft= fft( c( (1:ceiling(np/2)*0), Fr.spp/mean(Fr.spp,na.rm=T)-1,(1:ceiling(np/2)*0) ) ) 
	fkd.spp = fft( c((1:ceiling(np/2)*0), kd.spp, (1:ceiling(np/2)*0) )) 
	fkc.spp = fft( c((1:ceiling(np/2)*0), kc.spp, (1:ceiling(np/2)*0) )) 

	nr_tmp.f = lam1*eft
	Rs = fkd.spp/(1-(sr.spp + fkd.spp*(1-(alpha.rr*n.ave/D.spp)*fkc.spp)))
	nr_tmp.f = nr_tmp.f *sqrt(Re(Rs)^2 + Im(Rs)^2)

	#nr_tmp.f = 0	

	# for( b in seq(10,1,-1)) { 

	# 	Bu = (lam1*fkd.spp*(1-((alpha.rr*n.ave)/D.spp)*fkc.spp)+sr.spp)^(b-1)
	# 	Be = lam1*fkd.spp
	# 	nr_tmp.f = nr_tmp.f + Bu*Be*eft
	# }


	# #The iterative solution to the resident stationary distribution

	 	nr_tmp = Re( fft(nr_tmp.f,inverse=T))/(np+0.2*np)
		
		nr_tmp= c(nr_tmp[ceiling(np+np/2+1):(np*2)], nr_tmp[1:ceiling(np/2)])
	# 	#Some numerical adjustments -- move the minimum vale up from its current
	# 	# value which is less than 0, and smooth out spurious edges created by fft
	# 	nr_tmp = nr_tmp+min(abs(nr_tmp ))
	# 	nr_tmp[nr_tmp<0] = 0

		return( (nr_tmp-min(nr_tmp))/ mean(nr_tmp-min(nr_tmp) ))

}

#=============================================================================
#Function to calculate the invader's stationary distribution in the pairwise scenarios
#based on analytical approximations. 
#1D array of sites. 
#It takes as variables: 
#	Fr.spp 	The intrinsic range of invader
#	nr.spp The resident stationary distribution
#	sr.spp 	The survival of the invader.
#	alpha.ir Interspecific competetion of resident 
#	kd.spp 	Dispersal kernel
#	kc.spp 	Interspecific competition kernel
#=============================================================================

get.isd.pair = function(Fr.spp, nr.spp, sr.spp,alpha.ir,kd.spp,kc.spp ){
	np = length(Fr.spp)
	#Theoretical average of resident density
	#n.ave= mean(Fr.spp)/(1-sr.spp)-1
	n.ave = mean(nr.spp)
	#Theoretical scaling constant
	D.spp = 1+(alpha.ir*n.ave)
	#Average LDG excluding survival
	lam1 =mean(Fr.spp)/(D.spp)
	#Transformed environmental response padded with 0s for fft
	eft= fft( c( (1:ceiling(np/2)*0), Fr.spp/mean(Fr.spp,na.rm=T)-1,(1:ceiling(np/2)*0) ) ) 
	fkd.spp = fft( c((1:ceiling(np/2)*0), kd.spp, (1:ceiling(np/2)*0) )) 
	fkc.spp = fft( c((1:ceiling(np/2)*0), kc.spp, (1:ceiling(np/2)*0) )) 
	#Transform resident stationary distribution
	fnr = fft(c((1:ceiling(np/2)*0), nr.spp/n.ave-1, (1:ceiling(np/2)*0) ))

	# ir_tmp.f = lam1*(eft-(alpha.ir*n.ave/D.spp)*fkc.spp*fnr )
	# Rs = fkd.spp/(1-(sr.spp + lam1*fkd.spp) ) 
	# ir_tmp.f = ir_tmp.f *sqrt(Re(Rs)^2 + Im(Rs)^2)

	ir_tmp.f = lam1*(eft-((alpha.ir*n.ave)/D.spp)*fkc.spp*fnr )
	Rs = fkd.spp/(1-(sr.spp + lam1*fkd.spp) ) 
	ir_tmp.f = ir_tmp.f *sqrt(Re(Rs)^2 + Im(Rs)^2)
	#ir_tmp.f = ir_tmp.f *Rs



	#Second order:
	#ir_tmp.f = lam1*(1+eft-(1+eft)*(alpha.ir*n.ave/D.spp)*fkc.spp*fnr+ ((alpha.ir*n.ave/D.spp)*fkc.spp*fnr)^2 )
	#Rs = fkd.spp/(1-(sr.spp + lam1*fkd.spp*(eft+1 - (alpha.ir*n.ave/D.spp)*fkc.spp*fnr )) ) 
	#ir_tmp.f = ir_tmp.f *sqrt(Re(Rs)^2 + Im(Rs)^2)

	# ir_tmp.f = lam1*(1+eft-(1+eft)*(alpha.ir*n.ave/D.spp)*fkc.spp*fnr)
	# Rs = fkd.spp/(1-(sr.spp + lam1*fkd.spp*(eft+1) - (alpha.ir*n.ave/D.spp)*fkc.spp*fnr ) ) 
	# ir_tmp.f = ir_tmp.f *sqrt(Re(Rs)^2 + Im(Rs)^2)


	# ir_tmp.f = 0	
	# T=50

	# #The iterative solution to the invader stationary distribution
	# 	for( b in seq(1,T,1)) { 

	# 		Bu = (lam1*fkd.spp+sr.spp)^(b-1)
	# 		Be = lam1*fkd.spp
	# 		Bur = -lam1*fkd.spp*fkc.spp
	# 		ir_tmp.f = ir_tmp.f + Bu*(Be*eft+Bur*fnr)
	# 	}
		# for( b in seq(1,T,1)) { 
		# 	Bu = (lam1*fkd.spp+sr.spp)^(b-1)
		# 	ir_tmp.f = ir_tmp.f + Bu
		# }

		# Bu = (lam1*fkd.spp+sr.spp)
		# Be = lam1*fkd.spp
		# Bur = -lam1*fkd.spp*fkc.spp
		
		#ir_tmp.f= ir_tmp.f*(Be*eft+Bur*fnr)
		# #(alpha.ir*n.ave)/D.spp

		ir_tmp = Re( fft(ir_tmp.f,inverse=T))/(np+0.2*np)
		
		ir_tmp= c(ir_tmp[ceiling(np+np/2+1):(np*2)], ir_tmp[1:ceiling(np/2)])
		#Some numerical adjustments -- move the minimum vale up from its current
		# value which is less than 0, and smooth out spurious edges created by fft
		#ir_tmp = ir_tmp+min(abs(ir_tmp ))
		#ir_tmp[ir_tmp<0] = 0

		return(ir_tmp)

}

#=============================================================================
#Low-Density growth rates
#=============================================================================

#=======================#
#Multi-species routines
#=======================#

#=============================================================================
#This function is to get the fitness-density covariance for a single
#1D array of sites, for an arbitrary number of species. It requires 
#multi-species equilibriums calculated from e.g. simulations. 
#The analysis is based on Snyder 2008, Snyder and Chesson 2004,2003.
#It takes as variables: 
#	Fr.spp 	The intrinsic range of invader
#	sp.ids	An array indicating which species is the invader (entry 1),
#			and which species are residents (all subsequent entries).
#	nrs.spp The residents stationary distributions
#			-This should be a matrix with species as columns
#	sr.spp 	The survival of the invader.
#	kd.spp 	The invader dispersal kernel
#	fkd.yes If this is yes, then it is the Fourier transformed kernels being
#			passed to the function
#	alpha.ir Interspecific competetion with residents
#			- passed as the invader's row in competition matrix
#	a_rr 	Parameter for competition kernel
#	kc.spp 	Interspecific competition kernels with residents
#			-If this matrix is a single column then it is used
#			-used iteratively/identically across interactions	 
#=============================================================================

get.fd.cov = function (Fr.spp, sp.ids, nrs.spp, sr.spp, alpha.ir, a_rr, kd.spp, fkd.yes=FALSE, kc.spp){ 
	
	np=length(Fr.spp)
	inv.id = sp.ids[1] #Invader
	res.ids = sp.ids[2:length(sp.ids)] #Residents
	nspp = length (res.ids)

	#Check for the appropriate number of competition kernels
	if( ncol(kc.spp)<nspp | is.null(ncol(kc.spp)) == TRUE ) { kc.spp=matrix(kc.spp,np,(nspp+1))} 
	
	#Invader intrinsic fitness
	rmFr1=mean(Fr.spp) 
	ei1=(Fr.spp/rmFr1-1)

	#Resident competitive effects
	nrs.spp=nrs.spp[,res.ids]
	rmcc2=colMeans(nrs.spp,na.rm=T) 
	
	#Weight this by the equilibrium average alpha
	alpha.ir = alpha.ir[res.ids]
	w.eq = alpha.ir*rmcc2

	muj=nrs.spp/(matrix(rmcc2, np,nspp,byrow=T))-1
	uijmuj = matrix(0,np,nspp)
	for(a in 1:nspp) { 
			ns = res.ids[a]
			uijmuj[,a] = convolve(muj[,a],kc.spp[,ns])
			uijmuj[,a]= c(uijmuj[ceiling(np/2):(np),a], uijmuj[1:floor(np/2),a] )
			uijmuj[,a] = w.eq[a]*(uijmuj[,a]+abs(min(uijmuj[,a])))

	}
	uijmuj[uijmuj<0]=0

	D.spp = 1+sum(w.eq)

	e1=(ei1-(apply(uijmuj,1,sum)/D.spp) )
	s_zeta=spectrum(e1,spans=3,taper=0, plot=F, na.action=na.omit)
	ls1=length(s_zeta$spec)
	#Rebuild dispersal kernel to accomodate width
	xx_pad=matrix(seq(-ls1,ls1))

	lam2 = matrix(0,np,np)
			#lam2[ceiling(np/2),] = lam2[ceiling(np/2),]+t(as.matrix(lam.tmp))
			diag(lam2) = as.matrix(lam.tmp)
			fd2[,,sa]=mvfft(lam2)
			nr_disp[,sa] = rowSums( Re(mvfft( (fd2[,,sa]*fkd[,,sa]), inverse=T)/(np+1)),na.rm=T	)

	kd.tmp=a_rr/2*exp(-a_rr*abs(xx_pad))
	kd.tmp=kd.tmp/(sum(kd.tmp))
	fdc2=fft(fft(kd.tmp[(ls1+2):(ls1+1+ls1)])*fft(s_zeta$spec[1:ls1],inverse=T)/(1-fft(kd.tmp[(ls1+2):(ls1+1+ls1)])), inverse=T)/(4*ls1)
	#Then multiply by average igr
	agr1= rmFr1/(D.spp*2*3.14)
	return(Re(fdc2[1])*agr1)

}

#=============================================================================
#This function is to get the fitness-density covariance for a single
#1D array of sites, for an arbitrary number of species. This version is
#based on calculating the invader stationary distribution. It requires 
#multi-species equilibriums calculated from e.g. simulations. 
#The analysis is based on Snyder 2008, Snyder and Chesson 2004,2003.
#It takes as variables: 
#	Fr.spp 	The intrinsic range of invader
#	sp.ids	An array indicating which species is the invader (entry 1),
#			and which species are residents (all subsequent entries).
#	nrs.spp The residents stationary distributions
#			-This should be a matrix with species as columns
#	sr.spp 	The survival of the invader.
#	kd.spp 	The invader dispersal kernel
#	alpha.spp Matrix of competition coefficients
#	a_rr 	Parameter for competition kernel
#	kc.spp 	Interspecific competition kernels with residents
#			-If this matrix is a single column then it is used
#			-used iteratively/identically across interactions	 
#=============================================================================

get.fd.cov2 = function (Fr.spp, sp.ids, nrs.spp, sr.spp, alpha.spp, kd.spp, kc.spp){ 
	
	np=length(Fr.spp)
	inv.id = sp.ids[1] #Invader
	res.ids = sp.ids[2:length(sp.ids)] #Residents
	nspp = length (res.ids)
	kd.spp = kd.spp[,inv.id]
	
	#Check for the appropriate number of competition kernels
	if( ncol(kc.spp)<nspp | is.null(ncol(kc.spp)) == TRUE ) { kc.spp=matrix(kc.spp,np,(nspp+1))} 

		#Invader intrinsic fitness
	rmFr1=mean(Fr.spp) 
	ei1=(Fr.spp/rmFr1-1)

	#Resident competitive effects
	nrs.spp=nrs.spp[,res.ids]
	rmcc2=colMeans(nrs.spp,na.rm=T) 
	
	#Weight this by the equilibrium average alpha
	alpha.ir = alpha.spp[inv.id, res.ids]
	w.eq = alpha.ir*rmcc2

	muj=nrs.spp/(matrix(rmcc2, np,nspp,byrow=T))-1
	uijmuj = matrix(0,np,nspp)
	for(a in 1:nspp) { 
			ns = res.ids[a]
			uijmuj[,a] = convolve(muj[,a],kc.spp[,ns])
			uijmuj[,a]= c(uijmuj[ceiling(np/2):(np),a], uijmuj[1:floor(np/2),a] )
			uijmuj[,a] = w.eq[a]*(uijmuj[,a]+abs(min(uijmuj[,a])))

	}
	uijmuj[uijmuj<0]=0

	D.spp = 1+sum(w.eq)

	e1=(ei1-(apply(uijmuj,1,sum)/D.spp) )

	isd.tmp = get.isd(Fr.spp, nrs.spp, inv.id, alpha.spp, sr.spp, kc, kd[,inv.id])

	#Fitness-density covariance, standardized 
	#fdcv2 = rmFr1/(D.spp) * cov(e1,isd.tmp)
	fdcv2 = cov(e1,isd.tmp)

	return(fdcv2)
}

#######
#This block of functions is specifically for the FAST scenario
#It includes get.fast.igr, get.spread.rate, and 

#Get the IGR for a single time step
get.fast.igr = function (Fr.inv, nr.spp, sr,alphas,kd,kc,fkd.yes,n.inv){


	inv.one =pop_lg(Fr.inv, nr.spp, sr,alphas,kd,kc,fkd.yes)[,n.inv]
	gr1.fast=mean(inv.one)/mean(nr.spp[,n.inv])

	return(gr1.fast)

	}

#=======================#
#Pairwise routines
#=======================#

#This function is to get the fitness-density covariance for a single
#1D array of sites, on a pairwise basis. It is based on Snyder 2008, 
#Snyder and Chesson 2004,2003.
#It takes as variables: 
#	Fr.spp 	The intrinsic range of invader
#	nr.spp The resident stationary distribution
#	sr.spp 	The survival of the invader.
#	alpha.ir Interspecific competetion of resident 
#	a_rr 	Parameter for competition kernel
#	kd.spp 	Dispersal kernel
#	kc.spp 	Interspecific competition kernel

get.fd.cov.p = function (Fr.spp, nr.spp, sr.spp, alpha.ir, a_rr, kc.spp){ 
	np=length(Fr.spp)

	rmFr1=mean(Fr.spp) 
	ei1=(Fr.spp/rmFr1-1)

	rmcc2=mean(nr.spp) 
	muj=(nr.spp/rmcc2-1)
	uijmuj = convolve(muj,kc.spp)
	uijmuj= c(uijmuj[ceiling(np/2):(np)], uijmuj[1:floor(np/2)] )

	D.spp = 1+(alpha.ir*rmcc2)

	uijmuj[uijmuj<0]=0
	e1=(ei1-(alpha.ir*rmcc2)/D.spp*uijmuj)
	s_zeta=spectrum(e1,spans=3,taper=0, plot=F, na.action=na.omit)
	ls1=length(s_zeta$spec)
	#Rebuild dispersal kernel to accomodate width
	xx_pad=matrix(seq(-ls1,ls1))
	kd.tmp=a_rr/2*exp(-a_rr*abs(xx_pad))
	kd.tmp=kd.tmp/(sum(kd.tmp))
	fdc2=fft(fft(kd.tmp[(ls1+2):(ls1+1+ls1)])*fft(s_zeta$spec[1:ls1],inverse=T)/(1-fft(kd.tmp[(ls1+2):(ls1+1+ls1)])), inverse=T)/(2*ls1)
	#Then multiply by average igr
	#agr1=(gr1.spp-sr.spp)
	agr1= rmFr1/(D.spp*2*3.14)
	#agr1=1
	return(Re(fdc2[1])*agr1)


}

#Calculate fitness-density covariance for a single 1D array of sites. 
#This approach calculates the invader stationary distribution first using 
#the function get.isd.pair () 
#It takes as variables: 
#	Fr.spp 	The intrinsic range of invader
#	nr.spp The resident stationary distribution
#	sr.spp 	The survival of the invader.
#	alpha.ir Interspecific competetion of resident 
#	kd.spp 	Dispersal kernel
#	kc.spp 	Interspecific competition kernel

get.fd.cov.p2 = function (Fr.spp,nr.spp, sr.spp,alpha.ir,kd.spp,kc.spp ){

	isd.tmp = get.isd.pair (Fr.spp, nr.spp, sr.spp, alpha.ir, kd.spp, kc.spp)
	
	np=length(Fr.spp)

	#Standardized intrinsic fitness
	rmFr1=mean(Fr.spp) 
	ei1=(Fr.spp/rmFr1-1)

	#Standardized resident/competitive effect
	rmcc2=mean(nr.spp) 
	muj=nr.spp/rmcc2-1
	uijmuj = convolve(muj,kc.spp)
	uijmuj= c(uijmuj[ceiling(np/2):(np)], uijmuj[1:floor(np/2)] )

	#This form is for the Leslie-Gower model with D defined as in Appendix B
	D.spp = 1+(alpha.ir*rmcc2)
	e1=(ei1-(alpha.ir*rmcc2)/D.spp*uijmuj)
	
	#Fitness-density covariance, standardized 
	fdcv2 = rmFr1/(D.spp) * cov(e1,isd.tmp)
	#fdcv2 =  cov(e1,isd.tmp)

	return(fdcv2)


}




#=============================================================================
#Environmental distance and site-LGR impact
#=============================================================================

#=============================================================================
#Get the environmental distance for a site between its current and future
#values. This just takes the multi-variate Euclidean distance for whichever
#parameters have been used to define the abiotic conditions. 
#	env1  	current(past) environmental conditions
#	env2  	future(current) environmental conditions
#	env.ind columns with variables to use in the environmental distance
#=============================================================================
get_env_distance = function (env1,env2,env.ind){

	env1 = env1[,env.ind]
	env2 = env2[,env.ind]

	#Normalize all variables relative to the baseline environment
	env1_new = (env1 - matrix(colMeans(env1),dim(env1)[1],dim(env1)[2],byrow=T ))/matrix(sqrt(apply(env1,2,var)),dim(env1)[1],dim(env1)[2],byrow=T )
	env2_new = (env2 - matrix(colMeans(env2),dim(env1)[1],dim(env1)[2],byrow=T ))/matrix(sqrt(apply(env2,2,var)),dim(env2)[1],dim(env2)[2],byrow=T )

	# env1 = (env1)/matrix(apply(env1,2,var),dim(env1)[1],dim(env1)[2],byrow=T )
	# env2 = (env2)/matrix(apply(env2,2,var),dim(env2)[1],dim(env2)[2],byrow=T )


	env_distance = sqrt(rowSums( (env2_new-env1_new)^2) )

	return(env_distance)

}

#=============================================================================
#Get the closest analogue to each site based on environmental distance 
#with the set of all past environments. 
#	env1  	current(past) environmental conditions
#	env2  	future(current) environmental conditions
#	env.ind columns with variables to use in the environmental distance
#=============================================================================
get_analogue = function (env1,env2,env.ind){

	env1 = as.matrix(env1[,env.ind])
	env2 = as.matrix(env2[,env.ind])
	
	#Normalize all variables relative to the baseline environment
	env1_new = (env1 - matrix(colMeans(env1),dim(env1)[1],dim(env1)[2],byrow=T ))/matrix(sqrt(apply(env1,2,var)),dim(env1)[1],dim(env1)[2],byrow=T )
	env2_new = (env2 - matrix(colMeans(env2),dim(env1)[1],dim(env1)[2],byrow=T ))/matrix(sqrt(apply(env2,2,var)),dim(env2)[1],dim(env2)[2],byrow=T )


	#Normalize all variables by variance
	# env1_new = (env1 )/(matrix(apply(env1,2,var),dim(env1)[1],dim(env1)[2],byrow=T ))
	# env2_new = (env2 )/(matrix(apply(env2,2,var),dim(env2)[1],dim(env2)[2],byrow=T ))

	#Distance between a site in variable env2 and all sites in env1
	all_distances = matrix(0,dim(env2)[1],dim(env1)[1] )
	#Find the minimum distance between env1 and env2
	closest_env = matrix(0,dim(env2)[1],1 )
	#The site ID of the closest environmental conditions
	id_closest = matrix(0,dim(env2)[1],1 )	
	#Proportion of the distance that each variable in env.ind contributes to minimum distance
	#all_var_distance = array(matrix(0,dim(env2)[1],dim(env1)[1] ), dim=(dim(env2)[1],dim(env1)[1], length(env.ind)))
	all_var_distance =matrix(0,dim(env2)[1],length(env.ind) )
	

	#Get the distances between future site and all past sites 
	for ( d in 1:dim(env2)[1]){ 

		env2.tmp = matrix(env2_new[d,],dim(env2_new)[1],dim(env2_new)[2],byrow=T)
		var_distances = (env2.tmp-env1_new)^2
		env_distance = sqrt(rowSums( var_distances ) )
		all_distances [d,] = env_distance 
		closest_env[d] = min(env_distance)
		id_closest[d] = which(env_distance == closest_env[d] )
		all_var_distance[d,] = var_distances[id_closest[d],]/(sum(var_distances[id_closest[d],]))

	}
	
	analogue = list(all_distances,closest_env,id_closest,all_var_distance)
	return(analogue)

}

#=============================================================================
# site_impact_sim()
# Use the simulation measurement of LGR to find site impact. Code just
# iteratively removes each site and recalculates LGR. 
#	Uses function pop_lg()
#	ngenst  The total number of environmental-change timesteps to do the calc
# 			for	
#	sp.ids 	The numerical IDs of species (positions in the matrix), with invader
# 			listed first
#=============================================================================
site_impact_sim = function (Frs.spp,nr.spp, sr.spp,alpha.spp, kd.spp, kc.spp, fkd.yes=FALSE,ngenst,sp.ids){
	
	np = dim(as.matrix(Frs.spp))[1] #Space
	nspp = dim(as.matrix(Frs.spp))[2] #Number of species
	site.impacts = matrix(0,np,1)

	for (ss in 1:np){

		Frs.spp.n = Frs.spp[-ss,]
		nr.spp.n = nr.spp[-ss,]

		#If the dispersal kernel is site specific then it will have a dim()[3]
		#Otherwise, it is a normal (not a function of spatial location) kernel
		if (!is.na(dim(kd.spp)[3])){ 
			kd.spp.n = kd.spp[,-ss,]
		}else{
			kd.spp.n = kd.spp[-ss,]
		}
		#Simple numerical LDG: use the invader low-density equilibrium density, 
		#then let it invade against the community: 

		inv.one =pop_lg(Frs.spp.n, nr.spp.n, sr.spp, alpha.spp, kd.spp.n, kc.spp, fkd.yes )[,s]
		site.impacts[ss]=mean(inv.one)/mean(nr.spp.n[,sp.ids[1]])
		if(ss%%ceiling(np/100)==0) {print(ss)}
	}

	return(site.impacts)

}


# # source("./range_coexistence_functionsWALD.R")
# Frs.spp = Frs[t,,]
# nr.spp =nr[t,,]
# sr.spp =sr
# alpha.spp =alphas
# kd.spp = fkd.n
# kc.spp = fkc.n

# kd.spp=kdg
# kc.spp = kc.n
# nr.spp =nr2[t,,]

# Frs.spp=Frs.spp.n
# nr.spp = nr.spp.n
# kd.spp=kd.spp.n

# Frs.spp=Frs[t,,]
# nr.spp=nr[t,,]
# kd.spp=fkd
# kc.spp=fkc