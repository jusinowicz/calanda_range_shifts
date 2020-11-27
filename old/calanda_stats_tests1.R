################################
#
# This code is designed to test different methods for statistically 
# (back)inferring competition coefficients from population equations. 
# The goal is to decide between several basic approaches that could be 
# applied to empirical data. 
#
# The first section of code implements basic integro-difference equations to
# measure species invasion rates. The spatial model is based on Snyder 2003 and 
# Usinowicz and Levine 2018. 
# This version includes the following features: 
# 	1. Leslie Gower population model
#	2. Two spp
#	3. Competition modelled with spatial kernels
#


#==============================================================================
# Load these libraries
#==============================================================================

library(MASS)
library(mgcv)

#==============================================================================
# Function definitions
#==============================================================================
#R equivalent for matlab's meshgrid
meshgrid=function(a,b) {
  list(
       x=outer(b*0,a,FUN="+"),
       y=outer(b,a*0,FUN="+")
       )
} 

#==============================================================================
#Spatial parameters, species parameters, time-related parameters
#==============================================================================
ns=1000 #Lattice width/height -- This should be even
#1D
#xx=seq((-ns/2),(ns/2),by=0.1)
xx=matrix(seq((-ns/2),(ns/2)))
np=length(xx)
ngens=100 #Number of generations 
nspp=2 #number of species
inv.t = 10 #Invasion timestep

#Note: by 3000-5000 generations, populations has reached stable spatial distribution

####For the sampling of invasion rates
#Note, there are two ways to run this: The first is where random fluctuations
#happen within the population dynamics and samples are taken through a loop of 
#invasions. The second is much simpler, which is to run the deterministic model 
#and then add variation afterwards. These two approaches are different in that 
#they generate results where the fluctuations are either correlated or not. 

rnum = 1000
nsamp = rnum
Frr_samp = matrix(0,rnum,2)
inv1r = matrix(0, rnum,1)
inv2r = matrix(0, rnum,1)
Frr = matrix(0,rnum,2)
#==============================================================================


#==============================================================================
#Species parameters
#These will vary based on the model being used
#==============================================================================
Fr=matrix(c(1,1))   #Reproduction rates
### For temporally variable reproduction rates:
cM = matrix(c(0.1, 0, 0 ,0.1),2,2) #Covariance matrix
frn = exp(mvrnorm(ngens, Fr, cM))
#K= matrix(c(1,1))   #Carrying capacity
sr=matrix(c(0.0,0.0)) #Seed survival
gr=matrix(c(1,1))    # Constant germination rates
#Initialize each species fecundity. Could be done as a single array, but would rather keep it this 
#way for now. Then, each species has its own array with 1 less dimension than stc (for 1D). 
#lam1= array(0, dim=c(np,ngens))
#lam2= array(0, dim=c(np,ngens))
#constant lambdas
lam=matrix(c(5,5))

#Competition coefficients 
alphas=matrix( c(1,0.7,0.7,1),2,2)

#==============================================================================
#Parameters for the spatial model
#==============================================================================

# ###### Competition parameters
# #b_rr=1/(100*np) #Essentially global
# b_rr=2

# ###### Dispersal parameters
# #a_rr=1/(100*np) #essentially global
# a_rr=2
#==============================================================================

#==============================================================================
#Internal spatial variables 
#Shouldn't need to tune anything here, but could change competition/dispersal 
#kernel type. 
#==============================================================================
# 1D 
# #Make the full array of space-time coordinates (for 1D space)
# #meshgrid will produce an R data frame
# stc.temp=meshgrid(1:ngens,xx)
# #Convert the coordinates into an array for easier access
# stc=array(c(matrix(stc.temp$x,ngens,np,byrow=T),matrix(stc.temp$y,ngens,np,byrow=T)),dim=c(ngens,np,2)) 
# #Note: in 1d, things are a bit odd. stc.temp$x ends up as time, and $y as space. So in stc, stc[,,1] is time
# and stc[,,2] is space. stc[,,1] is the same as the column vector 1:ngens repeated over np columns. Then
# stc[,,2] is space as a row vector (xx) repeated over ngens rows. 

# #Competition kernel 
# kc= b_rr/2*exp(-b_rr*abs(xx))

# #Uniform spatial competition
# #kc=matrix(1,dim(xx)[1], dim(xx)[2]) 

# #Moore neighborhood. In this case, b_rr is the radius (b_rr in each direction from individual) 
# #kc=matrix(0, 1, np)
# #for( i in 0:(2*b_rr)){ kc[((ceiling(np/2)-b_rr)+i)]=1} #Find the middle site, go b_rr sites left, then start
# #making 1s until reaching a distance of b_rr to the right of middle site
# #kc[ceiling(np/2)]=0

#Normalize the kernel:
# kc=kc/(sum(kc)) #Normalize so that integral =1 
# fkc=fft(kc)#/(np+1) #Fourier transform
# #


# #Dispersal kernel 
# kd= a_rr/2*exp(-a_rr*abs(xx))
# #kd=matrix(1,dim(xx)[1], dim(xx)[2]) #Uniform dispersal
# kd=kd/(sum(kd))
# fkd=fft(kd)#/(np+1)
#==============================================================================


#==============================================================================
#Outermost loop: 
#==============================================================================

for( s in 1:rnum){ 

#The population matrix. Populations are of seeds. Columns represent space, rows are time. 
nr1=matrix(0,ngens,np)
nr2=matrix(0,ngens,np)
#nr[1,]= 10/2*exp(-abs(xx))#Initial Conditions
#nr[1,][nr[1,]<0.1]=0
nr1[1,np/2]= 0.01
nr2[1,np/2]= 0.01
#nr2[1,]=nr2[1,]+0.96/np

#nonspatial model
nrns1=matrix(0,ngens,1)
nrns1[1]=sum(nr1[1,])
nrns2=matrix(0,ngens,1)
nrns2[1]=sum(nr2[1,])


#==============================================================================
#Population dynamics (internal loop)
#==============================================================================

	for (n in 1:(ngens-1)) {

	#Seedlings germinate and produce seeds at rate Fr, weighted by competition
	if(n < inv.t) { 
		nr1=matrix(0,ngens,np)
		nr1[n,np/2]= 0.01
		#print(nr1[1,np/2])
		nrns1[n] =0.01
		}

	#==============================================================================
	#This section is for the spatial model. 
	#==============================================================================

	# #Germination fraction
	# #1.constant germination rates
	# fp1=fft(nr1[n,]) #/(np+1) #fft of population
	# fp2=fft(nr2[n,]) #/(np+1) #fft of population

	# #There are a total of 4 possible combinations for competition. All are the same now.
	# #That changes if the kernels differ for e.g. 1 from 2, and 2 from 1 (e.g. alpha_12 vs. alpha_21)
	# #As far as I can tell, literature always assumes 12=21 (e.g. all of Snyder papers) 
	# #Inverse fft of convolution of population and competition kernel -- gives competition 
	# #1 from itself
	# Cr11=Re(fft((1+alphas[1,1]*fp1*fkc),inverse=T)/(np+1)) 
	# #1 from 2
	# Cr12=Re(fft((1+alphas[1,2]*fp2*fkc),inverse=T)/(np+1)) 
	# #2 from itself
	# Cr22=Re(fft((1+alphas[2,2]*fp2*fkc),inverse=T)/(np+1)) 
	# #2 from 1
	# Cr21=Re(fft((1+alphas[2,1]*fp1*fkc),inverse=T)/(np+1)) 

	# #competition-weighted seed production
	# #1. Constant Fr
	# # lam1r1=t(Fr[1])/(c(Cr11[ceiling(np/2):(np)],Cr11[1:floor(np/2)])+c(Cr12[ceiling(np/2):(np)],Cr12[1:floor(np/2)])) 
	# # lam1r2=t(Fr[2])/(c(Cr22[ceiling(np/2):(np)],Cr22[1:floor(np/2)])+c(Cr21[ceiling(np/2):(np)],Cr21[1:floor(np/2)])) 

	# #2.fluctuating Fr
	# lam1r1=t(frn[n, 1])/(c(Cr11[ceiling(np/2):(np)],Cr11[1:floor(np/2)])+c(Cr12[ceiling(np/2):(np)],Cr12[1:floor(np/2)])) 
	# lam1r2=t(frn[n, 2])/(c(Cr22[ceiling(np/2):(np)],Cr22[1:floor(np/2)])+c(Cr21[ceiling(np/2):(np)],Cr21[1:floor(np/2)])) 


	# #Note: the form of Cr is necessary because the inverse fft puts things in the (standard) order of frequencies 1: np/2+1, -np/2: -1
	# #Piecing it together as above (and below after the dispersal step) is necessary. This is equivalent to using fftshift in matlab. 

	# lam1r1=lam1r1*(nr1[n,]>1e-4)
	# lam1r2=lam1r2*(nr2[n,]>1e-4)

	# #If both species have zero germination, then NAs are produced: 
	# lam1r1[is.na(lam1r1)]=0
	# lam1r2[is.na(lam1r2)]=0

	# #Seeds disperse
	# fd1=fft(lam1r1*nr1[n,])
	# fd2=fft(lam1r2*nr2[n,])

	# nr_disp1=Re(fft((fd1*fkd),inverse=T)/(np+1))
	# nr_disp1=c(nr_disp1[ceiling(np/2):(np)],nr_disp1[1:floor(np/2)])

	# nr_disp2=Re(fft((fd2*fkd),inverse=T)/(np+1))
	# nr_disp2=c(nr_disp2[ceiling(np/2):(np)],nr_disp2[1:floor(np/2)])

	# #Final population step
	# #1. constant Fr
	# nr1[n+1, ] = nr_disp1+nr1[n,]*sr[1]
	# nr2[n+1, ] = nr_disp2+nr2[n,]*sr[2]
	#==============================================================================

	#==============================================================================
	#Nonspatial model
	#1. constant Fr
	# nrns1[n+1,] = nrns1[n]*Fr[1]/(1+nrns1[n]*alphas[1,1]+nrns2[n]*alphas[1,2])+nrns1[n]*sr[1]
	# nrns2[n+1,] = nrns2[n]*Fr[2]/(1+nrns1[n]*alphas[2,1]+nrns2[n]*alphas[2,2])+nrns2[n]*sr[2]

	#2. fluctuating Fr
	nrns1[n+1,] = nrns1[n]*frn[n, 1]/(1+nrns1[n]*alphas[1,1]+nrns2[n]*alphas[1,2])+nrns1[n]*sr[1]
	nrns2[n+1,] = nrns2[n]*frn[n, 2]/(1+nrns1[n]*alphas[2,1]+nrns2[n]*alphas[2,2])+nrns2[n]*sr[2]
	#==============================================================================


	}

#Outer loop samples over multiple invasions
inv1r[s] = nrns1[inv.t]/nrns1[(inv.t-1)]
inv2r[s] = nrns2[inv.t]/nrns2[(inv.t-1)]
Frr[s,1] = frn[n, 1]
Frr[s,2] = frn[n, 2]

}

plot(rowMeans(nr1), t="l")
lines(rowMeans(nr2), col="red")
lines(nrns1,col="blue")

inv1=nrns1[inv.t]/nrns1[(inv.t-1)]
inv2=nrns2[inv.t]/nrns2[(inv.t-1)]

a22s = Fr[2]/(inv2-sr[2])-1
a12s = Fr[1]/(inv1-sr[1])-1
a12 =a12s/a22s

#==============================================================================
#Inferring the competiteion coefficients from invasion stats:
#==============================================================================
######With variance added afterwards: 
#Notes:
#	When the variance in Fr and the igr are independent, then 
#	bootstrapping and NLS works the best with small sample sizes.
#	The bootstrapping of the ratio is totally wrong.
#	Averaging and then taking the ratio (1) works fine, but 
#	becomes a worse and worse estimate with small sample sizes.
#
#	When the variance in Fr and igr are not independent -- because 
# 	variance is added to the demographic parameteres, making var(igr)
# 	a f(var(Fr)) -- then all 3 approaches work ok, but NLS still works
#	the best. 
#
#	Another cautionary tale from playing with this simulation: Be
#	careful about working with the log(data)!! Log-transforming the
#	the underlying Frs does not make sense, because the population-level
#	result (the inv1 and invr1s here) are a function of a lognormal
#	distrubtion. I.e., on exp( rand_norm(m,v)), and not on the rand_norm(m,v)
#	itself. 
# 
#==============================================================================

# Uncomment this if running scenario 2: Deterministic model, add 
# the variation afterwards. 
# inv1r = exp(inv1+rnorm(rnum,sd=1))
# inv2r = exp(inv2+rnorm(rnum,sd=1))
# Frr = matrix(0,rnum,2)
# Frr[,1] = exp(Fr[1]+rnorm(rnum,sd=1))
# Frr[,2] = exp(Fr[2]+rnorm(rnum,sd=1))

inv1r_samp = inv1r[ceiling(runif(nsamp)*rnum) ]
inv2r_samp = inv2r[ceiling(runif(nsamp)*rnum) ]
Frr_samp = matrix(0,nsamp,2)
Frr_samp[,1] = Frr[,1][ceiling(runif(nsamp)*rnum) ] 
Frr_samp[,2] = Frr[,2][ceiling(runif(nsamp)*rnum) ] 

#1: With means of variable
#inv1p=mean(log(inv1r_samp))
#inv2p=mean(log(inv2r_samp))

inv1p=mean((inv1r_samp))
inv2p=mean((inv2r_samp))

Fr1 = mean((Frr_samp[,1]))
Fr2 = mean((Frr_samp[,2]))

a22p = Fr2/(inv2p-sr[2])-1
a12p = Fr1/(inv1p-sr[1])-1
a12pp =a12p/a22p

#2: Boostrap the solution itself
bsamp=1000
# crdat1 = log(data.frame(cbind(inv1r_samp[ceiling(runif(bsamp)*nsamp)], 
# 	 		Frr_samp[,1] [ceiling(runif(bsamp)*nsamp)])))
# crdat2 = log(data.frame(cbind(inv2r_samp[ceiling(runif(bsamp)*nsamp)], 
# 	 		Frr_samp[,2] [ceiling(runif(bsamp)*nsamp)])))
crdat1 = (data.frame(cbind(inv1r_samp[ceiling(runif(bsamp)*nsamp)], 
	 		Frr_samp[,1] [ceiling(runif(bsamp)*nsamp)])))
crdat2 = (data.frame(cbind(inv2r_samp[ceiling(runif(bsamp)*nsamp)], 
	 		Frr_samp[,2] [ceiling(runif(bsamp)*nsamp)])))
a22b = mean(crdat2[,2]/(crdat2[,1]-sr[2]))-1
a12b = mean(crdat1[,2]/(crdat1[,1]-sr[1]))-1
a12bb = a12b/a22b

#Bootstrap but then fit using NLS: 
bsamp=1000
# crdat1 = log(data.frame(cbind(inv1r_samp[ceiling(runif(bsamp)*nsamp)], 
# 	 		Frr_samp[,1] [ceiling(runif(bsamp)*nsamp)])))
# colnames(crdat1) = c("l","R")
# crdat2 = log(data.frame(cbind(inv2r_samp[ceiling(runif(bsamp)*nsamp)], 
# 	 		Frr_samp[,2] [ceiling(runif(bsamp)*nsamp)])))
# colnames(crdat2) = c("l","R")

rdat1 = (data.frame(cbind(inv1r_samp[ceiling(runif(bsamp)*nsamp)], 
	 		Frr_samp[,1] [ceiling(runif(bsamp)*nsamp)])))
colnames(crdat1) = c("l","R")
crdat2 = (data.frame(cbind(inv2r_samp[ceiling(runif(bsamp)*nsamp)], 
	 		Frr_samp[,2] [ceiling(runif(bsamp)*nsamp)])))
colnames(crdat2) = c("l","R")

a22nls = nls(l ~ R/(1+cr),data=crdat2, start=list(cr=0.5))
	print(summary(a22nls))
a22n = a22nls$m$getPars()
	
a12nls = nls(l ~ R/(1+cr),data=crdat1, start=list(cr=0.5))
	print(summary(a12nls))
a12n = a12nls$m$getPars()	
a12nn = a12n/a22n

#==============================================================================
# GAMs and species distributions. 
# Some simulations and fitting procedures to help understand the best approach 
# to fitting GAMs to species ranges. 
#==============================================================================
#==============================================================================

## Fake some data... 
 
set.seed(0) 
n = 100 
mu1 = 1; sig1 = 0.5; 
x = runif(n)*10-2;x = sort(x); 
f = exp(-((x - mu1)^2/(2*sig1)));y = f+rnorm(100)*0.1;plot(x,y) 
dat = data.frame(x=x,y=y) 

## Create a spline basis and penalty, making sure there is a knot 
## at the constraint point, (0 here, but could be anywhere) 
knots = data.frame(x=seq(-2,8,length=11)) ## create knots 

## set up smoother... 
sm = smoothCon(s(x,k=length(knots$x),bs="cr"),dat,knots=knots)[[1]] 

## 3rd parameter is value of spline at knot location 0, 
## set it to 0 by dropping... 
X = sm$X[,-3]        ## spline basis 
S = sm$S[[1]][-3,-3] ## spline penalty 
off = y*0 + .6       ## offset term to force curve through (0, .6) 

off = y*0+0
## fit spline constrained through (0, .6)... 
b = gam(y ~ X - 1 + offset(off),paraPen=list(X=list(S))) 
lines(x,predict(b)) 

points(0,0.6,col="green")
## compare to unconstrained fit... 
b.u = gam(y ~ s(x,k=9),data=dat,knots=knots) 
lines(x,predict(b.u),col="red")

#This is most similar to what I'm dealing with with species ranges 
#where measurements have been made at spatial locations.
b.ud = gam(y ~ s(x,k=4),data=dat,knots=knots) 
lines(x,predict(b.ud),col="blue")


#########
#Multiple point constraints
# e.g. (-1,0) and (5,0)


### Use the formula for a Normal distribution to create a uni-humped 
### shaped distribution of points. 
set.seed(0) 
n = 100 
mu1 = 1; sig1 = 0.5; 
x = runif(n)*10-2;x = sort(x); 
f = exp(-((x - mu1)^2/(2*sig1)));y = f+rnorm(100)*0.1;plot(x,y) 
dat = data.frame(x=x,y=y) 

## Create a spline basis and penalty, making sure there is a knot 
## at the constraint point, (0 here, but could be anywhere) 
knots = data.frame(x=seq(-2,8,length=11)) ## create knots 
cp_x = c(2,8)
#cp_x = c(2)
cp_y = c(0.6,0.2)

## set up smoother... 
sm = smoothCon(s(x,k=length(knots$x),bs="cr"),dat,knots=knots)[[1]] 

## 3rd parameter is value of spline at knot location 0, 
## set it to 0 by dropping... 
X = sm$X[,-cp_x]        ## spline basis 
S = sm$S[[1]][-(cp_x),-(cp_x)] ## spline penalty 
off1 = y*0 + cp_y[1] ## offset term to force curve through (0, .6) 
off2 = y*0 + cp_y[2]

#This is to divide the domain into the two different offsets, based
#on where the peak of the Gaussian is in this problem.
peak1 = 1 
off=off1
off[x>=peak1] = off2[x>=peak1]
off[x<peak1] = off1[x<peak1]

## fit spline constrained through (0, .6)... 
b = gam(y ~ X - 1 +offset(off),paraPen=list(X=list(S))) 
lines(x,predict(b)) 
points(cp_x[1]-2,cp_y[1],col="green")
points(cp_x[2]-2,cp_y[2],col="green")

## compare to unconstrained fit... 
b.u = gam(y ~ s(x,k=11),data=dat,knots=knots) 
lines(x,predict(b.u),col="red")

### Now use Normal distribution and instead of assuming arbitrary points, 
### pick points that make sense to constrain the tails of the fit with 
### only 4 knots. 

### Use the formula for a Normal distribution to create a uni-humped 
### shaped distribution of points. 
set.seed(0) 
cp_x = c(1,7)
#cp_x = c(2)
cp_y = c(0,0)

n = 100 
mu1 = 1; sig1 = 0.5; 
x = runif(n)*10-2;x = sort(x); 
f = exp(-((x - mu1)^2/(2*sig1)));y = f+rnorm(100)*0.1;plot(x,y) 
y[ x< cp_x[1]] =NA
y[ x> cp_x[2]] =NA
dat = data.frame(x=x,y=y) 

## Create a spline basis and penalty, making sure there is a knot 
## at the constraint point, (0 here, but could be anywhere) 
#knots = data.frame(x=seq(-2,8,length=11)) ## create knots 
knots = data.frame(x=seq(-2,4,length=7)) ## create knots 

## set up smoother... 
sm = smoothCon(s(x,k=length(knots$x),bs="cr"),dat,knots=knots)[[1]] 

## 3rd parameter is value of spline at knot location 0, 
## set it to 0 by dropping... 
X = sm$X[,-cp_x]        ## spline basis 
S = sm$S[[1]][-(cp_x),-(cp_x)] ## spline penalty 
off1 = y*0 + cp_y[1] ## offset term to force curve through (0, .6) 
off2 = y*0 + cp_y[2]

#This is to divide the domain into the two different offsets, based
#on where the peak of the Gaussian is in this problem.
peak1 = 1 
off=off1
off[x>=peak1] = off2[x>=peak1]
off[x<peak1] = off1[x<peak1]

## fit spline constrained through (0, .6)... 
b = gam(y ~ X - 1 +offset(off),paraPen=list(X=list(S))) 
lines(x,predict(b)) 
points(cp_x[1]-2,cp_y[1],col="green")
points(cp_x[2]-2,cp_y[2],col="green")

## compare to unconstrained fit... 
b.u = gam(y ~ s(x,k=5),data=dat,knots=knots) 
lines(x,predict(b.u),col="red")


#==============================================================================
#
# A handy example from Simon Wood
# http://r.789695.n4.nabble.com/Use-pcls-in-quot-mgcv-quot-package-to-achieve-constrained-cubic-spline-td4660966.html
#==============================================================================

## Example constraining spline to pass through a 
## particular point (0,.6)... 

## Fake some data... 
 
set.seed(0) 
n <- 100 
x <- runif(n)*4-1;x <- sort(x); 
f <- exp(4*x)/(1+exp(4*x));y <- f+rnorm(100)*0.1;plot(x,y) 
dat <- data.frame(x=x,y=y) 

## Create a spline basis and penalty, making sure there is a knot 
## at the constraint point, (0 here, but could be anywhere) 
knots <- data.frame(x=seq(-1,3,length=9)) ## create knots 
## set up smoother... 
sm <- smoothCon(s(x,k=9,bs="cr"),dat,knots=knots)[[1]] 

## 3rd parameter is value of spline at knot location 0, 
## set it to 0 by dropping... 
X <- sm$X[,-3]        ## spline basis 
S <- sm$S[[1]][-3,-3] ## spline penalty 
off <- y*0 + .6       ## offset term to force curve through (0, .6) 

## fit spline constrained through (0, .6)... 
b <- gam(y ~ X - 1 + offset(off),paraPen=list(X=list(S))) 
lines(x,predict(b)) 

## compare to unconstrained fit... 

b.u <- gam(y ~ s(x,k=9),data=dat,knots=knots) 
lines(x,predict(b.u))


#==============================================================================
#Miscellaneous stats notes 
#==============================================================================
### A note about fitting GAMMs with mgcv vs. gamm4: 
### The following are equivalent: 
#mgcv and bs = "re"
rr_gam= gam(nr.inflor ~ s(gs_mean_temp, bs="cr",k=kn[sp])+s(year, bs="re",by = flyes),data=rr_dat)
rr_tmp=as.vector(predict.gam(rr_gam, newdata=data.frame(xx_new,flyes=0),type= "response"))

#mgcv and bs = "re", with "exclude" used in predict.gam
rr_gam= gam(nr.inflor ~ s(gs_mean_temp, bs="cr",k=kn[sp])+s(year, bs="re",by = flyes),data=rr_dat)
rr_tmp=as.vector(predict.gam(rr_gam, newdata=xx_new,type= "response",exclude = "s(year)"))

#gamm4 
rr_gam1= gamm4(nr.inflor ~ s(gs_mean_temp, k=kn[sp]),random = ~(1|year),data=rr_dat)
rr_tmp1=as.vector(predict.gam(rr_gam1$gam, newdata=xx_new,type= "response"))


#Ignore this: 
rr_gamA = gamm4(nr.inflor ~ 1,random = ~(1|year),data=rr_dat)
new.inflor = residuals(rr_gamA$mer)+rr_gamA$gam$coefficients
