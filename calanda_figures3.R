###########################################################
#Figure making
###########################################################
#
#==============================================================================
#	Figure 1: Just make some ranges for N species and their 
#	expansions for a picture. Should encompass 3 main scenarios: 
#	A. Range translocation, 
#	B. Range expansion
#	C. Translocation and expansion
#==============================================================================
library(abind) # To combine arrays along desired dimension
#==============================================================================
# Load Data
#==============================================================================

#This requires user input!

# variable.list=list("l1", "D", "var_mu_Us","cov_e_mu_Us",
# "cov_lam_vc", "cov_lam_vc2", "Elam1", "Elam2", "gr1.n", "gr1", "y.full", "w.eq.full", "l1.p",
# "y.p", "w.eq.p", "D.p", "var_mu_Us.p", "cov_e_mu_Us.p","cov_lam_vc.p", "cov_lam_vc.p2",
# "Elam1.p","Elam2.p", "gr1.n.p", "gr1.p","ldg.sim","ldgIMP.sim","covIMP.sim","ldg.sim.p","ldgIMP.sim.p",
# "covIMP.sim.p","Frs")

variable.list=list("l1", "D", "var_mu_Us","cov_e_mu_Us",
"cov_lam_vc", "cov_lam_vc2", "Elam1", "Elam2", "gr1.n", "gr1", "y.full", "w.eq.full", "l1.p",
"y.p", "w.eq.p", "D.p", "var_mu_Us.p", "cov_e_mu_Us.p","cov_lam_vc.p", "cov_lam_vc.p2",
"Elam1.p","Elam2.p", "gr1.n.p", "gr1.p","ldg.sim","ldgIMP.sim","covIMP.sim","ldg.sim.p","ldgIMP.sim.p",
"lded","lded.p","covIMP.sim.p","Frs")

#List of scenario file names 
#file.name.list= list("calanda_igr4_sdef2_2017_waldmean.var","calanda_igr4_sdef2_fut_2017_waldmean.var")
#file.name.list= list("calanda_igr6_multi_2017_waldmean.var","calanda_igr6_multi_fut_2017_waldmean.var")
#file.name.list= list("calanda_igr6_multiB_2017_waldmean.var","calanda_igr6_multi_futB_2017_waldmean.var")
#file.name.list= list("calanda_igr7_multiB_2017_waldmean.var","calanda_igr7_multiB_fut_2017_waldmean.var")
#file.name.list= list("calanda_igr8_multiB_2017_waldmean.var","calanda_igr8_multiB_fut_2017_waldmean.var")
#file.name.list= list("calanda_igr9_multiB_2017_waldmean.var","calanda_igr9_multiB_fut_2017_waldmean.var")
#file.name.list= list("calanda_igr9_multiC_2017_waldmean.var","calanda_igr9_multiC_fut_2017_waldmean.var")
#file.name.list= list("calanda_igr10_multi_2017_waldmean.var","calanda_igr10_multi_fut_2017_waldmean.var")
#file.name.list= list("calanda_igr11_multi_2017_waldmean.var","calanda_igr11_multi_fut_2017_waldmean.var")
file.name.list= list("calanda_igr12_multi_2017_waldmean.var","calanda_igr12_multi_fut_2017_waldmean.var")

#file.name.list= list("calanda_igrsites_multi_2017_waldmean.var","calanda_igrsites_multi_fut_2017_waldmean.var")
#file.name.list= list("calanda_igrsites2_multi_2017_waldmean.var","calanda_igrsites2_multi_fut_2017_waldmean.var")
#file.name.list= list("calanda_igrsites3_multi_2017_waldmean.var","calanda_igrsites3_multi_fut_2017_waldmean.var")
#In order to deal with each of the scenarios simultaneously, it is necessary
#to systematically rename variables based on their scenario. 
#Assosciated prefixes for renaming
prefix.list= list("start","shift")

#Rename the files from each scenario
var.length=length(variable.list)
nscen = length( prefix.list)


for (g in 1:nscen){
	load(file.name.list[[g]])
	for( f in 1:(var.length)) {      
		new.name=paste(prefix.list[[g]],variable.list[[f]],sep="")	
		assign(new.name,eval(as.name(variable.list[[f]])))

	}	
}

#When the files/scenarios actually represent shifts from i.c.s, 
#then the files can be combined into a single scenario. 


for( f in 1:(var.length-1)) { 
	tmp.var=NULL
	for (g in 1:(nscen)){   
		tmp.name=paste(prefix.list[[g]],variable.list[[f]],sep="")	
		tmp.var=rbind(tmp.var,eval(as.name(tmp.name)))	
	}	

	tmp.name2=paste("comb",variable.list[[f]],sep="")	
	assign(tmp.name2,tmp.var)

}

#Combine the intrinsic ranges
combFrs = abind(startFrs, shiftFrs, along=1)

prefix.list= list("comb")
nscen=1


#==============================================================================
# Figure: Intrinsic ranges, before and after climate change
#==============================================================================
source("./range_coexistence_functionsD.R")

nspp=3 #Number of species
ns=dim(Frs)[2] #Lattice width/height -- This should be even
ngens=1 #Number of generations for environmental change
iconfig=1 #NUmber of interattions of initial configuration of species

#Make these automatically:
xx=matrix(seq((-ns/2),(ns/2)))
np=length(xx)
ngenst=iconfig+ngens

#color.use = c("gray30","red", "darkgreen", "gray30","blue", "gray50","gray30","gray60" )
color.use = c( "darkgreen","red","blue" )


#fig.name = paste("calanda_rangeshifts_3sppC.pdf",sep="")
#pdf(file=fig.name, family='Helvetica', pointsize=16)

par(mfrow=c(1,1))
plot(combFrs[1,,1],t="l",xlab="Location", ylab="Fitness", cex.lab=1.7,cex.axis=1.2,ylim=c(0,30),
		col=color.use[1], xaxs="i",xaxt='n',yaxs="i",yaxt='n',bty="n",lwd = 3,lty=2)
abline(h=0)
for(j in 2:nspp) {
	lines( combFrs[1,,j],col = color.use[j],lwd=3,lty=2)
}

#plot(Frs[ngenst,,1],t="l",xlab="Location", ylab="Fitness", cex.lab=1.7,cex.axis=1.2,ylim=c(0,3),
#		col=color.use[1],xaxs="i",xaxt='n',yaxs="i",yaxt='n',bty="n",lwd=3)
abline(h=0)
for(j in 1:nspp) {
	lines( combFrs[ngenst,,j],col = color.use[j],lwd=3)
}

dev.off()


#Set this up 

#==============================================================================
# Figure: Style 2: Intrinsic ranges, before and after climate change
#==============================================================================

#Spatial scale: Assume that the gradient is from 500m to 2500m in units of 10m
elevations=c(1000,1400,1600,1800,2000)
el1=0
el2=3500

# el1=100
# el2=3000

# el1=700
# el2=3000

#Make some internal spatial variables:
xx1= matrix(seq(el1,el2,1))
xx2= matrix(xx1^2)
yx= matrix(2017,length(xx1),1)


xx=cbind(xx1,yx,xx2)
colnames(xx)=c('elevation','year','e2')
xx=data.frame(xx)
xx$elevation = as.matrix(xx$elevation)
xx$year = as.matrix(xx$year)
xx$e2 = as.matrix(xx$e2)

#=============================================================================

epoints = match( elevations,xx$elevation)

fig.name = paste("seeds_kriged_igrsites3.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

#par(mfrow=c(1,1),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
ylim=c(-1,70)
plot(xx$elevation, combFrs[1,,1], t="l", ylab="Seeds (per-capita, Kriged) ", xlab="Elevation", xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim)
lines(xx$elevation, combFrs[1,,1])
lines(xx$elevation, combFrs[1,,2],col="red")
lines(xx$elevation, combFrs[1,,3],col="blue")
lines(xx$elevation, combFrs[2,,1],lty=2)
lines(xx$elevation, combFrs[2,,2],col="red",lty=2)
lines(xx$elevation, combFrs[2,,3],col="blue",lty=2)

abline(v=c(elevations),lwd=10,col="grey80")
points(xx$elevation[epoints], combFrs[1,epoints,1])
points(xx$elevation[epoints], combFrs[1,epoints,2],col="red")
points(xx$elevation[epoints], combFrs[1,epoints,3],col="blue")
points(xx$elevation[epoints], combFrs[2,epoints,1])
points(xx$elevation[epoints], combFrs[2,epoints,2],col="red")
points(xx$elevation[epoints], combFrs[2,epoints,3],col="blue")

# f1 = flower_prob_act*flower_act
# points(xx$elevation[epoints], f1[,1],pch=3)
# points(xx$elevation[epoints], f1[,2],col="red",pch=3)
# points(xx$elevation[epoints], f1[,3],col="blue",pch=3)

dev.off()


################################################################
#Figure 3A: Non-Arch (linear) figures showing how the igr of each 
# scenario moves according to the benefit of either range width 
# or range overlap.
# This version for the 2SPP picture. 
# Components STANDARDIZED BY MEAN coefficients
# Figures are a SINGLE PLOT
################################################################
#Set these variables manually (hopefully you kept good notes)
iconfig=1
ngens=0
ngenst=iconfig+ngens
nspp=3
scale=10
#tol2=c(100,100,100,100)
tol2=c(5,5,5,5)

#Skipping columns
uu=c(1,2,3,5)
#Different limits

xlims = matrix(c(-0.0, 2.3, -0.0,2.3,-0.0,2.3, -0.0,2.3),4,2,byrow=T)
x.axx = matrix(c(0.2,2.3,0.4,0.2,2.3,0.4,0.2,2.3,0.4, 0.2,2.3,0.4 ),4,3,byrow=T)

#fig.name = paste("all4pairs_linear_widthOverlap_2Dgr4.pdf",sep="")
#pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

#Set this up as a panel of 4, or see below for individual files. 
par(mfrow=c(1,4),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
#This is a constant to scale the range overlap
min_overlap1=0
#Find the minimum value across the 3 normal scenarios
for( g in 1:3){ 
	for(u in 1:6){

		#Identify the variable names by scenario, then assign them to the right
		#variable for plotting
		cov_e_mu_Us.name=paste(prefix.list[[g]],variable.list[[18]] ,sep="")
		cov_lam_vc.name=paste(prefix.list[[g]],variable.list[[30]] ,sep="")
		l1.name=paste(prefix.list[[g]],variable.list[[13]] ,sep="")
		D.name = paste(prefix.list[[g]],variable.list[[16]] ,sep="")

		assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name)) )
		assign("cov_lam_vc",eval(as.name(cov_lam_vc.name)) )
		assign("l1",eval(as.name(l1.name)) )
		assign("D",eval(as.name(D.name)) )

		#Create a mean standardization coefficient to simplify presentation. 
		C1 =  median((l1[1,u]/(D[1,u])^3), na.rm=T)
		C2 =  median(l1[1,u]/(D[1,u])^2, na.rm=T )
		
		min_overlap1 = min(c (min_overlap1,(-C2*cov_e_mu_Us[tol2[g],u]+cov_lam_vc[tol2[g],u]),(-C2*cov_e_mu_Us[1,u]+cov_lam_vc[1,u])),na.rm=T)
	}
}

#Offset points in the y direction for readibility 	
j1 = jitter((1:4)*0, 0.5)

#Now make the plots: 
for( g in 1:4){ 

	color.use=list("darkgreen","darkgreen","red", "red","blue","blue")
	up.low=list("upper","lower")

	#Or just switch between devices for plotting
	#dev.set(g)
	#par(mfrow=c(1,4), mar=c(5,6,6,2))

	
	ylim=c(0,0.2)
	#xlim=c(0,0.8)
	xlim=xlims[g,]
	#plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Range width", xlab="Range overlap",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3 )
	if (g ==1){
	plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Benefit of Competitor Aggregation", xlab="",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ann=F,yaxt='n',xaxt='n' )
	axis(side=1, labels=T, at=seq(x.axx[g,1], x.axx[g,2],x.axx[g,3]))
	axis(side=2, labels=T, at=seq(0.1,0.6,0.2))

	} else { 
	plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="", xlab="",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ann=F,yaxt='n',xaxt='n')
	axis(side=1, labels=T, at=seq(x.axx[g,1], x.axx[g,2],x.axx[g,3]))

		} 

	for (ug in 1:4){ 
			u=uu[ug]
			
			#Identify the variable names by scenario, then assign them to the right
			#variable for plotting
			gr.name=paste(prefix.list[[g]],variable.list[[24]] ,sep="")
			var_mu_Us.name=paste(prefix.list[[g]],variable.list[[17]] ,sep="")
			cov_e_mu_Us.name=paste(prefix.list[[g]],variable.list[[18]] ,sep="")
			cov_lam_vc.name=paste(prefix.list[[g]],variable.list[[30]] ,sep="")
			#For scaling
			l1.name=paste(prefix.list[[g]],variable.list[[13]] ,sep="")
			w.eq.name = paste(prefix.list[[g]],variable.list[[15]] ,sep="")
			D.name = paste(prefix.list[[g]],variable.list[[16]] ,sep="")

			assign("gr",eval(as.name(gr.name)) )
			assign("var_mu_Us",eval(as.name(var_mu_Us.name)) )
			assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name)) )
			assign("cov_lam_vc",eval(as.name(cov_lam_vc.name)) )
			assign("l1",eval(as.name(l1.name)) )
			assign("D",eval(as.name(D.name)) )
			assign("cr_mean",D-1 )

			#This line represents persistence in the absence of spatial mechanisms
			#mlin.non1 =1-(l1[iconfig,u]/D[iconfig,u]+0.9)-min_overlap1
			mlin.non1 =0.376
			lin.non1= -(seq(xlim[1],xlim[2],by=0.1))+mlin.non1
			#lines((seq(xlim[1],xlim[2],by=0.1)),lin.non1,col=color.use[[u]],lty = 2,lwd=1.5)
			lines((seq(xlim[1],xlim[2],by=0.1)),lin.non1,lty = 2,lwd=1.5)

			#Create a mean standardization coefficient to simplify presentation. 
			#C1 =  median((l1[,u]*cr_mean[,u]^2/(D[,u])^3), na.rm=T)
			#C2 =  median(l1[,u]*cr_mean[,u]/(D[,u])^2, na.rm=T )

			C1 = (l1[1,u]/(D[1,u])^3)
			C2 = l1[1,u]/(D[1,u])^2


			#Spatial mechanisms add this much
			width_y1 = C1*var_mu_Us[iconfig,u] 
			#Strongest negative impact (this is only accurate if data include case of complete overlap)
			#min_overlap1 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
			overlap_x1 = -min_overlap1+((-C2*cov_e_mu_Us[iconfig,u])+cov_lam_vc[iconfig,u])+0.08
			#Point where lines intercept: x = (b1-b2)/2
			#non.x1 = (mlin.non1 - (width_y1 - overlap_x1) )/2
			#non.y1 = -non.x1+ mlin.non1 
			#Intersection through origin instead:
			non.x1 = -mlin.non1/(-1-(width_y1/overlap_x1) )
			non.y1 = -non.x1+ mlin.non1 
			#points(non.x1,non.y1,col=color.use[[u]], pch=16)
			points(overlap_x1,width_y1+j1[ug],col=color.use[[u]] )
			#segments(non.x1,non.y1,overlap_x1,width_y1,col=color.use[[u]],lty=2)			
			print(C2)

			#Now the later point
			#mlin.non2 =1-(l1[(tol2[g]),u]/D[(tol2[g]),u]+0.9)-min_overlap1
			mlin.non2 =0.376
			lin.non2= -(seq(xlim[1],xlim[2],by=0.1))+mlin.non2
			#lines((seq(xlim[1],xlim[2],by=0.1)),y=lin.non2,col=color.use[[u]],lty=2,lwd=1.5  )
			lines((seq(xlim[1],xlim[2],by=0.1)),y=lin.non2,lty=2,lwd=1.5  )

			#Spatial mechanisms add this much
			width_y2 = C1*var_mu_Us[(tol2[g]),u] 
			#Strongest negative impact (this is only accurate if data include case of complete overlap)
			#min_overlap2 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
			overlap_x2 = -min_overlap1+((-C2*cov_e_mu_Us[(tol2[g]),u])+cov_lam_vc[(tol2[g]),u])+0.08
			non.x2 = -mlin.non2/(-1-(width_y2/overlap_x2) )
			non.y2 = -non.x2+ mlin.non2 
			#points(non.x2,non.y2 )
			points(overlap_x2,width_y2+j1[ug],col=color.use[[u]])
			#segments(non.x2,non.y2,overlap_x2,width_y2,col=color.use[[u]],lty=2)			
			#print(overlap_x2)
			#segments(non.x1,non.y1,new.x2,new.y2,col=color.use[[u]],lty=2)
			arrows(overlap_x1,width_y1+j1[ug],overlap_x2,width_y2+j1[ug],col=color.use[[u]],length=0.1,lwd=2)
		} 
}

dev.off()


################################################################
#Figure 3B: Non-Arch (linear) figures showing how the igr of each 
# scenario moves according to the benefit of either range width 
# or range overlap.
# This version for the full 3spp picture. 
# Components STANDARDIZED BY MEAN coefficients
# Figures are a SINGLE PLOT
################################################################
#Set these variables manually (hopefully you kept good notes)
iconfig=1
ngens=2
ngenst=iconfig+ngens
nspp=3
scale=10
#tol2=c(100,100,100,100)
#tol2=c(5,5,5,5)
tol2=c(2,2,2,2)


#xlims = matrix(c(-0.0, 0.15, -0.0,2.3,-0.0,2.3, -0.0,2.3),4,2,byrow=T)
#x.axx = matrix(c(0.0,0.15,0.05,0.2,2.3,0.4,0.2,2.3,0.4, 0.2,2.3,0.4 ),4,3,byrow=T)


fig.name = paste("all3_calanda_widthOverlap_igr12.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

#Set this up as a panel of 4, or see below for individual files. 
par(mfrow=c(1,nscen),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
#This is a constant to scale the range overlap
min_overlap1=0
xlims=0
ylims=0

#Find the minimum value across the scenarios
for( g in 1:nscen){ 
	for(u in 1:nspp){
		#Identify the variable names by scenario, then assign them to the right
		#variable for plotting
		cov_e_mu_Us.name=paste(prefix.list[[g]],variable.list[[4]] ,sep="")
		cov_lam_vc.name=paste(prefix.list[[g]],variable.list[[27]] ,sep="")
		l1.name=paste(prefix.list[[g]],variable.list[[1]] ,sep="")
		D.name = paste(prefix.list[[g]],variable.list[[2]] ,sep="")
		var_mu_Us.name=paste(prefix.list[[g]],variable.list[[3]] ,sep="")

		assign("var_mu_Us",eval(as.name(var_mu_Us.name)) )
		assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name)) )
		assign("cov_lam_vc",eval(as.name(cov_lam_vc.name)) )
		assign("l1",eval(as.name(l1.name)) )
		assign("D",eval(as.name(D.name)) )

		#Create a mean standardization coefficient to simplify presentation. 
		#C1 =  median((l1[1,u]/(D[1,u])^3), na.rm=T)
		#C2 =  median(l1[1,u]/(D[1,u])^2, na.rm=T )

		C1 =  median((l1[1,u]/(D[1,u])), na.rm=T)
		C2 =  median(l1[1,u]/(D[1,u]), na.rm=T )

		min_overlap1 = min(c (min_overlap1,(-C2*cov_e_mu_Us[tol2[u],u]+cov_lam_vc[tol2[u],u]),(-C2*cov_e_mu_Us[1,u]+cov_lam_vc[1,u])),na.rm=T)
		xlims = max(c (xlims,(-C2*cov_e_mu_Us[tol2[u],u]+cov_lam_vc[tol2[u],u]),(-C2*cov_e_mu_Us[1,u]+cov_lam_vc[1,u])),na.rm=T)
		ylims = max(c (ylims,(var_mu_Us[tol2[u],u]),(var_mu_Us[tol2[u],u])),na.rm=T)

	}
}
xlims = matrix(c(0,round(-min_overlap1+xlims+0.2,1)),1,2)
ylims = matrix(c(0,round(ylims+0.1,1)),1,2)
x.axx = matrix( c(xlims[1],xlims[2],(( xlims[2]-xlims[1]) / 5) ),1,3)
y.axx = matrix( c(ylims[1],ylims[2],(( ylims[2]-ylims[1]) / 5) ),1,3)


#Now make the plots: 
for( g in 1:nscen){ 

	color.use=list("black","red","blue")
	up.low=list("upper","lower")

	#Or just switch between devices for plotting
	#dev.set(g)
	#par(mfrow=c(1,4), mar=c(5,6,6,2))

	ylim=ylims
	#xlim=c(0,0.8)
	xlim=xlims[g,]
	#plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Range width", xlab="Range overlap",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3 )
	if (g ==1){
	plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Benefit of Competitor Aggregation", xlab="",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ann=F,yaxt='n',xaxt='n' )
	axis(side=1, labels=T, at=seq(x.axx[g,1], x.axx[g,2],x.axx[g,3]))
	axis(side=2, labels=T, at=seq(y.axx[g,1], y.axx[g,2],y.axx[g,3]))

	} else { 
	plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="", xlab="",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ann=F,yaxt='n',xaxt='n')
	axis(side=1, labels=T, at=seq(x.axx[g,1], x.axx[g,2],x.axx[g,3]))

		} 
	#mlin.non1 =0.3229569
	#lines((seq(xlim[1],xlim[2],by=0.1)),lin.non1,lty=2,lwd=1.5  )

	for (u in 1:nspp){ 
			
			

			#Identify the variable names by scenario, then assign them to the right
			#variable for plotting
			gr.name=paste(prefix.list[[g]],variable.list[[10]] ,sep="")
			var_mu_Us.name=paste(prefix.list[[g]],variable.list[[3]] ,sep="")
			cov_e_mu_Us.name=paste(prefix.list[[g]],variable.list[[4]] ,sep="")
			cov_lam_vc.name=paste(prefix.list[[g]],variable.list[[27]] ,sep="")
			#For scaling
			l1.name=paste(prefix.list[[g]],variable.list[[1]] ,sep="")
			w.eq.name = paste(prefix.list[[g]],variable.list[[12]] ,sep="")
			D.name = paste(prefix.list[[g]],variable.list[[2]] ,sep="")

			assign("gr",eval(as.name(gr.name)) )
			assign("var_mu_Us",eval(as.name(var_mu_Us.name)) )
			assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name)) )
			assign("cov_lam_vc",eval(as.name(cov_lam_vc.name)) )
			assign("l1",eval(as.name(l1.name)) )
			assign("D",eval(as.name(D.name)) )
			assign("w.eq",eval(as.name(w.eq.name))  )

			#This line represents persistence in the absence of spatial mechanisms
			mlin.non1 =1-(l1[iconfig,u]/D[iconfig,u]+sr[u])-min_overlap1
			#mlin.non1 =0.3229569
			lin.non1= -(seq(xlim[1],xlim[2],by=0.05))+mlin.non1
			lines((seq(xlim[1],xlim[2],by=0.05)),lin.non1,col=color.use[[u]],lty=2,lwd=1.5  )


			#Create a mean standardization coefficient to simplify presentation. 
			#C1 =  median((l1[,u]/(D[,u])^3), na.rm=T)
			#C2 =  median(l1[,u]/(D[,u])^2, na.rm=T )

			#C1 = (l1[1,u]/(D[1,u])^3)
			#C2 = l1[1,u]/(D[1,u])^2

			C1 =  median((l1[1,u]/(D[1,u])), na.rm=T)
			C2 =  median(l1[1,u]/(D[1,u]), na.rm=T )


			#Spatial mechanisms add this much
			width_y1 = C1*var_mu_Us[iconfig,u] 
			#Strongest negative impact (this is only accurate if data include case of complete overlap)
			#min_overlap1 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
			overlap_x1 = -min_overlap1+((-C2*cov_e_mu_Us[iconfig,u])+cov_lam_vc[iconfig,u])
			#Point where lines intercept: x = (b1-b2)/2
			#non.x1 = (mlin.non1 - (width_y1 - overlap_x1) )/2
			#non.y1 = -non.x1+ mlin.non1 
			#Intersection through origin instead:
			non.x1 = -mlin.non1/(-1-(width_y1/overlap_x1) )
			non.y1 = -non.x1+ mlin.non1 
			#points(non.x1,non.y1,col=color.use[[u]], pch=16)
			points(overlap_x1,width_y1,col=color.use[[u]] )
			#segments(non.x1,non.y1,overlap_x1,width_y1,col=color.use[[u]],lty=2)			

			#Now the later point
			mlin.non2 =1-(l1[(tol2[g]),u]/D[(tol2[g]),u]+sr[u])-min_overlap1
			#mlin.non2 =0.3229569
			lin.non2= -(seq(xlim[1],xlim[2],by=0.1))+mlin.non2
			lines((seq(xlim[1],xlim[2],by=0.1)),lin.non2,col=color.use[[u]] ,lty=2,lwd=1.5  )
			#lines((seq(xlim[1],xlim[2],by=0.1)),lin.non2,lty=2,lwd=1.5  )

			#Spatial mechanisms add this much
			width_y2 = C1*var_mu_Us[(tol2[g]),u] 
			#Strongest negative impact (this is only accurate if data include case of complete overlap)
			#min_overlap2 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
			overlap_x2 = -min_overlap1+((-C2*cov_e_mu_Us[(tol2[g]),u])+cov_lam_vc[(tol2[g]),u])
			non.x2 = -mlin.non2/(-1-(width_y2/overlap_x2) )
			non.y2 = -non.x2+ mlin.non2 
			#points(non.x2,non.y2 )
			points(overlap_x2,width_y2,col=color.use[[u]])
			#segments(non.x2,non.y2,overlap_x2,width_y2,col=color.use[[u]],lty=2)			

			#segments(non.x1,non.y1,new.x2,new.y2,col=color.use[[u]],lty=2)
			arrows(overlap_x1,width_y1,overlap_x2,width_y2,col=color.use[[u]],length=0.1,lwd=2)
		} 
}

dev.off()

################################################################
#Figure ??: 
#Make a fancy multi-figure image showing the environmental novelty, 
#physical distance, and contribution of variables to the novelty.
################################################################

layout.matrix=matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)
layout(mat = layout.matrix,
       heights = c(5,1,1), # Heights of the two rows
       widths = c(10,1)) # Widths of the two columns

layout.show(6)

nlevel=64
par(oma = c(3,3,3,3))
par( mar = c(0.5,4,0,4) )
image(x=xx$elevation, y = 1:length(env.ind), var_dist, ylab="", yaxt='n',xaxt='n',col=viridis(nlevel))
axis(2, at=1:(length(env.ind)), labels = colnames(var_dist),cex.axis=0.9)

par( mar = c(0.5,4,0,4) )
image(x=xx$elevation, y = 1, env_analogue[[2]],ylab="", yaxt='n',xaxt='n',col=viridis(nlevel))
axis(2, at=1, labels = c("total"),cex.axis=0.9)

par( mar = c(0.5,4,0,4) )
image(x=xx$elevation, y = 1, an_dist, xlab="", ylab="",yaxt='n',col=viridis(nlevel))
axis(2, at=1, labels = c("distance"),cex.axis=0.9)

par( mar = c(4,1,4,1) )
image(1,(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), t(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), ylab="",xaxt='n',col=viridis(nlevel))
 
par( mar = c(1,1,1,1) )
image(1, (seq(min(env_analogue[[2]]),max(env_analogue[[2]]),max(env_analogue[[2]])/nlevel)), t(seq(min(env_analogue[[2]]),max(env_analogue[[2]]),max(env_analogue[[2]])/nlevel)), ylab="",xaxt='n',col=viridis(nlevel))
 
par( mar = c(1,1,1,1) )
image(1, (seq(min(an_dist),max(an_dist),max(an_dist)/nlevel)), t(seq(min(an_dist),max(an_dist),max(an_dist)/nlevel)), ylab="",xaxt='n',col=viridis(nlevel))
 
