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
library(viridis)
library(scales)
#==============================================================================
# Load Data
#==============================================================================

#This requires user input!

# variable.list=list("l1.all", "D.all", "var_mu_Us.all","cov_e_mu_Us.all",
# "cov_lam_vc.all", "cov_lam_vc2.all", "Elam1.all", "Elam2.all", "gr1.n.all", "gr1.all", "y.full.all", 
# "w.eq.full.all","ldg.sim.all", "ldgIMP.sim.all", "covIMP.sim.all", "Frs.all", "env_analogue_all")

variable.list=list("l1", "D", "var_mu_Us","cov_e_mu_Us",
"cov_lam_vc", "cov_lam_vc2", "Elam1", "Elam2", "gr1.n", "gr1", "y.full", 
"w.eq.full","ldg.sim", "ldgIMP.sim", "covIMP.sim", "Frs", "env_analogue_all",
"sim.impacts", "sim.impactsIMP","covIMP.impacts","lded","nf","is.final")



#List of scenario file names 
#file.name.list= list("calanda_igrsites12it2_multi_2017_waldmean.var")
#file.name.list= list("calanda_igrsites12it3_global_2017_waldmean.var")
file.name.list= list("calanda_igrsites12it3_multi_2017_waldmean.var")


#In order to deal with each of the scenarios simultaneously, it is necessary
#to systematically rename va"calanda_repeat1_A_multi_2017","calanda_repeat1_B_multi_2017",riables based on their scenario. 
#Assosciated prefixes for renaming

prefix.list= list("sc1")

#Rename the files from each scenario
var.length=length(variable.list)
nscen = length(file.name.list)


for (g in 1:nscen){
	load(file.name.list[[g]])
	for( f in 1:(var.length)) {      
		new.name=paste(prefix.list[[g]],variable.list[[f]],sep="")	
		assign(new.name,eval(as.name(variable.list[[f]])))

	}	
}

# #Take each file (which is a list) and turn it into an array
# for( f in 1:(var.length-1)) { 
# 	tmp.var=NULL
# 	for (g in 1:(nscen)){   
# 		tmp.name=paste(prefix.list[[g]],variable.list[[f]],sep="")	
# 		tmp.var=rbind(tmp.var,eval(as.name(tmp.name)))	
# 	}	

# 	tmp.name2=paste("comb",variable.list[[f]],sep="")	
# 	assign(tmp.name2,tmp.var)

# }

# #Combine the intrinsic ranges
# combFrs = abind(startFrs, shiftFrs, along=1)

# prefix.list= list("comb")
# nscen=1


#==============================================================================
# Figure: Intrinsic ranges, before and after climate change
#==============================================================================
source("./range_coexistence_functionsD.R")

nspp=3 #Number of species
ns=dim(sc1Frs)[2] #Lattice width/height -- This should be even
ngens=6#Number of generations for environmental change
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
plot(combFrs[1,,1],t="l",xlab="Location", ylab="Fitness", cex.lab=1.7,cex.axis=1.2,ylim=c(0,70),xlim=c(0,ns),
		col=color.use[1], xaxs="i",xaxt='n',yaxs="i",yaxt='n',bty="n",lwd = 3,lty=2)
abline(h=0)
for(j in 2:nspp) {
	lines( sc1Frs.all[[1]][1,,j],col = color.use[j],lwd=3,lty=2)
}

#plot(Frs[ngenst,,1],t="l",xlab="Location", ylab="Fitness", cex.lab=1.7,cex.axis=1.2,ylim=c(0,3),
#		col=color.use[1],xaxs="i",xaxt='n',yaxs="i",yaxt='n',bty="n",lwd=3)
abline(h=0)
for(j in 1:nspp) {
	lines( sc1Frs.all[[ngenst]][1,,j],col = color.use[j],lwd=3)
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

#
xx=cbind(xx1,yx,xx2)
colnames(xx)=c('elevation','year','e2')
xx=data.frame(xx)
xx$elevation = as.matrix(xx$elevation)
xx$year = as.matrix(xx$year)
xx$e2 = as.matrix(xx$e2)

#Temp increments
mstart = 0
mstop = 6
minc = 0.5
mtot= ceiling((mstop-mstart)/minc)
m.labels=as.character(seq(mstart,mstop,minc))
#=============================================================================

epoints = match( elevations,xx$elevation)

fig.name = paste("seeds_kriged_igrsites12it3.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

par(mfrow=c(2,1),mai= c( 0.0, 0.2, 0.0, 0.2), omi=c(0.5,0.75,0.5,0.75)) #,mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))

#####
# Panel 1: Intrinsic ranges
#####

#par(mfrow=c(1,1),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
ylim=c(-1,60)
x.txt = c(1900,1200,2700)
#Current: 
plot(xx$elevation, sc1Frs[1,,1], t="l", ylab="Seeds (per-capita, Kriged) ",  xaxt="n", xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim)
abline(v=c(elevations),lwd=10,col="grey80")
lines(xx$elevation,sc1Frs[1,,1])
lines(xx$elevation, sc1Frs[1,,2],col="red")
lines(xx$elevation, sc1Frs[1,,3],col="blue")
text(x=x.txt[1], y=sc1Frs[1,x.txt[1],1] , labels=m.labels[1])
text(x=x.txt[2], y=sc1Frs[1,x.txt[2],2] , labels=m.labels[1],col="red")
text(x=x.txt[3], y=sc1Frs[1,x.txt[3],3] , labels=m.labels[1],col="blue")

# points(xx$elevation[epoints], sc1Frs.all[[1]][1,epoints,1])
# points(xx$elevation[epoints], sc1Frs.all[[1]][1,epoints,2],col="red")
# points(xx$elevation[epoints], sc1Frs.all[[1]][1,epoints,3],col="blue")

#Future

for(n in seq(3,length(m.labels),2)){
lines(xx$elevation, sc1Frs[n,,1],lty=2)
text(x=x.txt[1], y=sc1Frs[n,x.txt[1],1] , labels=m.labels[(n)])
lines(xx$elevation, sc1Frs[n,,2],col="red",lty=2)
text(x=x.txt[2], y=sc1Frs[n,x.txt[2],2] , labels=m.labels[(n)],col="red")
lines(xx$elevation, sc1Frs[n,,3],col="blue",lty=2)
text(x=x.txt[3], y=sc1Frs[n,x.txt[3],3] , labels=m.labels[(n)],col="blue")
# points(xx$elevation[epoints], sc1Frs.all[[n]][1,epoints,1])
# points(xx$elevation[epoints], sc1Frs.all[[n]][1,epoints,2],col="red")
# points(xx$elevation[epoints], sc1Frs.all[[n]][1,epoints,3],col="blue")
}
# f1 = flower_prob_act*flower_act
# points(xx$elevation[epoints], f1[,1],pch=3)
# points(xx$elevation[epoints], f1[,2],col="red",pch=3)
# points(xx$elevation[epoints], f1[,3],col="blue",pch=3)

#####
# Panel 2: Realized ranges
#####

ylim=c(-1,150)
x.txt = c(1900,1200,2700)
#Current: 
plot(xx$elevation, sc1nf[1,,1], t="l", ylab="Seeds (per-capita, Kriged) ", xlab="Elevation", xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim)
abline(v=c(elevations),lwd=10,col="grey80")
lines(xx$elevation,sc1nf[1,,1])
lines(xx$elevation, sc1nf[1,,2],col="red")
lines(xx$elevation, sc1nf[1,,3],col="blue")
text(x=x.txt[1], y=sc1nf[1,x.txt[1],1] , labels=m.labels[1])
text(x=x.txt[2], y=sc1nf[1,x.txt[2],2] , labels=m.labels[1],col="red")
text(x=x.txt[3], y=sc1nf[1,x.txt[3],3] , labels=m.labels[1],col="blue")

# points(xx$elevation[epoints], sc1nf.all[[1]][1,epoints,1])
# points(xx$elevation[epoints], sc1nf.all[[1]][1,epoints,2],col="red")
# points(xx$elevation[epoints], sc1nf.all[[1]][1,epoints,3],col="blue")

#Future

for(n in seq(3,length(m.labels),2)){
lines(xx$elevation, sc1nf[n,,1],lty=2)
text(x=x.txt[1], y=sc1nf[n,x.txt[1],1] , labels=m.labels[(n)])
lines(xx$elevation, sc1nf[n,,2],col="red",lty=2)
text(x=x.txt[2], y=sc1nf[n,x.txt[2],2] , labels=m.labels[(n)],col="red")
lines(xx$elevation, sc1nf[n,,3],col="blue",lty=2)
text(x=x.txt[3], y=sc1nf[n,x.txt[3],3] , labels=m.labels[(n)],col="blue")
# points(xx$elevation[epoints], sc1nf.all[[n]][1,epoints,1])
# points(xx$elevation[epoints], sc1nf.all[[n]][1,epoints,2],col="red")
# points(xx$elevation[epoints], sc1nf.all[[n]][1,epoints,3],col="blue")
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
ngens=6
ngenst=iconfig+ngens
nspp=3
scale=10
#tol2=c(100,100,100,100)
#tol2=c(5,5,5,5)
tol2=c(2,2,2,2)

nscen=6


#xlims = matrix(c(-0.0, 0.15, -0.0,2.3,-0.0,2.3, -0.0,2.3),4,2,byrow=T)
#x.axx = matrix(c(0.0,0.15,0.05,0.2,2.3,0.4,0.2,2.3,0.4, 0.2,2.3,0.4 ),4,3,byrow=T)


fig.name = paste("calanda_widthOverlap_igrsites12it3.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

#Set this up as a panel of 4, or see below for individual files. 
par(mfrow=c(nscen/2,2),mai= c( 0.25, 0.2, 0.25, 0.2), omi=c(0.5,0.75,0.5,0.75)) #,mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
#This is a constant to scale the range overlap
min_overlap1=0

xlimsa=0
ylimsa=0

#Make different limits per plot
xlims = matrix(0,ngenst,2)
ylims = matrix(0,ngenst,2)
x.axx = matrix(0,ngenst,3)
y.axx = matrix(0,ngenst,3)

#Find the minimum value across the scenarios
for( g in 1:1){ 


	#Identify the variable names by scenario, then assign them to the right
	#variable for plotting
	cov_e_mu_Us.name=paste(prefix.list[[g]],variable.list[[4]] ,sep="")
	cov_lam_vc.name=paste(prefix.list[[g]],variable.list[[15]] ,sep="")
	l1.name=paste(prefix.list[[g]],variable.list[[1]] ,sep="")
	D.name = paste(prefix.list[[g]],variable.list[[2]] ,sep="")
	var_mu_Us.name=paste(prefix.list[[g]],variable.list[[3]] ,sep="")

	assign("var_mu_Us",eval(as.name(var_mu_Us.name)))
	assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name))  )
	assign("cov_lam_vc",eval(as.name(cov_lam_vc.name))  )
	assign("l1",eval(as.name(l1.name)) )
	assign("D",eval(as.name(D.name))  )

	for(s in seq(3,length(m.labels),2)) {
		
		min_overlap1=0
		xlimsa=0
		ylimsa=0

		for(u in 1:nspp){


			#Create a mean standardization coefficient to simplify presentation. 
			#C1 =  median((l1[1,u]/(D[1,u])^3), na.rm=T)
			#C2 =  median(l1[1,u]/(D[1,u])^2, na.rm=T )

			C1 =  median((l1[1,u]/(D[1,u])), na.rm=T)
			C2 =  median(l1[1,u]/(D[1,u]), na.rm=T )

			min_overlap1 = min(c (min_overlap1,(-C2*cov_e_mu_Us[s,u]+cov_lam_vc[s,u]),(-C2*cov_e_mu_Us[1,u]+cov_lam_vc[1,u])),na.rm=T)
			xlimsa = max(c (xlimsa,(-C2*cov_e_mu_Us[s,u]+cov_lam_vc[s,u]),(-C2*cov_e_mu_Us[1,u]+cov_lam_vc[1,u])),na.rm=T)
			ylimsa = max(c (ylimsa,(C1*var_mu_Us[s,u]),(C1*var_mu_Us[1,u])),na.rm=T)

		}
	
	xlims[s,2] = round(-min_overlap1+xlimsa+0.2,2)
	ylims[s,2] = round(ylimsa+0.1,1)



	}
	x.axx= round(matrix( c(xlims[,1],xlims[,2],(( xlims[,2]-xlims[,1]) / 3) ),ngenst,3),2)
	y.axx= round(matrix( c(ylims[,1],ylims[,2],(( ylims[,2]-ylims[,1]) / 3) ),ngenst,3),2)

}

#Uncomment this to make all of the limits the same
# xlims = matrix(c(0,round(-min_overlap1+xlims+0.2,1)),ngenst,2,byrow=T)
# ylims = matrix(c(0,round(ylims+0.1,1)),ngenst,2,byrow=T)
# x.axx = matrix( c(xlims[,1],xlims[,2],(( xlims[,2]-xlims[,1]) / 5) ),ngenst,3)
# y.axx = matrix( c(ylims[,1],ylims[,2],(( ylims[,2]-ylims[,1]) / 5) ),ngenst,3)

#Now make the plots: 

#Identify the variable names by scenario, then assign them to the right
#variable for plotting
gr.name=paste(prefix.list[[g]],variable.list[[13]] ,sep="")
cov_e_mu_Us.name=paste(prefix.list[[g]],variable.list[[4]] ,sep="")
cov_lam_vc.name=paste(prefix.list[[g]],variable.list[[15]] ,sep="")
var_mu_Us.name=paste(prefix.list[[g]],variable.list[[3]] ,sep="")

#For scaling
l1.name=paste(prefix.list[[g]],variable.list[[1]] ,sep="")
D.name = paste(prefix.list[[g]],variable.list[[2]] ,sep="")
w.eq.name = paste(prefix.list[[g]],variable.list[[12]] ,sep="")

assign("gr",eval(as.name(gr.name)) )
assign("var_mu_Us",eval(as.name(var_mu_Us.name)) )
assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name)))
assign("cov_lam_vc",eval(as.name(cov_lam_vc.name)) )
assign("l1",eval(as.name(l1.name)) )
assign("D",eval(as.name(D.name)) )
assign("w.eq",eval(as.name(w.eq.name))[[s]]  )


for( s in seq(3,length(m.labels),2)){ 

	color.use=list("black","red","blue")
	up.low=list("upper","lower")

	#Or just switch between devices for plotting
	#dev.set(g)
	#par(mfrow=c(1,4), mar=c(5,6,6,2))

	ylim=ylims[s,]
	#xlim=c(0,0.8)
	xlim=xlims[s,]
	#plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Range width", xlab="Range overlap",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3 )
	if ((s-1)%%3 == 0){
	plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Benefit of Competitor Aggregation", xlab="",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ann=F,yaxt='n',xaxt='n' )
	axis(side=1, labels=T, at=seq(x.axx[s,1], x.axx[s,2],x.axx[s,3]))
	axis(side=2, labels=T, at=seq(y.axx[s,1], y.axx[s,2],y.axx[s,3]))

	} else { 
	plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="", xlab="",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ann=F,yaxt='n',xaxt='n')
	axis(side=1, labels=T, at=seq(x.axx[s,1], x.axx[s,2],x.axx[s,3]))
	axis(side=2, labels=T, at=seq(y.axx[s,1], y.axx[s,2],y.axx[s,3]))
		} 
	#mlin.non1 =0.3229569
	#lines((seq(xlim[1],xlim[2],by=0.1)),lin.non1,lty=2,lwd=1.5  )

	for (u in 1:nspp){ 
			
			
			#This line represents persistence in the absence of spatial mechanisms
			mlin.non1 =1-(l1[iconfig,u]/D[iconfig,u]+sr[u])-min_overlap1
			#mlin.non1 =0.3229569
			lin.non1= -(seq(xlim[1],xlim[2],by=0.05))+mlin.non1
			lines((seq(xlim[1],xlim[2],by=0.05)),lin.non1,col=alpha(rgb(t(matrix(col2rgb(color.use[[u]]))),maxColorValue=255),0.4),lty=5,lwd=1.5  )


			#Create a mean standardization coefficient to simplify presentation. 
			#C1 =  median((l1[,u]/(D[,u])^3), na.rm=T)
			#C2 =  median(l1[,u]/(D[,u])^2, na.rm=T )

			#C1 = (l1[1,u]/(D[1,u])^3)
			#C2 = l1[1,u]/(D[1,u])^2

			C1 =  median((l1[iconfig,u]/(D[iconfig,u])), na.rm=T)
			C2 =  median(l1[iconfig,u]/(D[iconfig,u]), na.rm=T )


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
			mlin.non2 =1-(l1[s,u]/D[s,u]+sr[u])-min_overlap1
			#mlin.non2 =0.3229569
			lin.non2= -(seq(xlim[1],xlim[2],by=0.1))+mlin.non2
			lines((seq(xlim[1],xlim[2],by=0.1)),lin.non2,col=alpha(rgb(t(matrix(col2rgb(color.use[[u]]))),maxColorValue=255),0.4),lty=3,lwd=1.5  )
			#lines((seq(xlim[1],xlim[2],by=0.1)),lin.non2,lty=2,lwd=1.5  )

			#Spatial mechanisms add this much
			width_y2 = C1*var_mu_Us[s,u]
			#Strongest negative impact (this is only accurate if data include case of complete overlap)
			#min_overlap2 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
			overlap_x2 = -min_overlap1+((-C2*cov_e_mu_Us[s,u])+cov_lam_vc[s,u])
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
#Figure 3B: Non-Arch (linear) figures showing how the igr of each 
# scenario moves according to the benefit of either range width 
# or range overlap.
# This version for the full 3spp picture. 
# Components STANDARDIZED BY MEAN coefficients
# Figures are a SINGLE PLOT
# This gives a sequence of arrows where the initial and final points
# are between subsequent time increments, instead of from beginning 
# to the changed conditions. 
################################################################
#Set these variables manually (hopefully you kept good notes)
iconfig=1
ngens=6
ngenst=iconfig+ngens
nspp=3
scale=10
#tol2=c(100,100,100,100)
#tol2=c(5,5,5,5)
tol2=c(2,2,2,2)

nscen=6


#xlims = matrix(c(-0.0, 0.15, -0.0,2.3,-0.0,2.3, -0.0,2.3),4,2,byrow=T)
#x.axx = matrix(c(0.0,0.15,0.05,0.2,2.3,0.4,0.2,2.3,0.4, 0.2,2.3,0.4 ),4,3,byrow=T)


fig.name = paste("calanda_widthOverlap_repeat_igr3b_deltaP.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

#Set this up as a panel of 4, or see below for individual files. 
par(mfrow=c(nscen/2,2),mai= c( 1, 0.2, 0.2, 0.2), omi=c(0.5,0.75,0.5,0.75)) #,mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
#This is a constant to scale the range overlap
min_overlap1=0

xlimsa=0
ylimsa=0

#Make different limits per plot
xlims = matrix(0,ngenst,2)
ylims = matrix(0,ngenst,2)
x.axx = matrix(0,ngenst,3)
y.axx = matrix(0,ngenst,3)

#Find the minimum value across the scenarios
for( g in 1:1){ 


	#Identify the variable names by scenario, then assign them to the right
	#variable for plotting
	cov_e_mu_Us.name=paste(prefix.list[[g]],variable.list[[4]] ,sep="")
	cov_lam_vc.name=paste(prefix.list[[g]],variable.list[[15]] ,sep="")
	l1.name=paste(prefix.list[[g]],variable.list[[1]] ,sep="")
	D.name = paste(prefix.list[[g]],variable.list[[2]] ,sep="")
	var_mu_Us.name=paste(prefix.list[[g]],variable.list[[3]] ,sep="")

	assign("var_mu_Us",eval(as.name(var_mu_Us.name)))
	assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name))  )
	assign("cov_lam_vc",eval(as.name(cov_lam_vc.name))  )
	assign("l1",eval(as.name(l1.name)) )
	assign("D",eval(as.name(D.name))  )

	for(s in 1:ngenst) {
		for(u in 1:nspp){


			#Create a mean standardization coefficient to simplify presentation. 
			#C1 =  median((l1[1,u]/(D[1,u])^3), na.rm=T)
			#C2 =  median(l1[1,u]/(D[1,u])^2, na.rm=T )

			C1 =  median((l1[[1]][u]/(D[[1]][u])), na.rm=T)
			C2 =  median(l1[[1]][u]/(D[[1]][u]), na.rm=T )

			min_overlap1 = min(c (min_overlap1,(-C2*cov_e_mu_Us[[s]][u]+cov_lam_vc[[s]][u]),(-C2*cov_e_mu_Us[[1]][u]+cov_lam_vc[[1]][u])),na.rm=T)
			xlimsa = max(c (xlims,(-C2*cov_e_mu_Us[[s]][u]+cov_lam_vc[[s]][u]),(-C2*cov_e_mu_Us[[1]][u]+cov_lam_vc[[1]][u])),na.rm=T)
			ylimsa = max(c (ylims,(var_mu_Us[[s]][u]),(var_mu_Us[[s]][u])),na.rm=T)

		}
	
	xlims[s,2] = round(-min_overlap1+xlimsa+0.2,1)
	ylims[s,2] = round(ylimsa+0.1,1)

	}
	x.axx= matrix( c(xlims[,1],xlims[,2],(( xlims[,2]-xlims[,1]) / 5) ),ngenst,3)
	y.axx= matrix( c(ylims[,1],ylims[,2],(( ylims[,2]-ylims[,1]) / 5) ),ngenst,3)

}

#Uncomment this to make all of the limits the same
# xlims = matrix(c(0,round(-min_overlap1+xlims+0.2,1)),ngenst,2,byrow=T)
# ylims = matrix(c(0,round(ylims+0.1,1)),ngenst,2,byrow=T)
# x.axx = matrix( c(xlims[,1],xlims[,2],(( xlims[,2]-xlims[,1]) / 5) ),ngenst,3)
# y.axx = matrix( c(ylims[,1],ylims[,2],(( ylims[,2]-ylims[,1]) / 5) ),ngenst,3)

#Now make the plots: 

#Identify the variable names by scenario, then assign them to the right
#variable for plotting
gr.name=paste(prefix.list[[g]],variable.list[[13]] ,sep="")
cov_e_mu_Us.name=paste(prefix.list[[g]],variable.list[[4]] ,sep="")
cov_lam_vc.name=paste(prefix.list[[g]],variable.list[[15]] ,sep="")
var_mu_Us.name=paste(prefix.list[[g]],variable.list[[3]] ,sep="")

#For scaling
l1.name=paste(prefix.list[[g]],variable.list[[1]] ,sep="")
D.name = paste(prefix.list[[g]],variable.list[[2]] ,sep="")
w.eq.name = paste(prefix.list[[g]],variable.list[[12]] ,sep="")

assign("gr",eval(as.name(gr.name)) )
assign("var_mu_Us",eval(as.name(var_mu_Us.name)) )
assign("cov_e_mu_Us",eval(as.name(cov_e_mu_Us.name)))
assign("cov_lam_vc",eval(as.name(cov_lam_vc.name)) )
assign("l1",eval(as.name(l1.name)) )
assign("D",eval(as.name(D.name)) )
assign("w.eq",eval(as.name(w.eq.name))[[s]]  )


for( s in 2:(ngenst)){ 

	iconfig = s-1
	color.use=list("black","red","blue")
	up.low=list("upper","lower")

	#Or just switch between devices for plotting
	#dev.set(g)
	#par(mfrow=c(1,4), mar=c(5,6,6,2))

	ylim=ylims[s,]
	#xlim=c(0,0.8)
	xlim=xlims[s,]
	#plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Range width", xlab="Range overlap",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3 )
	if (s ==1){
	plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Benefit of Competitor Aggregation", xlab="",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ann=F,yaxt='n',xaxt='n' )
	axis(side=1, labels=T, at=seq(x.axx[s,1], x.axx[s,2],x.axx[s,3]))
	axis(side=2, labels=T, at=seq(y.axx[s,1], y.axx[s,2],y.axx[s,3]))

	} else { 
	plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="", xlab="",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ann=F,yaxt='n',xaxt='n')
	axis(side=1, labels=T, at=seq(x.axx[s,1], x.axx[s,2],x.axx[s,3]))
		} 
	#mlin.non1 =0.3229569
	#lines((seq(xlim[1],xlim[2],by=0.1)),lin.non1,lty=2,lwd=1.5  )

	for (u in 1:nspp){ 
			
			#This line represents persistence in the absence of spatial mechanisms
			mlin.non1 =1-(l1[[iconfig]][u]/D[[iconfig]][u]+sr[u])-min_overlap1
			#mlin.non1 =0.3229569
			lin.non1= -(seq(xlim[1],xlim[2],by=0.05))+mlin.non1
			lines((seq(xlim[1],xlim[2],by=0.05)),lin.non1,col=color.use[[u]],lty=2,lwd=1.5  )


			#Create a mean standardization coefficient to simplify presentation. 
			#C1 =  median((l1[,u]/(D[,u])^3), na.rm=T)
			#C2 =  median(l1[,u]/(D[,u])^2, na.rm=T )

			#C1 = (l1[1,u]/(D[1,u])^3)
			#C2 = l1[1,u]/(D[1,u])^2

			C1 =  median((l1[[1]][u]/(D[[1]][u])), na.rm=T)
			C2 =  median(l1[[1]][u]/(D[[1]][u]), na.rm=T )


			#Spatial mechanisms add this much
			width_y1 = C1*var_mu_Us[[iconfig]][u] 
			#Strongest negative impact (this is only accurate if data include case of complete overlap)
			#min_overlap1 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
			overlap_x1 = -min_overlap1+((-C2*cov_e_mu_Us[[iconfig]][u])+cov_lam_vc[[iconfig]][u])
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
			mlin.non2 =1-(l1[[(s)]][u]/D[[(s)]][u]+sr[u])-min_overlap1
			#mlin.non2 =0.3229569
			lin.non2= -(seq(xlim[1],xlim[2],by=0.1))+mlin.non2
			lines((seq(xlim[1],xlim[2],by=0.1)),lin.non2,col=color.use[[u]] ,lty=2,lwd=1.5  )
			#lines((seq(xlim[1],xlim[2],by=0.1)),lin.non2,lty=2,lwd=1.5  )

			#Spatial mechanisms add this much
			width_y2 = C1*var_mu_Us[[(s)]][u] 
			#Strongest negative impact (this is only accurate if data include case of complete overlap)
			#min_overlap2 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
			overlap_x2 = -min_overlap1+((-C2*cov_e_mu_Us[[(s)]][u])+cov_lam_vc[[(s)]][u])
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
###############################################################


#ldg.sim.a=t(matrix(unlist(sc1ldg.sim.all[1:13]),3,13))
env_analogue_a2 = matrix(0,3501,13)
for(a in 1:13){ env_analogue_a2[,a]=matrix(unlist(sc1env_analogue_all[[a]][[2]]),3501,1) }
env_analogue_a3 = matrix(0,3501,13)
for(a in 1:13){ env_analogue_a3[,a]=matrix(unlist(sc1env_analogue_all[[a]][[3]]),3501,1) }

c1 = colMeans(env_analogue_a2)
lmm1=lm(ldg.sim[,1]~c1)
lmm2=lm(ldg.sim[,2]~c1)
lmm3=lm(ldg.sim[,3]~c1)

fig.name = paste("calanda_noveltylm_ldg3.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

plot(c1, ldg.sim[,1],ylim=c(0,15),ylab="LDG",xlab ="average novelty")
points(c1, ldg.sim[,2],col="red")
points(c1, ldg.sim[,3],col="blue")
#lines(seq(mstart,mstop,minc),c1)

lines(c1,lmm1$coefficients[2]*c1+lmm1$coefficients[1],lty=3)
lines(c1,lmm2$coefficients[2]*c1+lmm2$coefficients[1],col="red",lty=3)
lines(c1,lmm3$coefficients[2]*c1+lmm3$coefficients[1],col="blue",lty=3)
#max(ldg.sim.a[,3])+1
dev.off()
################################################################
#Figure ??: 
#Make a fancy multi-figure image showing the environmental novelty, 
#physical distance, and contribution of variables to the novelty.
################################################################
env.ind = c(1,4,7,9,10)
#env_distance = get_env_distance(xx,xx_new,env.ind)

fig.name = paste("calanda_novel_igr1_5C.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

s=4
env.name=paste(prefix.list[[g]],variable.list[[17]] ,sep="")
assign("env_analogue",eval(as.name(env.name)) )
env_analogue = env_analogue [[s]]
var_dist = env_analogue[[4]]
colnames(var_dist)=colnames(xx)[env.ind]
an_dist = sqrt( (env_analogue[[3]] - xx_new$elevation)^2)

layout.matrix=matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)
layout(mat = layout.matrix,
       heights = c(5,1,1), # Heights of the two rows
       widths = c(10,1)) # Widths of the two columns

#layout.show(6)

nlevel=64
par(oma = c(3,3,3,3))
par( mar = c(0.5,4,0,4) )
image(x=xx$elevation, y = 1:length(env.ind), var_dist, ylab="", yaxt='n',xaxt='n',col=viridis(nlevel))
axis(2, at=1:(length(env.ind)), labels = colnames(var_dist),cex.axis=0.9)

par( mar = c(0.5,4,0,4) )
#image(x=xx$elevation, y = 1, env_analogue[[2]],ylab="", yaxt='n',xaxt='n',col=viridis(nlevel))
plot(x=xx$elevation,  env_analogue[[2]],ylab="total",xaxt='n',t="l")
#axis(2, at=1, labels = c("total"),cex.axis=0.9)

#This is how far a site is in actual space from its current location
par( mar = c(0.5,4,0,4) )
#image(x=xx$elevation, y = 1, an_dist, xlab="", ylab="",yaxt='n',col=viridis(nlevel))
plot(x=xx$elevation,an_dist,ylab="distance",t="l")
#lines(xx$elevation,xx$elevation,lty=3)
abline(h=0,lty=3)
#axis(2, at=1, labels = c("distance"),cex.axis=0.9)

#Color bar 1
par( mar = c(4,1,4,1) )
image(1,(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), t(seq(min(var_dist),max(var_dist),max(var_dist)/nlevel)), ylab="",xaxt='n',col=viridis(nlevel))

#Color bar 2
# par( mar = c(1,1,1,1) )
# image(1, (seq(min(env_analogue[[2]]),max(env_analogue[[2]]),max(env_analogue[[2]])/nlevel)), t(seq(min(env_analogue[[2]]),max(env_analogue[[2]]),max(env_analogue[[2]])/nlevel)), ylab="",xaxt='n',col=viridis(nlevel))

#Color bar 3
# par( mar = c(1,1,1,1) )
# image(1, (seq(min(an_dist),max(an_dist),max(an_dist)/nlevel)), t(seq(min(an_dist),max(an_dist),max(an_dist)/nlevel)), ylab="",xaxt='n',col=viridis(nlevel))
 
dev.off()


################################################################
#Figure ??: 
#The low-density equilibrium distriubtions of species
#################################################################

epoints = match( elevations,xx$elevation)

# fig.name = paste("ldeq_multi_igrsites3.pdf",sep="")
# pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

#par(mfrow=c(1,1),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
ylim=c(-1,70)
#Current: 
plot(xx$elevation, sc1Frs.all[[1]][1,,1], t="l", ylab="Seeds (per-capita, Kriged) ", xlab="Elevation", xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim)
lines(xx$elevation,sc1Frs.all[[1]][1,,1])
lines(xx$elevation, sc1Frs.all[[1]][1,,2],col="red")
lines(xx$elevation, sc1Frs.all[[1]][1,,3],col="blue")
abline(v=c(elevations),lwd=10,col="grey80")
points(xx$elevation[epoints], sc1Frs.all[[1]][1,epoints,1])
points(xx$elevation[epoints], sc1Frs.all[[1]][1,epoints,2],col="red")
points(xx$elevation[epoints], sc1Frs.all[[1]][1,epoints,3],col="blue")

#Future
x.txt = c(1900,1200,2700)
for(n in 2:ngenst){
lines(xx$elevation, sc1Frs.all[[n]][1,,1],lty=2)
text(x=x.txt[1], y=sc1Frs.all[[n]][1,x.txt[1],1] , labels=m.labels[(n)])
lines(xx$elevation, sc1Frs.all[[n]][1,,2],col="red",lty=2)
text(x=x.txt[2], y=sc1Frs.all[[n]][1,x.txt[2],2] , labels=m.labels[(n)],col="red")
lines(xx$elevation, sc1Frs.all[[n]][1,,3],col="blue",lty=2)
text(x=x.txt[3], y=sc1Frs.all[[n]][1,x.txt[3],3] , labels=m.labels[(n)],col="blue")
points(xx$elevation[epoints], sc1Frs.all[[n]][1,epoints,1])
points(xx$elevation[epoints], sc1Frs.all[[n]][1,epoints,2],col="red")
points(xx$elevation[epoints], sc1Frs.all[[n]][1,epoints,3],col="blue")
}
# f1 = flower_prob_act*flower_act
# points(xx$elevation[epoints], f1[,1],pch=3)
# points(xx$elevation[epoints], f1[,2],col="red",pch=3)
# points(xx$elevation[epoints], f1[,3],col="blue",pch=3)

dev.off()
