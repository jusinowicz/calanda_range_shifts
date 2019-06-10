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
library(cowplot)

#==============================================================================
# Load Data
#==============================================================================

#This requires user input!

variable.list=list("l1", "D", "var_mu_Us","cov_e_mu_Us",
"cov_lam_vc", "cov_lam_vc2", "Elam1", "Elam2", "gr1.n", "gr1", "y.full", 
"w.eq.full","ldg.sim", "ldgIMP.sim", "covIMP.sim", "Frs", "env_analogue_all",
"sim.impacts", "sim.impactsIMP","covIMP.impacts","lded","nf","is.final","Frs_ci")



#List of scenario file names 

# file.name.list=c("calanda_ccs26_temp_2017_waldmean.var", 
#   "calanda_ccs26_tsm_2017_waldmean.var", 
#   "calanda_ccs26B_tsm_2017_waldmean.var", 
#   "calanda_ccs85_tsm_2017_waldmean.var", 
#   "calanda_ccs85B_tsm_2017_waldmean.var")

# file.name.list=c("calanda_ccs26_temp_2017_waldmean.var", 
#   "calanda_ccs26_tsm_2017_waldmean.var", 
#   "calanda_ccs26B_tsm_2017_waldmean.var", 
#   "calanda_ccs85_tsmfull_2017_waldmean.var", 
#   "calanda_ccs85B_tsmfull_2017_waldmean.var")

file.name.list=c(
  "calanda_ccs85_tsmfull3_2017_waldmean.var", 
  "calanda_ccs85B_tsmfull3_2017_waldmean.var",  
  "calanda_ccs26_tsmfull3_2017_waldmean.var", 
  "calanda_ccs26B_tsmfull3_2017_waldmean.var")

# file.name.list=c("calanda_ccs26_temp_2017_waldmean.var", 
#   "calanda_ccs85_temp_2017_waldmean.var",
#   "calanda_ccs26_tsm_2017_waldmean.var", 
#   "calanda_ccs85_tsm_2017_waldmean.var",
#   "calanda_ccs26B_tsm_2017_waldmean.var", 
#   "calanda_ccs85B_tsm_2017_waldmean.var")


#In order to deal with each of the scenarios simultaneously, it is necessary
#to systematically rename va"calanda_repeat1_A_multi_2017","calanda_repeat1_B_multi_2017",riables based on their scenario. 
#Assosciated prefixes for renaming

#prefix.list= list("sc1")
#prefix.list= list("sc26temp","sc26tsm","sc26Btsm","sc85tsm","sc85Btsm")
prefix.list= list("sc85tsm","sc85Btsm","sc26tsm","sc26Btsm")

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

#The environmental novelty files:
load("env_all_ccs85full.var")
env_analogue_all85 = env_analogue_all
load("env_all_ccs26full.var")
env_analogue_all26 = env_analogue_all

#==============================================================================
# Figure: Intrinsic and realized ranges, before and after climate change
# panel with environmental novelty. 
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
# mstart = 0
# mstop = 6
# minc = 0.5
# mtot= ceiling((mstop-mstart)/minc)
# m.labels=as.character(seq(mstart,mstop,minc))

#Time increments
mstart = 1
mstop = 61
minc = 20
mtot= ceiling((mstop-mstart)/minc)
m.labels=as.character(seq(mstart,mstop,minc)-1)

#=============================================================================

epoints = match( elevations,xx$elevation)

# fig.name = paste("calanda_ranges_ipcc3.pdf",sep="")
# pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)
fig.name = paste("calanda_ranges_tsmipcc3.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

layout.matrix=matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)
layout(mat = layout.matrix,
       heights = c(5,5,1.5), # Heights of the two rows
       widths = c(10,10)) # Widths of the two columns

#layout.show(6)

#par(mfrow=c(2,1),mai= c( 0.0, 0.2, 0.0, 0.2), omi=c(0.5,0.75,0.5,0.75)) #,mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))

###Common figure properties
g=c(1,3) #Which data to use? E.g. Best-fit intrinsic ranges. 
nlevel = 64 #For viridis color scheme
ylimR1 =c(-1,70) # Common y limits for row 1 
ylimR2 =c(-1,280) # Common y limits for row 2  
xlimC1 =c(0,3500) # Common x limits for Column 1 
xlimC2 =c(0,3500) # Common x limits for Column 2
x.txt = c(1600,2000,1200) #X location of labels on ranges
par(oma = c(3,2,3,3) )

# sc1env_analogue_all = env_analogue_all85
# env_analogue_a2 = matrix(0,3501,7)
# for(a in 1:7){ env_analogue_a2[,a]=matrix(unlist(sc1env_analogue_all[[a]][[2]]),3501,1) }

# sc1env_analogue_all = env_analogue_all26
# env_analogue_a2 = matrix(0,3501,7)
# for(a in 1:7){ env_analogue_a2[,a]=matrix(unlist(sc1env_analogue_all[[a]][[2]]),3501,1) }


#####Column 1
###   Panel 1: Intrinsic range, IPCC 8.5

new.name=paste(prefix.list[[g[1]]],"Frs",sep="")  
assign("sc1Frs",eval(as.name(new.name)))
new.name2=paste(prefix.list[[g[1]]],"nf",sep="")  
assign("sc1nf",eval(as.name(new.name2)))


#par(mfrow=c(1,1),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))

#Current: 
par( mar = c(0.5,4,0,0.5) )
plot(xx$elevation, sc1Frs[1,,1], t="l", xlab="", ylab="Seeds (per-capita, Kriged) ",  xaxt="n", xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylimR1)
#plot(xx$elevation, sc1Frs[1,,1], t="l", xlab="Elevation", ylab="Seeds (per-capita, Kriged) ",   xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim)

abline(v=c(elevations),lwd=10,col="grey80")
lines(xx$elevation,sc1Frs[1,,1])
lines(xx$elevation, sc1Frs[1,,2],col="red")
lines(xx$elevation, sc1Frs[1,,3],col="blue")
text(x=x.txt[1], y=sc1Frs[1,x.txt[1],1], labels=m.labels[1])
text(x=x.txt[2], y=sc1Frs[1,x.txt[2],2] , labels=m.labels[1],col="red")
text(x=x.txt[3], y=sc1Frs[1,x.txt[3],3] , labels=m.labels[1],col="blue")

# polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[1]][[1]][1,]*6*8*1.5), sc1Frs_ci[[1]][[1]][2,]*6*8*1.5),col=alpha(rgb(t(matrix(col2rgb("grey20"))),maxColorValue=255),0.2), border = NA)
# polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[1]][[2]][1,]), sc1Frs_ci[[1]][[2]][2,]),col=alpha(rgb(t(matrix(col2rgb("red"))),maxColorValue=255),0.2), border = NA)
# polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[1]][[3]][1,]), sc1Frs_ci[[1]][[3]][2,]),col=alpha(rgb(t(matrix(col2rgb("blue"))),maxColorValue=255),0.2), border = NA)


# points(xx$elevation[epoints], sc1Frs.all[[1]][1,epoints,1])
# points(xx$elevation[epoints], sc1Frs.all[[1]][1,epoints,2],col="red")
# points(xx$elevation[epoints], sc1Frs.all[[1]][1,epoints,3],col="blue")

#Future

#for(n in seq(3,length(m.labels),2)){
for(n in 1:length(m.labels)) {
  lines(xx$elevation, sc1Frs[n,,1],lty=2)
  text(x=x.txt[1], y=sc1Frs[n,x.txt[1],1] , labels=m.labels[(n)])
  lines(xx$elevation, sc1Frs[n,,2],col="red",lty=2)
  text(x=x.txt[2], y=sc1Frs[n,x.txt[2],2] , labels=m.labels[(n)],col="red")
  lines(xx$elevation, sc1Frs[n,,3],col="blue",lty=2)
  text(x=x.txt[3], y=sc1Frs[n,x.txt[3],3] , labels=m.labels[(n)],col="blue")

  # polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[n]][[1]][1,]*6*8*1.5), sc1Frs_ci[[n]][[1]][2,]*6*8*1.5),col=alpha(rgb(t(matrix(col2rgb("grey80"))),maxColorValue=255),0.4), border = NA)
  # polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[n]][[2]][1,]), sc1Frs_ci[[n]][[2]][2,]),col=alpha(rgb(t(matrix(col2rgb("red"))),maxColorValue=255),0.4), border = NA)
  # polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[n]][[3]][1,]), sc1Frs_ci[[n]][[3]][2,]),col=alpha(rgb(t(matrix(col2rgb("blue"))),maxColorValue=255),0.4), border = NA)

  # points(xx$elevation[epoints], sc1Frs.all[[n]][1,epoints,1])
  # points(xx$elevation[epoints], sc1Frs.all[[n]][1,epoints,2],col="red")
  # points(xx$elevation[epoints], sc1Frs.all[[n]][1,epoints,3],col="blue")
}
# f1 = flower_prob_act*flower_act
# points(xx$elevation[epoints], f1[,1],pch=3)
# points(xx$elevation[epoints], f1[,2],col="red",pch=3)
# points(xx$elevation[epoints], f1[,3],col="blue",pch=3)


###  Panel 2: Realized ranges
a1=sc1nf[,,1]/sc1Frs[,,1]
a1[is.infinite(a1)]=0
a1 = apply(a1 ,1,max)
c1 = a1[2]/a1[1]

a2=sc1nf[,,2]/sc1Frs[,,2]
a2[is.infinite(a2)]=0
a2 = apply(a2 ,1,max)
c2 = a2[2]/a2[1]

a3=sc1nf[,,3]/sc1Frs[,,3]
a3[is.infinite(a3)]=0
a3 = apply(a3 ,1,max)
c3 = a3[6]/a3[1]

c1=1; c2=1; 
c3=1


#Current: 
par( mar = c(0.5,4,0,0.5) )
plot(xx$elevation, sc1nf[1,,1]*c1, t="l", ylab="Population density ", xlab="", xaxs="i", xaxt="n",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylimR2)
abline(v=c(elevations),lwd=10,col="grey80")
lines(xx$elevation,sc1nf[1,,1]*c1)
lines(xx$elevation, sc1nf[1,,2]*c2,col="red")
lines(xx$elevation, sc1nf[1,,3]*c3,col="blue")
text(x=x.txt[1], y=sc1nf[1,x.txt[1],1]*c1 , labels=m.labels[1])
text(x=x.txt[2], y=sc1nf[1,x.txt[2],2]*c2 , labels=m.labels[1],col="red")
text(x=x.txt[3], y=sc1nf[1,x.txt[3],3] *c3, labels=m.labels[1],col="blue")

# points(xx$elevation[epoints], sc1nf.all[[1]][1,epoints,1])
# points(xx$elevation[epoints], sc1nf.all[[1]][1,epoints,2],col="red")
# points(xx$elevation[epoints], sc1nf.all[[1]][1,epoints,3],col="blue")

#Future

#for(n in seq(3,length(m.labels),2)){
for(n in 2:length(m.labels)) {
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


###Environmental novelty
# new.name=paste(prefix.list[[g[1]]],"env_analogue_all",sep="")  
# assign("sc1env_analogue_all",eval(as.name(new.name)))
sc1env_analogue_all = env_analogue_all85
env_analogue_a2 = matrix(0,3501,7)
for(a in 1:7){ env_analogue_a2[,a]=matrix(unlist(sc1env_analogue_all[[a]][[2]]),3501,1) }
# env_analogue_a3 = matrix(0,3501,7)
# for(a in 1:7){ env_analogue_a3[,a]=matrix(unlist(sc1env_analogue_all[[a]][[3]]),3501,1) }

par( mar = c(2,4,0,0.5) )
image(xx$elevation, 2:7,env_analogue_a2[,2:7], ylab="", yaxt='n',xaxt='n',xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3, col=viridis(nlevel))
abline(v=c(elevations),lwd=10,col=alpha(rgb(t(matrix(col2rgb("grey80"))),maxColorValue=255),0.3))

par(new=TRUE)
contour(xx$elevation,2:7, env_analogue_a2[,2:7],ylab="",yaxt='n',xaxs="i",yaxs="i",cex.main=1.3, cex.lab=1.3)
axis(2, at=seq(3,7,2), labels = as.character(seq(1,5,2)),cex.axis=0.9)
par(new=FALSE)

#This assumes that env_analogue_all85 has the bigger range
zlimR3 = range(env_analogue_a2)


#####Column 2
###   Panel 1: Intrinsic range, IPCC 2.6

new.name=paste(prefix.list[[g[2]]],"Frs",sep="")  
assign("sc1Frs",eval(as.name(new.name)))
new.name2=paste(prefix.list[[g[2]]],"nf",sep="")  
assign("sc1nf",eval(as.name(new.name2)))


#Current: 
par( mar = c(0.5,0,0,4) )
plot(xx$elevation, sc1Frs[1,,1], t="l", xlab="", ylab=" ",  xaxt="n", yaxt='n',
  xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylimR1)
#plot(xx$elevation, sc1Frs[1,,1], t="l", xlab="Elevation", ylab="Seeds (per-capita, Kriged) ",   xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim)

abline(v=c(elevations),lwd=10,col="grey80")
lines(xx$elevation,sc1Frs[1,,1])
lines(xx$elevation, sc1Frs[1,,2],col="red")
lines(xx$elevation, sc1Frs[1,,3],col="blue")
text(x=x.txt[1], y=sc1Frs[1,x.txt[1],1], labels=m.labels[1])
text(x=x.txt[2], y=sc1Frs[1,x.txt[2],2] , labels=m.labels[1],col="red")
text(x=x.txt[3], y=sc1Frs[1,x.txt[3],3] , labels=m.labels[1],col="blue")

# polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[1]][[1]][1,]*6*8*1.5), sc1Frs_ci[[1]][[1]][2,]*6*8*1.5),col=alpha(rgb(t(matrix(col2rgb("grey20"))),maxColorValue=255),0.2), border = NA)
# polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[1]][[2]][1,]), sc1Frs_ci[[1]][[2]][2,]),col=alpha(rgb(t(matrix(col2rgb("red"))),maxColorValue=255),0.2), border = NA)
# polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[1]][[3]][1,]), sc1Frs_ci[[1]][[3]][2,]),col=alpha(rgb(t(matrix(col2rgb("blue"))),maxColorValue=255),0.2), border = NA)


# points(xx$elevation[epoints], sc1Frs.all[[1]][1,epoints,1])
# points(xx$elevation[epoints], sc1Frs.all[[1]][1,epoints,2],col="red")
# points(xx$elevation[epoints], sc1Frs.all[[1]][1,epoints,3],col="blue")

#Future

#for(n in seq(3,length(m.labels),2)){
for(n in 1:length(m.labels)) {
  lines(xx$elevation, sc1Frs[n,,1],lty=2)
  text(x=x.txt[1], y=sc1Frs[n,x.txt[1],1] , labels=m.labels[(n)])
  lines(xx$elevation, sc1Frs[n,,2],col="red",lty=2)
  text(x=x.txt[2], y=sc1Frs[n,x.txt[2],2] , labels=m.labels[(n)],col="red")
  lines(xx$elevation, sc1Frs[n,,3],col="blue",lty=2)
  text(x=x.txt[3], y=sc1Frs[n,x.txt[3],3] , labels=m.labels[(n)],col="blue")

  # polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[n]][[1]][1,]*6*8*1.5), sc1Frs_ci[[n]][[1]][2,]*6*8*1.5),col=alpha(rgb(t(matrix(col2rgb("grey80"))),maxColorValue=255),0.4), border = NA)
  # polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[n]][[2]][1,]), sc1Frs_ci[[n]][[2]][2,]),col=alpha(rgb(t(matrix(col2rgb("red"))),maxColorValue=255),0.4), border = NA)
  # polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[n]][[3]][1,]), sc1Frs_ci[[n]][[3]][2,]),col=alpha(rgb(t(matrix(col2rgb("blue"))),maxColorValue=255),0.4), border = NA)

  # points(xx$elevation[epoints], sc1Frs.all[[n]][1,epoints,1])
  # points(xx$elevation[epoints], sc1Frs.all[[n]][1,epoints,2],col="red")
  # points(xx$elevation[epoints], sc1Frs.all[[n]][1,epoints,3],col="blue")
}
# f1 = flower_prob_act*flower_act
# points(xx$elevation[epoints], f1[,1],pch=3)
# points(xx$elevation[epoints], f1[,2],col="red",pch=3)
# points(xx$elevation[epoints], f1[,3],col="blue",pch=3)


###  Panel 2: Realized ranges
a1=sc1nf[,,1]/sc1Frs[,,1]
a1[is.infinite(a1)]=0
a1 = apply(a1 ,1,max)
c1 = a1[7]/a1[1]

a2=sc1nf[,,2]/sc1Frs[,,2]
a2[is.infinite(a2)]=0
a2 = apply(a2 ,1,max)
c2 = a2[7]/a2[1]

a3=sc1nf[,,3]/sc1Frs[,,3]
a3[is.infinite(a3)]=0
a3 = apply(a3 ,1,max)
c3 = a3[7]/a3[1]

#c1=1; 
c2=1; 
c3=1

#Current: 
par( mar = c(0.5,0,0,4) )
plot(xx$elevation, sc1nf[1,,1]*c1, t="l", ylab="", xlab="", xaxs="i", xaxt="n", yaxt='n',
  yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylimR2)
abline(v=c(elevations),lwd=10,col="grey80")
lines(xx$elevation,sc1nf[1,,1]*c1)
lines(xx$elevation, sc1nf[1,,2]*c2,col="red")
lines(xx$elevation, sc1nf[1,,3],col="blue")
text(x=x.txt[1], y=sc1nf[1,x.txt[1],1]*c1 , labels=m.labels[1])
text(x=x.txt[2], y=sc1nf[1,x.txt[2],2]*c2 , labels=m.labels[1],col="red")
text(x=x.txt[3], y=sc1nf[1,x.txt[3],3] , labels=m.labels[1],col="blue")

# points(xx$elevation[epoints], sc1nf.all[[1]][1,epoints,1])
# points(xx$elevation[epoints], sc1nf.all[[1]][1,epoints,2],col="red")
# points(xx$elevation[epoints], sc1nf.all[[1]][1,epoints,3],col="blue")

#Future

#for(n in seq(3,length(m.labels),2)){
for(n in 2:length(m.labels)) {
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


###Environmental novelty
# new.name=paste(prefix.list[[g[1]]],"env_analogue_all",sep="")  
# assign("sc1env_analogue_all",eval(as.name(new.name)))

sc1env_analogue_all = env_analogue_all26
env_analogue_a2 = matrix(0,3501,7)
for(a in 1:7){ env_analogue_a2[,a]=matrix(unlist(sc1env_analogue_all[[a]][[2]]),3501,1) }
# env_analogue_a3 = matrix(0,3501,7)
# for(a in 1:7){ env_analogue_a3[,a]=matrix(unlist(sc1env_analogue_all[[a]][[3]]),3501,1) }


par( mar = c(2,0,0,4) )
image(xx$elevation, 2:7,env_analogue_a2[,2:7], zlim =zlimR3,  ylab="", yaxt='n',xaxt='n',xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3, col=viridis(nlevel))
abline(v=c(elevations),lwd=10,col=alpha(rgb(t(matrix(col2rgb("grey80"))),maxColorValue=255),0.3))

par(new=TRUE)
contour(xx$elevation,2:7, env_analogue_a2[,2:7],ylab="",yaxt='n',xaxs="i",yaxs="i",cex.main=1.3, cex.lab=1.3)
#axis(2, at=seq(3,7,2), labels = as.character(seq(1,5,2)),cex.axis=0.9)
par(new=FALSE)

dev.off()


#=============================================================================
#Figure 3B: Non-Arch (linear) figures showing how the igr of each 
# scenario moves according to the benefit of either range width 
# or range overlap.
# This version for the full 3spp picture. 
# Components STANDARDIZED BY MEAN coefficients
# Figures are a SINGLE PLOT
#=============================================================================

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


fig.name = paste("calanda_widthOverlap_full85B.pdf",sep="")
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
for( g in 2:2){ 


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

  for(s in 2:7) {
    
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


#for( s in seq(3,length(m.labels),2)){ 
for( s in 2:7){ 


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

#==============================================================================
# Figure: Compare the change in LDGs with climate change  across scenarios 
# and fits. 
#==============================================================================
#==============================================================================
# Proportionate change. 
#==============================================================================

fig.name = paste("calanda_igrs_bothmodels1.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Times', pointsize=16)

color.use=list("black","red","blue")
nspp=3
ylims=c(0.5,1.5)
pch85B = 16
pch26 = 2
pch26B = 17

plot(sc85tsmldg.sim[,1]/sc85tsmldg.sim[1,1], ylim=ylims, ylab="Decade", 
  xlab="Relative change in LGR", xaxs="i", cex.main=1.3,cex.lab=1.3)

#sc85
for( n in 1:nspp){
  points( sc85tsmldg.sim[,n]/sc85tsmldg.sim[1,n], col =color.use[[n]])
  lines(  sc85tsmldg.sim[,n]/sc85tsmldg.sim[1,n], col =color.use[[n]],lty = 3)
}
#sc85 B
for( n in 1:nspp){
  points( sc85Btsmldg.sim[,n]/sc85Btsmldg.sim[1,n], col =color.use[[n]],pch=pch85B)
  lines(  sc85Btsmldg.sim[,n]/sc85Btsmldg.sim[1,n], col =color.use[[n]],lty = 1)

}
#sc26
for( n in 1:nspp){
  points( sc26tsmldg.sim[,n]/sc26tsmldg.sim[1,n], col =color.use[[n]],pch=pch26)
  lines(  sc26tsmldg.sim[,n]/sc26tsmldg.sim[1,n], col =color.use[[n]],lty = 3)

}
#sc26B
for( n in 1:nspp){
  points( sc26Btsmldg.sim[,n]/sc26Btsmldg.sim[1,n], col =color.use[[n]],pch=pch26B)
  lines(  sc26Btsmldg.sim[,n]/sc26Btsmldg.sim[1,n], col =color.use[[n]],lty = 1)

}

dev.off()
#==============================================================================
# Absolute change. 
#==============================================================================

ylims=c(1,4)
#ylims=c(1,1.3)
plot(sc85tsmldg.sim[,1], ylim=ylims)
points(sc85Btsmldg.sim[,1], pch=16)
points(sc26tsmldg.sim[,1], pch=2)
points(sc26Btsmldg.sim[,1], pch=17)
points(sc85tsmldg.sim[,2], col="red")
points(sc85Btsmldg.sim[,2], col="red", pch=16)
points(sc26tsmldg.sim[,2], col="red", pch=2)
points(sc26Btsmldg.sim[,2], pch=17, col="red")
points(sc85tsmldg.sim[,3], col="blue")
points(sc85Btsmldg.sim[,3], pch=16, col="blue")
points(sc26tsmldg.sim[,3], pch=2, col="blue")
points(sc26Btsmldg.sim[,3], pch=17, col="blue")

#==============================================================================
# Figure: Population spread rates under an annual change in climate. 
# Load different data files for this figure: 
#==============================================================================
file.name.list=c(
  "calanda_ccs85B_tsmfull3_fast_2017_waldmean.var",  
  "calanda_ccs26B_tsmfull3_fast_2017_waldmean.var")
prefix.list= list("sc1","sc2")

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

#Climate change scenarios
ccs_temp1 = read.csv("rcp85_2019_2100.csv")[2:3]
ccs_temp2 = read.csv("rcp26_2019_2100.csv")[2:3]

#Make these relative increases
temp_field = mean(allB2$gs_mean_temp)
ccs_temp1 = ccs_temp1[,2]- temp_field
ccs_temp1=matrix(c(0,ccs_temp1))
ccs_temp2 = ccs_temp2[,2]- temp_field
ccs_temp2=matrix(c(0,ccs_temp2))

#Calculate the spread rates with a_rr for exponential kernel
spread_rates1 = sc1ldg.sim[(1:ngenst), ]
spread_rates2 = sc2ldg.sim[(1:ngenst), ]

for (q in 1:ngenst){ 
  for (q2 in 1:nspp){ 
  spread_rates1[q,q2] = get.spread.rate(sc1ldg.sim[q,q2],a_rr[q2],sr[q2])
  spread_rates2[q,q2] = get.spread.rate(sc2ldg.sim[q,q2],a_rr[q2],sr[q2])
  }
}

#==============================================================================
# Figure: Change in LGR etc for IPCC scenarios
#==============================================================================
color.use = c( "black","red","blue" )
nspp=3 #Number of species
ngens= 40
iconfig=1 #NUmber of interattions of initial configuration of species
ngenst=iconfig+ngens
nlevel=64

#fig.name = paste("calanda_ldg_ipcc2.pdf",sep="")
#pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

layout.matrix=matrix(c(1,2,5,3,4,6), nrow = 3, ncol = 2)
layout(mat = layout.matrix,
       heights = c(1,5,3), # Heights of the rows
       widths = c(5,5)) # Widths of the columns

#layout.show(6)

#Panel 1: Image of temp increase for scenario 1
par( mar = c(0,5,4,0.2) )
years1 = as.character(2019:(2019+ngens))
xlim = c(2019,(2019+ngens))
ylim=c(0.8,5.5)
zlim = range(c(ccs_temp1[1:ngens],ccs_temp2[1:ngens]))
#Color bar 1 for temp change
maxt1 = floor(max(ccs_temp1[1:ngens]))
mint1 = 0 
mincs1 = maxt1-mint1
image(as.matrix(1:ngens), 1, (as.matrix(ccs_temp1[1:ngens])), zlim=zlim, ylab="",xaxt='n',yaxt='n', col=viridis(nlevel))
#axis(3, at = as.matrix(1:ngens)[which(round(ccs_temp1[1:ngens],1) == seq(1,maxt1,1))], labels = round(ccs_temp1[which(round(ccs_temp1[1:ngens],1) == seq(1,maxt1,1))],1))
axis(3, at =seq(10,ngens,((xlim[2]-xlim[1])/4)), labels = round(ccs_temp1[1:ngenst][(2019:(2019+ngens)) %in% seq(xlim[1]+10,xlim[2],((xlim[2]-xlim[1])/4))],1))



#Panel 2: LDG for scenario 1
#par( mar = c(5,5,0,0.1) )
par( mar = c(0,5,0,0.2) )
plot(years1, sc1ldg.sim[(1:ngenst),1], t="l", xlab="", ylab="Low-density growth rate",  xaxt="n", yaxt="n",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim,col=color.use[1])
axis(1, at=seq(xlim[1]+10,xlim[2],((xlim[2]-xlim[1])/4)), labels = FALSE) #as.character(seq(xlim[1]+11,xlim[2]+1,((xlim[2]-xlim[1])/4))),cex.axis=1)
axis(2, at=seq(ylim[1],ylim[2],((ylim[2]-ylim[1])/5)), labels = as.character(round(seq(ylim[1],ylim[2],((ylim[2]-ylim[1])/5)),2)),cex.axis=1)

for(n in 2:3){
  lines(years1, sc1ldg.sim[(1:ngenst),n],col=color.use[n])
}
abline(h=1, col="grey50", lty=3)


#Panel 3: Image of temp increase for scenario 2
par( mar = c(0,0.2,4,5.5) )
#Color bar 2 for temp change
maxt2 = floor(max(ccs_temp2[1:ngens]))
mint2 = 0 
mincs2 = maxt2-mint2
image(as.matrix(1:ngens), 1, (as.matrix(ccs_temp2[1:ngens])), zlim=zlim,  ylab="",xaxt='n',yaxt='n', col=viridis(nlevel))
#axis(3, at = as.matrix(1:ngens)[which(round(ccs_temp2[1:ngens],1) == seq(1,maxt2,1))], labels = round(ccs_temp2[which(round(ccs_temp2[1:ngens],1) == seq(1,maxt2,1))],1))
axis(3, at =seq(10,ngens,((xlim[2]-xlim[1])/4)), labels = round(ccs_temp2[1:ngenst][(2019:(2019+ngens)) %in% seq(xlim[1]+10,xlim[2],((xlim[2]-xlim[1])/4))],1))


#Panel 4: LDG for scenario 2
par( mar = c(0,0.2,0,5.5) )
years1 = as.character(2019:(2019+ngens))
xlim = c(2019,(2019+ngens))
ylim=c(0.8,5.5)
plot(years1, sc2ldg.sim[(1:ngenst),1], t="l", xlab="", ylab="",  xaxt="n", yaxt="n",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim,col=color.use[1])
axis(1, at=seq(xlim[1]+10,xlim[2],((xlim[2]-xlim[1])/4)), labels = FALSE) #as.character(seq(xlim[1]+11,xlim[2]+1,((xlim[2]-xlim[1])/4))),cex.axis=1)
axis(2, at=seq(ylim[1],ylim[2],((ylim[2]-ylim[1])/5)), labels = FALSE) #as.character(seq(ylim[1],ylim[2],((ylim[2]-ylim[1])/5))),cex.axis=1)

for(n in 2:3){
  lines(years1, sc2ldg.sim[(1:ngenst),n],col=color.use[n])
}
abline(h=1, col="grey50", lty=3)



#Panel 5: Spread rates for scenario 1
par( mar = c(5,5,0.2,0.2) )
years1 = as.character(2019:(2019+ngens))
xlim = c(2019,(2019+ngens))
ylim=c(0.5,2)
plot(years1, spread_rates1[(1:ngenst),1], t="l", ylab="Population spread rate(m)",  xaxt="n", yaxt="n",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim,col=color.use[1])
axis(1, at=seq(xlim[1]+10,xlim[2],((xlim[2]-xlim[1])/4)), labels = as.character(seq(xlim[1]+11,xlim[2]+1,((xlim[2]-xlim[1])/4))),cex.axis=1)
axis(2, at=seq(ylim[1],ylim[2],((ylim[2]-ylim[1])/2)), labels = as.character(seq(ylim[1],ylim[2],((ylim[2]-ylim[1])/2))),cex.axis=1)

for(n in 2:3){
  lines(years1, spread_rates1[(1:ngenst),n],col=color.use[n])
}

#Panel 6: Spread rates for scenario 2
par( mar = c(5,0.2,0.2,5.5) )
years1 = as.character(2019:(2019+ngens))
xlim = c(2019,(2019+ngens))
ylim=c(0.5,2)
plot(years1, spread_rates2[(1:ngenst),1], t="l", ylab="",  xaxt="n", yaxt="n",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim,col=color.use[1])
axis(1, at=seq(xlim[1]+10,xlim[2],((xlim[2]-xlim[1])/4)), labels = as.character(seq(xlim[1]+11,xlim[2]+1,((xlim[2]-xlim[1])/4))),cex.axis=1)
axis(2, at=seq(ylim[1],ylim[2],((ylim[2]-ylim[1])/2)), labels = FALSE) #as.character(seq(ylim[1],ylim[2],((ylim[2]-ylim[1])/5))),cex.axis=1)

for(n in 2:3){
  lines(years1, spread_rates2[(1:ngenst),n],col=color.use[n])
}

dev.off()


#==============================================================================
#Visualizing environmental space
#==============================================================================

#==============================================================================
# Figure: Visulaize environmental space as a difference between bivariate bins 
# with hexbin
#==============================================================================
#Create two data frames of the data: 

df1  = data.frame( gs_mean_temp = allB2$gs_mean_temp, soil_moist = allB2$soil_moist)
df2  = data.frame( gs_mean_temp = xx$gs_mean_temp, soil_moist = xx$soil_moist)

#Make identical bins by constructing one hexbin object.  
#Use dplyr's bind_rows to keep a track of which data.frame the data came from 

bothDF = bind_rows(A = df1, B = df2, .id = "df")
bothHex = hexbin(bothDF$gs_mean_temp, bothDF$soil_moist, IDs = TRUE)

#Count the occurrences of each within each cell. 
#First, apply across the bins, constructing a table (needs to use factor to make sure all levels are shown; 
#not needed if your column is already a factor). 
#Then, simplify and construct a data.frame that is 
#then manipluated with mutate to calculate the difference in counts and 
#then joined back to a table that gives 
#the x and y values for each of the id's.

counts =
  hexTapply(bothHex, factor(bothDF$df), table) %>%
  simplify2array %>%
  t %>%
  data.frame() %>%
  mutate(id = as.numeric(row.names(.))
         , diff = A - B) %>%
  left_join(data.frame(id = bothHex@cell, hcell2xy(bothHex)))

counts %>%
  ggplot(aes(x = x, y = y
             , fill = diff)) +
  geom_hex(stat = "identity") +
  coord_equal() +
  scale_fill_gradient2()

ggplot(aes(df1$, y = y
             , fill = diff))


#==============================================================================
# Figure: Visulaize environmental space as two separate 2D bins
#==============================================================================
p1 = ggplot(xx, aes(x=gs_mean_temp, y=soil_moist) ) + geom_bin2d(bins=10)+ scale_fill_gradient(low="red", high="red", limits=c(0,500) )
p2 = ggplot(allB2, aes(x=gs_mean_temp, y=soil_moist) ) + geom_bin2d(bins=10)+ scale_fill_gradient(low="white", high="blue", limits=c(0,500) )

# combined plot using layer_data() to extract data
# from each plot above, & scale_identity to use the
# already calculated fill values
p.combined = ggplot(data = layer_data(p1), aes(x = x, y = y, fill = fill)) +
  geom_tile() +
  geom_tile(data = layer_data(p2)) +
  scale_fill_identity()

# optional: add legends to the combined plot
plot_grid(p.combined,
          plot_grid(get_legend(p2),
                    get_legend(p1),
                    ncol = 1),
          nrow = 1,
          rel_widths = c(1, 0.1)) # optimal relative widths may change, depending on plot dimensions

#==============================================================================
# Figure: Visulaize environmental space as two separate 2D bins. Layer the bins 
# for successive projections of the climate space
#==============================================================================

#Time increments
mstart = 1
mstop = 41
minc = 10
mtot= ceiling((mstop-mstart)/minc)
m.labels=as.character(seq(mstart,mstop,minc))

#Rename the files from each scenario
var.length=mtot
nscen = length(file.name.list)

#The combined call to ggplot for later: 
p_call = "ggplot(data = layer_data(p1), aes(x = x, y = y, fill = fill))+geom_tile() + "

for( f in 1:(var.length)) {      
	new.name=paste("p",f,sep="")	
	#assign(new.name,eval((sc1env_var[[f]])))
	assign(new.name,ggplot(sc1env_var[[f]],aes(x=gs_mean_temp, y=soil_moist))+ 
		geom_bin2d(bins=10)+ 
		scale_fill_gradient(low="grey80", high=viridis(var.length+1)[f], limits=c(0,500))  )
	p_call = paste(p_call, "geom_tile(data = layer_data(",print(new.name),")) +", sep="")
}
last.p = paste("p",(var.length+1),sep="")	
assign(last.p, ggplot(allB2, aes(x=gs_mean_temp, y=soil_moist) ) + geom_bin2d(bins=10)+ scale_fill_gradient(low="grey80", high="yellow", limits=c(0,500) ) )

p_call = paste( p_call, "geom_tile(data = layer_data(",print(last.p),")) + scale_fill_identity()+ geom_hline( yintercept=c( max(allB2$soil_moist),min(allB2$soil_moist) ),linetype=\"dashed\", color =\"grey80\") +	geom_vline( xintercept=c( max(allB2$gs_mean_temp,na.rm=T),min(allB2$gs_mean_temp,na.rm=T) ),linetype=\"dashed\", color =\"grey80\")
", sep="")

#Use p_call to generate the plotting object
# combined plot using layer_data() to extract data
# from each plot above, & scale_identity to use the
# already calculated fill values
assign("p.combined", eval(parse(text=p_call)) )


 fig.name = paste("envspace_ipcc26_temp_soilmoist.pdf",sep="")
 pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)


# optional: add legends to the combined plot
library(cowplot)
plot_grid(p.combined,
          plot_grid(get_legend(p5),
                    ncol = 1),
          nrow = 1,
          rel_widths = c(1, 0.1)) # optimal relative widths may change, depending on plot dimensions



#==============================================================================
#Misc. Not sure yet if these are useful
#==============================================================================


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