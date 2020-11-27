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
library(dplyr)

#==============================================================================
# Load Data
#==============================================================================

#This requires user input!


# file.name.list=c(
#   "calanda_ccs85_tsmfull3_2017_waldmean.var", 
#   "calanda_ccs26_tsmfull3_2017_waldmean.var")

# variable.list=list("l1", "D", "var_mu_Us","cov_e_mu_Us",
# "cov_lam_vc", "cov_lam_vc2", "Elam1", "Elam2", "gr1.n", "gr1", "y.full", 
# "w.eq.full","ldg.sim", "ldgIMP.sim", "covIMP.sim", "Frs", "Frs_post", "nf","is.final","Frs_ci")

file.name.list=c(
  "calanda_ccs85_tsmfull4ci_2017_waldmean.var", 
  "calanda_ccs26_tsmfull4ci2_2017_waldmean.var")

variable.list=list("l1", "D", "var_mu_Us","cov_e_mu_Us",
"cov_lam_vc", "cov_lam_vc2", "Elam1", "Elam2", "gr1.n", "gr1", "y.full", 
"w.eq.full","ldg.sim", "ldgIMP.sim", "covIMP.sim", "Frs", "nf","is.final","Frs_ci")


# variable.list=list("l1", "D", "var_mu_Us","cov_e_mu_Us",
# "cov_lam_vc", "cov_lam_vc2", "Elam1", "Elam2", "gr1.n", "gr1", "y.full", 
# "w.eq.full","ldg.sim", "ldgIMP.sim", "covIMP.sim", "Frs", "nf","is.final")

# file.name.list=c(
# "calanda_ccs85_tsmfull4_noci_2017_waldmean.var",
# "calanda_ccs26_tsmfull4_noci_2017_waldmean.var")


#Old scenario file names 

# variable.list=list("l1", "D", "var_mu_Us","cov_e_mu_Us",
# "cov_lam_vc", "cov_lam_vc2", "Elam1", "Elam2", "gr1.n", "gr1", "y.full", 
# "w.eq.full","ldg.sim", "ldgIMP.sim", "covIMP.sim", "Frs", "env_analogue_all",
# "sim.impacts", "sim.impactsIMP","covIMP.impacts","lded","nf","is.final","Frs_ci")


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

# file.name.list=c(
#   "calanda_ccs85_tsmfull3_2017_waldmean.var", 
#   "calanda_ccs85B_tsmfull3_2017_waldmean.var",  
#   "calanda_ccs26_tsmfull3_2017_waldmean.var", 
#   "calanda_ccs26B_tsmfull3_2017_waldmean.var")

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
prefix.list= list("sc85tsm","sc26tsm")

#Rename the files from each scenario
var.length=length(variable.list)
nscen = length(file.name.list)


for (g in 1:nscen){
	load(file.name.list[[g]])
	for( f in 1:(var.length)) {      
		new.name=paste(prefix.list[[g]],variable.list[[f]],sep="")	
		assign(new.name,eval(as.name(variable.list[[f]])))
    rm( list=paste(variable.list[[f]]) )
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
# Figure 1: Intrinsic and realized ranges, before and after climate change
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
# mstart = 1
# mstop = 61
# minc = 20
# mtot= ceiling((mstop-mstart)/minc)
# #m.labels=as.character(seq(mstart,mstop,minc)-1)
t_use = c(1,4,7)
m.labels=c("0","30","60")
epoints = match( elevations,xx$elevation)
mtot = length(m.labels)

#For the bar plots of LGR: 
igrs1 = read.csv(file="igrs1.csv")
igrs1=cbind(c(1:7,1:7),igrs1)
#Split by IPCC scenario
igrs1a = igrs1[1:7,] #8.5
igrs1b = igrs1[8:14,] #2.6

# color.use = c( "black","red","blue" )
# par(mfrow = c(1,2))
#plt_years = c (2:7)
#plt_years = c(1,seq(2,6,by=2))
plt_years = c(1,3,6)

myrange = c(min(cbind((igrs1a[plt_years,2:4]),(igrs1b[plt_years,2:4]))),
  max(cbind((igrs1a[plt_years,2:4]),(igrs1b[plt_years,2:4]))))

myrange[1] = 0
myrange[2] = myrange[2] + 2 

nspp = 3

###Make 95% confidence intervals for the average intrinsic fitness 
for (g in 1:2){ 
    new.name1=paste(prefix.list[[g]],"l1",sep="")  
    new.name2=paste(prefix.list[[g]],"l1_ci",sep="")  
    assign("l1_post_now",eval(as.name(new.name1)))
    li_ci_now = vector("list", mtot)

for (t in 1:7){
  li_ci_now[[t]] = vector("list",3) 

  for( s in 1:nspp){
      # t1 = data_frame(num = 1:5000) %>% 
      # group_by(num) %>% 
      # mutate(means = mean(sample(l1_post_now[,t,s], replace = TRUE)))

      t1=l1_post_now[,t,s]
      li_ci_now[[t]][[s]] = quantile(t1, c(0.05,0.95) )  #sd(t1$means)  #
  }

  assign(new.name2, li_ci_now)

  }
}

 # l1_ci_lower= matrix(unlist(l1_ci)[names(unlist(l1_ci)) == "5%"], 7,nspp,byrow=T) 
 # l1_ci_upper= matrix(unlist(l1_ci)[names(unlist(l1_ci)) == "95%"], 7,nspp,byrow=T) 


###Make 95% confidence intervals for the LGRs
igrs1a_ci = vector("list", 7 )
igrs1b_ci = vector("list", 7 )
igrs1a_med = igrs1a*0
igrs1b_med = igrs1b*0


for (t in 1:7){
  igrs1a_ci[[t]] = vector("list",3) 
  igrs1b_ci[[t]] = vector("list",3) 

  for( s in 1:nspp){
     # t1 = data_frame(num = 1:10000) %>% 
     #  group_by(num) %>% 
     #  mutate(means = mean(sample(sc85tsmgr1[,t,s], replace = TRUE)))

      #t1 = sc85tsmgr1[,t,s]
      t1 = sc85tsmldg.sim[,t,s]
      igrs1a_ci[[t]][[s]] = quantile(t1, c(0.05,0.95) ) # sd(t1$means) # 
      #igrs1a_med[t,(s+1)] = median (sc85tsmldg.sim[,t,s]  )
      igrs1a_med[t,(s+1)] = (sc85tsmldg.sim[1,t,s]  )


      # t2 = data_frame(num = 1:10000) %>% 
      # group_by(num) %>% 
      # mutate(means = mean(sample(sc26tsmgr1[,t,s], replace = TRUE)))

      #t2 = sc26tsmgr1[,t,s]
      t2 = sc26tsmldg.sim[,t,s]
      igrs1b_ci[[t]][[s]] = quantile(t2, c(0.05,0.95) ) #sd(t2$means) #
      #igrs1b_med[t,(s+1)] = median (sc26tsmldg.sim[,t,s]  )
      igrs1b_med[t,(s+1)] = (sc26tsmldg.sim[1,t,s]  )

  }
}


###Make 95% confidence intervals for the FRs 
# for (g in 1:2){ 
#     new.name1=paste(prefix.list[[g]],"Frs_post",sep="")  
#     new.name2=paste(prefix.list[[g]],"Frs_ci",sep="")  
#     assign("Frs_post_now",eval(as.name(new.name1)))
#     Frs_ci_now = vector("list", mtot)

#   for (t in 1:7){
#     Frs_ci_now[[t]] = vector("list", nspp)

#     for( s in 1:nspp){
#       Frs_ci_now[[t]][[s]] = apply(Frs_post_now[[t]][[s]], 1, quantile, c(0.05,0.95))

#     }
#   }

#    assign(new.name2, Frs_ci_now)
#    rm(Frs_post_now,Frs_ci_now)
# }


#=============================================================================
# fig.name = paste("calanda_ranges_ipcc3.pdf",sep="")
# pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)
fig.name = paste("calanda_ranges_tsmRCP5_bars.pdf",sep="")
pdf(file=fig.name, height=8, width=8, onefile=TRUE, family='Helvetica', pointsize=16)

#layout.matrix=matrix(c(1,2,3,4,5,6,7,8), nrow = 4, ncol = 2)
layout.matrix=matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)
layout(mat = layout.matrix,
       heights = c(1.5, 3.5,1.5), # Heights of the rows
       widths = c(10,10)) # Widths of the columns

#layout.show(6)

#par(mfrow=c(2,1),mai= c( 0.0, 0.2, 0.0, 0.2), omi=c(0.5,0.75,0.5,0.75)) #,mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))

###Common figure properties
g=c(1,2) #Which data to use? E.g. Best-fit intrinsic ranges. 
nlevel = 64 #For viridis color scheme
ylimR1 =c(-1,200) # Common y limits for row 1 
ylimR2 =c(-1,800) # Common y limits for row 2  
xlimC1 =c(0,3500) # Common x limits for Column 1 
xlimC2 =c(0,3500) # Common x limits for Column 2
x.txt = c(1600,2000,1200) #X location of labels on ranges
par(oma = c(3,2,3,3) )

#####Column 1
### Panel 1: Bars of LGR without competition
#Just the intrinsic growth rates
l1.name=paste(prefix.list[[g[1]]],variable.list[[1]] ,sep="")
assign("l1",eval(as.name(l1.name)) )
l1_ci.name=paste(prefix.list[[g[1]]],"l1_ci" ,sep="")
assign("l1_ci",eval(as.name(l1_ci.name)) )

ci_a = unlist(l1_ci,recursive=F)
l1_ci_a =  array(sapply(ci_a,rbind), dim=c(2,nspp,7) )
l1_ci_lower= matrix(sapply(ci_a,rbind)[1,],7,nspp,byrow=T)
l1_ci_upper= matrix(sapply(ci_a,rbind)[2,],7,nspp,byrow=T)

myrange1 = c(min(l1_ci_lower[plt_years,]),max(l1_ci_upper[plt_years,]) )
#if(myrange1[1]>0){myrange1[1] = 0 }
myrange1[1] = 0
myrange1[2] = myrange1[2] + 5 

##IPCC 8.5
#Full model
#When doing stacked: 
#myrange = c(min(colSums(igrs1am_n[,2:4])),max(colSums(igrs1am_p[,2:4])))
#myrange = c(min((igrs1am_n[plt_years,2:4])),max((igrs1am_p[,2:4])))
par( mar = c(0.5,4,0,0.5) )
a1 = barplot(l1[1,plt_years,],ylim=myrange1,beside = TRUE, xaxt='n',
  ylab = expression( paste("LGR",sep=" " )),cex.lab =1.3)
barplot(l1[1,plt_years,],add=TRUE,ylim=rev(myrange1),beside = TRUE,xaxt='n' )
text(x=a1[2,], y=myrange1[2]-1, labels=c("Dactylis","Alchemilla","Helianthemum"), pos=3, xpd=NA)


segments(a1, l1_ci_lower[plt_years,], a1, l1_ci_upper[plt_years,], lwd = 1, col = "grey50")


###   Panel 2: Intrinsic range, IPCC 8.5
new.name=paste(prefix.list[[g[1]]],"Frs",sep="")  
assign("sc1Frs",eval(as.name(new.name)))
new.name=paste(prefix.list[[g[1]]],"Frs_ci",sep="")  
assign("sc1Frs_ci",eval(as.name(new.name)))
new.name2=paste(prefix.list[[g[1]]],"nf",sep="")  
assign("sc1nf",eval(as.name(new.name2)))


#par(mfrow=c(1,1),mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))

#Current: 
par( mar = c(2,4,0,0.5) )
plot(xx$elevation, sc1Frs[1,,1], t="l", xlab="", ylab="Per-capita reproduction",  xaxt="n", xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylimR1)
#plot(xx$elevation, sc1Frs[1,,1], t="l", xlab="Elevation", ylab="Seeds (per-capita, Kriged) ",   xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim)

#abline(v=c(elevations),lwd=10,col="grey80")
lines(xx$elevation,sc1Frs[1,,1])
lines(xx$elevation, sc1Frs[1,,2],col="red")
lines(xx$elevation, sc1Frs[1,,3],col="blue")
text(x=x.txt[1], y=sc1Frs[1,x.txt[1],1], labels=m.labels[1])
text(x=x.txt[2], y=sc1Frs[1,x.txt[2],2] , labels=m.labels[1],col="red")
text(x=x.txt[3], y=sc1Frs[1,x.txt[3],3] , labels=m.labels[1],col="blue")

polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[1]][[1]][1,]), sc1Frs_ci[[1]][[1]][2,]),col=alpha(rgb(t(matrix(col2rgb("grey20"))),maxColorValue=255),0.1), border = NA)
polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[1]][[2]][1,]), sc1Frs_ci[[1]][[2]][2,]),col=alpha(rgb(t(matrix(col2rgb("red"))),maxColorValue=255),0.1), border = NA)
polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[1]][[3]][1,]), sc1Frs_ci[[1]][[3]][2,]),col=alpha(rgb(t(matrix(col2rgb("blue"))),maxColorValue=255),0.1), border = NA)

#Future

#for(n in seq(3,length(m.labels),2)){
for(nn in 1:length(m.labels)) {
  n = t_use[nn]
  lines(xx$elevation, sc1Frs[n,,1],lty=2)
  text(x=x.txt[1], y=sc1Frs[n,x.txt[1],1] , labels=m.labels[nn])
  lines(xx$elevation, sc1Frs[n,,2],col="red",lty=2)
  text(x=x.txt[2], y=sc1Frs[n,x.txt[2],2] , labels=m.labels[nn],col="red")
  lines(xx$elevation, sc1Frs[n,,3],col="blue",lty=2)
  text(x=x.txt[3], y=sc1Frs[n,x.txt[3],3] , labels=m.labels[nn],col="blue")

}

axis(1, at=seq(500,3000,500),cex.axis=1)

###Panel 3
#New stacked bar plot to replace table 3? Show only the spatial mechanisms, and 
#then the model without spatial mechanisms?  

ci_a = unlist(igrs1a_ci,recursive=F)
l1_ci_a =  array(sapply(ci_a,rbind), dim=c(2,nspp,7) )
igrs1a_ci_lower= matrix(sapply(ci_a,rbind)[1,],7,nspp,byrow=T)
igrs1a_ci_upper= matrix(sapply(ci_a,rbind)[2,],7,nspp,byrow=T)

myrange = c(min(igrs1a_ci_lower[plt_years,]),max(igrs1a_ci_upper[plt_years,]) )
myrange[1] = 0

par( mar = c(0.5,4,2,0.5) )
a1= barplot(as.matrix(igrs1a[plt_years,2:4]),ylim=myrange,beside = TRUE, xaxt='n',
  ylab = expression( paste("LGR",sep=" " )),cex.lab =1.3)
# a1= barplot(as.matrix(igrs1a_med[plt_years,2:4]),ylim=myrange,beside = TRUE, xaxt='n',
#   ylab = expression( paste("LGR",sep=" " )),cex.lab =1.3)


segments(a1, igrs1a_ci_lower[plt_years,], a1, igrs1a_ci_upper[plt_years,], lwd = 1, col = "grey50")



#####Column 2

### Panel 1: Bars of LGR without competition
#Just the intrinsic growth rates
l1.name=paste(prefix.list[[g[2]]],variable.list[[1]] ,sep="")
assign("l1",eval(as.name(l1.name)) )
l1_ci.name=paste(prefix.list[[g[2]]],"l1_ci" ,sep="")
assign("l1_ci",eval(as.name(l1_ci.name)) )


par( mar = c(0.5,0.5,0,4) )
a1=barplot(l1[1,plt_years,],ylim=myrange1,beside = TRUE, axes=F,
  ylab = expression( paste("LGR",sep=" " )),cex.lab =1.3)
barplot(l1[1,plt_years,],add=TRUE,ylim=rev(myrange1),beside = TRUE, axes=F )
Axis(side=2, labels=FALSE)
text(x=a1[2,], y=myrange1[2]-1, labels=c("Dactylis","Alchemilla","Helianthemum"), pos=3, xpd=NA)

ci_a = unlist(l1_ci,recursive=F)
l1_ci_a =  array(sapply(ci_a,rbind), dim=c(2,nspp,7) )
l1_ci_lower= matrix(sapply(ci_a,rbind)[1,],7,nspp,byrow=T)
l1_ci_upper= matrix(sapply(ci_a,rbind)[2,],7,nspp,byrow=T)

segments(a1, l1_ci_lower[plt_years,], a1, l1_ci_upper[plt_years,], lwd = 1, col = "grey50")


###   Panel 2: Intrinsic range, IPCC 2.6
new.name=paste(prefix.list[[g[2]]],"Frs",sep="")  
assign("sc1Frs",eval(as.name(new.name)))
new.name=paste(prefix.list[[g[2]]],"Frs_ci",sep="")  
assign("sc1Frs_ci",eval(as.name(new.name)))
new.name2=paste(prefix.list[[g[2]]],"nf",sep="")  
assign("sc1nf",eval(as.name(new.name2)))

#Current: 
par( mar = c(2,0.5,0,4) )
plot(xx$elevation, sc1Frs[1,,1], t="l", xlab="", ylab=" ",  xaxt="n", yaxt='n',
  xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylimR1)
#plot(xx$elevation, sc1Frs[1,,1], t="l", xlab="Elevation", ylab="Seeds (per-capita, Kriged) ",   xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3,ylim=ylim)

# abline(v=c(elevations),lwd=10,col="grey80")
lines(xx$elevation,sc1Frs[1,,1])
lines(xx$elevation, sc1Frs[1,,2],col="red")
lines(xx$elevation, sc1Frs[1,,3],col="blue")
text(x=x.txt[1], y=sc1Frs[1,x.txt[1],1], labels=m.labels[1])
text(x=x.txt[2], y=sc1Frs[1,x.txt[2],2] , labels=m.labels[1],col="red")
text(x=x.txt[3], y=sc1Frs[1,x.txt[3],3] , labels=m.labels[1],col="blue")

polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[1]][[1]][1,]), sc1Frs_ci[[1]][[1]][2,]),col=alpha(rgb(t(matrix(col2rgb("grey20"))),maxColorValue=255),0.1), border = NA)
polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[1]][[2]][1,]), sc1Frs_ci[[1]][[2]][2,]),col=alpha(rgb(t(matrix(col2rgb("red"))),maxColorValue=255),0.1), border = NA)
polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[1]][[3]][1,]), sc1Frs_ci[[1]][[3]][2,]),col=alpha(rgb(t(matrix(col2rgb("blue"))),maxColorValue=255),0.1), border = NA)

#Future

#for(n in seq(3,length(m.labels),2)){
for(nn in 1:length(m.labels)) {
  n = t_use[nn]
  lines(xx$elevation, sc1Frs[n,,1],lty=2)
  text(x=x.txt[1], y=sc1Frs[n,x.txt[1],1] , labels=m.labels[nn])
  lines(xx$elevation, sc1Frs[n,,2],col="red",lty=2)
  text(x=x.txt[2], y=sc1Frs[n,x.txt[2],2] , labels=m.labels[nn],col="red")
  lines(xx$elevation, sc1Frs[n,,3],col="blue",lty=2)
  text(x=x.txt[3], y=sc1Frs[n,x.txt[3],3] , labels=m.labels[nn],col="blue")

  # polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[n]][[1]][1,]), sc1Frs_ci[[n]][[1]][2,]),col=alpha(rgb(t(matrix(col2rgb("grey80"))),maxColorValue=255),0.1), border = NA)
  # polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[n]][[2]][1,]), sc1Frs_ci[[n]][[2]][2,]),col=alpha(rgb(t(matrix(col2rgb("red"))),maxColorValue=255),0.1), border = NA)
  # polygon(c(rev(xx$elevation), xx$elevation), c(rev(sc1Frs_ci[[n]][[3]][1,]), sc1Frs_ci[[n]][[3]][2,]),col=alpha(rgb(t(matrix(col2rgb("blue"))),maxColorValue=255),0.1), border = NA)

}

axis(1, at=seq(500,3000,500),cex.axis=1)


###Panel 3:
#New stacked bar plot to replace table 3? Show only the spatial mechanisms, and 
#then the model without spatial mechanisms?  
par( mar = c(0.5,0.5,2,4) )
#par( mar = c(1,4,2,0.5) )
#Full model

barplot(as.matrix(igrs1b[plt_years,2:4]),ylim=myrange,beside = TRUE,axes=F,xaxt="n",)
#barplot(as.matrix(igrs1b_med[plt_years,2:4]),ylim=myrange,beside = TRUE,axes=F,xaxt="n",)
Axis(side=2, labels=FALSE)

ci_a = unlist(igrs1b_ci,recursive=F)
l1_ci_a =  array(sapply(ci_a,rbind), dim=c(2,nspp,7) )
igrs1b_ci_lower= matrix(sapply(ci_a,rbind)[1,],7,nspp,byrow=T)
igrs1b_ci_upper= matrix(sapply(ci_a,rbind)[2,],7,nspp,byrow=T)

segments(a1, igrs1b_ci_lower[plt_years,], a1, igrs1b_ci_upper[plt_years,], lwd = 1, col = "grey50")

dev.off()

#=============================================================================
#Figure 2: Change in LGR related to dispersal and population spread rates
#
#=============================================================================
igrs1 = read.csv(file="igrs1.csv")
#Split by IPCC scenario
igrs1a = igrs1[1:7,]
igrs1b = igrs1[8:14,]
#Make tables: Limited, 
#IPCC 8.5
#igrs1am=as.matrix(igrs1a[,2:4])-as.matrix(igrs1a[,5:7])
igrs1am=as.matrix(igrs1a[,1:3])
igrs1am_n=igrs1am
igrs1am_p=igrs1am
igrs1am_p[igrs1am_p<0] = 0
igrs1am_n[igrs1am_n>0] = 0

##IPCC 2.6
#igrs1bm=as.matrix(igrs1b[,2:4])-as.matrix(igrs1b[,5:7])
igrs1bm=as.matrix(igrs1b[,1:3])
igrs1bm_n=igrs1bm
igrs1bm_p=igrs1bm
igrs1bm_p[igrs1bm_p<0] = 0
igrs1bm_n[igrs1bm_n>0] = 0

# color.use = c( "black","red","blue" )
# par(mfrow = c(1,2))
#plt_years = c (2:7)
#plt_years = seq(2,6,by=2)
myrange = c(min(cbind((igrs1am_n[plt_years,]),(igrs1bm_n[plt_years,]))),
  max(cbind((igrs1am_p[plt_years,]),(igrs1bm_p[plt_years,]))))
if(myrange[1]>0){myrange[1]=-0.10}
if(myrange[2]<0){myrange[2]=0.10}

# Load data of LGRs with annual increments and limited pop spread rates
#==============================================================================
file.name.list=c(
  "calanda_ccs85_tsmfull4_fast_2017_waldmean.var",  
  "calanda_ccs26_tsmfull4_fast_2017_waldmean.var")
prefix.list= list("sc1","sc2")

variable.list=list("l1", "D", "var_mu_Us","cov_e_mu_Us",
"cov_lam_vc", "cov_lam_vc2", "Elam1", "Elam2", "gr1.n", "gr1", "y.full", 
"w.eq.full","ldg.sim", "ldgIMP.sim", "covIMP.sim", "Frs", "nf")

#Rename the files from each scenario
var.length=length(variable.list)
nscen = length(file.name.list)


for (g in 1:nscen){
  load(file.name.list[[g]])
  for( f in 1:(var.length)) {      
    new.name=paste(prefix.list[[g]],variable.list[[f]],sep="")  
    assign(new.name,eval(as.name(variable.list[[f]])))
    rm( list=paste(variable.list[[f]]) )
  } 
}


t_use = c(1,4,7)
m.labels=c("20","40","60")
yrs_use1=c(21,41,61)
#yrs_use1=c(1,10,20,30,40,50,60)

epoints = match( elevations,xx$elevation)
mtot = length(m.labels)

ndgens2 = 40
yrs_use2=seq(10,ndgens2,length = 3)

#igrs1ld  = sc1ldg.sim[yrs_use1,]-as.matrix(igrs1a[plt_years,2:4])
igrs1ld  = sc1gr1[yrs_use1,]
igrs1ld_p =igrs1ld
igrs1ld_n = igrs1ld
igrs1ld_p[igrs1ld_p<0] = 0
igrs1ld_n[igrs1ld_n>0] = 0

#igrs2ld  =sc2ldg.sim[yrs_use2,]-as.matrix(igrs1b[plt_years,2:4])
igrs2ld  =sc2gr1[yrs_use1,]
igrs2ld_p =igrs2ld
igrs2ld_n = igrs2ld
igrs2ld_p[igrs2ld_p<0] = 0
igrs2ld_n[igrs2ld_n>0] = 0

myrange2 = c(min(cbind((igrs1ld_n),(igrs2ld_n))),
  max(cbind((igrs1ld_p),(igrs2ld_p))))
if(myrange2[1]>0){myrange2[1]=-0.10}
if(myrange2[2]<0){myrange2[2]=0.10}
#myrange2=c(-2,2)

myrange = c(min( c(myrange[1],myrange2[1])),max(c(myrange[2],myrange2[2]) ))

#==============================================================================

fig.name = paste("calanda_limdispersalLGRs3.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

#==============================================================================
#plt_years = c(seq(2,6,by=2))
plt_years = c(2,4,6)
#plt_years = 1:7

#Version with plot per species, per decade 
#Species 1 
igrsp1 = rbind(t(as.matrix(igrs1a[plt_years,c(4,1)])), t(igrs1ld)[1,])
#Species 2 
igrsp2 = rbind(t(as.matrix(igrs1a[plt_years,c(5,2)])), t(igrs1ld)[2,])
#Species 3
igrsp3 = rbind(t(as.matrix(igrs1a[plt_years,c(6,3)])), t(igrs1ld)[3,])

layout.matrix=matrix(c(1,2,3,4,5,6), nrow = 3, ncol = 2)
layout(mat = layout.matrix,
       heights = c(3.5, 3.5,3.5), # Heights of the two rows
       widths = c(10,10)) # Widths of the two columns


###IPCC 8.5
#par(mfrow=c(3,1),oma = c(3,2,3,3))
par( mar = c(0.2,4,2,4) )

a1= barplot(igrsp1,ylim=myrange,beside = TRUE,xaxt="n",
  ylab = "",cex.lab =1.3)
barplot(igrsp1,main = "RCP 8.5", add=TRUE,ylim=rev(myrange),beside = TRUE,xaxt="n")

par( mar = c(0.2,4,0.2,4) )
a1= barplot(igrsp2,ylim=myrange,beside = TRUE,xaxt="n",
  ylab = expression( paste("LGR",sep=" " )),cex.lab =1.3)
barplot(igrsp2,add=TRUE,ylim=rev(myrange),beside = TRUE,xaxt="n")

par( mar = c(2,4,0.2,4) )
a1= barplot(igrsp3,ylim=myrange,beside = TRUE,xaxt="n",
  ylab = "",cex.lab =1.3)
barplot(igrsp3,add=TRUE,ylim=rev(myrange),beside = TRUE,xaxt="n")
#  legend = c(igrs1a[plt_years+1,1]*10),)
text(x=a1[2,], y=myrange[1], labels=c("20","40","60"), pos=1, xpd=NA,cex=1.5)


###IPCC 2.6
#Species 1 
igrsp1 = rbind(t(as.matrix(igrs1b[plt_years,c(4,1)])), t(igrs2ld)[1,])
#Species 2 
igrsp2 = rbind(t(as.matrix(igrs1b[plt_years,c(5,2)])), t(igrs2ld)[2,])
#Species 3
igrsp3 = rbind(t(as.matrix(igrs1b[plt_years,c(6,3)])), t(igrs2ld)[3,])

#par(mfrow=c(3,1),oma = c(3,2,3,3))
par( mar = c(0.2,4,2,4) )

a1= barplot(igrsp1,ylim=myrange,beside = TRUE,xaxt="n",
  ylab = "",cex.lab =1.3)
barplot(igrsp1,main = "RCP 2.6", add=TRUE,ylim=rev(myrange),beside = TRUE,xaxt="n")

par( mar = c(0.2,4,0.2,4) )
a1= barplot(igrsp2,ylim=myrange,beside = TRUE,xaxt="n",
  ylab = expression( paste("LGR",sep=" " )),cex.lab =1.3)
barplot(igrsp2,add=TRUE,ylim=rev(myrange),beside = TRUE,xaxt="n")

par( mar = c(2,4,0.2,4) )
a1= barplot(igrsp3,ylim=myrange,beside = TRUE,xaxt="n",
  ylab = "",cex.lab =1.3)
barplot(igrsp3,add=TRUE,ylim=rev(myrange),beside = TRUE,xaxt="n")
#  legend = c(igrs1a[plt_years+1,1]*10),)
text(x=a1[2,], y=myrange[1], labels=c("20","40","60"), pos=1, xpd=NA,cex=1.5)


dev.off()

#==============================================================================
#Old version:

#IPCC 8.5, global dispersal
par(mfrow=c(2,2),oma = c(3,2,3,3))
par( mar = c(0.5,4,1,1) )
a1= barplot(igrs1am_p[plt_years,],ylim=myrange,beside = TRUE,xaxt="n",
  ylab = expression( paste(Delta," LGR",sep=" " )),cex.lab =1.3)
barplot(igrs1am_n[plt_years,],add=TRUE,ylim=rev(myrange),beside = TRUE,xaxt="n")
#  legend = c(igrs1a[plt_years+1,1]*10),)
text(x=a1[2,], y=myrange[2], labels=c("Dactylis","Alchemilla","Helianthemum"), pos=3, xpd=NA)

#IPCC 2.6, global dispersal
par( mar = c(0.5,1,1,4) )
barplot(igrs1bm_p[plt_years,],ylim=myrange,beside = TRUE,axes=F,xaxt="n",)
barplot(igrs1bm_n[plt_years,],add=TRUE,ylim=rev(myrange),beside = TRUE,axes=F,xaxt="n",)
#  legend = c(igrs1a[plt_years+1,1]*10),)
Axis(side=2, labels=FALSE)
text(x=a1[2,], y=myrange[2], labels=c("Dactylis","Alchemilla","Helianthemum"), pos=3, xpd=NA)

#IPCC 8.5, limited spread
par( mar = c(0.5,4,1,1) )
a1= barplot(igrs1ld_p,ylim=myrange,beside = TRUE,xaxt="n",
  ylab = expression( paste(Delta," LGR",sep=" " )),cex.lab =1.3)
barplot(igrs1ld_n,add=TRUE,ylim=rev(myrange),beside = TRUE,xaxt="n")
#  legend = c(igrs1a[plt_years+1,1]*10),)

#IPCC 2.6, limited spread
par( mar = c(0.5,1,1,4) )
barplot(igrs2ld_p,ylim=myrange,beside = TRUE,axes=F,xaxt="n",)
barplot(igrs2ld_n,add=TRUE,ylim=rev(myrange),beside = TRUE,axes=F,xaxt="n",)
#  legend = c(igrs1a[plt_years+1,1]*10),)
Axis(side=2, labels=FALSE)


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
sr = c(0.9108988, 0.8657372, 0.8426754)
d = 1
nscen=6

t_use = c(1,4,7)
m.labels=c("0","30","60")
epoints = match( elevations,xx$elevation)
mtot = length(m.labels)

#xlims = matrix(c(-0.0, 0.15, -0.0,2.3,-0.0,2.3, -0.0,2.3),4,2,byrow=T)
#x.axx = matrix(c(0.0,0.15,0.05,0.2,2.3,0.4,0.2,2.3,0.4, 0.2,2.3,0.4 ),4,3,byrow=T)


fig.name = paste("calanda_widthOverlap_full5.pdf",sep="")
pdf(file=fig.name, height=7.5, width=7.5, onefile=TRUE, family='Helvetica', pointsize=16)

#layout.matrix=matrix(c(1,2,3,4,5,6,7,8), nrow = 4, ncol = 2)
layout.matrix=matrix(c(1,2,3,4), nrow = 2, ncol = 2)
layout(mat = layout.matrix,
       heights = c(5,5), # Heights of the rows
       widths = c(5,5)) # Widths of the columns

#layout.show(4)

#par(mfrow=c(2,1),mai= c( 0.0, 0.2, 0.0, 0.2), omi=c(0.5,0.75,0.5,0.75)) #,mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))

#Set this up as a panel of 4, or see below for individual files. 
#par( mfcol=c(2,2),mai= c( 0.25, 0.2, 0.25, 0.2), omi=c(0.5,0.75,0.5,0.75)) #,mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
#par( mfcol=c(2,2),mar= c(1, 1, 0.5, 0.5) ) #, omi=c(0.5,0.75,0.5,0.75)) #,mai= c( 1, 0, 0.2, 0), omi=c(2,0.75,2,0.75))
par(oma = c(3,2,3,3) )


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
for( g in 1:2){ 


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

  for(si in 1:mtot) {
    
    s = t_use[si]
    # min_overlap1=0
    # xlimsa=0
    # ylimsa=0


    for(u in 1:nspp){


      #Create a mean standardization coefficient to simplify presentation. 
      #C1 =  median((l1[d,1,u]/(D[d,1,u])^3), na.rm=T)
      #C2 =  median(l1[d,1,u]/(D[d,1,u])^2, na.rm=T )

      C1 =  median((l1[d,1,u]/(D[d,1,u])), na.rm=T)
      C2 =  median(l1[d,1,u]/(D[d,1,u]), na.rm=T )

      #Spatial mechanisms add this much
      width_y1 = C1*var_mu_Us[,1,u] 
      width_y1_ci = quantile(width_y1, c(0.05,0.95))
      #Strongest negative impact (this is only accurate if data include case of complete overlap)
      #min_overlap1 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
      overlap_x1 = ((-C2*cov_e_mu_Us[,1,u])+cov_lam_vc[,1,u])
      overlap_x1_ci = quantile(overlap_x1 ,c(0.05,0.95) )

      #Spatial mechanisms add this much
      width_y2 = C1*var_mu_Us[,s,u] 
      width_y2_ci = quantile(width_y1, c(0.05,0.95))
      #Strongest negative impact (this is only accurate if data include case of complete overlap)
      #min_overlap1 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
      overlap_x2 = ((-C2*cov_e_mu_Us[,s,u])+cov_lam_vc[,s,u])
      overlap_x2_ci = quantile(overlap_x1 ,c(0.05,0.95) )



      min_overlap1 = min(c (min_overlap1, median(overlap_x1) ,median(overlap_x2) ),na.rm=T)
      xlimsa = max(c (xlimsa,overlap_x2_ci["95%"],overlap_x1_ci["95%"],na.rm=T))
      ylimsa = max(c (ylimsa,width_y2_ci["95%"],width_y2_ci["95%"],na.rm=T))
      #print(c (xlimsa,(-C2*cov_e_mu_Us[d,s,u]+cov_lam_vc[d,s,u]),(-C2*cov_e_mu_Us[d,1,u]+cov_lam_vc[d,1,u])))

    }
  }
}
  xlims[,2] = round(-min_overlap1+xlimsa+0.2,2)
  ylims[,2] = round(ylimsa+0.1,1)
  xlims[,1] = -1
  x.axx= round(matrix( c(xlims[,1],xlims[,2],(( xlims[,2]-xlims[,1]) / 3) ),ngenst,3),2)
  y.axx= round(matrix( c(ylims[,1],ylims[,2],(( ylims[,2]-ylims[,1]) / 3) ),ngenst,3),2)

#Uncomment this to make all of the limits the same
# xlims = matrix(c(0,round(-min_overlap1+xlims+0.2,1)),ngenst,2,byrow=T)
# ylims = matrix(c(0,round(ylims+0.1,1)),ngenst,2,byrow=T)
# x.axx = matrix( c(xlims[,1],xlims[,2],(( xlims[,2]-xlims[,1]) / 5) ),ngenst,3)
# y.axx = matrix( c(ylims[,1],ylims[,2],(( ylims[,2]-ylims[,1]) / 5) ),ngenst,3)

#Now make the plots: 

for( g in 1:2){ 
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
assign("w.eq",eval(as.name(w.eq.name)) )




#for( s in seq(3,length(m.labels),2)){ 
for( si in 2:mtot){ 

  s = t_use[si]
  color.use=list("black","red","blue")
  up.low=list("upper","lower")

  #Or just switch between devices for plotting
  #dev.set(g)
  #par(mfrow=c(1,4), mar=c(5,6,6,2))

  ylim=ylims[s,]
  #xlim=c(0,0.8)
  xlim=xlims[s,]
  #plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Range width", xlab="Range overlap",xaxs="i",yaxs="i",cex.main=1.3,cex.lab=1.3 )
  if ( g ==1) {   

    if ((si-1)%%(mtot-1) == 0){
      #Lower left
      par( mar = c(4,4,0,0.5) )
      plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="", 
      xlab="Benefit of Competitor Segregation",xaxs="i",yaxs="i",cex.main=1.2,cex.lab=1.2,yaxt='n',xaxt='n' )
      axis(side=1, labels=T, at=c(seq(x.axx[s,1], x.axx[s,2],x.axx[s,3])[-1],0))
      axis(side=2, labels=T, at=seq(y.axx[s,1], y.axx[s,2],y.axx[s,3]))
      abline(v = 0 )

    } else { 
      #Upper left
      par( mar = c(0.5,4,4,0.5) )
      plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="Benefit of Competitor Aggregation", 
      xlab="",xaxs="i",yaxs="i",cex.main=1.2,cex.lab=1.2,yaxt='n',xaxt='n')
      axis(side=1, labels=F, at=c(seq(x.axx[s,1], x.axx[s,2],x.axx[s,3])[-1],0))
      axis(side=2, labels=T, at=seq(y.axx[s,1], y.axx[s,2],y.axx[s,3]))
      abline(v = 0)
    } 
      
  } else {

    if ((si-1)%%(mtot-1) == 0){
      #Lower right
        par( mar = c(4,0.5,0,4) )
        plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="", 
        xlab="",xaxs="i",yaxs="i",cex.main=1.2,cex.lab=1.2,yaxt='n',xaxt='n' )
        axis(side=1, labels=T,at=c(seq(x.axx[s,1], x.axx[s,2],x.axx[s,3])[-1],0))
        axis(side=2, labels=F, at=seq(y.axx[s,1], y.axx[s,2],y.axx[s,3]))
        abline(v = 0)

      } else { 
        #Upper right
        par( mar = c(0.5,0.5,4,4) )
        plot( 0,0, type="n", xlim=xlim, ylim=ylim, ylab="", 
        xlab="",xaxs="i",yaxs="i",cex.main=1.2,cex.lab=1.2,yaxt='n',xaxt='n')
        axis(side=1, labels=F, at=c(seq(x.axx[s,1], x.axx[s,2],x.axx[s,3])[-1],0))
        axis(side=2, labels=F, at=seq(y.axx[s,1], y.axx[s,2],y.axx[s,3]))
       abline(v = 0)
      } 

  }


  #mlin.non1 =0.3229569
  #lines((seq(xlim[1],xlim[2],by=0.1)),lin.non1,lty=2,lwd=1.5  )

  for (u in 1:nspp){ 
      
      
      #This line represents persistence in the absence of spatial mechanisms
      # mlin.non1 =1-(l1[d,iconfig,u]/D[d,iconfig,u]+sr[u])-min_overlap1
      # #mlin.non1 =0.3229569
      # lin.non1= -(seq(xlim[1],xlim[2],by=0.05))+mlin.non1
      # lines((seq(xlim[1],xlim[2],by=0.05)),lin.non1,col=alpha(rgb(t(matrix(col2rgb(color.use[[u]]))),maxColorValue=255),0.4),lty=5,lwd=1.5  )


      #Create a mean standardization coefficient to simplify presentation. 
      #C1 =  median((l1[,u]/(D[,u])^3), na.rm=T)
      #C2 =  median(l1[,u]/(D[,u])^2, na.rm=T )

      #C1 = (l1[d,1,u]/(D[d,1,u])^3)
      #C2 = l1[d,1,u]/(D[d,1,u])^2

      C1 =  median((l1[d,iconfig,u]/(D[d,iconfig,u])), na.rm=T)
      C2 =  median(l1[d,iconfig,u]/(D[d,iconfig,u]), na.rm=T )


      #Spatial mechanisms add this much
      width_y1 = C1*var_mu_Us[,iconfig,u] 
      width_y1_ci = quantile(width_y1, c(0.05,0.95))
      #Strongest negative impact (this is only accurate if data include case of complete overlap)
      #min_overlap1 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
      overlap_x1 = -min_overlap1+((-C2*cov_e_mu_Us[,iconfig,u])+cov_lam_vc[,iconfig,u])
      overlap_x1_ci = quantile(overlap_x1 ,c(0.05,0.95) )

      #Point where lines intercept: x = (b1-b2)/2
      #non.x1 = (mlin.non1 - (width_y1 - overlap_x1) )/2
      #non.y1 = -non.x1+ mlin.non1 
      #Intersection through origin instead:
      # non.x1 = -mlin.non1/(-1-(width_y1/overlap_x1) )
      # non.y1 = -non.x1+ mlin.non1 
      #points(non.x1,non.y1,col=color.use[[u]], pch=16)
      #points(overlap_x1[d],width_y1[d],col=color.use[[u]] )
      points(median(overlap_x1),median(width_y1),col=color.use[[u]] )

      #segments(non.x1,non.y1,overlap_x1,width_y1,col=color.use[[u]],lty=2)     
      segments( overlap_x1_ci["5%"], median(width_y1), overlap_x1_ci["95%"], median(width_y1), lwd = 1,col = "grey50")
      segments(  median(overlap_x1), width_y1_ci["5%"], median(overlap_x1),width_y1_ci["95%"],  lwd = 1,col = "grey50")

      #Now the later point
      # mlin.non2 =1-(l1[d,s,u]/D[d,s,u]+sr[u])-min_overlap1
      # #mlin.non2 =0.3229569
      # lin.non2= -(seq(xlim[1],xlim[2],by=0.1))+mlin.non2
      # lines((seq(xlim[1],xlim[2],by=0.1)),lin.non2,col=alpha(rgb(t(matrix(col2rgb(color.use[[u]]))),maxColorValue=255),0.4),lty=3,lwd=1.5  )
      # #lines((seq(xlim[1],xlim[2],by=0.1)),lin.non2,lty=2,lwd=1.5  )

      #Spatial mechanisms add this much
      width_y2 = C1*var_mu_Us[,s,u]
      width_y2_ci = quantile(width_y2, c(0.05,0.95))

      #Strongest negative impact (this is only accurate if data include case of complete overlap)
      #min_overlap2 = min((-C2*cov_e_mu_Us[,u]+cov_lam_vc[,u]),na.rm=T)
      overlap_x2 = -min_overlap1+((-C2*cov_e_mu_Us[,s,u])+cov_lam_vc[,s,u])
      overlap_x2_ci = quantile(overlap_x2 ,c(0.05,0.95) )

      # non.x2 = -mlin.non2/(-1-(width_y2/overlap_x2) )
      # non.y2 = -non.x2+ mlin.non2 
      #points(non.x2,non.y2 )
      points(median(overlap_x2[]),median(width_y2[]),col=color.use[[u]])
      #segments(non.x2,non.y2,overlap_x2,width_y2,col=color.use[[u]],lty=2)     
      segments( overlap_x2_ci["5%"], median(width_y2), overlap_x2_ci["95%"], median(width_y2), lwd = 1,,col = "grey50")
      segments(  median(overlap_x2), width_y2_ci["5%"], median(overlap_x2),width_y2_ci["95%"],  lwd = 1,,col = "grey50")


      #segments(non.x1,non.y1,new.x2,new.y2,col=color.use[[u]],lty=2)
      arrows( median(overlap_x1), median(width_y1), median(overlap_x2), median(width_y2),col=color.use[[u]],length=0.1,lwd=2)
    } 
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
#New stacked bar plot to replace table 3? Show only the spatial mechanisms, and 
#then the model without spatial mechanisms?  
#==============================================================================
igrs1 = read.csv(file="igrs1.csv")
#Split by IPCC scenario
igrs1a = igrs1[1:7,]
igrs1b = igrs1[8:14,]
#Make tables: Limited, global dispersal, no spatial mechanisms
igrs1am=(as.matrix(igrs1a)-matrix(as.matrix(igrs1a[1,]),7,10,byrow=T))
igrs1am_n=igrs1am
igrs1am_p=igrs1am
igrs1am_p[igrs1am_p<0] = 0
igrs1am_n[igrs1am_n>0] = 0

color.use = c( "black","red","blue" )
par(mfrow = c(1,2))
#plt_years = c (2:7)
plt_years = seq(2,6,by=2)
myrange = c(min(cbind((igrs1am_n[plt_years,2:4]),(igrs1bm_n[plt_years,2:4]))),
  max(cbind((igrs1am_p[plt_years,2:4]),(igrs1bm_p[plt_years,2:4]))))
#myrange = c(-1,
#  max(cbind((igrs1am_p[plt_years,2:4]),(igrs1bm_p[plt_years,2:4]))))



##IPCC 8.5
#Full model
#When doing stacked: 
#myrange = c(min(colSums(igrs1am_n[,2:4])),max(colSums(igrs1am_p[,2:4])))
#myrange = c(min((igrs1am_n[plt_years,2:4])),max((igrs1am_p[,2:4])))
barplot(igrs1am_p[plt_years,2:4],ylim=myrange,beside = TRUE, 
  ylab = expression( paste(Delta," LGR",sep=" " )),cex.lab =1.5)
barplot(igrs1am_n[plt_years,2:4],add=TRUE,ylim=rev(myrange),beside = TRUE, )

#No spatial mechanisms
#myrange = c(min((igrs1am_n[plt_years,8:10])),max((igrs1am_p[,8:10])))
# barplot(igrs1am_p[plt_years,8:10],ylim=myrange,legend = c(igrs1a[plt_years,1]),beside = TRUE, )
# barplot(igrs1am_n[plt_years,8:10],add=TRUE,ylim=rev(myrange),beside = TRUE, )

##IPCC 2.6
igrs1bm=(as.matrix(igrs1b)-matrix(as.matrix(igrs1b[1,]),7,10,byrow=T))
igrs1bm_n=igrs1bm
igrs1bm_p=igrs1bm
igrs1bm_p[igrs1bm_p<0] = 0
igrs1bm_n[igrs1bm_n>0] = 0

#Full model
#myrange = c(min((igrs1bm_n[,2:4])),max((igrs1bm_p[plt_years,2:4])))
barplot(igrs1bm_p[plt_years,2:4],ylim=myrange,beside = TRUE,)
barplot(igrs1bm_n[plt_years,2:4],add=TRUE,ylim=rev(myrange),beside = TRUE,
  legend = c(igrs1a[plt_years+1,1]*10),)

#No spatial mechanisms
#myrange = c(min((igrs1bm_n[plt_years,8:10])),max((igrs1bm_p[,8:10])))
# barplot(igrs1bm_p[plt_years,8:10],ylim=myrange,legend = c(igrs1b[plt_years,1]),beside = TRUE,)
# barplot(igrs1bm_n[plt_years,8:10],add=TRUE,ylim=rev(myrange),beside = TRUE,)


#==============================================================================
#New stacked bar plot to replace table 3? 
#==============================================================================
igrs1 = read.csv(file="igrs1.csv")
#Split by IPCC scenario
igrs1a = igrs1[1:7,]
igrs1b = igrs1[8:14,]
#Make tables: Limited, global dispersal, no spatial mechanisms
igrs1am=(matrix(as.matrix(igrs1a[1,]),7,10,byrow=T)-as.matrix(igrs1a))
igrs1am_n=igrs1am
igrs1am_p=igrs1am
igrs1am_p[igrs1am_p<0] = 0
igrs1am_n[igrs1am_n>0] = 0

#par(mfrow = c(2,3))
par(mfrow = c(1,3))
#plt_years = c (2:7)
plt_years = c(3,7,by=2)

##IPCC 8.5
#Full model
#When doing stacked: 
#myrange = c(min(colSums(igrs1am_n[,2:4])),max(colSums(igrs1am_p[,2:4])))
myrange = c(min((igrs1am_n[plt_years,2:4])),max((igrs1am_p[,2:4])))
barplot(igrs1am_p[plt_years,2:4],ylim=myrange,beside = TRUE, ylab = expression( paste(Delta," LGR",sep=" " )),cex.lab =1.5)
barplot(igrs1am_n[plt_years,2:4],add=TRUE,ylim=rev(myrange),beside = TRUE, )

#Global dispersal
myrange = c(min((igrs1am_n[plt_years,5:7])),max((igrs1am_p[,5:7])))
barplot(igrs1am_p[plt_years,5:7],ylim=myrange,beside = TRUE, )
barplot(igrs1am_n[plt_years,5:7],add=TRUE,ylim=rev(myrange),beside = TRUE, )

#No spatial mechanisms
myrange = c(min((igrs1am_n[plt_years,8:10])),max((igrs1am_p[,8:10])))
barplot(igrs1am_p[plt_years,8:10],ylim=myrange,legend = c(igrs1a[plt_years,1]),beside = TRUE, )
barplot(igrs1am_n[plt_years,8:10],add=TRUE,ylim=rev(myrange),beside = TRUE, )

##IPCC 2.6
igrs1bm=(matrix(as.matrix(igrs1b[1,]),7,10,byrow=T)-as.matrix(igrs1b))
igrs1bm_n=igrs1bm
igrs1bm_p=igrs1bm
igrs1bm_p[igrs1bm_p<0] = 0
igrs1bm_n[igrs1bm_n>0] = 0

#Full model
myrange = c(min((igrs1bm_n[,2:4])),max((igrs1bm_p[plt_years,2:4])))
barplot(igrs1bm_p[plt_years,2:4],ylim=myrange,beside = TRUE,)
barplot(igrs1bm_n[plt_years,2:4],add=TRUE,ylim=rev(myrange),beside = TRUE,)

#Global dispersal
myrange = c(min((igrs1bm_n[plt_years,5:7])),max((igrs1bm_p[plt_years,5:7])))
barplot(igrs1bm_p[plt_years,5:7],ylim=myrange,beside = TRUE,)
barplot(igrs1bm_n[plt_years,5:7],add=TRUE,ylim=rev(myrange),beside = TRUE,)

#No spatial mechanisms
myrange = c(min((igrs1bm_n[plt_years,8:10])),max((igrs1bm_p[,8:10])))
barplot(igrs1bm_p[plt_years,8:10],ylim=myrange,legend = c(igrs1b[plt_years,1]),beside = TRUE,)
barplot(igrs1bm_n[plt_years,8:10],add=TRUE,ylim=rev(myrange),beside = TRUE,)


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