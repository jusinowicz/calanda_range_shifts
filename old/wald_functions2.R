#=============================================================================
# Calculate the WALD dispersal kernel 
#
# This is a collection of functions to calculate the WALD dispersal kernel 
# Most of these functions are to calculate the vertical and horizontal 
# properties of the wind profile, following certain assumptions about the 
# foliar properties of the canopy.
#
# This is linked to the vertical wind profile measured in the field through
# the average horizontal windspeed, u_mean. 
#
#=============================================================================


#=============================================================================
#Functions to calculate wind and foliar profiles (needed for the WALD kernel)
#=============================================================================
#=============================================================================
#Get u_mean
#=============================================================================
#This script is useful for concatenating multiple csvs without the headers
#for f in `ls *.csv`; do sed '2,$!d' $f >> nehobo.csv; done
#
#This is just a function to load data files and return each site's mean
#windspeed.  
#
# site_list     A vector with the names of the files to load
# headers       Does the file contain headers? 
# col_name      If headers = T, 
#               What is the name of the column with the wind speed? 
#               If false, what column number? 

 get.u_mean=function(site_files,headers, col_name){

  nsites= length(site_files)
  site_means=matrix(0,nsites,1)

  for(n in 1:nsites){
    if (headers == TRUE){ 
      stmp = read.csv(paste(site_files[n]), header=T)
      #Fix something that happens in the big csv files for unknown reasons
      #Convert to numeric
      wsp_tmp=as.numeric(levels(stmp[,6]))[stmp[,6]]
      #Text converted to NAs
      wsp_isna = which(is.na(wsp_tmp))
      #The correct value should be in the previous column
      wsp_tmp[wsp_isna] = as.numeric(as.character(stmp [wsp_isna,5]))
      
      #Now calculate the mean
      site_means[n] = mean(wsp_tmp)

    } 

    if (headers == FALSE){ 
      stmp = read.csv(paste(site_files[n]), header=F)
      #Fix something that happens in the big csv files for unknown reasons
      #Convert to numeric
      wsp_tmp=as.numeric(levels(stmp[,6]))[stmp[,6]]
      #Text converted to NAs
      wsp_isna = which(is.na(wsp_tmp))
      #The correct value should be in the previous column
      wsp_tmp[wsp_isna] = as.numeric(as.character(stmp [wsp_isna,5]))
      
      #Now calculate the mean
      site_means[n] = mean(wsp_tmp)

    }

 }
}
#=============================================================================
#Get beta_sum
#=============================================================================
get.beta_sum = function (a1,a2,a3,a4) {
  #The normalizing integral
  F2 = function(a1,a2,a3,a4){
      f2 = function(x){
         (a1-x)^a2*(a3+x)^a4
     }
  return(f2)
  }

  beta_sum=integrate(F2 (a1,a2,a3,a4),lower=0,upper=1)$value

  return(beta_sum)
}
#=============================================================================
#Get sigma_end
#=============================================================================
#This function takes a number of additional parameters over get.sigma_w, 
#corresponding to foliar canopy characteristics. These may be calculate from 
#field data, or taken from the literature. 
#
#Default values here approximate an open meadow, and parameters come from several
#sources as cited. 
#
#LAI= 2.8     LAI is for diverse meadows see e.g.Kull and Zobel 1991
#alpha=0.06   Fixed parameter, see i.e. Massman and Weil 1999 pg 9.
#Cd=0.2       Fixed parameter, see Massman and Weil 1999
#
#This function returns sigma_end for the flux profile
get.sigma_end = function (a1,a2,a3, h, LAI=2.8, alpha=0.06, Cd=0.2, alpha1=0.6 ) {

  #Start by calculating a(z) i.e. from Massman and Weil 1999, pg 9. 
  df = 0.01
  eta=seq(0,1,df)
  pv=0.1
  a4=(a2*pv+a2*a3)/(a1-pv) #solve for a4 at a certain peak value. 
  #Note, by setting this < h, we're saying that the top of the canopy is not the peak in
  #foliar density
  beta1 = ((a1-eta)^a2*(a3+eta)^a4)
  #beta_sum = sum(beta1)
  beta_max = max(beta1)
  #The distribution of foliage
  #az=LAI/h *beta1/beta_sum
  az=LAI/h *beta1/beta_max

  #plot(az,eta)

  #Use a(z) to generate canopy flux profiles for sigma_end
  sigma_end=matrix(0, length(eta),1)
  zindex=1
  for(z in 1:length(eta)){
    sig_tmp = Cd * az[1:z]
    sigma_end[z]= sum(sig_tmp)*df
  }

  #Use sigma_end to get ustar
  c1= 0.320; c2= 0.264; c3=15.1 #From Su et al. 2001 
  ustar = c1 - c2*exp(-c3*(sigma_end[h/df]))

  sigma=list(az = az, sigma_end=sigma_end,ustar=ustar)

  return(sigma)
}

#=============================================================================
#Get the vertical variance m_sigma_w
#=============================================================================
#This function takes a number of parameters corresponding to foliar canopy
#characteristics. These may be calculate from field data, or taken from the 
#literature. 
#
#Default values here approximate an open meadow, and parameters come from several
#sources as cited. 
#
#u_mean Mean horizontal windspeed at the canopy top (from field data)
#h=1          canopy height (or sensor height) in meters
#alpha1=0.6   Fixed parameter, see Massman and Weil 1999
#ustar=0.5    Fixed parameter, see Massman and Weil 1999
#Au, Av, Aw: Constants for the variance relations. These probably could be measured in field,
#but for now could come from Massman and Weil 1999 or Katul 2005. See appendix B
#of Katul 2005 or Massman and Weil 1999. 
#

get.m_sigma_w = function (a1,a2,a3,u_mean, h=1, alpha1=0.6, Au =2.1, Av = 1.8, Aw=1.1) {

  df = 0.01 #This is the scaling factor used throughout
  nu1=(Au^2+Av^2+Aw^2)^(1/2)
  nu3=(Au^2+Av^2+Aw^2)^(3/2)
  A2=3*nu1/alpha1^2
  nu2=nu3/6-(Aw^2/(2*nu1))

  #The empirical mean of the wind profile above canopy
  #u_mean = get.u_mean()

  #Get the canopy density profile, flux profile, and ustar
  #E.g. Su et al 2012, Massman 1999
  se=get.sigma_end(a1,a2,a3,h)
  ustar = se$ustar
  sigma_end = se$sigma_end
  az = se$az

  #Calculate the energy variance, sigma_e, which is the key to calculating 
  #all of the velocity variances. See Katul et al. 2005, Massman and Weil 1999.
  gamma_z= 1 - sigma_end/sigma_end[(h/df)]
  B1= (-9*ustar/u_mean)/(2*alpha1*nu1*(9/4-A2*(ustar^4/u_mean^4)))
  n=0.5*(ustar/(u_mean))^(-2) * sigma_end[(h/df)]
  sigma_e=ustar*(nu3*exp(-A2*sigma_end[(h/df)]*gamma_z)+B1*(exp(-3*n*gamma_z)-exp(-A2*sigma_end[(h/df)]*gamma_z) ))^(1/3)

  #Now each of these follow:
  #For 1D we just need sigma_w, vertical variance: 
  sigma_w=Aw*nu1*sigma_e

  #"Integrate" 
  m_sigma_w = Aw*nu1*(sum(sigma_e,na.rm=T)*df)

  return(m_sigma_w)
}
#=============================================================================
#Get the horizontal windspeed profile, m_u_mean
#=============================================================================
#
#u_mean       Mean horizontal windspeed at the canopy top (from field data)
#a1,a2,a3     Parameters for the Beta function: E.g.  Su et al. 2001
#h=1          canopy height (or sensor height) in meters

get.m_u_mean = function ( a1,a2,a3,u_mean,h=1){
  df = 0.01
 #Get the canopy density profile, flux profile, and ustar
  #E.g. Su et al 2012, Massman 1999
  se=get.sigma_end(a1,a2,a3,h)
  ustar = se$ustar
  sigma_end = se$sigma_end
  az = se$az

  gamma_z= 1 - sigma_end/sigma_end[(h/df)]
  n=0.5*(ustar/(u_mean))^(0.5) * sigma_end[(h/df)]

  #u_mean_z=u_mean*exp(-n*gamma_z)
  m_u_mean = u_mean*sum(exp(-n*gamma_z))*df


}



#=============================================================================
#Function to calculate the WALD kernel 
#=============================================================================
#This is the primary function for the WALD calculations. 
#It returns a list object with: 
#site.kernels[[1]][[1]]   the dispersal kernel
#site.kernels[[1]][[2]]   the mean dispersal distance
#site.kernels[[1]][[3]]   the 90th percentile distance 

#This function takes the parameters:
#
# u_mean     Mean horizontal windspeed at the canopy top (from field data)
# xx         the spatial extent
# heights    the release height in meters (i.e. height of flower) by site
# Vt         the seed terminal velocity 
# a1,a2,a3    Parameters for the Beta function: E.g.  Su et al. 2001

get.WALD.kernel = function (u_mean, xx, heights, Vt,a1,a2,a3 ){
  #Do this for each species at each sites. 
  #Lattice width/height. Should be even. Half for WALD, which does not seem to 
  #handle negative values. 
  np=floor(length(xx)/2)

  xx.f=matrix(seq(0,np),np+1,length(heights))

  #The variables that will be returned in a list: fitted kernel, 
  #the means, and 90% interval of the kernels 
  site.kernels =  NULL

    #Calculate the key kernel parameters 
    m_u_mean = get.m_u_mean( a1,a2,a3,u_mean)
    m_sigma_w = get.m_sigma_w(a1,a2,a3,u_mean)
    u_prime=(heights*m_u_mean)/ Vt
    lam_prime=(heights/m_sigma_w)^2

    #Make these matrixes the appropriate size (repeat across space)
    u_prime = matrix (u_prime, (np+1), length(heights),byrow=T)
    lam_prime = matrix (lam_prime, (np+1), length(heights),byrow=T)

    #Kernel
    kdw = (lam_prime/(2*pi*xx.f^3))^(1/2)*exp(- (lam_prime*(xx.f-u_prime)^2)/(2*u_prime^2*xx.f) )
    #Normalize the kernel
    kdw=kdw/matrix(apply(kdw, 2, sum, na.rm=T),(np+1), length(heights),byrow=T )

    #Calculate the distance within which 90% of seeds lie
    new_u=(sqrt(lam_prime/xx.f)*(-xx.f+u_prime))/(sqrt(2)*u_prime)
    new_u2=(sqrt(lam_prime/xx.f)*(xx.f+u_prime))/(sqrt(2)*u_prime)
    erf = function(x) 2 * pnorm(x * sqrt(2)) - 1
    erfc = function(x) 2 * pnorm(x * sqrt(2), lower = FALSE)
    #cdf1 = 0.5*(1+erf((xx-new_u)/sqrt(2) ))+exp(2*lam_prime/u_prime) *  0.5*(1+erf((xx-(-new_u))/sqrt(2)) )
    cdf1 = 0.5*erfc(new_u)+0.5*exp(2*lam_prime/u_prime)*erfc(new_u2)
   
    s.90= apply( (xx.f*cdf1[,] <= 0.9), 2, max)
  
  site.kernels[[1]][[1]] = kdw[2:(np+1),]
  site.kernels[[1]][[2]] = u_prime
  site.kernels[[1]][[3]] = s.90
  return(site.kernels)

}
