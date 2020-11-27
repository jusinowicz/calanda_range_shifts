####################################################
# Try fitting some models 
#
#=============================================================================
# Load these libraries
#=============================================================================
library(mgcv)
library(zoo)
library(ncdf4)
library(fields)
#=============================================================================
# Function definitions 
#=============================================================================


#=============================================================================
# Load data
#=============================================================================
# site_files=c("/home/jacob/Desktop/working/zurich/data/hobo_data/Aurella/arhobo.csv",
# 			"/home/jacob/Desktop/working/zurich/data/hobo_data/Neselboden/nshobo.csv",
# 			"/home/jacob/Desktop/working/zurich/data/hobo_data/Barenmos/bmhobo.csv",
# 			"/home/jacob/Desktop/working/zurich/data/hobo_data/Neusass/nehobo.csv",
# 			"/home/jacob/Desktop/working/zurich/data/hobo_data/Calanda/cahobo.csv")

site_files=c("/home/jacob/Desktop/working/zurich/data/hobo_data/arhobo.csv",
			"/home/jacob/Desktop/working/zurich/data/hobo_data/nshobo.csv",
			"/home/jacob/Desktop/working/zurich/data/hobo_data/bmhobo.csv",
			"/home/jacob/Desktop/working/zurich/data/hobo_data/nehobo.csv",
			"/home/jacob/Desktop/working/zurich/data/hobo_data/cahobo.csv")

data_names = c( "arhobo","nshobo", "bmhobo","nehobo","cahobo")


#=============================================================================
#Load each of the files and calculate the following per site after
#some data processing/formatting: 
#	1. ma_daily 	The daily average of temperature
#	2. rmin_daily 	The daily minimum of temperature
#	3. rmax_daily	The daily maximum of temperature 
#
#These also lead to the further calculation of: 
#	4. gdd_cum/_all	The cumulative number of Growing Degree Days
#=============================================================================


for (of in 1:length(site_files)){

	hb_tmp = read.csv(paste(site_files[of]),header=F)
	#Remove NAs (why are these here?)
	hb_tmp = hb_tmp[-which(is.na(hb_tmp[,3])),]

	#Procb_tmp=ess the dates:
	#dates=as.Date(as.character(hb_tmp[,6]), "%Y-%m-%d")
	#dates = as.Date(as.character(hb_tmp[,2]), "%m/%d/%y %I:%M:%S %p")
	dates=(strptime(hb_tmp[,2], format="%m/%d/%y %I:%M:%S %p"))
	dates2= as.Date(as.character(hb_tmp[,2]), "%m/%d/%y %I:%M:%S %p")

	#hb_tmp[,6] = dates
	day = (format(dates, "%d"))
	month = (format(dates, "%m"))
	year =(format(dates, "%Y"))

	#Remake and order the data frame
	hb_tmp = data.frame(dates2, dates,hb_tmp$V3,hb_tmp$V4)
	hb_tmp=hb_tmp[order(hb_tmp$dates),]


	nyears=unique(year)
	#nyears=nyears[-which(is.na(nyears))]

	nmonths = unique(month)
	#nmonths=nmonths[-which(is.na(nmonths))]

	ndays = unique(day)
	#ndays=ndays[-which(is.na(ndays))]

	#Caclulate a daily average
	ma_daily=aggregate(. ~ dates2, hb_tmp,FUN=mean)
	rmin_daily=aggregate(. ~ dates2, hb_tmp,FUN=min)
	rmax_daily=aggregate(. ~ dates2, hb_tmp,FUN=max)

	#Add columns for each (day.month.year)
	hb_tmp=data.frame(hb_tmp,day)
	hb_tmp=data.frame(hb_tmp,month)
	hb_tmp=data.frame(hb_tmp,year)

	gdd_all = NULL
	gdd_per = NULL
	gdd1 = NULL
	
	for(yr in 1:length(nyears)){
		
		ma_tmp = subset(ma_daily, format(ma_daily[,1], "%Y") ==nyears[yr])
		hb_match = subset(hb_tmp, format(hb_tmp[,1], "%Y") ==nyears[yr])

		#Growing degree days. 
		#Use the MA to calculate GDD with the assumption that 
		#the baseline temp is 10C and the max is 30
		tmin=10
		tmax=30
		
		ma_tmp[,3][ma_tmp[,3]<tmin] = tmin
		ma_tmp[,3][ma_tmp[,3]>tmax] = tmax

		gdd_tmp = ma_tmp[,3]-tmin
		gdd_tmp[is.na(gdd_tmp)]=0

		#Total gdd 
		gdd=cumsum(gdd_tmp)
		gdd1 = rbind(gdd1, matrix(gdd,length(gdd),1))

		#Set up a lookup table to match each quantity to the appropriate
		#hourly entry in the full data set.
		gdd_tab = data.frame(unique(ma_tmp[,1]),gdd)
		gdd_sum = gdd_tab[match(hb_match[,1],gdd_tab[,1]),]
		gdd_all = rbind(gdd_all, matrix(gdd_sum[,2], length(gdd_sum[,2]),1))

		#Each day's contribution to GDD:
		gdd_per = rbind(gdd_per, matrix(gdd_tmp, length(gdd_tmp),1))
	
	}

	#Add the degree days to the main data object
	hb_tmp=data.frame(hb_tmp,gdd_all)

	#Assign each new data frame to its appropriate site
	assign(data_names[of], eval(as.name("hb_tmp")))
	
	#Create a second data frame that has the daily measures
	daily_tmp = data.frame(unique(dates2),rmin_daily[,3],ma_daily[,3],rmax_daily[,3],
				matrix(gdd_per, length(gdd_per),1), gdd1, rmin_daily[,4],ma_daily[,4],rmax_daily[,4])
	colnames(daily_tmp) = c("date","temp_min","temp_mean","temp_max","daily_gdd","cum_gdd","light_min","light_mean","light_max")
	assign( paste(data_names[of],"_daily",sep=""), eval(as.name("daily_tmp")))
}


#Recursively save each of the data frames as an R variable
for (of in 1:length(site_files)){
	save( list=data_names[of],file=paste(data_names[of],".var",sep=""))
	save( list=paste(data_names[of],"_daily",sep=""),file=paste(data_names[of],"_daily",".var",sep=""))
}

for (of in 1:length(site_files)){
	load(file=paste(data_names[of],".var",sep=""))
	load(file=paste(data_names[of],"_daily",".var",sep=""))
}


#pdf(file="degree_days_all.pdf", height=7.5, width=7.5,family='Helvetica', pointsize=16)
png(file="degree_days_all.png", height=7.5, width=7.5, units="in", res=600, family='Helvetica', pointsize=16)

 par(mfrow=c(1,1), mar=c(5,6,6,2))
 plot(arhobo$dates2, arhobo$gdd_all, ylab="Cumulative GDD", xlab="Date")
 points(nshobo$dates2, nshobo$gdd_all, col="red")
 points(bmhobo$dates2, bmhobo$gdd_all, col="blue")
 points(nehobo$dates2, nehobo$gdd_all, col="green")
 points(cahobo$dates2, cahobo$gdd_all, col="orange")

dev.off()

#=============================================================================
# Use the data on daily average, min, max, and GDD to find predictors
# of intrinsic growth rates for each species. 
#=============================================================================
#Standardize all of the dates to one another. Compare AR and NS to get 909 dates
#dates_use = unique(arhobo_daily[,1])[unique(arhobo_daily[,1])%in%unique(nshobo_daily[,1])]
elevations = c(600,1000,1400,1600,1800,2000,2100,2700)
ne= 1

dates_use =  unique(nshobo_daily[,1])[format(nshobo_daily[,1], "%Y")>2014 & format(nshobo_daily[,1], "%Y")< 2018]

arhobo_sub = arhobo_daily[arhobo_daily[,1]%in%dates_use,]
nshobo_sub = nshobo_daily[nshobo_daily[,1]%in%dates_use,]
bmhobo_sub = bmhobo_daily[bmhobo_daily[,1]%in%dates_use,]
nehobo_sub = nehobo_daily[nehobo_daily[,1]%in%dates_use,]
cahobo_sub = cahobo_daily[cahobo_daily[,1]%in%dates_use,]

#Combine into one giant dataframe with elevations 
arhobo_sub = data.frame(arhobo_sub, elevations=matrix(elevations[ne+1],dim(arhobo_sub)[1]))
nshobo_sub = data.frame(nshobo_sub, elevations=matrix(elevations[ne+2],dim(nshobo_sub)[1]))
bmhobo_sub = data.frame(bmhobo_sub, elevations=matrix(elevations[ne+3],dim(bmhobo_sub)[1]))
nehobo_sub = data.frame(nehobo_sub, elevations=matrix(elevations[ne+4],dim(nehobo_sub)[1]))
cahobo_sub = data.frame(cahobo_sub, elevations=matrix(elevations[ne+5],dim(cahobo_sub)[1]))

all_sub = rbind(arhobo_sub,nshobo_sub,bmhobo_sub,nehobo_sub,cahobo_sub)

arhobo=cbind(arhobo, elevations = matrix(1000,dim(arhobo)[1],1))
nshobo=cbind(nshobo, elevations = matrix(1400,dim(nshobo)[1],1))
bmhobo=cbind(bmhobo, elevations = matrix(1600,dim(bmhobo)[1],1))
nehobo=cbind(nehobo, elevations = matrix(1800,dim(nehobo)[1],1))
cahobo=cbind(cahobo, elevations = matrix(2000,dim(cahobo)[1],1))

all_sub2 = rbind(arhobo,nshobo,bmhobo,nehobo,cahobo)

#=============================================================================
#Soil moisture data 
#Potential evapotranspiration data
#This comes from the Climatology Lab at the University of Idaho
#http://www.climatologylab.org/terraclimate.html
#
#The functions below are taken from the example R code at 
#http://www.climatologylab.org/uploads/2/2/1/3/22133936/read_terraclimate.r
#
#=============================================================================

#Added a loop to go over each site: 
#Real:
#lat_long = data.frame(lat = c(46.874, 46.869, 46.878, 46.887, 46.888), lon = c(9.509, 9.490,9.494, 9.495, 9.489))
#This set of coordinates produces a gradient 
#lat_long = data.frame(lat = c(46.874, 46.869, 46.878, 46.887, 46.95), lon = c(9.509, 9.490,9.494, 9.495, 9.489))

#This set produces a longer gradient 
#lat_long = data.frame(lat = c(46.83788,46.874, 46.869, 46.878, 46.887, 46.95,46.90979,46.88997), lon = c(9.48069,9.509, 9.490,9.494, 9.495, 9.489,9.48682,9.47481))

#This set produces a longer gradient, but starts on a peak, goes across the valley, then up again
#lat_long = data.frame(lat = c(46.66, 46.7, 46.74, 46.78, 46.81, 46.85, 46.88  ,46.91,46.94), lon = c(9.490,9.490, 9.490,9.490, 9.490, 9.490,9.490,9.490,9.490))

#This set combines these
#lat_long = data.frame(lat = c(46.85, 46.88  ,46.91), lon = c(9.490,9.490, 9.490,9.490, 9.490, 9.490,9.490,9.490,9.490))
#elevations = c(562,1800,2100)
#Best currently: 
lat_long = data.frame(lat = c(46.85, 46.874, 46.869, 46.878, 46.887, 46.95,46.91,46.89750), lon = c(9.509, 9.509, 9.490,9.494, 9.495, 9.489,9.490,9.46858))
elevations = c(600,1000,1400,1600,1800,2000,2100,2700)


moist_yearly = (expand.grid(elevations, years))
moist_yearly=data.frame(cbind(moist_yearly,matrix(0,dim(moist_yearly[1]),1)))
colnames(moist_yearly) = c('elevations','years','soil_moist')

#Added a loop to go over each year: 

for (yy in 1:3) {
	
	year = years[yy]

	# enter in variable you want to download see: http://thredds.northwestknowledge.net:8080/thredds/terraclimate_aggregated.html
	var="soil"

	#baseurlagg <- paste0(paste0("http://thredds.northwestknowledge.net:8080/thredds/dodsC/agg_terraclimate_",var),"_1958_CurrentYear_GLOBE.nc")
	basefile = paste("soil_moisture/TerraClimate_", var, "_", year,".nc", sep="")

	nc = nc_open(basefile)
	
	lat =ncvar_get(nc, varid="lat")
	lon =ncvar_get(nc, varid="lon")
	#time =ncvar_get(nc, varid="time")
	soil = ncvar_get(nc, varid="soil")

	for( ll in 1:(length(elevations))) {

		# enter in longitude, latitude here
		x=c(unlist(lat_long[ll,]))

		# time_d = as.Date(time, format="%d", origin=as.Date( "1899-01-01" ) )
		# time_years = format(time_d, "%Y")
		# time_months = format(time_d, "%m")
		# time_ym = format(time_d, "%Y-%m")

		flat = match(abs(lat - x[1]) < 1/48, 1)
		latindex = which(flat %in% 1)
		flon = match(abs(lon - x[2]) < 1/48, 1)
		lonindex = which(flon %in% 1)
		start =c(lonindex, latindex, 1)
		count = c(1, 1, -1)

		# #To plot an area:
		lat_range <- seq(latindex-10, latindex+10)
		lon_range <- seq(lonindex-10, lonindex+10)
		#image.plot(lon_range, lat_range,soil[lon_range,c(max(lat_range):min(lat_range)),1])
		#mat1 = apply(soil[lon_range,c(min(lat_range):max(lat_range)),1], 2, rev)
		
		if(ll==1){
			image.plot(lon_range, lat_range,soil[lon_range,c(min(lat_range):max(lat_range)),1])
		}
		#abline(h=latindex)
		#abline(v=lonindex)
		#points(lonindex,latindex )
		text(lonindex,latindex, ll)
		# # #image.plot((soil[,c(4320:1),1]))
		# grid <- expand.grid(lon=lon_range, lat=lat_range)
		# levelplot(soil ~ lon * lat, data=grid, at=cutpts, cuts=11, pretty=T, col.regions=(rev(brewer.pal(10,"RdBu"))))




		# read in the full period of record using aggregated files
		data <- as.numeric(ncvar_get(nc, varid = var,start = start, count))
		moist_yearly[elevations %in% elevations[ll], ][years %in% years[yy],]$soil_moist = mean(data,na.rm=T)

	}

}

########################################################################################################################
#PET
########################################################################################################################

#Added a loop to go over each site: 
#lat_long = data.frame(lat = c(46.874, 46.869, 46.878, 46.887, 46.888), lon = c(9.509, 9.490,9.494, 9.495, 9.489))

#This set of coordinates produces a gradient 
#lat_long = data.frame(lat = c(46.874, 46.869, 46.878, 46.887, 46.95), lon = c(9.509, 9.490,9.494, 9.495, 9.489))
#Best currently: 
lat_long = data.frame(lat = c(46.85, 46.874, 46.869, 46.878, 46.887, 46.95,46.91,46.89750), lon = c(9.509, 9.509, 9.490,9.494, 9.495, 9.489,9.490,9.46858))
elevations = c(600,1000,1400,1600,1800,2000,2100,2700)

 
pet_yearly = (expand.grid(elevations, years))
pet_yearly=data.frame(cbind(pet_yearly,matrix(0,dim(pet_yearly[1]),1)))
colnames(pet_yearly) = c('elevations','years','soil_pet')

#Added a loop to go over each year: 

for (yy in 1:3) {
	
	year = years[yy]

	# enter in variable you want to download see: http://thredds.northwestknowledge.net:8080/thredds/terraclimate_aggregated.html
	var="pet"

	basefile = paste("soil_moisture/TerraClimate_", var, "_", year,".nc", sep="")
	nc = nc_open(basefile)
	lat =ncvar_get(nc, varid="lat")
	lon =ncvar_get(nc, varid="lon")
	#time =ncvar_get(nc, varid="time")
	pet = ncvar_get(nc, varid="pet")


	for( ll in 1:5) {

		# enter in longitude, latitude here
		x=c(unlist(lat_long[ll,]))



		# time_d = as.Date(time, format="%d", origin=as.Date( "1899-01-01" ) )
		# time_years = format(time_d, "%Y")
		# time_months = format(time_d, "%m")
		# time_ym = format(time_d, "%Y-%m")

		flat = match(abs(lat - x[1]) < 1/48, 1)
		latindex = which(flat %in% 1)
		flon = match(abs(lon - x[2]) < 1/48, 1)
		lonindex = which(flon %in% 1)
		start =c(lonindex, latindex, 1)
		count = c(1, 1, -1)

		# # #To plot an area:
		lat_range <- seq(latindex-10, latindex+10)
		lon_range <- seq(lonindex-10, lonindex+10)
		#image.plot(lon_range, lat_range,pet[lon_range,c(max(lat_range):min(lat_range)),6])
		image.plot(lon_range, lat_range,pet[lon_range,lat_range,6])
		abline(h=latindex)
		abline(v=lonindex)
		# #image.plot((pet[,c(4320:1),1]))
		



		# read in the full period of record using aggregated files
		data <- as.numeric(ncvar_get(nc, varid = var,start = start, count))
		pet_yearly[elevations %in% elevations[ll], ][years %in% years[yy],]$soil_pet = mean(data,na.rm=T)

	}

}

########################################################################################################################
#Climate water deficit
########################################################################################################################

#Added a loop to go over each site: 
#lat_long = data.frame(lat = c(46.874, 46.869, 46.878, 46.887, 46.888), lon = c(9.509, 9.490,9.494, 9.495, 9.489))

#This set of coordinates produces a gradient 
#lat_long = data.frame(lat = c(46.874, 46.869, 46.878, 46.887, 46.95), lon = c(9.509, 9.490,9.494, 9.495, 9.489))
#Best currently: 
lat_long = data.frame(lat = c(46.85, 46.874, 46.869, 46.878, 46.887, 46.95,46.91,46.89750), lon = c(9.509, 9.509, 9.490,9.494, 9.495, 9.489,9.490,9.46858))
elevations = c(600,1000,1400,1600,1800,2000,2100,2700)

def_yearly = (expand.grid(elevations, years))
def_yearly=data.frame(cbind(def_yearly,matrix(0,dim(def_yearly[1]),1)))
colnames(def_yearly) = c('elevations','years','soil_def')

#Added a loop to go over each year: 

for (yy in 1:3) {
	
	year = years[yy]

	# enter in variable you want to download see: http://thredds.northwestknowledge.net:8080/thredds/terraclimate_aggregated.html
	var="def"

	basefile = paste("soil_moisture/TerraClimate_", var, "_", year,".nc", sep="")
	nc = nc_open(basefile)
	lat =ncvar_get(nc, varid="lat")
	lon =ncvar_get(nc, varid="lon")
	#time =ncvar_get(nc, varid="time")
	def = ncvar_get(nc, varid="def")


	for( ll in 1:5) {

		# enter in longitude, latitude here
		x=c(unlist(lat_long[ll,]))



		# time_d = as.Date(time, format="%d", origin=as.Date( "1899-01-01" ) )
		# time_years = format(time_d, "%Y")
		# time_months = format(time_d, "%m")
		# time_ym = format(time_d, "%Y-%m")

		flat = match(abs(lat - x[1]) < 1/48, 1)
		latindex = which(flat %in% 1)
		flon = match(abs(lon - x[2]) < 1/48, 1)
		lonindex = which(flon %in% 1)
		start =c(lonindex, latindex, 1)
		count = c(1, 1, -1)

		# # #To plot an area:
		lat_range <- seq(latindex-10, latindex+10)
		lon_range <- seq(lonindex-10, lonindex+10)
		#image.plot(lon_range, lat_range,def[lon_range,c(max(lat_range):min(lat_range)),6])
		image.plot(lon_range, lat_range,def[lon_range,lat_range,6])
		abline(h=latindex)
		abline(v=lonindex)
		# #image.plot((def[,c(4320:1),1]))
		



		# read in the full period of record using aggregated files
		data <- as.numeric(ncvar_get(nc, varid = var,start = start, count))
		def_yearly[elevations %in% elevations[ll], ][years %in% years[yy],]$soil_def = mean(data,na.rm=T)

	}

}


#=============================================================================
#Look at the relationship between each of these and elevation per year
#Make data objects that can be used with the model fitting. These include
#a yearly per-elevation mean/min/max summary statistic for:
#	min temp
#	max temp
#	mean temp
#	max light
#	mean light
#	max gdd
#	growing season length
#	soil moisture deficit
#	potential evapotranspiration
#	each of these but for the growing season only, where the growing season
#	is determined as the width of time between when gdd start and then threshold
#	14 total new variables here. 
#=============================================================================

#This is the main data object, which will be used as a lookup table later to 
#match entries in the larger data set

#Note: This is not currently set up to actually include the soil moisture
#This happens in a separate step below. 

nyears = c("2015","2016","2017")
new_length = length(nyears)*length(elevations)
env_vars = c("an_min_temp", "an_max_temp", "an_mean_temp", 
			"an_max_light", "an_mean_light","max_gdd","gsl",
			"gs_min_temp", "gs_max_temp", "gs_mean_temp", 
			"gs_max_light", "gs_mean_light")
env_yearly = data.frame(expand.grid(elevations, nyears))
for(n in 1:length(env_vars)){
	env_yearly = cbind(env_yearly, matrix(0,dim(env_yearly)[1],1))
}
colnames(env_yearly) = c("elevations", "year", env_vars)

for (y in 1:length(nyears)){
	all_sub_year = subset(all_sub, format(all_sub[,1], "%Y") == nyears[y])
	all_sub_year2 = subset(all_sub2, format(all_sub2[,1], "%Y") == nyears[y])
	

	for ( ee in 1:length(elevations)){
		all_sub_yel = all_sub_year[all_sub_year$elevations == elevations[ee],]
		all_sub_yel2 = all_sub_year2[all_sub_year2$elevations == elevations[ee],]
		index1 = (y-1)*length(elevations)+ee
		if(dim(all_sub_yel)[1] == 0){
			env_yearly[index1,3:14] = NA
		}else{ 
			env_yearly[index1,3] = min(all_sub_yel$temp_min )
			env_yearly[index1,4] = max(all_sub_yel$temp_max )
			env_yearly[index1,5] = mean(all_sub_yel$temp_mean )
			env_yearly[index1,6] = max(all_sub_yel$light_max )
			env_yearly[index1,7] = mean(all_sub_yel$light_mean )
			#Growing season
			mgdd = max(all_sub_yel2$gdd_all )
			s_start = max(which(all_sub_yel2$gdd_all == 0))
			s_end = min(which(all_sub_yel2$gdd_all == mgdd))
			if(is.infinite(s_start) ) { 
				s_start = 1
				s_end = length(all_sub_yel2$gdd_all)
			}
			d1 = as.Date(as.character(  all_sub_yel2[s_start,1]), "%Y-%m-%d")
			d2 = as.Date(as.character(  all_sub_yel2[s_end,1]), "%Y-%m-%d")
			env_yearly[index1,8] = mgdd
			env_yearly[index1,9] = as.numeric(d2-d1)
			#Variables over growing season only 
			dd1 = which(all_sub_yel[,1] == d1)
			dd2 = which(all_sub_yel[,1] == d2)
			env_yearly[index1,10] = min(all_sub_yel$temp_min[dd1:dd2] )
			env_yearly[index1,11] = max(all_sub_yel$temp_max[dd1:dd2] )
			env_yearly[index1,12] = mean(all_sub_yel$temp_mean[dd1:dd2] )
			env_yearly[index1,13] = max(all_sub_yel$light_max[dd1:dd2] )
			env_yearly[index1,14] = mean(all_sub_yel$light_mean[dd1:dd2] )
		}
	}
}

#Add the soil moisture column, PET, and climate water deficit
env_yearly=cbind(env_yearly, moist_yearly$soil_moist)
env_yearly=cbind(env_yearly, pet_yearly$soil_pet)
env_yearly=cbind(env_yearly, def_yearly$soil_def)

colnames(env_yearly)[15] = "soil_moist"
colnames(env_yearly)[16] = "soil_pet"
colnames(env_yearly)[17] = "soil_def"
