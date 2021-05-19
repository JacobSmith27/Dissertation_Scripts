library(ncdf4)
library(fields)
library(dplyr)
library(tidyr)
library(lubridate)
library(raster)
library(zoo)
library(ggplot2)
library(gridExtra)

# Reading SSP126
Drivers_ssp126 = nc_open("Kenya_DRIVERS_OBS_MOHC_ssp126_climate_change_2001_2100.nc")

# Reading SSP245
Drivers_ssp245 = nc_open("Kenya_DRIVERS_OBS_MOHC_ssp245_climate_change_2001_2100.nc")

# Reading SSP370
Drivers_ssp370 = nc_open("Kenya_DRIVERS_OBS_MOHC_ssp370_climate_change_2001_2100.nc")

# Reading SSP585
Drivers_ssp585 = nc_open("Kenya_DRIVERS_OBS_MOHC_ssp585_climate_change_2001_2100.nc")

# Reading SSP434
Drivers_ssp434 = nc_open("Kenya_DRIVERS_OBS_MOHC_ssp434_climate_change_2001_2100.nc")

# Drivers data
# Driver data only has 3 dimensions not 4 like the rest as it doesnt have the 7 different estimates.

# Mean precipitation (kg.m-2.s-1)
pr_ssp126 = ncvar_get(Drivers_ssp126,"pr",count = c(19,27,1200)) 
pr_ssp245 = ncvar_get(Drivers_ssp245,"pr",count = c(19,27,1200))
pr_ssp370 = ncvar_get(Drivers_ssp370,"pr",count = c(19,27,1200))
pr_ssp585 = ncvar_get(Drivers_ssp585,"pr",count = c(19,27,1200))

# CO2 conc (ppm)
co2_ssp126 = ncvar_get(Drivers_ssp126,"co2",count = c(19,27,1200))
co2_ssp245 = ncvar_get(Drivers_ssp245,"co2",count = c(19,27,1200))
co2_ssp370 = ncvar_get(Drivers_ssp370,"co2",count = c(19,27,1200))
co2_ssp434 = ncvar_get(Drivers_ssp434,"co2",count = c(19,27,1200))
co2_ssp585 = ncvar_get(Drivers_ssp585,"co2",count = c(19,27,1200))

# Mean vapour pressure deficit (Pa)
vpd_ssp126 = ncvar_get(Drivers_ssp126,"vpd",count = c(19,27,1200))
vpd_ssp245 = ncvar_get(Drivers_ssp245,"vpd",count = c(19,27,1200))
vpd_ssp370 = ncvar_get(Drivers_ssp370,"vpd",count = c(19,27,1200))
vpd_ssp585 = ncvar_get(Drivers_ssp585,"vpd",count = c(19,27,1200))

# Mean daily minimum near surface air temp (C)
tas_min_ssp126 = ncvar_get(Drivers_ssp126,"tas_min",count = c(19,27,1200))
tas_min_ssp245 = ncvar_get(Drivers_ssp245,"tas_min",count = c(19,27,1200))
tas_min_ssp370 = ncvar_get(Drivers_ssp370,"tas_min",count = c(19,27,1200))
tas_min_ssp585 = ncvar_get(Drivers_ssp585,"tas_min",count = c(19,27,1200))

# Mean daily maximum near surface air temp (C)
tas_max_ssp126 = ncvar_get(Drivers_ssp126,"tas_max",count = c(19,27,1200))
tas_max_ssp245 = ncvar_get(Drivers_ssp245,"tas_max",count = c(19,27,1200))
tas_max_ssp370 = ncvar_get(Drivers_ssp370,"tas_max",count = c(19,27,1200))
tas_max_ssp585 = ncvar_get(Drivers_ssp585,"tas_max",count = c(19,27,1200))

# Mean downwelling short wave rad(MJ.m-2.d-1)
rsds_ssp126 = ncvar_get(Drivers_ssp126,"rsds",count = c(19,27,1200))
rsds_ssp245 = ncvar_get(Drivers_ssp245,"rsds",count = c(19,27,1200))
rsds_ssp370 = ncvar_get(Drivers_ssp370,"rsds",count = c(19,27,1200))
rsds_ssp585 = ncvar_get(Drivers_ssp585,"rsds",count = c(19,27,1200))

# Mean wind speed (m.s-1)
wsp_ssp126 = ncvar_get(Drivers_ssp126,"wsp",count = c(19,27,1200))
wsp_ssp245 = ncvar_get(Drivers_ssp245,"wsp",count = c(19,27,1200))
wsp_ssp370 = ncvar_get(Drivers_ssp370,"wsp",count = c(19,27,1200))
wsp_ssp585 = ncvar_get(Drivers_ssp585,"wsp",count = c(19,27,1200))

# LAI observed
laiobs_ssp126 = ncvar_get(Drivers_ssp126,"LAI_OBS",count = c(19,27,228))

# Get lat and lon
lon <- ncvar_get(cstock_ssp126, "lon")
lat <- ncvar_get(cstock_ssp126, "lat")
lon <- lon[1:24]
lat <- lat[1:25]

# Get quantile variable to identify 1-7 
quantiles = ncvar_get(cstock_ssp126,"quantile")
print(quantiles)
# 1 = 0.025
# 2 = 0.050
# 3 = 0.250
# 4 = 0.500
# 5 = 0.750
# 6 = 0.950
# 7 = 0.975

# When done with reading in variables remember to close the file
nc_close(Drivers_ssp126);nc_close(Drivers_ssp245);nc_close(Drivers_ssp370);nc_close(Drivers_ssp585);nc_close(Drivers_ssp434)

# Converting drivers into a dataframe



# Precipitation
par(mfrow=c(3,2))
image.plot(lon, lat, pr_ssp126[,,1],
main = "Dom in 2019 (gC.m-2)",
xlab = "Longitude",
ylab = "Latitude",
legend.line = 4,
zlim=c(0,60000),
col = viridis(200))
plot(Kenya_shp, add = T, border = "Black")
image.plot(pr_ssp126[,,1], zlim= c(0,0.0002),legend.lab='Pr_2001')
image.plot(pr_ssp126[,,1200], zlim= c(0,0.0002),legend.lab='Pr_126')
image.plot(pr_ssp245[,,1200], zlim= c(0,0.0002),legend.lab='Pr_245')
image.plot(pr_ssp370[,,1200], zlim= c(0,0.0002),legend.lab='Pr_370')
image.plot(pr_ssp585[,,1200], zlim= c(0,0.0002),legend.lab='Pr_585')


x <- (2001:2100)

# Converting values from monthly to annual to remove some noise
# How many time steps per year?
steps_per_year = 12

# Converting ssp126 co2 to annual
co2_ssp126_annual = array(NA, dim=c(dim(co2_ssp126)[1],dim(co2_ssp126)[2],dim(co2_ssp126)[3]/steps_per_year))
for (i in seq(1,dim(co2_ssp126)[1])) {
     for (j in seq(1,dim(co2_ssp126)[2])) {
          co2_ssp126_annual[i,j,] = rollapply(co2_ssp126[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# Converting ssp245 co2 to annual
co2_ssp245_annual = array(NA, dim=c(dim(co2_ssp245)[1],dim(co2_ssp245)[2],dim(co2_ssp245)[3]/steps_per_year))
for (i in seq(1,dim(co2_ssp245)[1])) {
     for (j in seq(1,dim(co2_ssp245)[2])) {
          co2_ssp245_annual[i,j,] = rollapply(co2_ssp245[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# Converting ssp370 co2 to annual
co2_ssp370_annual = array(NA, dim=c(dim(co2_ssp370)[1],dim(co2_ssp370)[2],dim(co2_ssp370)[3]/steps_per_year))
for (i in seq(1,dim(co2_ssp370)[1])) {
     for (j in seq(1,dim(co2_ssp370)[2])) {
          co2_ssp370_annual[i,j,] = rollapply(co2_ssp370[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# Converting ssp585 co2 to annual
co2_ssp585_annual = array(NA, dim=c(dim(co2_ssp585)[1],dim(co2_ssp585)[2],dim(co2_ssp585)[3]/steps_per_year))
for (i in seq(1,dim(co2_ssp585)[1])) {
     for (j in seq(1,dim(co2_ssp585)[2])) {
          co2_ssp585_annual[i,j,] = rollapply(co2_ssp585[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# CO2 conc plot
plot(x,apply(co2_ssp126_annual,3,mean, na.rm=TRUE), type="l", ylim=c(300,1200),xlab="Year", ylab="Atmospheric CO2 concentration (ppm)", col="goldenrod3")
lines(x,apply(co2_ssp245_annual,3,mean, na.rm=TRUE), type="l", col="darkgreen")
lines(x,apply(co2_ssp370_annual,3,mean, na.rm=TRUE), type="l", col ="blue")
lines(x,apply(co2_ssp585_annual,3,mean, na.rm=TRUE), type="l", col ="firebrick")
legend(x = "top",          # Position
       legend = c("SSP126", "SSP245", "SSP370", "SSP585"), 
	bty="n",
       col = c("goldenrod3", "darkgreen", "blue", "firebrick"), 
       lwd = 2)                


# Converting values from monthly to 5 years to remove some noise
# How many time steps per year?
steps_per_5years = 60

# Converting ssp126 pr to annual
pr_ssp126_annual = array(NA, dim=c(dim(pr_ssp126)[1],dim(pr_ssp126)[2],dim(pr_ssp126)[3]/steps_per_5years))
for (i in seq(1,dim(pr_ssp126)[1])) {
     for (j in seq(1,dim(pr_ssp126)[2])) {
          pr_ssp126_annual[i,j,] = rollapply(pr_ssp126[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp245 pr to annual
pr_ssp245_annual = array(NA, dim=c(dim(pr_ssp245)[1],dim(pr_ssp245)[2],dim(pr_ssp245)[3]/steps_per_5years))
for (i in seq(1,dim(pr_ssp245)[1])) {
     for (j in seq(1,dim(pr_ssp245)[2])) {
          pr_ssp245_annual[i,j,] = rollapply(pr_ssp245[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp370 pr to annual
pr_ssp370_annual = array(NA, dim=c(dim(pr_ssp370)[1],dim(pr_ssp370)[2],dim(pr_ssp370)[3]/steps_per_5years))
for (i in seq(1,dim(pr_ssp370)[1])) {
     for (j in seq(1,dim(pr_ssp370)[2])) {
          pr_ssp370_annual[i,j,] = rollapply(pr_ssp370[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp585 pr to annual
pr_ssp585_annual = array(NA, dim=c(dim(pr_ssp585)[1],dim(pr_ssp585)[2],dim(pr_ssp585)[3]/steps_per_5years))
for (i in seq(1,dim(pr_ssp585)[1])) {
     for (j in seq(1,dim(pr_ssp585)[2])) {
          pr_ssp585_annual[i,j,] = rollapply(pr_ssp585[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# VPD
# Converting ssp126 vpd to annual
vpd_ssp126_annual = array(NA, dim=c(dim(vpd_ssp126)[1],dim(vpd_ssp126)[2],dim(vpd_ssp126)[3]/steps_per_5years))
for (i in seq(1,dim(vpd_ssp126)[1])) {
     for (j in seq(1,dim(vpd_ssp126)[2])) {
          vpd_ssp126_annual[i,j,] = rollapply(vpd_ssp126[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp245 vpd to annual
vpd_ssp245_annual = array(NA, dim=c(dim(vpd_ssp245)[1],dim(vpd_ssp245)[2],dim(vpd_ssp245)[3]/steps_per_5years))
for (i in seq(1,dim(vpd_ssp245)[1])) {
     for (j in seq(1,dim(vpd_ssp245)[2])) {
          vpd_ssp245_annual[i,j,] = rollapply(vpd_ssp245[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp370 vpd to annual
vpd_ssp370_annual = array(NA, dim=c(dim(vpd_ssp370)[1],dim(vpd_ssp370)[2],dim(vpd_ssp370)[3]/steps_per_5years))
for (i in seq(1,dim(vpd_ssp370)[1])) {
     for (j in seq(1,dim(vpd_ssp370)[2])) {
          vpd_ssp370_annual[i,j,] = rollapply(vpd_ssp370[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp585 vpd to annual
vpd_ssp585_annual = array(NA, dim=c(dim(vpd_ssp585)[1],dim(vpd_ssp585)[2],dim(vpd_ssp585)[3]/steps_per_5years))
for (i in seq(1,dim(vpd_ssp585)[1])) {
     for (j in seq(1,dim(vpd_ssp585)[2])) {
          vpd_ssp585_annual[i,j,] = rollapply(vpd_ssp585[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# rsds
# Converting ssp126 rsds to annual
rsds_ssp126_annual = array(NA, dim=c(dim(rsds_ssp126)[1],dim(rsds_ssp126)[2],dim(rsds_ssp126)[3]/steps_per_5years))
for (i in seq(1,dim(rsds_ssp126)[1])) {
     for (j in seq(1,dim(rsds_ssp126)[2])) {
          rsds_ssp126_annual[i,j,] = rollapply(rsds_ssp126[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp245 rsds to annual
rsds_ssp245_annual = array(NA, dim=c(dim(rsds_ssp245)[1],dim(rsds_ssp245)[2],dim(rsds_ssp245)[3]/steps_per_5years))
for (i in seq(1,dim(rsds_ssp245)[1])) {
     for (j in seq(1,dim(rsds_ssp245)[2])) {
          rsds_ssp245_annual[i,j,] = rollapply(rsds_ssp245[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp370 rsds to annual
rsds_ssp370_annual = array(NA, dim=c(dim(rsds_ssp370)[1],dim(rsds_ssp370)[2],dim(rsds_ssp370)[3]/steps_per_5years))
for (i in seq(1,dim(rsds_ssp370)[1])) {
     for (j in seq(1,dim(rsds_ssp370)[2])) {
          rsds_ssp370_annual[i,j,] = rollapply(rsds_ssp370[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp585 rsds to annual
rsds_ssp585_annual = array(NA, dim=c(dim(rsds_ssp585)[1],dim(rsds_ssp585)[2],dim(rsds_ssp585)[3]/steps_per_5years))
for (i in seq(1,dim(rsds_ssp585)[1])) {
     for (j in seq(1,dim(rsds_ssp585)[2])) {
          rsds_ssp585_annual[i,j,] = rollapply(rsds_ssp585[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# tas_min
# Converting ssp126 tas_min to annual
tas_min_ssp126_annual = array(NA, dim=c(dim(tas_min_ssp126)[1],dim(tas_min_ssp126)[2],dim(tas_min_ssp126)[3]/steps_per_5years))
for (i in seq(1,dim(tas_min_ssp126)[1])) {
     for (j in seq(1,dim(tas_min_ssp126)[2])) {
          tas_min_ssp126_annual[i,j,] = rollapply(tas_min_ssp126[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp245 tas_min to annual
tas_min_ssp245_annual = array(NA, dim=c(dim(tas_min_ssp245)[1],dim(tas_min_ssp245)[2],dim(tas_min_ssp245)[3]/steps_per_5years))
for (i in seq(1,dim(tas_min_ssp245)[1])) {
     for (j in seq(1,dim(tas_min_ssp245)[2])) {
          tas_min_ssp245_annual[i,j,] = rollapply(tas_min_ssp245[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp370 tas_min to annual
tas_min_ssp370_annual = array(NA, dim=c(dim(tas_min_ssp370)[1],dim(tas_min_ssp370)[2],dim(tas_min_ssp370)[3]/steps_per_5years))
for (i in seq(1,dim(tas_min_ssp370)[1])) {
     for (j in seq(1,dim(tas_min_ssp370)[2])) {
          tas_min_ssp370_annual[i,j,] = rollapply(tas_min_ssp370[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp585 tas_min to annual
tas_min_ssp585_annual = array(NA, dim=c(dim(tas_min_ssp585)[1],dim(tas_min_ssp585)[2],dim(tas_min_ssp585)[3]/steps_per_5years))
for (i in seq(1,dim(tas_min_ssp585)[1])) {
     for (j in seq(1,dim(tas_min_ssp585)[2])) {
          tas_min_ssp585_annual[i,j,] = rollapply(tas_min_ssp585[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# tas_max
# Converting ssp126 tas_max to annual
tas_max_ssp126_annual = array(NA, dim=c(dim(tas_max_ssp126)[1],dim(tas_max_ssp126)[2],dim(tas_max_ssp126)[3]/steps_per_5years))
for (i in seq(1,dim(tas_max_ssp126)[1])) {
     for (j in seq(1,dim(tas_max_ssp126)[2])) {
          tas_max_ssp126_annual[i,j,] = rollapply(tas_max_ssp126[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp245 tas_max to annual
tas_max_ssp245_annual = array(NA, dim=c(dim(tas_max_ssp245)[1],dim(tas_max_ssp245)[2],dim(tas_max_ssp245)[3]/steps_per_5years))
for (i in seq(1,dim(tas_max_ssp245)[1])) {
     for (j in seq(1,dim(tas_max_ssp245)[2])) {
          tas_max_ssp245_annual[i,j,] = rollapply(tas_max_ssp245[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp370 tas_max to annual
tas_max_ssp370_annual = array(NA, dim=c(dim(tas_max_ssp370)[1],dim(tas_max_ssp370)[2],dim(tas_max_ssp370)[3]/steps_per_5years))
for (i in seq(1,dim(tas_max_ssp370)[1])) {
     for (j in seq(1,dim(tas_max_ssp370)[2])) {
          tas_max_ssp370_annual[i,j,] = rollapply(tas_max_ssp370[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp585 tas_max to annual
tas_max_ssp585_annual = array(NA, dim=c(dim(tas_max_ssp585)[1],dim(tas_max_ssp585)[2],dim(tas_max_ssp585)[3]/steps_per_5years))
for (i in seq(1,dim(tas_max_ssp585)[1])) {
     for (j in seq(1,dim(tas_max_ssp585)[2])) {
          tas_max_ssp585_annual[i,j,] = rollapply(tas_max_ssp585[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# wsp
# Converting ssp126 wsp to annual
wsp_ssp126_annual = array(NA, dim=c(dim(wsp_ssp126)[1],dim(wsp_ssp126)[2],dim(wsp_ssp126)[3]/steps_per_5years))
for (i in seq(1,dim(wsp_ssp126)[1])) {
     for (j in seq(1,dim(wsp_ssp126)[2])) {
          wsp_ssp126_annual[i,j,] = rollapply(wsp_ssp126[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp245 wsp to annual
wsp_ssp245_annual = array(NA, dim=c(dim(wsp_ssp245)[1],dim(wsp_ssp245)[2],dim(wsp_ssp245)[3]/steps_per_5years))
for (i in seq(1,dim(wsp_ssp245)[1])) {
     for (j in seq(1,dim(wsp_ssp245)[2])) {
          wsp_ssp245_annual[i,j,] = rollapply(wsp_ssp245[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp370 wsp to annual
wsp_ssp370_annual = array(NA, dim=c(dim(wsp_ssp370)[1],dim(wsp_ssp370)[2],dim(wsp_ssp370)[3]/steps_per_5years))
for (i in seq(1,dim(wsp_ssp370)[1])) {
     for (j in seq(1,dim(wsp_ssp370)[2])) {
          wsp_ssp370_annual[i,j,] = rollapply(wsp_ssp370[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

# Converting ssp585 wsp to annual
wsp_ssp585_annual = array(NA, dim=c(dim(wsp_ssp585)[1],dim(wsp_ssp585)[2],dim(wsp_ssp585)[3]/steps_per_5years))
for (i in seq(1,dim(wsp_ssp585)[1])) {
     for (j in seq(1,dim(wsp_ssp585)[2])) {
          wsp_ssp585_annual[i,j,] = rollapply(wsp_ssp585[i,j,], by = steps_per_5years, width = steps_per_5years, FUN=mean, na.rm=TRUE) } } 

x1 <- seq(2001, 2100, by = 5)

# Unit conversion
# Creating functions for conversion
pr_conv <- function(x){ # Converts each pixel from kg.m-2.s-1 to mm.d-1
	x*86400
	}

pr_ssp126_annual <- apply(pr_ssp126_annual, c(1,2),pr_conv)
pr_ssp245_annual <- apply(pr_ssp245_annual, c(1,2),pr_conv)
pr_ssp370_annual <- apply(pr_ssp370_annual, c(1,2),pr_conv)
pr_ssp585_annual <- apply(pr_ssp585_annual, c(1,2),pr_conv)

par(mfrow=c(3,2))
# Precipitation plot
plot(x1,apply(pr_ssp126_annual,1,mean, na.rm=TRUE), type="l",ylim=c(2.5,6.5),xlab="Year", ylab="Precipitation (mm.d-1)",cex.lab=1.3,  col="goldenrod3")
lines(x1,apply(pr_ssp245_annual,1,mean, na.rm=TRUE), type="l", col="darkgreen")
lines(x1,apply(pr_ssp370_annual,1,mean, na.rm=TRUE), type="l", col ="blue")
lines(x1,apply(pr_ssp585_annual,1,mean, na.rm=TRUE), type="l", col ="firebrick")

# rsds plot
plot(x1,apply(rsds_ssp126_annual,3,mean, na.rm=TRUE), type="l",ylim=c(18.5,20),xlab="Year", ylab="Downwelling SWR (MJ.m-1.d-1)",cex.lab=1.3,  col="goldenrod3")
lines(x1,apply(rsds_ssp245_annual,3,mean, na.rm=TRUE), type="l", col="darkgreen")
lines(x1,apply(rsds_ssp370_annual,3,mean, na.rm=TRUE), type="l", col ="blue")
lines(x1,apply(rsds_ssp585_annual,3,mean, na.rm=TRUE), type="l", col ="firebrick")

# vpd plot
plot(x1,apply(vpd_ssp126_annual,3,mean, na.rm=TRUE), type="l", ylim=c(800,1400), xlab="Year", ylab="Vapour pressure deficit (Pa)",cex.lab=1.3,  col="goldenrod3")
lines(x1,apply(vpd_ssp245_annual,3,mean, na.rm=TRUE), type="l", col="darkgreen")
lines(x1,apply(vpd_ssp370_annual,3,mean, na.rm=TRUE), type="l", col ="blue")
lines(x1,apply(vpd_ssp585_annual,3,mean, na.rm=TRUE), type="l", col ="firebrick")

# wsp plot
plot(x1,apply(wsp_ssp126_annual,3,mean, na.rm=TRUE), type="l",ylim=c(1.8,2.8),xlab="Year", ylab="Wind speed(m.s-1)",cex.lab=1.3,  col="goldenrod3")
lines(x1,apply(wsp_ssp245_annual,3,mean, na.rm=TRUE), type="l", col="darkgreen")
lines(x1,apply(wsp_ssp370_annual,3,mean, na.rm=TRUE), type="l", col ="blue")
lines(x1,apply(wsp_ssp585_annual,3,mean, na.rm=TRUE), type="l", col ="firebrick")

# tas_min plot
plot(x1,apply(tas_min_ssp126_annual,3,mean, na.rm=TRUE), type="l", ylim=c(13,21), xlab="Year", ylab="Minimum near surface air temperature (C)", cex.lab=1.3, col="goldenrod3")
lines(x1,apply(tas_min_ssp245_annual,3,mean, na.rm=TRUE), type="l", col="darkgreen")
lines(x1,apply(tas_min_ssp370_annual,3,mean, na.rm=TRUE), type="l", col ="blue")
lines(x1,apply(tas_min_ssp585_annual,3,mean, na.rm=TRUE), type="l", col ="firebrick")

# tas_max plot
plot(x1,apply(tas_max_ssp126_annual,3,mean, na.rm=TRUE), type="l", ylim=c(26.5,34), xlab="Year", ylab="Maximum near surface air temperature (C)", cex.lab=1.3, col="goldenrod3")
lines(x1,apply(tas_max_ssp245_annual,3,mean, na.rm=TRUE), type="l", col="darkgreen")
lines(x1,apply(tas_max_ssp370_annual,3,mean, na.rm=TRUE), type="l", col ="blue")
lines(x1,apply(tas_max_ssp585_annual,3,mean, na.rm=TRUE), type="l", col ="firebrick")

