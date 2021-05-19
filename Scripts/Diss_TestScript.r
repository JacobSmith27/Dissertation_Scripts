#Load libraries
library(ncdf4)
library(fields)
library(dplyr)
library(lubridate)
library(ggplot2)
library(raster)
library(zoo)
library(kableExtra)
# install.packages("kableExtra")

# Read in SSP126
cstock_ssp126 = nc_open("Kenya_CSTOCK_MOHC_ssp126_climate_change_2001_2100.nc")
cflux_ssp126 = nc_open("Kenya_CFLUX_MOHC_ssp126_climate_change_2001_2100.nc")
Drivers_ssp126 = nc_open("Kenya_DRIVERS_OBS_MOHC_ssp126_climate_change_2001_2100.nc")

# Read in SSP245
cstock_ssp245 = nc_open("Kenya_CSTOCK_MOHC_ssp245_climate_change_2001_2100.nc")
cflux_ssp245 = nc_open("Kenya_CFLUX_MOHC_ssp245_climate_change_2001_2100.nc")
Drivers_ssp245 = nc_open("Kenya_DRIVERS_OBS_MOHC_ssp245_climate_change_2001_2100.nc")

# Read in SSP370
cstock_ssp370 = nc_open("Kenya_CSTOCK_MOHC_ssp370_climate_change_2001_2100.nc")
cflux_ssp370 = nc_open("Kenya_CFLUX_MOHC_ssp370_climate_change_2001_2100.nc")
Drivers_ssp370 = nc_open("Kenya_DRIVERS_OBS_MOHC_ssp370_climate_change_2001_2100.nc")

# Read in SSP434
Drivers_ssp434 = nc_open("Kenya_DRIVERS_OBS_MOHC_ssp434_climate_change_2001_2100.nc")

# Read in SSP585
cstock_ssp585 = nc_open("Kenya_CSTOCK_MOHC_ssp585_climate_change_2001_2100.nc")
cflux_ssp585 = nc_open("Kenya_CFLUX_MOHC_ssp585_climate_change_2001_2100.nc")
Drivers_ssp585 = nc_open("Kenya_DRIVERS_OBS_MOHC_ssp585_climate_change_2001_2100.nc")

# Reading in from STOCKS
lai_ssp126 = ncvar_get(cstock_ssp126,"lai_ensemble",count = c(19,27,7,1200))
lai_ssp245 = ncvar_get(cstock_ssp245,"lai_ensemble",count = c(19,27,7,1200))
lai_ssp370 = ncvar_get(cstock_ssp370,"lai_ensemble",count = c(19,27,7,1200))
lai_ssp585 = ncvar_get(cstock_ssp585,"lai_ensemble",count = c(19,27,7,1200))

bio_ssp126 = ncvar_get(cstock_ssp126,"cVeg_ensemble",count = c(19,27,7,1200))
bio_ssp245 = ncvar_get(cstock_ssp245,"cVeg_ensemble",count = c(19,27,7,1200))
bio_ssp370 = ncvar_get(cstock_ssp370,"cVeg_ensemble",count = c(19,27,7,1200))
bio_ssp585 = ncvar_get(cstock_ssp585,"cVeg_ensemble",count = c(19,27,7,1200))

dom_ssp126 = ncvar_get(cstock_ssp126,"cDOM_ensemble",count = c(19,27,7,1200))
dom_ssp245 = ncvar_get(cstock_ssp245,"cDOM_ensemble",count = c(19,27,7,1200))
dom_ssp370 = ncvar_get(cstock_ssp370,"cDOM_ensemble",count = c(19,27,7,1200))
dom_ssp585 = ncvar_get(cstock_ssp585,"cDOM_ensemble",count = c(19,27,7,1200))

cTotal_ssp126 = ncvar_get(cstock_ssp126,"cTotal",count = c(19,27,7,1200))
cTotal_ssp245 = ncvar_get(cstock_ssp245,"cTotal",count = c(19,27,7,1200))
cTotal_ssp370 = ncvar_get(cstock_ssp370,"cTotal",count = c(19,27,7,1200))
cTotal_ssp585 = ncvar_get(cstock_ssp585,"cTotal",count = c(19,27,7,1200))

# Reading in from FLUX
# Gross primary productivity (gC.m-2.d-1)
gpp_ssp126 = ncvar_get(cflux_ssp126,"gpp_ensemble",count = c(19,27,7,1200))
gpp_ssp245 = ncvar_get(cflux_ssp245,"gpp_ensemble",count = c(19,27,7,1200))
gpp_ssp370 = ncvar_get(cflux_ssp370,"gpp_ensemble",count = c(19,27,7,1200))
gpp_ssp585 = ncvar_get(cflux_ssp585,"gpp_ensemble",count = c(19,27,7,1200))

# Autotrophic respiration (gC.m-2.d-1)
ra_ssp126 = ncvar_get(cflux_ssp126,"ra_ensemble",count = c(19,27,7,1200))
ra_ssp245 = ncvar_get(cflux_ssp245,"ra_ensemble",count = c(19,27,7,1200))
ra_ssp370 = ncvar_get(cflux_ssp370,"ra_ensemble",count = c(19,27,7,1200))
ra_ssp585 = ncvar_get(cflux_ssp585,"ra_ensemble",count = c(19,27,7,1200))

# Net primary productivity (gC.m-2.d-1)
npp_ssp126 = ncvar_get(cflux_ssp126,"npp_ensemble",count = c(19,27,7,1200))
npp_ssp245 = ncvar_get(cflux_ssp245,"npp_ensemble",count = c(19,27,7,1200))
npp_ssp370 = ncvar_get(cflux_ssp370,"npp_ensemble",count = c(19,27,7,1200))
npp_ssp585 = ncvar_get(cflux_ssp585,"npp_ensemble",count = c(19,27,7,1200))

# Heterotrophic respiration (gC.m-2.d-1)
rh_ssp126 = ncvar_get(cflux_ssp126,"rh_ensemble",count = c(19,27,7,1200))
rh_ssp245 = ncvar_get(cflux_ssp245,"rh_ensemble",count = c(19,27,7,1200))
rh_ssp370 = ncvar_get(cflux_ssp370,"rh_ensemble",count = c(19,27,7,1200))
rh_ssp585 = ncvar_get(cflux_ssp585,"rh_ensemble",count = c(19,27,7,1200))

# Net ecosystem exchange (gC.m-2.d-1)
nee_ssp126 = ncvar_get(cflux_ssp126,"nee_ensemble",count = c(19,27,7,1200))
nee_ssp245 = ncvar_get(cflux_ssp245,"nee_ensemble",count = c(19,27,7,1200))
nee_ssp370 = ncvar_get(cflux_ssp370,"nee_ensemble",count = c(19,27,7,1200))
nee_ssp585 = ncvar_get(cflux_ssp585,"nee_ensemble",count = c(19,27,7,1200))

# Fire losses (gC.m-2.d-1)
fire_ssp126 = ncvar_get(cflux_ssp126,"fFire_ensemble",count = c(19,27,7,1200))
fire_ssp245 = ncvar_get(cflux_ssp245,"fFire_ensemble",count = c(19,27,7,1200))
fire_ssp370 = ncvar_get(cflux_ssp370,"fFire_ensemble",count = c(19,27,7,1200))
fire_ssp585 = ncvar_get(cflux_ssp585,"fFire_ensemble",count = c(19,27,7,1200))

# Net biome exchange (NEE+Fire) (gC.m-2.d-1)
nbe_ssp126 = ncvar_get(cflux_ssp126,"nbe_ensemble",count = c(19,27,7,1200))
nbe_ssp245 = ncvar_get(cflux_ssp245,"nbe_ensemble",count = c(19,27,7,1200))
nbe_ssp370 = ncvar_get(cflux_ssp370,"nbe_ensemble",count = c(19,27,7,1200))
nbe_ssp585 = ncvar_get(cflux_ssp585,"nbe_ensemble",count = c(19,27,7,1200))

# Forest harvet (gC.m-2.d-1)
fLuc_ssp126 = ncvar_get(cflux_ssp126,"fLuc_ensemble",count = c(19,27,7,1200))
fLuc_ssp245 = ncvar_get(cflux_ssp245,"fLuc_ensemble",count = c(19,27,7,1200))
fLuc_ssp370 = ncvar_get(cflux_ssp370,"fLuc_ensemble",count = c(19,27,7,1200))
fLuc_ssp585 = ncvar_get(cflux_ssp585,"fLuc_ensemble",count = c(19,27,7,1200))

# Internal:Ambiant CO2 ratio
CiCa_ssp126 = ncvar_get(cflux_ssp126,"CiCa_Ensemble",count = c(19,27,7,1200))
CiCa_ssp245 = ncvar_get(cflux_ssp245,"CiCa_Ensemble",count = c(19,27,7,1200))
CiCa_ssp370 = ncvar_get(cflux_ssp370,"CiCa_Ensemble",count = c(19,27,7,1200))
CiCa_ssp585 = ncvar_get(cflux_ssp585,"CiCa_Ensemble",count = c(19,27,7,1200))

# Net biome productivity (-NEE-fire-fluv)(gC.m-2.d-1)
nbp_ssp126 = ncvar_get(cflux_ssp126,"nbp_ensemble",count = c(19,27,7,1200))
nbp_ssp245 = ncvar_get(cflux_ssp245,"nbp_ensemble",count = c(19,27,7,1200))
nbp_ssp370 = ncvar_get(cflux_ssp370,"nbp_ensemble",count = c(19,27,7,1200))
nbp_ssp585 = ncvar_get(cflux_ssp585,"nbp_ensemble",count = c(19,27,7,1200))

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
nc_close(cflux_ssp126) ; nc_close(cstock_ssp126); nc_close(Drivers_ssp126)
nc_close(cflux_ssp245) ; nc_close(cstock_ssp245); nc_close(Drivers_ssp245)
nc_close(cflux_ssp370) ; nc_close(cstock_ssp370); nc_close(Drivers_ssp370)
nc_close(cflux_ssp585) ; nc_close(cstock_ssp585); nc_close(Drivers_ssp585)
nc_close(Drivers_ssp434)

# Plotting Stocks data

# Plot map lai
par(mfrow=c(3,2))
image.plot(lai_ssp126[,,4,1], zlim=c(0,20), legend.lab='LAI_2001')
image.plot(lai_ssp126[,,4,1200], zlim=c(0,20),legend.lab='LAI_126')
image.plot(lai_ssp245[,,4,1200], zlim=c(0,20),legend.lab='LAI_245')
image.plot(lai_ssp370[,,4,1200], zlim=c(0,20),legend.lab='LAI_370')
image.plot(lai_ssp585[,,4,1200], zlim=c(0,20),legend.lab='LAI_585')
# 370 much lower than the others, why?

# Plotting better maps

# Slicing through to just get the medians (4)
lai_ssp126_4 <- lai_ssp126[,,4,]
lai_ssp245_4 <- lai_ssp245[,,4,]
lai_ssp370_4 <- lai_ssp370[,,4,]
lai_ssp585_4 <- lai_ssp585[,,4,]

# Slicing through to just get 0.250(4)
lai_ssp126_1 <- lai_ssp126[,,1,]
lai_ssp245_1 <- lai_ssp245[,,1,]
lai_ssp370_1 <- lai_ssp370[,,1,]
lai_ssp585_1 <- lai_ssp585[,,1,]

# Slicing through to just get 0.750 (4)
lai_ssp126_7 <- lai_ssp126[,,7,]
lai_ssp245_7 <- lai_ssp245[,,7,]
lai_ssp370_7 <- lai_ssp370[,,7,]
lai_ssp585_7 <- lai_ssp585[,,7,]

# Plotting maps with more detail
# Slicing to get just the required time references 
lai_ssp126_4_1 <- lai_ssp126_4[,,1]
lai_ssp126_4_2100 <- lai_ssp126_4[,,1200]
lai_ssp245_4_2100 <- lai_ssp245_4[,,1200]
lai_ssp370_4_2100 <- lai_ssp370_4[,,1200]
lai_ssp585_4_2100 <- lai_ssp585_4[,,1200]

# Conditional replacement of LAI > 10
# Use whcih() to find out where LAIs are more than 10. 
which(lai_ssp126_4_1 >= 10, arr.ind=TRUE)
 
par(mfrow=c(2,3))
image.plot(lon, lat, lai_ssp126_4_1,
main = "LAI in 2001",
xlab = "Longitude",
ylab = "Latitude",
legend.lab = "Leaf Area Index",
legend.line = 2.5,
zlim=c(0,10))

image.plot(lon, lat, lai_ssp126_4_2100,
main = "LAI in 2100 with SSP126",
xlab = "Longitude",
ylab = "Latitude",
legend.lab = "Leaf Area Index",
legend.line = 2.5, 
zlim=c(0,10))

image.plot(lon, lat, lai_ssp245_4_2100,
main = "LAI in 2100 with SSP245",
xlab = "Longitude",
ylab = "Latitude",
legend.lab = "Leaf Area Index",
legend.line = 2.5, 
zlim=c(0,10))

image.plot(lon, lat, lai_ssp370_4_2100,
main = "LAI in 2100 with SSP370",
xlab = "Longitude",
ylab = "Latitude",
legend.lab = "Leaf Area Index",
legend.line = 2.5, 
zlim=c(0,10))

image.plot(lon, lat, lai_ssp585_4_2100,
main = "LAI in 2100 with SSP585",
xlab = "Longitude",
ylab = "Latitude",
legend.lab = "Leaf Area Index",
legend.line = 2.5, 
zlim=c(0,10))

# PLotting change in LAI over whole time period
image.plot(apply(lai_ssp126[,,4,],3, range, na.rm=TRUE)) # did not work

par(mfrow=c(2,2))
plot(apply(lai_ssp126[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,14))
lines(apply(lai_ssp126[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(lai_ssp126[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(lai_ssp245[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,14))
lines(apply(lai_ssp245[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(lai_ssp245[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(lai_ssp370[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,14))
lines(apply(lai_ssp370[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(lai_ssp370[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(lai_ssp585[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,14))
lines(apply(lai_ssp585[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(lai_ssp585[,,7,],3,mean, na.rm=TRUE), col="red")

# Try to change LAI into anual time steps not monthly by averaging

lai_ssp126_annual <- apply(lai_ssp126[,,4,],3,mean, na.rm=TRUE)

# New solution using raster package and zoo package. DIDNT WORK

stack_lai_126 <- stack(cstock_ssp126[19,27,4,1200], varname ="lai_ensemble")
lai_ssp126_annual <- as.Date(lai_ssp126$dim$time$vals/12, origin='2001-01-01')

# Solution from LUKE using zoo package
# Assuming LAI is in an array with dimensions (long,lat,time)
# In this example LAI is stored in a variable called lai_ssp126_4

# How many time steps per year?
steps_per_year = 12

# create new object which will contain the aggregated to mean annual rather than monthly
lai_ssp126_4_annual = array(NA, dim=c(dim(lai_ssp126_4)[1],dim(lai_ssp126_4)[2],dim(lai_ssp126_4)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp126_4)[1])) {
     for (j in seq(1,dim(lai_ssp126_4)[2])) {
          lai_ssp126_annual[i,j,] = rollapply(lai_ssp126_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

lai_ssp126_1_annual = array(NA, dim=c(dim(lai_ssp126_1)[1],dim(lai_ssp126_4)[2],dim(lai_ssp126_1)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp126_1)[1])) {
     for (j in seq(1,dim(lai_ssp126_1)[2])) {
          lai_ssp126_annual[i,j,] = rollapply(lai_ssp126_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

lai_ssp126_1_annual = array(NA, dim=c(dim(lai_ssp126_7)[1],dim(lai_ssp126_7)[2],dim(lai_ssp126_7)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp126_7)[1])) {
     for (j in seq(1,dim(lai_ssp126_7)[2])) {
          lai_ssp126_annual[i,j,] = rollapply(lai_ssp126_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# Plot the annual averages.
plot(apply(lai_ssp126_4_annual,3,mean, na.rm=TRUE), type="l", ylim=c(0,14), xlab="Year since 2000", ylab="Mean annual LAI values under the ssp126")
lines(apply(lai_ssp126_1_annual,3,mean, na.rm=TRUE), col="blue")
lines(apply(lai_ssp126_1_annual,3,mean, na.rm=TRUE), col="red")

# Plot map veg Carbon stock

range(bio_ssp585, na.rm=TRUE)
# Range is 34889

par(mfrow=c(3,2))
image.plot(bio_ssp126[,,4,1], zlim=c(0,25000),legend.lab='Cveg_2001')
image.plot(bio_ssp126[,,4,1200], zlim=c(0,25000),legend.lab='Cveg_126')
image.plot(bio_ssp245[,,4,1200], zlim=c(0,25000),legend.lab='Cveg_245')
image.plot(bio_ssp370[,,4,1200], zlim=c(0,25000),legend.lab='Cveg_370')
image.plot(bio_ssp585[,,4,1200], zlim=c(0,25000),legend.lab='Cveg_585')
# Some pixels could have some interesting dynamics and could look at how these evolve over time even at a regional level
# Pixels up by turkana area have excesively high veg carbon stocks ~40,000 gC.m-2. around 400 tonnes per ha.
# This is almost 4 times the carbon density of the most carbon dense system at current in the country.
# What fileting could be don? This is the middle estimate?
# could i filter using most desnse estimates that mat has from paper? this is the denssest in the country at current so surely things would need to be less dense?

# Mapping the most conservative posibility
par(mfrow=c(3,2))
image.plot(bio_ssp126[,,1,1], zlim=c(0,13000))
image.plot(bio_ssp126[,,1,1200], zlim=c(0,13000))
image.plot(bio_ssp245[,,1,1200], zlim=c(0,13000))
image.plot(bio_ssp370[,,1,1200], zlim=c(0,13000))
image.plot(bio_ssp585[,,1,1200], zlim=c(0,13000))

# Making better maps
# Slicing through to just get the medians (4)
bio_ssp126_4 <- bio_ssp126[,,4,]
bio_ssp245_4 <- bio_ssp245[,,4,]
bio_ssp370_4 <- bio_ssp370[,,4,]
bio_ssp585_4 <- bio_ssp585[,,4,]

# Slicing to get just the required time references 
bio_ssp126_4_1 <- bio_ssp126_4[,,1]
bio_ssp126_4_1200 <- bio_ssp126_4[,,1200]
bio_ssp245_4_1200 <- bio_ssp245_4[,,1200]
bio_ssp370_4_1200 <- bio_ssp370_4[,,1200]
bio_ssp585_4_1200 <- bio_ssp585_4[,,1200]

par(mfrow=c(2,3))
image.plot(lon, lat, bio_ssp126_4_1,
main = "Carbon in live biomass in 2001",
xlab = "Longitude",
ylab = "Latitude",
legend.lab = "Carbon in live biomass (gC.m-2)",
legend.line = 4,
zlim=c(0,13000))

image.plot(lon, lat, bio_ssp126_4_1200,
main = "Carbon in live biomass in 2100 with SSP126",
xlab = "Longitude",
ylab = "Latitude",
legend.lab = "Carbon in live biomass (gC.m-2)",
legend.line = 4,
zlim=c(0,13000))

image.plot(lon, lat, bio_ssp245_4_1200,
main = "Carbon in live biomass in 2100 with SSP245",
xlab = "Longitude",
ylab = "Latitude",
legend.lab = "Carbon in live biomass (gC.m-2)",
legend.line = 4,
zlim=c(0,13000))

image.plot(lon, lat, bio_ssp370_4_1200,
main = "Carbon in live biomass in 2100 with SSP370",
xlab = "Longitude",
ylab = "Latitude",
legend.lab = "Carbon in live biomass (gC.m-2)",
legend.line = 4,
zlim=c(0,13000))

image.plot(lon, lat, bio_ssp585_4_1200,
main = "Carbon in live biomass in 2100 with SSP585",
xlab = "Longitude",
ylab = "Latitude",
legend.lab = "Carbon in live biomass (gC.m-2)",
legend.line = 4,
zlim=c(0,13000))


par(mfrow=c(2,2))
plot(apply(bio_ssp126[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,9000))
lines(apply(bio_ssp126[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(bio_ssp126[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(bio_ssp245[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,9000))
lines(apply(bio_ssp245[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(bio_ssp245[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(bio_ssp370[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,9000))
lines(apply(bio_ssp370[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(bio_ssp370[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(bio_ssp585[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,9000))
lines(apply(bio_ssp585[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(bio_ssp585[,,7,],3,mean, na.rm=TRUE), col="red")
# An increase in 370 at 600 time step then decrease. peak isnt at final time step.

# Plot map Total Carbon stock
par(mfrow=c(3,2))
image.plot(cTotal_ssp126[,,4,1], zlim=c(0,120000),legend.lab='Ctot_2001')
image.plot(cTotal_ssp126[,,4,1200], zlim=c(0,120000),legend.lab='Ctot_126')
image.plot(cTotal_ssp245[,,4,1200], zlim=c(0,120000),legend.lab='Ctot_245')
image.plot(cTotal_ssp370[,,4,1200], zlim=c(0,120000),legend.lab='Ctot_370')
image.plot(cTotal_ssp585[,,4,1200], zlim=c(0,120000),legend.lab='Ctot_585')

par(mfrow=c(2,2))
plot(apply(cTotal_ssp126[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,70000))
lines(apply(cTotal_ssp126[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(cTotal_ssp126[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(cTotal_ssp245[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,70000))
lines(apply(cTotal_ssp245[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(cTotal_ssp245[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(cTotal_ssp370[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,70000))
lines(apply(cTotal_ssp370[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(cTotal_ssp370[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(cTotal_ssp585[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,70000))
lines(apply(cTotal_ssp585[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(cTotal_ssp585[,,7,],3,mean, na.rm=TRUE), col="red")
# increases and slows in the highest 7. lowest estimate seems lower than startiing time step. 

# Plotting for FLUXES
# More duificult than the stock data as constantly changes so if we plot maps its just the result for that timestep. 
# Avg the last 5 years of observed data (2014 -2019) vs the las 5 year of the the porjections

# Plotting gpp
# Slice for quantiles
gpp_ssp126_4 <- gpp_ssp126[,,4,]
gpp_ssp245_4 <- gpp_ssp245[,,4,]
gpp_ssp370_4 <- gpp_ssp370[,,4,]
gpp_ssp585_4 <- gpp_ssp585[,,4,]

# For 50% confidence intervals
gpp_ssp126_3 <- gpp_ssp126[,,3,]
gpp_ssp245_3 <- gpp_ssp245[,,3,]
gpp_ssp370_3 <- gpp_ssp370[,,3,]
gpp_ssp585_3 <- gpp_ssp585[,,3,]

gpp_ssp126_5 <- gpp_ssp126[,,5,]
gpp_ssp245_5 <- gpp_ssp245[,,5,]
gpp_ssp370_5 <- gpp_ssp370[,,5,]
gpp_ssp585_5 <- gpp_ssp585[,,5,]

# For 90% confidence intervals
gpp_ssp126_1 <- gpp_ssp126[,,1,]
gpp_ssp245_1 <- gpp_ssp245[,,1,]
gpp_ssp370_1 <- gpp_ssp370[,,1,]
gpp_ssp585_1 <- gpp_ssp585[,,1,]

gpp_ssp126_7 <- gpp_ssp126[,,7,]
gpp_ssp245_7 <- gpp_ssp245[,,7,]
gpp_ssp370_7 <- gpp_ssp370[,,7,]
gpp_ssp585_7 <- gpp_ssp585[,,7,]



# create new object which will contain the aggregated to mean annual rather than monthly
# gpp ssp126
gpp_ssp126_4_annual = array(NA, dim=c(dim(gpp_ssp126_4)[1],dim(gpp_ssp126_4)[2],dim(gpp_ssp126_4)[3]/steps_per_year))
for (i in seq(1,dim(gpp_ssp126_4)[1])) {
     for (j in seq(1,dim(gpp_ssp126_4)[2])) {
          gpp_ssp126_4_annual[i,j,] = rollapply(gpp_ssp126_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

gpp_ssp126_1_annual = array(NA, dim=c(dim(gpp_ssp126_1)[1],dim(gpp_ssp126_1)[2],dim(gpp_ssp126_1)[3]/steps_per_year))
for (i in seq(1,dim(gpp_ssp126_1)[1])) {
     for (j in seq(1,dim(gpp_ssp126_1)[2])) {
          gpp_ssp126_1_annual[i,j,] = rollapply(gpp_ssp126_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

gpp_ssp126_7_annual = array(NA, dim=c(dim(gpp_ssp126_7)[1],dim(gpp_ssp126_7)[2],dim(gpp_ssp126_7)[3]/steps_per_year))
for (i in seq(1,dim(gpp_ssp126_7)[1])) {
     for (j in seq(1,dim(gpp_ssp126_7)[2])) {
          gpp_ssp126_7_annual[i,j,] = rollapply(gpp_ssp126_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# gpp ssp245
gpp_ssp245_4_annual = array(NA, dim=c(dim(gpp_ssp245_4)[1],dim(gpp_ssp245_4)[2],dim(gpp_ssp245_4)[3]/steps_per_year))
for (i in seq(1,dim(gpp_ssp245_4)[1])) {
     for (j in seq(1,dim(gpp_ssp245_4)[2])) {
          gpp_ssp245_4_annual[i,j,] = rollapply(gpp_ssp245_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

gpp_ssp245_1_annual = array(NA, dim=c(dim(gpp_ssp245_1)[1],dim(gpp_ssp245_1)[2],dim(gpp_ssp245_1)[3]/steps_per_year))
for (i in seq(1,dim(gpp_ssp245_1)[1])) {
     for (j in seq(1,dim(gpp_ssp245_1)[2])) {
          gpp_ssp245_1_annual[i,j,] = rollapply(gpp_ssp245_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

gpp_ssp245_7_annual = array(NA, dim=c(dim(gpp_ssp245_7)[1],dim(gpp_ssp245_7)[2],dim(gpp_ssp245_7)[3]/steps_per_year))
for (i in seq(1,dim(gpp_ssp245_7)[1])) {
     for (j in seq(1,dim(gpp_ssp245_7)[2])) {
          gpp_ssp245_7_annual[i,j,] = rollapply(gpp_ssp245_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# Plot the annual averages trough time
par(mfrow=c(2,2)
plot(apply(gpp_ssp126_4_annual,3,mean, na.rm=TRUE), type="l", ylim=c(0,12), xlab="Year since 2000", ylab="Mean annual GPP values under the ssp126")
lines(apply(gpp_ssp126_1_annual,3,mean, na.rm=TRUE), col="blue")
lines(apply(gpp_ssp126_7_annual,3,mean, na.rm=TRUE), col="red")

plot(apply(gpp_ssp245_4_annual,3,mean, na.rm=TRUE), type="l", ylim=c(0,12), xlab="Year since 2000", ylab="Mean annual GPP values under the ssp245")
lines(apply(gpp_ssp245_1_annual,3,mean, na.rm=TRUE), col="blue")
lines(apply(gpp_ssp245_7_annual,3,mean, na.rm=TRUE), col="red")

# Plotting gpp through time with 50% confidecne intervals
par(mfrow=c(2,2))
plot(apply(gpp_ssp126[,,4,],3,mean, na.rm=TRUE), type="l", ylim = c(0,12))
lines(apply(gpp_ssp126[,,3,],3,mean, na.rm=TRUE), col="blue")
lines(apply(gpp_ssp126[,,5,],3,mean, na.rm=TRUE), col="red")

plot(apply(gpp_ssp245[,,4,],3,mean, na.rm=TRUE), type="l", ylim = c(0,12))
lines(apply(gpp_ssp245[,,3,],3,mean, na.rm=TRUE), col="blue")
lines(apply(gpp_ssp245[,,5,],3,mean, na.rm=TRUE), col="red")

plot(apply(gpp_ssp370[,,4,],3,mean, na.rm=TRUE), type="l", ylim = c(0,12))
lines(apply(gpp_ssp370[,,3,],3,mean, na.rm=TRUE), col="blue")
lines(apply(gpp_ssp370[,,5,],3,mean, na.rm=TRUE), col="red")

plot(apply(gpp_ssp585[,,4,],3,mean, na.rm=TRUE), type="l", ylim = c(0,12))
lines(apply(gpp_ssp585[,,3,],3,mean, na.rm=TRUE), col="blue")
lines(apply(gpp_ssp585[,,5,],3,mean, na.rm=TRUE), col="red")

gpp_ssp126_avgfirst <- apply(fire_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
gpp_ssp126_avglast <- apply(gpp_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
gpp_ssp245_avglast <- apply(gpp_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
gpp_ssp370_avglast <- apply(gpp_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
gpp_ssp585_avglast <- apply(gpp_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)

# Plotting fire flux avgs over last 5 years
fire_ssp126_avglast <- apply(fire_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
fire_ssp245_avglast <- apply(fire_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
fire_ssp370_avglast <- apply(fire_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
fire_ssp585_avglast <- apply(fire_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
par(mfrow= c(2,2))
image.plot(fire_ssp126_avglast, zlim =c(0,1.4))
image.plot(fire_ssp245_avglast, zlim =c(0,1.4))
image.plot(fire_ssp370_avglast, zlim =c(0,1.4))
image.plot(fire_ssp585_avglast, zlim =c(0,1.4))
# Fire does increase but seems to be concentrated in that one reagion. This is arround the serengeti/mara so fire is expected.

# Plotting fire flux avgs over first 5 years
fire_ssp126_avgfirst <- apply(fire_ssp126[,,4,1:60], c(1,2), mean, na.rm=TRUE)
image.plot(fire_ssp126_avgfirst, zlim =c(0,1.4))
# Much lower than the final 5 years
image.plot(diif

par(mfrow=c(2,2))
plot(apply(fire_ssp126[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,0.5))
lines(apply(fire_ssp126[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(fire_ssp126[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(fire_ssp245[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,0.5))
lines(apply(fire_ssp245[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(fire_ssp245[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(fire_ssp370[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,0.5))
lines(apply(fire_ssp370[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(fire_ssp370[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(fire_ssp585[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,0.5))
lines(apply(fire_ssp585[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(fire_ssp585[,,7,],3,mean, na.rm=TRUE), col="red")
# Fire regime increases in 245, is much lower in 370 and climbs much higher in 585.
# Not the case anymore, it increases with higher ssp 

# Plotting NEE
nee_ssp126_avglast <- apply(nee_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nee_ssp245_avglast <- apply(nee_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nee_ssp370_avglast <- apply(nee_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nee_ssp585_avglast <- apply(nee_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
# Plotting NEE over the first 5 years
nee_ssp126_avgfirst <- apply(nee_ssp126[,,4,1:60], c(1,2), mean, na.rm=TRUE)

par(mfrow=c(3,2))
image.plot(nee_ssp126_avgfirst, zlim =c(-4.5,1),legend.lab='NEE_2001')
image.plot(nee_ssp126_avglast, zlim =c(-4.5,1),legend.lab='NEE_126')
image.plot(nee_ssp245_avglast, zlim =c(-4.5,1),legend.lab='NEE_245')
image.plot(nee_ssp370_avglast, zlim =c(-4.5,1),legend.lab='NEE_370') 
image.plot(nee_ssp585_avglast, zlim =c(-4.5,1),legend.lab='NEE_585')
# General decrease in NEE up ssps



par(mfrow=c(2,2))
plot(apply(nee_ssp126[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(-4,4))
lines(apply(nee_ssp126[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(nee_ssp126[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(nee_ssp245[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(-4,4))
lines(apply(nee_ssp245[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(nee_ssp245[,,7,],3,mean, na.rm=TRUE), col="red")
# at high estimationg 7 there is limited nee, all slightly +ve

plot(apply(nee_ssp370[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(-4,4))
lines(apply(nee_ssp370[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(nee_ssp370[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(nee_ssp585[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(-4,4))
lines(apply(nee_ssp585[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(nee_ssp585[,,7,],3,mean, na.rm=TRUE), col="red")
# high estimates are much more in the positive

# Plotting Forest harvest
fLuc_ssp126_avglast <- apply(fLuc_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
fLuc_ssp245_avglast <- apply(fLuc_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
fLuc_ssp370_avglast <- apply(fLuc_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
fLuc_ssp585_avglast <- apply(fLuc_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
fLuc_ssp126_avgfirst <- apply(fLuc_ssp126[,,4,1:60], c(1,2), mean, na.rm=TRUE)
par(mfrow=c(3,2))
image.plot(fLuc_ssp126_avgfirst, zlim =c(0,1.4))
image.plot(fLuc_ssp126_avglast, zlim =c(0,1.4))
image.plot(fLuc_ssp245_avglast, zlim =c(0,1.4))
image.plot(fLuc_ssp370_avglast, zlim =c(0,1.4))
image.plot(fLuc_ssp585_avglast, zlim =c(0,1.4))
#by 2100 all have very very low forest harvest

par(mfrow=c(2,2))
plot(apply(fLuc_ssp126[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,0.6))
lines(apply(fLuc_ssp126[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(fLuc_ssp126[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(fLuc_ssp245[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,0.6))
lines(apply(fLuc_ssp245[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(fLuc_ssp245[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(fLuc_ssp370[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,0.6))
lines(apply(fLuc_ssp370[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(fLuc_ssp370[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(fLuc_ssp585[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,0.6))
lines(apply(fLuc_ssp585[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(fLuc_ssp585[,,7,],3,mean, na.rm=TRUE), col="red")

# Plotting Net Biome Productivity
par(mfrow=c(3,2))
image.plot(nbp_ssp126[,,4,1], zlim=c(0,12))
image.plot(nbp_ssp126[,,4,1200], zlim=c(0,12))
image.plot(nbp_ssp245[,,4,1200], zlim=c(0,12))
image.plot(nbp_ssp370[,,4,1200], zlim=c(0,12))
image.plot(nbp_ssp585[,,4,1200], zlim=c(0,12))

# Plotting the last 5 years average nbp
nbp_ssp126_avglast <- apply(nbp_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nbp_ssp245_avglast <- apply(nbp_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nbp_ssp370_avglast <- apply(nbp_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nbp_ssp585_avglast <- apply(nbp_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nbp_ssp126_avgfirst <- apply(nbp_ssp126[,,4,1:6], c(1,2), mean, na.rm=TRUE)
par(mfrow=c(3,2))
image.plot(nbp_ssp126_avgfirst, zlim=c(-1,5))
image.plot(nbp_ssp126_avglast, zlim=c(-1,5))
image.plot(nbp_ssp245_avglast, zlim=c(-1,5))
image.plot(nbp_ssp370_avglast, zlim=c(-1,5))
image.plot(nbp_ssp585_avglast, zlim=c(-1,5))
# General increase in NBP in different ssp scenarios

par(mfrow=c(2,2))
plot(apply(nbp_ssp126[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(-3,7))
lines(apply(nbp_ssp126[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(nbp_ssp126[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(nbp_ssp245[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(-3,7))
lines(apply(nbp_ssp245[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(nbp_ssp245[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(nbp_ssp370[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(-3,7))
lines(apply(nbp_ssp370[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(nbp_ssp370[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(nbp_ssp585[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(-3,7))
lines(apply(nbp_ssp585[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(nbp_ssp585[,,7,],3,mean, na.rm=TRUE), col="red")

# Plotting Internal:Ambiant CO2 ratio
par(mfrow=c(3,2))
image.plot(CiCa_ssp126[,,4,1])
image.plot(CiCa_ssp126[,,4,1200])
image.plot(CiCa_ssp245[,,4,1200])
image.plot(CiCa_ssp370[,,4,1200])
image.plot(CiCa_ssp585[,,4,1200])
#by 2100 ssp

par(mfrow=c(2,2))
plot(apply(CiCa_ssp126[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,1))
lines(apply(CiCa_ssp126[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(CiCa_ssp126[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(CiCa_ssp245[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,1))
lines(apply(CiCa_ssp245[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(CiCa_ssp245[,,7,],3,mean, na.rm=TRUE), col="red")
# ssp245 the ratio decreases but then begins to increase again. 

plot(apply(CiCa_ssp370[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,1))
lines(apply(CiCa_ssp370[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(CiCa_ssp370[,,7,],3,mean, na.rm=TRUE), col="red")

plot(apply(CiCa_ssp585[,,4,],3,mean, na.rm=TRUE), type="l", ylim=c(0,1))
lines(apply(CiCa_ssp585[,,1,],3,mean, na.rm=TRUE), col="blue")
lines(apply(CiCa_ssp585[,,7,],3,mean, na.rm=TRUE), col="red")
# ssp585 the ratio deacreases then stabalises 

# Plotting Drivers 
# Different as its only 3 dimensions doesnt have the 7 different scenarios

# Precipitation
par(mfrow=c(3,2))
image.plot(pr_ssp126[,,1], zlim= c(0,0.0002),legend.lab='Pr_2001')
image.plot(pr_ssp126[,,1200], zlim= c(0,0.0002),legend.lab='Pr_126')
image.plot(pr_ssp245[,,1200], zlim= c(0,0.0002),legend.lab='Pr_245')
image.plot(pr_ssp370[,,1200], zlim= c(0,0.0002),legend.lab='Pr_370')
image.plot(pr_ssp585[,,1200], zlim= c(0,0.0002),legend.lab='Pr_585')

par(mfrow=c(2,2))
plot(apply(pr_ssp126,3,mean, na.rm=TRUE), type="l", ylim=c(0,0.00018))
plot(apply(pr_ssp245,3,mean, na.rm=TRUE), type="l", ylim=c(0,0.00018))
plot(apply(pr_ssp370,3,mean, na.rm=TRUE), type="l", ylim=c(0,0.00018))
plot(apply(pr_ssp585,3,mean, na.rm=TRUE), type="l", ylim=c(0,0.00018))

# CO2 conc plots
par(mfrow=c(3,2))
plot(apply(co2_ssp126,3,mean, na.rm=TRUE), type="l", ylim=c(300,1200))
plot(apply(co2_ssp245,3,mean, na.rm=TRUE), type="l", ylim=c(300,1200))
plot(apply(co2_ssp370,3,mean, na.rm=TRUE), type="l", ylim=c(300,1200))
plot(apply(co2_ssp585,3,mean, na.rm=TRUE), type="l", ylim=c(300,1200))

# Observed LAI
par(mfrow=c(2,2))
plot(apply(laiobs_ssp126,3,mean, na.rm=TRUE), type="l")
image.plot(laiobs_ssp126[,,1])# LAI at very start
image.plot(laiobs_ssp126[,,228])# LAI at end of obs

# Making data frames and tables to read off final values.
# Get averages of first and last 5 years of gpp
gpp_ssp126_4_avgfirst <- apply(gpp_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
gpp_ssp126_4_avglast <- apply(gpp_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
gpp_ssp245_4_avglast <- apply(gpp_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
gpp_ssp370_4_avglast <- apply(gpp_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
gpp_ssp585_4_avglast <- apply(gpp_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
par(mfrow= c(3,2))
image.plot(gpp_ssp126_avgfirst, zlim =c(0,16))
image.plot(gpp_ssp126_avglast, zlim =c(0,16))
image.plot(gpp_ssp245_avglast, zlim =c(0,16))
image.plot(gpp_ssp370_avglast, zlim =c(0,16))
image.plot(gpp_ssp585_avglast, zlim =c(0,16))

gpp_ssp126_4 <-gpp_ssp126[,,4,]

# Calculating the area of a pixel 

###
## Function to determine the area of a pixel at any location in meter squared
###
resolution <- 0.25

## function to convert degrees to radians
deg2rad<-function(degree_in) {return(degree_in*(pi/180))}
rad2deg<-function(radian_in) {return(radian_in*(180/pi))}

calc_pixel_area<-function(lat,lon,resolution) {

    # resolution in degrees
    # lat (-90,90) degrees
    # lon (-180,180) degrees

    # mean earth radius (m)
    R = 6371e3

    pixel_area = R**2 * ( deg2rad(lon+resolution*0.5)-deg2rad(lon-resolution*0.5) ) * (sin(deg2rad(lat+resolution*0.5))-sin(deg2rad(lat-resolution*0.5)))

    # return to user in meters
    return(pixel_area)

}
# pixel area = 770886443 m.2

# Creating functions for conversion
MtC.yr <- function(x){ # Converts each pixel from gC.m-2.yr-1 to MtC.yr-1
	x*365.25*pixel_area/1000000000000 
	}

tC.ha.yr <- function(x){
	x*365.25*10000/1000000
	}

plot(apply(gpp_ssp126[,,4,168:288],3,mean, na.rm=TRUE), type="l")

gpp_ssp126_4_MtC <- apply(gpp_ssp126_4_avgfirst,c(1,2), MtC.yr)
sum(gpp_ssp126_4_MtC, na.rm=TRUE) # returns 209.0653 MtC.yr-1

# GPP tC.ha-1.yr-1
# With scaling
gpp_ssp126_4_tC <- apply(gpp_ssp126_4[,,168:288], c(1,2), scale) # scaling still not working 
gpp_ssp126_4_tC <- apply(gpp_ssp126_4_tC,c(2,3), mean, na.rm=TRUE)
gpp_ssp126_4_tC <- apply(gpp_ssp126_4_tC,c(1,2), tC.ha.yr)
mean(gpp_ssp126_4_tC, na.rm=TRUE) # returns -1569e-17

# Without scaling
gpp_ssp126_4_tC <- apply(gpp_ssp126_4_avgfirst,c(1,2), mean, na.rm=TRUE)
gpp_ssp126_4_tC <- apply(gpp_ssp126_4_tC,c(1,2), tC.ha.yr)
mean(gpp_ssp126_4_tC, na.rm=TRUE) # returns 8.144179




