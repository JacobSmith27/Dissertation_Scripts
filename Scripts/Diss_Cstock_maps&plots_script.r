#Load libraries
library(ncdf4)
library(fields)
library(dplyr)
library(ggplot2)
library(raster)
library(zoo)
library(rgdal)
library(viridis)
library(RColourBrewer)

# STOCKS
# Read in SSP126
cstock_ssp126 = nc_open("Kenya_CSTOCK_MOHC_ssp126_climate_change_2001_2100.nc")

# Read in SSP245
cstock_ssp245 = nc_open("Kenya_CSTOCK_MOHC_ssp245_climate_change_2001_2100.nc")

# Read in SSP370
cstock_ssp370 = nc_open("Kenya_CSTOCK_MOHC_ssp370_climate_change_2001_2100.nc")

# Read in SSP585
cstock_ssp585 = nc_open("Kenya_CSTOCK_MOHC_ssp585_climate_change_2001_2100.nc")


# Reading in from STOCKS
cur_ssp126 = ncvar_get(cstock_ssp126,"cVeg_ensemble")
cur_ssp245 = ncvar_get(cstock_ssp245,"cVeg_ensemble")
cur_ssp370 = ncvar_get(cstock_ssp370,"cVeg_ensemble")
cur_ssp585 = ncvar_get(cstock_ssp585,"cVeg_ensemble")

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

# Get lat and lon
lon <- ncvar_get(cstock_ssp126, "lon")
lat <- ncvar_get(cstock_ssp126, "lat")
lon1 <- lon
lat1 <- lat
lon <- lon[3:19]
lat <- lat[5:27]

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
nc_close(cstock_ssp126); nc_close(cstock_ssp245); nc_close(cstock_ssp370); nc_close(cstock_ssp585)

# Read in Kenyas shape file 
Kenya_shp <- shapefile("/home/s1635288/Desktop/Kenya_shapefiles/County.shp")

# Calculating the area of a pixel 
## Function to determine the area of a pixel at any location in meter squared
###
resolution <- 0.25
lat2 <- -5.1750000
lon2 <- 33.12500

## function to convert degrees to radians
deg2rad<-function(degree_in) {return(degree_in*(pi/180))}
rad2deg<-function(radian_in) {return(radian_in*(180/pi))}

calc_pixel_area<-function(lat,lon,resolution) {

    # resolution in degrees
    # lat (-90,90) degrees
    # lon (-180,180) degrees

    # mean earth radius (m)
    R = 6371e3

    pixel_area = R**2 * (deg2rad(lon+resolution*0.5)-deg2rad(lon-resolution*0.5) ) * (sin(deg2rad(lat+resolution*0.5))-sin(deg2rad(lat-resolution*0.5)))

    # return to user in meters
    return(pixel_area)

}

pixel_area <- calc_pixel_area(lat,lon,0.25)
# pixel area = 769618947 m.2

# Creating functions for conversion
MtC.ha <- function(x){ # Converts each pixel from gC.m-2 to MtC.ha-1
	x* 10000/1000000000000
}

# Unit conversion functions
tC.ha <- function(x){ # Converts each pixel from gC.m-2 to tC.ha-1
	x/100
	}

# Creating necessary maps
# Live biomass 
# Making better maps
# Slicing through to just get the medians (4) and the required time references 
bio_ssp126_4_1 <- apply(bio_ssp126[3:19,5:275:27,4,168:228], c(1,2), mean, na.rm=TRUE)
bio_ssp126_4_1200 <- apply(bio_ssp126[3:19,5:27,4,1140:1200], c(1,2), mean, na.rm=TRUE)
bio_ssp245_4_1200 <- apply(bio_ssp245[3:19,5:27,4,1140:1200], c(1,2), mean, na.rm=TRUE)
bio_ssp370_4_1200 <- apply(bio_ssp370[3:19,5:27,4,1140:1200], c(1,2), mean, na.rm=TRUE)
bio_ssp585_4_1200 <- apply(bio_ssp585[3:19,5:27,4,1140:1200], c(1,2), mean, na.rm=TRUE)

bio_ssp126_4_1 <- apply(bio_ssp126_4_1,c(1,2), tC.ha)
bio_ssp126_4_1200 <- apply(bio_ssp126_4_1200,c(1,2), tC.ha)
bio_ssp245_4_1200 <- apply(bio_ssp245_4_1200,c(1,2), tC.ha)
bio_ssp370_4_1200 <- apply(bio_ssp370_4_1200,c(1,2), tC.ha)
bio_ssp585_4_1200 <- apply(bio_ssp585_4_1200,c(1,2), tC.ha)

par(mfrow=c(2,3))
image.plot(lon, lat, bio_ssp126_4_1,
main = "Current(2014-2019) ",
xlab = "Longitude",
ylab = "Latitude",
cex.main=1.6,
legend.line = 4,
zlim=c(0,100),
col = viridis(200))
plot(Kenya_shp, add = T, border = "Black")

image.plot(lon, lat, bio_ssp126_4_1200,
main = "SSP126",
xlab = "Longitude",
ylab = "Latitude",
cex.main=1.6,
legend.line = 4,
zlim=c(0,100),
col = viridis(200))
plot(Kenya_shp, add = T, border = "Black")

image.plot(lon, lat, bio_ssp245_4_1200,
main = "SSP245",
xlab = "Longitude",
ylab = "Latitude",
cex.main=1.6,
legend.line = 4,
zlim=c(0,100),
col = viridis(200))
plot(Kenya_shp, add = T, border = "Black")

image.plot(lon, lat, bio_ssp370_4_1200,
main = "SSP370",
xlab = "Longitude",
ylab = "Latitude",
cex.main=1.6,
legend.line = 4,
zlim=c(0,100),
col = viridis(200))
plot(Kenya_shp, add = T, border = "Black")

image.plot(lon, lat, bio_ssp585_4_1200,
main = "SSP585",
xlab = "Longitude",
ylab = "Latitude",
cex.main=1.6,
legend.line = 4,
zlim=c(0,100),
col = viridis(200))
plot(Kenya_shp, add = T, border = "Black")

# Dead biomass (DOM)
# Slicing through to just get the medians (4) and the required time references 
dom_ssp126_4_1 <- apply(dom_ssp126[3:19,5:27,4,168:228], c(1,2), mean, na.rm=TRUE)
dom_ssp126_4_1200 <- apply(dom_ssp126[3:19,5:27,4,1140:1200], c(1,2), mean, na.rm=TRUE)
dom_ssp245_4_1200 <- apply(dom_ssp245[3:19,5:27,4,1140:1200], c(1,2), mean, na.rm=TRUE)
dom_ssp370_4_1200 <- apply(dom_ssp370[3:19,5:27,4,1140:1200], c(1,2), mean, na.rm=TRUE)
dom_ssp585_4_1200 <- apply(dom_ssp585[3:19,5:27,4,1140:1200], c(1,2), mean, na.rm=TRUE)

dom_ssp126_4_1 <- apply(dom_ssp126_4_1,c(1,2), tC.ha)
dom_ssp126_4_1200 <- apply(dom_ssp126_4_1200,c(1,2), tC.ha)
dom_ssp245_4_1200 <- apply(dom_ssp245_4_1200,c(1,2), tC.ha)
dom_ssp370_4_1200 <- apply(dom_ssp370_4_1200,c(1,2), tC.ha)
dom_ssp585_4_1200 <- apply(dom_ssp585_4_1200,c(1,2), tC.ha)

par(mfrow=c(2,3))
image.plot(lon, lat, dom_ssp126_4_1,
main = "Current (2014-2019)",
xlab = "Longitude",
ylab = "Latitude",cex.main=1.6,
legend.line = 4,
zlim=c(0,600),
col = viridis(200))
plot(Kenya_shp, add = T, border = "Black")

image.plot(lon, lat, dom_ssp126_4_1200,
main = "SSP126",
xlab = "Longitude",
ylab = "Latitude",cex.main=1.6,
legend.line = 4,
zlim=c(0,600),
col = viridis(200))
plot(Kenya_shp, add = T, border = "Black")

image.plot(lon, lat, dom_ssp245_4_1200,
main = "SSP245",
xlab = "Longitude",
ylab = "Latitude",cex.main=1.6,
legend.line = 4,
zlim=c(0,600),
col = viridis(200))
plot(Kenya_shp, add = T, border = "Black")

image.plot(lon, lat, dom_ssp370_4_1200,
main = "SSP370",
xlab = "Longitude",
ylab = "Latitude",cex.main=1.6,
legend.line = 4,
zlim=c(0,600),
col = viridis(200))
plot(Kenya_shp, add = T, border = "Black")

image.plot(lon, lat, dom_ssp585_4_1200,
main = "SSP585",
xlab = "Longitude",
ylab = "Latitude",cex.main=1.6,
legend.line = 4,
zlim=c(0,600),
col = viridis(200))
plot(Kenya_shp, add = T, border = "Black")

# Creating through time plots

# LAI
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

# Solution from LUKE using zoo package
# Assuming LAI is in an array with dimensions (long,lat,time)
# In this example LAI is stored in a variable called lai_ssp126_4
# How many time steps per year?
steps_per_year = 12

# create new object which will contain the aggregated to mean annual rather than monthly
# ssp126
lai_ssp126_4_annual = array(NA, dim=c(dim(lai_ssp126_4)[1],dim(lai_ssp126_4)[2],dim(lai_ssp126_4)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp126_4)[1])) {
     for (j in seq(1,dim(lai_ssp126_4)[2])) {
          lai_ssp126_4_annual[i,j,] = rollapply(lai_ssp126_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

lai_ssp126_1_annual = array(NA, dim=c(dim(lai_ssp126_1)[1],dim(lai_ssp126_4)[2],dim(lai_ssp126_1)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp126_1)[1])) {
     for (j in seq(1,dim(lai_ssp126_1)[2])) {
          lai_ssp126_1_annual[i,j,] = rollapply(lai_ssp126_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

lai_ssp126_7_annual = array(NA, dim=c(dim(lai_ssp126_7)[1],dim(lai_ssp126_7)[2],dim(lai_ssp126_7)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp126_7)[1])) {
     for (j in seq(1,dim(lai_ssp126_7)[2])) {
          lai_ssp126_7_annual[i,j,] = rollapply(lai_ssp126_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# ssp245
lai_ssp245_4_annual = array(NA, dim=c(dim(lai_ssp245_4)[1],dim(lai_ssp245_4)[2],dim(lai_ssp245_4)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp245_4)[1])) {
     for (j in seq(1,dim(lai_ssp245_4)[2])) {
          lai_ssp245_4_annual[i,j,] = rollapply(lai_ssp245_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

lai_ssp245_1_annual = array(NA, dim=c(dim(lai_ssp245_1)[1],dim(lai_ssp245_4)[2],dim(lai_ssp245_1)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp245_1)[1])) {
     for (j in seq(1,dim(lai_ssp245_1)[2])) {
          lai_ssp245_1_annual[i,j,] = rollapply(lai_ssp245_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

lai_ssp245_7_annual = array(NA, dim=c(dim(lai_ssp245_7)[1],dim(lai_ssp245_7)[2],dim(lai_ssp245_7)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp245_7)[1])) {
     for (j in seq(1,dim(lai_ssp245_7)[2])) {
          lai_ssp245_7_annual[i,j,] = rollapply(lai_ssp245_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# ssp370
lai_ssp370_4_annual = array(NA, dim=c(dim(lai_ssp370_4)[1],dim(lai_ssp370_4)[2],dim(lai_ssp370_4)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp370_4)[1])) {
     for (j in seq(1,dim(lai_ssp370_4)[2])) {
          lai_ssp370_4_annual[i,j,] = rollapply(lai_ssp370_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

lai_ssp370_1_annual = array(NA, dim=c(dim(lai_ssp370_1)[1],dim(lai_ssp370_4)[2],dim(lai_ssp370_1)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp370_1)[1])) {
     for (j in seq(1,dim(lai_ssp370_1)[2])) {
          lai_ssp370_1_annual[i,j,] = rollapply(lai_ssp370_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

lai_ssp370_7_annual = array(NA, dim=c(dim(lai_ssp370_7)[1],dim(lai_ssp370_7)[2],dim(lai_ssp370_7)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp370_7)[1])) {
     for (j in seq(1,dim(lai_ssp370_7)[2])) {
          lai_ssp370_7_annual[i,j,] = rollapply(lai_ssp370_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# ssp585
lai_ssp585_4_annual = array(NA, dim=c(dim(lai_ssp585_4)[1],dim(lai_ssp585_4)[2],dim(lai_ssp585_4)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp585_4)[1])) {
     for (j in seq(1,dim(lai_ssp585_4)[2])) {
          lai_ssp585_4_annual[i,j,] = rollapply(lai_ssp585_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

lai_ssp585_1_annual = array(NA, dim=c(dim(lai_ssp585_1)[1],dim(lai_ssp585_4)[2],dim(lai_ssp585_1)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp585_1)[1])) {
     for (j in seq(1,dim(lai_ssp585_1)[2])) {
          lai_ssp585_1_annual[i,j,] = rollapply(lai_ssp585_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

lai_ssp585_7_annual = array(NA, dim=c(dim(lai_ssp585_7)[1],dim(lai_ssp585_7)[2],dim(lai_ssp585_7)[3]/steps_per_year))
for (i in seq(1,dim(lai_ssp585_7)[1])) {
     for (j in seq(1,dim(lai_ssp585_7)[2])) {
          lai_ssp585_7_annual[i,j,] = rollapply(lai_ssp585_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# Plot the annual averages.
par(mfrow=c(2,2))
plot(apply(lai_ssp126_4_annual,1,mean, na.rm=TRUE), type="l", ylim=c(0,12), xlab="Year since 2000", ylab="Mean annual LAI values under the ssp126")
lines(x,apply(lai_ssp126_1_annual,1,mean, na.rm=TRUE), col="blue")
lines(x,apply(lai_ssp126_7_annual,1,mean, na.rm=TRUE), col="red")

plot(apply(lai_ssp245_4_annual,1,mean, na.rm=TRUE), type="l", ylim=c(0,12), xlab="Year since 2000", ylab="Mean annual LAI values under the ssp245")
lines(x,apply(lai_ssp245_1_annual,1,mean, na.rm=TRUE), col="blue")
lines(x,apply(lai_ssp245_7_annual,1,mean, na.rm=TRUE), col="red")

plot(apply(lai_ssp370_4_annual,1,mean, na.rm=TRUE), type="l", ylim=c(0,12), xlab="Year since 2000", ylab="Mean annual LAI values under the ssp370")
lines(x,apply(lai_ssp370_1_annual,1,mean, na.rm=TRUE), col="blue")
lines(x,apply(lai_ssp370_7_annual,1,mean, na.rm=TRUE), col="red")

plot(apply(lai_ssp585_4_annual,1,mean, na.rm=TRUE), type="l", ylim=c(0,12), xlab="Year since 2000", ylab="Mean annual LAI values under the ssp585")
lines(x,apply(lai_ssp585_1_annual,1,mean, na.rm=TRUE), col="blue")
lines(x,apply(lai_ssp585_7_annual,1,mean, na.rm=TRUE), col="red")

plot(apply(lai_ssp126_4_annual,1,mean, na.rm=TRUE), type="l", ylim=c(0, 9000), xlab="Years since 2000", ylab="Carbon in living laimass (gC.m-2) ", col="goldenrod3", lwd=3)
lines(x,apply(lai_ssp126_1_annual,1,mean, na.rm=TRUE), type="l", lty=2, col="goldenrod3")
lines(x,apply(lai_ssp126_7_annual,1,mean, na.rm=TRUE), type="l", lty=2, col="goldenrod3")
lines(x,apply(lai_ssp245_4_annual,1,mean, na.rm=TRUE), type="l", lwd=3, col="darkgreen")
lines(x,apply(lai_ssp245_1_annual,1,mean, na.rm=TRUE), type="l", lty=2, col="darkgreen")
lines(x,apply(lai_ssp245_7_annual,1,mean, na.rm=TRUE), type="l", lty=2, col="darkgreen")
lines(x,apply(lai_ssp370_4_annual,1,mean, na.rm=TRUE), type="l", lwd=3, col ="blue")
lines(x,apply(lai_ssp370_1_annual,1,mean, na.rm=TRUE), type="l", lty=2, col ="blue")
lines(x,apply(lai_ssp370_7_annual,1,mean, na.rm=TRUE), type="l", lty=2, col ="blue")
lines(x,apply(lai_ssp585_4_annual,1,mean, na.rm=TRUE), type="l", lwd=3, col ="firebrick")
lines(x,apply(lai_ssp585_1_annual,1,mean, na.rm=TRUE), type="l", lty=2, col ="firebrick")
lines(x,apply(lai_ssp585_7_annual,1,mean, na.rm=TRUE), type="l", lty=2, col ="firebrick")
legend(x = "top",          # Position
       legend = c("SSP126", "SSP245", "SSP370", "SSP585"), 
	bty="n",
       col = c("goldenrod3", "darkgreen", "blue", "firebrick"), 
       lwd = 2)   

# cVeg
# Slicing through to just get the medians (4)
bio_ssp126_4 <- bio_ssp126[,,4,]
bio_ssp245_4 <- bio_ssp245[,,4,]
bio_ssp370_4 <- bio_ssp370[,,4,]
bio_ssp585_4 <- bio_ssp585[,,4,]

# Slicing through to just get 0.250(4)
bio_ssp126_1 <- bio_ssp126[,,1,]
bio_ssp245_1 <- bio_ssp245[,,1,]
bio_ssp370_1 <- bio_ssp370[,,1,]
bio_ssp585_1 <- bio_ssp585[,,1,]

# Slicing through to just get 0.750 (4)
bio_ssp126_7 <- bio_ssp126[,,7,]
bio_ssp245_7 <- bio_ssp245[,,7,]
bio_ssp370_7 <- bio_ssp370[,,7,]
bio_ssp585_7 <- bio_ssp585[,,7,]

# How many time steps per year?
steps_per_year = 12

# create new object which will contain the aggregated to mean annual rather than monthly
# ssp126
bio_ssp126_4_annual = array(NA, dim=c(dim(bio_ssp126_4)[1],dim(bio_ssp126_4)[2],dim(bio_ssp126_4)[3]/steps_per_year))
for (i in seq(1,dim(bio_ssp126_4)[1])) {
     for (j in seq(1,dim(bio_ssp126_4)[2])) {
          bio_ssp126_4_annual[i,j,] = rollapply(bio_ssp126_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

bio_ssp126_1_annual = array(NA, dim=c(dim(bio_ssp126_1)[1],dim(bio_ssp126_4)[2],dim(bio_ssp126_1)[3]/steps_per_year))
for (i in seq(1,dim(bio_ssp126_1)[1])) {
     for (j in seq(1,dim(bio_ssp126_1)[2])) {
          bio_ssp126_1_annual[i,j,] = rollapply(bio_ssp126_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

bio_ssp126_7_annual = array(NA, dim=c(dim(bio_ssp126_7)[1],dim(bio_ssp126_7)[2],dim(bio_ssp126_7)[3]/steps_per_year))
for (i in seq(1,dim(bio_ssp126_7)[1])) {
     for (j in seq(1,dim(bio_ssp126_7)[2])) {
          bio_ssp126_7_annual[i,j,] = rollapply(bio_ssp126_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# ssp245
bio_ssp245_4_annual = array(NA, dim=c(dim(bio_ssp245_4)[1],dim(bio_ssp245_4)[2],dim(bio_ssp245_4)[3]/steps_per_year))
for (i in seq(1,dim(bio_ssp245_4)[1])) {
     for (j in seq(1,dim(bio_ssp245_4)[2])) {
          bio_ssp245_4_annual[i,j,] = rollapply(bio_ssp245_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

bio_ssp245_1_annual = array(NA, dim=c(dim(bio_ssp245_1)[1],dim(bio_ssp245_4)[2],dim(bio_ssp245_1)[3]/steps_per_year))
for (i in seq(1,dim(bio_ssp245_1)[1])) {
     for (j in seq(1,dim(bio_ssp245_1)[2])) {
          bio_ssp245_1_annual[i,j,] = rollapply(bio_ssp245_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

bio_ssp245_7_annual = array(NA, dim=c(dim(bio_ssp245_7)[1],dim(bio_ssp245_7)[2],dim(bio_ssp245_7)[3]/steps_per_year))
for (i in seq(1,dim(bio_ssp245_7)[1])) {
     for (j in seq(1,dim(bio_ssp245_7)[2])) {
          bio_ssp245_7_annual[i,j,] = rollapply(bio_ssp245_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# ssp370
bio_ssp370_4_annual = array(NA, dim=c(dim(bio_ssp370_4)[1],dim(bio_ssp370_4)[2],dim(bio_ssp370_4)[3]/steps_per_year))
for (i in seq(1,dim(bio_ssp370_4)[1])) {
     for (j in seq(1,dim(bio_ssp370_4)[2])) {
          bio_ssp370_4_annual[i,j,] = rollapply(bio_ssp370_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

bio_ssp370_1_annual = array(NA, dim=c(dim(bio_ssp370_1)[1],dim(bio_ssp370_4)[2],dim(bio_ssp370_1)[3]/steps_per_year))
for (i in seq(1,dim(bio_ssp370_1)[1])) {
     for (j in seq(1,dim(bio_ssp370_1)[2])) {
          bio_ssp370_1_annual[i,j,] = rollapply(bio_ssp370_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

bio_ssp370_7_annual = array(NA, dim=c(dim(bio_ssp370_7)[1],dim(bio_ssp370_7)[2],dim(bio_ssp370_7)[3]/steps_per_year))
for (i in seq(1,dim(bio_ssp370_7)[1])) {
     for (j in seq(1,dim(bio_ssp370_7)[2])) {
          bio_ssp370_7_annual[i,j,] = rollapply(bio_ssp370_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# ssp585
bio_ssp585_4_annual = array(NA, dim=c(dim(bio_ssp585_4)[1],dim(bio_ssp585_4)[2],dim(bio_ssp585_4)[3]/steps_per_year))
for (i in seq(1,dim(bio_ssp585_4)[1])) {
     for (j in seq(1,dim(bio_ssp585_4)[2])) {
          bio_ssp585_4_annual[i,j,] = rollapply(bio_ssp585_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

bio_ssp585_1_annual = array(NA, dim=c(dim(bio_ssp585_1)[1],dim(bio_ssp585_4)[2],dim(bio_ssp585_1)[3]/steps_per_year))
for (i in seq(1,dim(bio_ssp585_1)[1])) {
     for (j in seq(1,dim(bio_ssp585_1)[2])) {
          bio_ssp585_1_annual[i,j,] = rollapply(bio_ssp585_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

bio_ssp585_7_annual = array(NA, dim=c(dim(bio_ssp585_7)[1],dim(bio_ssp585_7)[2],dim(bio_ssp585_7)[3]/steps_per_year))
for (i in seq(1,dim(bio_ssp585_7)[1])) {
     for (j in seq(1,dim(bio_ssp585_7)[2])) {
          bio_ssp585_7_annual[i,j,] = rollapply(bio_ssp585_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

bio_ssp126_4_annual <- apply(bio_ssp126_4_annual,c(1,2), tC.ha)
bio_ssp126_1_annual <- apply(bio_ssp126_1_annual,c(1,2), tC.ha)
bio_ssp126_7_annual <- apply(bio_ssp126_7_annual,c(1,2), tC.ha)
bio_ssp245_4_annual <- apply(bio_ssp245_4_annual,c(1,2), tC.ha)
bio_ssp245_1_annual <- apply(bio_ssp245_1_annual,c(1,2), tC.ha)
bio_ssp245_7_annual <- apply(bio_ssp245_7_annual,c(1,2), tC.ha)
bio_ssp370_4_annual <- apply(bio_ssp370_4_annual,c(1,2), tC.ha)
bio_ssp370_1_annual <- apply(bio_ssp370_1_annual,c(1,2), tC.ha)
bio_ssp370_7_annual <- apply(bio_ssp370_7_annual,c(1,2), tC.ha)
bio_ssp585_4_annual <- apply(bio_ssp585_4_annual,c(1,2), tC.ha)
bio_ssp585_1_annual <- apply(bio_ssp585_1_annual,c(1,2), tC.ha)
bio_ssp585_7_annual <- apply(bio_ssp585_7_annual,c(1,2), tC.ha)

# Creat axis to plot against
x <- c(2001:2100)

# Plotting
# Plotting Living biomass
plot(x,apply(bio_ssp126_4_annual,1,mean, na.rm=TRUE), type="l", ylim=c(0, 90), xlab="Year", ylab="Carbon in living biomass (tC/ha) ", col="goldenrod3", lwd=3)
lines(x,apply(bio_ssp126_1_annual,1,mean, na.rm=TRUE), type="l", lty=3, col="goldenrod3")
lines(x,apply(bio_ssp126_7_annual,1,mean, na.rm=TRUE), type="l", lty=3, col="goldenrod3")
lines(x,apply(bio_ssp245_4_annual,1,mean, na.rm=TRUE), type="l", lwd=3, col="darkgreen")
lines(x,apply(bio_ssp245_1_annual,1,mean, na.rm=TRUE), type="l", lty=3, col="darkgreen")
lines(x,apply(bio_ssp245_7_annual,1,mean, na.rm=TRUE), type="l", lty=3, col="darkgreen")
lines(x,apply(bio_ssp370_4_annual,1,mean, na.rm=TRUE), type="l", lwd=3, col ="blue")
lines(x,apply(bio_ssp370_1_annual,1,mean, na.rm=TRUE), type="l", lty=3, col ="blue")
lines(x,apply(bio_ssp370_7_annual,1,mean, na.rm=TRUE), type="l", lty=3, col ="blue")
lines(x,apply(bio_ssp585_4_annual,1,mean, na.rm=TRUE), type="l", lwd=3, col ="firebrick")
lines(x,apply(bio_ssp585_1_annual,1,mean, na.rm=TRUE), type="l", lty=3, col ="firebrick")
lines(x,apply(bio_ssp585_7_annual,1,mean, na.rm=TRUE), type="l", lty=3, col ="firebrick")
legend(x = "top",          # Position
       legend = c("SSP126", "SSP245", "SSP370", "SSP585"), 
	bty="n",
       col = c("goldenrod3", "darkgreen", "blue", "firebrick"), 
       lwd = 2)   

# DOM
# Slicing through to just get the medians (4)
dom_ssp126_4 <- dom_ssp126[,,4,]
dom_ssp245_4 <- dom_ssp245[,,4,]
dom_ssp370_4 <- dom_ssp370[,,4,]
dom_ssp585_4 <- dom_ssp585[,,4,]

# Slicing through to just get 0.250(4)
dom_ssp126_1 <- dom_ssp126[,,1,]
dom_ssp245_1 <- dom_ssp245[,,1,]
dom_ssp370_1 <- dom_ssp370[,,1,]
dom_ssp585_1 <- dom_ssp585[,,1,]

# Slicing through to just get 0.750 (4)
dom_ssp126_7 <- dom_ssp126[,,7,]
dom_ssp245_7 <- dom_ssp245[,,7,]
dom_ssp370_7 <- dom_ssp370[,,7,]
dom_ssp585_7 <- dom_ssp585[,,7,]

# How many time steps per year?
steps_per_year = 12

# create new object which will contain the aggregated to mean annual rather than monthly
# ssp126
dom_ssp126_4_annual = array(NA, dim=c(dim(dom_ssp126_4)[1],dim(dom_ssp126_4)[2],dim(dom_ssp126_4)[3]/steps_per_year))
for (i in seq(1,dim(dom_ssp126_4)[1])) {
     for (j in seq(1,dim(dom_ssp126_4)[2])) {
          dom_ssp126_4_annual[i,j,] = rollapply(dom_ssp126_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

dom_ssp126_1_annual = array(NA, dim=c(dim(dom_ssp126_1)[1],dim(dom_ssp126_4)[2],dim(dom_ssp126_1)[3]/steps_per_year))
for (i in seq(1,dim(dom_ssp126_1)[1])) {
     for (j in seq(1,dim(dom_ssp126_1)[2])) {
          dom_ssp126_1_annual[i,j,] = rollapply(dom_ssp126_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

dom_ssp126_7_annual = array(NA, dim=c(dim(dom_ssp126_7)[1],dim(dom_ssp126_7)[2],dim(dom_ssp126_7)[3]/steps_per_year))
for (i in seq(1,dim(dom_ssp126_7)[1])) {
     for (j in seq(1,dim(dom_ssp126_7)[2])) {
          dom_ssp126_7_annual[i,j,] = rollapply(dom_ssp126_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# ssp245
dom_ssp245_4_annual = array(NA, dim=c(dim(dom_ssp245_4)[1],dim(dom_ssp245_4)[2],dim(dom_ssp245_4)[3]/steps_per_year))
for (i in seq(1,dim(dom_ssp245_4)[1])) {
     for (j in seq(1,dim(dom_ssp245_4)[2])) {
          dom_ssp245_4_annual[i,j,] = rollapply(dom_ssp245_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

dom_ssp245_1_annual = array(NA, dim=c(dim(dom_ssp245_1)[1],dim(dom_ssp245_4)[2],dim(dom_ssp245_1)[3]/steps_per_year))
for (i in seq(1,dim(dom_ssp245_1)[1])) {
     for (j in seq(1,dim(dom_ssp245_1)[2])) {
          dom_ssp245_1_annual[i,j,] = rollapply(dom_ssp245_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

dom_ssp245_7_annual = array(NA, dim=c(dim(dom_ssp245_7)[1],dim(dom_ssp245_7)[2],dim(dom_ssp245_7)[3]/steps_per_year))
for (i in seq(1,dim(dom_ssp245_7)[1])) {
     for (j in seq(1,dim(dom_ssp245_7)[2])) {
          dom_ssp245_7_annual[i,j,] = rollapply(dom_ssp245_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# ssp370
dom_ssp370_4_annual = array(NA, dim=c(dim(dom_ssp370_4)[1],dim(dom_ssp370_4)[2],dim(dom_ssp370_4)[3]/steps_per_year))
for (i in seq(1,dim(dom_ssp370_4)[1])) {
     for (j in seq(1,dim(dom_ssp370_4)[2])) {
          dom_ssp370_4_annual[i,j,] = rollapply(dom_ssp370_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

dom_ssp370_1_annual = array(NA, dim=c(dim(dom_ssp370_1)[1],dim(dom_ssp370_4)[2],dim(dom_ssp370_1)[3]/steps_per_year))
for (i in seq(1,dim(dom_ssp370_1)[1])) {
     for (j in seq(1,dim(dom_ssp370_1)[2])) {
          dom_ssp370_1_annual[i,j,] = rollapply(dom_ssp370_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

dom_ssp370_7_annual = array(NA, dim=c(dim(dom_ssp370_7)[1],dim(dom_ssp370_7)[2],dim(dom_ssp370_7)[3]/steps_per_year))
for (i in seq(1,dim(dom_ssp370_7)[1])) {
     for (j in seq(1,dim(dom_ssp370_7)[2])) {
          dom_ssp370_7_annual[i,j,] = rollapply(dom_ssp370_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

# ssp585
dom_ssp585_4_annual = array(NA, dim=c(dim(dom_ssp585_4)[1],dim(dom_ssp585_4)[2],dim(dom_ssp585_4)[3]/steps_per_year))
for (i in seq(1,dim(dom_ssp585_4)[1])) {
     for (j in seq(1,dim(dom_ssp585_4)[2])) {
          dom_ssp585_4_annual[i,j,] = rollapply(dom_ssp585_4[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

dom_ssp585_1_annual = array(NA, dim=c(dim(dom_ssp585_1)[1],dim(dom_ssp585_4)[2],dim(dom_ssp585_1)[3]/steps_per_year))
for (i in seq(1,dim(dom_ssp585_1)[1])) {
     for (j in seq(1,dim(dom_ssp585_1)[2])) {
          dom_ssp585_1_annual[i,j,] = rollapply(dom_ssp585_1[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

dom_ssp585_7_annual = array(NA, dim=c(dim(dom_ssp585_7)[1],dim(dom_ssp585_7)[2],dim(dom_ssp585_7)[3]/steps_per_year))
for (i in seq(1,dim(dom_ssp585_7)[1])) {
     for (j in seq(1,dim(dom_ssp585_7)[2])) {
          dom_ssp585_7_annual[i,j,] = rollapply(dom_ssp585_7[i,j,], by = steps_per_year, width = steps_per_year, FUN=mean, na.rm=TRUE) } } 

dom_ssp126_4_annual <- apply(dom_ssp126_4_annual,c(1,2), tC.ha)
dom_ssp126_1_annual <- apply(dom_ssp126_1_annual,c(1,2), tC.ha)
dom_ssp126_7_annual <- apply(dom_ssp126_7_annual,c(1,2), tC.ha)
dom_ssp245_4_annual <- apply(dom_ssp245_4_annual,c(1,2), tC.ha)
dom_ssp245_1_annual <- apply(dom_ssp245_1_annual,c(1,2), tC.ha)
dom_ssp245_7_annual <- apply(dom_ssp245_7_annual,c(1,2), tC.ha)
dom_ssp370_4_annual <- apply(dom_ssp370_4_annual,c(1,2), tC.ha)
dom_ssp370_1_annual <- apply(dom_ssp370_1_annual,c(1,2), tC.ha)
dom_ssp370_7_annual <- apply(dom_ssp370_7_annual,c(1,2), tC.ha)
dom_ssp585_4_annual <- apply(dom_ssp585_4_annual,c(1,2), tC.ha)
dom_ssp585_1_annual <- apply(dom_ssp585_1_annual,c(1,2), tC.ha)
dom_ssp585_7_annual <- apply(dom_ssp585_7_annual,c(1,2), tC.ha)

# Plotting
# Plotting Living dommass
plot(x,apply(dom_ssp126_4_annual,1,mean, na.rm=TRUE), type="l", ylim=c(0,650), xlab="Year", ylab="Carbon in dead organic matter (tC/ha) ", col="goldenrod3", lwd=3)
lines(x,apply(dom_ssp126_1_annual,1,mean, na.rm=TRUE), type="l", lty=3, col="goldenrod3")
lines(x,apply(dom_ssp126_7_annual,1,mean, na.rm=TRUE), type="l", lty=3, col="goldenrod3")
lines(x,apply(dom_ssp245_4_annual,1,mean, na.rm=TRUE), type="l", lwd=3, col="darkgreen")
lines(x,apply(dom_ssp245_1_annual,1,mean, na.rm=TRUE), type="l", lty=3, col="darkgreen")
lines(x,apply(dom_ssp245_7_annual,1,mean, na.rm=TRUE), type="l", lty=3, col="darkgreen")
lines(x,apply(dom_ssp370_4_annual,1,mean, na.rm=TRUE), type="l", lwd=3, col ="blue")
lines(x,apply(dom_ssp370_1_annual,1,mean, na.rm=TRUE), type="l", lty=3, col ="blue")
lines(x,apply(dom_ssp370_7_annual,1,mean, na.rm=TRUE), type="l", lty=3, col ="blue")
lines(x,apply(dom_ssp585_4_annual,1,mean, na.rm=TRUE), type="l", lwd=3, col ="firebrick")
lines(x,apply(dom_ssp585_1_annual,1,mean, na.rm=TRUE), type="l", lty=3, col ="firebrick")
lines(x,apply(dom_ssp585_7_annual,1,mean, na.rm=TRUE), type="l", lty=3, col ="firebrick")
legend(x = "top",          # Position
       legend = c("SSP126", "SSP245", "SSP370", "SSP585"), 
	bty="n",
       col = c("goldenrod3", "darkgreen", "blue", "firebrick"), 
       lwd = 2)   
# Calculating total AGBC

cur_ssp126 <- apply(dom_ssp585_7_annual,c(1,2), tC.ha)

# Current stocks
cur_ssp126_4 <- apply(cur_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
cur_ssp126_4 <- apply(cur_ssp126_4,c(1,2), tC.ha)
image.plot(lon1, lat1, cur_ssp126_4,
xlab = "Longitude",
ylab = "Latitude",cex.main=1.6,
legend.line = 4,
col = viridis(200))
plot(Kenya_shp, add = T, border = "Black")

