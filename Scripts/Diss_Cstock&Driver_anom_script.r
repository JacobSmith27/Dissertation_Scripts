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

# DRIVERS
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

# FLux for Fluc
# Read in SSP126
cflux_ssp126 = nc_open("Kenya_CFLUX_MOHC_ssp126_climate_change_2001_2100.nc")

# Read in SSP245
cflux_ssp245 = nc_open("Kenya_CFLUX_MOHC_ssp245_climate_change_2001_2100.nc")

# Read in SSP370
cflux_ssp370 = nc_open("Kenya_CFLUX_MOHC_ssp370_climate_change_2001_2100.nc")

# Read in SSP585
cflux_ssp585 = nc_open("Kenya_CFLUX_MOHC_ssp585_climate_change_2001_2100.nc")

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

# Drivers data
# Driver data only has 3 dimensions not 4 like the rest as it doesnt have the 7 different estimates.
# Mean precipitation (kg.m-2.s-1)
pr_ssp126 = ncvar_get(Drivers_ssp126,"pr",count = c(19,27,1200)) 
pr_ssp245 = ncvar_get(Drivers_ssp245,"pr",count = c(19,27,1200))
pr_ssp370 = ncvar_get(Drivers_ssp370,"pr",count = c(19,27,1200))
pr_ssp585 = ncvar_get(Drivers_ssp585,"pr",count = c(19,27,1200))

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

# Forest loss fraction
flf_ssp126 = ncvar_get(Drivers_ssp126,"forest_loss_fraction",count = c(19,27,1200))
flf_ssp245 = ncvar_get(Drivers_ssp245,"forest_loss_fraction",count = c(19,27,1200))
flf_ssp370 = ncvar_get(Drivers_ssp370,"forest_loss_fraction",count = c(19,27,1200))
flf_ssp585 = ncvar_get(Drivers_ssp585,"forest_loss_fraction",count = c(19,27,1200))

# Forest loss fraction
bf_ssp126 = ncvar_get(Drivers_ssp126,"burnt_fraction",count = c(19,27,1200))
bf_ssp245 = ncvar_get(Drivers_ssp245,"burnt_fraction",count = c(19,27,1200))
bf_ssp370 = ncvar_get(Drivers_ssp370,"burnt_fraction",count = c(19,27,1200))
bf_ssp585 = ncvar_get(Drivers_ssp585,"burnt_fraction",count = c(19,27,1200))

# Reading in for fluc
# Forest harvet (gC.m-2.d-1)
fLuc_ssp126 = ncvar_get(cflux_ssp126,"fLuc_ensemble",count = c(19,27,7,1200))
fLuc_ssp245 = ncvar_get(cflux_ssp245,"fLuc_ensemble",count = c(19,27,7,1200))
fLuc_ssp370 = ncvar_get(cflux_ssp370,"fLuc_ensemble",count = c(19,27,7,1200))
fLuc_ssp585 = ncvar_get(cflux_ssp585,"fLuc_ensemble",count = c(19,27,7,1200))

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
nc_close(Drivers_ssp126);nc_close(Drivers_ssp245);
nc_close(Drivers_ssp370);nc_close(Drivers_ssp585);nc_close(Drivers_ssp434)
nc_close(cflux_ssp126) ; nc_close(cflux_ssp245) ; nc_close(cflux_ssp370) ; nc_close(cflux_ssp585) 


# Converting drivers into a dataframe

# Read in Kenyas shape file 
Kenya_shp <- shapefile("~/Kenya_shapefiles/County.shp")


# Creating unit convertion functions
tC.ha <- function(x){ # Converts each pixel from gC.m-2 to tC.ha-1
	x/100
	}

pr_conv <- function(x){ # Converts each pixel from kg.m-2.s-1 to mm.d-1
	x*86400
	}

tC.ha.yr <- function(x){ # Converts each pixel from gC.m-2.d-1 to tC.ha-1.yr-1
	x*365.25/1000
	}

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

# cVeg
bio_ssp126_anom <- bio_ssp126_4_1200 - bio_ssp126_4_1
bio_ssp245_anom <- bio_ssp245_4_1200 - bio_ssp126_4_1
bio_ssp370_anom <- bio_ssp370_4_1200 - bio_ssp126_4_1
bio_ssp585_anom <- bio_ssp585_4_1200 - bio_ssp126_4_1

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

# DOM
dom_ssp126_anom <- dom_ssp126_4_1200 - dom_ssp126_4_1
dom_ssp245_anom <- dom_ssp245_4_1200 - dom_ssp126_4_1
dom_ssp370_anom <- dom_ssp370_4_1200 - dom_ssp126_4_1
dom_ssp585_anom <- dom_ssp585_4_1200 - dom_ssp126_4_1

# pr
# Slicing through to just getthe required time references 
pr_ssp126_1 <- apply(pr_ssp126[3:19,5:27,168:228], c(1,2), mean, na.rm=TRUE)
pr_ssp126_1200 <- apply(pr_ssp126[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
pr_ssp245_1200 <- apply(pr_ssp245[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
pr_ssp370_1200 <- apply(pr_ssp370[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
pr_ssp585_1200 <- apply(pr_ssp585[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)

pr_ssp126_1 <- apply(pr_ssp126_1,c(1,2), pr_conv)
pr_ssp126_1200 <- apply(pr_ssp126_1200,c(1,2), pr_conv)
pr_ssp245_1200 <- apply(pr_ssp245_1200,c(1,2), pr_conv)
pr_ssp370_1200 <- apply(pr_ssp370_1200,c(1,2), pr_conv)
pr_ssp585_1200 <- apply(pr_ssp585_1200,c(1,2), pr_conv)

# pr
pr_ssp126_anom <- pr_ssp126_1200 - pr_ssp126_1
pr_ssp245_anom <- pr_ssp245_1200 - pr_ssp126_1
pr_ssp370_anom <- pr_ssp370_1200 - pr_ssp126_1
pr_ssp585_anom <- pr_ssp585_1200 - pr_ssp126_1

# vpd
# Slicing through to just get the required time references 
vpd_ssp126_1 <- apply(vpd_ssp126[3:19,5:27,168:228], c(1,2), mean, na.rm=TRUE)
vpd_ssp126_1200 <- apply(vpd_ssp126[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
vpd_ssp245_1200 <- apply(vpd_ssp245[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
vpd_ssp370_1200 <- apply(vpd_ssp370[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
vpd_ssp585_1200 <- apply(vpd_ssp585[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)

# vpd
vpd_ssp126_anom <- vpd_ssp126_1200 - vpd_ssp126_1
vpd_ssp245_anom <- vpd_ssp245_1200 - vpd_ssp126_1
vpd_ssp370_anom <- vpd_ssp370_1200 - vpd_ssp126_1
vpd_ssp585_anom <- vpd_ssp585_1200 - vpd_ssp126_1

# tas_min
# Slicing through to just get the required time references 
tas_min_ssp126_1 <- apply(tas_min_ssp126[3:19,5:27,168:228], c(1,2), mean, na.rm=TRUE)
tas_min_ssp126_1200 <- apply(tas_min_ssp126[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
tas_min_ssp245_1200 <- apply(tas_min_ssp245[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
tas_min_ssp370_1200 <- apply(tas_min_ssp370[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
tas_min_ssp585_1200 <- apply(tas_min_ssp585[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)

# tas_max
# Slicing through to just getthe required time references 
tas_max_ssp126_1 <- apply(tas_max_ssp126[3:19,5:27,168:228], c(1,2), mean, na.rm=TRUE)
tas_max_ssp126_1200 <- apply(tas_max_ssp126[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
tas_max_ssp245_1200 <- apply(tas_max_ssp245[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
tas_max_ssp370_1200 <- apply(tas_max_ssp370[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
tas_max_ssp585_1200 <- apply(tas_max_ssp585[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)

# tas_mean
# calculating mean temp
tas_mean_ssp126_1 <- ((tas_max_ssp126_1+tas_min_ssp126_1)/2)
tas_mean_ssp126_1200 <- ((tas_max_ssp126_1200+tas_min_ssp126_1200)/2)
tas_mean_ssp245_1200 <- ((tas_max_ssp245_1200+tas_min_ssp245_1200)/2)
tas_mean_ssp370_1200 <- ((tas_max_ssp370_1200+tas_min_ssp370_1200)/2)
tas_mean_ssp585_1200 <- ((tas_max_ssp585_1200+tas_min_ssp585_1200)/2)

# tas_mean_anom
tas_mean_ssp126_anom <- tas_mean_ssp126_1200 - tas_mean_ssp126_1
tas_mean_ssp245_anom <- tas_mean_ssp245_1200 - tas_mean_ssp126_1
tas_mean_ssp370_anom <- tas_mean_ssp370_1200 - tas_mean_ssp126_1
tas_mean_ssp585_anom <- tas_mean_ssp585_1200 - tas_mean_ssp126_1

# rsds
# Slicing through to just get the required time references 
rsds_ssp126_1 <- apply(rsds_ssp126[3:19,5:27,168:228], c(1,2), mean, na.rm=TRUE)
rsds_ssp126_1200 <- apply(rsds_ssp126[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
rsds_ssp245_1200 <- apply(rsds_ssp245[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
rsds_ssp370_1200 <- apply(rsds_ssp370[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
rsds_ssp585_1200 <- apply(rsds_ssp585[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)

# rsds
rsds_ssp126_anom <- rsds_ssp126_1200 - rsds_ssp126_1
rsds_ssp245_anom <- rsds_ssp245_1200 - rsds_ssp126_1
rsds_ssp370_anom <- rsds_ssp370_1200 - rsds_ssp126_1
rsds_ssp585_anom <- rsds_ssp585_1200 - rsds_ssp126_1

# flf
# Slicing through to just get the required time references 
flf_ssp126_1 <- apply(flf_ssp126[3:19,5:27,168:228], c(1,2), mean, na.rm=TRUE)
flf_ssp126_1200 <- apply(flf_ssp126[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
flf_ssp245_1200 <- apply(flf_ssp245[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
flf_ssp370_1200 <- apply(flf_ssp370[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
flf_ssp585_1200 <- apply(flf_ssp585[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)

# flf
flf_ssp126_anom <- flf_ssp126_1200 - flf_ssp126_1
flf_ssp245_anom <- flf_ssp245_1200 - flf_ssp126_1
flf_ssp370_anom <- flf_ssp370_1200 - flf_ssp126_1
flf_ssp585_anom <- flf_ssp585_1200 - flf_ssp126_1

# bf
# Slicing through to just get the required time references 
bf_ssp126_1 <- apply(bf_ssp126[3:19,5:27,168:228], c(1,2), mean, na.rm=TRUE)
bf_ssp126_1200 <- apply(bf_ssp126[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
bf_ssp245_1200 <- apply(bf_ssp245[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
bf_ssp370_1200 <- apply(bf_ssp370[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)
bf_ssp585_1200 <- apply(bf_ssp585[3:19,5:27,1140:1200], c(1,2), mean, na.rm=TRUE)

# bf
bf_ssp126_anom <- bf_ssp126_1200 - bf_ssp126_1
bf_ssp245_anom <- bf_ssp245_1200 - bf_ssp126_1
bf_ssp370_anom <- bf_ssp370_1200 - bf_ssp126_1
bf_ssp585_anom <- bf_ssp585_1200 - bf_ssp126_1

# fLuc
# Slicing through to just get the required time references 
fLuc_ssp126_1 <- apply(fLuc_ssp126[3:19,5:27,4,168:228], c(1,2), mean, na.rm=TRUE)
fLuc_ssp126_1200 <- apply(fLuc_ssp126[3:19,5:27,4,1080:1140], c(1,2), mean, na.rm=TRUE)
fLuc_ssp245_1200 <- apply(fLuc_ssp245[3:19,5:27,4,1080:1140], c(1,2), mean, na.rm=TRUE)
fLuc_ssp370_1200 <- apply(fLuc_ssp370[3:19,5:27,4,1080:1140], c(1,2), mean, na.rm=TRUE)
fLuc_ssp585_1200 <- apply(fLuc_ssp585[3:19,5:27,4,1080:1140], c(1,2), mean, na.rm=TRUE)

fLuc_ssp126_1 <- apply(fLuc_ssp126_1,c(1,2), tC.ha.yr)
fLuc_ssp126_1200 <- apply(fLuc_ssp126_1200,c(1,2), tC.ha.yr)
fLuc_ssp245_1200 <- apply(fLuc_ssp245_1200,c(1,2), tC.ha.yr)
fLuc_ssp370_1200 <- apply(fLuc_ssp370_1200,c(1,2), tC.ha.yr)
fLuc_ssp585_1200 <- apply(fLuc_ssp585_1200,c(1,2), tC.ha.yr)

# fLuc
fLuc_ssp126_anom <- fLuc_ssp126_1200 - fLuc_ssp126_1
fLuc_ssp245_anom <- fLuc_ssp245_1200 - fLuc_ssp126_1
fLuc_ssp370_anom <- fLuc_ssp370_1200 - fLuc_ssp126_1
fLuc_ssp585_anom <- fLuc_ssp585_1200 - fLuc_ssp126_1

# PLotting cVeg anomally
par(mfrow=c(2,2))
image.plot(lon, lat, bio_ssp126_anom,
main = "SSP126",
xlab = "Longitude",
ylab = "Latitude",cex.main=1.6,
legend.line = 4,
zlim=c(-35,90),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

image.plot(lon, lat, bio_ssp245_anom,
main = "SSP245",
xlab = "Longitude",
ylab = "Latitude",cex.main=1.6,
legend.line = 4,
zlim=c(-35,90),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

image.plot(lon, lat, bio_ssp370_anom,
main = "SSP370",
xlab = "Longitude",
ylab = "Latitude",cex.main=1.6,
legend.line = 4,
zlim=c(-35,90),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

image.plot(lon, lat, bio_ssp585_anom,
main = "SSP585",
xlab = "Longitude",
ylab = "Latitude",cex.main=1.6,
legend.line = 4,
zlim=c(-35,90),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# Plotting DOM anommaly
par(mfrow=c(2,2))
image.plot(lon, lat, dom_ssp126_anom,
main = "SSP126",
xlab = "Longitude",
ylab = "Latitude",cex.main=1.6,
legend.line = 4,
zlim=c(-300,400),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

image.plot(lon, lat, dom_ssp245_anom,
main = "SSP245",
xlab = "Longitude",
ylab = "Latitude",cex.main=1.6,
legend.line = 4,
zlim=c(-300,400),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

image.plot(lon, lat, dom_ssp370_anom,
main = "SSP370",
xlab = "Longitude",
ylab = "Latitude",cex.main=1.6,
legend.line = 4,
zlim=c(-300,400),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

image.plot(lon, lat, dom_ssp585_anom,
main = "SSP585",
xlab = "Longitude",
ylab = "Latitude",cex.main=1.6,
legend.line = 4,
zlim=c(-300,400),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# bio anom vs climate
par(mfrow=c(5,4),mar=c(3,5.5,2,3.2))
# bio 126
image.plot(lon, lat, bio_ssp126_anom,
main = "SSP126",
xaxt = "n", xlab = "",
ylab = "Latitude",cex.main=1.5,
legend.line = 4,
zlim=c(-35,90),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")
mtext("Bio (tC/ha)", side = 2, line=4)

# bio 245
image.plot(lon, lat, bio_ssp245_anom,
main = "SSP245",
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 4,
zlim=c(-35,90),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# bio370
image.plot(lon, lat, bio_ssp370_anom,
main = "SSP370",
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 4,
zlim=c(-35,90),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# bio 585
image.plot(lon, lat, bio_ssp585_anom,
main = "SSP585",
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 4,
zlim=c(-35,90),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# pr 126
image.plot(lon, lat, pr_ssp126_anom,
xaxt = "n", xlab = "",
ylab = "Latitude",cex.main=1.5,
legend.line = 1.6,
zlim=c(0,3.5),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")
mtext("Precipitation (mm/day)", side = 2, line=4)

# pr 245
image.plot(lon, lat, pr_ssp245_anom,
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(0,3.5),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# pr 370
image.plot(lon, lat, pr_ssp370_anom,
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(0,3.5),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# pr 585
image.plot(lon, lat, pr_ssp585_anom,
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(0,3.5),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

#tas_mean 126
image.plot(lon, lat, tas_mean_ssp126_anom,
xaxt = "n", xlab = "",
ylab = "Latitude",cex.main=1.5,
legend.line = 1.6,
zlim=c(0,7),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")
mtext("Mean air temperature (c)", side = 2, line=4)

# tas_mean 245
image.plot(lon, lat, tas_mean_ssp245_anom,
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(0,7),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# tas_mean 370
image.plot(lon, lat, tas_mean_ssp370_anom,
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(0,7),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# tas_mean 585
image.plot(lon, lat, tas_mean_ssp585_anom,
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(0,7),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

#vpd 126
image.plot(lon, lat, vpd_ssp126_anom,
xaxt = "n", xlab = "",
ylab = "Latitude",cex.main=1.5,
legend.line = 1.6,
zlim=c(-150,750),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")
mtext("Mean VPD (pa)", side = 2, line=4)

# vpd 245
image.plot(lon, lat, vpd_ssp245_anom,
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(-150,750),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# vpd 370
image.plot(lon, lat, vpd_ssp370_anom,
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(-150,750),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# vpd 585
image.plot(lon, lat, vpd_ssp585_anom,
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(-150,750),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

#rsds 126
image.plot(lon, lat, rsds_ssp126_anom,
xlab = "Longitude",
ylab = "Latitude",cex.main=1.5,
legend.line = 1.6,
zlim=c(-0.5,1.5),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")
mtext("Mean SWR (MJ/m2/day)", side = 2, line=4)

# rsds 245
image.plot(lon, lat, rsds_ssp245_anom,
xlab = "Longitude",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(-0.5,1.5),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# rsds 370
image.plot(lon, lat, rsds_ssp370_anom,
xlab = "Longitude",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(-0.5,1.5),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# rsds 585
image.plot(lon, lat, rsds_ssp585_anom,
xlab = "Longitude",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(-0.5,1.5),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# Other drivers map

# bio anom vs climate
par(mfrow=c(3,4),mar=c(4.1,5.2,2,4.5))
# bio 126
image.plot(lon, lat, bio_ssp126_anom,
main = "SSP126",
xaxt = "n", xlab = "",
ylab = "Latitude",cex.main=1.5,
legend.line = 4,
zlim=c(-35,90),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")
mtext("Bio (tC/ha)", side = 2, line=4)

# bio 245
image.plot(lon, lat, bio_ssp245_anom,
main = "SSP245",
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 4,
zlim=c(-35,90),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# bio370
image.plot(lon, lat, bio_ssp370_anom,
main = "SSP370",
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 4,
zlim=c(-35,90),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# bio 585
image.plot(lon, lat, bio_ssp585_anom,
main = "SSP585",
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 4,
zlim=c(-35,90),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# fLuc 126
image.plot(lon, lat, fLuc_ssp126_anom,
xaxt = "n", xlab = "",
ylab = "Latitude",cex.main=1.5,
legend.line = 1.6,
zlim=c(-0.141,0.6),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")
mtext("Forrest Harvest (tC/ha/yr)", side = 2, line=4)

# fLuc 245
image.plot(lon, lat, fLuc_ssp245_anom,
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(-0.141,0.6),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# fLuc 370
image.plot(lon, lat, fLuc_ssp370_anom,
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(-0.141,0.6),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# fLuc 585
image.plot(lon, lat, fLuc_ssp585_anom,
xaxt = "n", xlab = "",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(-0.141,0.6),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

#bf 126
image.plot(lon, lat, bf_ssp126_anom,
xlab = "Longitude",
ylab = "Latitude",cex.main=1.5,
legend.line = 1.6,
zlim=c(-0.0004,0.018),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")
mtext("Burned fraction", side = 2, line=4)

# bf 245
image.plot(lon, lat, bf_ssp245_anom,
xlab = "Longitude",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(-0.0003,0.018),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# bf 370
image.plot(lon, lat, bf_ssp370_anom,
xlab = "Longitude",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(-0.0003,0.018),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

# bf 585
image.plot(lon, lat, bf_ssp585_anom,
xlab = "Longitude",
yaxt = "n", ylab = "",cex.main=1.5,
legend.line = 1.6,
zlim=c(-0.0003,0.018),
col = plasma(200))
plot(Kenya_shp, add = T, border = "Black")

