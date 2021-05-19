#Load libraries
library(ncdf4)
library(fields)
library(dplyr)
library(lubridate)
library(raster)
library(zoo)
library(kableExtra)

# Read in SSP126
cflux_ssp126 = nc_open("Kenya_CFLUX_MOHC_ssp126_climate_change_2001_2100.nc")

# Read in SSP245
cflux_ssp245 = nc_open("Kenya_CFLUX_MOHC_ssp245_climate_change_2001_2100.nc")

# Read in SSP370
cflux_ssp370 = nc_open("Kenya_CFLUX_MOHC_ssp370_climate_change_2001_2100.nc")

# Read in SSP585
cflux_ssp585 = nc_open("Kenya_CFLUX_MOHC_ssp585_climate_change_2001_2100.nc")

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

# Net biome productivity (-NEE-fire-fluv)(gC.m-2.d-1)
nbp_ssp126 = ncvar_get(cflux_ssp126,"nbp_ensemble",count = c(19,27,7,1200))
nbp_ssp245 = ncvar_get(cflux_ssp245,"nbp_ensemble",count = c(19,27,7,1200))
nbp_ssp370 = ncvar_get(cflux_ssp370,"nbp_ensemble",count = c(19,27,7,1200))
nbp_ssp585 = ncvar_get(cflux_ssp585,"nbp_ensemble",count = c(19,27,7,1200))

# Get lat and lon
lon <- ncvar_get(cflux_ssp126, "lon")
lat <- ncvar_get(cflux_ssp126, "lat")
lon <- lon[1:19]
lat <- lat[1:27]

# Get quantile variable to identify 1-7 
quantiles = ncvar_get(cflux_ssp126,"quantile")
print(quantiles)
# 1 = 0.025
# 2 = 0.050
# 3 = 0.250
# 4 = 0.500
# 5 = 0.750
# 6 = 0.950
# 7 = 0.975

# When done with reading in variables remember to close the file
nc_close(cflux_ssp126) ; nc_close(cflux_ssp245) ; nc_close(cflux_ssp370) ; nc_close(cflux_ssp585) 

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

    pixel_area = R**2 * (deg2rad(lon+resolution*0.5)-deg2rad(lon-resolution*0.5) ) * (sin(deg2rad(lat+resolution*0.5))-sin(deg2rad(lat-resolution*0.5)))

    # return to user in meters
    return(pixel_area)

}

pixel_area <- calc_pixel_area(lat,lon,0.25)
# pixel area = 769618947 m.2

# Creating functions for conversion
MtC.yr <- function(x){ # Converts each pixel from gC.m-2.d-1 to MtC.yr-1
	x*365.25*pixel_area/1000000000000 
	}

tC.ha.yr <- function(x){ # Converts each pixel from cG.m-1.d-1 to tC.ha.yr
	x*365.25/100
	}

# Calculating FLUX values
# Avg the last 5 years of observed data (2014 -2019) vs the las 5 year of the the porjections
# gpp median
gpp_ssp126_4_avgfirst <- apply(gpp_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
gpp_ssp126_4_avglast <- apply(gpp_ssp126[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
gpp_ssp245_4_avglast <- apply(gpp_ssp245[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
gpp_ssp370_4_avglast <- apply(gpp_ssp370[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
gpp_ssp585_4_avglast <- apply(gpp_ssp585[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)

# gpp 25%
gpp_ssp126_3_avgfirst <- apply(gpp_ssp126[,,3,168:228], c(1,2), mean, na.rm=TRUE)
gpp_ssp126_3_avglast <- apply(gpp_ssp126[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
gpp_ssp245_3_avglast <- apply(gpp_ssp245[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
gpp_ssp370_3_avglast <- apply(gpp_ssp370[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
gpp_ssp585_3_avglast <- apply(gpp_ssp585[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)

# gpp 75%
gpp_ssp126_5_avgfirst <- apply(gpp_ssp126[,,5,168:228], c(1,2), mean, na.rm=TRUE)
gpp_ssp126_5_avglast <- apply(gpp_ssp126[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
gpp_ssp245_5_avglast <- apply(gpp_ssp245[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
gpp_ssp370_5_avglast <- apply(gpp_ssp370[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
gpp_ssp585_5_avglast <- apply(gpp_ssp585[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)

# ra median
ra_ssp126_4_avgfirst <- apply(ra_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
ra_ssp126_4_avglast <- apply(ra_ssp126[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
ra_ssp245_4_avglast <- apply(ra_ssp245[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
ra_ssp370_4_avglast <- apply(ra_ssp370[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
ra_ssp585_4_avglast <- apply(ra_ssp585[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)

# ra 25%
ra_ssp126_3_avgfirst <- apply(ra_ssp126[,,3,168:228], c(1,2), mean, na.rm=TRUE)
ra_ssp126_3_avglast <- apply(ra_ssp126[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
ra_ssp245_3_avglast <- apply(ra_ssp245[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
ra_ssp370_3_avglast <- apply(ra_ssp370[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
ra_ssp585_3_avglast <- apply(ra_ssp585[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)

# ra 75%
ra_ssp126_5_avgfirst <- apply(ra_ssp126[,,5,168:228], c(1,2), mean, na.rm=TRUE)
ra_ssp126_5_avglast <- apply(ra_ssp126[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
ra_ssp245_5_avglast <- apply(ra_ssp245[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
ra_ssp370_5_avglast <- apply(ra_ssp370[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
ra_ssp585_5_avglast <- apply(ra_ssp585[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)

# npp median

npp_ssp126_4_avgfirst <- apply(npp_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
npp_ssp126_4_avglast <- apply(npp_ssp126[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
npp_ssp245_4_avglast <- apply(npp_ssp245[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
npp_ssp370_4_avglast <- apply(npp_ssp370[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
npp_ssp585_4_avglast <- apply(npp_ssp585[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)

# npp 25%
npp_ssp126_3_avgfirst <- apply(npp_ssp126[,,3,168:228], c(1,2), mean, na.rm=TRUE)
npp_ssp126_3_avglast <- apply(npp_ssp126[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
npp_ssp245_3_avglast <- apply(npp_ssp245[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
npp_ssp370_3_avglast <- apply(npp_ssp370[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
npp_ssp585_3_avglast <- apply(npp_ssp585[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)

# npp 75%
npp_ssp126_5_avgfirst <- apply(npp_ssp126[,,5,168:228], c(1,2), mean, na.rm=TRUE)
npp_ssp126_5_avglast <- apply(npp_ssp126[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
npp_ssp245_5_avglast <- apply(npp_ssp245[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
npp_ssp370_5_avglast <- apply(npp_ssp370[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
npp_ssp585_5_avglast <- apply(npp_ssp585[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)

# rh median
rh_ssp126_4_avgfirst <- apply(rh_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
rh_ssp126_4_avglast <- apply(rh_ssp126[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
rh_ssp245_4_avglast <- apply(rh_ssp245[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
rh_ssp370_4_avglast <- apply(rh_ssp370[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
rh_ssp585_4_avglast <- apply(rh_ssp585[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)

# rh 25%
rh_ssp126_3_avgfirst <- apply(rh_ssp126[,,3,168:228], c(1,2), mean, na.rm=TRUE)
rh_ssp126_3_avglast <- apply(rh_ssp126[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
rh_ssp245_3_avglast <- apply(rh_ssp245[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
rh_ssp370_3_avglast <- apply(rh_ssp370[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
rh_ssp585_3_avglast <- apply(rh_ssp585[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)

# rh 75%
rh_ssp126_5_avgfirst <- apply(rh_ssp126[,,5,168:228], c(1,2), mean, na.rm=TRUE)
rh_ssp126_5_avglast <- apply(rh_ssp126[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
rh_ssp245_5_avglast <- apply(rh_ssp245[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
rh_ssp370_5_avglast <- apply(rh_ssp370[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
rh_ssp585_5_avglast <- apply(rh_ssp585[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)

# nee median
nee_ssp126_4_avgfirst <- apply(nee_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
nee_ssp126_4_avglast <- apply(nee_ssp126[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
nee_ssp245_4_avglast <- apply(nee_ssp245[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
nee_ssp370_4_avglast <- apply(nee_ssp370[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
nee_ssp585_4_avglast <- apply(nee_ssp585[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)

# nee 25%
nee_ssp126_3_avgfirst <- apply(nee_ssp126[,,3,168:228], c(1,2), mean, na.rm=TRUE)
nee_ssp126_3_avglast <- apply(nee_ssp126[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
nee_ssp245_3_avglast <- apply(nee_ssp245[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
nee_ssp370_3_avglast <- apply(nee_ssp370[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
nee_ssp585_3_avglast <- apply(nee_ssp585[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)

# nee 75%
nee_ssp126_5_avgfirst <- apply(nee_ssp126[,,5,168:228], c(1,2), mean, na.rm=TRUE)
nee_ssp126_5_avglast <- apply(nee_ssp126[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
nee_ssp245_5_avglast <- apply(nee_ssp245[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
nee_ssp370_5_avglast <- apply(nee_ssp370[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
nee_ssp585_5_avglast <- apply(nee_ssp585[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)

# fire median
fire_ssp126_4_avgfirst <- apply(fire_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
fire_ssp126_4_avglast <- apply(fire_ssp126[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
fire_ssp245_4_avglast <- apply(fire_ssp245[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
fire_ssp370_4_avglast <- apply(fire_ssp370[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
fire_ssp585_4_avglast <- apply(fire_ssp585[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)

# fire 25%
fire_ssp126_3_avgfirst <- apply(fire_ssp126[,,3,168:228], c(1,2), mean, na.rm=TRUE)
fire_ssp126_3_avglast <- apply(fire_ssp126[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
fire_ssp245_3_avglast <- apply(fire_ssp245[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
fire_ssp370_3_avglast <- apply(fire_ssp370[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
fire_ssp585_3_avglast <- apply(fire_ssp585[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)

# fire 75%
fire_ssp126_5_avgfirst <- apply(fire_ssp126[,,5,168:228], c(1,2), mean, na.rm=TRUE)
fire_ssp126_5_avglast <- apply(fire_ssp126[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
fire_ssp245_5_avglast <- apply(fire_ssp245[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
fire_ssp370_5_avglast <- apply(fire_ssp370[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
fire_ssp585_5_avglast <- apply(fire_ssp585[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)

# nbe median
nbe_ssp126_4_avgfirst <- apply(nbe_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
nbe_ssp126_4_avglast <- apply(nbe_ssp126[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
nbe_ssp245_4_avglast <- apply(nbe_ssp245[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
nbe_ssp370_4_avglast <- apply(nbe_ssp370[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
nbe_ssp585_4_avglast <- apply(nbe_ssp585[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)

# nbe 25%
nbe_ssp126_3_avgfirst <- apply(nbe_ssp126[,,3,168:228], c(1,2), mean, na.rm=TRUE)
nbe_ssp126_3_avglast <- apply(nbe_ssp126[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
nbe_ssp245_3_avglast <- apply(nbe_ssp245[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
nbe_ssp370_3_avglast <- apply(nbe_ssp370[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
nbe_ssp585_3_avglast <- apply(nbe_ssp585[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)

# nbe 75%
nbe_ssp126_5_avgfirst <- apply(nbe_ssp126[,,5,168:228], c(1,2), mean, na.rm=TRUE)
nbe_ssp126_5_avglast <- apply(nbe_ssp126[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
nbe_ssp245_5_avglast <- apply(nbe_ssp245[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
nbe_ssp370_5_avglast <- apply(nbe_ssp370[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
nbe_ssp585_5_avglast <- apply(nbe_ssp585[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)

# fLuc median
fLuc_ssp126_4_avgfirst <- apply(fLuc_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
fLuc_ssp126_4_avglast <- apply(fLuc_ssp126[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
fLuc_ssp245_4_avglast <- apply(fLuc_ssp245[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
fLuc_ssp370_4_avglast <- apply(fLuc_ssp370[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
fLuc_ssp585_4_avglast <- apply(fLuc_ssp585[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)

# fLuc 25%
fLuc_ssp126_3_avgfirst <- apply(fLuc_ssp126[,,3,168:228], c(1,2), mean, na.rm=TRUE)
fLuc_ssp126_3_avglast <- apply(fLuc_ssp126[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
fLuc_ssp245_3_avglast <- apply(fLuc_ssp245[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
fLuc_ssp370_3_avglast <- apply(fLuc_ssp370[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
fLuc_ssp585_3_avglast <- apply(fLuc_ssp585[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)

# fLuc 75%
fLuc_ssp126_5_avgfirst <- apply(fLuc_ssp126[,,5,168:228], c(1,2), mean, na.rm=TRUE)
fLuc_ssp126_5_avglast <- apply(fLuc_ssp126[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
fLuc_ssp245_5_avglast <- apply(fLuc_ssp245[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
fLuc_ssp370_5_avglast <- apply(fLuc_ssp370[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
fLuc_ssp585_5_avglast <- apply(fLuc_ssp585[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)

# nbp median
nbp_ssp126_4_avgfirst <- apply(nbp_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
nbp_ssp126_4_avglast <- apply(nbp_ssp126[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
nbp_ssp245_4_avglast <- apply(nbp_ssp245[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
nbp_ssp370_4_avglast <- apply(nbp_ssp370[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)
nbp_ssp585_4_avglast <- apply(nbp_ssp585[,,4,1128:1188], c(1,2), mean, na.rm=TRUE)

# nbp 25%
nbp_ssp126_3_avgfirst <- apply(nbp_ssp126[,,3,168:228], c(1,2), mean, na.rm=TRUE)
nbp_ssp126_3_avglast <- apply(nbp_ssp126[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
nbp_ssp245_3_avglast <- apply(nbp_ssp245[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
nbp_ssp370_3_avglast <- apply(nbp_ssp370[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)
nbp_ssp585_3_avglast <- apply(nbp_ssp585[,,3,1128:1188], c(1,2), mean, na.rm=TRUE)

# nbp 75%
nbp_ssp126_5_avgfirst <- apply(nbp_ssp126[,,5,168:228], c(1,2), mean, na.rm=TRUE)
nbp_ssp126_5_avglast <- apply(nbp_ssp126[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
nbp_ssp245_5_avglast <- apply(nbp_ssp245[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
nbp_ssp370_5_avglast <- apply(nbp_ssp370[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)
nbp_ssp585_5_avglast <- apply(nbp_ssp585[,,5,1128:1188], c(1,2), mean, na.rm=TRUE)

# GPP
plot(apply(gpp_ssp126[,,4,168:228],3,mean, na.rm=TRUE), type="l")

# MtC.yr-1
# 2014-2019
gpp_ssp126_4_avgfirst_MtC <- apply(gpp_ssp126_4_avgfirst,c(1,2), MtC.yr)
sum(gpp_ssp126_4_avgfirst_MtC, na.rm=TRUE) # returns 197.3466 MtC.yr-1

gpp_ssp126_3_avgfirst_MtC <- apply(gpp_ssp126_3_avgfirst,c(1,2), MtC.yr)
sum(gpp_ssp126_3_avgfirst_MtC, na.rm=TRUE) # returns 153.6578 MtC.yr-1

gpp_ssp126_5_avgfirst_MtC <- apply(gpp_ssp126_5_avgfirst,c(1,2), MtC.yr)
sum(gpp_ssp126_5_avgfirst_MtC, na.rm=TRUE) # returns 243.8919 MtC.yr-1

# SSP126
gpp_ssp126_4_MtC <- apply(gpp_ssp126_4_avglast,c(1,2), MtC.yr)
sum(gpp_ssp126_4_MtC, na.rm=TRUE) # returns 317.6131 MtC.yr-1

gpp_ssp126_3_MtC <- apply(gpp_ssp126_3_avglast,c(1,2), MtC.yr)
sum(gpp_ssp126_3_MtC, na.rm=TRUE) # returns 256.5847 MtC.yr-1

gpp_ssp126_5_MtC <- apply(gpp_ssp126_5_avglast,c(1,2), MtC.yr)
sum(gpp_ssp126_5_MtC, na.rm=TRUE) # returns 383.8129 MtC.yr-1

# SSP245
gpp_ssp245_4_MtC <- apply(gpp_ssp245_4_avglast,c(1,2), MtC.yr)
sum(gpp_ssp245_4_MtC, na.rm=TRUE) # returns 492.7542 MtC.yr-1

gpp_ssp245_3_MtC <- apply(gpp_ssp245_3_avglast,c(1,2), MtC.yr)
sum(gpp_ssp245_3_MtC, na.rm=TRUE) # returns 410.6502 MtC.yr-1

gpp_ssp245_5_MtC <- apply(gpp_ssp245_5_avglast,c(1,2), MtC.yr)
sum(gpp_ssp245_5_MtC, na.rm=TRUE) # returns 575.9434 MtC.yr-1

# SSP370
gpp_ssp370_4_MtC <- apply(gpp_ssp370_4_avglast,c(1,2), MtC.yr)
sum(gpp_ssp370_4_MtC, na.rm=TRUE) # returns 425.5089 MtC.yr-1

gpp_ssp370_3_MtC <- apply(gpp_ssp370_3_avglast,c(1,2), MtC.yr)
sum(gpp_ssp370_3_MtC, na.rm=TRUE) # returns 358.9473 MtC.yr-1

gpp_ssp370_5_MtC <- apply(gpp_ssp370_5_avglast,c(1,2), MtC.yr)
sum(gpp_ssp370_5_MtC, na.rm=TRUE) # returns 490.4101 MtC.yr-1

# SSP585
gpp_ssp585_4_MtC <- apply(gpp_ssp585_4_avglast,c(1,2), MtC.yr)
sum(gpp_ssp585_4_MtC, na.rm=TRUE) # returns 796.6199 MtC.yr-1

gpp_ssp585_3_MtC <- apply(gpp_ssp585_3_avglast,c(1,2), MtC.yr)
sum(gpp_ssp585_3_MtC, na.rm=TRUE) # returns 710.004 MtC.yr-1

gpp_ssp585_5_MtC <- apply(gpp_ssp585_5_avglast,c(1,2), MtC.yr)
sum(gpp_ssp585_5_MtC, na.rm=TRUE) # returns 895.7877 MtC.yr-1

# Converting the units of all matricies from gC.m-2.d-1
# GPP 
# 2014-2019
gpp_ssp126_4_avgfirst_tC <- apply(gpp_ssp126_4_avgfirst,c(1,2), tC.ha.yr)
gpp_ssp126_3_avgfirst_tC <- apply(gpp_ssp126_3_avgfirst,c(1,2), tC.ha.yr)
gpp_ssp126_5_avgfirst_tC <- apply(gpp_ssp126_5_avgfirst,c(1,2), tC.ha.yr)

# SSP126
gpp_ssp126_4_tC <- apply(gpp_ssp126_4_avglast,c(1,2), tC.ha.yr)
gpp_ssp126_3_tC <- apply(gpp_ssp126_3_avglast,c(1,2), tC.ha.yr)
gpp_ssp126_5_tC <- apply(gpp_ssp126_5_avglast,c(1,2), tC.ha.yr)

# SSP245
gpp_ssp245_4_tC <- apply(gpp_ssp245_4_avglast,c(1,2), tC.ha.yr)
median(gpp_ssp245_4_tC, na.rm=TRUE) # 24.467

gpp_ssp245_3_tC <- apply(gpp_ssp245_3_avglast,c(1,2), tC.ha.yr)
median(gpp_ssp245_3_tC, na.rm=TRUE) # 21.0587
gpp_ssp245_5_tC <- apply(gpp_ssp245_5_avglast,c(1,2), tC.ha.yr)
median(gpp_ssp245_5_tC, na.rm=TRUE) # 27.845

# SSP370
gpp_ssp370_4_tC <- apply(gpp_ssp370_4_avglast,c(1,2), tC.ha.yr)
median(gpp_ssp370_4_tC, na.rm=TRUE) # 20.273

gpp_ssp370_3_tC <- apply(gpp_ssp370_3_avglast,c(1,2), tC.ha.yr)
median(gpp_ssp370_3_tC, na.rm=TRUE) # 12.8413

gpp_ssp370_5_tC <- apply(gpp_ssp370_5_avglast,c(1,2), tC.ha.yr)
median(gpp_ssp370_5_tC, na.rm=TRUE) # 24.262

# SSP585
gpp_ssp585_4_tC <- apply(gpp_ssp585_4_avglast,c(1,2), tC.ha.yr)
median(gpp_ssp585_4_tC, na.rm=TRUE) # 39.442

gpp_ssp585_3_tC <- apply(gpp_ssp585_3_avglast,c(1,2), tC.ha.yr)
median(gpp_ssp585_3_tC, na.rm=TRUE) # 35.640

gpp_ssp585_5_tC <- apply(gpp_ssp585_5_avglast,c(1,2), tC.ha.yr)
median(gpp_ssp585_5_tC, na.rm=TRUE) # 43.279

# Ra
# 2014-2019
ra_ssp126_4_avgfirst_tC <- apply(ra_ssp126_4_avgfirst,c(1,2), tC.ha.yr)
median(ra_ssp126_4_avgfirst_tC, na.rm=TRUE) # returns 2.750 tC.ha-1.yr-1

ra_ssp126_3_avgfirst_tC <- apply(ra_ssp126_3_avgfirst,c(1,2), tC.ha.yr)
median(ra_ssp126_3_avgfirst_tC, na.rm=TRUE) # 1.941 tC.ha-1.yr-1

ra_ssp126_5_avgfirst_tC <- apply(ra_ssp126_5_avgfirst,c(1,2), tC.ha.yr)
median(ra_ssp126_5_avgfirst_tC, na.rm=TRUE) # 3.862 tC.ha-1.yr-1

# SSP126
ra_ssp126_4_tC <- apply(ra_ssp126_4_avglast,c(1,2), tC.ha.yr)
median(ra_ssp126_4_tC, na.rm=TRUE) # 6.537

ra_ssp126_3_tC <- apply(ra_ssp126_3_avglast,c(1,2), tC.ha.yr)
median(ra_ssp126_3_tC, na.rm=TRUE) # 4.863

ra_ssp126_5_tC <- apply(ra_ssp126_5_avglast,c(1,2), tC.ha.yr)
median(ra_ssp126_5_tC, na.rm=TRUE) # 8.536

# SSP245
ra_ssp245_4_tC <- apply(ra_ssp245_4_avglast,c(1,2), tC.ha.yr)
median(ra_ssp245_4_tC, na.rm=TRUE) # 10.240

ra_ssp245_3_tC <- apply(ra_ssp245_3_avglast,c(1,2), tC.ha.yr)
median(ra_ssp245_3_tC, na.rm=TRUE) # 7.916

ra_ssp245_5_tC <- apply(ra_ssp245_5_avglast,c(1,2), tC.ha.yr)
median(ra_ssp245_5_tC, na.rm=TRUE) # 13.102

# SSP370
ra_ssp370_4_tC <- apply(ra_ssp370_4_avglast,c(1,2), tC.ha.yr)
median(ra_ssp370_4_tC, na.rm=TRUE) # 8.070

ra_ssp370_3_tC <- apply(ra_ssp370_3_avglast,c(1,2), tC.ha.yr)
median(ra_ssp370_3_tC, na.rm=TRUE) # 4.691

ra_ssp370_5_tC <- apply(ra_ssp370_5_avglast,c(1,2), tC.ha.yr)
median(ra_ssp370_5_tC, na.rm=TRUE) # 9.923

# SSP585
ra_ssp585_4_tC <- apply(ra_ssp585_4_avglast,c(1,2), tC.ha.yr)
median(ra_ssp585_4_tC, na.rm=TRUE) # 16.783

ra_ssp585_3_tC <- apply(ra_ssp585_3_avglast,c(1,2), tC.ha.yr)
median(ra_ssp585_3_tC, na.rm=TRUE) # 13.155

ra_ssp585_5_tC <- apply(ra_ssp585_5_avglast,c(1,2), tC.ha.yr)
median(ra_ssp585_5_tC, na.rm=TRUE) # 20.624

# NPP
# 2014-2019
npp_ssp126_4_avgfirst_tC <- apply(npp_ssp126_4_avgfirst,c(1,2), tC.ha.yr)
median(npp_ssp126_4_avgfirst_tC, na.rm=TRUE) # returns 9.566 tC.ha-1.yr-1

npp_ssp126_3_avgfirst_tC <- apply(npp_ssp126_3_avgfirst,c(1,2), tC.ha.yr)
median(npp_ssp126_3_avgfirst_tC, na.rm=TRUE) # 6.997 tC.ha-1.yr-1

npp_ssp126_5_avgfirst_tC <- apply(npp_ssp126_5_avgfirst,c(1,2), tC.ha.yr)
median(npp_ssp126_5_avgfirst_tC, na.rm=TRUE) # 12.496 tC.ha-1.yr-1

# SSP126
npp_ssp126_4_tC <- apply(npp_ssp126_4_avglast,c(1,2), tC.ha.yr)
median(npp_ssp126_4_tC, na.rm=TRUE) # 21.815

npp_ssp126_3_tC <- apply(npp_ssp126_3_avglast,c(1,2), tC.ha.yr)
median(npp_ssp126_3_tC, na.rm=TRUE) # 17.632

npp_ssp126_5_tC <- apply(npp_ssp126_5_avglast,c(1,2), tC.ha.yr)
median(npp_ssp126_5_tC, na.rm=TRUE) # 26.503

# SSP245
npp_ssp245_4_tC <- apply(npp_ssp245_4_avglast,c(1,2), tC.ha.yr)
median(npp_ssp245_4_tC, na.rm=TRUE) # 34.911

npp_ssp245_3_tC <- apply(npp_ssp245_3_avglast,c(1,2), tC.ha.yr)
median(npp_ssp245_3_tC, na.rm=TRUE) # 29.304

npp_ssp245_5_tC <- apply(npp_ssp245_5_avglast,c(1,2), tC.ha.yr)
median(npp_ssp245_5_tC, na.rm=TRUE) # 40.4762

# SSP370
npp_ssp370_4_tC <- apply(npp_ssp370_4_avglast,c(1,2), tC.ha.yr)
median(npp_ssp370_4_tC, na.rm=TRUE) # 28.487

npp_ssp370_3_tC <- apply(npp_ssp370_3_avglast,c(1,2), tC.ha.yr)
median(npp_ssp370_3_tC, na.rm=TRUE) # 18.214

npp_ssp370_5_tC <- apply(npp_ssp370_5_avglast,c(1,2), tC.ha.yr)
median(npp_ssp370_5_tC, na.rm=TRUE) # 34.332

# SSP585
npp_ssp585_4_tC <- apply(npp_ssp585_4_avglast,c(1,2), tC.ha.yr)
median(npp_ssp585_4_tC, na.rm=TRUE) # 56.151

npp_ssp585_3_tC <- apply(npp_ssp585_3_avglast,c(1,2), tC.ha.yr)
median(npp_ssp585_3_tC, na.rm=TRUE) # 49.797

npp_ssp585_5_tC <- apply(npp_ssp585_5_avglast,c(1,2), tC.ha.yr)
median(npp_ssp585_5_tC, na.rm=TRUE) # 62.962

# Ra
# 2014-2019
rh_ssp126_4_avgfirst_tC <- apply(rh_ssp126_4_avgfirst,c(1,2), tC.ha.yr)
median(rh_ssp126_4_avgfirst_tC, na.rm=TRUE) # returns 3.366 tC.ha-1.yr-1

rh_ssp126_3_avgfirst_tC <- apply(rh_ssp126_3_avgfirst,c(1,2), tC.ha.yr)
median(rh_ssp126_3_avgfirst_tC, na.rm=TRUE) # 2.414 tC.ha-1.yr-1

rh_ssp126_5_avgfirst_tC <- apply(rh_ssp126_5_avgfirst,c(1,2), tC.ha.yr)
median(rh_ssp126_5_avgfirst_tC, na.rm=TRUE) # 4.545 tC.ha-1.yr-1

# SSP126
rh_ssp126_4_tC <- apply(rh_ssp126_4_avglast,c(1,2), tC.ha.yr)
median(rh_ssp126_4_tC, na.rm=TRUE) # 6.777

rh_ssp126_3_tC <- apply(rh_ssp126_3_avglast,c(1,2), tC.ha.yr)
median(rh_ssp126_3_tC, na.rm=TRUE) # 5.035

rh_ssp126_5_tC <- apply(rh_ssp126_5_avglast,c(1,2), tC.ha.yr)
median(rh_ssp126_5_tC, na.rm=TRUE) # 8.764

# SSP245
rh_ssp245_4_tC <- apply(rh_ssp245_4_avglast,c(1,2), tC.ha.yr)
median(rh_ssp245_4_tC, na.rm=TRUE) # 9.748

rh_ssp245_3_tC <- apply(rh_ssp245_3_avglast,c(1,2), tC.ha.yr)
median(rh_ssp245_3_tC, na.rm=TRUE) # 7.335

rh_ssp245_5_tC <- apply(rh_ssp245_5_avglast,c(1,2), tC.ha.yr)
median(rh_ssp245_5_tC, na.rm=TRUE) # 12.569

# SSP370
rh_ssp370_4_tC <- apply(rh_ssp370_4_avglast,c(1,2), tC.ha.yr)
median(rh_ssp370_4_tC, na.rm=TRUE) # 7.356

rh_ssp370_3_tC <- apply(rh_ssp370_3_avglast,c(1,2), tC.ha.yr)
median(rh_ssp370_3_tC, na.rm=TRUE) # 4.753

rh_ssp370_5_tC <- apply(rh_ssp370_5_avglast,c(1,2), tC.ha.yr)
median(rh_ssp370_5_tC, na.rm=TRUE) # 9.658

# SSP585
rh_ssp585_4_tC <- apply(rh_ssp585_4_avglast,c(1,2), tC.ha.yr)
median(rh_ssp585_4_tC, na.rm=TRUE) # 15.194

rh_ssp585_3_tC <- apply(rh_ssp585_3_avglast,c(1,2), tC.ha.yr)
median(rh_ssp585_3_tC, na.rm=TRUE) # 11.536

rh_ssp585_5_tC <- apply(rh_ssp585_5_avglast,c(1,2), tC.ha.yr)
median(rh_ssp585_5_tC, na.rm=TRUE) # 19.275

# nee
# 2014-2019
nee_ssp126_4_avgfirst_tC <- apply(nee_ssp126_4_avgfirst,c(1,2), tC.ha.yr)
median(nee_ssp126_4_avgfirst_tC, na.rm=TRUE) # returns -0.124 tC.ha-1.yr-1

nee_ssp126_3_avgfirst_tC <- apply(nee_ssp126_3_avgfirst,c(1,2), tC.ha.yr)
median(nee_ssp126_3_avgfirst_tC, na.rm=TRUE) # -0.842 tC.ha-1.yr-1

nee_ssp126_5_avgfirst_tC <- apply(nee_ssp126_5_avgfirst,c(1,2), tC.ha.yr)
median(nee_ssp126_5_avgfirst_tC, na.rm=TRUE) # 0.617 tC.ha-1.yr-1

# SSP126
nee_ssp126_4_tC <- apply(nee_ssp126_4_avglast,c(1,2), tC.ha.yr)
median(nee_ssp126_4_tC, na.rm=TRUE) # -0.950

nee_ssp126_3_tC <- apply(nee_ssp126_3_avglast,c(1,2), tC.ha.yr)
median(nee_ssp126_3_tC, na.rm=TRUE) # -1.871

nee_ssp126_5_tC <- apply(nee_ssp126_5_avglast,c(1,2), tC.ha.yr)
median(nee_ssp126_5_tC, na.rm=TRUE) # -0.0582

# SSP245
nee_ssp245_4_tC <- apply(nee_ssp245_4_avglast,c(1,2), tC.ha.yr)
median(nee_ssp245_4_tC, na.rm=TRUE) # -2.758

nee_ssp245_3_tC <- apply(nee_ssp245_3_avglast,c(1,2), tC.ha.yr)
median(nee_ssp245_3_tC, na.rm=TRUE) # -4.128

nee_ssp245_5_tC <- apply(nee_ssp245_5_avglast,c(1,2), tC.ha.yr)
median(nee_ssp245_5_tC, na.rm=TRUE) # -1.347

# SSP370
nee_ssp370_4_tC <- apply(nee_ssp370_4_avglast,c(1,2), tC.ha.yr)
median(nee_ssp370_4_tC, na.rm=TRUE) # -2.316

nee_ssp370_3_tC <- apply(nee_ssp370_3_avglast,c(1,2), tC.ha.yr)
median(nee_ssp370_3_tC, na.rm=TRUE) # -3.794

nee_ssp370_5_tC <- apply(nee_ssp370_5_avglast,c(1,2), tC.ha.yr)
median(nee_ssp370_5_tC, na.rm=TRUE) # -0.520

# SSP585
nee_ssp585_4_tC <- apply(nee_ssp585_4_avglast,c(1,2), tC.ha.yr)
median(nee_ssp585_4_tC, na.rm=TRUE) # -5.343

nee_ssp585_3_tC <- apply(nee_ssp585_3_avglast,c(1,2), tC.ha.yr)
median(nee_ssp585_3_tC, na.rm=TRUE) # -7.521

nee_ssp585_5_tC <- apply(nee_ssp585_5_avglast,c(1,2), tC.ha.yr)
median(nee_ssp585_5_tC, na.rm=TRUE) # -3.013

# Fire
# 2014-2019
fire_ssp126_4_avgfirst_tC <- apply(fire_ssp126_4_avgfirst,c(1,2), tC.ha.yr)
median(fire_ssp126_4_avgfirst_tC, na.rm=TRUE) # returns 0.00864 tC.ha-1.yr-1

fire_ssp126_3_avgfirst_tC <- apply(fire_ssp126_3_avgfirst,c(1,2), tC.ha.yr)
median(fire_ssp126_3_avgfirst_tC, na.rm=TRUE) # 0.00748 tC.ha-1.yr-1

fire_ssp126_5_avgfirst_tC <- apply(fire_ssp126_5_avgfirst,c(1,2), tC.ha.yr)
median(fire_ssp126_5_avgfirst_tC, na.rm=TRUE) # 0.0104 tC.ha-1.yr-1

# SSP126
fire_ssp126_4_tC <- apply(fire_ssp126_4_avglast,c(1,2), tC.ha.yr)
median(fire_ssp126_4_tC, na.rm=TRUE) # 0.0307

fire_ssp126_3_tC <- apply(fire_ssp126_3_avglast,c(1,2), tC.ha.yr)
median(fire_ssp126_3_tC, na.rm=TRUE) # 0.0230

fire_ssp126_5_tC <- apply(fire_ssp126_5_avglast,c(1,2), tC.ha.yr)
median(fire_ssp126_5_tC, na.rm=TRUE) # 0.0402

# SSP245
fire_ssp245_4_tC <- apply(fire_ssp245_4_avglast,c(1,2), tC.ha.yr)
median(fire_ssp245_4_tC, na.rm=TRUE) # 0.0437

fire_ssp245_3_tC <- apply(fire_ssp245_3_avglast,c(1,2), tC.ha.yr)
median(fire_ssp245_3_tC, na.rm=TRUE) # 0.0329

fire_ssp245_5_tC <- apply(fire_ssp245_5_avglast,c(1,2), tC.ha.yr)
median(fire_ssp245_5_tC, na.rm=TRUE) # 0.0563

# SSP370
fire_ssp370_4_tC <- apply(fire_ssp370_4_avglast,c(1,2), tC.ha.yr)
median(fire_ssp370_4_tC, na.rm=TRUE) # 0.0218

fire_ssp370_3_tC <- apply(fire_ssp370_3_avglast,c(1,2), tC.ha.yr)
median(fire_ssp370_3_tC, na.rm=TRUE) # 0.0145

fire_ssp370_5_tC <- apply(fire_ssp370_5_avglast,c(1,2), tC.ha.yr)
median(fire_ssp370_5_tC, na.rm=TRUE) # 0.0294

# SSP585
fire_ssp585_4_tC <- apply(fire_ssp585_4_avglast,c(1,2), tC.ha.yr)
median(fire_ssp585_4_tC, na.rm=TRUE) # 0.0619

fire_ssp585_3_tC <- apply(fire_ssp585_3_avglast,c(1,2), tC.ha.yr)
median(fire_ssp585_3_tC, na.rm=TRUE) # 0.0462

fire_ssp585_5_tC <- apply(fire_ssp585_5_avglast,c(1,2), tC.ha.yr)
median(fire_ssp585_5_tC, na.rm=TRUE) # 0.0797

# NBE
# 2014-2019
nbe_ssp126_4_avgfirst_tC <- apply(nbe_ssp126_4_avgfirst,c(1,2), tC.ha.yr)
median(nbe_ssp126_4_avgfirst_tC, na.rm=TRUE) # returns -0.0969 tC.ha-1.yr-1

nbe_ssp126_3_avgfirst_tC <- apply(nbe_ssp126_3_avgfirst,c(1,2), tC.ha.yr)
median(nbe_ssp126_3_avgfirst_tC, na.rm=TRUE) # -0.813 tC.ha-1.yr-1

nbe_ssp126_5_avgfirst_tC <- apply(nbe_ssp126_5_avgfirst,c(1,2), tC.ha.yr)
median(nbe_ssp126_5_avgfirst_tC, na.rm=TRUE) # 0.638 tC.ha-1.yr-1

# SSP126
nbe_ssp126_4_tC <- apply(nbe_ssp126_4_avglast,c(1,2), tC.ha.yr)
median(nbe_ssp126_4_tC, na.rm=TRUE) # -0.796

nbe_ssp126_3_tC <- apply(nbe_ssp126_3_avglast,c(1,2), tC.ha.yr)
median(nbe_ssp126_3_tC, na.rm=TRUE) # -1.679

nbe_ssp126_5_tC <- apply(nbe_ssp126_5_avglast,c(1,2), tC.ha.yr)
median(nbe_ssp126_5_tC, na.rm=TRUE) # 0.0514

# SSP245
nbe_ssp245_4_tC <- apply(nbe_ssp245_4_avglast,c(1,2), tC.ha.yr)
median(nbe_ssp245_4_tC, na.rm=TRUE) # -2.595

nbe_ssp245_3_tC <- apply(nbe_ssp245_3_avglast,c(1,2), tC.ha.yr)
median(nbe_ssp245_3_tC, na.rm=TRUE) # -3.926

nbe_ssp245_5_tC <- apply(nbe_ssp245_5_avglast,c(1,2), tC.ha.yr)
median(nbe_ssp245_5_tC, na.rm=TRUE) # -1.197

# SSP370
nbe_ssp370_4_tC <- apply(nbe_ssp370_4_avglast,c(1,2), tC.ha.yr)
median(nbe_ssp370_4_tC, na.rm=TRUE) # -.2.257

nbe_ssp370_3_tC <- apply(nbe_ssp370_3_avglast,c(1,2), tC.ha.yr)
median(nbe_ssp370_3_tC, na.rm=TRUE) # -3.725

nbe_ssp370_5_tC <- apply(nbe_ssp370_5_avglast,c(1,2), tC.ha.yr)
median(nbe_ssp370_5_tC, na.rm=TRUE) # -0.340

# SSP585
nbe_ssp585_4_tC <- apply(nbe_ssp585_4_avglast,c(1,2), tC.ha.yr)
median(nbe_ssp585_4_tC, na.rm=TRUE) # -5.143

nbe_ssp585_3_tC <- apply(nbe_ssp585_3_avglast,c(1,2), tC.ha.yr)
median(nbe_ssp585_3_tC, na.rm=TRUE) # -7.230

nbe_ssp585_5_tC <- apply(nbe_ssp585_5_avglast,c(1,2), tC.ha.yr)
median(nbe_ssp585_5_tC, na.rm=TRUE) # -2.777

# fLuc
# 2014-2019
fLuc_ssp126_4_avgfirst_tC <- apply(fLuc_ssp126_4_avgfirst,c(1,2), tC.ha.yr)
median(fLuc_ssp126_4_avgfirst_tC, na.rm=TRUE) # returns 0 tC.ha-1.yr-1

fLuc_ssp126_3_avgfirst_tC <- apply(fLuc_ssp126_3_avgfirst,c(1,2), tC.ha.yr)
median(fLuc_ssp126_3_avgfirst_tC, na.rm=TRUE) # 0 tC.ha-1.yr-1

fLuc_ssp126_5_avgfirst_tC <- apply(fLuc_ssp126_5_avgfirst,c(1,2), tC.ha.yr)
median(fLuc_ssp126_5_avgfirst_tC, na.rm=TRUE) # 0 tC.ha-1.yr-1

# SSP126
fLuc_ssp126_4_tC <- apply(fLuc_ssp126_4_avglast,c(1,2), tC.ha.yr)
median(fLuc_ssp126_4_tC, na.rm=TRUE) # 0.0840

fLuc_ssp126_3_tC <- apply(fLuc_ssp126_3_avglast,c(1,2), tC.ha.yr)
median(fLuc_ssp126_3_tC, na.rm=TRUE) # 0.0614

fLuc_ssp126_5_tC <- apply(fLuc_ssp126_5_avglast,c(1,2), tC.ha.yr)
median(fLuc_ssp126_5_tC, na.rm=TRUE) # 0.122

# SSP245
fLuc_ssp245_4_tC <- apply(fLuc_ssp245_4_avglast,c(1,2), tC.ha.yr)
median(fLuc_ssp245_4_tC, na.rm=TRUE) # 0.248

fLuc_ssp245_3_tC <- apply(fLuc_ssp245_3_avglast,c(1,2), tC.ha.yr)
median(fLuc_ssp245_3_tC, na.rm=TRUE) # 0.167

fLuc_ssp245_5_tC <- apply(fLuc_ssp245_5_avglast,c(1,2), tC.ha.yr)
median(fLuc_ssp245_5_tC, na.rm=TRUE) # 0.338

# SSP370
fLuc_ssp370_4_tC <- apply(fLuc_ssp370_4_avglast,c(1,2), tC.ha.yr)
median(fLuc_ssp370_4_tC, na.rm=TRUE) # 4.15e-5

fLuc_ssp370_3_tC <- apply(fLuc_ssp370_3_avglast,c(1,2), tC.ha.yr)
median(fLuc_ssp370_3_tC, na.rm=TRUE) # 1.81e-6

fLuc_ssp370_5_tC <- apply(fLuc_ssp370_5_avglast,c(1,2), tC.ha.yr)
median(fLuc_ssp370_5_tC, na.rm=TRUE) # 0.00154

# SSP585
fLuc_ssp585_4_tC <- apply(fLuc_ssp585_4_avglast,c(1,2), tC.ha.yr)
median(fLuc_ssp585_4_tC, na.rm=TRUE) # 0.439

fLuc_ssp585_3_tC <- apply(fLuc_ssp585_3_avglast,c(1,2), tC.ha.yr)
median(fLuc_ssp585_3_tC, na.rm=TRUE) # 0.366

fLuc_ssp585_5_tC <- apply(fLuc_ssp585_5_avglast,c(1,2), tC.ha.yr)
median(fLuc_ssp585_5_tC, na.rm=TRUE) # 0.561

# nbp
# 2014-2019
nbp_ssp126_4_avgfirst_tC <- apply(nbp_ssp126_4_avgfirst,c(1,2), tC.ha.yr)
median(nbp_ssp126_4_avgfirst_tC, na.rm=TRUE) # returns 0.0927 tC.ha-1.yr-1

nbp_ssp126_3_avgfirst_tC <- apply(nbp_ssp126_3_avgfirst,c(1,2), tC.ha.yr)
median(nbp_ssp126_3_avgfirst_tC, na.rm=TRUE) # -0.641 tC.ha-1.yr-1

nbp_ssp126_5_avgfirst_tC <- apply(nbp_ssp126_5_avgfirst,c(1,2), tC.ha.yr)
median(nbp_ssp126_5_avgfirst_tC, na.rm=TRUE) # 0.803 tC.ha-1.yr-1

# SSP126
nbp_ssp126_4_tC <- apply(nbp_ssp126_4_avglast,c(1,2), tC.ha.yr)
median(nbp_ssp126_4_tC, na.rm=TRUE) # 0.682

nbp_ssp126_3_tC <- apply(nbp_ssp126_3_avglast,c(1,2), tC.ha.yr)
median(nbp_ssp126_3_tC, na.rm=TRUE) # -0.196

nbp_ssp126_5_tC <- apply(nbp_ssp126_5_avglast,c(1,2), tC.ha.yr)
median(nbp_ssp126_5_tC, na.rm=TRUE) # 1.517

# SSP245
nbp_ssp245_4_tC <- apply(nbp_ssp245_4_avglast,c(1,2), tC.ha.yr)
median(nbp_ssp245_4_tC, na.rm=TRUE) # 2.183

nbp_ssp245_3_tC <- apply(nbp_ssp245_3_avglast,c(1,2), tC.ha.yr)
median(nbp_ssp245_3_tC, na.rm=TRUE) # 0.766

nbp_ssp245_5_tC <- apply(nbp_ssp245_5_avglast,c(1,2), tC.ha.yr)
median(nbp_ssp245_5_tC, na.rm=TRUE) # 3.492

# SSP370
nbp_ssp370_4_tC <- apply(nbp_ssp370_4_avglast,c(1,2), tC.ha.yr)
median(nbp_ssp370_4_tC, na.rm=TRUE) # 0.889

nbp_ssp370_3_tC <- apply(nbp_ssp370_3_avglast,c(1,2), tC.ha.yr)
median(nbp_ssp370_3_tC, na.rm=TRUE) # -0.340

nbp_ssp370_5_tC <- apply(nbp_ssp370_5_avglast,c(1,2), tC.ha.yr)
median(nbp_ssp370_5_tC, na.rm=TRUE) # 2.068

# SSP585
nbp_ssp585_4_tC <- apply(nbp_ssp585_4_avglast,c(1,2), tC.ha.yr)
median(nbp_ssp585_4_tC, na.rm=TRUE) # 4.177

nbp_ssp585_3_tC <- apply(nbp_ssp585_3_avglast,c(1,2), tC.ha.yr)
median(nbp_ssp585_3_tC, na.rm=TRUE) # 1.894

nbp_ssp585_5_tC <- apply(nbp_ssp585_5_avglast,c(1,2), tC.ha.yr)
median(nbp_ssp585_5_tC, na.rm=TRUE) # 6.348

# Creating Collums containing the medians for each of the Carbonf luxes from all scenarios
c1 <- c("GPP","Ra","NPP","Rh","NEE","Fire","NBE","FLuc","NBP")

currentmedian <- c(median(gpp_ssp126_4_avgfirst_tC, na.rm=TRUE), median(ra_ssp126_4_avgfirst_tC, na.rm=TRUE),median(npp_ssp126_4_avgfirst_tC, na.rm=TRUE),median(rh_ssp126_4_avgfirst_tC, na.rm=TRUE),median(nee_ssp126_4_avgfirst_tC, na.rm=TRUE),median(fire_ssp126_4_avgfirst_tC, na.rm=TRUE),median(nbe_ssp126_4_avgfirst_tC, na.rm=TRUE),median(fLuc_ssp126_4_avgfirst_tC, na.rm=TRUE),median(nbp_ssp126_4_avgfirst_tC, na.rm=TRUE))

current25 <- c(median(gpp_ssp126_3_avgfirst_tC, na.rm=TRUE), median(ra_ssp126_3_avgfirst_tC, na.rm=TRUE),median(npp_ssp126_3_avgfirst_tC, na.rm=TRUE),median(rh_ssp126_3_avgfirst_tC, na.rm=TRUE),median(nee_ssp126_3_avgfirst_tC, na.rm=TRUE),median(fire_ssp126_3_avgfirst_tC, na.rm=TRUE),median(nbe_ssp126_3_avgfirst_tC, na.rm=TRUE),median(fLuc_ssp126_3_avgfirst_tC, na.rm=TRUE),median(nbp_ssp126_3_avgfirst_tC, na.rm=TRUE))

current75 <- c(median(gpp_ssp126_5_avgfirst_tC, na.rm=TRUE), median(ra_ssp126_5_avgfirst_tC, na.rm=TRUE),median(npp_ssp126_5_avgfirst_tC, na.rm=TRUE),median(rh_ssp126_5_avgfirst_tC, na.rm=TRUE),median(nee_ssp126_5_avgfirst_tC, na.rm=TRUE),median(fire_ssp126_5_avgfirst_tC, na.rm=TRUE),median(nbe_ssp126_5_avgfirst_tC, na.rm=TRUE),median(fLuc_ssp126_5_avgfirst_tC, na.rm=TRUE),median(nbp_ssp126_5_avgfirst_tC, na.rm=TRUE))
 
SSP126median <- c(median(gpp_ssp126_4_tC, na.rm=TRUE), median(ra_ssp126_4_tC, na.rm=TRUE),median(npp_ssp126_4_tC, na.rm=TRUE),median(rh_ssp126_4_tC, na.rm=TRUE),median(nee_ssp126_4_tC, na.rm=TRUE),median(fire_ssp126_4_tC, na.rm=TRUE),median(nbe_ssp126_4_tC, na.rm=TRUE),median(fLuc_ssp126_4_tC, na.rm=TRUE),median(nbp_ssp126_4_tC, na.rm=TRUE))

SSP12625 <- c(median(gpp_ssp126_3_tC, na.rm=TRUE), median(ra_ssp126_3_tC, na.rm=TRUE),median(npp_ssp126_3_tC, na.rm=TRUE),median(rh_ssp126_3_tC, na.rm=TRUE),median(nee_ssp126_3_tC, na.rm=TRUE),median(fire_ssp126_3_tC, na.rm=TRUE),median(nbe_ssp126_3_tC, na.rm=TRUE),median(fLuc_ssp126_3_tC, na.rm=TRUE),median(nbp_ssp126_3_tC, na.rm=TRUE))

SSP12675 <- c(median(gpp_ssp126_5_tC, na.rm=TRUE), median(ra_ssp126_5_tC, na.rm=TRUE),median(npp_ssp126_5_tC, na.rm=TRUE),median(rh_ssp126_5_tC, na.rm=TRUE),median(nee_ssp126_5_tC, na.rm=TRUE),median(fire_ssp126_5_tC, na.rm=TRUE),median(nbe_ssp126_5_tC, na.rm=TRUE),median(fLuc_ssp126_5_tC, na.rm=TRUE),median(nbp_ssp126_5_tC, na.rm=TRUE))

SSP245median <- c(median(gpp_ssp245_4_tC, na.rm=TRUE), median(ra_ssp245_4_tC, na.rm=TRUE),median(npp_ssp245_4_tC, na.rm=TRUE),median(rh_ssp245_4_tC, na.rm=TRUE),median(nee_ssp245_4_tC, na.rm=TRUE),median(fire_ssp245_4_tC, na.rm=TRUE),median(nbe_ssp245_4_tC, na.rm=TRUE),median(fLuc_ssp245_4_tC, na.rm=TRUE),median(nbp_ssp245_4_tC, na.rm=TRUE))

SSP24525 <- c(median(gpp_ssp245_3_tC, na.rm=TRUE),median(ra_ssp245_3_tC, na.rm=TRUE),median(npp_ssp245_3_tC, na.rm=TRUE),median(rh_ssp245_3_tC, na.rm=TRUE),median(nee_ssp245_3_tC, na.rm=TRUE),median(fire_ssp245_3_tC, na.rm=TRUE),median(nbe_ssp245_3_tC, na.rm=TRUE),median(fLuc_ssp245_3_tC, na.rm=TRUE),median(nbp_ssp245_3_tC, na.rm=TRUE))

SSP24575 <- c(median(gpp_ssp245_5_tC, na.rm=TRUE), median(ra_ssp245_5_tC, na.rm=TRUE),median(npp_ssp245_5_tC, na.rm=TRUE),median(rh_ssp245_5_tC, na.rm=TRUE),median(nee_ssp245_5_tC, na.rm=TRUE),median(fire_ssp245_5_tC, na.rm=TRUE),median(nbe_ssp245_5_tC, na.rm=TRUE),median(fLuc_ssp245_5_tC, na.rm=TRUE),median(nbp_ssp245_5_tC, na.rm=TRUE))

SSP370median <- c(median(gpp_ssp370_4_tC, na.rm=TRUE), median(ra_ssp370_4_tC, na.rm=TRUE),median(npp_ssp370_4_tC, na.rm=TRUE),median(rh_ssp370_4_tC, na.rm=TRUE),median(nee_ssp370_4_tC, na.rm=TRUE),median(fire_ssp370_4_tC, na.rm=TRUE),median(nbe_ssp370_4_tC, na.rm=TRUE),median(fLuc_ssp370_4_tC, na.rm=TRUE),median(nbp_ssp370_4_tC, na.rm=TRUE))

SSP37025 <- c(median(gpp_ssp370_3_tC, na.rm=TRUE), median(ra_ssp370_3_tC, na.rm=TRUE),median(npp_ssp370_3_tC, na.rm=TRUE),median(rh_ssp370_3_tC, na.rm=TRUE),median(nee_ssp370_3_tC, na.rm=TRUE),median(fire_ssp370_3_tC, na.rm=TRUE),median(nbe_ssp370_3_tC, na.rm=TRUE),median(fLuc_ssp370_3_tC, na.rm=TRUE),median(nbp_ssp370_3_tC, na.rm=TRUE))

SSP37075 <- c(median(gpp_ssp370_5_tC, na.rm=TRUE), median(ra_ssp370_5_tC, na.rm=TRUE),median(npp_ssp370_5_tC, na.rm=TRUE),median(rh_ssp370_5_tC, na.rm=TRUE),median(nee_ssp370_5_tC, na.rm=TRUE),median(fire_ssp370_5_tC, na.rm=TRUE),median(nbe_ssp370_5_tC, na.rm=TRUE),median(fLuc_ssp370_5_tC, na.rm=TRUE),median(nbp_ssp370_5_tC, na.rm=TRUE))

SSP585median <- c(median(gpp_ssp585_4_tC, na.rm=TRUE), median(ra_ssp585_4_tC, na.rm=TRUE),median(npp_ssp585_4_tC, na.rm=TRUE),median(rh_ssp585_4_tC, na.rm=TRUE),median(nee_ssp585_4_tC, na.rm=TRUE),median(fire_ssp585_4_tC, na.rm=TRUE),median(nbe_ssp585_4_tC, na.rm=TRUE),median(fLuc_ssp585_4_tC, na.rm=TRUE),median(nbp_ssp585_4_tC, na.rm=TRUE))

SSP58525 <- c(median(gpp_ssp585_3_tC, na.rm=TRUE), median(ra_ssp585_3_tC, na.rm=TRUE),median(npp_ssp585_3_tC, na.rm=TRUE),median(rh_ssp585_3_tC, na.rm=TRUE),median(nee_ssp585_3_tC, na.rm=TRUE),median(fire_ssp585_3_tC, na.rm=TRUE),median(nbe_ssp585_3_tC, na.rm=TRUE),median(fLuc_ssp585_3_tC, na.rm=TRUE),median(nbp_ssp585_3_tC, na.rm=TRUE))

SSP58575 <- c(median(gpp_ssp585_5_tC, na.rm=TRUE), median(ra_ssp585_5_tC, na.rm=TRUE),median(npp_ssp585_5_tC, na.rm=TRUE),median(rh_ssp585_5_tC, na.rm=TRUE),median(nee_ssp585_5_tC, na.rm=TRUE),median(fire_ssp585_5_tC, na.rm=TRUE),median(nbe_ssp585_5_tC, na.rm=TRUE),median(fLuc_ssp585_5_tC, na.rm=TRUE),median(nbp_ssp585_5_tC, na.rm=TRUE))

# Changing all values to 2 significant figures
currentmedian <- signif(currentmedian,2)
current25 <- signif(current25,2)
current75 <- signif(current75,2)
SSP126median <-signif(SSP126median,2)
SSP12625 <- signif(SSP12625,2)
SSP12675 <- signif(SSP12675,2)
SSP245median <-signif(SSP245median,2)
SSP24525 <- signif(SSP24525,2)
SSP24575 <- signif(SSP24575,2)
SSP370median <-signif(SSP370median,2)
SSP37025 <- signif(SSP37025,2)
SSP37075 <- signif(SSP37075,2)
SSP585median <-signif(SSP585median,2)
SSP58525 <- signif(SSP58525,2)
SSP58575 <- signif(SSP58575,2)

x <- c(")",")",")",")",")",")",")",")",")")

FLUXES <- data.frame(c1,currentmedian,current25,current75,
SSP126median,SSP12625,SSP12675,
SSP245median,SSP24525,SSP24575,
SSP370median,SSP37025,SSP37075,
SSP585median,SSP58525,SSP58575,x)

# Colapse collums together
FLUXES$Current <- paste(FLUXES$currentmedian, FLUXES$current25, sep="(")
FLUXES$Current <- paste(FLUXES$Current, FLUXES$current75, sep="/")
FLUXES$Current <- paste(FLUXES$Current, FLUXES$x, sep="")

FLUXES$SSP126 <- paste(FLUXES$SSP126median, FLUXES$SSP12625, sep="(")
FLUXES$SSP126 <- paste(FLUXES$SSP126, FLUXES$SSP12675, sep="/")
FLUXES$SSP126 <- paste(FLUXES$SSP126, FLUXES$x, sep="")

FLUXES$SSP245 <- paste(FLUXES$SSP245median, FLUXES$SSP24525, sep="(")
FLUXES$SSP245 <- paste(FLUXES$SSP245, FLUXES$SSP24575, sep="/")
FLUXES$SSP245 <- paste(FLUXES$SSP245, FLUXES$x, sep="")

FLUXES$SSP370 <- paste(FLUXES$SSP370median, FLUXES$SSP37025, sep="(")
FLUXES$SSP370 <- paste(FLUXES$SSP370, FLUXES$SSP37075, sep="/")
FLUXES$SSP370 <- paste(FLUXES$SSP370, FLUXES$x, sep="")

FLUXES$SSP585 <- paste(FLUXES$SSP585median, FLUXES$SSP58525, sep="(")
FLUXES$SSP585 <- paste(FLUXES$SSP585, FLUXES$SSP58575, sep="/")
FLUXES$SSP585 <- paste(FLUXES$SSP585, FLUXES$x, sep="")

# Delete unecessary collumns 
FLUXES <- subset (FLUXES, select = -c(currentmedian:x)) 

kbl(FLUXES) %>%
  kable_classic() 

