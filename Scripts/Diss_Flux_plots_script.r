library(ncdf4)
library(fields)
library(dplyr)
library(lubridate)
library(raster)
library(zoo)
library(ggplot2)

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

# When done with reading in variables remember to close the file
nc_close(cflux_ssp126) ; nc_close(cflux_ssp245) ; nc_close(cflux_ssp370) ; nc_close(cflux_ssp585) 

# Slicing to get median and 90% confidence intervals.
# Calculating FLUX values
# Avg the last 5 Years of observed data (2014 -2019) vs the las 5 Year of the the porjections
# gpp median
gpp_ssp126_4 <- gpp_ssp126[,,4,]
gpp_ssp126_4 <- gpp_ssp126[,,4,]
gpp_ssp245_4 <- gpp_ssp245[,,4,]
gpp_ssp370_4 <- gpp_ssp370[,,4,]
gpp_ssp585_4 <- gpp_ssp585[,,4,]

# gpp 25%
gpp_ssp126_1 <- gpp_ssp126[,,1,]
gpp_ssp126_1 <- gpp_ssp126[,,1,]
gpp_ssp245_1 <- gpp_ssp245[,,1,]
gpp_ssp370_1 <- gpp_ssp370[,,1,]
gpp_ssp585_1 <- gpp_ssp585[,,1,]

# gpp 75%
gpp_ssp126_7 <- gpp_ssp126[,,7,]
gpp_ssp126_7 <- gpp_ssp126[,,7,]
gpp_ssp245_7 <- gpp_ssp245[,,7,]
gpp_ssp370_7 <- gpp_ssp370[,,7,]
gpp_ssp585_7 <- gpp_ssp585[,,7,]

# ra median
ra_ssp126_4 <- ra_ssp126[,,4,]
ra_ssp126_4 <- ra_ssp126[,,4,]
ra_ssp245_4 <- ra_ssp245[,,4,]
ra_ssp370_4 <- ra_ssp370[,,4,]
ra_ssp585_4 <- ra_ssp585[,,4,]

# ra 25%
ra_ssp126_1 <- ra_ssp126[,,1,]
ra_ssp126_1 <- ra_ssp126[,,1,]
ra_ssp245_1 <- ra_ssp245[,,1,]
ra_ssp370_1 <- ra_ssp370[,,1,]
ra_ssp585_1 <- ra_ssp585[,,1,]

# ra 75%
ra_ssp126_7 <- ra_ssp126[,,7,]
ra_ssp126_7 <- ra_ssp126[,,7,]
ra_ssp245_7 <- ra_ssp245[,,7,]
ra_ssp370_7 <- ra_ssp370[,,7,]
ra_ssp585_7 <- ra_ssp585[,,7,]

# npp median

npp_ssp126_4 <- npp_ssp126[,,4,]
npp_ssp126_4 <- npp_ssp126[,,4,]
npp_ssp245_4 <- npp_ssp245[,,4,]
npp_ssp370_4 <- npp_ssp370[,,4,]
npp_ssp585_4 <- npp_ssp585[,,4,]

# npp 25%
npp_ssp126_1 <- npp_ssp126[,,1,]
npp_ssp126_1 <- npp_ssp126[,,1,]
npp_ssp245_1 <- npp_ssp245[,,1,]
npp_ssp370_1 <- npp_ssp370[,,1,]
npp_ssp585_1 <- npp_ssp585[,,1,]

# npp 75%
npp_ssp126_7 <- npp_ssp126[,,7,]
npp_ssp126_7 <- npp_ssp126[,,7,]
npp_ssp245_7 <- npp_ssp245[,,7,]
npp_ssp370_7 <- npp_ssp370[,,7,]
npp_ssp585_7 <- npp_ssp585[,,7,]

# rh median
rh_ssp126_4 <- rh_ssp126[,,4,]
rh_ssp126_4 <- rh_ssp126[,,4,]
rh_ssp245_4 <- rh_ssp245[,,4,]
rh_ssp370_4 <- rh_ssp370[,,4,]
rh_ssp585_4 <- rh_ssp585[,,4,]

# rh 25%
rh_ssp126_1 <- rh_ssp126[,,1,]
rh_ssp126_1 <- rh_ssp126[,,1,]
rh_ssp245_1 <- rh_ssp245[,,1,]
rh_ssp370_1 <- rh_ssp370[,,1,]
rh_ssp585_1 <- rh_ssp585[,,1,]

# rh 75%
rh_ssp126_7 <- rh_ssp126[,,7,]
rh_ssp126_7 <- rh_ssp126[,,7,]
rh_ssp245_7 <- rh_ssp245[,,7,]
rh_ssp370_7 <- rh_ssp370[,,7,]
rh_ssp585_7 <- rh_ssp585[,,7,]

# nee median
nee_ssp126_4 <- nee_ssp126[,,4,]
nee_ssp126_4 <- nee_ssp126[,,4,]
nee_ssp245_4 <- nee_ssp245[,,4,]
nee_ssp370_4 <- nee_ssp370[,,4,]
nee_ssp585_4 <- nee_ssp585[,,4,]

# nee 25%
nee_ssp126_1 <- nee_ssp126[,,1,]
nee_ssp126_1 <- nee_ssp126[,,1,]
nee_ssp245_1 <- nee_ssp245[,,1,]
nee_ssp370_1 <- nee_ssp370[,,1,]
nee_ssp585_1 <- nee_ssp585[,,1,]

# nee 75%
nee_ssp126_7 <- nee_ssp126[,,7,]
nee_ssp126_7 <- nee_ssp126[,,7,]
nee_ssp245_7 <- nee_ssp245[,,7,]
nee_ssp370_7 <- nee_ssp370[,,7,]
nee_ssp585_7 <- nee_ssp585[,,7,]

# fire median
fire_ssp126_4 <- fire_ssp126[,,4,]
fire_ssp126_4 <- fire_ssp126[,,4,]
fire_ssp245_4 <- fire_ssp245[,,4,]
fire_ssp370_4 <- fire_ssp370[,,4,]
fire_ssp585_4 <- fire_ssp585[,,4,]

# fire 25%
fire_ssp126_1 <- fire_ssp126[,,1,]
fire_ssp126_1 <- fire_ssp126[,,1,]
fire_ssp245_1 <- fire_ssp245[,,1,]
fire_ssp370_1 <- fire_ssp370[,,1,]
fire_ssp585_1 <- fire_ssp585[,,1,]

# fire 75%
fire_ssp126_7 <- fire_ssp126[,,7,]
fire_ssp126_7 <- fire_ssp126[,,7,]
fire_ssp245_7 <- fire_ssp245[,,7,]
fire_ssp370_7 <- fire_ssp370[,,7,]
fire_ssp585_7 <- fire_ssp585[,,7,]

# nbe median
nbe_ssp126_4 <- nbe_ssp126[,,4,]
nbe_ssp126_4 <- nbe_ssp126[,,4,]
nbe_ssp245_4 <- nbe_ssp245[,,4,]
nbe_ssp370_4 <- nbe_ssp370[,,4,]
nbe_ssp585_4 <- nbe_ssp585[,,4,]

# nbe 25%
nbe_ssp126_1 <- nbe_ssp126[,,1,]
nbe_ssp126_1 <- nbe_ssp126[,,1,]
nbe_ssp245_1 <- nbe_ssp245[,,1,]
nbe_ssp370_1 <- nbe_ssp370[,,1,]
nbe_ssp585_1 <- nbe_ssp585[,,1,]

# nbe 75%
nbe_ssp126_7 <- nbe_ssp126[,,7,]
nbe_ssp126_7 <- nbe_ssp126[,,7,]
nbe_ssp245_7 <- nbe_ssp245[,,7,]
nbe_ssp370_7 <- nbe_ssp370[,,7,]
nbe_ssp585_7 <- nbe_ssp585[,,7,]

# fLuc median
fLuc_ssp126_4 <- fLuc_ssp126[,,4,]
fLuc_ssp126_4 <- fLuc_ssp126[,,4,]
fLuc_ssp245_4 <- fLuc_ssp245[,,4,]
fLuc_ssp370_4 <- fLuc_ssp370[,,4,]
fLuc_ssp585_4 <- fLuc_ssp585[,,4,]

# fLuc 25%
fLuc_ssp126_1 <- fLuc_ssp126[,,1,]
fLuc_ssp126_1 <- fLuc_ssp126[,,1,]
fLuc_ssp245_1 <- fLuc_ssp245[,,1,]
fLuc_ssp370_1 <- fLuc_ssp370[,,1,]
fLuc_ssp585_1 <- fLuc_ssp585[,,1,]

# fLuc 75%
fLuc_ssp126_7 <- fLuc_ssp126[,,7,]
fLuc_ssp126_7 <- fLuc_ssp126[,,7,]
fLuc_ssp245_7 <- fLuc_ssp245[,,7,]
fLuc_ssp370_7 <- fLuc_ssp370[,,7,]
fLuc_ssp585_7 <- fLuc_ssp585[,,7,]

# nbp median
nbp_ssp126_4 <- nbp_ssp126[,,4,]
nbp_ssp126_4 <- nbp_ssp126[,,4,]
nbp_ssp245_4 <- nbp_ssp245[,,4,]
nbp_ssp370_4 <- nbp_ssp370[,,4,]
nbp_ssp585_4 <- nbp_ssp585[,,4,]

# nbp 25%
nbp_ssp126_1 <- nbp_ssp126[,,1,]
nbp_ssp126_1 <- nbp_ssp126[,,1,]
nbp_ssp245_1 <- nbp_ssp245[,,1,]
nbp_ssp370_1 <- nbp_ssp370[,,1,]
nbp_ssp585_1 <- nbp_ssp585[,,1,]

# nbp 75%
nbp_ssp126_7 <- nbp_ssp126[,,7,]
nbp_ssp126_7 <- nbp_ssp126[,,7,]
nbp_ssp245_7 <- nbp_ssp245[,,7,]
nbp_ssp370_7 <- nbp_ssp370[,,7,]
nbp_ssp585_7 <- nbp_ssp585[,,7,]

# How many time steps per Year?
steps_per_Year = 12

# create new object which will contain the aggregated to mean annual rather than monthly

# gpp
# Converting gpp to annual ssp126
gpp_ssp126_4_annual = array(NA, dim=c(dim(gpp_ssp126_4)[1],dim(gpp_ssp126_4)[2],dim(gpp_ssp126_4)[3]/steps_per_Year))
for (i in seq(1,dim(gpp_ssp126_4)[1])) {
     for (j in seq(1,dim(gpp_ssp126_4)[2])) {
          gpp_ssp126_4_annual[i,j,] = rollapply(gpp_ssp126_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

gpp_ssp126_1_annual = array(NA, dim=c(dim(gpp_ssp126_1)[1],dim(gpp_ssp126_4)[2],dim(gpp_ssp126_1)[3]/steps_per_Year))
for (i in seq(1,dim(gpp_ssp126_1)[1])) {
     for (j in seq(1,dim(gpp_ssp126_1)[2])) {
          gpp_ssp126_1_annual[i,j,] = rollapply(gpp_ssp126_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

gpp_ssp126_7_annual = array(NA, dim=c(dim(gpp_ssp126_1)[1],dim(gpp_ssp126_7)[2],dim(gpp_ssp126_7)[3]/steps_per_Year))
for (i in seq(1,dim(gpp_ssp126_7)[1])) {
     for (j in seq(1,dim(gpp_ssp126_7)[2])) {
          gpp_ssp126_7_annual[i,j,] = rollapply(gpp_ssp126_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting gpp to annual ssp245
gpp_ssp245_4_annual = array(NA, dim=c(dim(gpp_ssp245_4)[1],dim(gpp_ssp245_4)[2],dim(gpp_ssp245_4)[3]/steps_per_Year))
for (i in seq(1,dim(gpp_ssp245_4)[1])) {
     for (j in seq(1,dim(gpp_ssp245_4)[2])) {
          gpp_ssp245_4_annual[i,j,] = rollapply(gpp_ssp245_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

gpp_ssp245_1_annual = array(NA, dim=c(dim(gpp_ssp245_1)[1],dim(gpp_ssp245_4)[2],dim(gpp_ssp245_1)[3]/steps_per_Year))
for (i in seq(1,dim(gpp_ssp245_1)[1])) {
     for (j in seq(1,dim(gpp_ssp245_1)[2])) {
          gpp_ssp245_1_annual[i,j,] = rollapply(gpp_ssp245_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

gpp_ssp245_7_annual = array(NA, dim=c(dim(gpp_ssp245_7)[1],dim(gpp_ssp245_7)[2],dim(gpp_ssp245_7)[3]/steps_per_Year))
for (i in seq(1,dim(gpp_ssp245_1)[1])) {
     for (j in seq(1,dim(gpp_ssp245_1)[2])) {
          gpp_ssp245_7_annual[i,j,] = rollapply(gpp_ssp245_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting gpp to annual ssp370
gpp_ssp370_4_annual = array(NA, dim=c(dim(gpp_ssp370_4)[1],dim(gpp_ssp370_4)[2],dim(gpp_ssp370_4)[3]/steps_per_Year))
for (i in seq(1,dim(gpp_ssp370_4)[1])) {
     for (j in seq(1,dim(gpp_ssp370_4)[2])) {
          gpp_ssp370_4_annual[i,j,] = rollapply(gpp_ssp370_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

gpp_ssp370_1_annual = array(NA, dim=c(dim(gpp_ssp370_1)[1],dim(gpp_ssp370_4)[2],dim(gpp_ssp370_1)[3]/steps_per_Year))
for (i in seq(1,dim(gpp_ssp370_1)[1])) {
     for (j in seq(1,dim(gpp_ssp370_1)[2])) {
          gpp_ssp370_1_annual[i,j,] = rollapply(gpp_ssp370_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

gpp_ssp370_7_annual = array(NA, dim=c(dim(gpp_ssp370_7)[1],dim(gpp_ssp370_7)[2],dim(gpp_ssp370_7)[3]/steps_per_Year))
for (i in seq(1,dim(gpp_ssp370_7)[1])) {
     for (j in seq(1,dim(gpp_ssp370_7)[2])) {
          gpp_ssp370_7_annual[i,j,] = rollapply(gpp_ssp370_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting gpp to annual ssp585
gpp_ssp585_4_annual = array(NA, dim=c(dim(gpp_ssp585_4)[1],dim(gpp_ssp585_4)[2],dim(gpp_ssp585_4)[3]/steps_per_Year))
for (i in seq(1,dim(gpp_ssp585_4)[1])) {
     for (j in seq(1,dim(gpp_ssp585_4)[2])) {
          gpp_ssp585_4_annual[i,j,] = rollapply(gpp_ssp585_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

gpp_ssp585_1_annual = array(NA, dim=c(dim(gpp_ssp585_1)[1],dim(gpp_ssp585_4)[2],dim(gpp_ssp585_1)[3]/steps_per_Year))
for (i in seq(1,dim(gpp_ssp585_1)[1])) {
     for (j in seq(1,dim(gpp_ssp585_1)[2])) {
          gpp_ssp585_1_annual[i,j,] = rollapply(gpp_ssp585_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

gpp_ssp585_7_annual = array(NA, dim=c(dim(gpp_ssp585_7)[1],dim(gpp_ssp585_7)[2],dim(gpp_ssp585_7)[3]/steps_per_Year))
for (i in seq(1,dim(gpp_ssp585_7)[1])) {
     for (j in seq(1,dim(gpp_ssp585_7)[2])) {
          gpp_ssp585_7_annual[i,j,] = rollapply(gpp_ssp585_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Ra 
# Converting ra to annual ssp126
ra_ssp126_4_annual = array(NA, dim=c(dim(ra_ssp126_4)[1],dim(ra_ssp126_4)[2],dim(ra_ssp126_4)[3]/steps_per_Year))
for (i in seq(1,dim(ra_ssp126_4)[1])) {
     for (j in seq(1,dim(ra_ssp126_4)[2])) {
          ra_ssp126_4_annual[i,j,] = rollapply(ra_ssp126_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

ra_ssp126_1_annual = array(NA, dim=c(dim(ra_ssp126_1)[1],dim(ra_ssp126_4)[2],dim(ra_ssp126_1)[3]/steps_per_Year))
for (i in seq(1,dim(ra_ssp126_1)[1])) {
     for (j in seq(1,dim(ra_ssp126_1)[2])) {
          ra_ssp126_1_annual[i,j,] = rollapply(ra_ssp126_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

ra_ssp126_7_annual = array(NA, dim=c(dim(ra_ssp126_1)[1],dim(ra_ssp126_7)[2],dim(ra_ssp126_7)[3]/steps_per_Year))
for (i in seq(1,dim(ra_ssp126_7)[1])) {
     for (j in seq(1,dim(ra_ssp126_7)[2])) {
          ra_ssp126_7_annual[i,j,] = rollapply(ra_ssp126_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting ra to annual ssp245
ra_ssp245_4_annual = array(NA, dim=c(dim(ra_ssp245_4)[1],dim(ra_ssp245_4)[2],dim(ra_ssp245_4)[3]/steps_per_Year))
for (i in seq(1,dim(ra_ssp245_4)[1])) {
     for (j in seq(1,dim(ra_ssp245_4)[2])) {
          ra_ssp245_4_annual[i,j,] = rollapply(ra_ssp245_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

ra_ssp245_1_annual = array(NA, dim=c(dim(ra_ssp245_1)[1],dim(ra_ssp245_4)[2],dim(ra_ssp245_1)[3]/steps_per_Year))
for (i in seq(1,dim(ra_ssp245_1)[1])) {
     for (j in seq(1,dim(ra_ssp245_1)[2])) {
          ra_ssp245_1_annual[i,j,] = rollapply(ra_ssp245_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

ra_ssp245_7_annual = array(NA, dim=c(dim(ra_ssp245_7)[1],dim(ra_ssp245_7)[2],dim(ra_ssp245_7)[3]/steps_per_Year))
for (i in seq(1,dim(ra_ssp245_1)[1])) {
     for (j in seq(1,dim(ra_ssp245_1)[2])) {
          ra_ssp245_7_annual[i,j,] = rollapply(ra_ssp245_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting ra to annual ssp370
ra_ssp370_4_annual = array(NA, dim=c(dim(ra_ssp370_4)[1],dim(ra_ssp370_4)[2],dim(ra_ssp370_4)[3]/steps_per_Year))
for (i in seq(1,dim(ra_ssp370_4)[1])) {
     for (j in seq(1,dim(ra_ssp370_4)[2])) {
          ra_ssp370_4_annual[i,j,] = rollapply(ra_ssp370_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

ra_ssp370_1_annual = array(NA, dim=c(dim(ra_ssp370_1)[1],dim(ra_ssp370_4)[2],dim(ra_ssp370_1)[3]/steps_per_Year))
for (i in seq(1,dim(ra_ssp370_1)[1])) {
     for (j in seq(1,dim(ra_ssp370_1)[2])) {
          ra_ssp370_1_annual[i,j,] = rollapply(ra_ssp370_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

ra_ssp370_7_annual = array(NA, dim=c(dim(ra_ssp370_7)[1],dim(ra_ssp370_7)[2],dim(ra_ssp370_7)[3]/steps_per_Year))
for (i in seq(1,dim(ra_ssp370_7)[1])) {
     for (j in seq(1,dim(ra_ssp370_7)[2])) {
          ra_ssp370_7_annual[i,j,] = rollapply(ra_ssp370_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting ra to annual ssp585
ra_ssp585_4_annual = array(NA, dim=c(dim(ra_ssp585_4)[1],dim(ra_ssp585_4)[2],dim(ra_ssp585_4)[3]/steps_per_Year))
for (i in seq(1,dim(ra_ssp585_4)[1])) {
     for (j in seq(1,dim(ra_ssp585_4)[2])) {
          ra_ssp585_4_annual[i,j,] = rollapply(ra_ssp585_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

ra_ssp585_1_annual = array(NA, dim=c(dim(ra_ssp585_1)[1],dim(ra_ssp585_4)[2],dim(ra_ssp585_1)[3]/steps_per_Year))
for (i in seq(1,dim(ra_ssp585_1)[1])) {
     for (j in seq(1,dim(ra_ssp585_1)[2])) {
          ra_ssp585_1_annual[i,j,] = rollapply(ra_ssp585_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

ra_ssp585_7_annual = array(NA, dim=c(dim(ra_ssp585_7)[1],dim(ra_ssp585_7)[2],dim(ra_ssp585_7)[3]/steps_per_Year))
for (i in seq(1,dim(ra_ssp585_7)[1])) {
     for (j in seq(1,dim(ra_ssp585_7)[2])) {
          ra_ssp585_7_annual[i,j,] = rollapply(ra_ssp585_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

#  npp
# Converting npp to annual ssp126
npp_ssp126_4_annual = array(NA, dim=c(dim(npp_ssp126_4)[1],dim(npp_ssp126_4)[2],dim(npp_ssp126_4)[3]/steps_per_Year))
for (i in seq(1,dim(npp_ssp126_4)[1])) {
     for (j in seq(1,dim(npp_ssp126_4)[2])) {
          npp_ssp126_4_annual[i,j,] = rollapply(npp_ssp126_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

npp_ssp126_1_annual = array(NA, dim=c(dim(npp_ssp126_1)[1],dim(npp_ssp126_4)[2],dim(npp_ssp126_1)[3]/steps_per_Year))
for (i in seq(1,dim(npp_ssp126_1)[1])) {
     for (j in seq(1,dim(npp_ssp126_1)[2])) {
          npp_ssp126_1_annual[i,j,] = rollapply(npp_ssp126_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

npp_ssp126_7_annual = array(NA, dim=c(dim(npp_ssp126_1)[1],dim(npp_ssp126_7)[2],dim(npp_ssp126_7)[3]/steps_per_Year))
for (i in seq(1,dim(npp_ssp126_7)[1])) {
     for (j in seq(1,dim(npp_ssp126_7)[2])) {
          npp_ssp126_7_annual[i,j,] = rollapply(npp_ssp126_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting npp to annual ssp245
npp_ssp245_4_annual = array(NA, dim=c(dim(npp_ssp245_4)[1],dim(npp_ssp245_4)[2],dim(npp_ssp245_4)[3]/steps_per_Year))
for (i in seq(1,dim(npp_ssp245_4)[1])) {
     for (j in seq(1,dim(npp_ssp245_4)[2])) {
          npp_ssp245_4_annual[i,j,] = rollapply(npp_ssp245_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

npp_ssp245_1_annual = array(NA, dim=c(dim(npp_ssp245_1)[1],dim(npp_ssp245_4)[2],dim(npp_ssp245_1)[3]/steps_per_Year))
for (i in seq(1,dim(npp_ssp245_1)[1])) {
     for (j in seq(1,dim(npp_ssp245_1)[2])) {
          npp_ssp245_1_annual[i,j,] = rollapply(npp_ssp245_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

npp_ssp245_7_annual = array(NA, dim=c(dim(npp_ssp245_7)[1],dim(npp_ssp245_7)[2],dim(npp_ssp245_7)[3]/steps_per_Year))
for (i in seq(1,dim(npp_ssp245_1)[1])) {
     for (j in seq(1,dim(npp_ssp245_1)[2])) {
          npp_ssp245_7_annual[i,j,] = rollapply(npp_ssp245_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting npp to annual ssp370
npp_ssp370_4_annual = array(NA, dim=c(dim(npp_ssp370_4)[1],dim(npp_ssp370_4)[2],dim(npp_ssp370_4)[3]/steps_per_Year))
for (i in seq(1,dim(npp_ssp370_4)[1])) {
     for (j in seq(1,dim(npp_ssp370_4)[2])) {
          npp_ssp370_4_annual[i,j,] = rollapply(npp_ssp370_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

npp_ssp370_1_annual = array(NA, dim=c(dim(npp_ssp370_1)[1],dim(npp_ssp370_4)[2],dim(npp_ssp370_1)[3]/steps_per_Year))
for (i in seq(1,dim(npp_ssp370_1)[1])) {
     for (j in seq(1,dim(npp_ssp370_1)[2])) {
          npp_ssp370_1_annual[i,j,] = rollapply(npp_ssp370_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

npp_ssp370_7_annual = array(NA, dim=c(dim(npp_ssp370_7)[1],dim(npp_ssp370_7)[2],dim(npp_ssp370_7)[3]/steps_per_Year))
for (i in seq(1,dim(npp_ssp370_7)[1])) {
     for (j in seq(1,dim(npp_ssp370_7)[2])) {
          npp_ssp370_7_annual[i,j,] = rollapply(npp_ssp370_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting npp to annual ssp585
npp_ssp585_4_annual = array(NA, dim=c(dim(npp_ssp585_4)[1],dim(npp_ssp585_4)[2],dim(npp_ssp585_4)[3]/steps_per_Year))
for (i in seq(1,dim(npp_ssp585_4)[1])) {
     for (j in seq(1,dim(npp_ssp585_4)[2])) {
          npp_ssp585_4_annual[i,j,] = rollapply(npp_ssp585_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

npp_ssp585_1_annual = array(NA, dim=c(dim(npp_ssp585_1)[1],dim(npp_ssp585_4)[2],dim(npp_ssp585_1)[3]/steps_per_Year))
for (i in seq(1,dim(npp_ssp585_1)[1])) {
     for (j in seq(1,dim(npp_ssp585_1)[2])) {
          npp_ssp585_1_annual[i,j,] = rollapply(npp_ssp585_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

npp_ssp585_7_annual = array(NA, dim=c(dim(npp_ssp585_7)[1],dim(npp_ssp585_7)[2],dim(npp_ssp585_7)[3]/steps_per_Year))
for (i in seq(1,dim(npp_ssp585_7)[1])) {
     for (j in seq(1,dim(npp_ssp585_7)[2])) {
          npp_ssp585_7_annual[i,j,] = rollapply(npp_ssp585_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# rh
# Converting rh to annual ssp126
rh_ssp126_4_annual = array(NA, dim=c(dim(rh_ssp126_4)[1],dim(rh_ssp126_4)[2],dim(rh_ssp126_4)[3]/steps_per_Year))
for (i in seq(1,dim(rh_ssp126_4)[1])) {
     for (j in seq(1,dim(rh_ssp126_4)[2])) {
          rh_ssp126_4_annual[i,j,] = rollapply(rh_ssp126_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

rh_ssp126_1_annual = array(NA, dim=c(dim(rh_ssp126_1)[1],dim(rh_ssp126_4)[2],dim(rh_ssp126_1)[3]/steps_per_Year))
for (i in seq(1,dim(rh_ssp126_1)[1])) {
     for (j in seq(1,dim(rh_ssp126_1)[2])) {
          rh_ssp126_1_annual[i,j,] = rollapply(rh_ssp126_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

rh_ssp126_7_annual = array(NA, dim=c(dim(rh_ssp126_1)[1],dim(rh_ssp126_7)[2],dim(rh_ssp126_7)[3]/steps_per_Year))
for (i in seq(1,dim(rh_ssp126_7)[1])) {
     for (j in seq(1,dim(rh_ssp126_7)[2])) {
          rh_ssp126_7_annual[i,j,] = rollapply(rh_ssp126_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting rh to annual ssp245
rh_ssp245_4_annual = array(NA, dim=c(dim(rh_ssp245_4)[1],dim(rh_ssp245_4)[2],dim(rh_ssp245_4)[3]/steps_per_Year))
for (i in seq(1,dim(rh_ssp245_4)[1])) {
     for (j in seq(1,dim(rh_ssp245_4)[2])) {
          rh_ssp245_4_annual[i,j,] = rollapply(rh_ssp245_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

rh_ssp245_1_annual = array(NA, dim=c(dim(rh_ssp245_1)[1],dim(rh_ssp245_4)[2],dim(rh_ssp245_1)[3]/steps_per_Year))
for (i in seq(1,dim(rh_ssp245_1)[1])) {
     for (j in seq(1,dim(rh_ssp245_1)[2])) {
          rh_ssp245_1_annual[i,j,] = rollapply(rh_ssp245_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

rh_ssp245_7_annual = array(NA, dim=c(dim(rh_ssp245_7)[1],dim(rh_ssp245_7)[2],dim(rh_ssp245_7)[3]/steps_per_Year))
for (i in seq(1,dim(rh_ssp245_1)[1])) {
     for (j in seq(1,dim(rh_ssp245_1)[2])) {
          rh_ssp245_7_annual[i,j,] = rollapply(rh_ssp245_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting rh to annual ssp370
rh_ssp370_4_annual = array(NA, dim=c(dim(rh_ssp370_4)[1],dim(rh_ssp370_4)[2],dim(rh_ssp370_4)[3]/steps_per_Year))
for (i in seq(1,dim(rh_ssp370_4)[1])) {
     for (j in seq(1,dim(rh_ssp370_4)[2])) {
          rh_ssp370_4_annual[i,j,] = rollapply(rh_ssp370_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

rh_ssp370_1_annual = array(NA, dim=c(dim(rh_ssp370_1)[1],dim(rh_ssp370_4)[2],dim(rh_ssp370_1)[3]/steps_per_Year))
for (i in seq(1,dim(rh_ssp370_1)[1])) {
     for (j in seq(1,dim(rh_ssp370_1)[2])) {
          rh_ssp370_1_annual[i,j,] = rollapply(rh_ssp370_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

rh_ssp370_7_annual = array(NA, dim=c(dim(rh_ssp370_7)[1],dim(rh_ssp370_7)[2],dim(rh_ssp370_7)[3]/steps_per_Year))
for (i in seq(1,dim(rh_ssp370_7)[1])) {
     for (j in seq(1,dim(rh_ssp370_7)[2])) {
          rh_ssp370_7_annual[i,j,] = rollapply(rh_ssp370_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting rh to annual ssp585
rh_ssp585_4_annual = array(NA, dim=c(dim(rh_ssp585_4)[1],dim(rh_ssp585_4)[2],dim(rh_ssp585_4)[3]/steps_per_Year))
for (i in seq(1,dim(rh_ssp585_4)[1])) {
     for (j in seq(1,dim(rh_ssp585_4)[2])) {
          rh_ssp585_4_annual[i,j,] = rollapply(rh_ssp585_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

rh_ssp585_1_annual = array(NA, dim=c(dim(rh_ssp585_1)[1],dim(rh_ssp585_4)[2],dim(rh_ssp585_1)[3]/steps_per_Year))
for (i in seq(1,dim(rh_ssp585_1)[1])) {
     for (j in seq(1,dim(rh_ssp585_1)[2])) {
          rh_ssp585_1_annual[i,j,] = rollapply(rh_ssp585_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

rh_ssp585_7_annual = array(NA, dim=c(dim(rh_ssp585_7)[1],dim(rh_ssp585_7)[2],dim(rh_ssp585_7)[3]/steps_per_Year))
for (i in seq(1,dim(rh_ssp585_7)[1])) {
     for (j in seq(1,dim(rh_ssp585_7)[2])) {
          rh_ssp585_7_annual[i,j,] = rollapply(rh_ssp585_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# nee
# Converting nee to annual ssp126
nee_ssp126_4_annual = array(NA, dim=c(dim(nee_ssp126_4)[1],dim(nee_ssp126_4)[2],dim(nee_ssp126_4)[3]/steps_per_Year))
for (i in seq(1,dim(nee_ssp126_4)[1])) {
     for (j in seq(1,dim(nee_ssp126_4)[2])) {
          nee_ssp126_4_annual[i,j,] = rollapply(nee_ssp126_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nee_ssp126_1_annual = array(NA, dim=c(dim(nee_ssp126_1)[1],dim(nee_ssp126_4)[2],dim(nee_ssp126_1)[3]/steps_per_Year))
for (i in seq(1,dim(nee_ssp126_1)[1])) {
     for (j in seq(1,dim(nee_ssp126_1)[2])) {
          nee_ssp126_1_annual[i,j,] = rollapply(nee_ssp126_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nee_ssp126_7_annual = array(NA, dim=c(dim(nee_ssp126_1)[1],dim(nee_ssp126_7)[2],dim(nee_ssp126_7)[3]/steps_per_Year))
for (i in seq(1,dim(nee_ssp126_7)[1])) {
     for (j in seq(1,dim(nee_ssp126_7)[2])) {
          nee_ssp126_7_annual[i,j,] = rollapply(nee_ssp126_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting nee to annual ssp245
nee_ssp245_4_annual = array(NA, dim=c(dim(nee_ssp245_4)[1],dim(nee_ssp245_4)[2],dim(nee_ssp245_4)[3]/steps_per_Year))
for (i in seq(1,dim(nee_ssp245_4)[1])) {
     for (j in seq(1,dim(nee_ssp245_4)[2])) {
          nee_ssp245_4_annual[i,j,] = rollapply(nee_ssp245_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nee_ssp245_1_annual = array(NA, dim=c(dim(nee_ssp245_1)[1],dim(nee_ssp245_4)[2],dim(nee_ssp245_1)[3]/steps_per_Year))
for (i in seq(1,dim(nee_ssp245_1)[1])) {
     for (j in seq(1,dim(nee_ssp245_1)[2])) {
          nee_ssp245_1_annual[i,j,] = rollapply(nee_ssp245_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nee_ssp245_7_annual = array(NA, dim=c(dim(nee_ssp245_7)[1],dim(nee_ssp245_7)[2],dim(nee_ssp245_7)[3]/steps_per_Year))
for (i in seq(1,dim(nee_ssp245_1)[1])) {
     for (j in seq(1,dim(nee_ssp245_1)[2])) {
          nee_ssp245_7_annual[i,j,] = rollapply(nee_ssp245_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting nee to annual ssp370
nee_ssp370_4_annual = array(NA, dim=c(dim(nee_ssp370_4)[1],dim(nee_ssp370_4)[2],dim(nee_ssp370_4)[3]/steps_per_Year))
for (i in seq(1,dim(nee_ssp370_4)[1])) {
     for (j in seq(1,dim(nee_ssp370_4)[2])) {
          nee_ssp370_4_annual[i,j,] = rollapply(nee_ssp370_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nee_ssp370_1_annual = array(NA, dim=c(dim(nee_ssp370_1)[1],dim(nee_ssp370_4)[2],dim(nee_ssp370_1)[3]/steps_per_Year))
for (i in seq(1,dim(nee_ssp370_1)[1])) {
     for (j in seq(1,dim(nee_ssp370_1)[2])) {
          nee_ssp370_1_annual[i,j,] = rollapply(nee_ssp370_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nee_ssp370_7_annual = array(NA, dim=c(dim(nee_ssp370_7)[1],dim(nee_ssp370_7)[2],dim(nee_ssp370_7)[3]/steps_per_Year))
for (i in seq(1,dim(nee_ssp370_7)[1])) {
     for (j in seq(1,dim(nee_ssp370_7)[2])) {
          nee_ssp370_7_annual[i,j,] = rollapply(nee_ssp370_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting nee to annual ssp585
nee_ssp585_4_annual = array(NA, dim=c(dim(nee_ssp585_4)[1],dim(nee_ssp585_4)[2],dim(nee_ssp585_4)[3]/steps_per_Year))
for (i in seq(1,dim(nee_ssp585_4)[1])) {
     for (j in seq(1,dim(nee_ssp585_4)[2])) {
          nee_ssp585_4_annual[i,j,] = rollapply(nee_ssp585_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nee_ssp585_1_annual = array(NA, dim=c(dim(nee_ssp585_1)[1],dim(nee_ssp585_4)[2],dim(nee_ssp585_1)[3]/steps_per_Year))
for (i in seq(1,dim(nee_ssp585_1)[1])) {
     for (j in seq(1,dim(nee_ssp585_1)[2])) {
          nee_ssp585_1_annual[i,j,] = rollapply(nee_ssp585_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nee_ssp585_7_annual = array(NA, dim=c(dim(nee_ssp585_7)[1],dim(nee_ssp585_7)[2],dim(nee_ssp585_7)[3]/steps_per_Year))
for (i in seq(1,dim(nee_ssp585_7)[1])) {
     for (j in seq(1,dim(nee_ssp585_7)[2])) {
          nee_ssp585_7_annual[i,j,] = rollapply(nee_ssp585_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# fire
# Converting fire to annual ssp126
fire_ssp126_4_annual = array(NA, dim=c(dim(fire_ssp126_4)[1],dim(fire_ssp126_4)[2],dim(fire_ssp126_4)[3]/steps_per_Year))
for (i in seq(1,dim(fire_ssp126_4)[1])) {
     for (j in seq(1,dim(fire_ssp126_4)[2])) {
          fire_ssp126_4_annual[i,j,] = rollapply(fire_ssp126_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fire_ssp126_1_annual = array(NA, dim=c(dim(fire_ssp126_1)[1],dim(fire_ssp126_4)[2],dim(fire_ssp126_1)[3]/steps_per_Year))
for (i in seq(1,dim(fire_ssp126_1)[1])) {
     for (j in seq(1,dim(fire_ssp126_1)[2])) {
          fire_ssp126_1_annual[i,j,] = rollapply(fire_ssp126_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fire_ssp126_7_annual = array(NA, dim=c(dim(fire_ssp126_1)[1],dim(fire_ssp126_7)[2],dim(fire_ssp126_7)[3]/steps_per_Year))
for (i in seq(1,dim(fire_ssp126_7)[1])) {
     for (j in seq(1,dim(fire_ssp126_7)[2])) {
          fire_ssp126_7_annual[i,j,] = rollapply(fire_ssp126_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting fire to annual ssp245
fire_ssp245_4_annual = array(NA, dim=c(dim(fire_ssp245_4)[1],dim(fire_ssp245_4)[2],dim(fire_ssp245_4)[3]/steps_per_Year))
for (i in seq(1,dim(fire_ssp245_4)[1])) {
     for (j in seq(1,dim(fire_ssp245_4)[2])) {
          fire_ssp245_4_annual[i,j,] = rollapply(fire_ssp245_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fire_ssp245_1_annual = array(NA, dim=c(dim(fire_ssp245_1)[1],dim(fire_ssp245_4)[2],dim(fire_ssp245_1)[3]/steps_per_Year))
for (i in seq(1,dim(fire_ssp245_1)[1])) {
     for (j in seq(1,dim(fire_ssp245_1)[2])) {
          fire_ssp245_1_annual[i,j,] = rollapply(fire_ssp245_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fire_ssp245_7_annual = array(NA, dim=c(dim(fire_ssp245_7)[1],dim(fire_ssp245_7)[2],dim(fire_ssp245_7)[3]/steps_per_Year))
for (i in seq(1,dim(fire_ssp245_1)[1])) {
     for (j in seq(1,dim(fire_ssp245_1)[2])) {
          fire_ssp245_7_annual[i,j,] = rollapply(fire_ssp245_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting fire to annual ssp370
fire_ssp370_4_annual = array(NA, dim=c(dim(fire_ssp370_4)[1],dim(fire_ssp370_4)[2],dim(fire_ssp370_4)[3]/steps_per_Year))
for (i in seq(1,dim(fire_ssp370_4)[1])) {
     for (j in seq(1,dim(fire_ssp370_4)[2])) {
          fire_ssp370_4_annual[i,j,] = rollapply(fire_ssp370_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fire_ssp370_1_annual = array(NA, dim=c(dim(fire_ssp370_1)[1],dim(fire_ssp370_4)[2],dim(fire_ssp370_1)[3]/steps_per_Year))
for (i in seq(1,dim(fire_ssp370_1)[1])) {
     for (j in seq(1,dim(fire_ssp370_1)[2])) {
          fire_ssp370_1_annual[i,j,] = rollapply(fire_ssp370_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fire_ssp370_7_annual = array(NA, dim=c(dim(fire_ssp370_7)[1],dim(fire_ssp370_7)[2],dim(fire_ssp370_7)[3]/steps_per_Year))
for (i in seq(1,dim(fire_ssp370_7)[1])) {
     for (j in seq(1,dim(fire_ssp370_7)[2])) {
          fire_ssp370_7_annual[i,j,] = rollapply(fire_ssp370_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting fire to annual ssp585
fire_ssp585_4_annual = array(NA, dim=c(dim(fire_ssp585_4)[1],dim(fire_ssp585_4)[2],dim(fire_ssp585_4)[3]/steps_per_Year))
for (i in seq(1,dim(fire_ssp585_4)[1])) {
     for (j in seq(1,dim(fire_ssp585_4)[2])) {
          fire_ssp585_4_annual[i,j,] = rollapply(fire_ssp585_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fire_ssp585_1_annual = array(NA, dim=c(dim(fire_ssp585_1)[1],dim(fire_ssp585_4)[2],dim(fire_ssp585_1)[3]/steps_per_Year))
for (i in seq(1,dim(fire_ssp585_1)[1])) {
     for (j in seq(1,dim(fire_ssp585_1)[2])) {
          fire_ssp585_1_annual[i,j,] = rollapply(fire_ssp585_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fire_ssp585_7_annual = array(NA, dim=c(dim(fire_ssp585_7)[1],dim(fire_ssp585_7)[2],dim(fire_ssp585_7)[3]/steps_per_Year))
for (i in seq(1,dim(fire_ssp585_7)[1])) {
     for (j in seq(1,dim(fire_ssp585_7)[2])) {
          fire_ssp585_7_annual[i,j,] = rollapply(fire_ssp585_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# nbe
# Converting nbe to annual ssp126
nbe_ssp126_4_annual = array(NA, dim=c(dim(nbe_ssp126_4)[1],dim(nbe_ssp126_4)[2],dim(nbe_ssp126_4)[3]/steps_per_Year))
for (i in seq(1,dim(nbe_ssp126_4)[1])) {
     for (j in seq(1,dim(nbe_ssp126_4)[2])) {
          nbe_ssp126_4_annual[i,j,] = rollapply(nbe_ssp126_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbe_ssp126_1_annual = array(NA, dim=c(dim(nbe_ssp126_1)[1],dim(nbe_ssp126_4)[2],dim(nbe_ssp126_1)[3]/steps_per_Year))
for (i in seq(1,dim(nbe_ssp126_1)[1])) {
     for (j in seq(1,dim(nbe_ssp126_1)[2])) {
          nbe_ssp126_1_annual[i,j,] = rollapply(nbe_ssp126_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbe_ssp126_7_annual = array(NA, dim=c(dim(nbe_ssp126_1)[1],dim(nbe_ssp126_7)[2],dim(nbe_ssp126_7)[3]/steps_per_Year))
for (i in seq(1,dim(nbe_ssp126_7)[1])) {
     for (j in seq(1,dim(nbe_ssp126_7)[2])) {
          nbe_ssp126_7_annual[i,j,] = rollapply(nbe_ssp126_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting nbe to annual ssp245
nbe_ssp245_4_annual = array(NA, dim=c(dim(nbe_ssp245_4)[1],dim(nbe_ssp245_4)[2],dim(nbe_ssp245_4)[3]/steps_per_Year))
for (i in seq(1,dim(nbe_ssp245_4)[1])) {
     for (j in seq(1,dim(nbe_ssp245_4)[2])) {
          nbe_ssp245_4_annual[i,j,] = rollapply(nbe_ssp245_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbe_ssp245_1_annual = array(NA, dim=c(dim(nbe_ssp245_1)[1],dim(nbe_ssp245_4)[2],dim(nbe_ssp245_1)[3]/steps_per_Year))
for (i in seq(1,dim(nbe_ssp245_1)[1])) {
     for (j in seq(1,dim(nbe_ssp245_1)[2])) {
          nbe_ssp245_1_annual[i,j,] = rollapply(nbe_ssp245_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbe_ssp245_7_annual = array(NA, dim=c(dim(nbe_ssp245_7)[1],dim(nbe_ssp245_7)[2],dim(nbe_ssp245_7)[3]/steps_per_Year))
for (i in seq(1,dim(nbe_ssp245_1)[1])) {
     for (j in seq(1,dim(nbe_ssp245_1)[2])) {
          nbe_ssp245_7_annual[i,j,] = rollapply(nbe_ssp245_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting nbe to annual ssp370
nbe_ssp370_4_annual = array(NA, dim=c(dim(nbe_ssp370_4)[1],dim(nbe_ssp370_4)[2],dim(nbe_ssp370_4)[3]/steps_per_Year))
for (i in seq(1,dim(nbe_ssp370_4)[1])) {
     for (j in seq(1,dim(nbe_ssp370_4)[2])) {
          nbe_ssp370_4_annual[i,j,] = rollapply(nbe_ssp370_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbe_ssp370_1_annual = array(NA, dim=c(dim(nbe_ssp370_1)[1],dim(nbe_ssp370_4)[2],dim(nbe_ssp370_1)[3]/steps_per_Year))
for (i in seq(1,dim(nbe_ssp370_1)[1])) {
     for (j in seq(1,dim(nbe_ssp370_1)[2])) {
          nbe_ssp370_1_annual[i,j,] = rollapply(nbe_ssp370_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbe_ssp370_7_annual = array(NA, dim=c(dim(nbe_ssp370_7)[1],dim(nbe_ssp370_7)[2],dim(nbe_ssp370_7)[3]/steps_per_Year))
for (i in seq(1,dim(nbe_ssp370_7)[1])) {
     for (j in seq(1,dim(nbe_ssp370_7)[2])) {
          nbe_ssp370_7_annual[i,j,] = rollapply(nbe_ssp370_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting nbe to annual ssp585
nbe_ssp585_4_annual = array(NA, dim=c(dim(nbe_ssp585_4)[1],dim(nbe_ssp585_4)[2],dim(nbe_ssp585_4)[3]/steps_per_Year))
for (i in seq(1,dim(nbe_ssp585_4)[1])) {
     for (j in seq(1,dim(nbe_ssp585_4)[2])) {
          nbe_ssp585_4_annual[i,j,] = rollapply(nbe_ssp585_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbe_ssp585_1_annual = array(NA, dim=c(dim(nbe_ssp585_1)[1],dim(nbe_ssp585_4)[2],dim(nbe_ssp585_1)[3]/steps_per_Year))
for (i in seq(1,dim(nbe_ssp585_1)[1])) {
     for (j in seq(1,dim(nbe_ssp585_1)[2])) {
          nbe_ssp585_1_annual[i,j,] = rollapply(nbe_ssp585_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbe_ssp585_7_annual = array(NA, dim=c(dim(nbe_ssp585_7)[1],dim(nbe_ssp585_7)[2],dim(nbe_ssp585_7)[3]/steps_per_Year))
for (i in seq(1,dim(nbe_ssp585_7)[1])) {
     for (j in seq(1,dim(nbe_ssp585_7)[2])) {
          nbe_ssp585_7_annual[i,j,] = rollapply(nbe_ssp585_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# fLuc
# Converting fLuc to annual ssp126
fLuc_ssp126_4_annual = array(NA, dim=c(dim(fLuc_ssp126_4)[1],dim(fLuc_ssp126_4)[2],dim(fLuc_ssp126_4)[3]/steps_per_Year))
for (i in seq(1,dim(fLuc_ssp126_4)[1])) {
     for (j in seq(1,dim(fLuc_ssp126_4)[2])) {
          fLuc_ssp126_4_annual[i,j,] = rollapply(fLuc_ssp126_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fLuc_ssp126_1_annual = array(NA, dim=c(dim(fLuc_ssp126_1)[1],dim(fLuc_ssp126_4)[2],dim(fLuc_ssp126_1)[3]/steps_per_Year))
for (i in seq(1,dim(fLuc_ssp126_1)[1])) {
     for (j in seq(1,dim(fLuc_ssp126_1)[2])) {
          fLuc_ssp126_1_annual[i,j,] = rollapply(fLuc_ssp126_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fLuc_ssp126_7_annual = array(NA, dim=c(dim(fLuc_ssp126_1)[1],dim(fLuc_ssp126_7)[2],dim(fLuc_ssp126_7)[3]/steps_per_Year))
for (i in seq(1,dim(fLuc_ssp126_7)[1])) {
     for (j in seq(1,dim(fLuc_ssp126_7)[2])) {
          fLuc_ssp126_7_annual[i,j,] = rollapply(fLuc_ssp126_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting fLuc to annual ssp245
fLuc_ssp245_4_annual = array(NA, dim=c(dim(fLuc_ssp245_4)[1],dim(fLuc_ssp245_4)[2],dim(fLuc_ssp245_4)[3]/steps_per_Year))
for (i in seq(1,dim(fLuc_ssp245_4)[1])) {
     for (j in seq(1,dim(fLuc_ssp245_4)[2])) {
          fLuc_ssp245_4_annual[i,j,] = rollapply(fLuc_ssp245_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fLuc_ssp245_1_annual = array(NA, dim=c(dim(fLuc_ssp245_1)[1],dim(fLuc_ssp245_4)[2],dim(fLuc_ssp245_1)[3]/steps_per_Year))
for (i in seq(1,dim(fLuc_ssp245_1)[1])) {
     for (j in seq(1,dim(fLuc_ssp245_1)[2])) {
          fLuc_ssp245_1_annual[i,j,] = rollapply(fLuc_ssp245_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fLuc_ssp245_7_annual = array(NA, dim=c(dim(fLuc_ssp245_7)[1],dim(fLuc_ssp245_7)[2],dim(fLuc_ssp245_7)[3]/steps_per_Year))
for (i in seq(1,dim(fLuc_ssp245_1)[1])) {
     for (j in seq(1,dim(fLuc_ssp245_1)[2])) {
          fLuc_ssp245_7_annual[i,j,] = rollapply(fLuc_ssp245_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting fLuc to annual ssp370
fLuc_ssp370_4_annual = array(NA, dim=c(dim(fLuc_ssp370_4)[1],dim(fLuc_ssp370_4)[2],dim(fLuc_ssp370_4)[3]/steps_per_Year))
for (i in seq(1,dim(fLuc_ssp370_4)[1])) {
     for (j in seq(1,dim(fLuc_ssp370_4)[2])) {
          fLuc_ssp370_4_annual[i,j,] = rollapply(fLuc_ssp370_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fLuc_ssp370_1_annual = array(NA, dim=c(dim(fLuc_ssp370_1)[1],dim(fLuc_ssp370_4)[2],dim(fLuc_ssp370_1)[3]/steps_per_Year))
for (i in seq(1,dim(fLuc_ssp370_1)[1])) {
     for (j in seq(1,dim(fLuc_ssp370_1)[2])) {
          fLuc_ssp370_1_annual[i,j,] = rollapply(fLuc_ssp370_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fLuc_ssp370_7_annual = array(NA, dim=c(dim(fLuc_ssp370_7)[1],dim(fLuc_ssp370_7)[2],dim(fLuc_ssp370_7)[3]/steps_per_Year))
for (i in seq(1,dim(fLuc_ssp370_7)[1])) {
     for (j in seq(1,dim(fLuc_ssp370_7)[2])) {
          fLuc_ssp370_7_annual[i,j,] = rollapply(fLuc_ssp370_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting fLuc to annual ssp585
fLuc_ssp585_4_annual = array(NA, dim=c(dim(fLuc_ssp585_4)[1],dim(fLuc_ssp585_4)[2],dim(fLuc_ssp585_4)[3]/steps_per_Year))
for (i in seq(1,dim(fLuc_ssp585_4)[1])) {
     for (j in seq(1,dim(fLuc_ssp585_4)[2])) {
          fLuc_ssp585_4_annual[i,j,] = rollapply(fLuc_ssp585_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fLuc_ssp585_1_annual = array(NA, dim=c(dim(fLuc_ssp585_1)[1],dim(fLuc_ssp585_4)[2],dim(fLuc_ssp585_1)[3]/steps_per_Year))
for (i in seq(1,dim(fLuc_ssp585_1)[1])) {
     for (j in seq(1,dim(fLuc_ssp585_1)[2])) {
          fLuc_ssp585_1_annual[i,j,] = rollapply(fLuc_ssp585_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

fLuc_ssp585_7_annual = array(NA, dim=c(dim(fLuc_ssp585_7)[1],dim(fLuc_ssp585_7)[2],dim(fLuc_ssp585_7)[3]/steps_per_Year))
for (i in seq(1,dim(fLuc_ssp585_7)[1])) {
     for (j in seq(1,dim(fLuc_ssp585_7)[2])) {
          fLuc_ssp585_7_annual[i,j,] = rollapply(fLuc_ssp585_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# nbp
# Converting nbp to annual ssp126
nbp_ssp126_4_annual = array(NA, dim=c(dim(nbp_ssp126_4)[1],dim(nbp_ssp126_4)[2],dim(nbp_ssp126_4)[3]/steps_per_Year))
for (i in seq(1,dim(nbp_ssp126_4)[1])) {
     for (j in seq(1,dim(nbp_ssp126_4)[2])) {
          nbp_ssp126_4_annual[i,j,] = rollapply(nbp_ssp126_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbp_ssp126_1_annual = array(NA, dim=c(dim(nbp_ssp126_1)[1],dim(nbp_ssp126_4)[2],dim(nbp_ssp126_1)[3]/steps_per_Year))
for (i in seq(1,dim(nbp_ssp126_1)[1])) {
     for (j in seq(1,dim(nbp_ssp126_1)[2])) {
          nbp_ssp126_1_annual[i,j,] = rollapply(nbp_ssp126_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbp_ssp126_7_annual = array(NA, dim=c(dim(nbp_ssp126_1)[1],dim(nbp_ssp126_7)[2],dim(nbp_ssp126_7)[3]/steps_per_Year))
for (i in seq(1,dim(nbp_ssp126_7)[1])) {
     for (j in seq(1,dim(nbp_ssp126_7)[2])) {
          nbp_ssp126_7_annual[i,j,] = rollapply(nbp_ssp126_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting nbp to annual ssp245
nbp_ssp245_4_annual = array(NA, dim=c(dim(nbp_ssp245_4)[1],dim(nbp_ssp245_4)[2],dim(nbp_ssp245_4)[3]/steps_per_Year))
for (i in seq(1,dim(nbp_ssp245_4)[1])) {
     for (j in seq(1,dim(nbp_ssp245_4)[2])) {
          nbp_ssp245_4_annual[i,j,] = rollapply(nbp_ssp245_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbp_ssp245_1_annual = array(NA, dim=c(dim(nbp_ssp245_1)[1],dim(nbp_ssp245_4)[2],dim(nbp_ssp245_1)[3]/steps_per_Year))
for (i in seq(1,dim(nbp_ssp245_1)[1])) {
     for (j in seq(1,dim(nbp_ssp245_1)[2])) {
          nbp_ssp245_1_annual[i,j,] = rollapply(nbp_ssp245_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbp_ssp245_7_annual = array(NA, dim=c(dim(nbp_ssp245_7)[1],dim(nbp_ssp245_7)[2],dim(nbp_ssp245_7)[3]/steps_per_Year))
for (i in seq(1,dim(nbp_ssp245_1)[1])) {
     for (j in seq(1,dim(nbp_ssp245_1)[2])) {
          nbp_ssp245_7_annual[i,j,] = rollapply(nbp_ssp245_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting nbp to annual ssp370
nbp_ssp370_4_annual = array(NA, dim=c(dim(nbp_ssp370_4)[1],dim(nbp_ssp370_4)[2],dim(nbp_ssp370_4)[3]/steps_per_Year))
for (i in seq(1,dim(nbp_ssp370_4)[1])) {
     for (j in seq(1,dim(nbp_ssp370_4)[2])) {
          nbp_ssp370_4_annual[i,j,] = rollapply(nbp_ssp370_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbp_ssp370_1_annual = array(NA, dim=c(dim(nbp_ssp370_1)[1],dim(nbp_ssp370_4)[2],dim(nbp_ssp370_1)[3]/steps_per_Year))
for (i in seq(1,dim(nbp_ssp370_1)[1])) {
     for (j in seq(1,dim(nbp_ssp370_1)[2])) {
          nbp_ssp370_1_annual[i,j,] = rollapply(nbp_ssp370_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbp_ssp370_7_annual = array(NA, dim=c(dim(nbp_ssp370_7)[1],dim(nbp_ssp370_7)[2],dim(nbp_ssp370_7)[3]/steps_per_Year))
for (i in seq(1,dim(nbp_ssp370_7)[1])) {
     for (j in seq(1,dim(nbp_ssp370_7)[2])) {
          nbp_ssp370_7_annual[i,j,] = rollapply(nbp_ssp370_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

# Converting nbp to annual ssp585
nbp_ssp585_4_annual = array(NA, dim=c(dim(nbp_ssp585_4)[1],dim(nbp_ssp585_4)[2],dim(nbp_ssp585_4)[3]/steps_per_Year))
for (i in seq(1,dim(nbp_ssp585_4)[1])) {
     for (j in seq(1,dim(nbp_ssp585_4)[2])) {
          nbp_ssp585_4_annual[i,j,] = rollapply(nbp_ssp585_4[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbp_ssp585_1_annual = array(NA, dim=c(dim(nbp_ssp585_1)[1],dim(nbp_ssp585_4)[2],dim(nbp_ssp585_1)[3]/steps_per_Year))
for (i in seq(1,dim(nbp_ssp585_1)[1])) {
     for (j in seq(1,dim(nbp_ssp585_1)[2])) {
          nbp_ssp585_1_annual[i,j,] = rollapply(nbp_ssp585_1[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

nbp_ssp585_7_annual = array(NA, dim=c(dim(nbp_ssp585_7)[1],dim(nbp_ssp585_7)[2],dim(nbp_ssp585_7)[3]/steps_per_Year))
for (i in seq(1,dim(nbp_ssp585_7)[1])) {
     for (j in seq(1,dim(nbp_ssp585_7)[2])) {
          nbp_ssp585_7_annual[i,j,] = rollapply(nbp_ssp585_7[i,j,], by = steps_per_Year, width = steps_per_Year, FUN=mean, na.rm=TRUE) } } 

#Creating X axis to plot against
x <- c(2001:2100)
x1 <- c(2090:2100)



par(mfrow=c(3,3))
# Plotting gpp
plot(x,apply(gpp_ssp126_4_annual,3,mean, na.rm=TRUE), type="l", ylim=c(2,14), xlab="Year",cex.lab=1.2, main="GPP", ylab="",col="goldenrod3", lwd=4)
lines(x1,apply(gpp_ssp126_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x1,apply(gpp_ssp126_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x,apply(gpp_ssp245_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col="darkgreen")
lines(x1,apply(gpp_ssp245_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x1,apply(gpp_ssp245_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x,apply(gpp_ssp370_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="blue")
lines(x1,apply(gpp_ssp370_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x1,apply(gpp_ssp370_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x,apply(gpp_ssp585_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="firebrick")
lines(x1,apply(gpp_ssp585_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")
lines(x1,apply(gpp_ssp585_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")

# Plotting ra
plot(x,apply(ra_ssp126_4_annual,3,mean, na.rm=TRUE), type="l", ylim=c(0.5,8), xlab="Year",cex.lab=1.2, main="Ra",ylab="", col="goldenrod3", lwd=4)
lines(x1,apply(ra_ssp126_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x1,apply(ra_ssp126_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x,apply(ra_ssp245_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col="darkgreen")
lines(x1,apply(ra_ssp245_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x1,apply(ra_ssp245_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x,apply(ra_ssp370_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="blue")
lines(x1,apply(ra_ssp370_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x1,apply(ra_ssp370_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x,apply(ra_ssp585_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="firebrick")
lines(x1,apply(ra_ssp585_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")
lines(x1,apply(ra_ssp585_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")

# Plotting npp
plot(x,apply(npp_ssp126_4_annual,3,mean, na.rm=TRUE), type="l", ylim=c(0.5,9), xlab="Year",cex.lab=1.2, main="NPP", ylab="",col="goldenrod3", lwd=4)
lines(x1,apply(npp_ssp126_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x1,apply(npp_ssp126_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x,apply(npp_ssp245_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col="darkgreen")
lines(x1,apply(npp_ssp245_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x1,apply(npp_ssp245_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x,apply(npp_ssp370_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="blue")
lines(x1,apply(npp_ssp370_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x1,apply(npp_ssp370_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x,apply(npp_ssp585_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="firebrick")
lines(x1,apply(npp_ssp585_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")
lines(x1,apply(npp_ssp585_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")

# Plotting rh
plot(x,apply(rh_ssp126_4_annual,3,mean, na.rm=TRUE), type="l", ylim=c(0.5,8), xlab="Year",cex.lab=1.2, main="Rh",ylab="", col="goldenrod3", lwd=4)
lines(x1,apply(rh_ssp126_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x1,apply(rh_ssp126_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x,apply(rh_ssp245_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col="darkgreen")
lines(x1,apply(rh_ssp245_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x1,apply(rh_ssp245_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x,apply(rh_ssp370_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="blue")
lines(x1,apply(rh_ssp370_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x1,apply(rh_ssp370_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x,apply(rh_ssp585_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="firebrick")
lines(x1,apply(rh_ssp585_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")
lines(x1,apply(rh_ssp585_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")

# Plotting nee
plot(x,apply(nee_ssp126_4_annual,3,mean, na.rm=TRUE), type="l", ylim=c(-3.5,1), xlab="Year",cex.lab=1.2, main="NEE", ylab="",col="goldenrod3", lwd=4)
lines(x1,apply(nee_ssp126_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x1,apply(nee_ssp126_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x,apply(nee_ssp245_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col="darkgreen")
lines(x1,apply(nee_ssp245_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x1,apply(nee_ssp245_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x,apply(nee_ssp370_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="blue")
lines(x1,apply(nee_ssp370_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x1,apply(nee_ssp370_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x,apply(nee_ssp585_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="firebrick")
lines(x1,apply(nee_ssp585_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")
lines(x1,apply(nee_ssp585_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")

# Plotting fire
plot(x,apply(fire_ssp126_4_annual,3,mean, na.rm=TRUE), type="l", ylim=c(0,0.1), xlab="Year",cex.lab=1.2, main="Fire",ylab="", col="goldenrod3", lwd=2)
lines(x1,apply(fire_ssp126_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=1, col="goldenrod3")
lines(x1,apply(fire_ssp126_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=1, col="goldenrod3")
lines(x,apply(fire_ssp245_4_annual,3,mean, na.rm=TRUE), type="l", lwd=2, col="darkgreen")
lines(x1,apply(fire_ssp245_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=1, col="darkgreen")
lines(x1,apply(fire_ssp245_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=1, col="darkgreen")
lines(x,apply(fire_ssp370_4_annual,3,mean, na.rm=TRUE), type="l", lwd=2, col ="blue")
lines(x1,apply(fire_ssp370_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=1, col ="blue")
lines(x1,apply(fire_ssp370_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=1, col ="blue")
lines(x,apply(fire_ssp585_4_annual,3,mean, na.rm=TRUE), type="l", lwd=2, col ="firebrick")
lines(x1,apply(fire_ssp585_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=1, col ="firebrick")
lines(x1,apply(fire_ssp585_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=1, col ="firebrick")

# Plotting nbe
plot(x,apply(nbe_ssp126_4_annual,3,mean, na.rm=TRUE), type="l", ylim=c(-3.5,1), xlab="Year",cex.lab=1.2, main="NBE", ylab="",col="goldenrod3", lwd=4)
lines(x1,apply(nbe_ssp126_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x1,apply(nbe_ssp126_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x,apply(nbe_ssp245_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col="darkgreen")
lines(x1,apply(nbe_ssp245_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x1,apply(nbe_ssp245_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x,apply(nbe_ssp370_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="blue")
lines(x1,apply(nbe_ssp370_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x1,apply(nbe_ssp370_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x,apply(nbe_ssp585_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="firebrick")
lines(x1,apply(nbe_ssp585_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")
lines(x1,apply(nbe_ssp585_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")

# Plotting fLuc
plot(x,apply(fLuc_ssp126_4_annual,3,mean, na.rm=TRUE), type="l", ylim=c(0,0.6), xlab="Year",cex.lab=1.2, main="FLuc",ylab="", col="goldenrod3", lwd=4)
lines(x1,apply(fLuc_ssp126_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x1,apply(fLuc_ssp126_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x,apply(fLuc_ssp245_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col="darkgreen")
lines(x1,apply(fLuc_ssp245_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x1,apply(fLuc_ssp245_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x,apply(fLuc_ssp370_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="blue")
lines(x1,apply(fLuc_ssp370_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x1,apply(fLuc_ssp370_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x,apply(fLuc_ssp585_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="firebrick")
lines(x1,apply(fLuc_ssp585_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")
lines(x1,apply(fLuc_ssp585_7_annual[,,90:100],3,mean, na.rm=TRUE), xlim=c(90,100),type="l", lty=3,lwd=3, col ="firebrick")

# Plotting nbp
plot(x, apply(nbp_ssp126_4_annual,3,mean, na.rm=TRUE), type="l", ylim=c(-1,3.5), xlab="Year",cex.lab=1.2, main="NBP",ylab="", col="goldenrod3", lwd=4)
lines(x1, apply(nbp_ssp126_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x1,apply(nbp_ssp126_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="goldenrod3")
lines(x,apply(nbp_ssp245_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col="darkgreen")
lines(x1,apply(nbp_ssp245_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x1,apply(nbp_ssp245_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col="darkgreen")
lines(x,apply(nbp_ssp370_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="blue")
lines(x1,apply(nbp_ssp370_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x1,apply(nbp_ssp370_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="blue")
lines(x,apply(nbp_ssp585_4_annual,3,mean, na.rm=TRUE), type="l", lwd=4, col ="firebrick")
lines(x1,apply(nbp_ssp585_1_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")
lines(x1,apply(nbp_ssp585_7_annual[,,90:100],3,mean, na.rm=TRUE), type="l", lty=3,lwd=3, col ="firebrick")


