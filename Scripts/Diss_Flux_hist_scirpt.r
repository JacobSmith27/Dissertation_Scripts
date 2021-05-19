library(ncdf4)
library(fields)
library(dplyr)
library(tidyr)
library(lubridate)
library(raster)
library(zoo)
library(ggplot2)
library(gridExtra)

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

# GPP

# Cutting sections
# Avg the last 5 years of observed data (2014 -2019) vs the las 5 year of the the porjections
# Gpp median
gpp_ssp126_4_avgfirst <- apply(gpp_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
gpp_ssp126_4_avglast <- apply(gpp_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
gpp_ssp245_4_avglast <- apply(gpp_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
gpp_ssp370_4_avglast <- apply(gpp_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
gpp_ssp585_4_avglast <- apply(gpp_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)

# creating a single dataframe
gpp_ssp126_4_avgfirst.df <- as.data.frame(as.table(gpp_ssp126_4_avgfirst,na.rm=TRUE))
gpp_ssp126_4_avglast.df <- as.data.frame(as.table(gpp_ssp126_4_avglast,na.rm=TRUE))
gpp_ssp245_4_avglast.df <- as.data.frame(as.table(gpp_ssp245_4_avglast,na.rm=TRUE))
gpp_ssp370_4_avglast.df <- as.data.frame(as.table(gpp_ssp370_4_avglast,na.rm=TRUE))
gpp_ssp585_4_avglast.df <- as.data.frame(as.table(gpp_ssp585_4_avglast,na.rm=TRUE))

gpp.df <- merge(x=gpp_ssp126_4_avgfirst.df, y=gpp_ssp126_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(current=Freq.x, ssp126=Freq.y)
gpp.df= merge(x=gpp.df, y=gpp_ssp245_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp245=Freq)
gpp.df= merge(x=gpp.df, y=gpp_ssp370_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp370=Freq)
gpp.df= merge(x=gpp.df, y=gpp_ssp585_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp585=Freq)

gpp.df <- gather(gpp.df, scenario, gpp, current:ssp585, factor_key=TRUE, na.rm=TRUE)

# Ra

# Cutting sections
# Avg the last 5 years of observed data (2014 -2019) vs the las 5 year of the the porjections
# ra median
ra_ssp126_4_avgfirst <- apply(ra_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
ra_ssp126_4_avglast <- apply(ra_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
ra_ssp245_4_avglast <- apply(ra_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
ra_ssp370_4_avglast <- apply(ra_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
ra_ssp585_4_avglast <- apply(ra_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)

# creating a single dataframe
ra_ssp126_4_avgfirst.df <- as.data.frame(as.table(ra_ssp126_4_avgfirst,na.rm=TRUE))
ra_ssp126_4_avglast.df <- as.data.frame(as.table(ra_ssp126_4_avglast,na.rm=TRUE))
ra_ssp245_4_avglast.df <- as.data.frame(as.table(ra_ssp245_4_avglast,na.rm=TRUE))
ra_ssp370_4_avglast.df <- as.data.frame(as.table(ra_ssp370_4_avglast,na.rm=TRUE))
ra_ssp585_4_avglast.df <- as.data.frame(as.table(ra_ssp585_4_avglast,na.rm=TRUE))

ra.df <- merge(x=ra_ssp126_4_avgfirst.df, y=ra_ssp126_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(current=Freq.x, ssp126=Freq.y)
ra.df= merge(x=ra.df, y=ra_ssp245_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp245=Freq)
ra.df= merge(x=ra.df, y=ra_ssp370_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp370=Freq)
ra.df= merge(x=ra.df, y=ra_ssp585_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp585=Freq)

ra.df <- gather(ra.df, scenario, ra, current:ssp585, factor_key=TRUE, na.rm=TRUE)
 
# npp

# Cutting sections
# Avg the last 5 years of observed data (2014 -2019) vs the las 5 year of the the porjections
# npp median
npp_ssp126_4_avgfirst <- apply(npp_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
npp_ssp126_4_avglast <- apply(npp_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
npp_ssp245_4_avglast <- apply(npp_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
npp_ssp370_4_avglast <- apply(npp_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
npp_ssp585_4_avglast <- apply(npp_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)

# creating a single dataframe
npp_ssp126_4_avgfirst.df <- as.data.frame(as.table(npp_ssp126_4_avgfirst,na.rm=TRUE))
npp_ssp126_4_avglast.df <- as.data.frame(as.table(npp_ssp126_4_avglast,na.rm=TRUE))
npp_ssp245_4_avglast.df <- as.data.frame(as.table(npp_ssp245_4_avglast,na.rm=TRUE))
npp_ssp370_4_avglast.df <- as.data.frame(as.table(npp_ssp370_4_avglast,na.rm=TRUE))
npp_ssp585_4_avglast.df <- as.data.frame(as.table(npp_ssp585_4_avglast,na.rm=TRUE))

npp.df <- merge(x=npp_ssp126_4_avgfirst.df, y=npp_ssp126_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(current=Freq.x, ssp126=Freq.y)
npp.df= merge(x=npp.df, y=npp_ssp245_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp245=Freq)
npp.df= merge(x=npp.df, y=npp_ssp370_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp370=Freq)
npp.df= merge(x=npp.df, y=npp_ssp585_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp585=Freq)

npp.df <- gather(npp.df, scenario, npp, current:ssp585, factor_key=TRUE, na.rm=TRUE)

# Plotting npp hist
# rh

# Cutting sections
# Avg the last 5 years of observed data (2014 -2019) vs the las 5 year of the the porjections
# rh median
rh_ssp126_4_avgfirst <- apply(rh_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
rh_ssp126_4_avglast <- apply(rh_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
rh_ssp245_4_avglast <- apply(rh_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
rh_ssp370_4_avglast <- apply(rh_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
rh_ssp585_4_avglast <- apply(rh_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)

# creating a single dataframe
rh_ssp126_4_avgfirst.df <- as.data.frame(as.table(rh_ssp126_4_avgfirst,na.rm=TRUE))
rh_ssp126_4_avglast.df <- as.data.frame(as.table(rh_ssp126_4_avglast,na.rm=TRUE))
rh_ssp245_4_avglast.df <- as.data.frame(as.table(rh_ssp245_4_avglast,na.rm=TRUE))
rh_ssp370_4_avglast.df <- as.data.frame(as.table(rh_ssp370_4_avglast,na.rm=TRUE))
rh_ssp585_4_avglast.df <- as.data.frame(as.table(rh_ssp585_4_avglast,na.rm=TRUE))

rh.df <- merge(x=rh_ssp126_4_avgfirst.df, y=rh_ssp126_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(current=Freq.x, ssp126=Freq.y)
rh.df= merge(x=rh.df, y=rh_ssp245_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp245=Freq)
rh.df= merge(x=rh.df, y=rh_ssp370_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp370=Freq)
rh.df= merge(x=rh.df, y=rh_ssp585_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp585=Freq)

rh.df <- gather(rh.df, scenario, rh, current:ssp585, factor_key=TRUE, na.rm=TRUE)

# nee

# Cutting sections
# Avg the last 5 years of observed data (2014 -2019) vs the las 5 year of the the porjections
# nee median
nee_ssp126_4_avgfirst <- apply(nee_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
nee_ssp126_4_avglast <- apply(nee_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nee_ssp245_4_avglast <- apply(nee_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nee_ssp370_4_avglast <- apply(nee_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nee_ssp585_4_avglast <- apply(nee_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)

# creating a single dataframe
nee_ssp126_4_avgfirst.df <- as.data.frame(as.table(nee_ssp126_4_avgfirst,na.rm=TRUE))
nee_ssp126_4_avglast.df <- as.data.frame(as.table(nee_ssp126_4_avglast,na.rm=TRUE))
nee_ssp245_4_avglast.df <- as.data.frame(as.table(nee_ssp245_4_avglast,na.rm=TRUE))
nee_ssp370_4_avglast.df <- as.data.frame(as.table(nee_ssp370_4_avglast,na.rm=TRUE))
nee_ssp585_4_avglast.df <- as.data.frame(as.table(nee_ssp585_4_avglast,na.rm=TRUE))

nee.df <- merge(x=nee_ssp126_4_avgfirst.df, y=nee_ssp126_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(current=Freq.x, ssp126=Freq.y)
nee.df= merge(x=nee.df, y=nee_ssp245_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp245=Freq)
nee.df= merge(x=nee.df, y=nee_ssp370_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp370=Freq)
nee.df= merge(x=nee.df, y=nee_ssp585_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp585=Freq)

nee.df <- gather(nee.df, scenario, nee, current:ssp585, factor_key=TRUE, na.rm=TRUE)

# fire

# Cutting sections
# Avg the last 5 years of observed data (2014 -2019) vs the las 5 year of the the porjections
# fire median
fire_ssp126_4_avgfirst <- apply(fire_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
fire_ssp126_4_avglast <- apply(fire_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
fire_ssp245_4_avglast <- apply(fire_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
fire_ssp370_4_avglast <- apply(fire_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
fire_ssp585_4_avglast <- apply(fire_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)

# creating a single dataframe
fire_ssp126_4_avgfirst.df <- as.data.frame(as.table(fire_ssp126_4_avgfirst,na.rm=TRUE))
fire_ssp126_4_avglast.df <- as.data.frame(as.table(fire_ssp126_4_avglast,na.rm=TRUE))
fire_ssp245_4_avglast.df <- as.data.frame(as.table(fire_ssp245_4_avglast,na.rm=TRUE))
fire_ssp370_4_avglast.df <- as.data.frame(as.table(fire_ssp370_4_avglast,na.rm=TRUE))
fire_ssp585_4_avglast.df <- as.data.frame(as.table(fire_ssp585_4_avglast,na.rm=TRUE))

fire.df <- merge(x=fire_ssp126_4_avgfirst.df, y=fire_ssp126_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(current=Freq.x, ssp126=Freq.y)
fire.df= merge(x=fire.df, y=fire_ssp245_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp245=Freq)
fire.df= merge(x=fire.df, y=fire_ssp370_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp370=Freq)
fire.df= merge(x=fire.df, y=fire_ssp585_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp585=Freq)

fire.df <- gather(fire.df, scenario, fire, current:ssp585, factor_key=TRUE, na.rm=TRUE)


# nbe

# Cutting sections
# Avg the last 5 years of observed data (2014 -2019) vs the las 5 year of the the porjections
# nbe median
nbe_ssp126_4_avgfirst <- apply(nbe_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
nbe_ssp126_4_avglast <- apply(nbe_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nbe_ssp245_4_avglast <- apply(nbe_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nbe_ssp370_4_avglast <- apply(nbe_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nbe_ssp585_4_avglast <- apply(nbe_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)

# creating a single dataframe
nbe_ssp126_4_avgfirst.df <- as.data.frame(as.table(nbe_ssp126_4_avgfirst,na.rm=TRUE))
nbe_ssp126_4_avglast.df <- as.data.frame(as.table(nbe_ssp126_4_avglast,na.rm=TRUE))
nbe_ssp245_4_avglast.df <- as.data.frame(as.table(nbe_ssp245_4_avglast,na.rm=TRUE))
nbe_ssp370_4_avglast.df <- as.data.frame(as.table(nbe_ssp370_4_avglast,na.rm=TRUE))
nbe_ssp585_4_avglast.df <- as.data.frame(as.table(nbe_ssp585_4_avglast,na.rm=TRUE))

nbe.df <- merge(x=nbe_ssp126_4_avgfirst.df, y=nbe_ssp126_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(current=Freq.x, ssp126=Freq.y)
nbe.df= merge(x=nbe.df, y=nbe_ssp245_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp245=Freq)
nbe.df= merge(x=nbe.df, y=nbe_ssp370_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp370=Freq)
nbe.df= merge(x=nbe.df, y=nbe_ssp585_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp585=Freq)

nbe.df <- gather(nbe.df, scenario, nbe, current:ssp585, factor_key=TRUE, na.rm=TRUE)

# fLuc

# Cutting sections
# Avg the last 5 years of observed data (2014 -2019) vs the las 5 year of the the porjections
# fLuc median
fLuc_ssp126_4_avgfirst <- apply(fLuc_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
fLuc_ssp126_4_avglast <- apply(fLuc_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
fLuc_ssp245_4_avglast <- apply(fLuc_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
fLuc_ssp370_4_avglast <- apply(fLuc_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
fLuc_ssp585_4_avglast <- apply(fLuc_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)

# creating a single dataframe
fLuc_ssp126_4_avgfirst.df <- as.data.frame(as.table(fLuc_ssp126_4_avgfirst,na.rm=TRUE))
fLuc_ssp126_4_avglast.df <- as.data.frame(as.table(fLuc_ssp126_4_avglast,na.rm=TRUE))
fLuc_ssp245_4_avglast.df <- as.data.frame(as.table(fLuc_ssp245_4_avglast,na.rm=TRUE))
fLuc_ssp370_4_avglast.df <- as.data.frame(as.table(fLuc_ssp370_4_avglast,na.rm=TRUE))
fLuc_ssp585_4_avglast.df <- as.data.frame(as.table(fLuc_ssp585_4_avglast,na.rm=TRUE))

fLuc.df <- merge(x=fLuc_ssp126_4_avgfirst.df, y=fLuc_ssp126_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(current=Freq.x, ssp126=Freq.y)
fLuc.df= merge(x=fLuc.df, y=fLuc_ssp245_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp245=Freq)
fLuc.df= merge(x=fLuc.df, y=fLuc_ssp370_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp370=Freq)
fLuc.df= merge(x=fLuc.df, y=fLuc_ssp585_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp585=Freq)

fLuc.df <- gather(fLuc.df, scenario, fLuc, current:ssp585, factor_key=TRUE, na.rm=TRUE)

# nbp

# Cutting sections
# Avg the last 5 years of observed data (2014 -2019) vs the las 5 year of the the porjections
# nbp median
nbp_ssp126_4_avgfirst <- apply(nbp_ssp126[,,4,168:228], c(1,2), mean, na.rm=TRUE)
nbp_ssp126_4_avglast <- apply(nbp_ssp126[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nbp_ssp245_4_avglast <- apply(nbp_ssp245[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nbp_ssp370_4_avglast <- apply(nbp_ssp370[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)
nbp_ssp585_4_avglast <- apply(nbp_ssp585[,,4,1140:1200], c(1,2), mean, na.rm=TRUE)

# creating a single dataframe
nbp_ssp126_4_avgfirst.df <- as.data.frame(as.table(nbp_ssp126_4_avgfirst,na.rm=TRUE))
nbp_ssp126_4_avglast.df <- as.data.frame(as.table(nbp_ssp126_4_avglast,na.rm=TRUE))
nbp_ssp245_4_avglast.df <- as.data.frame(as.table(nbp_ssp245_4_avglast,na.rm=TRUE))
nbp_ssp370_4_avglast.df <- as.data.frame(as.table(nbp_ssp370_4_avglast,na.rm=TRUE))
nbp_ssp585_4_avglast.df <- as.data.frame(as.table(nbp_ssp585_4_avglast,na.rm=TRUE))

nbp.df <- merge(x=nbp_ssp126_4_avgfirst.df, y=nbp_ssp126_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(current=Freq.x, ssp126=Freq.y)
nbp.df= merge(x=nbp.df, y=nbp_ssp245_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp245=Freq)
nbp.df= merge(x=nbp.df, y=nbp_ssp370_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp370=Freq)
nbp.df= merge(x=nbp.df, y=nbp_ssp585_4_avglast.df, by=c('Var1','Var2'))%>%
	rename(ssp585=Freq)

nbp.df <- gather(nbp.df, scenario, nbp, current:ssp585, factor_key=TRUE, na.rm=TRUE)


# Plotting all density plots

# Plotting gpp denisty
gpp_den <- ggplot(gpp.df) +                
  geom_density(aes(x = gpp, fill=scenario), alpha=0.2)+     
    theme_bw() +                                                 
  ylab("") +                                            
  xlab("\nGPP")  +                      
  theme(axis.text = element_text(size = 12),                   
        axis.title.x = element_text(size = 14, face = "plain"), 
        panel.grid = element_blank(),                           
        plot.margin = unit(c(1,1,1,1), units = , "cm"))

# Plotting ra denisty
ra_den <- ggplot(ra.df) +                
  geom_density(aes(x = ra, fill=scenario), alpha=0.2)+
    theme_bw() +                                                 
  ylab("") +                                            
  xlab("\nRa")  +                      
  theme(axis.text = element_text(size = 12),                   
        axis.title.x = element_text(size = 14, face = "plain"), 
        panel.grid = element_blank(),                           
        plot.margin = unit(c(1,1,1,1), units = , "cm"))

# Plotting npp denisty
npp_den <- ggplot(npp.df) +                
  geom_density(aes(x = npp, fill=scenario), alpha=0.2)+     
    theme_bw() +                                                 
  ylab("") +                                            
  xlab("\nNPP")  +                      
  theme(axis.text = element_text(size = 12),                   
        axis.title.x = element_text(size = 14, face = "plain"), 
        panel.grid = element_blank(),                           
        plot.margin = unit(c(1,1,1,1), units = , "cm"))

# Plotting rh denisty
rh_den <- ggplot(rh.df) +                
  geom_density(aes(x = rh, fill=scenario), alpha=0.2)+  
    theme_bw() +                                                 
  ylab("") +                                            
  xlab("\nRh")  +                      
  theme(axis.text = element_text(size = 12),                   
        axis.title.x = element_text(size = 14, face = "plain"), 
        panel.grid = element_blank(),                           
        plot.margin = unit(c(1,1,1,1), units = , "cm"))

# Plotting nee denisty
nee_den <- ggplot(nee.df) +                
  geom_density(aes(x = nee, fill=scenario), alpha=0.2)+  
    theme_bw() +                                                 
  ylab("") +                                            
  xlab("\nNEE")  +                      
  theme(axis.text = element_text(size = 12),                   
        axis.title.x = element_text(size = 14, face = "plain"), 
        panel.grid = element_blank(),                           
        plot.margin = unit(c(1,1,1,1), units = , "cm"))

# Plotting fire denisty
fire_den <- ggplot(fire.df) +                
  geom_density(aes(x = fire, fill=scenario), alpha=0.2)+
  scale_x_log10()+    
    theme_bw() +                                                 
  ylab("") +                                            
  xlab("\nFire")  +                      
  theme(axis.text = element_text(size = 12),                   
        axis.title.x = element_text(size = 14, face = "plain"), 
        panel.grid = element_blank(),                           
        plot.margin = unit(c(1,1,1,1), units = , "cm"))

# Plotting nbe denisty
nbe_den <- ggplot(nbe.df) +                
  geom_density(aes(x = nbe, fill=scenario), alpha=0.2)+     
    theme_bw() +                                                 
  ylab("") +                                            
  xlab("\nNBE")  +                      
  theme(axis.text = element_text(size = 12),                   
        axis.title.x = element_text(size = 14, face = "plain"), 
        panel.grid = element_blank(),                           
        plot.margin = unit(c(1,1,1,1), units = , "cm"))

# Plotting fLuc denisty
fLuc_den<- ggplot(fLuc.df) +                
  geom_density(aes(x = fLuc, fill=scenario), alpha=0.2)+  
  scale_x_log10()+    
    theme_bw() +                                                 
  ylab("") +                                            
  xlab("\nFLuc")  +                      
  theme(axis.text = element_text(size = 12),                   
        axis.title.x = element_text(size = 14, face = "plain"), 
        panel.grid = element_blank(),                           
        plot.margin = unit(c(1,1,1,1), units = , "cm"))

# Plotting nbp denisty
nbp_den <- ggplot(nbp.df) +                
  geom_density(aes(x = nbp, fill=scenario), alpha=0.2)+      
    theme_bw() +                                                 
  ylab("") +                                            
  xlab("\nNBP")  +                      
  theme(axis.text = element_text(size = 12),                   
        axis.title.x = element_text(size = 14, face = "plain"), 
        panel.grid = element_blank(),                           
        plot.margin = unit(c(1,1,1,1), units = , "cm"))

grid.arrange(gpp_den, ra_den, npp_den, rh_den, nee_den, fire_den, nbe_den,fLuc_den,nbp_den,
                   # labels = c("A", "B", "C","D","E","F","G","H","I"),
                    ncol = 3, nrow = 3)
