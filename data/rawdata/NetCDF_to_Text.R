#______________________________________________________________________________#
#Script for converting NetCDF4 files to Dataframe. 
#Reason - Issues with loading the package in Habanero. The Columbia Cluster



#______________________________________________________________________________#
#Setting Working Directory
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
library(ncdf4) 


#______________________________________________________________________________#
#Wind Power
name <- c("ERCOT_daily_WP_0_5_deg.nc")
nc_data <- nc_open(name)
print(nc_data) 
lat_lon_index <- ncvar_get(nc_data, "lat_lon_index")
t <- ncvar_get(nc_data, "time")
WP_CF <- ncvar_get(nc_data, "WP") 
WP_CF <- t(WP_CF) #Shifting to a correct format. 
nc_close(nc_data) #Closing the Netcdf file. 
write.table(WP_CF, "ERCOT_Wind_Power_Daily.txt", sep=" ")


#______________________________________________________________________________#
#Solar
name <- c("ERCOT_daily_RAD_0_5_deg.nc")
nc_data <- nc_open(name)
print(nc_data) 
lat_lon_index <- ncvar_get(nc_data, "lat_lon_index")
t <- ncvar_get(nc_data, "time")
ssrd <- ncvar_get(nc_data, "ssrd") 
ssrd <- t(ssrd) #Shifting to a correct format. 
nc_close(nc_data) #Closing the Netcdf file. 
write.table(ssrd, "ERCOT_Solar_Rad_Daily.txt", sep=" ")



#______________________________________________________________________________#
#Wind Power - Extended
name <- c("ERCOT_daily_WP_CF_0_5_deg_1950_2020.nc")
nc_data <- nc_open(name)
print(nc_data) 
lat_lon_index <- ncvar_get(nc_data, "lat_lon_index")
t <- ncvar_get(nc_data, "time")
WP_CF <- ncvar_get(nc_data, "WP CF") 
WP_CF <- t(WP_CF) #Shifting to a correct format. 
nc_close(nc_data) #Closing the Netcdf file. 
write.table(WP_CF, "ERCOT_Wind_Power_Daily_Extended.txt", sep=" ")


#______________________________________________________________________________#
#Solar - Extended
name <- c("ERCOT_daily_SP_CF_0_5_deg_1950_2020.nc")
nc_data <- nc_open(name)
print(nc_data) 
lat_lon_index <- ncvar_get(nc_data, "lat_lon_index")
t <- ncvar_get(nc_data, "time")
ssrd <- ncvar_get(nc_data, "SP CF") 
ssrd <- t(ssrd) #Shifting to a correct format. 
nc_close(nc_data) #Closing the Netcdf file. 
write.table(ssrd, "ERCOT_Solar_Rad_Daily_Extended.txt", sep=" ")
