## To match-up in situ stations with nearest L4 satellite data for each day
## The output is a list for each sensor (e.g. 'lst_MODISA_L3/4') of data frames for each day

## The climatologies can be determined from the extracted nearest values using decompose

## A series of data frame are created for each satellite product
## Individual station data frames are created and saved as csv for the entire time series

## Requires the csv files produced by 'extractSatSSTlevel[3C/4]_csv.R'

## Requires 'getInsituRecords.R' to be run first (or load 'OutPut01...Rdata')
## Requires 'list_dates', 'row_stations' 'col_time' from 'getInsituRecords.R'
## Calls 'callSatDataList.R'
## Calls the set of 'callSENSORSatVal.R' scripts

#################################################################################
# WARNING: Check that rounding 23:46 doesn't end up as 00:00
# anything beyond 11:45 should be rounded to 11:30
# code contains subsetting for 'station' only data (line 46)
#################################################################################
library(parallel)
library(Imap)
library(ncdf4)
library(tools)
library(xts)
library(maptools) # also loads sp package
library(data.table)

#setwd("~/R_projects-UWC/ExtractTS/")
setwd("~/R_projects-CS/ame-temporalbabe/ExtractSSTanomaly/")

# FUNCTION: get the date from the filename
getstr = function(mystring, initial.character, final.character) {
  sub(sprintf(".*%s(.*)%s.*", initial.character, final.character), "\\1", mystring)
}
# SOURCE: function
source("~/R/myFunctions/func_getNearestPlus.R")

# LOAD: get station information
#load("~/R_projects-UWC/ExtractTS/Output/OutPut01_InsituRecords.RData")
#load("~/R_projects-UWC/ExtractTS/Output/all_samples_stations.Rdata") # created by 'bind_all_stations'
#my.data <- all_samples

#---------------------------------------------------------------------------------------------
## THIS IS TO GET THE NEAREST FOR 'STATION' ONLY
#my.data <- subset(my.data, rownames(my.data)=="station")
station.loc <- read.table("~/R/myFunctions/station_locations.Rdata",stringsAsFactors = F)
#---------------------------------------------------------------------------------------------

# Get the coordinates
#source("func_perpCont.R")
#sourceURL <-("/media/robert/Seagate Backup Plus Drive/SST/GHRSST/")
sourceURL <-("/media/robert/Seagate Backup Plus Drive/Backup Satellite Data/SST/GHRSST/")
#sourceURL <- ("/media/robert/KELP-HDD-Portable/SST/GHRSST/")

# create directory strings for file locations
product <- c(
  "AVHRR_L3",
  "AVHRR-OI",
  "CMC",
  "G1SST",
  "K10",
  "MUR",
  "MW-IR-v1",
  "MW-IR-v4",
  "ODYSSEA",
  "OSTIA")

sourceCSV <- c(
  "L3C/AVHRR_L3/csv",
  "L4/AVHRR-OI/csv",
  "L4/CMC/csv",
  "L4/G1SST/csv",
  "L4/K10/csv",
  "L4/MUR/csv",
  "L4/MW_IR_v1.0/csv",
  "L4/MW_IR_v4.0/csv",
  "L4/ODYSSEA_SAF/csv",
  "L4/OSTIA/csv/")                            # directories where csv files kept

saveTS <- c(
  "L3C/AVHRR_L3/ts/",
  "L4/AVHRR-OI/ts/",
  "L4/CMC/ts/",
  "L4/G1SST/ts/",
  "L4/K10/ts/",
  "L4/MUR/ts/",
  "L4/MW_IR_v1.0/ts/",
  "L4/MW_IR_v4.0/ts/",
  "L4/ODYSSEA_SAF/ts/",
  "L4/OSTIA/ts/")                             # directories to save the time series

sensorRes <- c(4,25,20,1,10,1,9,9,10,5)       # vector of product pixel resolution

dataSources <-data.frame(product,sourceCSV,saveTS,sensorRes,stringsAsFactors = FALSE)
rm(product,sourceCSV,saveTS,sensorRes)

for (h in 5:5){#:nrow(dataSources)) {            # for each satellite product
  
  # create empty list for all the days data
  mylist <- sapply(station.loc[,3],function(x) NULL)
  
  productIndex <- h
  product <- dataSources[productIndex,1]
  sourceData <- dataSources[productIndex, 2]
  saveTS <- dataSources[productIndex, 3]
  resol <- dataSources[productIndex,4]
  
  #setwd(paste0("/media/robert/KELP-HDD-Portable/SST/GHRSST/L4/",product,"/"))

  # get a list of all the csv files
  fileURL <- list.files(paste(sourceURL, sourceData, sep = ''), full.names = TRUE)
  fileList <- list.files(paste(sourceURL, sourceData, sep = ''), recursive = FALSE)
  
  # open the climatology
  # open.file3 <- paste0("Output/",product,"_clim.grd")
  # r.clim <- brick(open.file3)
  
  ########## FOR EACH DAY 'i'
  for (i in 1:length(fileList)){
  
    # Get a list of files for this date for each sensor
    cat(paste("\nStarting on day", i,"\n"))
    
    # get the date from the filename
    if (product=="AVHRR-OI"){
      dateStr <- getstr(fileList[i],"_",".csv")
      temp_date <- dateStr
      #temp_date <- as.Date(dateStr,"%Y%m%d")
    } else {
      dateStr <- getstr(fileList[i],"_","_")
      temp_date <- dateStr
      #temp_date <- as.Date(dateStr,"%Y%m%d")
    }
    
    # # check if the data file is for night or day
    # if (grepl("day",fileList[i])){
    #   don <- "day"
    # } else if (grepl("night",fileList[i])){
    #   don <- "night"
    # } else {
    #   don <- "0"
    # }
    # 
    # get the month value for the climatology
    # mymonth <- as.numeric(format(as.Date(temp_date, format = "%Y%m%d"),"%m"))
    
    # get the corresponding month raster layer and convert to data.table
    # myanom <- r.clim[[mymonth]]
    # myClim <- as.data.frame(rasterToPoints(myanom))
    # colnames(myClim) <- c("lon","lat","clim")
    
    # open the csv file
    myData <- data.table::fread(fileURL[i])
    
    # for each station get the nearest data points dependent on resolution
    temp.df <- mcmapply(get_5nearest, station.loc[,1], station.loc[,2],station.loc[,3],
                        mc.cores = 5L,SIMPLIFY=FALSE)
    
    # append to mylist
    mylist <- mapply(rbind, mylist, temp.df, SIMPLIFY=FALSE)
    
    rm(dateStr,temp_date, temp.df)
    
  } # for each day 'i'
  
  cat("Writing csv file\n")
  
  # create each station as data frame from list and save as csv
  for (k in 1:length(mylist)){
    #myFrame <- do.call(rbind, c(mylist[[k]], make.row.names=FALSE))     # combine all the daily data into dataframe
    myFrame <- mylist[[k]]
    #myFrame <- as.data.frame(do.call(rbind, myFrame))

    filename1 <- paste("/media/robert/KELP-HDD-Portable/SST/GHRSST/",
                       saveTS,product,"_",
                       station.loc[k,3],
                       "_SST_timeseries_5nearest.csv",sep = '')
    
    data.table::fwrite(myFrame, file = filename1)
    rm(myFrame, filename1)
    gc()
  } # for each 'k'
} # for each 'h'
rm(h,i,k,mylist,sourceURL,dataSources,saveTS,product,resol,productIndex,fileURL,fileList,my.data)
