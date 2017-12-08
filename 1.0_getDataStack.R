## To load all the L4 satellite data as raster stack

library(raster)
library(rasterVis)
library(rgdal)
library(sp)
library(parallel)
library(maptools)
library(lubridate)
library(spacetime)

source("~/R/myFunctions/func_fileList.R")
source("~/R/myFunctions/func_csv2pts.R")
source("~/R/myFunctions/func_startenddates.R")
source("~/R/myFunctions/func_getstring.R")
source("~/R/myFunctions/func_allDates.R")

######################################################################################################
# select the data product
# ----------------------- SLP
#baseURL <- "/media/robert/Seagate Backup Plus Drive"
#type <- "SLP"
#source <- "NOAA"
#level <- "L4"
#product <- "NCEP"

# get a single list of all the data and the start and end dates of the sequence using local functions
#list.files <- fileList(baseURL,type,source,level,product)
#dates <- startdateend(list.files)

#plot.title1 <- paste(product,level,"Climatologies\n ",dates[[1]],"to",dates[[2]])
#plot.title2 <- paste(product,level,"Std. Deviation\n ",dates[[1]],"to",dates[[2]])
#save.file <- "Output/SLPstack.Rdata"

# ------------------------- SST
#baseURL <- "/media/robert/KELP-HDD-Portable/"
baseURL <- "/media/robert/Seagate Backup Plus Drive/Backup Satellite Data/"
type <- "SST"
source <- "GHRSST"
level <- "L4"
#product <- "AVHRR-OI"
#product <- "CMC"
#product <- "K10"
product <- "MUR"

# get a single list of all the data and the start and end dates of the sequence using local functions
list.files <- fileList(baseURL,type,source,level,product)   # 'fileList' is a function
all.dates <- alldates(list.files)  
dates <- startdateend(list.files)

# convert csv to raster layer and stack layers using 'csv2pts' and 'pts2stack'
setwd("/media/robert/KELP-HDD-Portable/SST/GHRSST/L4/MUR/tempfiles")
pts <- csv2pts(list.files)

if (product=="AVHRR-OI"){                             # for AVHRR
  r.list <- lapply(pts, function(f) {rasterFromXYZ(f)})
  
} else if (product=="CMC") {                          # for CMC create a grid as cells are irregular
  r <- raster(xmn=min(pts[[1]][,1]),xmx=max(pts[[1]][,1]),ymn=min(pts[[1]][,2]),ymx=max(pts[[1]][,2]),
              resolution=0.25)
  r.list <- mclapply(pts, function(f) rasterize(f[, 1:2], r, f[,3], fun=mean), mc.cores = 5)
} else if (product=="K10") {                          # for CMC create a grid as cells are irregular
  r <- raster(xmn=min(pts[[1]][,1]),xmx=max(pts[[1]][,1]),ymn=min(pts[[1]][,2]),ymx=max(pts[[1]][,2]),
              resolution=0.1)
  r.list <- mclapply(pts, function(f) rasterize(f[, 1:2], r, f[,3], fun=mean), mc.cores = 5)
}

#r.stack <- stack(r.list)
#names(r.list) <- as.character(all.dates)
r.brick <- brick(r.list)
rm(pts,r.list)
gc()

# get the dates of missing data
idx <- seq(as.Date(dates[[1]]), as.Date(dates[[2]]), 'day') # sequence from first to last date                         # dates from file names
temp <- is.na(match(idx,all.dates,nomatch=NA))
miss.dates <- idx[which(temp == TRUE)]
miss.dates.loc <- which(temp == TRUE)

filename.temp <- paste0("Output/",product,"_missing.Rdata")
save(miss.dates,miss.dates.loc,file = filename.temp)
rm(temp,miss.dates,miss.dates.loc,filename.temp)

# change '99' value to NA
#r.stack[r.stack==99] <- NA
r.brick[r.brick==99] <- NA

names(r.brick) <- as.character(all.dates)
rm(all.dates)
#names(SISmm) <- month.abb

# write the data stack
save.file <- paste0("Output/",product,"_data.grd")
#writeRaster(r.stack, save.file, "raster", overwrite=T)
writeRaster(r.brick, save.file, "raster", overwrite=T)
#rm(r.stack)
