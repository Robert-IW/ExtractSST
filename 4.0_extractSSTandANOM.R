## Opens the sst  (AVHRR-OI_data.grd) and anomaly (AVHRR-OI_anom.grd) raster bricks created by '3.0_getSSTanom.R'
## Uses the 2, 1 and 0.5 deg polygons created by OAFlux '1.0_createBB.R' to extract:
##      mean SST and anomaly
##      std dev SST and anomaly

library(raster)
library(rasterVis)
library(xts)
library(maptools)
library(grid)
library(data.table)
library(parallel)

#setwd("~/R_projects-CS/ame-temporalbabe/CS-ExtractFeatures/")
doGraphic <- 0
doSST <- 0
doANOM <- 0
doCLIM <- 1

source("~/R/myFunctions/func_allDates.R")
source("~/R/myFunctions/func_fileList.R")
source("~/R/myFunctions/func_getstring.R")

# ------------------------- SST
baseURL <- "/media/robert/KELP-HDD-Portable/"
type <- "SST"
source <- "GHRSST"
level <- "L4"
product <- "AVHRR-OI"; resol <- 25
#product <- "CMC"; resol <- 25
#product <- "K10"; resol <- 10
#product <- "MUR"; resol <- 1

crs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

list.files <- fileList(baseURL,type,source,level,product)
all.dates <- alldates(list.files)
#all.doy <- as.numeric(strftime(all.dates, format = "%j"))

open.file1 <- paste0("Output/",product,"_data.grd")
open.file2 <- paste0("Output/",product,"_anom.grd")
open.file3 <- paste0("Output/",product,"_clim.grd")

r.sst <- brick(open.file1)
r.anom <- brick(open.file2)
r.clim <- brick(open.file3)

all.dates <- names(r.sst)
all.dates <- gsub("X","",all.dates)
all.dates <- gsub("\\.","",all.dates)

# convert date to numeric month of the year
all.months <- as.numeric(format(as.Date(all.dates, format = "%Y%m%d"),"%m"))

# get station locations
station.loc <- read.table("~/R/myFunctions/station_locations.Rdata",stringsAsFactors = F)
station.pts <- SpatialPointsDataFrame(station.loc[,1:2],station.loc,proj4string = crs)

for (h in 1:nrow(station.loc)){
  
  stat.name <- station.loc$station[h]
  cat(paste0("Starting bbox for ",stat.name,'\n'))
  
  st <- station.pts@data[h,1:2]
  coordinates(st) <- ~lon+lat
  
  # load the spatial polygons for each station
  load(paste0("~/R/myFunctions/Trans_Polys/",stat.name,"_OAfluxPolygons.Rdata"))
  
  # create a box for the nearest
  if (resol<2){
    resol <- resol+5
  }
  latEdge <- 2*resol* (1/110.574)
  lonEdge <-2*resol* (1/(111.320*cos(stLat*pi/180)))
  
  sub.df <- subset(myData, lat <= stLat+latEdge & lat >= stLat-latEdge 
                   & lon <= stLon+lonEdge & lon >= stLon-lonEdge)
  
  # ---------------------------------------------------------------------------------
  if (doSST){
    cat("Extracting SST data\n")
    
    tasks <- list(sst.median2 = function() as.vector(raster::extract(r.sst, bb.2sp, fun=median, na.rm=T)),
                  sst.median1 = function() as.vector(raster::extract(r.sst, bb.1sp, fun=median, na.rm=T)),
                  sst.median05 = function() as.vector(raster::extract(r.sst, bb.05sp, fun=median, na.rm=T)),
                  sst.stdDev2 = function() as.vector(raster::extract(r.sst, bb.2sp, fun=sd, na.rm=T)),
                  sst.stdDev1 = function() as.vector(raster::extract(r.sst, bb.1sp, fun=sd, na.rm=T)),
                  sst.stdDev05 = function() as.vector(raster::extract(r.sst, bb.05sp, fun=sd, na.rm=T)))
    
    out <- mclapply( 
      tasks, 
      function(f) f(), 
      mc.cores = 6
    )
    
    df.sst <- data.frame("date"=all.dates,"lon"=station.pts[[1]][h],"lat"=station.pts[[2]][h],
                         "med2deg"=out[[1]],"std2deg"=out[[4]],
                         "med1deg"=out[[2]],"std1deg"=out[[5]],
                         "med05deg"=out[[3]],"std05deg"=out[[6]])
    
    rm(out, tasks)
    
    cat("Writing Files\n")
    filename1 <- paste0("/media/robert/KELP-HDD-Portable/SST/GHRSST/L4/",product,"/extractCSV/",
                        stat.name,"_SST_",product,".csv")
    fwrite(df.sst, file.path = filename1, col.names=T)
  }
  # --------------------------------------------------------------------
  if (doANOM){
    cat("Extracting ANOM data\n")
    
    tasks <-  list(anom.median2 = function() as.vector(raster::extract(r.anom, bb.2sp, fun=median, na.rm=T)),
                   anom.median1 = function() as.vector(raster::extract(r.anom, bb.1sp, fun=median, na.rm=T)),
                   anom.median05 = function() as.vector(raster::extract(r.anom, bb.05sp, fun=median, na.rm=T)),
                   anom.stdDev2 = function() as.vector(raster::extract(r.anom, bb.2sp, fun=sd, na.rm=T)),
                   anom.stdDev1 = function() as.vector(raster::extract(r.anom, bb.1sp, fun=sd, na.rm=T)),
                   anom.stdDev05 = function() as.vector(raster::extract(r.anom, bb.05sp, fun=sd, na.rm=T)))
    
    out <- mclapply( 
      tasks, 
      function(f) f(), 
      mc.cores = 6
    )
    
    # create the climatology based on the all.dates months
    df.anom <- data.frame("date"=all.dates,"lon"=station.pts[[1]][h],"lat"=station.pts[[2]][h],
                          "median2deg"=out[[1]],"std2deg"=out[[4]],
                          "median1deg"=out[[2]],"std1deg"=out[[5]],
                          "median05deg"=out[[3]],"std05deg"=out[[6]])
    
    rm(out)
    
    cat("Writing Files\n")
    filename2 <- paste0("/media/robert/KELP-HDD-Portable/SST/GHRSST/L4/",product,"/extractCSV/",
                        stat.name,"_ANOM_",product,".csv")
    fwrite(df.anom, file.path = filename2, col.names=T)
  }
  # --------------------------------------------------------------------
  if (doCLIM){
    cat("Extracting CLIM data\n")
    
    # this extracts the block climatology for each month 1-12
    tasks <-  list(clim.median2 = function() as.vector(raster::extract(r.clim, bb.2sp, fun=median, na.rm=T)),
                   clim.median1 = function() as.vector(raster::extract(r.clim, bb.1sp, fun=median, na.rm=T)),
                   clim.median05 = function() as.vector(raster::extract(r.clim, bb.05sp, fun=median, na.rm=T)),
                   clim.stdDev2 = function() as.vector(raster::extract(r.clim, bb.2sp, fun=sd, na.rm=T)),
                   clim.stdDev1 = function() as.vector(raster::extract(r.clim, bb.1sp, fun=sd, na.rm=T)),
                   clim.stdDev05 = function() as.vector(raster::extract(r.clim, bb.05sp, fun=sd, na.rm=T)))
    
    out <- mclapply( 
      tasks, 
      function(f) f(), 
      mc.cores = 6
    )
    
    median2deg <-  vector()
    median1deg <-  vector()
    median05deg <-  vector()
    std2deg <-  vector()
    std1deg <-  vector()
    std05deg <-  vector()
    
    for (t in 1:length(all.months)){
      median2deg[t] <- out[[1]][[all.months[t]]]
      median1deg[t] <- out[[2]][[all.months[t]]]
      median05deg[t] <- out[[3]][[all.months[t]]]
      std2deg[t] <- out[[4]][[all.months[t]]]
      std1deg[t] <- out[[5]][[all.months[t]]]
      std05deg[t] <- out[[6]][[all.months[t]]]
    }
    
    df.clim <- data.frame("date"=all.dates,"lon"=station.pts[[1]][h],"lat"=station.pts[[2]][h],
                          "median2deg"=median2deg,"std2deg"=std2deg,
                          "median1deg"=median1deg,"std1deg"=std1deg,
                          "median05deg"=median05deg,"std05deg"=std05deg)
    rm(out)
    
    cat("Writing Files\n")
    filename3 <- paste0("/media/robert/KELP-HDD-Portable/SST/GHRSST/L4/",product,"/extractCSV/",
                        stat.name,"_CLIM_",product,".csv")
    
    fwrite(df.clim, file.path = filename3, col.names=T)
  }
  
  rm(df.sst, df.anom)
  gc()
  
}