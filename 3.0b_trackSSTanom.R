## Identify the anomaly extend and direction
## Track the growth/decay of the anomaly

library(raster)
library(rasterVis)
library(xts)
library(maptools)
library(gridExtra)
library(grid)

setwd("~/R_projects-CS/ame-temporalbabe/CS-ExtractFeatures/")

source("func_allDates.R")
source("func_fileList.R")
source("func_getstring.R")

# setup bathymetry map
bathy <- raster("/home/robert/R/x86_64-pc-linux-gnu-library/3.3/bathy-GEBCO/gebco-southernAfrica_1min.nc")
locbb <- matrix(c(10,-40,40,-20),nrow=2,ncol=2)
x <- extent(locbb)
bathy.cont <- crop(bathy,x)
bathy.loc <- crop(bathy,x)      # crop bathy map
bathy.land <- crop(bathy,x)
rm(bathy,locbb,x)

borders <- readShapePoly("SAfica-borders.shp")
borders <- SpatialPolygons(borders@polygons)
minlat <- -40;minlon <- 10;maxlat <- -20;maxlon <- 40

ncol <- 20
mapThemeQ10 <- rasterTheme(region = colorRampPalette(c("white","blue"))(ncol))
mapThemeQ90 <- rasterTheme(region = colorRampPalette(c("white","red"))(ncol))
mapThemeBathy <- colorRampPalette(c("blue","white"),alpha=T)
at <- c(6000,0,-200,-500,-1000,-2000,-3000,-4000,-5000,-6000,-7000)

crs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")
station.loc <- read.table("~/R_projects-CS/ame-temporalbabe/CS-ExtractFeatures/station_locations.Rdata")
station.pts <- SpatialPointsDataFrame(station.loc[,1:2],station.loc,proj4string = crs)

# ------------------------- SST
baseURL <- "/media/robert/Seagate Backup Plus Drive"
type <- "SST"
source <- "GHRSST"
level <- "L4"
product <- "AVHRR-OI"

list.files <- fileList(baseURL,type,source,level,product)
all.dates <- alldates(list.files)
all.doy <- as.numeric(strftime(all.dates, format = "%j"))

# open existing raster data
open.file1 <- paste0("Output/",product,"_anom.grd")
open.file2 <- paste0("Output/",product,"_100kmMask.grd")
r.anom <- brick(open.file1)
bathy.mask <- raster(open.file2)
rm(open.file1, open.file2)

# get the 10% and 90% quantiles
Q1090 <- calc(r.anom, fun = function(x) {quantile(x,probs = c(.1,.9),na.rm=TRUE)})

# plot the upper and lower quantiles ---------------------------------------------
# note that this does not yet print the legend title

filename1 <- paste0("Output/",product,"_AnomalyQuant.png")
png(filename1,width = 8.5,height = 5,units = 'in',res = 300)

leg.title <- gpar(fontsize=10,cex=1,lineheight=0.4,fontface="italic")

plot.title1 <- paste0(product," Anomaly Quantile 10%")
p1 <- levelplot(Q1090$layer.1, par.settings=mapThemeQ10,margin = F,
                  main = plot.title1, xlab="Latitude",ylab="Longitude") +
  layer(sp.lines(borders, col = "black", lwd = 1))
trellis.focus("legend",side = "right",clipp.off = TRUE,highlight = FALSE)
grid.text(expression(atop(paste("SST (",degree,"C)"),"below clim.")), 0.2, -0.1,
          rot = 0, gp=leg.title)
trellis.unfocus()

plot.title2 <- paste0(product," Anomaly Quantile 90%")
p2 <- levelplot(Q1090$layer.2, par.settings=mapThemeQ90, margin = F,
                  main = plot.title2, xlab="Latitude",ylab="Longitude") +
  layer(sp.lines(borders, col = "black", lwd = 1))
trellis.focus("legend",side = "right",clip.off = TRUE,highlight = FALSE)
grid.text(expression(atop(paste("SST (",degree,"C)"),"above clim.")), 0.2, -0.1,
          rot = 0, gp=leg.title)
trellis.unfocus()

print(p1, split=c(1,1,2,1))
print(p2, split=c(2,1,2,1), newpage=FALSE)

dev.off()

# if bathy.mask does not exist ---------------------------------------------------
if (!exists("bathy.mask")){
  # extract the coastal quantiles from the raster layer
  buf.dist <- 100000                                  # coastal boudary distance metres
  bathy.land[bathy.land >= 0] = NA                    # set land values to NA
  bathy.loc[bathy.loc > 5] <- NA                      # select 5 m height as coastal limits
  bathy.loc[bathy.loc < -5] <- NA
  bathy.buf <- buffer(bathy.loc, width=buf.dist, doEdge=TRUE)
  bathy.mask.temp <- mask(bathy.buf,bathy.land)       # remove the inland buffer distance
  
  Q1090@crs <- crs                                    # set the projection to bathymetry map
  bathy.mask <- projectRaster(bathy.mask.temp,Q1090)  # resample bathymetry to data resolution
  
  # save the mask image and mask raster layer
  filename1 <- paste0("Output/",product,"_100kmMask.grd")
  writeRaster(bathy.mask,filename1,"raster", overwrite=T)
}

# print mask map ------------------------------------------------------------------
bathy.mask.plot <- bathy.mask
bathy.mask.plot[is.na(bathy.mask.plot)] <- 0

filename1 <- paste0("Output/",product,"_100kmMask.png")
png(filename1,width = 8.5,height = 5,units = 'in',res = 300)

obj <- levelplot(bathy.mask.plot, margin=F, col.regions=c("white","yellow"), colorkey=F, 
          xlab="Longitude",ylab="Latitude", contour=T,
          main=paste(product," 100 km Coastal Boundary Mask")) +
contourplot(bathy.cont,region=T,col.regions=mapThemeBathy, alpha.regions=.3,
            at=at,margin = F, labels = F) +
  layer(sp.lines(borders, col = "black", lwd = 1))

print(obj)
dev.off()

rm(bathy.mask.plot)

# 2.0 extract 10 and 90 % anomaly limits --------------------------------------------
Q1090coast <- mask(Q1090$layer.1,bathy.mask)
names(Q1090coast) <- "Qlower"
Q1090coast$Qupper <- mask(Q1090$layer.2,bathy.mask)

r.anom.low.list <- list()
r.anom.hi.list <- list()

# build stack with
for (i in 1:nlayers(r.anom)){
  cat("Now on Day ",i,"\n")
  r.anom.low.list[i] <- r.anom[[i]] < Q1090coast$Qlower
  r.anom.hi.list[i] <- r.anom[[i]] > Q1090coast$Qupper
}

r.anom.low <- stack(r.anom.low.list)
r.anom.hi <- stack(r.anom.hi.list)

station.idx <- c(1,2,4,7,13,23,27,31,33,41,49,53,62,82,94,105,108)

# plot cold anomalies -----------------------------------------------------------------------
for (i in 367:nlayers(r.anom)){
  cat("Printing ",i," of ",nlayers(r.anom),"\n")
  
  filename1 <- paste0("Output/",strftime(all.dates[i],format="%Y%m%d"),".png")
  png(filename1,width = 8.5,height = 5,units = 'in',res = 300)
  
  plot.title1 <- paste0(product," < 10 % Anomaly\n",strftime(all.dates[i],format="%d %B %Y"))
  plot.title2 <- paste0(product," > 90 % Anomaly\n",strftime(all.dates[i],format="%d %B %Y"))
  
  obj1 <- levelplot(r.anom.low[[i]], margin=F, at=c(-0.5,0.5,1.5),
                   col.regions=c("old lace","blue"), colorkey=F, 
                   xlab="Longitude",ylab="Latitude", contour=F,
                   main=plot.title1) +
    layer(sp.lines(borders, col = "black", lwd = 1))
    #layer(sp.points(station.pts[station.idx,], col = "red", pch = 9))
  
  obj2 <- levelplot(r.anom.hi[[i]], margin=F, at=c(-0.5,0.5,1.5),
                   col.regions=c("old lace","darkorange1"), colorkey=F, 
                   xlab="Longitude",ylab=NULL, contour=F,
                   main=plot.title2) +
    layer(sp.lines(borders, col = "black", lwd = 1))
    #layer(sp.points(station.pts, col = "red", pch = 9))
  
  print(obj1, split=c(1,1,2,1))
  print(obj2, split=c(2,1,2,1), newpage=FALSE)
  
  dev.off()
}



