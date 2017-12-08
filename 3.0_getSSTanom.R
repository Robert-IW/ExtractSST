## Opens the data and climatology raster bricks created by '1.0_getDataStack.R'
## Creates and saves or opens the anomaly stack '*_anom.grd'
## What does it output and why?

library(raster)
library(rasterVis)
library(xts)
library(maptools)
library(grid)

#setwd("~/R_projects-CS/ame-temporalbabe/CS-ExtractFeatures/")
doGraphic <- 0

source("~/R/myFunctions/func_allDates.R")
source("~/R/myFunctions/func_fileList.R")
source("~/R/myFunctions/func_getstring.R")

# ------------------------- SST
baseURL <- "/media/robert/KELP-HDD-Portable/"
type <- "SST"
source <- "GHRSST"
level <- "L4"
#product <- "AVHRR-OI"
#product <- "CMC"
#product <- "K10"
product <- "MUR"

list.files <- fileList(baseURL,type,source,level,product)
all.dates <- alldates(list.files)
all.doy <- as.numeric(strftime(all.dates, format = "%j"))

open.file1 <- paste0("Output/",product,"_clim.grd")
open.file2 <- paste0("Output/",product,"_data.grd")
# open.file3 <- paste0("Output/",product,"_anom.grd")

r.clim <- brick(open.file1)
r.brick <- brick(open.file2)
# if (exists(open.file3)){
#   r.anom <- brick(open.file3)
# }

# interpolate day-of-the-year climatology for smooth anomalies
date.seq <- seq.Date(as.Date("1971-12-01"),as.Date("1973-01-01"), by = "day")
mon.loc <- axTicksByTime(date.seq,"months",format.labels="%b")
data.loc <- as.vector(mon.loc[1:14])

if (product=="AVHRR-OI" | product=="CMC"){
  r <- raster(nrows=nrow(r.clim), ncols=ncol(r.clim),resolution=c(0.25,0.25),ext=r.clim@extent)
} else if (product=="K10"){
  r <- raster(nrows=nrow(r.clim), ncols=ncol(r.clim),resolution=c(0.1,0.1),ext=r.clim@extent)
}

values(r) <- NA
x <- sapply(1:397, function(...) r)
rm(r)

x[data.loc] <- c(r.clim[[12]],r.clim[[1]],r.clim[[2]],
                 r.clim[[3]],r.clim[[4]],r.clim[[5]],
                 r.clim[[6]],r.clim[[7]],r.clim[[8]],
                 r.clim[[9]],r.clim[[10]],r.clim[[11]],
                 r.clim[[12]],r.clim[[1]])

s <- stack(x)
r.interp <- approxNA(s)
rm(x,s)

if (doGraphic){
  
  borders <- readShapePoly("~/R/myFunctions/SAfrica-borders.shp")
  borders <- SpatialPolygons(borders@polygons)
  minlat <- -40;minlon <- 10;maxlat <- -20;maxlon <- 40
  
  # plot each day average
  mincol <- 10;maxcol <- 32;brk <- 0.25;ncol <- (maxcol-mincol)/brk
  mapTheme1 <- rasterTheme(region = rev(colorRampPalette(c("red","yellow","springgreen","royalblue"))(ncol)))
  leg.title <- gpar(fontsize=10,cex=1,lineheight=0.4,fontface="italic")
  
  for (i in 1:nlayers(r.interp)){
    plot.title1 <- paste0(product," Daily Climatology\n",strftime(date.seq[i],format="%d %B"))
    
    filename1 <- paste0("Output/",i,".png")
    png(filename1,width = 8.5,height = 5,units = 'in',res = 300)
    
    obj <- levelplot(r.interp[[i]], at=seq(mincol,maxcol,length=ncol), par.settings=mapTheme1,
                     colorkey=list(labels=list(at = c(10,15,20,25,30))),
                     margin = F,main = plot.title1, 
                     xlab="Latitude",ylab="Longitude") +
      layer(sp.lines(borders, col = "black", lwd = 1))
    print(obj)
    trellis.focus("legend",side = "right",clip.off = TRUE,highlight = FALSE)
    grid.text(expression(paste("SST (",~degree,"C)", sep = "")), 1.4, 0.50, rot = 90,gp=leg.title)
    trellis.unfocus()
    dev.off()
  }
}

# subset data from 01 Jan to 31 Dec
r.clim.daily <- subset(r.interp,32:397)

# get the anomalies
r <- raster(nrows=nrow(r.clim), ncols=ncol(r.clim),resolution=c(0.25,0.25),ext=r.clim@extent)
values(r) <- NA
x <- sapply(1:nlayers(r.brick), function(...) r)
rm(r)

for (i in 1:nlayers(r.brick)){
  x[i] <- r.brick[[i]] - r.clim.daily[[all.doy[i]]]
}

r.anom <- stack(x)
names(r.anom) <- names(r.brick)
save.file1 <- paste0("Output/",product,"_anom.grd")
writeRaster(r.anom, save.file1, "raster", overwrite=T)

# plot each day average
mincol <- -10;maxcol <- 10;brk <- 0.25;ncol <- (maxcol-mincol)/brk
mapTheme1 <- rasterTheme(region = colorRampPalette(c("blue","white","red"))(ncol))

if (doGraphic){
  
  borders <- readShapePoly("~/R/myFunctions/SAfrica-borders.shp")
  borders <- SpatialPolygons(borders@polygons)
  minlat <- -40;minlon <- 10;maxlat <- -20;maxlon <- 40
  
  for (i in 1:nlayers(r.anom)){
    
    plot.title1 <- paste0(product," Daily Anomaly\n",strftime(all.dates[i],format="%d %B %Y"))
    
    filename1 <- paste0("Output/",strftime(all.dates[i],format="%Y%m%d"),".png")
    png(filename1,width = 8.5,height = 5,units = 'in',res = 300)
    
    obj <- levelplot(r.anom[[i]], at=seq(mincol,maxcol,length=ncol), par.settings=mapTheme1,
                     colorkey=list(labels=list(at = c(-10,-5,0,5,10))),
                     margin = F,main = plot.title1, 
                     xlab="Latitude",ylab="Longitude") +
      layer(sp.lines(borders, col = "black", lwd = 1))
    print(obj)          # this is required to run the loop
    trellis.focus("legend",side = "right",clipp.off = TRUE,highlight = FALSE)
    grid.text(expression(paste("Anomaly (",~degree,"C)", sep = "")), 1.4, 0.50, rot = 90)
    trellis.unfocus()
    dev.off()
  }
}
