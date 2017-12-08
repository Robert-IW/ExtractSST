## Creates a monthly climatology and std dev from L4 data

## A detrended climatology is considered but not yet applied (10-05-2017)

library(raster)
library(rasterVis)
library(rgdal)
library(sp)
library(doParallel)
library(maptools)
library(lubridate)
library(spacetime)
library(mgcv)
library(foreach)
library(rgdal)

source("~/R/myFunctions/func_gls.R")
source("~/R/myFunctions/func_nonlingamm.R")
source("~/R/myFunctions/func_allDates.R")
source("~/R/myFunctions/func_fileList.R")
source("~/R/myFunctions/func_getstring.R")

setwd("~/R_projects-CS/ame-temporalbabe/ExtractSSTanomaly/")

#bathy <- raster("/home/robert/R/x86_64-pc-linux-gnu-library/3.3/bathy-GEBCO/gebco-southernAfrica.nc")
crs <- CRS("+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0")

# set up plotting
borders <- readShapePoly("~/R/myFunctions/SAfrica-borders.shp", proj4string = crs)
borders <- SpatialPolygons(borders@polygons)

mapTheme1 <- rasterTheme(region = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)))
mapTheme2 <- rasterTheme(region = colorRampPalette(brewer.pal(9, "YlGn"))(100))
mapTheme3 <- rasterTheme(region = colorRampPalette(c("green3","white","darkorange1"))(50))
mapTheme4 <- rasterTheme(region = colorRampPalette(c("white","blueviolet"))(50))
color.pal <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(256))

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

# ------------------------- SST
baseURL <- "/media/robert/KELP-HDD-Portable/"
type <- "SST"
source <- "GHRSST"
level <- "L4"
#product <- "AVHRR-OI"
#product <- "CMC"
#product <- "K10"
product <- "MUR"

# load the raw data
list.files <- fileList(baseURL,type,source,level,product)
all.dates <- alldates(list.files)
rm(list.files)

# set temoerture range of climatology
clim.1 <- 10
clim.2 <- 32

open.file1 <- paste0("Output/",product,"_data.grd")
r.stack <- brick(open.file1)
rm(open.file1)

############################################### Monthly CLIMATOLOGY
# produce the monthly climatolgy by month means and detrended month means

cat("Starting on Climatologies\n")

# visualize linear model slope coefficient and least sqaures error
times <- 1:nlayers(r.stack)
fun1 <- function(x) { if (is.na(x[1])){ NA } else { lm(x ~ times)$coefficients[2]}}
fun2 <- function(x) { if (is.na(x[1])){ NA } else { summary(lm(x ~ times))$sigma}}
r.slope <- calc(r.stack, fun1)
r.slope.err <- calc(r.stack,fun2)

save.file1 <- paste0("Output/",product,"_slope.grd")
save.file2 <- paste0("Output/",product,"_slopeErr.grd")

writeRaster(r.slope, save.file1, "raster", overwrite=T)
writeRaster(r.slope.err, save.file2, "raster",overwrite=T)
rm(save.file1,save.file2)

# plot the slope and errors -----------------------------------------------------
# prepare colors to centre white on zero for slope plot
cat("Plotting regression graphs\n")

filename1 <- paste0("Output/",product,"_LinearRegress.png")
png(filename1,width = 8.5,height = 5,units = 'in',res = 300)

paletteLength <- 50
myBreaks <- c(seq(minValue(r.slope), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(maxValue(r.slope)/paletteLength, maxValue(r.slope), length.out=floor(paletteLength/2)))

plot.title1 <- paste0(product," Linear Regression\nSlope Coefficient")
obj1 <- levelplot(r.slope, at=myBreaks, par.settings = mapTheme3,
                  margin = F,main = plot.title1,
                  xlab="Latitude",ylab="Longitude") +
  layer(sp.lines(borders, col = "black", lwd = 1))

plot.title2 <- paste0(product," Linear Regression\nResidual Std. Err.")
obj2 <- levelplot(r.slope.err, par.settings = mapTheme4,
                  margin = F,main = plot.title2,
                  xlab="Latitude",ylab="Longitude") +
  layer(sp.lines(borders, col = "black", lwd = 1))

print(obj1, split=c(1,1,2,1))
print(obj2, split=c(2,1,2,1), newpage=FALSE)

dev.off()

# ------------------------------------------------------------------------------------
cat("Creating climatology graphs\n")

plot.title1 <- paste(product,level,"Climatologies\n ",all.dates[[1]],"to",all.dates[[length(all.dates)]])
plot.title2 <- paste(product,level,"Std. Deviation\n ",all.dates[[1]],"to",all.dates[[length(all.dates)]])

idxM <- as.integer(strftime(all.dates,"%m"))                  # get a list of months
idxDOY <- as.integer(strftime(all.dates,"%j"))
idxY <- as.integer(strftime(all.dates,"%Y"))
idxNo <- 1:nlayers(r.stack)
rm(all.dates)

# create data frame from the raster brick
temp <- data.frame(rasterToPoints(r.stack))
r.df <- list()

cat("Beginning parallel processing\n")
cl<-makeCluster(4)
registerDoParallel(cl)

r.df <- foreach(i=1:nrow(temp), .inorder=T) %dopar% {
  #cat("Now on ",i," of ",nrow(temp),"\n")
    r.df.temp <- data.frame("d.o.y"=idxDOY,"year"=idxY,"id"=idxNo,"sst"=t(temp[i,3:ncol(temp)]))
    names(r.df.temp)[4] <- "sst"
    r.df.temp
}
stopCluster(cl)

cat("Saving data\n")
save.file3 <- paste0("Output/",product,"_listDFcells.Rdata")
save(r.df, file = save.file3)
rm(save.file3, r.df,  temp)
gc()

# before detrending
r.clim.temp <- stackApply(r.stack, idxM, median, na.rm=TRUE)    # by means of months

# # after detrending
# cl=makeCluster(6)
# clusterEvalQ(cl,c(library(mgcv)))     # Load package on all instances of R on all cores
# 
# r.gamm <- parLapply(cl,r.df,function(i) {
#   temp <- gamm(sst ~ s(d.o.y, bs = "cc", k = 365) + s(id, bs = "cr"),   # k is the no. of unique values
#        correlation = corARMA(form = ~ 1 | year, p=1),                 # note p=3 better, p=2 error
#        method = "REML",
#        data = i, na.action = na.omit)
#   return(temp$lme)
# })
# stopCluster(cl)

# get order of months produced by stackApply and re-order months 1-12
if (idxM[1] != 1){
  temp.ord <- c(seq(idxM[1],12,1),seq(1,idxM[1]-1,1))
  idxJan <- which(temp.ord == 1)
  data.ord <- c(seq(idxJan,12,1),seq(1,idxJan-1,1))
} else {
  data.ord <- seq(1,12,1)
}

r.clim <- subset(r.clim.temp, order(data.ord))
names(r.clim) <- month.name
rm(r.clim.temp)

filename.png1 <- paste0("Output/",product,"_climato.png")

png(filename.png1, width = 8.5,height = 5,units = 'in',res = 300)

obj <- levelplot(r.clim, at=seq(clim.1,clim.2,length=50), par.settings = mapTheme1,
                 margin = F,main = plot.title1,
                 xlab="Latitude",ylab="Longitude") +
  layer(sp.lines(borders, col = "black", lwd = 1))# +

print(obj)
dev.off()
rm(obj)
gc()

############################################### Monthly STD DEVIATION
# produce the standard deviation
cat("Estimating Standard Deviation\n")
r.sd.temp <- stackApply(r.stack, idxM, sd, na.rm=TRUE)

r.sd <- subset(r.sd.temp, order(data.ord))
names(r.sd) <- month.name
rm(r.sd.temp)

filename.png2 <- paste0("Output/",product,"_stddev.png")

png(filename.png2,width = 8.5,height = 5,units = 'in',res = 300)

obj <- levelplot(r.sd, par.settings = mapTheme2, 
                 margin = F,main = plot.title2,
                 xlab="Latitude",ylab="Longitude") +
  layer(sp.lines(borders, col = "black", lwd = 1))# +

print(obj)
dev.off()
rm(obj)
gc()

# save the raster stack
cat("Saving the raster stack\n")
save.file1 <- paste0("Output/",product,"_clim.grd")
save.file2 <- paste0("Output/",product,"_std.grd")

writeRaster(r.clim, save.file1, "raster", overwrite=T)
rm(r.clim)
writeRaster(r.sd, save.file2, "raster", overwrite=T)
rm(r.sd)

rm(mapTheme1,mapTheme2,mapTheme3,mapTheme4)