####! To do:
# consider what species to lump
#     potentially based on whether species ID gets coarser in older years (to genus?)
#     plus based on ability to distinguish from air
# figure out the layer names of the flight path and mortality layers in older survey years GDBs
#     so they can be opened in a loop
# 
# Buffer climate data and master raster extent to a little further outside CA since flight path goes outside
# 
#Consider that due to reprojection, some distribution layer cells will have very low values for species because only a tiny fraction of the cell might contain the species

### Set your working directory to the "Aerial mortality" folder location
setwd("~/Aerial mortality")

#load libraries
library(sp)
library(rgdal)
library(raster)
library(ncdf4) # this is for TopoWx; you may need to install this following instructions here: http://cirrus.ucsd.edu/~pierce/ncdf/
library(latticeExtra)
library(gdalUtils)
library(spatial.tools)
library(plyr)
library(reshape2)
library(lme4)
library(ggplot2)
library(rgeos)
library(reshape2)

###### Utility functions #######
fn.bin <- function(x){ifelse(x >0,1,0)}

#############################################################
#### Grid-level variables (superceding year and species) ####
#############################################################


#### Create master rasters ####

# Open CA shapefile and project in Albers (meter units)
albers.proj <- CRS("+proj=aea +lat_1=34 +lat_2=40.5 +lat_0=0 +lon_0=-120 +x_0=0 +y_0=-4000000 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")
calif <- readOGR("California shapefile","statep010")
calif <- spTransform(calif,albers.proj) #project to albers (meter units)

# Define resolutions
res.fine <- 100 #spatial resolution for rasterization of flight and mortality data (in meters)
res.coarse <- 3500 #resolution for aggregation of mortality data and extraction of climate data for statistical analysis (should be somewhat finer than climate raster to avoid smoothing it out too much)
agg.factor <- round(res.coarse/res.fine) #aggregation factor for conversion from fine to coarse

# Create raster templates
master.fine <- raster(calif,res=res.fine) #define extent and coord system based on calif shapefile (shoudl be in albers--meters)
master.coarse <- aggregate(master.fine,agg.factor,FUN=mode) #coarse grid aligns exactly with fine grid

# # output sample
# values(master.fine) <- runif(length(values(master.fine)))
# writeRaster(master.fine,"master_fine_sample.tif")


# Extract data needed to re-create the same grid
master.fine.res <- res(master.fine)
master.fine.extent <- extent(master.fine)[c(1,3,2,4)]



##############################
####### Mortality loop #######
##############################

###################
#### Year loop ####
###################
year.set <- c(2009:2015)
######################
#### Species loop ####
######################
#Select species that have mortality data in 2014 and 2015
# file <- list.files(path="Survey GDBs",pattern=as.character(2015),full.names=TRUE)
# layers <- ogrListLayers(file)
# mortality.index <- grep("ADS",layers,ignore.case=TRUE)
# mort.polygon <- readOGR(dsn=file,layer=layers[mortality.index])
# focal.sp <- unique(c(mort.polygon$HOST1,mort.polygon$HOST2,mort.polygon$HOST3))
# focal.sp <- unique(c(focal.sp,746,41,312,17)) # add grand fir, bigleaf maple, port-orford cedar, quaking aspen (present in 2014 but not 2015)



## Make species list based on intersection of species found in all survey layers
mort.2009 <- readOGR("Survey shapefiles","ADS_2009")
mort.2010 <- readOGR("Survey shapefiles","ADS_2010")
mort.2011 <- readOGR("Survey shapefiles","ADS_2011")
mort.2012 <- readOGR("Survey shapefiles","ADS_2012")
mort.2013 <- readOGR("Survey shapefiles","ADS_2013")
mort.2014 <- readOGR("Survey shapefiles","ADS_2014")
mort.2015 <- readOGR("Survey shapefiles","ADS_2015")

sp.2009 <- unique(c(mort.2009$HOST1,mort.2009$HOST2,mort.2009$HOST3))
sp.2010 <- unique(c(mort.2010$HOST1,mort.2010$HOST2,mort.2010$HOST3))
sp.2011 <- unique(c(mort.2011$HOST1,mort.2011$HOST2,mort.2011$HOST3))
sp.2012 <- unique(c(mort.2012$HOST1,mort.2012$HOST2,mort.2012$HOST3))
sp.2013 <- unique(c(mort.2013$HOST1,mort.2013$HOST2,mort.2013$HOST3))
sp.2014 <- unique(c(mort.2014$HOST1,mort.2014$HOST2,mort.2014$HOST3))
sp.2015 <- unique(c(mort.2015$HOST1,mort.2015$HOST2,mort.2015$HOST3))

focal.sp <- intersect(sp.2009,(intersect(sp.2010,intersect(sp.2011,intersect(sp.2012,intersect(sp.2013,intersect(sp.2014,sp.2015)))))))


#Only include species that we have raster layers for distribution
file.distr <- list.files(path="Raster Live Basal Area Synced")
file.distr <- unique(substr(file.distr,2,nchar(file.distr)-4))
focal.sp <- as.list(focal.sp[focal.sp>=10 & focal.sp < 998])
focal.sp <- focal.sp[focal.sp %in% as.numeric(file.distr)]

species.set <- focal.sp

## OPTIONAL: add genera to species.set if necessary for consistency with older years



for(i in 1:length(year.set)){
  
  #define year of interest; and aerial survey file for that year
  year <- year.set[i]
  mort.file <- paste("ADS_",year,sep="")
  flight.file <- paste("Flown_area_",year,sep="")
  
  ## Open corresponding shapefiles
  mort.polygon <- readOGR("Survey shapefiles",mort.file)
  flight.polygon <- readOGR("Survey shapefiles",flight.file)
  
  ## Reproject to Albers (units in meters and compatible with master rasters)
  mort.polygon <- spTransform(mort.polygon,albers.proj)
  flight.polygon <- spTransform(flight.polygon,albers.proj)

  #Add column to flight shapefile that indicates that it was flown
  flight.polygon$FLOWN1 <- 1
  
  ## Eliminate polygons where mortality was due to fire
  fire <- 30000:30005
  mort.polygon <- mort.polygon[!((mort.polygon$DCA1 %in% fire) | (mort.polygon$DCA2 %in% fire) | (mort.polygon$DCA3 %in% fire)),]
  
  # Sum the three TPA columns whever they are positive
  tpa.cols <- cbind(mort.polygon$TPA1,mort.polygon$TPA2,mort.polygon$TPA3)
  tpa.cols[tpa.cols <= 0] <- NA
  mort.polygon$TPA.tot <- rowSums(tpa.cols,na.rm=TRUE)
  
  ## Eliminate polygons where total TPA is 1 or less
  mort.polygon <- mort.polygon[mort.polygon$TPA.tot > 1,]
  
  ## Eliminate polygons where the number of trees is 3 or less
  mort.polygon <- mort.polygon[mort.polygon$NO_TREES1 > 3 , ]

  #when multiple host species present, divide TPA by the number of focal sp.
  host.cols <- cbind((mort.polygon$HOST1), (mort.polygon$HOST2), (mort.polygon$HOST3))
  host.cols[host.cols < 1] <- NA
  n.host <- rowSums(host.cols > 0,na.rm=TRUE)
  count.dup <- function(x) {  sum(duplicated(x,incomparables=NA))  }
  n.dup <- apply(host.cols,1,count.dup)
  n.host.unique <- n.host - n.dup
  mort.polygon$TPA.split <- mort.polygon$TPA.tot/n.host.unique

  
  #write flight path polygon to shape file for rasterization
  writeOGR(flight.polygon, getwd(),"flightpolygon_test",driver="ESRI Shapefile",overwrite=TRUE)
  
  
  #######writeOGR(mort.polygon, getwd(),"mortpolygon_test2",driver="ESRI Shapefile",overwrite=TRUE)
  
  #rasterize it to the master.fine grid
  flightraster <- gdal_rasterize("flightpolygon_test.shp","flightrasterfocal_test.tif",
                                 a="FLOWN1", at=TRUE,tr=master.fine.res, te=master.fine.extent,
                                 l="flightpolygon_test",a_nodata=NA,verbose=TRUE,output_Raster=TRUE)
  
  #Aggregate flight raster to the scale for statistical analysis
  #! This sets cells all coarse cells that are partially outside the flight path to NA
  flightraster.coarse <- aggregate(flightraster,agg.factor,FUN=mode,na.rm=FALSE)
  
  
  
  for(j in 1:length(species.set)){

    #Define tree species groups 
    focal.sp <- species.set[[j]]
    
    ## Convert polygon TPA to TPA of focal species only (i.e., set TPA 0 if polygon does not contain any of focal species)
    # see if focal species is present
    focal.sp.present <- (mort.polygon$HOST1 %in% focal.sp) | (mort.polygon$HOST2 %in% focal.sp) | (mort.polygon$HOST3 %in% focal.sp)
    

    
    
    #only count TPAs when at least one of the mortality species is one of the focal species; otherwise set TPA to 0
    mort.polygon$tpa.focal <- ifelse(focal.sp.present,mort.polygon$TPA.split,0)
    
    # eliminate polygons with none of the current focal species
    if(sum(mort.polygon$tpa.focal) > 0) {
      mort.polygon.focal <- mort.polygon[mort.polygon$tpa.focal > 0,]
    } else {
      mort.polygon.focal <- mort.polygon[0,]
    }
    
    #### Rasterize mortality polygon; define unobserved cells, i.e. not in flight path, as NA ####
    #Using actual gdal code to rasterize; couldn't get rasterize function to generate sufficient sensitivity
    writeOGR(mort.polygon.focal, getwd(),"mortpolygonfocal_test", driver="ESRI Shapefile",overwrite=TRUE) #write mortality polygon of focal species to shape file for next step (rasterization)
    
    
    #! checkme: with if multiple polygons per raster cell? and does AT work as we think?
    mortraster <- gdal_rasterize("mortpolygonfocal_test.shp","mortrasterfocal_test.tif",
                                 a="tpa_fcl", at=TRUE,tr=master.fine.res, te=master.fine.extent,
                                 l="mortpolygonfocal_test",verbose=TRUE,output_Raster=TRUE)
    
    mortraster[is.na(mortraster)] <- 0 #this is a temp fix because we can't get the a_nodata to work above...
    
    #Aggregate to the scale for statistical analysis
    mortraster.coarse <- aggregate(mortraster,agg.factor,FUN=mean,na.rm=FALSE)
    
    #Set aggregated mortality cells that are at all outside flight path to NA
    mort.flight <- mortraster.coarse * flightraster.coarse
    

    ##########################################
    #### Prepare for statistical analysis ####
    ##########################################

    #stack rasters, name layers, and write to external file
    raster_stack <- stack(mort.flight)
    spgrp.text <- sprintf("%02d",j) # add leading zero

    layer.names <- paste("Y",year,".spgrp",spgrp.text,".",c("mort.tpa"),sep="")
    names(raster_stack) <- layer.names
    
    writeRaster(raster_stack, file=paste0("RasterOutput/Y",year,"_","spgrp",j,".grd",sep=""),overwrite=T)
    cat("\rFinished Year",year,", species group",j)
  }
}


###########################################
############ Species distr loop ###########
###########################################



for(i in 1:length(species.set)) {
  #Define tree species groups 
  focal.sp <- species.set[[i]]
  
  # Species distribution by live basal area raster #
  #check that species of interest have already been synced to master.coarse raster; if not run script at bottom of page
  # this code can open multiple rasters (for spp groups) and sums them

  distr.stack <- stack()
  
  for(j in 1:length(focal.sp)) {
    
    filename <- paste("Raster Live Basal Area Synced/s",focal.sp[j],".grd",sep="")
    
    # try to open it; if not, do not return an error; just skip it but print a warning
    # the purpose of the "try" component is to allow the user to include species groups that include genera,
    #       and since many genera do not exist in the live basal area rasters, they would throw an error
    
    sp.distr <- try( stack(filename), silent=TRUE )
    if(class(sp.distr) == "try-error") { # if it returned an error
      warning("Could not find the following species distribution raster: ",filename," Skipped it.",immediate.=TRUE)
      next() 
    } # otherwise it opened successfully
    
    distr.stack <- stack(distr.stack,sp.distr)
                         
  }

  if(nlayers(distr.stack) > 1) distr.stack <- sum(distr.stack)
  
  spgrp.text <- sprintf("%02d",i) # add leading zero
  layer.names <- paste("spgrp",spgrp.text,".distr",sep="")
  names(distr.stack) <- layer.names
  
  writeRaster(distr.stack,paste("RasterOutput/spgrp",spgrp.text,"_distr.grd",sep=""),overwrite=TRUE)

}


###########################################
############### Climate data ##############
###########################################

#### Open and prep climate layers ####

tmax.ann <- brick("Climate layers/PRISM_tmax_wateryear/PRISM_tmax_wateryear.grd")
ppt.ann <- brick("Climate layers/PRISM_ppt_wateryear/PRISM_ppt_wateryear.grd")
climate.rast <- stack(tmax.ann,ppt.ann)

# interpolate them down to the grid used for statistical analysis (master.coarse)
climate.rast <- spatial_sync_raster(climate.rast,master.coarse,method="bilinear")

# Write to folder
writeRaster(climate.rast,"RasterOutput/climate.grd",overwrite=TRUE)


##############################
#### Raster to data frame ####
##############################


#Read files from RasterOutput folder
RasterOutput.files <- list.files(path="RasterOutput",pattern=c(".grd"),full.names=TRUE)

#Pull them all together into a stack
raster.output <- stack(RasterOutput.files)

# Convert raster stack to data.frame
d <- as.data.frame(raster.output) # NAs in the mortality columns mean outside flightpath

# add in cell position
cells <- 1:ncell(master.coarse)
pos <- rowColFromCell(master.coarse,cells)
d$pos.x <- pos[,1]
d$pos.y <- pos[,2]

# Remove cells that are outside the boundaries of the climate layer
d <- d[!is.na(d$tmax.ann.2010),]

# Set mort to NA when distr = 0 or NA, unless mort > 0 (make sure doesn't bias)
# Remove cells that are outside the boundaries of the distribution layer

dist.col <- list()
mort.col <- list()

for(i in 1:length(species.set)) {
  spgrp.text <- sprintf("%02d",i) # add leading zero
  search <- paste("spgrp",spgrp.text,".distr",sep="")
  dist.col[[i]] <- grep(search,names(d))  
  
  
  spgrp.text <- sprintf("%02d",i) # add leading zero
  search <- paste("spgrp",spgrp.text,".mort",sep="")
  mort.col[[i]] <- grep(search,names(d))
  
  total.mort <- rowSums(d[,mort.col[[i]]],na.rm=TRUE)
  
  dist.present <- ifelse( (total.mort > 0) | (d[,dist.col[[i]] ] > 0) ,1,NA)
  
  d[, mort.col[[i]] ] <- d[, mort.col[[i]] ] * dist.present    

}

distr.cols <- grep("distr",names(d))
d$allsp.liveba <- rowSums(d[,distr.cols],na.rm=TRUE)


##Melt data frame so that all mortality variables are in a single column
mort.cols <- grep("mort.tpa",names(d))
d.melt <- melt(d,measure.vars=mort.cols,variable.name="ID",value.name="mort.tpa")
d.melt$mort.bin <- fn.bin(d.melt$mort.tpa)

#remove rows where mortality is NA (outside of flight path or species distribution)
d.melt <- d.melt[!is.na(d.melt$mort.tpa),]

#Add new columns
d.melt$year <- substr(d.melt$ID,2,5) #add column for year
d.melt$spgroup <- substr(d.melt$ID,7,13) #add column for species group


#################################################
#### Derive climate variables for data frame ####
#################################################

#Calculate climate normal and standard deviation
normal.years <- 1981:2015
tmax.names <- paste("tmax.ann.",normal.years,sep="")
ppt.names <- paste("ppt.ann.",normal.years,sep="")
tmax.cols <- grep("tmax.ann.",names(d.melt))
ppt.cols <- grep("ppt.ann",names(d.melt))

d.melt$Tnorm <- rowMeans(d.melt[,tmax.names],na.rm=TRUE) #for all years 1981-2014
d.melt$Pnorm <- rowMeans(d.melt[,ppt.names],na.rm=TRUE) #for all years 1981-2014

d.melt$Tsd <- apply(d.melt[,tmax.names],1,sd)
d.melt$Psd <- apply(d.melt[,ppt.names],1,sd)

#Calculate current and lagged temp/precip
v.year <- 2005:2015 #

for(year in v.year) {
  
  tmax.colname <- paste("tmax.ann.",year,sep="")
  ppt.colname <- paste("ppt.ann.",year,sep="")
  tmax.col <- grep(tmax.colname,names(d.melt))
  ppt.col <- grep(ppt.colname,names(d.melt))

  d.melt[d.melt$year==year,"T0"] <- d.melt[d.melt$year==year,tmax.col]
  d.melt[d.melt$year==year,"T1"] <- d.melt[d.melt$year==year,tmax.col-1]
  d.melt[d.melt$year==year,"T2"] <- d.melt[d.melt$year==year,tmax.col-2]
  d.melt[d.melt$year==year,"T3"] <- d.melt[d.melt$year==year,tmax.col-3]
  d.melt[d.melt$year==year,"T4"] <- d.melt[d.melt$year==year,tmax.col-4]
  d.melt[d.melt$year==year,"T5"] <- d.melt[d.melt$year==year,tmax.col-5]
  
  d.melt[d.melt$year==year,"P0"] <- d.melt[d.melt$year==year,ppt.col]
  d.melt[d.melt$year==year,"P1"] <- d.melt[d.melt$year==year,ppt.col-1]
  d.melt[d.melt$year==year,"P2"] <- d.melt[d.melt$year==year,ppt.col-2]
  d.melt[d.melt$year==year,"P3"] <- d.melt[d.melt$year==year,ppt.col-3]
  d.melt[d.melt$year==year,"P4"] <- d.melt[d.melt$year==year,ppt.col-4]
  d.melt[d.melt$year==year,"P5"] <- d.melt[d.melt$year==year,ppt.col-5]
}



#Calculate Temp z.score
d.melt$Tz0 <- (d.melt$T0-d.melt$Tnorm)/d.melt$Tsd
d.melt$Tz1 <- (d.melt$T1-d.melt$Tnorm)/d.melt$Tsd
d.melt$Tz2 <- (d.melt$T2-d.melt$Tnorm)/d.melt$Tsd
d.melt$Tz3 <- (d.melt$T3-d.melt$Tnorm)/d.melt$Tsd
d.melt$Tz4 <- (d.melt$T4-d.melt$Tnorm)/d.melt$Tsd
d.melt$Tz5 <- (d.melt$T5-d.melt$Tnorm)/d.melt$Tsd

#Calculate Precip z.score
d.melt$Pz0 <- (d.melt$P0-d.melt$Pnorm)/d.melt$Psd
d.melt$Pz1 <- (d.melt$P1-d.melt$Pnorm)/d.melt$Psd
d.melt$Pz2 <- (d.melt$P2-d.melt$Pnorm)/d.melt$Psd
d.melt$Pz3 <- (d.melt$P3-d.melt$Pnorm)/d.melt$Psd
d.melt$Pz4 <- (d.melt$P4-d.melt$Pnorm)/d.melt$Psd
d.melt$Pz5 <- (d.melt$P5-d.melt$Pnorm)/d.melt$Psd

#Remove raw climate data from dataframe
d.melt <- d.melt[,-c(tmax.cols,ppt.cols)]


## extract BA of the focal sp.
get.colnums <- function(x,colnames) {
  grep(x,colnames)
}
distr.colname <- paste(d.melt$spgroup,".distr",sep="")
distr.col.numbers <- sapply(distr.colname,get.colnums,names(d.melt))
d.melt$spgrp.focal.liveba <- as.numeric(d.melt[cbind((1:nrow(d.melt)),distr.col.numbers)])

# remove raw species live ba columns
distr.cols <- grep("distr",names(d.melt))
d.melt <- d.melt[,-distr.cols]  #<------THIS TABLE IS READY FOR ANALYSIS

# #Check tables
# names(d.melt)
# head(d.melt)

#Write table to external file
write.table(d.melt,file="df_export_18spp.txt",sep=",",col.names=TRUE)


# if starting with data analysis, start here
d.melt <- read.table("df_export.txt",sep=",",header=TRUE)


##########Plot data#############
#plot boxplot by species group and year
d.mort <- d.melt[d.melt$mort.tpa > 0,]
boxplot(d.mort$mort.tpa~d.mort$spgroup+d.mort$year,na.rm=TRUE)


###########################
####Gamma hurdle model#####
###########################

d.melt$pixelID <- paste0(d.melt$pos.x,"_",d.melt$pos.y)

#New model
m1 <- glm(mort.bin ~ Tnorm*Pnorm+spgrp.focal.liveba+allsp.liveba, data=d.melt[d.melt$spgroup=="spgrp04",], family=binomial(link = logit))
m2 <- glm(mort.bin ~ Tz0*Pz0+spgrp.focal.liveba+allsp.liveba, data=d.melt[d.melt$spgroup=="spgrp04",], family=binomial(link = logit))


AIC(m1,m2)
summary(m1)
summary(m2)

#New model
m1 <- glm(mort.bin ~ Tnorm + spgroup, data=d.melt, family=binomial(link = logit))
m2 <- glm(mort.bin ~ Tnorm + Pnorm, data=d.melt, family=binomial(link = logit))
m3 <- glm(mort.bin ~ Tnorm*Pnorm, data=d.melt, family=binomial(link = logit))
m4 <- glm(mort.bin ~ Tnorm*Pnorm + spgroup, data=d.melt, family=binomial(link = logit))
m5 <- glm(mort.bin ~ Tz0 + spgroup, data=d.melt, family=binomial(link = logit))
m6 <- glm(mort.bin ~ Tz0 + Pz0 + spgroup, data=d.melt, family=binomial(link = logit))
m7 <- glm(mort.bin ~ Tz0*Pz0 + spgroup, data=d.melt, family=binomial(link = logit))
m8 <- glm(mort.bin ~ Tz0*Pz0 + Tnorm*Pnorm + spgroup, data=d.melt, family=binomial(link = logit))
m9 <- glm(mort.bin ~ Tz0*Pz0 + Tnorm + Pnorm + spgroup, data=d.melt, family=binomial(link = logit))
m10 <- glm(mort.bin ~ Pz0 + I(Tz0*Pz0) + Tnorm + Pnorm + spgroup, data=d.melt, family=binomial(link = logit))
m11 <- glm(mort.bin ~ Pz0 + I(Tz0*Pz0) + Pz1*Tz1 + Tnorm + Pnorm + spgroup, data=d.melt, family=binomial(link = logit))
m12 <- glm(mort.bin ~ Pz0 + I(Tz0*Pz0) + Pz1 + I(Pz1*Tz1) + Tnorm + Pnorm + spgroup, data=d.melt, family=binomial(link = logit))
m13 <- glm(mort.bin ~ I(Tz0*Pz0) + Pz1 + I(Pz1*Tz1) + Tnorm + Pnorm + spgroup, data=d.melt, family=binomial(link = logit))
m14 <- glm(mort.bin ~ I(Tz0*Pz0) + Pz1 + I(Pz1*Tz1) + Pz2*Tz2 + Tnorm + Pnorm + spgroup, data=d.melt, family=binomial(link = logit))
m15 <- glm(mort.bin ~ Tz0*Pz0 + Pz1*Tz1 + Pz2*Tz2 + Tnorm*Pnorm + spgroup, data=d.melt, family=binomial(link = logit))

BIC(m1,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,m14,m15)
summary(m15)


m1b <- glm(mort.tpa ~ Tnorm*Pnorm + spgroup, data=d.melt[d.melt$mort.tpa > 0,], family=Gamma(link = log))
m2b <- glm(mort.tpa ~ Tnorm*Pnorm + Tz0*Pz0 + spgroup, data=d.melt[d.melt$mort.tpa > 0,], family=Gamma(link = log))
m3b <- glm(mort.tpa ~ Tnorm*Pnorm + Tz0*Pz0 + Tz1 + Pz1 + spgroup, data=d.melt[d.melt$mort.tpa > 0,], family=Gamma(link = log))
m4b <- glm(mort.tpa ~ Tz1*Pz1 + spgroup, data=d.melt[d.melt$mort.tpa > 0,], family=Gamma(link = log))
m5b <- glm(mort.tpa ~ Pz1 + I(Tz1*Pz1) + spgroup, data=d.melt[d.melt$mort.tpa > 0,], family=Gamma(link = log))
m6b <- glm(mort.tpa ~ Pz1 + I(Tz1*Pz1), data=d.melt[d.melt$mort.tpa > 0,], family=Gamma(link = log))
m7b <- glm(mort.tpa ~ Pz1*Tz1, data=d.melt[d.melt$mort.tpa > 0,], family=Gamma(link = log))

m8b <- glmer(mort.tpa ~ Pz1*Tz1 + (Pz1|spgroup), data=d.melt[d.melt$mort.tpa > 0,], family=Gamma(link = log))

m1c <- glm(mort.tpa ~ Pz1+Tz1+Pz2+Tz2, data=d.melt[d.melt$mort.tpa > 0 & 
            d.melt$spgroup=="spgrp01",], family=Gamma(link = log))

m2c <- glm(mort.tpa ~ Pz1+Tz1+Pz2+Tz2, data=d.melt[d.melt$mort.tpa > 0 & 
           d.melt$spgroup=="spgrp02",], family=Gamma(link = log))
summary(m1c)
summary(m2c)
summary(m8b)
ranef(m8b)

BIC(m1b,m2b,m3b,m4b,m5b,m6b,m7b)
summary(m7b)



#Logit model for estimating probability that pixel will contain mortality given ppt and tmax.aug
m1 <- glm(mort.bin ~ ppt*tmax, data = df, family = binomial(link = logit))
m2 <- glm(mort.bin ~ ppt, data = df[complete.cases(df$ppt) & complete.cases(df$tmax),], family = binomial(link = logit))
m3 <- glm(mort.bin ~ tmax, data = df[complete.cases(df$ppt) & complete.cases(df$tmax),], family = binomial(link = logit))
m4 <- glm(mort.bin ~ tmax + I(ppt*tmax), data = df[complete.cases(df$ppt) & complete.cases(df$tmax),], family = binomial(link = logit))
m5 <- glm(mort.bin ~ ppt + I(ppt*tmax), data = df[complete.cases(df$ppt) & complete.cases(df$tmax),], family = binomial(link = logit))

#compare models using AIC
AIC(m1,m2,m3,m4,m5)

#summarize best model
summary(m4)

#Gamma model for predicting TPA dead 
m6 <- glm(mort.tpa ~ ppt, data=subset(df, mort.bin == 1), family = Gamma(link = log))
m7 <- glm(mort.tpa ~ tmax, data=subset(df, mort.bin == 1), family = Gamma(link = log))
m8 <- glm(mort.tpa ~ ppt*tmax, data = subset(df, mort.bin == 1), family = Gamma(link = log))
m9 <- glm(mort.tpa ~ ppt + I(ppt*tmax), data = subset(df, mort.bin == 1), family = Gamma(link = log))
m10 <- glm(mort.tpa ~ tmax + I(ppt*tmax), data = subset(df, mort.bin == 1), family = Gamma(link = log))

#compare models using AIC
AIC(m6,m7,m8,m9,m10)

#summarize best model
summary(m8)

#run logit mortality probability simulation
mortprob.sim <- summary(m1)$coeff[1]+summary(m1)$coeff[2]*d.melt$Tnorm+summary(m1)$coeff[3]*d.melt$Pnorm+summary(m1)$coeff[4]*d.melt$Tnorm*d.melt$Pnorm
d.melt$mortprob.sim <- exp(mortprob.sim)/(1+exp(mortprob.sim))

#plot tpa mortality model
morttpa.sim <- summary(m8)$coeff[1]+summary(m8)$coeff[2]*df$ppt+summary(m8)$coeff[3]*df$tmax+summary(m8)$coeff[4]*df$ppt*df$tmax
df$morttpa.sim <- exp(morttpa.sim)

#Calculate conditional probabilities of tpa mortality
df$morttpa.condprob <- df$morttpa.sim*df$mortprob.sim

#Plot simulations
ggplot(d.melt) + geom_point(aes(Tnorm,Pnorm,colour=mortprob.sim)) + scale_colour_gradientn(colours=rainbow(3)) + theme_bw() #mortality probability simulation
ggplot(df) + geom_point(aes(tmax,ppt,colour=morttpa.sim)) + scale_colour_gradientn(colours=rainbow(3)) + theme_bw() #mortality tpa simulation
ggplot(df) + geom_point(aes(tmax,ppt,colour=morttpa.condprob)) + scale_colour_gradientn(colours=rainbow(3)) + theme_bw() #conditional mortality tpa simulation

#bootstrapped confidence intervals for gamma hurdle model (run function below first)
#not currently working
library(boot)
b <- boot(df, hurdle_fn, R = 500,parallel="multicore")
b.ci <- boot.ci(b, type = "bca",parallel="multicore")
print(b.ci)


################################
#One-time use functions and code
################################

#bootstrap function
hurdle_fn <- function(data, i) {
  dat_boot <- data[i, ]
  m1 <- glm(mort.bin ~ ppt*tmax, data = dat_boot,
            family = binomial(link = logit))
  m2 <- glm(mort.tpa~ppt*tmax, data = subset(dat_boot, mort.bin == 1),
            family = Gamma(link = log))
  bin_coef <- plogis(coef(m1)[[1]])
  gamma_coef <- exp(coef(m2)[[1]])
  exp(log(bin_coef) + log(gamma_coef))
}

#############
# species distribution file sync loop
#############

file <- list.files(path="Survey GDBs",pattern=as.character(2015),full.names=TRUE)
layers <- ogrListLayers(file)
mortality.index <- grep("ADS",layers,ignore.case=TRUE)
mort.polygon <- readOGR(dsn=file,layer=layers[mortality.index])
focal.sp <- unique(c(mort.polygon$HOST1,mort.polygon$HOST2,mort.polygon$HOST3))
focal.sp <- unique(c(focal.sp,746,41,312,17)) # add grand fir, bigleaf maple, port-orford cedar, quaking aspen (present in 2014 but not 2015)
focal.sp <- focal.sp[focal.sp>=10 & focal.sp < 998]

# eliminate genus-level
genus <- c(10,800,100,510) #fir, oak, pine, eucalyptus
focal.sp <- focal.sp[!(focal.sp %in% genus)]

#read all files that correspond to species in aerial survey data
file <- list.files(path="LiveBasalAreaRasters",full.names=TRUE)
file <- file[as.numeric(substr(file,23,as.numeric(nchar(file[1]))-3)) %in% focal.sp]


for(i in 1:length(file)){
  s <- raster(file[i])
  
  #compute aggregation factor
  agg.factor.sp <- round(res.coarse/(res(s)[1]))
  
  #aggregate raster
  s <- aggregate(s,fact=agg.factor.sp,fun=mean)
  
  #reproject and align grids
  s <- spatial_sync_raster(s,master.coarse,method="bilinear")
  
  #write to file
  sp.name <- substr(file[i],22,nchar(file[i])-4) #22
  writeRaster(s,paste("Raster Live Basal Area Synced/",sp.name,sep=""),overwrite=TRUE) 
}











#Only run this section of script if the raster hasn't been generated for a species
#Open raster files
s15 <- raster("Raster Live Basal Area/s15.tif")
s81 <- raster("Raster Live Basal Area/s81.tif")
s116 <- raster("Raster Live Basal Area/s116.tif")
s122 <- raster("Raster Live Basal Area/s122.tif")
s202 <- raster("Raster Live Basal Area/s202.tif")


#! Need to aggregate first; and after that, syncing will result in some "bleeding"--saying it is present some places where it is not.

#Spatially sync live basal area raster to master.coarse raster
s15 <- spatial_sync_raster(s15,master.coarse,method="ngb")
s81 <- spatial_sync_raster(s81,master.coarse,method="ngb")
s116 <- spatial_sync_raster(s116,master.coarse,method="ngb")
s122 <- spatial_sync_raster(s122,master.coarse,method="ngb") 
s202 <- spatial_sync_raster(s202,master.coarse,method="ngb")

#write synced species distributions to raster files
writeRaster(s15,"Raster Live Basal Area Synced/s15",format="GTiff",overwrite=TRUE) 
writeRaster(s81,"Raster Live Basal Area Synced/s81",format="GTiff",overwrite=TRUE)
writeRaster(s116,"Raster Live Basal Area Synced/s116",format="GTiff",overwrite=TRUE)
writeRaster(s122,"Raster Live Basal Area Synced/s122",format="GTiff",overwrite=TRUE)
writeRaster(s202,"Raster Live Basal Area Synced/s202",format="GTiff",overwrite=TRUE)

