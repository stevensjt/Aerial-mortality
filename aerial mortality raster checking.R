# load an example raster to get res, extent, etc


df.to.raster <- function(d,raster.template,spgroup,year,data.colname)  {
  
  #wipe the raster clean
  newraster <- raster.template * NA
  
  #look up cell number based on x and y position (rows and cols)
  d$cellnum <- cellFromRowCol(mortraster.coarse,d$pos.x,d$pos.y)
  
  # thin to desired spgroup and year
  d.thin <- d[d$spgroup == spgroup,]
  d.thin <- d.thin[d.thin$year == year,]
  
  # extract values to raster
  newraster[d.thin$cellnum] <- d.thin[,data.colname]
  
  return(newraster)
  
}


mortraster.coarse <- raster("RasterOutput/Y2015_spgrp1.grd")
d <- read.table("df_export_18spp.txt",header=TRUE,sep=",")


## change to dead tpa
a <- df.to.raster(d,mortraster.coarse,"spgrp01",2011,"P0")

plot(a)




## change to dead tpa
writeRaster(a,"Raster checking/precip_pipo_2011_0.tif")


##open and write climate rasters etc
c <- stack("Climate layers/PRISM_ppt_wateryear/PRISM_ppt_wateryear.grd")
c <- c[["ppt.ann.2011"]]
writeRaster(c,"Raster checking/ppt_2011_orig.tif",overwrite=TRUE)

distr <- stack("RasterOutput/spgrp10_distr.grd")
writeRaster(distr,"Raster checking/spgrp10_distr_orig.tif")
















# load the data frame to rasterize








newraster2 <- setValues(newraster,values=d.thin$mort.bin,index=d.thin$cellnum)


newraster[]








# Extract data needed to re-create the same grid
master.coarse.res <- res(master.fine)
master.coarse.extent <- extent(master.fine)[c(1,3,2,4)]



master.fine.res
agg.factor
master.fine.extent

master.fine.coarse <- master.fine.res/agg.factor



