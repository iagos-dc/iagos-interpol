args = commandArgs(trailingOnly=TRUE)

nlon     =as.integer(args[1])
nlat     =as.integer(args[2])
lon_dom_1=as.numeric(args[3])
lon_dom_2=as.numeric(args[4])
lat_dom_1=as.numeric(args[5])
lat_dom_2=as.numeric(args[6])
outfile  =args[7]

print(paste(nlon, nlat, lon_dom_1, lon_dom_2, lat_dom_1, lat_dom_2))
## Deriving the horizontal resolution.
lon_resol=(lon_dom_2-lon_dom_1)/nlon
lat_resol=round((lat_dom_2-lat_dom_1)/nlat,digits=6)
## Deriving the domain boundaries, in terms of gridcell coordinates.
lon1=lon_dom_1+lon_resol/2
lon2=lon_dom_2-lon_resol/2
lat1=lat_dom_1+lat_resol/2
lat2=lat_dom_2-lat_resol/2
mytext=paste(lon_resol,lat_resol,lon1,lon2,lat1,lat2)

con <- file(outfile,open="w")
writeLines(mytext, con = con)
close(con)
