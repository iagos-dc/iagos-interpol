args = commandArgs(trailingOnly=TRUE)

lon1_region=as.numeric(args[1]) ## Region western edge.
lon2_region=as.numeric(args[2]) ## Region eastern edge.
lat1_region=as.numeric(args[3]) ## Region southern edge.
lat2_region=as.numeric(args[4]) ## Region northern edge.
lon1=as.numeric(args[5]) ## First gridcell longitude in the whole domain.
lon2=as.numeric(args[6]) ## Last gridcell longitude in the whole domain.
lat1=as.numeric(args[7]) ## First gridcell latitude in the whole domain.
lat2=as.numeric(args[8]) ## Last gridcell latitude in the whole domain.
lon_resol=as.numeric(args[9]) ## Longitude resolution in degrees.
lat_resol=as.numeric(args[10]) ## Latitude resolution in degrees.
model=args[11]
outfile=args[12]

## The idea is to find new edges for the current region,
## that have to fit the model grid.
## I could have chosen to widen the regions, but it would cause a superimposition of the regions that share an edge.
## The best choice would be the closest gridcell edges.

## We must be careful with the gridcell definition. INCA locates its gridcells with their center.
## Thus, by assuming the gridcells to be located like MOCAGE, we make an error shifting the grid by lon_resol/2 and lat_resol/2.
## In that case, we must shift the grid back.
if (model %in% c('INCA', 'EMAC', 'Oslo-CTM3', 'GEOS-Chem', 'MOZART3', 'UKESM1.1')){
    shift_grid_lon=-lon_resol/2
    shift_grid_lat=-lat_resol/2
}else if (model %in% c('MOCAGE')){
    shift_grid_lon=0
    shift_grid_lat=0
}

## Initialization.
x=lon1+shift_grid_lon ; ilon=1
y=lat1+shift_grid_lat ; ilat=1

## 1/ Region western limit.
while (x < lon1_region-lon_resol){ ## The loop stops when the region edge is inside the current gridcell.
    x=x+lon_resol
    ilon=ilon+1
}
x1_west=x ## Western closest gridcell edge.
x1_east=x+lon_resol ## Eastern closest gridcell edge.

## Now we chose the closest one between the eastern and the western size of the regional edge.
if (abs(x1_east-lon1_region) <= abs(x1_west-lon1_region)){
    x1=x1_east
    ilon1=ilon+1
}else{
    x1=x1_west
    ilon1=ilon
}
## 2/ Region eastern limit.
while (x < lon2_region-lon_resol){
    x=x+lon_resol
    ilon=ilon+1
}
x2_west=x ## Western closest gridcell edge.
x2_east=x+lon_resol ## Eastern closest gridcell edge.

## Now we chose the closest one between the eastern and the western size of the regional edge.
## It is not exactly the same than for the region western limit, since the gridcells are defined with their western edge.
if (abs(x2_east-lon2_region) <= abs(x2_west-lon2_region)){
    x2=x2_east ## Physical limits. Remain unchanged.
    ilon2=ilon ## Last gridcell coordinate. Must be adjusted compared to ilon1.
}else{
    x2=x2_west
    ilon2=ilon-1
}

## 3/ Region southern limit.
while (y < lat1_region-lat_resol){ ## The loop stops when the region edge is inside the current gridcell.
    y=y+lat_resol
    ilat=ilat+1
}
y1_south=y ## Southern closest gridcell edge.
y1_north=y+lat_resol ## Northern closest gridcell edge.

## Now we chose the closest one between the northern and the southern size of the regional edge.
if (abs(y1_north-lat1_region) <= abs(y1_south-lat1_region)){
    y1=y1_north
    ilat1=ilat+1
}else{
    y1=y1_south
    ilat1=ilat
}
## 4/ Region northern limit.
while (y < lat2_region-lat_resol){
    y=y+lat_resol
    ilat=ilat+1
}
y2_south=y ## Southern closest gridcell edge.
y2_north=y+lat_resol ## Northern closest gridcell edge.

## Now we chose the closest one between the northern and the southern size of the regional edge.
## It is not exactly the same than for the region southern limit, since the gridcells are defined with their southern edge.
if (abs(y2_north-lat2_region) <= abs(y2_south-lat2_region)){
    y2=y2_north ## Physical limits. Remain unchanged.
    ilat2=ilat ## Last gridcell coordinate. Must be adjusted compared to ilat1.
}else{
    y2=y2_south
    ilat2=ilat-1
}

mytext=paste(ilon1,ilon2,ilat1,ilat2,x1,x2,y1,y2)
con <- file(outfile,open="w")
writeLines(mytext, con = con)
close(con)
