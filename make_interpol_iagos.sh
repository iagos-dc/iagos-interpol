#!/bin/bash


## First of all: please indicate here the paths for the NetCDF libraries.
NETCDF_INCDIR="-I$NETCDF_INCDIR -I$NETCDFFORTRAN_INCDIR -I$FORTRAN_INTEL_INCDIR"
NETCDF_LIBDIR="-L$NETCDF_LIBDIR -L$NETCDFFORTRAN_LIBDIR -L$FORTRAN_INTEL_LIBDIR"
NETCDF_LIB="-lnetcdf -lnetcdff -lstdc++"


########## 1/ Main commands ##########
### This first bin is useful if one wants to derive some new model variables (e.g. O3/CO ratio).
### It also extracts pressure and PV, only if it has not been done yet with extract_model_variables.sh.
bin_extract_last_variables=TRUE
### When all the model variables are extracted, the matter is, if needed, to interpolate them
### onto another grid. It is useful in case of a multimodel assessment: first we put all the models
### on a same grid. Then, the projection of the IAGOS data can be done on one single grid.
bin_change_grid=TRUE
### The first main part of this program interpolates the NetCDF IAGOS data onto a chosen model grid,
### and computes monthly means for several measured variables (then in NetCFD format).
### It reads each IAGOS file and makes a reverse linear interpolation of each observation onto the grid.
### Next, it computes monthly means on each grid cell, using spatial weighting factors based on distances.
### It yields one output file for each month and each layer, both for observations and model.
### Almost everything is made by the Fortran routine interpol_iagos.f90 (and associated subroutines).
### This part is activated when the following logical bin_interpol_iagos is TRUE.
### Duration: - a few hours for the monthly-resolved files, with two layers and two chemical species, on 20 years.
###           - several days for the daily-resolved files (same conditions).
bin_interpol_iagos=TRUE
### The second main part of the script generates regional time series.
### First, it extracts the data representing the selected regions only.
### Then it computes monthly time series (for each season) by averaging the gridded regional fields.
### It yields one output file for each region and each layer, both for observations and model.
### The extraction relies on an R routine that adapts the regions coordinates to the grid resolution,
### and on a few nco commands. So you need to have the R and nco modules loaded.
### The time series calculation is made by the Fortran routine regional_average_iagos.f90 (and associated subroutines).
### This part is activated when the following logical bin_regional_time_series is TRUE.
### Duration: approximately 3h with two layers and two chemical species, on 20 years.
bin_regional_time_series=TRUE
### The third main part of the script generates seasonal and yearly climatologies on each vertical grid level.
### In order to keep an equilibrated sampling between the four seasons, it adjusts
### the first and the last dates set by the user (a little bit later in the file) by filtering out the first and last years if they are not fully sampled (at least one month per season), fixing winter as the first season
### of the year and fall as the last season.
### The yearly means are computed by averaging the seasonal means together.
### It yields one output file for each season and also for the whole year, both for observations and model.
### Only nco commands are used here.
### This part is activated when the following logical bin_seasonal_climatology is TRUE.
bin_seasonal_climatology=TRUE
### Averaging vertically, following a given filtering process.
bin_average_layers=TRUE
### Averaging vertically and zonally over predefined longitude bands, following a given filtering process.
bin_zonal_cross_sections=TRUE

### Options: bin_mask_NOy if you want to process the individual NOy species output once again, this time with the IAGOS NOy mask.
bin_mask_NOy=TRUE
bin_recompile=TRUE ## TRUE if the user wants to compile the fortran routines.
bin_colissimo=TRUE ## Copies the final output files into an easy-to-find repertory.

### Start and end of the comparison period.
### Please note that the IAGOS measurements start in Aug. 1994.
yyyymm1=199408 ; yyyymm2=202212
### One very practical option if you run this code on your local machine, before leaving your office,
### especially before starting your week-end. Set 'yes' to the following binary
### if you want your computer to shut down when the work is done.
shutdown_when_run_completed=no
### Next: tell me where the Pack_code directory is located. Its path will be called base_dir.
base_dir=../../../
### Also: tell me where all the data are stored. This has to be the main data repertory.
### (in my case - see below, another Pack_code folder exists in a storing space)
data_dir=$CCCSCRATCHDIR/Pack_code/data


### Name of the model, the configuration and the name of the simulation to be compared to the observations.

model=your_model_name
config=your_config_name
experiment=your_run_name
sci_program=your_scientific_program ## Helpful if you have several models involved in a same experiment. If you don't have one, you can answer "idk".

### Define here the configuration of your model output as input for this script.
nlon_init=144
nlat_init=143
nlev_init=39

### Define here the coordinate system in your model output. If at least one of the following three binaries is TRUE,
### then the adjust_model_output_files.sh script is called to readapt the coordinate system.
lon_axis_0_360=TRUE_or_FALSE ## If TRUE, then the axis origin will be shifted to the 0 meridian.
lat_axis_90N_90S=TRUE_or_FALSE ## If TRUE (i.e. if latitude is decreasing), then the axis will be reversed.
P_axis_bottom_to_top=TRUE_or_FALSE ## If TRUE (i.e. if pressure is decreasing), then the axis will be reversed.

### Define here the vertical extension of the layers, for the vertical boundaries. No impact if the bin_change_grid binary is activated.
## First: the vertical domain for the interpol_iagos.f90 routine.
upper_level_UTLS=21 ; lower_level_UTLS=27
## Second: the vertical domain for the climatology calculations.
upper_level_UTLS_climato=21 ; lower_level_UTLS_climato=27
## Third: the vertical domain for the tropical zonal cross sections.
upper_level_UTLS_zonal=21 ; lower_level_UTLS_zonal=25
## Fourth: the vertical domain for the free troposphere, if one wants to include the ascent/descent phases (LTO).
upper_level_FT=21 ; lower_level_FT=$nlev_init

### Pick the time resolution you need to separate the upper troposphere and the lower stratosphere, 
### based on the potential vorticity (PV) fields. In other words, it corresponds to the PV time resolution.
### The monthly resolution is adapted to the CCMI experiments because the model outputs are stored with a monthly resolution. You can guess it is less voluminous.
### With the daily resolution, you need daily fields concatenated into monthly files.
time_resolution='daily'
### If the PV is not available in the model files, then set extract_PV to 'no'.
extract_PV='yes'
## If the pressure at the full levels are not available, then set extract_P to 'no'. 
## The program will interpolate with regard to hybrid pressures instead.
extract_P='yes'

### Tropopause dynamical definition, for extratropics.
## Choose the value for the tropopause potential vorticity (PV).
PV_tropopause=2 # unit: potential vorticity units (PVU)
### Pressure constant P_0, in case we have to calculate the 3D pressure from the sigma-hybrid coefficients,
### using the equation: P_i = P_0*A_i + P_S*B_i
P_0=100000 ## In Pa
### In case we modify the model's grid, this integer tells which version of the VERTICAL_INTERP.f90 routine we use.
### Version 0 is straightforward if the final grid is much more resolved than the initial grid.
version_vertical_interp=0

#############################################
########### User's shopping list. ###########
#############################################

### List of regions where you plan to compute mean values, for the regional analyse.
### Their acronyms are defined like in Cohen et al. (2018). Their coordinates are
### defined a little bit later.
list_regions="WNAm EUS NAt Eur WMed MidE Sib NEAs"
### List of layers where you plan to compute mean values. 
### IMPORTANT: Make sure they exist in the regional_average_iagos.f90 routine.
## Their characteristics are defined a little bit later too.
## UT for upper troposphere
## LS for lower stratosphere
## LFT for lower-free-troposphere, typically 3-6 km
## MFT for middle-free-troposphere, typically 6-9 km
## UBL for upper boundary layer
## all_layers for computing mean values without dependence on the layer.
## UTLS is quite the same, but with the UT as ground level.
list_layers="UT LS UTLS"
list_layers_climatologies="UT LS UTLS"
list_layers_zonal_cross_sections="UT" ## Not really useful to have the zonal cross sections in the LS, since it is not sampled in the tropics.
list_seasons_def_tropics=Boreal' 'Tropical_Atl' 'Tropical_Afr' 'South_Asia
# list_summers_def=individual
id_job=${experiment}_${config}_${model}_${yyyymm1}_${yyyymm2}
function define_seasons(){
    case ${1} in
	'Boreal' ) list_seasons=DJF' 'MAM' 'JJA' 'SON
		   seasons_start_indexes='03_06_09_12' ;; ## Needs the 1st month of every season, whatever the order.
	'Tropical_Atl' ) list_seasons=DJF' 'MAMJ' 'JA' 'SON 
	    seasons_start_indexes='03_07_09_12' ;;
	'Tropical_Afr' ) list_seasons=DJFM' 'AM' 'JJASO' 'N
	    seasons_start_indexes='04_06_11_12' ;;
	'South_Asia' ) list_seasons=DJF' 'MAM' 'JJAS' 'ON
	    seasons_start_indexes='03_06_10_12' ;;
	'individual' ) list_seasons=Dec' 'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov
	    seasons_start_indexes='01_02_03_04_05_06_07_08_09_10_11_12' ;; ## Remains on the boreal model, but has no impact.
    esac
    echo $list_seasons
}
function define_zonal_bands(){
    zonal_region=${1}
    tmpdir=${2}
    model=${3}
    case $zonal_region in
	'Tropical_Atl' )
	    lon01=-60 ; lon02=-15 ;;
	'Tropical_Afr' )
	    lon01=-5 ; lon02=30 ;; ## As in Sauvage et al., 2007, GRL
	'South_Asia' )
	    lon01=60 ; lon02=90 ;;
    esac
## Computes lon and lat resolutions, and the coordinates for the 1st and last gridcells in each dimension.
    txtfile_resol=$tmpdir/resol_${config}.txt
    Rscript $pdir/compute_resolution.R $nlon $nlat $lon_dom_1 $lon_dom_2 $lat_dom_1 $lat_dom_2 $txtfile_resol
    lon_resol=$(awk '{printf "%8.6f",$1}' ${txtfile_resol})
    lat_resol=$(awk '{printf "%8.6f",$2}' ${txtfile_resol})
    lon1=$(awk '{printf "%9.6f",$3}' ${txtfile_resol})
    lon2=$(awk '{printf "%9.6f",$4}' ${txtfile_resol})
    lat1=$(awk '{printf "%9.6f",$5}' ${txtfile_resol})
    lat2=$(awk '{printf "%9.6f",$6}' ${txtfile_resol})
    ## Now adapts the longitude bands coordinates to the grid resolution, and finds out the corresponding grid indexes.
    tmpfile_indexes=$tmpdir/indexes.txt
    Rscript $pdir/compute_grid_indexes.R $lon01 $lon02 $lat_dom_1 $lat_dom_2 $lon1 $lon2 $lat1 $lat2 $lon_resol $lat_resol $model $tmpfile_indexes
    ### First set of outputs from the R routine: the indexes locating the new regions coordinates in the model grid.
    ilon01=$(awk '{print $1 }' ${tmpfile_indexes})
    ilon02=$(awk '{print $2 }' ${tmpfile_indexes})
    klat01=$(awk '{print $3 }' ${tmpfile_indexes})
    klat02=$(awk '{print $4 }' ${tmpfile_indexes})
    ### Second set of outputs from the R routine: the corresponding coordinates in degrees East and North.
    lon01_grid=$(awk '{printf "%8.5f",$5 }' ${tmpfile_indexes})
    lon02_grid=$(awk '{printf "%8.5f",$6 }' ${tmpfile_indexes})
    lat01_grid=$(awk '{printf "%8.5f",$7 }' ${tmpfile_indexes})
    lat02_grid=$(awk '{printf "%8.5f",$8 }' ${tmpfile_indexes})
    echo $zonal_region $ilon01 $ilon02 $klat01 $klat02 $lon01_grid $lon02_grid $lat01_grid $lat02_grid >> $tmpdir/zonal_regions_def_${config}.txt
    echo $ilon01 $ilon02 $klat01 $klat02
}
function west_bound(){
    zonal_region=${1}
    ilon1=$(awk '{print $1}' `define_zonal_bands $zonal_region $tmpdir`)
    echo $ilon1
}
function east_bound(){
    zonal_region=${1}
    ilon2=$(awk '{print $2}' `define_zonal_bands $zonal_region $tmpdir`)
    echo $ilon2
}
##### Complete list of observed variables. #####
## Here is the list of the observed variables you can ask to the program.
list_observed_variables="Ozone Carbon_monoxide Water_vapour Relative_humidity_liquid \
Relative_humidity_ice Cloud_liquid_water_particles Nitrogen_oxides Nitrogen_monoxide \
Nitrogen_dioxide Nitrogen_reactive_species Pressure Temperature"
################################################
##### Complete list of derived variables. ######
## See the following list:
list_derived_variables="Potential_temperature Ozone_to_CO_ratio"
## Caution: for each derived variable in the following list, you need to have its ingredients
## in the previous list.
################################################
##### Complete list of modelled variables. #####
## Here is the list of the modelled variables you can ask to the program.
list_model_variables="Ozone Carbon_monoxide Water_vapour Nitrogen_oxides Nitrogen_monoxide Nitrogen_dioxide Nitrogen_reactive_species Temperature Potential_temperature Ozone_to_CO_ratio Stratospheric_ozone_ratio Inert_strat_ozone Inert_strat_ozone_ratio Peroxyacetyl_nitrate Nitric_acid Stratospheric_ozone Production_rate_ozone Loss_rate_ozone Net_production_ozone NOy_ACACIA"
################################################
## List of packages contributing to the analysed IAGOS database.
## IMPORTANT: the double IAGOS-CARIBIC is important if you want to include the CARIBIC measurements.
## The first one stands for the CARIBIC measurements before it was merged with MOZAIC, and the second one after.
## Do not worry if they have the same name, it is the case in the IAGOS database too. Probably a mistake, but anyway.
list_packages_IAGOS="IAGOS-MOZAIC IAGOS-CORE IAGOS-CARIBIC IAGOS-CARIBIC"
## List of dimensions to be read from the IAGOS files. You do not have to change it.
listcoord="Time_UTC Latitude Longitude Baro_altitude Pressure" ## IAGOS files, version CNRS 8.0
## List of dimensions to be read from the simulation output files. You do not have to change it.
list_model_dimensions="Time Longitude Latitude Vertical_grid_level" # Do not touch!
#############################################
########### End of shopping list. ###########
#############################################

###### Initial definition for the regions. 
## These coordinates will be automatically readjusted to the model resolution.
## Fortunately, you do not have to care about it.
lon01_WNAm=-125    ; lon02_WNAm=-105     ; lat01_WNAm=40     ; lat02_WNAm=60  # Western North America
lon01_EUS=-90      ; lon02_EUS=-60       ; lat01_EUS=35      ; lat02_EUS=50   # Eastern North America
lon01_NAt=-50      ; lon02_NAt=-20       ; lat01_NAt=50      ; lat02_NAt=60   # Northern Atlantic
lon01_Eur=-15      ; lon02_Eur=15        ; lat01_Eur=45      ; lat02_Eur=55   # Europe
lon01_WMed=-5      ; lon02_WMed=15       ; lat01_WMed=35     ; lat02_WMed=45  # Western Mediterranean basin
lon01_MidE=25      ; lon02_MidE=55       ; lat01_MidE=35     ; lat02_MidE=45  # Proche-Moyen Orient
lon01_Sib=40       ; lon02_Sib=120       ; lat01_Sib=50      ; lat02_Sib=65   # Siberia
lon01_NEAs=105     ; lon02_NEAs=145      ; lat01_NEAs=30     ; lat02_NEAs=50  # Northeastern Asia

###### Vertical axis: sets the different layers properties.
## The most important is: the minimum PV value for the LS, and its max value for the UT.
## The pressure specification is set only to avoid sampling the free troposphere too deeply.
## Please note that the PV boundaries are for absolute PV values,
## because PV is negative in South Hemisphere.
## Also note that pressures are in Pa, not hPa.

## UTLS
lower_PV_UTLS=-1                 ; upper_PV_UTLS=100 # -1 just to be sure that the prefilled 0 value is included.
lower_P_UTLS=40000               ; upper_P_UTLS=0
lower_P_UTLS_climato=40000       ; upper_P_UTLS_climato=0
lower_P_UTLS_zonal=30000         ; upper_P_UTLS_zonal=0
## LS
lower_PV_LS=$((PV_tropopause+1)) ; upper_PV_LS=100 ## In PVU.
lower_P_LS=40000                 ; upper_P_LS=0 ## In Pa.
lower_P_LS_climato=40000         ; upper_P_LS_climato=0
lower_P_LS_zonal=30000           ; upper_P_LS_zonal=0
## TPL (transition layer)
lower_PV_TPL=$PV_tropopause      ; upper_PV_TPL=$((PV_tropopause+1))
lower_P_TPL=40000                ; upper_P_TPL=0
lower_P_TPL_climato=30000        ; upper_P_TPL_climato=0
## UT
lower_PV_UT=0                    ; upper_PV_UT=$PV_tropopause
lower_P_UT=40000                 ; upper_P_UT=0 ## The upper limit is controlled by the potential vorticity.
lower_P_UT_climato=40000         ; upper_P_UT_climato=0
lower_P_UT_zonal=30000           ; upper_P_UT_zonal=0
## FT (free troposphere)
lower_PV_FT=0                    ; upper_PV_FT=$PV_tropopause
lower_P_FT=70000                 ; upper_P_FT=0 ## In Pa.
## MFT (mid-free troposphere)
lower_PV_MFT=-100                ; upper_PV_MFT=100
lower_P_MFT=47000                ; upper_P_MFT=30000 ## In Pa.
## LFT (lower free troposphere)
lower_PV_LFT=-100                ; upper_PV_LFT=100
lower_P_LFT=70000                ; upper_P_LFT=47000 ## In Pa.
## UBL (upper boundary layer)
lower_PV_UBL=-100                ; upper_PV_UBL=100
lower_P_UBL=95000                ; upper_P_UBL=70000 ## In Pa.
## Without distinction between the layers.
lower_PV_all_layers=-1           ; upper_PV_all_layers=100 # -1 just to be sure that the prefilled 0 value is included.
lower_P_all_layers=1000000       ; upper_P_all_layers=0

## Chemical conditions:
### As seen before, the separation between UT and LS relies on the modelled PV field.
### Still, this one does not always catches all the tropopause features
### as stratospheric intrusions. Relying on it only would bring a non realistic
### stratospheric influence in the UT, and vice versa. Moreover, our aim consists in
### assessing the modelled chemical fields in the different layers AS DEFINED IN THE MODELS,
### not assessing the modelled tropopause features. Thus, we need to use the observations
### in order to complement the PV-based selection.
### I chose observed ozone values because I know it quite well, but it is also possible
### to use the observed potential temperature. You could also use the O3/CO ratio, but
### you would be limited to the CO measurement period, i.e. since Dec. 2001.
### Here, you can fix the threshold above which observed ozone cannot represent the UT and the TPL.
### Just make sure you do not cut off any part of the realistic ozone variability in the current layer.
ozone_max_UT=140 ## in ppb.
ozone_max_TPL=225 ## in ppb.
### Same for the minimum ozone threshold in the LS.
ozone_min_LS=60 ## in ppb.
## These other ones (not very important) complete the conditions fixed just above.
ozone_min_UT=0          ; ozone_max_LS=10000
ozone_min_TPL=0
ozone_min_MFT=0         ; ozone_max_MFT=10000
ozone_min_LFT=0         ; ozone_max_LFT=10000
ozone_min_UBL=0         ; ozone_max_UBL=10000
ozone_min_UTLS=0        ; ozone_max_UTLS=10000
ozone_min_all_layers=0  ; ozone_max_all_layers=10000
###########
## One filtering parameter for the regional time series calculation:
## N_obs_min_regional is the minimum threshold for the equivalent amount of data into one given gridcell during a month.
## The gridcells with a sampling below this one are not taken into account during the regional average calculation.
## Make sure it is strictly greater than 1, else the standard deviations computation will probably generate bugs.
N_obs_min_regional=2
## One filtering parameter for the seasonal climatologies calculation:
## N_obs_min_climato is the reference threshold for the IAGOS sampling, once each monthly grid cell has been averaged
## over the whole time period. For a climatological grid cell at a given latitude, the effective sampling threshold
## is obtained by multiplying N_obs_min_climato with the ratio cos(lat)/<cos(lat)>, in order to account for the size of the grid cells.
N_obs_min_climato=140
## The following line takes into account the fact that N_obs_min_climato must be adapted to the mean size of the gridcells. Thus, we have to give precision on the meaning of the amount chosen above: which resolution would it correspond to?
ref_grid=$config

## Instructions concerning the new grid, in case the user wants to regrid the model output.
## If not, then these values have no impact, as they are redefined later in the routine.
## 1/ Horizontal dimensions for the new grid.
nlon_ref=191 ; nlat_ref=144
## 2/ Vertical dimensions for the new grid. If it is cst_P_UTLS, then one has to define the
## new pressure boundaries and resolution here, in hPa.
top_level_hPa=140 ; down_level_hPa=440 ; delta_P_hPa=20
if [ "$bin_change_grid" = TRUE ] ; then
    ref_grid=cst_P_UTLS
fi

### If quicktest_interpol_IAGOS is turned on 1, then the script will only read the first IAGOS file indicated in each monthly list.
### Although it is a fast test, it spreads over the whole period.
quicktest_interpol_IAGOS=0

######### End of user's instructions #########
### What is coded next are automatic commands. The user does not have to edit it.

### The directory where the extrapolated IAGOS data are stored.
pack_code_dir=$base_dir/Pack_code
### The directory where the files will be ready for analytic treatments.
comparison_dir=$data_dir/comparison_IAGOS_models ; mkdir -p $comparison_dir
if [ "$sci_program" != "idk" ] ; then
    outdir_comparison_iagos_model=$comparison_dir/IAGOS_${sci_program}_models/IAGOS_$model
    mkdir -p $comparison_dir/IAGOS_${sci_program}_models
else
    outdir_comparison_iagos_model=$comparison_dir/IAGOS_$model
fi
mkdir -p $outdir_comparison_iagos_model
global_outdir=$outdir_comparison_iagos_model/global_NetCDF ; mkdir -p $global_outdir
global_outdir_config=$global_outdir/$ref_grid ; mkdir -p $global_outdir_config
global_outdir_exp=$global_outdir_config/$experiment ; mkdir -p $global_outdir_exp
global_outdir_resol=$global_outdir_exp/$time_resolution ; mkdir -p $global_outdir_resol
outdir_interpol_iagos=$global_outdir_resol/global_monthly_maps ; mkdir -p $outdir_interpol_iagos ## The gridded data sets, in monthly means. Used for both global climatologies and regional time series.
for layer in $list_layers ; do
    mkdir -p $outdir_interpol_iagos/$layer
done
### Output directory for climatologies.
outdir_seasonal_clim=$global_outdir_resol/seasonal_climatologies ; mkdir -p $outdir_seasonal_clim
outdir_level_by_level=$outdir_seasonal_clim/level_by_level ; mkdir -p $outdir_level_by_level
### Current and parent directories
# pdir=$base_dir/Pack_code/routines/grid_IAGOS_model
pdir=`pwd`
parent_dir=$base_dir/Pack_code/routines
tmpdir=$pdir/tmp/tmp_${id_job} 
mkdir -p $pdir/tmp ; mkdir -p $tmpdir
tmpdir_there=$data_dir/tmp/job_${id_job} ## Voluminous files: on the scratch dir.
mkdir -p $data_dir/tmp ; mkdir -p $tmpdir_there
exe_dir=$pdir/executables ; mkdir -p $exe_dir



### The directory where the regional output files are created.
### They will be the input files for the analysis R routine.
regional_outdir=$outdir_comparison_iagos_model/regional_NetCDF ; mkdir -p $regional_outdir
regional_outdir_config=$regional_outdir/$ref_grid ; mkdir -p $regional_outdir_config
regional_outdir_exp=$regional_outdir_config/$experiment ; mkdir -p $regional_outdir_exp
regional_outdir_resol=$regional_outdir_exp/$time_resolution ; mkdir -p $regional_outdir_resol
dir_regional_time_series=$regional_outdir_resol/time_series ; mkdir -p $dir_regional_time_series
dir_regions_config_experiment_resol=$dir_regional_time_series ## YannC: temporary, as long as I do it in parallel.

if [ "$bin_regional_time_series" = TRUE ] ; then
    echo "dir_regional_time_series = "$dir_regional_time_series
fi

###### Configuration of the geommetry for the current experiment.
## Domain definition (in degrees north and east respectively). Here, global scale.
lat_dom_1=-90 ; lat_dom_2=90
lon_dom_1=-180 ; lon_dom_2=180
horiz_resol=${nlon_init}_${nlat_init}
if  [ "$model" = MOCAGE ] || [ "$model" = EMAC ] || [ "$model" = Oslo-CTM3 ] || [ "$model" = MOZART3 ] || [ "$model" = UKESM1.1 ] ; then
    nlat_eff=$nlat_init
elif [ "$model" = INCA ] || [ "$model" = CMIP6_INCA ] ; then
    nlat_eff=$((nlat_init+1)) ## This has no incidence on the resolution. It just allows to create an output file with nlat+1 meridian cells, as in the original INCA file.
else
    echo 'UNKNOWN MODEL: STOPPING THE EXECUTION'
    exit
fi

if [ "$ref_grid" = cst_P_UTLS ] ; then
    file_cst_P_levels=$tmpdir/cst_P_levels.txt
    file_horiz_grid=$pdir/ref_horizontal_grid_${nlon_ref}x${nlat_ref}.nc
    ref_nlon_filter=144 ; ref_nlat_filter=143 ## Dimensions for the grid taken as reference for the filter (initially INCA). Note that the _ref suffix does not mean the same as in the previous variables.
    nlat_eff=$nlat_ref
    ## This grid starts at 80 hPa, and is 10 hPa resolved at the cruise altitudes.
    ## Since we make a linear interpolation for each output level from the closest input levels, 
    ## this high resolution allows to account for all the vertical grid levels from all the participating models.
    top_level_hPa=140 ; down_level_hPa=440 ; delta_P_hPa=20
    P_levels_hPa=`seq $top_level_hPa $delta_P_hPa $down_level_hPa`
    nlev_ref=0 ; P_levels_Pa=""
    for lev in $P_levels_hPa ; do nlev_ref=$((nlev_ref+1)) ; P_levels_Pa=$P_levels_Pa" "$((lev*100)); done
    lon_dom_ref_1=$lon_dom_1 ; lon_dom_ref_2=$lon_dom_2
    lat_dom_ref_1=$lat_dom_1 ; lat_dom_ref_2=$lat_dom_2
    echo $P_levels_Pa > $file_cst_P_levels
    # Initialising the vertical level indexes.
    upper_level_UTLS=1 ; lower_level_UTLS=1
    upper_level_FT=1 ; lower_level_FT=1
    upper_level_UTLS_climato=1 ; lower_level_UTLS_climato=1
    upper_level_UTLS_zonal=1 ; lower_level_UTLS_zonal=1
    # Then incrementing.
    for P_lev in $P_levels_Pa ; do
	if [ $P_lev -lt $upper_P_UTLS ] ; then
	    upper_level_UTLS=$((upper_level_UTLS+1))
	fi
	if [ $P_lev -lt $lower_P_UTLS ] ; then
	    lower_level_UTLS=$((lower_level_UTLS+1))
	fi
	if [ $P_lev -lt $upper_P_FT ] ; then
	    upper_level_FT=$((upper_level_FT+1))
	fi
	if [ $P_lev -lt $lower_P_FT ] ; then
	    lower_level_FT=$((lower_level_FT+1))
	fi
	if [ $P_lev -lt $upper_P_UTLS_climato ] ; then
	    upper_level_UTLS_climato=$((upper_level_UTLS_climato+1))
	fi
	if [ $P_lev -lt $lower_P_UTLS_climato ] ; then
	    lower_level_UTLS_climato=$((lower_level_UTLS_climato+1))
	fi
	if [ $P_lev -lt $upper_P_UTLS_zonal ] ; then
	    upper_level_UTLS_zonal=$((upper_level_UTLS_zonal+1))
	fi
	if [ $P_lev -lt $lower_P_UTLS_zonal ] ; then
	    lower_level_UTLS_zonal=$((lower_level_UTLS_zonal+1))
	fi
    done
else
    file_cst_P_levels=NONE
    ref_nlon=144 ; ref_nlat=143
    nlon_ref=144 ; nlat_ref=143
    nlev_ref=$nlev_init
    ref_nlon_filter=$nlon_ref; ref_nlat_filter=$nlat_ref
fi

## There, we adjust the threshold N_obs_min_climato with respect to the current resolution, compared to the chosen reference resolution.
ncells=$((nlon_ref*nlat_ref)) ; ncells_0=$((ref_nlon_filter*ref_nlat_filter))
N_obs_min_climato=$(echo "scale=2 ; $N_obs_min_climato * $ncells_0 / $ncells" | bc)
echo "ncells, ncells_0, N_obs_min_climato=" $ncells $ncells_0 $N_obs_min_climato
## Computes lon and lat resolutions, and the coordinates for the 1st and last gridcells in each dimension.
txtfile_resol_init=$tmpdir/resol_init_${config}.txt
txtfile_resol_ref=$tmpdir/resol_ref_${config}.txt

# $rscript_dir/Rscript $pdir/compute_resolution.R $nlon_init $nlat_init $lon_dom_1 $lon_dom_2 $lat_dom_1 $lat_dom_2 $txtfile_resol_init
# $rscript_dir/Rscript $pdir/compute_resolution.R $nlon_ref $nlat_ref $lon_dom_1 $lon_dom_2 $lat_dom_1 $lat_dom_2 $txtfile_resol_ref
Rscript $pdir/compute_resolution.R $nlon_init $nlat_init $lon_dom_1 $lon_dom_2 $lat_dom_1 $lat_dom_2 $txtfile_resol_init
Rscript $pdir/compute_resolution.R $nlon_ref $nlat_ref $lon_dom_1 $lon_dom_2 $lat_dom_1 $lat_dom_2 $txtfile_resol_ref
## The next line prepares input parameters for the analysing routines in .../routines/comparison_IAGOS_models/. This script is not concerned.
mkdir -p $parent_dir/comparison_IAGOS_models/ ; cp $txtfile_resol_ref $parent_dir/comparison_IAGOS_models/
lon_resol_init=$(awk '{printf "%8.6f",$1}' ${txtfile_resol_init})
lat_resol_init=$(awk '{printf "%8.6f",$2}' ${txtfile_resol_init})
lon1_init=$(awk '{printf "%9.6f",$3}' ${txtfile_resol_init})
lon2_init=$(awk '{printf "%9.6f",$4}' ${txtfile_resol_init})
lat1_init=$(awk '{printf "%9.6f",$5}' ${txtfile_resol_init})
lat2_init=$(awk '{printf "%9.6f",$6}' ${txtfile_resol_init})

lon_resol_ref=$(awk '{printf "%8.6f",$1}' ${txtfile_resol_ref})
lat_resol_ref=$(awk '{printf "%8.6f",$2}' ${txtfile_resol_ref})
lon1_ref=$(awk '{printf "%9.6f",$3}' ${txtfile_resol_ref})
lon2_ref=$(awk '{printf "%9.6f",$4}' ${txtfile_resol_ref})
lat1_ref=$(awk '{printf "%9.6f",$5}' ${txtfile_resol_ref})
lat2_ref=$(awk '{printf "%9.6f",$6}' ${txtfile_resol_ref})

###### Vertical axis: sets the different layers properties.
## First defining variables for the prefix.
lower_PV_char='lower_PV_'   ; upper_PV_char='upper_PV_'
lower_P_char='lower_P_'     ; upper_P_char='upper_P_'
ozone_min_char='ozone_min_' ; ozone_max_char='ozone_max_'
## Then evaluating the corresponding limits, depending on the layer.
list_lower_PV_layer=""
list_upper_PV_layer=""
list_lower_P_layer=""
list_upper_P_layer=""
list_ozone_min_layer=""
list_ozone_max_layer=""
for layer in $list_layers ; do
    ## Defining the name.
    eval lower_PV_layer=\$$lower_PV_char$layer
    eval upper_PV_layer=\$$upper_PV_char$layer
    eval lower_P_layer=\$$lower_P_char$layer
    eval upper_P_layer=\$$upper_P_char$layer
    eval ozone_min_layer=\$$ozone_min_char$layer
    eval ozone_max_layer=\$$ozone_max_char$layer
    list_lower_PV_layer=${list_lower_PV_layer}" "$lower_PV_layer
    list_upper_PV_layer=${list_upper_PV_layer}" "$upper_PV_layer
    list_lower_P_layer=${list_lower_P_layer}" "$lower_P_layer
    list_upper_P_layer=${list_upper_P_layer}" "$upper_P_layer
    list_ozone_min_layer=${list_ozone_min_layer}" "$ozone_min_layer
    list_ozone_max_layer=${list_ozone_max_layer}" "$ozone_max_layer
done
nlayers=0 ; for layer in $list_layers ; do nlayers=$((nlayers+1)) ; done

######## Source files ########
### Dynamical fields repertory
ef_values_dir=$CCCSCRATCHDIR/Pack_code/data/reference_values/$model ## Useful for reading the lists of IAGOS files to be read.
hybrid_coef_file_input=./hybrid_coefs_$nlev_init.nc ## A_i and B_i coefficients.
model_output_dir=$CCCSCRATCHDIR/Pack_code/data/models_output/$model/$config/$experiment/$time_resolution ## Here are the pressure fields.
if [ "$ref_grid" = cst_P_UTLS ] ; then
    interp_output_dir=$CCCSCRATCHDIR/Pack_code/data/models_output/$model/$config/$experiment/$time_resolution/$ref_grid
    mkdir -p $interp_output_dir
fi
### IAGOS data repertory
iagos_data_dir=$CCCSCRATCHDIR/Pack_code/data/IAGOS_database
iagos_data_dir_format=$iagos_data_dir/NetCDF
##############################
### Lists of the IAGOS files to be read. There is one list per month.
list_iagos_dir=$tmpdir/list_iagos_files ; mkdir -p $list_iagos_dir
### Lists of input variables for the fortran routines.
list_file_extract=$pdir/namelist_extract
list_file_interp=$pdir/namelist_interp
list_file_stock=$tmpdir/namelist_stock

echo -e "\e[00;32m Experiment:                    \e[01;32m $experiment \e[00m"
echo -e "\e[00;32m Configuration:                 \e[01;32m $config \e[00m"
echo -e "\e[00;32m Model:                         \e[01;32m $model \e[00m"
echo -e "\e[00;32m Period:                        \e[01;32m $yyyymm1 \e[00;32m until \e[01;32m $yyyymm2 \e[00m"
echo -e "\e[00;32m Tropopause PV:                 \e[01;32m $PV_tropopause PVU \e[00m"
echo -e "\e[00;32m Time resolution for tropopause:\e[01;32m $time_resolution \e[00m"
echo -e "\e[00;32m Layers:                        \e[01;32m $list_layers \e[00m"
if [ "$bin_seasonal_climatology" = TRUE ] | [ "$bin_average_layers" = TRUE ] ; then
    echo -e "\e[00;32m Layers for climatologies:  \e[01;32m     $list_layers_climatologies \e[00m"
fi


echo -e "\e[00;37m=============================================================\e[00m"
echo -e "\e[00;37m==================== START GRIDDING IAGOS ===================\e[00m"
echo -e "\e[00;37m=============================================================\e[00m"

### Here a chrono starts, in order to compute how much time the run takes.
function convertsec(){ # Convert seconds for chrono routines.
    echo $((${1}/86400))':'`date +%H:%M:%S -d "0000-01-01 $(($1 % 86400)) seconds"`
}
start_run=`date +%s`

#### Generates the list of months from yyyymm1 to yyyymm2.
yyyymm=$yyyymm1
list_months=$yyyymm
nmonths=1
while [ $yyyymm -lt $yyyymm2 ] ; do
    yyyymm=`date -d "${yyyymm}01 1 months" +%Y%m`
    list_months=${list_months}" "${yyyymm}
    nmonths=$((nmonths+1))
done

#### Generates the list of years from yyyy1 to yyyy2.
yyyy1=${yyyymm1:0:4} ; yyyy2=${yyyymm2:0:4}
yyyy=$yyyy1
listyears=$yyyy
while [ $yyyy -lt $yyyy2 ] ; do
    yyyy=`date -d "${yyyy}0101 1 years" +%Y`
    listyears=${listyears}' '$yyyy
done

## Creating the yearly directories.
for layer in $list_layers ; do
    for yyyy in $listyears ; do
	mkdir -p $outdir_interpol_iagos/$layer/$yyyy
    done
done

## Extracting the PV and the pressure fields.
## The following two loops are an instruction for the fortran routines.
if [ "$extract_PV" = 'yes' ] ; then
    logical_PV_fortran=.TRUE.
else
    logical_PV_fortran=.FALSE.
fi
if [ "$extract_P" = "yes" ] ; then
    logical_P_fortran=.TRUE.
else
    logical_P_fortran=.FALSE.
fi
## Also generating a set of variables if these ones are absent in the model files.
model_var_header=$(head -1 $pdir/table_variable_names_in_models.txt)
icolmodel=1
current_column_in_header=$(echo $model_var_header | awk -v icolmodel=$icolmodel '{print $icolmodel}')
## Finds the right column in the variables tabular for the input name of the variable (the one in the initial model file).
## Note that if the program is CMIP6 (that must be indicated in the model name), then we copy the 'CMIP6' column. Useful because it avoids to repeat n_models_in_CMIP6 times the same column in the Table_variable files.
if [ "$sci_program" = CMIP6 ] ; then
    column_to_find=CMIP6
elif  [ "$sci_program" = ACACIA ] & [ "$model" != INCA ] ; then
    column_to_find=ACACIA
else
    column_to_find=$model
fi
while [ "$current_column_in_header" != "$column_to_find" ] ; do
    icolmodel=$((icolmodel+1))
    current_column_in_header=$(echo $model_var_header | awk -v icolmodel=$icolmodel '{print $icolmodel}')
done
icol_output_models=1
current_column_in_header=$(echo $model_var_header | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
## Finds the right column in the variables tabular for the output name of the variable.
while [ "$current_column_in_header" != 'Output_name' ] ; do
    icol_output_models=$((icol_output_models+1))
    current_column_in_header=$(echo $model_var_header | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
done
## These next lines do not describe an exhaustive list of the variables. 
## It only accounts for the variables that are exclusively involved in some special operations.
## For example, for deriving the O3/CO ratio, or the missing NOx variable in case the other two exist.
lon_dimname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Longitude' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
lat_dimname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Latitude' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
lev_dimname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Vertical_grid_level' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
O3_to_CO_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Ozone_to_CO_ratio' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
O3_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Ozone' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
CO_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Carbon_monoxide' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
NO_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Nitrogen_monoxide' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
NO2_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Nitrogen_dioxide' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
NOx_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Nitrogen_oxides' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
NOy_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Nitrogen_reactive_species' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
PAN_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Peroxyacetyl_nitrate' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
HNO3_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Nitric_acid' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
O3prod_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Production_rate_ozone' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
O3loss_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Loss_rate_ozone' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
O3net_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Net_production_ozone' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
O3S_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Stratospheric_ozone' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
O3S_to_O3_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Stratospheric_ozone_ratio' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
O3I_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Inert_strat_ozone' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
O3I_to_O3_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Inert_strat_ozone_ratio' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
NOx_to_NOy_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'NOx_to_NOy_ratio' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
PAN_to_NOy_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'PAN_to_NOy_ratio' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
HNO3_to_NOy_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'HNO3_to_NOy_ratio' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
O3T_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Tropospheric_ozone' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
## Of course, the PV and pressure names.
PV_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Potential_vorticity' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
P_S_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Surface_pressure' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
P_varname=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Pressure' | awk -v icolmodel=$icolmodel '{print $icolmodel}')
### O3 variables, differing by their origin layer.
O3_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Ozone' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
O3_to_CO_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Ozone_to_CO_ratio' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
O3S_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Stratospheric_ozone' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
O3T_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Tropospheric_ozone' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
O3S_to_O3_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Stratospheric_ozone_ratio' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
O3I_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Inert_strat_ozone' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
O3I_to_O3_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Inert_strat_ozone_ratio' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
NOx_to_NOy_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'NOx_to_NOy_ratio' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
PAN_to_NOy_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'PAN_to_NOy_ratio' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
HNO3_to_NOy_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'HNO3_to_NOy_ratio' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
### Rate variables concerning ozone. Production, loss and net production.
O3prod_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Production_rate_ozone' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
O3loss_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Loss_rate_ozone' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
O3net_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Net_production_ozone' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
### Time. Has to be the same as in monthly IAGOS files.
time_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Time' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
### Spatial coordinates. Same comment as previously.
lon_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Longitude' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
lat_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Latitude' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
level_varname_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Vertical_grid_level' | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')

## Extracting some model variables, and calculating new model variables.
if [ "$bin_extract_last_variables" = TRUE ] ; then
    for yyyymm in $list_months ; do
	yyyy=${yyyymm:0:4}
	model_file=$model_output_dir/$yyyy/${experiment}_${yyyymm}.nc
	model_PV_file=$model_output_dir/$yyyy/PV_${yyyymm}.nc
	model_P_S_file=$model_output_dir/$yyyy/P_S_${yyyymm}.nc
	model_P_file=$model_output_dir/$yyyy/P_${yyyymm}.nc
	if [ ! -f "$model_PV_file" ] && [ "$extract_PV" = "yes" ] ; then
	    echo -e "\e[00;37m Extracting PV: \e[00m" $yyyymm
    	    ncks -h -O -v $PV_varname $model_file -o $model_PV_file
	    ncks -h -O -x -v $PV_varname $model_file -o $model_file
	fi
	## Extracting surface pressure P_S. It is useful to reproduce the vertical grid levels (sigma-P), thus
	## only if there is not any pressure variable in the model output (hence the "$extract_P" = "no" condition).
	if [ ! -f "$model_P_S_file" ] && [ "$extract_P" = "no" ]; then
	    echo -e "\e[00;37m Extracting surface pressure: \e[00m" $yyyymm
	    ncks -h -O -v $P_S_varname $model_file -o $model_P_S_file
	    ncks -h -O -x -v $P_S_varname $model_file -o $model_file
	fi
	if [ ! -f "$model_P_file" ] && [ "$extract_P" = "yes" ] ; then
	    echo -e "\e[00;37m Extracting pressure: \e[00m" $yyyymm
	    ncks -h -O -v $P_varname $model_file -o $model_P_file
	    ncks -h -O -x -v $P_varname $model_file -o $model_file
	fi
	o3toco_exists=$(ncdump -h $model_file | grep -E -w $O3_to_CO_varname | wc -l)
	o3toco_queried=$(echo $list_model_variables | grep "Ozone_to_CO_ratio" | wc -l)
	if [ $o3toco_exists -lt 1 ] && [ $o3toco_queried -ge 1 ] ; then
    	    ozone_exists=$(ncdump -h $model_file | grep -E -w $O3_varname | wc -l)
    	    CO_exists=$(ncdump -h $model_file | grep -E -w $CO_varname | wc -l)
    	    if [ $ozone_exists -ge 1 ] && [ $CO_exists -ge 1 ] ; then
		echo -e "\e[00;37m Computing O3/CO: \e[00m" $yyyymm
		ncap2 -h -A -s "${O3_to_CO_varname}=${O3_varname}/${CO_varname}" $model_file $model_file
    		echo -e '\e[00;31m O3-to-CO ratio calculated in the model file for \e[00m' $yyyymm
    	    fi
	fi
	## Now, deriving the NOx from the sum between NO and NO2.
	NOx_exists=$(ncdump -h $model_file | grep -E -w $NOx_varname | wc -l)
	NOx_queried=$(echo $list_model_variables | grep 'Nitrogen_oxides' | wc -l)
	if [ $NOx_exists -lt 1 ] && [ $NOx_queried -eq 1 ] ; then
    	    NO_exists=$(ncdump -h $model_file | grep -E -w $NO_varname | wc -l)
    	    NO2_exists=$(ncdump -h $model_file | grep -E -w $NO2_varname | wc -l)
    	    if [ $NO_exists -ge 1 ] && [ $NO2_exists -ge 1 ] ; then
		echo -e "\e[00;37m Computing NOx: \e[00m" $yyyymm
		ncap2 -h -A -s "${NOx_varname}=${NO_varname}+${NO2_varname}" $model_file $model_file
    		echo -e '\e[00;31m NOx calculated in the model file for \e[00m' $yyyymm
    	    fi
	fi
	## Derives the net ozone production in the model files, 
	## if not already done and if the Prod_O3 and Loss_O3 variables are available.
	O3net_queried=$(echo $list_model_variables | grep 'Net_production_ozone' | wc -l)
	if [ $O3net_queried -eq 1 ] ; then
	    O3net_exists=$(ncdump -h $model_file | grep -E -w $O3net_varname | wc -l)
	    if [ $O3net_exists -lt 1 ] ; then
    		O3prod_exists=$(ncdump -h $model_file | grep -E -w $O3prod_varname | wc -l)
    		O3loss_exists=$(ncdump -h $model_file | grep -E -w $O3loss_varname | wc -l)
    		if [ $O3prod_exists -ge 1 ] && [ $O3loss_exists -ge 1 ] ; then
		    echo -e "\e[00;37m Computing net O3 production: \e[00m" $yyyymm
		    ncap2 -h -A -s "${O3net_varname}=${O3prod_varname}-${O3loss_varname}" $model_file $model_file
    		    echo -e '\e[00;31m Net_production_ozone calculated in the model file for \e[00m' $yyyymm
    		fi
	    fi
	fi
	## Derives tropospheric-origin ozone in the model files, 
	## if not already done and if the O3 and O3S variables are available.
	O3T_queried=$(echo $list_model_variables | grep 'Tropospheric_ozone' | wc -l)
	if [ $O3T_queried -eq 1 ] ; then
	    O3T_exists=$(ncdump -h $model_file | grep -E -w $O3T_varname | wc -l)
	    if [ $O3T_exists -lt 1 ] ; then
    		O3_exists=$(ncdump -h $model_file | grep -E -w $O3_varname | wc -l)
    		O3S_exists=$(ncdump -h $model_file | grep -E -w $O3S_varname | wc -l)
    		if [ $O3_exists -ge 1 ] && [ $O3S_exists -ge 1 ] ; then
		    echo -e "\e[00;37m Computing tropospheric ozone: \e[00m" $yyyymm
		    ncap2 -h -A -s "${O3T_varname}=${O3_varname}-${O3S_varname}" $model_file $model_file
    		    echo -e '\e[00;31m Tropospheric_ozone calculated in the model file for \e[00m' $yyyymm
    		fi
	    fi
	fi
	## Derives the ratio between stratospheric-origin ozone and total ozone in the model files (O3S/O3), 
	## if not already done and if the O3 and O3S variables are available.
	O3S_to_O3_queried=$(echo $list_model_variables | grep 'Stratospheric_ozone_ratio' | wc -l)
	if [ $O3S_to_O3_queried -eq 1 ] ; then
	    O3S_to_O3_exists=$(ncdump -h $model_file | grep -E -w $O3S_to_O3_varname | wc -l)
	    if [ $O3S_to_O3_exists -lt 1 ] ; then
    		O3_exists=$(ncdump -h $model_file | grep -E -w $O3_varname | wc -l)
    		O3S_exists=$(ncdump -h $model_file | grep -E -w $O3S_varname | wc -l)
    		if [ $O3_exists -ge 1 ] && [ $O3S_exists -ge 1 ] ; then
		    ncap2 -h -A -s "${O3S_to_O3_varname}=${O3S_varname}/${O3_varname}" $model_file $model_file
    		    echo -e '\e[00;31m Stratospheric_ozone_ratio calculated in the model file for \e[00m' $yyyymm
    		fi
	    fi
	fi
	## Derives the ratio between inert stratospheric-origin ozone and total ozone in the model files (O3I/O3), 
	## if not already done and if the O3 and O3I variables are available.
	O3I_to_O3_queried=$(echo $list_model_variables | grep 'Inert_strat_ozone_ratio' | wc -l)
	if [ $O3I_to_O3_queried -eq 1 ] ; then
	    O3I_to_O3_exists=$(ncdump -h $model_file | grep -E -w $O3I_to_O3_varname | wc -l)
	    if [ $O3I_to_O3_exists -lt 1 ] ; then
    		O3_exists=$(ncdump -h $model_file | grep -E -w $O3_varname | wc -l)
    		O3I_exists=$(ncdump -h $model_file | grep -E -w $O3I_varname | wc -l)
    		if [ $O3_exists -ge 1 ] && [ $O3I_exists -ge 1 ] ; then
		    echo -e "\e[00;37m Computing inert-stratospheric ozone ratio: \e[00m" $yyyymm
		    ncap2 -h -A -s "${O3I_to_O3_varname}=${O3I_varname}/${O3_varname}" $model_file $model_file
    		    echo -e '\e[00;31m Inert_strat_ozone_ratio calculated in the model file for \e[00m' $yyyymm
    		fi
	    fi
	fi
	## Derives the ratio between HNO3 or PAN and NOy in the model files, 
	## if not already done and if the PAN, HNO3 and NOy variables are available.
	NOy_ratios_queried=$(echo $list_model_variables | grep 'NOx_to_NOy_ratio' | wc -l) ## By consistency, if there is one of the main three NOy, the other one should be included too.
	if [ $NOy_ratios_queried -eq 1 ] ; then
	    NOx_to_NOy_exists=$(ncdump -h $model_file | grep -E -w $NOx_to_NOy_varname | wc -l)
	    PAN_to_NOy_exists=$(ncdump -h $model_file | grep -E -w $PAN_to_NOy_varname | wc -l)
	    HNO3_to_NOy_exists=$(ncdump -h $model_file | grep -E -w $HNO3_to_NOy_varname | wc -l)
	    if [ $NOx_to_NOy_exists -lt 1 ] || [ $PAN_to_NOy_exists -lt 1 ] || [ $HNO3_to_NOy_exists -lt 1 ] ; then
    		NOy_exists=$(ncdump -h $model_file | grep -E -w $NOy_varname | wc -l)
		NOx_exists=$(ncdump -h $model_file | grep -E -w $NOx_varname | wc -l)
    		PAN_exists=$(ncdump -h $model_file | grep -E -w $PAN_varname | wc -l)
		HNO3_exists=$(ncdump -h $model_file | grep -E -w $HNO3_varname | wc -l)
    		if [ $NOy_exists -ge 1 ] && [ $NOx_exists -ge 1 ] && [ $PAN_exists -ge 1 ] && [ $HNO3_exists -ge 1 ] ; then
		    echo -e "\e[00;37m Computing individual NOy ratios: \e[00m" $yyyymm
		    ncap2 -h -A -s "${NOy_varname}=${NOx_varname}+${HNO3_varname}+${PAN_varname}" $model_file $model_file
		    ncap2 -h -A -s "${NOx_to_NOy_varname}=${NOx_varname}/${NOy_varname}" $model_file $model_file
		    ncap2 -h -A -s "${PAN_to_NOy_varname}=${PAN_varname}/${NOy_varname}" $model_file $model_file
		    ncap2 -h -A -s "${HNO3_to_NOy_varname}=${HNO3_varname}/${NOy_varname}" $model_file $model_file
    		    echo -e '\e[00;31m NOx_to_NOy_ratio, PAN_to_NOy_ratio and HNO3_to_NOy_ratio calculated in the model file for \e[00m' $yyyymm
    		fi
	    fi
	fi
    done
fi

if [ "$lon_axis_0_360" = TRUE ] || [ "$lat_axis_90N_90S" = TRUE ] || [ "$P_axis_bottom_to_top" = TRUE ] ; then
    ## The latitude and vertical coordinates I used are indexed northwards and downwards.
    ## It corresponds to the MOCAGE configuration, the model I first used.
    ## The LMDZ-OR-INCA configuration, on the contrary, has a reverse indexing for latitude and vertical levels.
    ## For consistency, the script called right below aims at imposing the MOCAGE configuration to the INCA outputs.
    ## As the previous steps concerning model variables extraction, it is coded in a way that ensures the modifications 
    ## do not happen twice on the same files, so you do not have to deactivate the call of that subscript 
    ## when you launch the main script again.
    shift_lon_axis=$lon_axis_0_360
    reverse_lat_axis=$lat_axis_90N_90S
    reverse_P_axis=$P_axis_bottom_to_top
    echo "Since we work with ${model}: here we call adjust_model_output_files.sh for the files in:"
    echo $model_output_dir
    $pdir/adjust_model_output_files.sh $model $config $experiment $sci_program $yyyymm1 $yyyymm2 $model_output_dir $lon_dimname $lat_dimname $lev_dimname $extract_PV $extract_P $shift_lon_axis $reverse_lat_axis $reverse_P_axis
fi

npackages=0 ; for package in $list_packages_IAGOS ; do npackages=$((npackages+1)) ; done
ncoord=0 ; for coord in $listcoord ; do ncoord=$((ncoord+1)) ; done
ndim_model=0 ; for dim in $list_model_dimensions ; do ndim_model=$((ndim_model+1)) ; done

### Output names for all the selected variables.
## First, finds the right column in the IAGOS table, as for the models.
IAGOS_var_header=$(head -1 ${pdir}/table_variable_names_in_IAGOS_files.txt)
icol_output_IAGOS=1
current_column_in_header=$(echo $IAGOS_var_header | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
while [ "$current_column_in_header" != "Output_name" ] ; do
    icol_output_IAGOS=$((icol_output_IAGOS+1))
    current_column_in_header=$(echo $IAGOS_var_header | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
done
listvar_IAGOS_output_name=''
for var in $list_observed_variables ; do
    var_output=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w $var | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
    listvar_IAGOS_output_name=$listvar_IAGOS_output_name' '$var_output
done
for var in $list_derived_variables ; do
    var_output=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w $var | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
    listvar_IAGOS_output_name=$listvar_IAGOS_output_name' '$var_output
done
### Creating a variable name for all the IAGOS variables, referring to their output version.
O3_varname_output_IAGOS=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w 'Ozone' | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
P_varname_output_IAGOS=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w 'Pressure' | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
T_varname_output_IAGOS=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w 'Temperature' | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
Theta_varname_output_IAGOS=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w 'Potential_temperature' | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
u_varname_output_IAGOS=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w 'Zonal_wind_speed' | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
v_varname_output_IAGOS=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w 'Meridional_wind_speed' | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
CO_varname_output_IAGOS=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w 'Carbon_monoxide' | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
O3toCO_varname_output_IAGOS=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w 'Ozone_to_CO_ratio' | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
H2O_varname_output_IAGOS=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w 'Water_vapour' | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
RHL_varname_output_IAGOS=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w 'Relative_humidity_liquid' | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
RHI_varname_output_IAGOS=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w 'Relative_humidity_ice' | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
Cloud_varname_output_IAGOS=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w 'Cloud_liquid_water_particles' | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')
NOy_varname_output_IAGOS=$(cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w 'Nitrogen_reactive_species' | awk -v icol_output_IAGOS=$icol_output_IAGOS '{print $icol_output_IAGOS}')

echo "icol_output_models=" $icol_output_models
listvar_model_output_name='' ; n_model_variables=0
for var in $list_model_variables ; do
    var_output=$(cat $pdir/table_variable_names_in_models.txt | grep -E -w $var | awk -v icol_output_models=$icol_output_models '{print $icol_output_models}')
    listvar_model_output_name=$listvar_model_output_name' '$var_output
    n_model_variables=$((n_model_variables+1))
done

n_observed_variables=0
list_observed_variables_nc="" ## The variable names as they will be written in the NetCDF output files.
## Generating the temporary variable_names table, reduced to the variables required by the user.
head -1 $pdir/table_variable_names_in_IAGOS_files.txt > $tmpdir/table_variable_names_in_IAGOS_files.txt
head -1 $pdir/table_variable_units_in_IAGOS_files.txt > $tmpdir/table_variable_units_in_IAGOS_files.txt

for coord in $listcoord ; do
    cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w $coord >> $tmpdir/table_variable_names_in_IAGOS_files.txt
    cat $pdir/table_variable_units_in_IAGOS_files.txt | grep -E -w $coord >> $tmpdir/table_variable_units_in_IAGOS_files.txt
done
for var in $list_observed_variables ; do
    n_observed_variables=$((n_observed_variables+1)) ## Counting the amount of observed variables.
    var_nc=$(grep -E -w $var $pdir/table_variable_names_in_IAGOS_files.txt | awk '{print $6}') ## Repairing the output name of the current variable.
    list_observed_variables_nc=$list_observed_variables_nc" "$var_nc ## Implementing the list of output variables.
    ## Creates a copy of the reference table for the variable names, selecting only the required variables.
    ## Also generating the temporary variable_units table, restrained too.
    cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w $var >> $tmpdir/table_variable_names_in_IAGOS_files.txt
    cat $pdir/table_variable_units_in_IAGOS_files.txt | grep -E -w $var >> $tmpdir/table_variable_units_in_IAGOS_files.txt    
done
list_derived_variables_nc=""
n_derived_variables=0
for var in $list_derived_variables ; do
    n_derived_variables=$((n_derived_variables+1))
    var_nc=$(grep -E -w $var $pdir/table_variable_names_in_IAGOS_files.txt | awk '{print $6}') ## Repairing the output name of the current variable.
    list_derived_variables_nc=$list_derived_variables_nc" "$var_nc
    cat $pdir/table_variable_names_in_IAGOS_files.txt | grep -E -w $var >> $tmpdir/table_variable_names_in_IAGOS_files.txt
    cat $pdir/table_variable_units_in_IAGOS_files.txt | grep -E -w $var >> $tmpdir/table_variable_units_in_IAGOS_files.txt
done
n_IAGOS_variables=$((n_observed_variables+n_derived_variables))

# Generating the list of variables for the model now. n_model_variables already exists.
n_model_dimensions=0
list_model_variables_nc="" ## The variable names as they will be written in the NetCDF output files.
## Generating the temporary variable_names table, restrained to the variables required by the user.
## Also generating the temporary variable_units table, restrained too.
model_var_header=$(head -1 $pdir/table_variable_names_in_models.txt)
head -1 $pdir/table_variable_names_in_models.txt > $tmpdir/table_variable_names_in_models.txt
n_col_header=$(echo $model_var_header | wc -w) ; n_models_registered=$((n_col_header-2))
model_units_header=$(head -1 $pdir/table_variable_units_in_models.txt)
head -1 $pdir/table_variable_units_in_models.txt > $tmpdir/table_variable_units_in_models.txt

## In case we are in the CMIP6 program (or any other one requiring homogeneous variable naming):
## in the temporary table, we replace the CMIP6 term by $model, which is named CMIP6_xx (xx = original model name).
if [ "$sci_program" = "CMIP6" ] || [ "$sci_program" = "ACACIA" ] ; then
    model_var_header_tmpfile=''
    ## 1/ Variable names.
    for column_title in $model_var_header ; do
    	if [ "$column_title" = "$sci_program" ] ; then
	    if [ "$model" = "INCA" ] ; then
    		column_title_tmpfile=$model
	    else
		column_title_tmpfile=$sci_program
	    fi
    	else
    	    column_title_tmpfile=$column_title
    	fi
    	model_var_header_tmpfile=$model_var_header_tmpfile" "$column_title_tmpfile
    done
    echo $model_var_header_tmpfile > $tmpdir/table_variable_names_in_models.txt
    ## 2/ Variable units.
    ## This time, CMIP6 is a suffix, not a prefix. 
    model_units_header_tmpfile=''
    for column_title in $model_units_header ; do
	nchar=$(echo -n $column_title | wc -m)
	nchar_sci_program=$(echo -n $sci_program | wc -m)
	first_letters=${sci_program:0:$((nchar-nchar_sci_program))}
	last_letters=${sci_program:$((nchar-nchar_sci_program)):$nchar_sci_program}
	# if [ $last_5_letters = CMIP6 ] || [ $last_6_letters = ACACIA ] ; then
	if [ "$last_letters" = CMIP6 ] ; then # Not with ACACIA, since the variables names must be homogenized.
	    column_title_tmpfile=${first_letters}${model}
	else
	    column_title_tmpfile=$column_title
	fi
	model_units_header_tmpfile=$model_units_header_tmpfile" "$column_title_tmpfile
    done
    echo $model_units_header_tmpfile > $tmpdir/table_variable_units_in_models.txt
fi
for dim in $list_model_dimensions ; do
    cat $pdir/table_variable_names_in_models.txt | grep -E -w $dim >> $tmpdir/table_variable_names_in_models.txt
    cat $pdir/table_variable_units_in_models.txt | grep -E -w $dim >> $tmpdir/table_variable_units_in_models.txt
    n_model_dimensions=$((n_model_dimensions+1))
done
for var in $list_model_variables ; do
    cat $pdir/table_variable_names_in_models.txt | grep -E -w $var >> $tmpdir/table_variable_names_in_models.txt
    cat $pdir/table_variable_units_in_models.txt | grep -E -w $var >> $tmpdir/table_variable_units_in_models.txt
done
## A last row for surface pressure. It is apart from both dimensions and variables, 
## cause it is only a temporary variable during the code execution. Not needed later.
cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Surface_pressure' >> $tmpdir/table_variable_names_in_models.txt
cat $pdir/table_variable_units_in_models.txt | grep -E -w 'Surface_pressure' >> $tmpdir/table_variable_units_in_models.txt
cat $pdir/table_variable_names_in_models.txt | grep -E -w 'Pressure' >> $tmpdir/table_variable_names_in_models.txt
cat $pdir/table_variable_units_in_models.txt | grep -E -w 'Pressure' >> $tmpdir/table_variable_units_in_models.txt
## The next line prepares input parameters for the analysing routines in .../routines/comparison_IAGOS_models/. This script is not concerned.
cp $tmpdir/table_variable_*.txt $parent_dir/comparison_IAGOS_models/

## This section interpolates the model grid onto another grid. 
## The grid called cst_P_UTLS has constant pressure levels with a focus on the cruise altitudes.
## No level in the deep stratosphere is treated here.
if [ "$bin_change_grid" = TRUE ] ; then
    echo -e "\e[00;31m Starting the interpolation of the model onto another grid \e[00m"
    ## Computes lon and lat resolutions, and the coordinates for the 1st and last gridcells in each dimension.
    txtfile_resol_ref=$tmpdir/resol_ref_${config}.txt
    # $rscript_dir/Rscript $pdir/compute_resolution.R $nlon_ref $nlat_ref $lon_dom_ref_1 $lon_dom_ref_2 $lat_dom_ref_1 $lat_dom_ref_2 $txtfile_resol_ref
    Rscript $pdir/compute_resolution.R $nlon_ref $nlat_ref $lon_dom_ref_1 $lon_dom_ref_2 $lat_dom_ref_1 $lat_dom_ref_2 $txtfile_resol_ref
    lon_resol_ref=$(awk '{printf "%8.6f",$1}' ${txtfile_resol_ref})
    lat_resol_ref=$(awk '{printf "%8.6f",$2}' ${txtfile_resol_ref})
    lon1_ref=$(awk '{printf "%9.6f",$3}' ${txtfile_resol_ref})
    lon2_ref=$(awk '{printf "%9.6f",$4}' ${txtfile_resol_ref})
    lat1_ref=$(awk '{printf "%9.6f",$5}' ${txtfile_resol_ref})
    lat2_ref=$(awk '{printf "%9.6f",$6}' ${txtfile_resol_ref})
    # exe_file=$pdir/exe_model_interp_${id_job}
    exe_file=$exe_dir/exe_model_interp
    check_file=$tmpdir/check_model_interp_${id_job} ; rm -f $check_file
    if [ "$bin_recompile" = TRUE ] ; then
	## Compilation, only once.
	ifort -g -DasuX $pdir/global_var_interpol_mod.f90 $pdir/functions_mod.f90 $pdir/log_P_interp.f90 $pdir/vertical_interp.f90 $pdir/eumetsat_ropp/ropp_io/ncdf/ncdf.f90 -o $exe_file $NETCDF_INCDIR $NETCDF_LIBDIR $NETCDF_LIB -free
	status=$?
	if [ $status -ne 0 ] ; then
	    echo "log_P_interp.f90: Compilation failed."
	    exit
	fi
    fi
    export FORT_FMT_NO_WRAP_MARGIN=true
    for yyyymm in $list_months ; do
	yyyy=${yyyymm:0:4}
	interp_output_dir_yyyy=$interp_output_dir/$yyyy ; mkdir -p $interp_output_dir_yyyy
	## Input files.
	interp_input_dir=$CCCSCRATCHDIR/Pack_code/data/models_output/$model/$config/$experiment/$time_resolution/$yyyy
	vmr_file_input=$interp_input_dir/${experiment}_${yyyymm}.nc
	pv_file_input=$interp_input_dir/PV_${yyyymm}.nc
	p_file_input=$interp_input_dir/P_${yyyymm}.nc
	p_0_file_input=$interp_input_dir/P_S_${yyyymm}.nc
	# Now writing the instructions into the namelist_interp file, to be read by interpol_iagos.f90.
	cat > $list_file_interp << EOF
$yyyymm $yyyy
$nlon_init $nlat_init $nlev_init
$lon_resol_init $lat_resol_init
$lon_dom_1 $lon_dom_2 $lat_dom_1 $lat_dom_2
$nlon_ref $nlat_ref $nlev_ref
$lon_resol_ref $lat_resol_ref
$lon1_ref $lon2_ref $lat1_ref $lat2_ref
$P_levels_Pa
$logical_PV_fortran
$logical_P_fortran
$ndim_model
$list_model_dimensions
$n_model_variables
$list_model_variables
$n_models_registered
$version_vertical_interp
'${ref_grid}'
'${model}'
'${experiment}'
'${tmpdir}/'
'${tmpdir_there}/'
'${ref_values_dir}/'
'${interp_input_dir}/'
'${vmr_file_input}'
'${pv_file_input}'
'${p_file_input}'
'${p_0_file_input}'
'${hybrid_coef_file_input}'
'${interp_output_dir_yyyy}/'
EOF
	echo -e "\e[00;31m List_file available as \e[00m" $list_file_interp
	echo -e "\e[00;31m This file represents the list of technical instructions asked to the program, depending on the user requirements. \e[00m"
	$exe_file
	exe_status=$?
	wait & # Waiting for the routine to finish.
	if [ $exe_status -eq 0 ] ; then
	    echo -e "\e[00;31m log_P_interp.f90 terminated from the $model model grid, month $yyyymm \e[00m"
	    # rm -f $exe_file
	else
	    echo -e "\e[00;31m Error during the log_P_interp.f90 execution, month $yyyymm. \e[00m"
	    exit
	fi
	cdo remapcon,$file_horiz_grid $interp_output_dir_yyyy/tmp_${experiment}_${yyyymm}.nc $interp_output_dir_yyyy/${experiment}_${yyyymm}.nc
	if [ "$extract_PV" == yes ] ; then
	    cdo remapcon,$file_horiz_grid $interp_output_dir_yyyy/tmp_PV_${yyyymm}.nc $interp_output_dir_yyyy/PV_${yyyymm}.nc
	    rm $interp_output_dir_yyyy/tmp_PV_${yyyymm}.nc
	fi
	cdo remapcon,$file_horiz_grid $interp_output_dir_yyyy/tmp_P_${yyyymm}.nc $interp_output_dir_yyyy/P_${yyyymm}.nc
	rm $interp_output_dir_yyyy/tmp_${experiment}_${yyyymm}.nc $interp_output_dir_yyyy/tmp_P_${yyyymm}.nc
    done
    # Then stops the chrono.
    end_run=`date +%s`
    time_spent=`convertsec $(($end_run-$start_run))`
    echo -e "\e[00;31m | END OF THE MODEL INTERPOLATION \e[00m (duration: ${time_spent})"
    echo duration: $time_spent > $tmpdir/duration_log_P_interp.txt
fi ##### Ends the inter-model interpolation

if [ "$bin_mask_NOy" = TRUE ] ; then
    ## Changing the tmp table_variable_names_in_models file. Before adding the mask_NOy variables,
    ## we have to reorganize the tmp table in order to keep the pressure and surface pressure at the end.
    ## Initially, the last two lines are the pressure and the surface pressure.
    nlines=$(cat $tmpdir/table_variable_names_in_models.txt | wc -l)
    head -$((nlines-2)) $tmpdir/table_variable_names_in_models.txt > $tmpdir/cp_table_variable_names_in_models.txt
    nlines=$(cat $tmpdir/table_variable_units_in_models.txt | wc -l)
    head -$((nlines-2)) $tmpdir/table_variable_units_in_models.txt > $tmpdir/cp_table_variable_units_in_models.txt
    listvar_mask_NOy="NOx_mask_NOy NO_mask_NOy NO2_mask_NOy HNO3_mask_NOy PAN_mask_NOy O3_mask_NOy"
    n_model_variables_mask_NOy=0
    for var in $listvar_mask_NOy ; do
	n_model_variables_mask_NOy=$((n_model_variables_mask_NOy+1))
	cat $pdir/table_variable_names_in_models.txt | grep -E -w $var >> $tmpdir/cp_table_variable_names_in_models.txt
	cat $pdir/table_variable_units_in_models.txt | grep -E -w $var >> $tmpdir/cp_table_variable_units_in_models.txt
    done
    ## Next: we add the mask_NOy variables into the model variable list.
    if [ $n_model_variables_mask_NOy -ne 0 ] ; then
	list_model_variables=$list_model_variables" "$listvar_mask_NOy
	listvar_model_output_name=$listvar_model_output_name" "$listvar_mask_NOy
    fi
    n_model_variables=$((n_model_variables+n_model_variables_mask_NOy))
    ## Last: we add the pressure and surface pressure rows, and replace the previous tmp table.
    tail -2 $tmpdir/table_variable_names_in_models.txt >> $tmpdir/cp_table_variable_names_in_models.txt
    tail -2 $tmpdir/table_variable_units_in_models.txt >> $tmpdir/cp_table_variable_units_in_models.txt
    mv $tmpdir/cp_table_variable_names_in_models.txt $tmpdir/table_variable_names_in_models.txt
    mv $tmpdir/cp_table_variable_units_in_models.txt $tmpdir/table_variable_units_in_models.txt
fi


if [ "$bin_interpol_iagos" = TRUE ] ; then
    ## Now we update the logical_P_fortran binary: there is no need of looking for a 3D pressure variable
    ## if the extract_P variable is false, or if the model output have been interpolated onto constant 
    ## pressure levels before.
    if [ "$extract_P" = "yes" ] & [ "$ref_grid" != "cst_P_UTLS" ] ; then
    	logical_P_fortran=.TRUE.
    else
    	logical_P_fortran=.FALSE.
    fi
    exe_file=$exe_dir/exe_iagos_interpolation_${id_job}
    rm -f $exe_file
    check_file=$tmpdir/check_interpol_iagos_${id_job}
    rm -f $check_file
    ifort -g -DasuX $pdir/global_var_interpol_mod.f90 $pdir/functions_mod.f90 $pdir/interpol_iagos.f90 $pdir/mean_values_interpol.f90 $pdir/eumetsat_ropp/ropp_io/ncdf/ncdf.f90 -o $exe_file $NETCDF_INCDIR $NETCDF_LIBDIR $NETCDF_LIB -free
    export FORT_FMT_NO_WRAP_MARGIN=true
    reinitialize_list_files=1
    if [ $reinitialize_list_files -eq 1 ] ; then
	rm -rf $list_iagos_dir
	mkdir -p $list_iagos_dir
    fi
    #### 1/ Month by month, lists the corresponding IAGOS files and stocks the list into a file.
    #### 2/ Also counts the amount of IAGOS files for each month.
    for yyyymm in $list_months ; do
	yyyy=${yyyymm:0:4}
	# 0/ Creates the list of days in the current month.
	yyyymmdd=${yyyymm}01
        this_month=${yyyymmdd:0:6}
	listdays=$yyyymmdd
	ndays=1
	## Computes the last date of the current month.
	next_month=`date -d "${yyyymm}01 1 months" +%Y%m`
	last_day=`date -d "${next_month}01 -1 days" +%Y%m%d`
	while [ $yyyymmdd -lt $last_day ] ; do
	    yyyymmdd=`date -d "$yyyymmdd 1 days" +%Y%m%d`
	    listdays=$listdays" "$yyyymmdd
	    this_month=${yyyymmdd:0:6}
	    ndays=$((ndays+1))
	done
	# 1/
	Tdir=$iagos_data_dir_format/$yyyy/$yyyymm
	list_iagos_files=$list_iagos_dir/$yyyy/$yyyymm
	list_iagos_nfiles_daily=$list_iagos_dir/$yyyy/list_nfiles_daily_${yyyymm}
	mkdir -p $list_iagos_dir/$yyyy
	## If $quicktest_interpol_IAGOS is set at 1, the script just picks the first file
	## from each month data in order to accelerate the tests while sounding the whole period.
	if [ $quicktest_interpol_IAGOS -eq 1 ] ; then
	    Tdir=$(ls $Tdir | head -1)
	    echo $Tdir > $list_iagos_files
	    echo 1 > $list_iagos_nfiles_daily
	    nfiles=1 ; nfiles_daily=1
	else
	    ## daily resolution
	    ## There is no need to reproduce all that was made for the monthly resolution.
	    ## The only thing we have to do here consists in assigning to each day its amount of data files.
	    for yyyymmdd in $listdays ; do
		nfiles_daily=$(ls $Tdir/*${yyyymmdd}* | wc -l)
		if [ ${yyyymmdd:6:2} -eq 01 ] ; then
		    list_nfiles_daily=$nfiles_daily
		else
		    list_nfiles_daily=$list_nfiles_daily" "$nfiles_daily
		fi
	    done
	    echo $list_nfiles_daily > $list_iagos_nfiles_daily
	    ## monthly resolution
	    ls $Tdir > $list_iagos_files
	    nfiles=$(ls $Tdir | wc -l)
	fi
	# 2/
	if [ $yyyymm -eq $yyyymm1 ] ; then
	    list_nfiles=$nfiles
	else
	    list_nfiles=$list_nfiles" "$nfiles
	fi
    done
    ## Here: do we want the model files with the original grid, or the new one, if it exists?
    init_model_output_dir=$model_output_dir
    if [ "$ref_grid" = cst_P_UTLS ] ; then
	model_output_dir=$interp_output_dir
	logical_cst_P_fortran=.TRUE.
	else
	logical_cst_P_fortran=.FALSE.
    fi
    # Another relevant step, useful for reading the day of the current IAGOS file, 
    # thus for computing the first and last sampled days in the current month:
    # detects what preceeds the date in the IAGOS file name.
    # Then, it tells to the fortran routine what this prefix is.
    # It is quite important: in case a change happens in the file names between two versions of the database,
    # the user will not have to change the code.
    # We assume that the prefix is the same in the whole data set. Thus, the only condition is to avoid mixing
    # two different versions of the data set.
    list_iagos_files=$(ls $list_iagos_dir/$yyyy2/$yyyymm2) # The last month is sufficient.
    filename=$(head -1 $list_iagos_files)
    i_first_char=0
    while [ ${filename:i_first_char:6} != $yyyymm2 ] ; do i_first_char=$((i_first_char+1)) ; done
    filename_prefix_iagos=${filename:0:i_first_char}
    filename_suffix_iagos=_L2_3.1.0
    echo -e "\e[00;31m List of IAGOS files generated in \e[00m" $list_iagos_dir
    # Now writing the instructions into the namelist_extract file, to be read by interpol_iagos.f90.
    cat > $list_file_extract <<EOF
$nlev_ref
$lat_resol_ref $lon_resol_ref
$nlat_ref $nlon_ref
$lat1_ref $lat2_ref
$lon1_ref $lon2_ref
$yyyymm1 $yyyymm2
$nmonths
$list_months
$list_nfiles
$ncoord
$listcoord
$n_observed_variables
$list_observed_variables
$list_observed_variables_nc
$n_derived_variables
$list_derived_variables
$list_derived_variables_nc
$logical_PV_fortran
$PV_tropopause
$logical_P_fortran
$P_0
$logical_cst_P_fortran
$nlayers
$list_layers
$list_lower_PV_layer
$list_upper_PV_layer
$list_lower_P_layer
$list_upper_P_layer
$list_ozone_min_layer
$list_ozone_max_layer
$n_model_dimensions
$list_model_dimensions
$n_model_variables
$list_model_variables
$n_models_registered
'${ref_values_dir}/'
'${model_output_dir}/'
'${init_model_output_dir}/'
'${iagos_data_dir_format}/'
'${tmpdir}/'
'${tmpdir_there}/'
'${file_cst_P_levels}'
'${outdir_interpol_iagos}/'
'${model}'
'${experiment}'
'${sci_program}'
'${time_resolution}'
'${filename_prefix_iagos}'
'${filename_suffix_iagos}'
$quicktest_interpol_IAGOS
$npackages
$list_packages_IAGOS
EOF

    echo -e "\e[00;31m List_file available as \e[00m" $list_file_extract
    echo -e "\e[00;31m This file represents the list of technical instructions asked to the program, depending on the user requirements. \e[00m"
    $exe_file
    exe_status=$?
    wait & # Waiting for the routine to finish.
    if [ $exe_status -eq 0 ] ; then
	echo -e "\e[00;31m interpol_iagos.f90 terminated on the $model model grid \e[00m"
	# rm -rf $list_iagos_dir
	rm -f $exe_file
    else
	echo -e "\e[00;31m Error during the interpol_iagos.f90 execution. \e[00m"
	# Then stops the chrono.
	end_run=`date +%s`
	time_spent=`convertsec $(($end_run-$start_run))`
	echo -e "\e[00;31m | END OF TRE PROGRAM \e[00m (duration: ${time_spent})"
	echo duration: $time_spent > $tmpdir/duration.txt
	exit
    fi
fi ##### Ends the IAGOS interpolation

## The pressure as a variable is not useful anymore, since we have it in attribute for each variable.
listvar_IAGOS_without_pressure=""
n_IAGOS_variables=0
for var in $listvar_IAGOS_output_name ; do
    if [ "$var" != "$P_varname_output_IAGOS" ] ; then
	listvar_IAGOS_without_pressure=$listvar_IAGOS_without_pressure" "$var
	n_IAGOS_variables=$((n_IAGOS_variables+1))
    fi
done
listvar_IAGOS_output_name=$listvar_IAGOS_without_pressure

if [ "$bin_regional_time_series" = TRUE ] ; then
    echo -e "\e[00;31m | Preparing the regional files. \e[00m"
    ############ Regional pre-treatment #############
    ### The next line is the initialisation of the variable names referring to the regions boundaries.
    ### It has to be followed by the name of a region.
    lon01char='lon01_' ; lon02char='lon02_' ; lat01char='lat01_' ; lat02char='lat02_'
    rm $tmpdir/regions_def_${config}.txt
    rm $parent_dir/comparison_IAGOS_models/regions_def_${config}.txt
    for region in $list_regions ; do
	tmpdir_region=$tmpdir_there/$region ; mkdir -p $tmpdir_region
	outdir_region=$dir_regions_config_experiment_resol/$region ; mkdir -p $outdir_region
	for layer in $list_layers ; do
	    mkdir -p $outdir_region/$layer
	    mkdir -p $tmpdir_region/$layer
	    for yyyy in $listyears ; do
		mkdir -p $tmpdir_region/$layer/$yyyy
		mkdir -p $outdir_region/$layer/$yyyy
	    done
	done
	bin_neutral_layers=$(echo $list_layers | grep -e 'all_layers' -e 'UTLS' test.txt | wc -l)
	if [ $bin_neutral_layers -eq 1 ] ; then
	    if [ $(echo $list_layers | grep 'all_layers' | wc -l) -eq 1 ] ; then
		pv_outdir=$outdir_region/all_layers
	    elif [ $(echo $list_layers | grep 'UTLS' | wc -l) -eq 1 ] ; then
		pv_outdir=$outdir_region/UTLS
	    else
		echo -e "\e[00;31m Needing a neutral layer for PV or P fields. Stopping the execution. \e[00m"
		exit
	    fi
	fi
	### Defining the limits of the current region.
	eval lon01=\$$lon01char$region
	eval lon02=\$$lon02char$region
	eval lat01=\$$lat01char$region
	eval lat02=\$$lat02char$region
	echo -e "\e[00;31m $region >>> $lon01 $lon02 $lat01 $lat02 \e[00m"
	### Now they have to be slightly modified to match with the model grid. This is what the R routine compute_grid_indexes.R does.
	tmpfile_indexes=$tmpdir/indexes.txt
	# $rscript_dir/Rscript $pdir/compute_grid_indexes.R $lon01 $lon02 $lat01 $lat02 $lon1_ref $lon2_ref $lat1_ref $lat2_ref $lon_resol_ref $lat_resol_ref $model $tmpfile_indexes
	Rscript $pdir/compute_grid_indexes.R $lon01 $lon02 $lat01 $lat02 $lon1_ref $lon2_ref $lat1_ref $lat2_ref $lon_resol_ref $lat_resol_ref $model $tmpfile_indexes
	### First set of outputs from the R routine: the indexes locating the new regions coordinates in the model grid.
	ilon01=$(awk '{print $1 }' ${tmpfile_indexes})
	ilon02=$(awk '{print $2 }' ${tmpfile_indexes})
	klat01=$(awk '{print $3 }' ${tmpfile_indexes})
	klat02=$(awk '{print $4 }' ${tmpfile_indexes})
	### Second set of outputs from the R routine: the corresponding coordinates in degrees East and North.
	lon01_grid=$(awk '{printf "%8.5f",$5 }' ${tmpfile_indexes})
	lon02_grid=$(awk '{printf "%8.5f",$6 }' ${tmpfile_indexes})
	lat01_grid=$(awk '{printf "%8.5f",$7 }' ${tmpfile_indexes})
	lat02_grid=$(awk '{printf "%8.5f",$8 }' ${tmpfile_indexes})
	echo $region $ilon01 $ilon02 $klat01 $klat02 $lon01_grid $lon02_grid $lat01_grid $lat02_grid >> $tmpdir/regions_def_${config}.txt
	## The next line prepares input parameters for the post-treatment routines in .../comparison_IAGOS_models/. 
	## The current script is not concerned.
	echo $region $ilon01 $ilon02 $klat01 $klat02 $lon01_grid $lon02_grid $lat01_grid $lat02_grid >> $parent_dir/comparison_IAGOS_models/regions_def_${config}.txt
	echo  -e "\e[00;32m | $region  >>>  $ilon01 $ilon02 $klat01 $klat02 $lon01_grid $lon02_grid $lat01_grid $lat02_grid \e[00m"
	for yyyymm in $list_months ; do
	    yyyy=${yyyymm:0:4}
	    input_file_model_P_S=$model_output_dir/$yyyy/P_S_${yyyymm}.nc
	    input_file_model_P=$model_output_dir/$yyyy/P_${yyyymm}.nc
	    input_file_model_pv=$model_output_dir/$yyyy/PV_${yyyymm}.nc
	    outfile_regional_model_P=$pv_outdir/$yyyy/P_${yyyymm}.nc
	    outfile_regional_model_P_S=$pv_outdir/$yyyy/P_S_${yyyymm}.nc
	    outfile_regional_model_pv=$pv_outdir/$yyyy/PV_${yyyymm}.nc
	    ## First preparing the PV file.
	    for layer in $list_layers ; do
		input_file_IAGOS=$outdir_interpol_iagos/$layer/$yyyy/IAGOS_${yyyymm}.nc
		input_file_model=$outdir_interpol_iagos/$layer/$yyyy/${experiment}_${yyyymm}.nc
		outfile_regional_IAGOS=$outdir_region/$layer/$yyyy/IAGOS_${yyyymm}.nc
		outfile_regional_model=$outdir_region/$layer/$yyyy/${experiment}_${yyyymm}.nc
		tmpfile_regional_IAGOS_CO_variables=$tmpdir_region/$layer/$yyyy/IAGOS_CO_variables.nc
		tmpfile_regional_model_CO_variables=$tmpdir_region/$layer/$yyyy/${model}_${experiment}_CO_variables.nc
		
		# Extracting the data points at the selected locations
		ncks -F -O -h -d lon,$ilon01,$ilon02,1 -d lat,$klat01,$klat02,1 $input_file_IAGOS -o $outfile_regional_IAGOS
		ncks -F -O -h -d lon,$ilon01,$ilon02,1 -d lat,$klat01,$klat02,1 $input_file_model -o $outfile_regional_model
		echo  -e "\e[00;36m | gridded regional file: \e[00m${outfile_regional_IAGOS}"
		echo  -e "\e[00;36m | gridded regional file: \e[00m${outfile_regional_model}"	
	    done
	done
    done
    echo  -e "\e[00;31m | Now computing regional time series. \e[00m" 
    exe_file=$exe_dir/exe_regional_average_IAGOS_${id_job}
    rm -f $exe_file
    ifort -g -DasuX $pdir/global_var_interpol_mod.f90 $pdir/functions_mod.f90 $pdir/regional_average_IAGOS.f90 $pdir/eumetsat_ropp/ropp_io/ncdf/ncdf.f90 -o $exe_file $NETCDF_INCDIR $NETCDF_LIBDIR $NETCDF_LIB -free
    export FORT_FMT_NO_WRAP_MARGIN=true
    for layer in $list_layers ; do
	outdir_time_series=$dir_regional_time_series/$layer ; mkdir -p $outdir_time_series
	for region in $list_regions ; do
	    echo -e "\e[00;32m | $layer $region \e[00m"  
	    cat > $pdir/namelist_regional <<EOF
$experiment
$yyyymm1 $yyyymm2
$nmonths
$list_months
$upper_level_UTLS $lower_level_UTLS
$upper_level_FT $lower_level_FT
$lon_varname_output $lat_varname_output $level_varname_output
$n_IAGOS_variables
$listvar_IAGOS_output_name
$n_model_variables
$listvar_model_output_name
$N_obs_min_regional
'${tmpdir_there}/'
$region
$layer
'${dir_regions_config_experiment_resol}/'
'${outdir_time_series}/'
EOF
	    # valgrind --leak-check=full --track-origins=yes $exe_file # One debug command.
	    $exe_file
	    exe_status=$?
	    wait & # Waiting for the routine to finish.
	    if [ $exe_status -eq 0 ] ; then # Stops the script if an error occurred during the routine execution.
		echo -e "\e[00;31m regional_average_iagos.f90 terminated on the $model model grid \e[00m"
	    else
		echo -e "\e[00;31m Error during the regional_average_iagos.f90 execution. \e[00m"
		exit
	    fi
	done
    done
    ## Cleaning the regional temporary files.
    for region in $list_regions ; do
	tmpdir_region=$tmpdir_there/$region ; rm -rf $tmpdir_region
	tmpdir_region_2=$dir_regions_config_experiment_resol/$region ; rm -rf $tmpdir_region_2
    done
fi

## Preliminary task 1: readjusting the first and the last month in a way we equilibrate
## the samples between the seasons. We start with DJF and end with SON.
yyyymm1_clim=${yyyy1}12 # Default: we readjust the first year to December in the same year.
case ${yyyymm1:4:2} in
    12 | 01 | 02) yyyymm1_clim=$yyyymm1 ;; # If winter is sampled during the first year, even partially, then it is ok.
esac
yyyymm2_clim=$((yyyy2-1))11 # Default: we readjust the last year to November in the previous year.
case ${yyyymm2:4:2} in
    09 | 10 | 11) yyyymm2_clim=$yyyymm2 ;; # If fall is sampled during the last year, even partially, then it is ok.
    12) yyyymm2_clim=${yyyy2}11 ;; # If the last month is December, then we readjust it to November in the same year.
esac
outdir_time_period=$outdir_level_by_level/${yyyymm1_clim}_${yyyymm2_clim}
mkdir -p $outdir_level_by_level ; mkdir -p $outdir_time_period

if [ "$bin_seasonal_climatology" = TRUE ] ; then
    echo -e "\e[00;31m | Preparing the seasonal files. \e[00m"
    ########## Locations of different directories ##########
    dir=$base_dir/Pack_code/routines/comparison_IAGOS_models
    data_dir=$base_dir/Pack_code/data
    echo -e "\e[00;32m First and last dates readjusted for interseasonal balance:\e[01;32m $yyyymm1_clim $yyyymm2_clim \e[00m"
    
    ## Preliminary task: preparing the variable list that will be used in the next loop.
    ## We only need mean values, SD and the amounts of data which name ends by _N_obs.
    ## Treats all the variables that are present in the input files with the suffix _Mean.
    ## Once for IAGOS, once for the model.
    ## Given that the list of variables is not supposed to change from one month to another
    ## (I am trusting you on this point), sounding the first month only is sufficient.
    for layer in $list_layers_climatologies ; do
	echo 'Layer: ' $layer
	monthly_IAGOS=$outdir_interpol_iagos/$layer/$yyyy1/IAGOS_${yyyymm1_clim}.nc
	monthly_model=$outdir_interpol_iagos/$layer/$yyyy1/${experiment}_${yyyymm1_clim}.nc
	## First, we generate a text containing the variable names once per line.
	grep_lines_IAGOS=$(ncdump -h $monthly_IAGOS | grep Mean | grep float)
	grep_lines_model=$(ncdump -h $monthly_model | grep Mean | grep float)
	## Second, we only keep the words containing the variable names.
	coarse_substr_array_IAGOS='' ; coarse_substr_array_model=''
	for line in $grep_lines_IAGOS ; do
	    coarse_substr_array_IAGOS=$coarse_substr_array_IAGOS' '$(echo $line | grep Mean)
	done
	for line in $grep_lines_model ; do
	    coarse_substr_array_model=$coarse_substr_array_model' '$(echo $line | grep Mean)
	done
	## Now, we can start the research of the variable name in each character string we extracted.
	## The information we have is that the variable name preceeds the substring _Mean.
	listvar_yyyymm1_clim_IAGOS='' ; listvar_yyyymm1_clim_model=''
	for line in $coarse_substr_array_IAGOS ; do
	    i_first_char=0
	    while [ "${line:i_first_char:5}" != "_Mean" ] ; do i_first_char=$((i_first_char+1)) ; done
	    if [ "${line:0:i_first_char}" != "$P_varname_output_IAGOS" ] ; then
		listvar_yyyymm1_clim_IAGOS=$listvar_yyyymm1_clim_IAGOS' '${line:0:i_first_char}
	    fi
	done
	for line in $coarse_substr_array_model ; do
	    i_first_char=0
	    while [ "${line:i_first_char:5}" != "_Mean" ] ; do i_first_char=$((i_first_char+1)) ; done
	    listvar_yyyymm1_clim_model=$listvar_yyyymm1_clim_model' '${line:0:i_first_char}
	done
	echo 'Variables in the first monthly IAGOS file:' $listvar_yyyymm1_clim_IAGOS
	echo 'Variables in the first monthly model file:' $listvar_yyyymm1_clim_model
	
	## Next step: preparing the list of coupled names (var_stat) for the relevant stats:
	## Mean, SD, and N_obs. The goal is to make a list with separating commas,
	## for the nco commands.
	listvarstat_Mean_IAGOS='' ; listvarstat_SD_IAGOS='' 
	listvarstat_Min_IAGOS=''  ; listvarstat_Max_IAGOS='' 
	listvarstat_N_obs_IAGOS=''; listvarstat_N_diff_IAGOS=''
	listvarstat_N_flights_IAGOS=''; listvarstat_P_IAGOS=''
	listvarstat_Mean_model='' ; listvarstat_SD_model=''
	listvarstat_Min_model=''  ; listvarstat_Max_model=''
	for var in $listvar_yyyymm1_clim_IAGOS ; do
	    listvarstat_Mean_IAGOS=$listvarstat_Mean_IAGOS' '${var}_Mean
	    listvarstat_SD_IAGOS=$listvarstat_SD_IAGOS' '${var}_SD
	    listvarstat_Min_IAGOS=$listvarstat_Min_IAGOS' '${var}_Min
	    listvarstat_Max_IAGOS=$listvarstat_Max_IAGOS' '${var}_Max
	    listvarstat_N_diff_IAGOS=$listvarstat_N_diff_IAGOS' '${var}_N_diff_updown
	    listvarstat_P_IAGOS=$listvarstat_P_IAGOS' '${var}_Pressure
	    listvarstat_N_obs_IAGOS=$listvarstat_N_obs_IAGOS' '${var}_N_obs
	    listvarstat_N_flights_IAGOS=$listvarstat_N_flights_IAGOS' '${var}_N_flights
	done
	for var in $listvar_yyyymm1_clim_model ; do
	    listvarstat_Mean_model=$listvarstat_Mean_model' '${var}_Mean
	    listvarstat_SD_model=$listvarstat_SD_model' '${var}_SD
	    listvarstat_Min_model=$listvarstat_Min_model' '${var}_Min
	    listvarstat_Max_model=$listvarstat_Max_model' '${var}_Max
	done
	listvarstat_Mean_IAGOS_nc=$(echo $listvarstat_Mean_IAGOS | sed s/' '/'\,'/g)
	listvarstat_SD_IAGOS_nc=$(echo $listvarstat_SD_IAGOS | sed s/' '/'\,'/g)
	listvarstat_Min_IAGOS_nc=$(echo $listvarstat_Min_IAGOS | sed s/' '/'\,'/g)
	listvarstat_Max_IAGOS_nc=$(echo $listvarstat_Max_IAGOS | sed s/' '/'\,'/g)
	listvarstat_N_obs_IAGOS_nc=$(echo $listvarstat_N_obs_IAGOS | sed s/' '/'\,'/g)
	listvarstat_N_diff_IAGOS_nc=$(echo $listvarstat_N_diff_IAGOS | sed s/' '/'\,'/g)
	listvarstat_N_flights_IAGOS_nc=$(echo $listvarstat_N_flights_IAGOS | sed s/' '/'\,'/g)
	listvarstat_P_IAGOS_nc=$(echo $listvarstat_P_IAGOS | sed s/' '/'\,'/g)
	listvarstat_Mean_model_nc=$(echo $listvarstat_Mean_model | sed s/' '/'\,'/g)
	listvarstat_SD_model_nc=$(echo $listvarstat_SD_model | sed s/' '/'\,'/g)
	listvarstat_Min_model_nc=$(echo $listvarstat_Min_model | sed s/' '/'\,'/g)
	listvarstat_Max_model_nc=$(echo $listvarstat_Max_model | sed s/' '/'\,'/g)

	echo '######################'
	echo 'Separating the seasons'
	echo '######################'
	rm -rf $tmpdir_there/monthly_*
	rm -rf $tmpdir_there/seasonal_*
	## Now defining the different seasonal averagings.
	## The following list corresponds to different yearly subdivisions into seasons, for the climatologies.
	## Each one is representative of at least one region in the world. The "summer" definition differs from one to another,
	## thus it is sufficient to determine which one we are talking about.
	## JJA is the boreal summer, classically. The whole year is subdivided into DJF, MAM, JJA and SON.
	## The next three definitions correspond to northmost and stabilized position of ITCZ in the tropics.
	## JAS corresponds to Atlantic coasts, JJASO to Africa (see Lannuque et al., 2020) and JJAS to South Asia.
	if [ "$layer" = "UTLS" ] || [ "$layer" = "LS" ] || [ "$layer" = "TPL" ] ; then
	    list_seasons_def=Boreal
	else
	    list_seasons_def=Boreal' '$list_seasons_def_tropics
	fi
	for seasons_def in $list_seasons_def ; do
	# for seasons_def in Tropical_Afr ; do
	    case $seasons_def in
		'Boreal' ) seasons_start_indexes='03_06_09_12' ;; ## Needs the 1st month of every season, whatever the order.
		'Tropical_Atl' ) seasons_start_indexes='03_07_09_12' ;;
		'Tropical_Afr' ) seasons_start_indexes='04_06_11_12' ;;
		'South_Asia' ) seasons_start_indexes='03_06_10_12' ;;
		'individual' ) list_seasons=Dec' 'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov
		    seasons_start_indexes='01_02_03_04_05_06_07_08_09_10_11_12' ;; ## Remains on the boreal model, but has no impact.
	    esac
	    list_seasons=`define_seasons $seasons_def`
	    for season in $list_seasons ; do
	    # for season in JJASO ; do
		case $season in
		    ##### Extratropical climate
		    'DJF' ) list_mm=12' '01' '02 ;;
		    'MAM' ) list_mm=03' '04' '05 ;;
		    'JJA' ) list_mm=06' '07' '08 ;;
		    'SON' ) list_mm=09' '10' '11 ;;
		    ##### Tropical climates
		    ## Atlantic coasts
		    ## DJF and SON already registered.
		    'MAMJ' ) list_mm=03' '04' '05' '06 ;;
		    'JA' ) list_mm=07' '08 ;;
		    ## Africa (Lannuque et al., 2021)
		    'DJFM' ) list_mm=12' '01' '02' '03 ;;
		    'AM' ) list_mm=04' '05 ;;
		    'JJASO' ) list_mm=06' '07' '08' '09' '10 ;;
		    'N' ) list_mm=11 ;;
		    ## Southeast Asia
		    ## DJF and MAM already registered.
		    'JJAS' ) list_mm=06' '07' '08' '09 ;;
		    'ON' ) list_mm=10' '11 ;;
		    ##### In case one wants to look at the individual months,
		    ## e.g. to test other definitions for seasons.
		    'Dec' ) list_mm=12 ;;
		    'Jan' ) list_mm=01 ;;
		    'Feb' ) list_mm=02 ;;
		    'Mar' ) list_mm=03 ;;
		    'Apr' ) list_mm=04 ;;
		    'May' ) list_mm=05 ;;
		    'Jun' ) list_mm=06 ;;
		    'Jul' ) list_mm=07 ;;
		    'Aug' ) list_mm=08 ;;
		    'Sep' ) list_mm=09 ;;
		    'Oct' ) list_mm=10 ;;
		    'Nov' ) list_mm=11 ;;
		esac
		echo 'Season: ' $season $list_mm
		# 0/ Cleaning the previous temporary files created by this loop.
		tmpdir_there_clim=${tmpdir_there}_$layer
		mkdir -p $tmpdir_there_clim
		rm -rf $tmpdir_there_clim/monthly_IAGOS_${layer}_${season}_*.nc
		rm -rf $tmpdir_there_clim/monthly_${experiment}_${layer}_${season}_*.nc

		## The two following conditions avoid useless operations: the seasons other than
		## DJF, ..., SON are defined in order to fit with the tropics, 
		## thus have nothing to do with the LS or the TPL in matter of IAGOS observations.
		if [ "$season" = DJF ] || [ "$season" = MAM ] || [ "$season" = JJA ] || [ "$season" = SON ] ; then
		    bin_extratropical_season=TRUE
		fi
		bin_tropospheric_layer=TRUE
		if [ "$layer" = LS ] || [ "$layer" = TPL ] ; then
		    bin_tropospheric_layer=FALSE
		fi
		if [ "$bin_tropospheric_layer" ] || [ "$bin_extratropical_season" ] ; then
		    for yyyy in $listyears ; do
			echo $yyyy
			# 1/ Creating a monthly-means file for each season and each year
			for mm in $list_mm ; do
			    yyyymm=${yyyy}$mm
			    ## Defining the year corresponding to the current season,
			    ## i.e. the previous one if the month is December for a DJF-winter,
			    ## and the next year if the month is January for an ONDJ-fall.
			    yyyy_seas=$yyyy
			    ## Defining input files
			    monthly_IAGOS=$outdir_interpol_iagos/$layer/$yyyy/IAGOS_${yyyymm}.nc
			    monthly_model=$outdir_interpol_iagos/$layer/$yyyy/${experiment}_${yyyymm}.nc

			    # 2/ Applying a mask with respect to the N_diff_updown metric.
			    echo -e "\e[00;31m Applying a mask with respect to the N_diff_updown metric. \e[00m"
			    # 2.1/ To the model grid cells.
			    tmpfile_IAGOS_masked=$tmpdir_there_clim/masked_IAGOS_${layer}_${yyyymm}.nc
			    tmpfile_model_masked=$tmpdir_there_clim/masked_${experiment}_${layer}_${yyyymm}.nc
			    cp $monthly_IAGOS $tmpfile_IAGOS_masked
			    cp $monthly_model $tmpfile_model_masked

			    ## Defining temporary files, outputs from step 1/
			    tmpfile_clim_IAGOS=$tmpdir_there_clim/monthly_IAGOS_${layer}_${season}_${yyyy_seas}.nc
			    tmpfile_clim_model=$tmpdir_there_clim/monthly_${experiment}_${layer}_${season}_${yyyy_seas}.nc
			    ## The following if loop just copies the file if the current month is
			    ## the first month of a season (3, 6, 9 or 12), as a seasonal initialization.
			    ## It also initializes the output file if it does not exist yet.
			    ## Else, it concatenates with the file corresponding to the previous seasons.
			    ## The additional 'm' is just a trick to avoid manipulating a 0-length string character.
			    if [ -f "$tmpfile_IAGOS_masked" ] ; then
				if [ $(echo $seasons_start_indexes | grep $mm)m != 'm' -o ! -f "$tmpfile_clim_IAGOS" ] ; then
				    cp $tmpfile_IAGOS_masked $tmpfile_clim_IAGOS
				    cp $tmpfile_model_masked $tmpfile_clim_model
				else
				    ncrcat -h -O $tmpfile_clim_IAGOS $tmpfile_IAGOS_masked $tmpfile_clim_IAGOS
				    ncrcat -h -O $tmpfile_clim_model $tmpfile_model_masked $tmpfile_clim_model
				fi
			    else
				echo $tmpfile_IAGOS_masked 'does not exist'
			    fi
			done
			echo -e "\e[00;31m Monthly files concatenated in: \e[00m" $tmpfile_clim_IAGOS
			echo -e "\e[00;31m and in: \e[00m" $tmpfile_clim_model
			# 2/ Concatenating $season files through the years. This step and the average can take a long time.
			## Defining the new temporary files, outputs from step 2/.
			tmpfile_concat_IAGOS=$tmpdir_there_clim/monthly_IAGOS_${layer}_${season}_${yyyymm1_clim}_${yyyymm2_clim}.nc
			tmpfile_concat_model=$tmpdir_there_clim/monthly_${experiment}_${layer}_${season}_${yyyymm1_clim}_${yyyymm2_clim}.nc
			if [ $yyyy -eq $yyyy1 -o ! -f "$tmpfile_concat_IAGOS" ] ; then
			    # Initializing the seasonal file.
			    cp $tmpdir_there_clim/monthly_IAGOS_${layer}_${season}_${yyyy}.nc $tmpfile_concat_IAGOS
			    cp $tmpdir_there_clim/monthly_${experiment}_${layer}_${season}_${yyyy}.nc $tmpfile_concat_model
			else
			    # Incrementing the seasonal file.
			    ncrcat -O $tmpfile_concat_IAGOS $tmpdir_there_clim/monthly_IAGOS_${layer}_${season}_${yyyy}.nc $tmpfile_concat_IAGOS
			    ncrcat -O $tmpfile_concat_model $tmpdir_there_clim/monthly_${experiment}_${layer}_${season}_${yyyy}.nc $tmpfile_concat_model
			fi
			echo -e "\e[00;31m $yyyy added to the seasonal files: \e[00m" $tmpfile_concat_IAGOS
			echo -e "\e[00;31m and:                               \e[00m" $tmpfile_concat_model
		    done
		    
		    
		    tmpfile_averaged_IAGOS=$tmpdir_there_clim/averaged_IAGOS_${layer}_${season}_${yyyymm1_clim}_${yyyymm2_clim}.nc
		    tmpfile_averaged_model=$tmpdir_there_clim/averaged_${experiment}_${layer}_${season}_${yyyymm1_clim}_${yyyymm2_clim}.nc
		    # 3/ Averaging all over the months, for each of the 4 seasons. Summing up for the N_obs stat.
		    echo -e "\e[00;32m Averaging into: \e[00m" $tmpfile_averaged_IAGOS
		    echo -e "\e[00;32m and in:         \e[00m" $tmpfile_averaged_model
		    # echo $listvarstat_Mean_IAGOS_nc
		    
		    ncra -h -O -v $listvarstat_Mean_IAGOS_nc $tmpfile_concat_IAGOS -o $tmpfile_averaged_IAGOS
		    # ncra -h -A -v $listvarstat_SD_IAGOS_nc $tmpfile_concat_IAGOS -o $tmpfile_averaged_IAGOS
		    # ncra -h -A -v $listvarstat_Min_IAGOS_nc $tmpfile_concat_IAGOS -o $tmpfile_averaged_IAGOS
		    # ncra -h -A -v $listvarstat_Max_IAGOS_nc $tmpfile_concat_IAGOS -o $tmpfile_averaged_IAGOS
		    ncra -h -A -v $listvarstat_N_obs_IAGOS_nc -y ttl $tmpfile_concat_IAGOS -o $tmpfile_averaged_IAGOS
		    ncra -h -A -v $listvarstat_N_diff_IAGOS_nc $tmpfile_concat_IAGOS -o $tmpfile_averaged_IAGOS
		    ncra -h -A -v $listvarstat_N_flights_IAGOS_nc -y ttl $tmpfile_concat_IAGOS -o $tmpfile_averaged_IAGOS
		    ncra -h -A -v $listvarstat_P_IAGOS_nc $tmpfile_concat_IAGOS -o $tmpfile_averaged_IAGOS

		    ncra -h -O -v $listvarstat_Mean_model_nc $tmpfile_concat_model -o $tmpfile_averaged_model
		    outfile_seasonal_clim_IAGOS=$outdir_time_period/IAGOS_${layer}_${season}_${yyyymm1_clim}_${yyyymm2_clim}.nc
		    outfile_seasonal_clim_model=$outdir_time_period/${experiment}_${layer}_${season}_${yyyymm1_clim}_${yyyymm2_clim}.nc	

		    mv $tmpfile_averaged_IAGOS $outfile_seasonal_clim_IAGOS
		    mv $tmpfile_averaged_model $outfile_seasonal_clim_model

		    # 5/ Concatenating the averages over the 4 seasons, then deriving the yearly mean values.
		    if [ "$seasons_def" = "Boreal" ] ; then
			echo '########################################################'
			echo 'Preparing yearly means: concatenating the seasonal files'
			echo '########################################################'
			tmpfile_all_seasons_gathered_IAGOS=$tmpdir_there_clim/all_seasons_IAGOS_${layer}_${yyyymm1_clim}_${yyyymm2_clim}.nc
			tmpfile_all_seasons_gathered_model=$tmpdir_there_clim/all_seasons_${experiment}_${layer}_${yyyymm1_clim}_${yyyymm2_clim}.nc
			if [ "$season" = DJF ] ; then
			    rm -rf $outdir_time_period/all_seasons_IAGOS_${layer}_${yyyymm1_clim}_${yyyymm2_clim}.nc $outdir_time_period/all_seasons_${experiment}_${layer}_${yyyymm1_clim}_${yyyymm2_clim}.nc
			    cp $outfile_seasonal_clim_IAGOS $tmpfile_all_seasons_gathered_IAGOS
			    cp $outfile_seasonal_clim_model $tmpfile_all_seasons_gathered_model
			else
			    ncrcat -O $tmpfile_all_seasons_gathered_IAGOS $outfile_seasonal_clim_IAGOS $tmpfile_all_seasons_gathered_IAGOS
			    ncrcat -O $tmpfile_all_seasons_gathered_model $outfile_seasonal_clim_model $tmpfile_all_seasons_gathered_model
			fi
			if [ "$season" = SON ] ; then
			    echo '#############################################################'
			    echo 'Computing yearly climatologies by averaging the seasonal ones'
			    echo '#############################################################'
			    ANN_outfile_IAGOS=$outdir_time_period/IAGOS_${layer}_ANN_${yyyymm1_clim}_${yyyymm2_clim}.nc
			    ANN_outfile_model=$outdir_time_period/${experiment}_${layer}_ANN_${yyyymm1_clim}_${yyyymm2_clim}.nc
			    ## At this step, all the var_stats in the tmpfiles are relevant. We can take them all without formulating lists again.
			    ## Only two exceptions: the amount of measurements (_N_obs suffix) and of flights (_N_flights suffix), that have to be summed up instead of averaged.
			    ncra -O -x -v $listvarstat_N_obs_IAGOS_nc,$listvarstat_N_flights_IAGOS_nc $tmpfile_all_seasons_gathered_IAGOS -o $ANN_outfile_IAGOS
			    ncra -O $tmpfile_all_seasons_gathered_model -o $ANN_outfile_model
			    ncra -A -v $listvarstat_N_obs_IAGOS_nc,$listvarstat_N_flights_IAGOS_nc -y ttl $tmpfile_all_seasons_gathered_IAGOS -o $ANN_outfile_IAGOS
			fi
		    fi ## if [ $seasons_def = 'Boreal' ] ; then
		fi ## if [ $bin_tropospheric_layer ] || [ $bin_extratropical_season ] ; then
	    done ## seasons
	done ## seasons_def
	rm -rf $tmpdir_there_clim/*.nc
	rm -rf $tmpdir_there_clim/*.tmp
	echo -e "\e[00;31m Yearly and seasonal files concatenated in: \e[00m" $outdir_time_period
    done ## layer
fi

outdir_seasonal_clim=$global_outdir_resol/seasonal_climatologies ## Should already exist, no need to create it then.
yyyymm1_clim=${yyyy1}12 # Default: we readjust the first year to December in the same year.
case ${yyyymm1:4:2} in
    12 | 01 | 02) yyyymm1_clim=$yyyymm1 ;; # If winter is sampled during the first year, even partially, then it is ok.
esac
yyyymm2_clim=$((yyyy2-1))11 # Default: we readjust the last year to November in the previous year.
case ${yyyymm2:4:2} in
    09 | 10 | 11) yyyymm2_clim=$yyyymm2 ;; # If fall is sampled during the last year, even partially, then it is ok.
    12) yyyymm2_clim=${yyyy2}11 ;; # If the last month is December, then we readjust it to November in the same year.
esac
if [ "$ref_grid" != "$config" ] ; then
    dir_level_by_level=$outdir_seasonal_clim/level_by_level/${yyyymm1_clim}_${yyyymm2_clim} ## Should already exist, no need to create it then.
else
    dir_level_by_level=$outdir_seasonal_clim/level_by_level/${yyyymm1_clim}_${yyyymm2_clim} ## Should already exist, no need to create it then.
fi

if [ "$bin_average_layers" = TRUE ] ; then
    list_seasons_def='Boreal'
    ### Output directory.
    outdir_layers_columns=$outdir_seasonal_clim/layers_columns ; mkdir -p $outdir_layers_columns
    echo -e "\e[00;32m First and last dates readjusted for interseasonal balance:\e[01;32m $yyyymm1_clim $yyyymm2_clim \e[00m"
    outdir_time_period=$outdir_layers_columns/${yyyymm1_clim}_${yyyymm2_clim} ; mkdir -p $outdir_time_period
    if [ "$ref_grid" != "$config" ] ; then
	outdir_columns_files=$outdir_time_period/resol_$ref_grid ; mkdir -p $outdir_columns_files
    else
	### Input directory for this section.
	outdir_columns_files=$outdir_time_period
    fi
    ### Computing the duration (in months) of the measurement period for each variable.
    listvar_nmonths=""
    for var in $listvar_IAGOS_output_name ; do
	yyyymm1_clim_var=$yyyymm1_clim ; yyyymm2_clim_var=$yyyymm2_clim
	case "$var" in
	    "$CO_varname_output_IAGOS" | "$O3toCO_varname_output_IAGOS")
		if [ $yyyymm1_clim -le 200112 ] ; then
		    yyyymm1_clim_var=200112
		else
		    yyyymm1_clim_var=$yyyymm1_clim
		fi
		;;
	    "$NOy_varname_output_IAGOS")
		if [ $yyyymm1_clim -le 200104 ] ; then
		    yyyymm1_clim_var=200104
		else
		    yyyymm1_clim_var=$yyyymm1_clim
		fi
		;;
	esac
	## Now counting the month between the 1st and the last month, for each IAGOS variable.
	yyyymm=$yyyymm1_clim_var ; nmonths_var=1
	while [ $yyyymm -lt $yyyymm2_clim_var ] ; do
	    yyyymm=`date -d "${yyyymm}01 1 months" +%Y%m`
	    nmonths_var=$((nmonths_var+1))
	done
	listvar_nmonths=$listvar_nmonths" "$nmonths_var
	if [ "$var" = "$O3_varname_output_IAGOS" ] ; then
	    nmonths_clim=$nmonths_var
	fi
    done
    list_seasons=''
    for seasons_def in $list_seasons_def ; do
	list_seasons=$list_seasons' '`define_seasons $seasons_def`
	if [ "$seasons_def" = Boreal ] ; then
	    list_seasons=$list_seasons' 'ANN
	fi
    done
    for season in $list_seasons ; do
    	for layer in $list_layers_climatologies ; do
	    input_file_IAGOS=$dir_level_by_level/IAGOS_${layer}_${season}_${yyyymm1_clim}_${yyyymm2_clim}.nc
	    input_file_model=$dir_level_by_level/${experiment}_${layer}_${season}_${yyyymm1_clim}_${yyyymm2_clim}.nc
	    tmpdir_there_clim=${tmpdir_there}_$layer
	    mkdir -p $tmpdir_there_clim
	    echo  -e "\e[00;31m | Now computing layered-seasonal climatologies. \e[00m"
	    exe_file=$exe_dir/exe_column_average_IAGOS_${id_job}
	    rm -f $exe_file
	    ifort -g -DasuX $pdir/global_var_interpol_mod.f90 $pdir/functions_mod.f90 $pdir/column_average_IAGOS.f90 $pdir/eumetsat_ropp/ropp_io/ncdf/ncdf.f90 -o $exe_file $NETCDF_INCDIR $NETCDF_LIBDIR $NETCDF_LIB -free
	    export FORT_FMT_NO_WRAP_MARGIN=true
	    compil_status=$?
	    if [ $compil_status -ne 0 ] ; then
		echo -e "\e[00;31m column_average_IAGOS.f90 compilation failed. Stops the program. \e[00m"
		# Then stops the chrono.
		end_run=`date +%s`
		time_spent=`convertsec $(($end_run-$start_run))`
		echo -e "\e[00;31m | END OF TRE PROGRAM \e[00m (duration: ${time_spent})"
		echo duration: $time_spent > $tmpdir/duration.txt
		exit
	    fi
	    echo  -e "\e[00;32m | $layer $season \e[00m"
	    cat > $pdir/namelist_columns <<EOF
$experiment
$yyyymm1_clim $yyyymm2_clim
$upper_level_UTLS_climato $lower_level_UTLS_climato
$upper_level_FT $lower_level_FT
$lon_varname_output $lat_varname_output $level_varname_output
$nlat_eff
$n_IAGOS_variables
$listvar_IAGOS_output_name
$listvar_nmonths
$nmonths_clim
$n_model_variables
$listvar_model_output_name
$N_obs_min_climato
'${tmpdir_there_clim}/'
$season
$layer
'${outdir_level_by_level}/'
'${outdir_columns_files}/'
'${input_file_IAGOS}'
'${input_file_model}'
EOF
	    # valgrind --leak-check=full --track-origins=yes $exe_file # One debug command.
	    $exe_file
	    exe_status=$?
	    wait & # Waiting for the routine to finish.
	    if [ $exe_status -eq 0 ] ; then # Stops the script if an error occurred during the routine execution.
		echo -e "\e[00;31m column_average_iagos.f90 terminated on the $model model grid \e[00m"
	    else
		echo -e "\e[00;31m Error during the column_average_iagos.f90 execution. \e[00m"
		exit
	    fi
	    ## Renaming the NOy variables into NOy_ACACIA variables, for the R post-treatment after.
	    if [ "$sci_program" == ACACIA ] ; then
		file_IAGOS=$outdir_columns_files/IAGOS_${layer}_${season}_${yyyymm1_clim}_${yyyymm2_clim}.nc
		file_model=$outdir_columns_files/${experiment}_${layer}_${season}_${yyyymm1_clim}_${yyyymm2_clim}.nc
		## Obs
		for stat in Mean SD Min Max N_diff_updown Pressure N_flights N_obs N_gridcells ; do
		    ncrename -v .NOy_${stat},NOy_ACACIA_${stat} $file_IAGOS
		done
		## Model
		for stat in Mean SD Min Max ; do
		    if [ "$model" == "INCA" ] ; then
			ncrename -v .NOy_${stat},NOy_ACACIA_${stat} $file_model
		    fi
		done
		if [ "$model" == "INCA" ] ; then ## Also adjusting the name of the model files into ACACIA as well, instead of Rescaled_Nudged.
		    mv $file_model $outdir_columns_files/ACACIA_${layer}_${season}_${yyyymm1_clim}_${yyyymm2_clim}.nc
		fi
	    fi
	done
    done
fi

if [ "$bin_zonal_cross_sections" = TRUE ] ; then
    ### Output directory.
    outdir_seasonal_clim=$global_outdir_resol/seasonal_climatologies ## Should already exist, no need to create it then.
    outdir_zonal_cross_sections=$outdir_seasonal_clim/zonal_cross_sections ; mkdir -p $outdir_zonal_cross_sections
    yyyymm1_clim=${yyyy1}12 # Default: we readjust the first year to December in the same year.
    case ${yyyymm1:4:2} in
	12 | 01 | 02) yyyymm1_clim=$yyyymm1 ;; # If winter is sampled during the first year, even partially, then it is ok.
    esac
    yyyymm2_clim=$((yyyy2-1))11 # Default: we readjust the last year to November in the previous year.
    case ${yyyymm2:4:2} in
	09 | 10 | 11) yyyymm2_clim=$yyyymm2 ;; # If fall is sampled during the last year, even partially, then it is ok.
	12) yyyymm2_clim=${yyyy2}11 ;; # If the last month is December, then we readjust it to November in the same year.
    esac
    echo -e "\e[00;32m First and last dates readjusted for interseasonal balance:\e[01;32m $yyyymm1_clim $yyyymm2_clim \e[00m"
    outdir_time_period=$outdir_zonal_cross_sections/${yyyymm1_clim}_${yyyymm2_clim} ; mkdir -p $outdir_time_period
    if [ "$ref_grid" != "$config" ] ; then
	outdir_zonal_files=$outdir_time_period/resol_$ref_grid ; mkdir -p $outdir_zonal_files
    else
	outdir_zonal_files=$outdir_time_period
    fi
    ### Computing the duration (in months) of the measurement period for each variable.
    listvar_nmonths=""
    for var in $listvar_IAGOS_output_name ; do
	yyyymm1_clim_var=$yyyymm1_clim ; yyyymm2_clim_var=$yyyymm2_clim
	case $var in
	    $CO_varname_output_IAGOS | $O3toCO_varname_output_IAGOS)
		if [ $yyyymm1_clim -le 200112 ] ; then
		    yyyymm1_clim_var=200112
		else
		    yyyymm1_clim_var=$yyyymm1_clim
		fi
		;;
	    $H2O_varname_output_IAGOS | $RHL_varname_output_IAGOS | $RHI_varname_output_IAGOS)
		if [ $yyyymm2_clim -ge 201312 ] ; then
		    yyyymm2_clim_var=201312
		else
		    yyyymm2_clim_var=$yyyymm2_clim
		fi
		;;
	    $NOy_varname_output_IAGOS)
		if [ $yyyymm1_clim -le 200104 ] ; then
		    yyyymm1_clim_var=200104
		else
		    yyyymm1_clim_var=$yyyymm1_clim
		fi
		;;
	esac
	## Now counting the month between the 1st and the last month, for each IAGOS variable.
	yyyymm=$yyyymm1_clim_var ; nmonths_var=1
	while [ $yyyymm -lt $yyyymm2_clim_var ] ; do
	    yyyymm=`date -d "${yyyymm}01 1 months" +%Y%m`
	    nmonths_var=$((nmonths_var+1))
	done
	listvar_nmonths=${listvar_nmonths}" "$nmonths_var
	if [ "$var" = "$O3_varname_output_IAGOS" ] ; then
	    nmonths_clim=$nmonths_var
	fi
    done
    for seasons_def in $list_seasons_def_tropics ; do
	## First finds out the coordinates for the zonal band.
	case $seasons_def in
	    'Tropical_Atl' )
		lon01=-60 ; lon02=-15 ;; ## As in Sauvage et al., 2007, GRL
	    'Tropical_Afr' )
		lon01=-5 ; lon02=30 ;; ## As in Sauvage et al., 2007, GRL
	    'South_Asia' )
		lon01=60 ; lon02=90 ;;
	esac
	## Then finds out the seasons definitions.
	list_seasons=`define_seasons $seasons_def`
	for season in $list_seasons ; do
    	    for layer in $list_layers_zonal_cross_sections ; do
		# for season in DJF ; do
		# 	for layer in UTLS ; do
		input_file_IAGOS=$dir_level_by_level/IAGOS_${layer}_${season}_${yyyymm1_clim}_${yyyymm2_clim}.nc
		input_file_model=$dir_level_by_level/${experiment}_${layer}_${season}_${yyyymm1_clim}_${yyyymm2_clim}.nc
		tmpdir_there_clim=${tmpdir_there}_$layer
		mkdir -p $tmpdir_there_clim
		echo  -e "\e[00;31m | Now computing mean zonal cross sections. \e[00m" 
		exe_file=$exe_dir/exe_zonal_${id_job}
		rm -f $exe_file
		ifort -g -DasuX $pdir/global_var_interpol_mod.f90 $pdir/functions_mod.f90 $pdir/zonal_cross_sections_IAGOS.f90 $pdir/eumetsat_ropp/ropp_io/ncdf/ncdf.f90 -o $exe_file $NETCDF_INCDIR $NETCDF_LIBDIR $NETCDF_LIB -free
		export FORT_FMT_NO_WRAP_MARGIN=true
		compil_status=$?
		if [ $compil_status -ne 0 ] ; then
		    echo -e "\e[00;31m zonal_cross_sections_IAGOS.f90 compilation failed. Stops the program. \e[00m"
		    # Then stops the chrono.
		    end_run=`date +%s`
		    time_spent=`convertsec $(($end_run-$start_run))`
		    echo -e "\e[00;31m | END OF TRE PROGRAM \e[00m (duration: ${time_spent})"
		    echo duration: $time_spent > $tmpdir/duration.txt
		    exit
		fi
		echo  -e "\e[00;32m | $layer $season \e[00m"  
		cat > $pdir/namelist_zonal <<EOF
$experiment $model
$yyyymm1_clim $yyyymm2_clim
$upper_level_UTLS_zonal $lower_level_UTLS_zonal
$upper_level_FT $lower_level_FT
$lon01 $lon02
$lon_varname_output $lat_varname_output $level_varname_output
$nlat_eff
$n_IAGOS_variables
$listvar_IAGOS_output_name
$listvar_nmonths
$nmonths_clim
$n_model_variables
$listvar_model_output_name
$N_obs_min_climato
'${tmpdir_there_clim}/'
$season
$layer
'${outdir_level_by_level}/'
'${outdir_zonal_files}/'
'${input_file_IAGOS}'
'${input_file_model}'
EOF
		# valgrind --leak-check=full --track-origins=yes $exe_file # One debug command.
		$exe_file
		exe_status=$?
		wait & # Waiting for the routine to finish.
		if [ $exe_status -eq 0 ] ; then # Stops the script if an error occurred during the routine execution.
		    echo -e "\e[00;31m zonal_cross_sections_iagos.f90 terminated on the $model model grid \e[00m"
		else
		    echo -e "\e[00;31m Error during the zonal_cross_sections_iagos.f90 execution. \e[00m"
		    exit
		fi
		## Renaming the NOy variables into NOy_ACACIA variables, for the R post-treatment after.
		if [ "$sci_program" == ACACIA ] ; then
		    file_IAGOS=$outdir_zonal_files/IAGOS_${layer}_${season}_${lon01}E_${lon02}E_${yyyymm1_clim}_${yyyymm2_clim}.nc
		    file_model=$outdir_zonal_files/${experiment}_${layer}_${season}_${lon01}E_${lon02}E_${yyyymm1_clim}_${yyyymm2_clim}.nc
		    ## Obs
		    for stat in Mean SD Min Max N_diff_updown Pressure N_flights N_obs N_gridcells Decile_1 Quartile_1 Quartile_3 Decile_9 ; do
			ncrename -v .NOy_${stat},NOy_ACACIA_${stat} $file_IAGOS
		    done
		    ## Model
		    for stat in Mean SD Min Max Decile_1 Quartile_1 Quartile_3 Decile_9 ; do
			if [ "$model" == "INCA" ] ; then
			    ncrename -v .NOy_${stat},NOy_ACACIA_${stat} $file_model
			fi
		    done
		    if [ "$model" == "INCA" ] ; then ## Also adjusting the name of the model files into ACACIA as well, instead of Rescaled_Nudged.
			mv $file_model $outdir_zonal_files/ACACIA_${layer}_${season}_${lon01}E_${lon02}E_${yyyymm1_clim}_${yyyymm2_clim}.nc
		    fi
		fi
	    done
	done
    done
fi
if [ "$bin_colissimo" == TRUE ] ; then
    outdir_easy=$CCCSCRATCHDIR/IAGOS_$sci_program/$model
    mkdir -p $CCCSCRATCHDIR/IAGOS_$sci_program
    mkdir -p $outdir_easy
    if [ "$bin_average_layers" == TRUE ] ; then
	 cp -r $outdir_layers_columns $outdir_easy
    fi
    if [ "$bin_zonal_cross_sections" == TRUE ] ; then
	cp -r $outdir_zonal_cross_sections $outdir_easy
    fi
fi

# Here stops the chrono.
end_run=`date +%s`
time_spent=`convertsec $(($end_run-$start_run))`
echo -e "\e[00;31m | END OF TRE PROGRAM \e[00m (duration: ${time_spent})"
echo duration: $time_spent > $tmpdir/duration.txt

if [ "$shutdown_when_run_completed" = "yes" ] ; then
    shutdown
else
    exit 0
fi
