! This program aims at converting the IAGOS NetCDF files downloaded from
! the IAGOS platform into monthly gridded NetCDF files.
! It can be used as a blackbox: you do not have anything to change here.
! The parameters that you chose to process must have been specified in
! the main shell script make_interpol_iagos.sh.
! The current version is used in order to make the IAGOS data directly
! comparable to the chosen model simulation outputs.

PROGRAM INTERPOL_IAGOS
  USE NETCDF
  USE GLOBAL_VAR_INTERPOL_MOD
  USE FUNCTIONS_MOD
  IMPLICIT NONE
  ! Checks the existence of the pressure files. If they do not exist, then we cannot do anything.
  LOGICAL pres_file_exists
  ! Depending on the layers asked by the user, logic_PV decides whether or not the program reads the PV fields.
  ! logic_P decided whether the pressure is read from the model 3D field or calculated as sigma-P hybrid levels,
  ! thus only depending on the surface pressure field and the vertical coefficients (A_i) and (B_i).
  LOGICAL logic_PV, logic_P, logic_cst_P
  ! The parts of the IAGOS filename preceeding and following the flight number.
  CHARACTER(LEN=50) filename_prefix_IAGOS, filename_suffix_IAGOS
  ! Index pointing to a given IAGOS file in each monthly list of files.
  INTEGER iFile
  ! Dates as integers, then as characters.
  ! Their lengths have to be the same as the variable in INT2CHAR subroutine (here: 20).
  INTEGER iyyyymm, iyyyy, iyyyymm1, iyyyymm2
  CHARACTER(LEN=30) cyyyymm, cyyyymm1, cyyyy, cnlev, cnlat, cnlon
  ! Period length of the time series, in months.
  INTEGER nmonths
  ! Indexes repairing in the IAGOS files, and in the lists of variables and dimensions (also for the model outputs). ilay is the index for the atmospheric layer. istat is the index for the metrics (mean value, SD,...).
  INTEGER icol_unit_IAGOS_input, icol_iagos_output, icol_unit_iagos_output
  INTEGER ivarmod, imodel, icol_model, icol_model_output, icol_unit_model_output, icol_model_full_name
  INTEGER icol_program
  INTEGER idimmod_time, idimmod_lon, idimmod_lat, idimmod_lev
  INTEGER istat, ilay
  ! Integers relative to the model grid dimension.
  ! Amount of days during the current month. Set at 1 if the time resolution is monthly.
  INTEGER ntime
  ! l,k,i and imth are used respectively for the height, latitude, longitude and time dimensions.
  ! itime is generally for the day into the current month.
  ! imth_this_data should be created to account for flights that take place at nighttime, at the end of the current month.
  INTEGER l, k, i, imth, itime
  ! Reals corresponding to the domain borders and to the grid horizontal resolution.
  REAL lat1, lat2, lon1, lon2, dlat, dlon
  ! The following shift_lon and shift_lat are a correction for the grid,
  ! in case the gridcells are defined by their center coordinates (as for the INCA model) instead of their
  ! southwestern border (as for MOCAGE).
  REAL shift_lon, shift_lat
  ! Tropopause definition using PV for extratropics.
  REAL PV_tropopause, P_0
!!!!!!! Dynamic arrays !!!!!!!!
  ! Averaging step: indexes repairing the gridcells validated for being averaged.
  INTEGER, ALLOCATABLE :: ilon_ok(:), klat_ok(:), llev_ok(:)
  ! Array stocking the modelled PV field. Important for the dynamical UT/LS separation.
  REAL, ALLOCATABLE :: PV(:,:,:,:)
  ! List of the IAGOS data files for the current month.
  CHARACTER*200,ALLOCATABLE :: list_data_files(:)
  ! list_months is the monthly time dimension. It contains every month
  ! in the period iyyyymm1 - iyyyymm2, in the yyyymm format.
  ! Based on the same dimension, list_nfiles tells how many IAGOS files are contained in each month.
  INTEGER,ALLOCATABLE :: list_months(:), list_nfiles(:)
  ! Lists of dimensions and variables. The suffix _nc is for the variable names as outputs from this routine.
  CHARACTER*30,ALLOCATABLE :: list_coord(:), list_variables(:), list_variables_nc(:)
!!!! Temporary lists, quickly deallocated. They are used only to make the code organization clearer.
  ! list_obs_variables(_nc) refers to the list of observed variables required by the user.
  CHARACTER*30,ALLOCATABLE :: list_obs_variables(:)
  CHARACTER*30,ALLOCATABLE :: list_obs_variables_nc(:) ! The names given to the variables in the output NetCDF file.
  ! list_der_variables is the list of variables (asked by the user) derived from the observed ones.
  CHARACTER*30,ALLOCATABLE :: list_der_variables(:)
  CHARACTER*30,ALLOCATABLE :: list_der_variables_nc(:)
!!!! Permanent lists.
  ! Lists of variables and dimensions as queried by the user (standard name).
  CHARACTER*30,ALLOCATABLE :: list_dim_model(:), list_var_model(:)
  ! Lists of variables and dimensions as named in the model output files (i.e. as inputs for this routine).
  CHARACTER*30,ALLOCATABLE :: model_dimensions(:), model_variables(:), model_var_output(:)
  ! Pressure, excluded from the model lists cause it is not a variable we want in output for this routine.
  CHARACTER*30                model_P, model_P_S
  ! Same list, gathering the variables and all their metrics in the format: variable_stat. Size: nvar_model*nstats_model.
  CHARACTER*30,ALLOCATABLE :: list_var_model_with_stats(:)
  ! Corresponding units in the input and output model files. Below, corresponding input-output conversion factors.
  CHARACTER*30,ALLOCATABLE :: model_units_output(:), model_units_dim_output(:), model_units_var_output(:)
  REAL,ALLOCATABLE :: model_conv_factor(:)
  ! Also for units in IAGOS files.
  CHARACTER*30,ALLOCATABLE :: IAGOS_units_output(:), IAGOS_units_input(:), IAGOS_units_output_with_stats(:)
  REAL,ALLOCATABLE :: IAGOS_conv_factor(:)
  ! list_variables_with_stats gathers both lists of variables and their respective statistics (mean, SD, min, max, data amounts,...).
  ! It corresponds to the names given to all the fields saved in the NetCDF files.
  CHARACTER*30,ALLOCATABLE :: list_variables_with_stats(:)
  ! list_var_tot is only used for checking that the program treats the variables correctly.
  ! It gathers every IAGOS column, in the following order: dimensions, observed variables and derived variables.
  CHARACTER*30,ALLOCATABLE :: list_var_tot(:)
  ! Names of the different IAGOS packages, and more generally the header of the table gathering the variable names under these IAGOS packages.
  CHARACTER*30,ALLOCATABLE :: package_names(:), header_packages(:)
  ! Packages asked in the make file.
  CHARACTER*30,ALLOCATABLE :: list_packages(:)
  ! Permanent tabular gathering all the possible names for each required variable.
  CHARACTER*30,ALLOCATABLE :: tab_var_names(:,:)
  ! Names of the different models which variables have been registered into a reference tabular.
  ! The size of model_names is n_models_registered. +2 for header models.
  CHARACTER*30,ALLOCATABLE :: model_names(:), header_models(:)
  ! Headers for the two other tables, gathering the names of the units in input (depending on the data set)
  ! and the units that are to be given to the outputs. Consequently, they also contain the appropriate
  ! conversion factors.
  CHARACTER*30,ALLOCATABLE :: header_units_IAGOS(:), header_units_models(:)
  ! Permanent tabular gathering all the possible names for each required variable.
  CHARACTER*30,ALLOCATABLE :: tab_model_names(:,:)
  ! Permanent tabulars gathering the units in input and output files in IAGOS and models, 
  ! with the corresponding conversion factors.
  CHARACTER*30,ALLOCATABLE :: tab_var_units_IAGOS(:,:), tab_var_units_models(:,:)
  ! Threshold values for each variable. Below it, the data point is filtered out.
  REAL,ALLOCATABLE :: thres_val(:)
!!!!!!!!! Layers properties. !!!!!!!!
  ! Names.
  CHARACTER*30,ALLOCATABLE :: list_layers(:)
  ! Boundaries, in matter of potential vorticity (PV), pressures and ozone volume mixing ratio.
  ! The last one is used to filter out the observations in the stratospheric air masses detected
  ! by the model PV field as tropospheric, and vice-versa.
  ! The pressures are given in Pa.
  REAL,ALLOCATABLE :: PV_down_layers(:)  , PV_up_layers(:)
  REAL,ALLOCATABLE :: P_down_layers(:)   , P_up_layers(:)
  REAL,ALLOCATABLE :: ozone_max_layers(:), ozone_min_layers(:)
!!!!!!!!! Gridded data. !!!!!!!!!!
  ! Gridded dimensions. The _netcdf suffix is for the dimension as it will be written in this routine output.
  REAL,ALLOCATABLE :: rlat(:), rlon(:), rlat_netcdf(:), rlon_netcdf(:)
  ! Arrays needed to build the pressure grid. A and B coefficients, and surface pressure for the current month.
  ! The dimensions for P_surf are longitude, latitude and the day in the current month.
  REAL,ALLOCATABLE :: A_inter(:), B_inter(:), P_surf(:,:,:)
  ! Half and full vertical grid levels.
  REAL,ALLOCATABLE :: P_edge(:,:,:,:), P_full(:,:,:,:)
  REAL,ALLOCATABLE :: P_mid(:,:,:,:)
  ! The next array is, if precised, the predefined pressure levels. In this case, P_mid will have their values.
  REAL,ALLOCATABLE :: cst_P_levels(:)
  ! 5-dimension arrays, for observations (monthly or daily) averages and their corresponding amounts of data.
  ! rcount_daily is the amount of data with a daily resolution. It allows to apply a IAGOS mask on the model
  ! fields during the monthly averaging step, in case the model outputs have a daily resolution.
  REAL,ALLOCATABLE :: mean_5D(:,:,:,:,:), rcount(:,:,:,:,:), rcount_daily(:,:,:,:,:,:), P_5D(:,:,:,:,:)
  ! N_flights_5D counts the amount of flights involved in the monthly mean.
  REAL,ALLOCATABLE :: N_flights_5D(:,:,:,:,:)
  ! Same for the minimum and maximum observed values.
  REAL,ALLOCATABLE :: min_5D(:,:,:,:,:), max_5D(:,:,:,:,:)
  ! Same, for standard deviation (SD).
  REAL,ALLOCATABLE :: SD_5D(:,:,:,:,:)
  ! day1_5D and dayN_5D are the first and last days of valid measurements during the current month.
  ! Both spread between 1 and 31.
  INTEGER,ALLOCATABLE :: day1_5D(:,:,:,:,:), dayN_5D(:,:,:,:,:)
  ! n_diff_updown represents the weights summed up during each month (like rcount, but with positive or 
  ! negative values depending on the sign of the pressure difference between the obs and the vertical grid level.
  ! It is normalized to rcount, thus it spreads from -1 up to 1. For a given grid cell, 
  ! the more n_diff_updown is close to 0, the more the sampling is representative of the grid cell center
  ! (on the vertical axis).
  REAL,ALLOCATABLE :: n_diff_updown(:,:,:,:,:)
  ! The previous arrays (mean, SD, min, max, rcount, day1, dayN) gathered into one, to be written in the NetCDF output files.
  REAL,ALLOCATABLE :: mean_5D_with_stats(:,:,:,:,:)
  ! ############## Interpolation variables ###############
  ! 5-dimension arrays for model output fields (var_model_daily) and for layered model fields (mean_5D_model).
  ! Their dimensions are not the same. The first one is daily resolved, the second one is layer-resolved and
  ! is composed by monthly means.
  REAL,ALLOCATABLE :: var_model_daily(:,:,:,:,:), mean_5D_model(:,:,:,:,:), mean_5D_model_with_stats(:,:,:,:,:)
  ! Same for min and max.
  REAL,ALLOCATABLE :: min_5D_model(:,:,:,:,:), max_5D_model(:,:,:,:,:)
  ! Same, for SD.
  REAL,ALLOCATABLE :: SD_5D_model(:,:,:,:,:)

  tmpdir=TRIM(pdir)
  CALL FILL_STR_CHAR(tmpdir,6,LEN(tmpdir),tmpdir)
  ! Reads all the input parameters asked in make_interpol_iagos.sh
  OPEN(14, FILE=TRIM(tmpdir)//'namelist_extract', FORM='formatted')
  READ(14,*) nlev
  WRITE(*,*) 'nlev=', nlev
  READ(14,*) dlat, dlon
  WRITE(*,*) 'dlat, dlon=', dlat, dlon
  READ(14,*) nlat, nlon
  WRITE(*,*) 'nlat, nlon=', nlat, nlon
  READ(14,*) lat1, lat2
  WRITE(*,*) 'lat1, lat2=', lat1, lat2
  READ(14,*) lon1, lon2
  WRITE(*,*) 'lon1, lon2=', lon1, lon2
  READ(14,*) iyyyymm1, iyyyymm2
  WRITE(*,*) 'iyyyymm1, iyyyymm2=', iyyyymm1, iyyyymm2
  READ(14,*) nmonths
  WRITE(*,*) 'nmonths=', nmonths
  ALLOCATE(list_months(nmonths))
  READ(14,*) list_months(1:nmonths)
  WRITE(*,*) 'list_months(1:nmonths)=',list_months(1:nmonths)
  ALLOCATE(list_nfiles(nmonths))
  READ(14,*) list_nfiles(1:nmonths)
  WRITE(*,*) 'list_nfiles(1:nmonths)=',list_nfiles(1:nmonths)
  READ(14,*) ncoord
  WRITE(*,*) 'ncoord=',ncoord
  ALLOCATE(list_coord(ncoord))
  READ(14,*) list_coord(1:ncoord)
  WRITE(*,*) 'list_coord(1:ncoord)=',list_coord(1:ncoord)
  READ(14,*) n_obs_variables
  WRITE(*,*) 'n_obs_variables=',n_obs_variables
  ALLOCATE(list_obs_variables(n_obs_variables)) ! Variables in their full name.
  ALLOCATE(list_obs_variables_nc(n_obs_variables))
  READ(14,*) list_obs_variables(1:n_obs_variables)
  WRITE(*,*) 'list_obs_variables(1:n_obs_variables)=',list_obs_variables(1:n_obs_variables)
  READ(14,*) list_obs_variables_nc(1:n_obs_variables)
  WRITE(*,*) 'list_obs_variables_nc(1:n_obs_variables)=',list_obs_variables_nc(1:n_obs_variables)
  READ(14,*) n_der_variables
  WRITE(*,*) 'n_der_variables=',n_der_variables
  ALLOCATE(list_der_variables(n_der_variables)) ! Variables not in the IAGOS files
  ALLOCATE(list_der_variables_nc(n_der_variables))
  READ(14,*) list_der_variables(1:n_der_variables)! (it means 'derived variables')
  WRITE(*,*) 'list_der_variables(1:n_der_variables)=',list_der_variables(1:n_der_variables)
  READ(14,*) list_der_variables_nc(1:n_der_variables)
  WRITE(*,*) 'list_der_variables_nc(1:n_der_variables)=',list_der_variables_nc(1:n_der_variables)
  READ(14,*) logic_PV
  WRITE(*,*) 'logic_PV=',logic_PV
  READ(14,*) PV_tropopause
  WRITE(*,*) 'PV_tropopause=',PV_tropopause
  READ(14,*) logic_P
  WRITE(*,*) 'logic_P=',logic_P
  READ(14,*) P_0
  WRITE(*,*) 'P_0=', P_0
  READ(14,*) logic_cst_P
  WRITE(*,*) 'logic_cst_P=',logic_cst_P
  READ(14,*) nlayers
  WRITE(*,*) 'nlayers=',nlayers
  ALLOCATE(list_layers(nlayers))
  ALLOCATE(PV_down_layers(nlayers))
  ALLOCATE(PV_up_layers(nlayers))
  ALLOCATE(P_down_layers(nlayers))
  ALLOCATE(P_up_layers(nlayers))
  ALLOCATE(ozone_min_layers(nlayers))
  ALLOCATE(ozone_max_layers(nlayers))
  READ(14,*) list_layers(1:nlayers)
  READ(14,*) PV_down_layers(1:nlayers)
  READ(14,*) PV_up_layers(1:nlayers)
  READ(14,*) P_down_layers(1:nlayers)
  READ(14,*) P_up_layers(1:nlayers)
  READ(14,*) ozone_min_layers(1:nlayers)
  READ(14,*) ozone_max_layers(1:nlayers)
  READ(14,*) ndim_model
  ALLOCATE(list_dim_model(ndim_model))
  READ(14,*) list_dim_model(1:ndim_model)
  READ(14,*) nvar_model
  ALLOCATE(list_var_model(nvar_model))
  READ(14,*) list_var_model(1:nvar_model)
  READ(14,*) n_models_registered
  READ(14,*) ref_dir
  READ(14,*) model_output_dir
  READ(14,*) init_model_output_dir
  READ(14,*) data_dir_cst
  READ(14,*) tmpdir
  READ(14,*) tmpdir_there
  READ(14,*) filename_cst_P_levels
  OPEN(unitFile_check, FILE=TRIM(tmpdir_there)//'check_interpol_iagos',   &
       FORM='formatted')
  WRITE(*,*) 'list_layers(1:nlayers)=',list_layers(1:nlayers)
  WRITE(*,*) 'PV_down_layers(1:nlayers)=',PV_down_layers(1:nlayers)
  WRITE(*,*) 'PV_up_layers(1:nlayers)=',PV_up_layers(1:nlayers)
  WRITE(*,*) 'P_down_layers(1:nlayers)=',P_down_layers(1:nlayers)
  WRITE(*,*) 'P_up_layers(1:nlayers)=',P_up_layers(1:nlayers)
  WRITE(*,*) 'P_up_layers(1:nlayers)=',P_up_layers(1:nlayers)
  WRITE(*,*) 'ozone_min_layers(1:nlayers)=',ozone_min_layers(1:nlayers)
  WRITE(*,*) 'ndim_model=',ndim_model
  WRITE(*,*) 'ozone_max_layers(1:nlayers)=',ozone_max_layers(1:nlayers)
  WRITE(*,*) 'list_dim_model(1:ndim_model)=',list_dim_model(1:ndim_model)
  WRITE(*,*) 'nvar_model=',nvar_model
  WRITE(*,*) 'list_var_model(1:nvar_model)=',list_var_model(1:nvar_model)
  WRITE(*,*) 'n_models_registered=',n_models_registered
  WRITE(*,*) 'ref_dir=',ref_dir
  WRITE(*,*) 'model_output_dir=',model_output_dir
  WRITE(*,*) 'init_model_output_dir=',init_model_output_dir
  WRITE(*,*) 'data_dir_cst=',data_dir_cst
  WRITE(*,*) 'tmpdir=',tmpdir
  WRITE(*,*) 'tmpdir_there=',tmpdir_there
  WRITE(*,*) 'filename_cst_P_levels=',TRIM(filename_cst_P_levels)
  READ(14,*) outdir_NC
  WRITE(*,*) 'outdir_NC=',outdir_NC
  READ(14,*) cmodel
  WRITE(*,*) 'cmodel=',cmodel
  READ(14,*) experiment
  WRITE(*,*) 'experiment=',experiment
  READ(14,*) sci_program
  WRITE(*,*) 'sci_program=',sci_program
  READ(14,*) time_resol
  WRITE(*,*) 'time_resol=',time_resol
  READ(14,*) filename_prefix_IAGOS
  WRITE(*,*) 'filename_prefix_IAGOS=',filename_prefix_IAGOS
  READ(14,*) filename_suffix_IAGOS
  WRITE(*,*) 'filename_suffix_IAGOS=',filename_suffix_IAGOS
  READ(14,*) quicktest
  WRITE(*,*) 'quicktest=',quicktest
  READ(14,*) npackages
  WRITE(*,*) 'npackages=',npackages
  ALLOCATE(list_packages(npackages))
  READ(14,*) list_packages
  WRITE(*,*) 'list_packages=',list_packages
  CLOSE(14)
  nvariables= n_der_variables + n_obs_variables
  WRITE(*,*) 'nvariables=',nvariables
  nvar_tot=nvariables + ncoord
  nvariables_with_stats = nvariables*nstats_obs
  nvar_iagos = n_obs_variables + ncoord ! The amount of columns picked in one IAGOS file.
  ALLOCATE(list_variables(nvariables),list_variables_nc(nvariables))
  ALLOCATE(list_var_tot(nvar_tot))
  DO ivar=1,ncoord
     list_var_tot(ivar)=list_coord(ivar)
  END DO
  DO ivar=1,n_obs_variables
     list_variables(ivar)=list_obs_variables(ivar)
     list_variables_nc(ivar)=list_obs_variables_nc(ivar)
     list_var_tot(ivar+ncoord)=list_variables(ivar)
  END DO
  DEALLOCATE(list_obs_variables)
  DEALLOCATE(list_obs_variables_nc)
  DO ivar=1,n_der_variables
     list_variables(ivar+n_obs_variables)=list_der_variables(ivar)
     list_variables_nc(ivar+n_obs_variables)=list_der_variables_nc(ivar)
     list_var_tot(ivar+n_obs_variables+ncoord)= &
          list_variables(ivar+n_obs_variables)
  END DO
  ! Converting the dimension sizes into characters.
  call INT2CHAR(nlev,tmpdir,cnlev)
  call INT2CHAR(nlat,tmpdir,cnlat)
  call INT2CHAR(nlon,tmpdir,cnlon)
  call INT2CHAR(iyyyymm1,tmpdir,cyyyymm1) ! Just in order to keep the time origin in memory, for the NetCDF time coordinate
  ALLOCATE(A_inter(nlev),B_inter(nlev))
  ! A difference must be made between the grids that define a gridcell by its center and by its southwestern border.
  IF (TRIM(cmodel).eq.'MOCAGE') THEN
     shift_lon=0
     shift_lat=0
  ELSE
     ! Modifying the gridded longitudes to target the correct grid cells.
     ! Note that we need to be very careful when crossing the 180E meridian.
     shift_lon=-dlon/2
     shift_lat=-dlat/2
  END IF
  ! shift_lon=0
  ! shift_lat=0
  ALLOCATE(rlon(nlon), rlat(nlat), rlon_netcdf(nlon), rlat_netcdf(nlat))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Assigning individual names to coordinates indexes.
!!! For example, icolUTC must correspond to the index of the UTC_time dimension in the list.
  ! Time_UTC
  idimUTC= which_one_str(c_Time_UTC, list_coord, ncoord)
  WRITE(unitFile_check,*) 'Time_UTC:::::::::', list_coord(idimUTC)
  ! Latitude
  idimlat= which_one_str(c_Latitude, list_coord, ncoord)
  WRITE(unitFile_check,*) 'Latitude:::::::::', list_coord(idimlat)
  ! Longitude
  idimlon= which_one_str(c_Longitude, list_coord, ncoord)
  WRITE(unitFile_check,*) 'Longitude:::::::::', list_coord(idimlon)
  ! Baro_altitude
  idimalt= which_one_str(c_Baro_altitude, list_coord, ncoord)
  WRITE(unitFile_check,*) 'Baro_altitude:::::::::', list_coord(idimalt)
  ! Pressure
  idimP= which_one_str(c_Pressure, list_coord, ncoord)
  WRITE(unitFile_check,*) 'Pressure:::::::::', list_coord(idimP)
  ! Ozone (we expect the ozone column to be contained in every file)
  IF (any(list_variables=='Ozone')) THEN
     ivarO3= which_one_str(c_Ozone, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Ozone:::::::::', list_variables(ivarO3)
  END IF
  ! Carbon_monoxide
  IF (any(list_variables=='Carbon_monoxide')) THEN
     ivarCO= which_one_str(c_Carbon_monoxide, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Carbon_monoxide:::::::::', list_variables(ivarCO)
  END IF
  ! Water_vapour
  IF (any(list_variables=='Water_vapour')) THEN
     ivarH2O= which_one_str(c_Water_vapour, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Water_vapour:::::::::', list_variables(ivarH2O)
  END IF  
  ! Relative_humidity_liquid
  IF (any(list_variables=='Relative_humidity_liquid')) THEN
     ivarRHL= which_one_str(c_Relative_humidity_liquid, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Relative_humidity_liquid:::::::::', list_variables(ivarRHL)
  END IF
  ! Relative_humidity_ice
  IF (any(list_variables=='Relative_humidity_ice')) THEN
     ivarRHI= which_one_str(c_Relative_humidity_ice, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Relative_humidity_ice:::::::::', list_variables(ivarRHI)
  END IF
  ! Cloud_liquid_water_particles
  IF (any(list_variables=='Cloud_liquid_water_particles')) THEN
     ivarCloud= which_one_str(c_Cloud_liquid_water_particles, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Cloud_liquid_water_particles:::::::::', list_variables(ivarCloud)
  END IF
  ! Nitrogen_reactive_species
  IF (any(list_variables=='Nitrogen_reactive_species')) THEN
     ivarNOy= which_one_str(c_Nitrogen_reactive_species, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Nitrogen_reactive_species:::::::::', list_variables(ivarNOy)
  END IF  
  ! Nitrogen_oxides
  IF (any(list_variables=='Nitrogen_oxides')) THEN
     ivarNOx= which_one_str(c_Nitrogen_oxides, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Nitrogen_oxides:::::::::', list_variables(ivarNOx)
  END IF
  ! Nitrogen_monoxide
  IF (any(list_variables=='Nitrogen_monoxide')) THEN
     ivarNO= which_one_str(c_Nitrogen_monoxide, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Nitrogen_monoxide:::::::::', list_variables(ivarNO)
  END IF
  ! Nitrogen_dioxide
  IF (any(list_variables=='Nitrogen_dioxide')) THEN
     ivarNO2= which_one_str(c_Nitrogen_dioxide, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Nitrogen_dioxide:::::::::', list_variables(ivarNO2)
  END IF
  ! Pressure (as a variable)
  IF (any(list_variables=='Pressure')) THEN
     ivarP= which_one_str(c_Pressure, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Pressure:::::::::', list_variables(ivarP)
  END IF
  ! Zonal wind speed (in m.s-1)
  IF (any(list_variables=='Zonal_wind_speed')) THEN
     ivarU= which_one_str(c_Zonal_wind_speed, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Zonal_wind_speed:::::::::', list_variables(ivarU)
  END IF
  ! Meridional wind speed (in m.s-1)
  IF (any(list_variables=='Meridional_wind_speed')) THEN
     ivarV= which_one_str(c_Meridional_wind_speed, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Meridional_wind_speed:::::::::', list_variables(ivarV)
  END IF
  ! Temperature (in Celsius degree)
  IF (any(list_variables=='Temperature')) THEN
     ivarT= which_one_str(c_Temperature, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Temperature:::::::::', list_variables(ivarT)
  END IF
  ! Theta (in Kelvin degree)
  IF (any(list_variables=='Potential_temperature')) THEN
     ivarTheta= which_one_str(c_Potential_temperature, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Potential_temperature:::::::::', list_variables(ivarTheta)
  END IF
  ! O3-to-CO ratio
  IF (any(list_variables=='Ozone_to_CO_ratio')) THEN
     ivarO3toCO= which_one_str(c_Ozone_to_CO_ratio, list_variables, nvariables)
     WRITE(unitFile_check,*) 'Ozone_to_CO_ratio:::::::::', list_variables(ivarO3toCO)
 END IF
 ! Model water vapour
  IF (any(list_var_model=='Water_vapour')) THEN
     ivarmod_H2O= which_one_str(c_Water_vapour, list_var_model, nvar_model)
     WRITE(unitFile_check,*) 'Model water vapour:::::::::', list_var_model(ivarmod_H2O)
 END IF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Filtering minimum values
  ALLOCATE(thres_val(nvariables))
  DO ivar=1,nvariables
     thres_val(ivar)=0
  END DO
  ! Exceptions: T (temperature in Celsius degree), and zonal and meridional wind speeds.
  thres_val(ivarT)=-97
  thres_val(ivarU)=-97
  thres_val(ivarV)=-97
!!!!!!!!!!! Reading the 4 reference tables (IAGOS and model variable names and their respective units),
!!!!!!!!!!! with their corresponding conversion factors.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IAGOS variable names.
  WRITE(unitFile_check,*) '############# IAGOS variable names ################'
  var_names_file=TRIM(tmpdir)//'table_variable_names_in_IAGOS_files.txt'
  OPEN(14,FILE=TRIM(var_names_file),FORM='formatted')
  ALLOCATE(header_packages(npackages+2))
  ALLOCATE(package_names(npackages))
  READ(14,*) header_packages
  WRITE(unitFile_check,*) header_packages
  ALLOCATE(tab_var_names(nvar_tot,npackages+2))
  DO ivar=1,nvar_tot
     READ(14,*) tab_var_names(ivar,:)
     WRITE(unitFile_check,*) tab_var_names(ivar,:)
  END DO
  package_names=header_packages(2:npackages+1)
  CLOSE(14)
  WRITE(unitFile_check,*) 'package_names=', package_names
  WRITE(unitFile_check,*) '###################################################'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! IAGOS variable units.
  WRITE(unitFile_check,*) '############# IAGOS variable units ################'
  var_units_file=TRIM(tmpdir)//'table_variable_units_in_IAGOS_files.txt'
  OPEN(14,FILE=TRIM(var_units_file),FORM='formatted')
  ALLOCATE(header_units_IAGOS(4))
  READ(14,*) header_units_IAGOS(:)
  WRITE(unitFile_check,*) header_units_IAGOS(:)
  ALLOCATE(tab_var_units_IAGOS(nvar_tot,4))
  DO ivar=1,nvar_tot
     READ(14,*) tab_var_units_IAGOS(ivar,:)
     WRITE(unitFile_check,*) tab_var_units_IAGOS(ivar,:)
  END DO
  CLOSE(14)
  WRITE(unitFile_check,*) '###################################################'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Model variable names.
  var_names_file=TRIM(tmpdir)//'table_variable_names_in_models.txt'
  OPEN(14,FILE=TRIM(var_names_file),FORM='formatted')
  ALLOCATE(header_models(n_models_registered+2))
  ALLOCATE(model_names(n_models_registered))
  READ(14,*) header_models
  WRITE(unitFile_check,*) 'n_models_registered=', n_models_registered
  WRITE(unitFile_check,*) 'ndim_model=',ndim_model
  WRITE(unitFile_check,*) 'nvar_model=',nvar_model
  ALLOCATE(tab_model_names(nvar_model+ndim_model+2,n_models_registered+2))
  WRITE(unitFile_check,*) '############# Models variable names ###############'
  WRITE(unitFile_check,*) header_models(:)
  DO ivarmod=1,(nvar_model+ndim_model+2) ! +2 stands for surface pressure and pressure. Neither accounted for in variables nor dimensions.
     READ(14,*) tab_model_names(ivarmod,:)
     WRITE(*,*) tab_model_names(ivarmod,:)
     WRITE(unitFile_check,*) tab_model_names(ivarmod,:)
  END DO
  model_names=header_models(2:n_models_registered+1)
  WRITE(unitFile_check,*) model_names
  CLOSE(14)
  WRITE(unitFile_check,*) '###################################################'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Model variable units.
  var_units_file=TRIM(tmpdir)//'table_variable_units_in_models.txt'
  OPEN(14,FILE=TRIM(var_units_file),FORM='formatted')
  ALLOCATE(header_units_models(2*n_models_registered+2)) ! 2 columns per model, and 2 additional columns.
  READ(14,*) header_units_models(:)
  ALLOCATE(tab_var_units_models(nvar_model+ndim_model+2,2*n_models_registered+2))
  WRITE(unitFile_check,*) '############# Models variable units ###############'
  WRITE(unitFile_check,*) header_units_models(:)
  DO ivarmod=1,(nvar_model+ndim_model+2) ! +2 stands for surface pressure and pressure. Neither accounted for in variables nor dimensions.
     READ(14,*) tab_var_units_models(ivarmod,:)
     WRITE(*,*) tab_var_units_models(ivarmod,:)
     WRITE(unitFile_check,*) tab_var_units_models(ivarmod,:)
  END DO
  WRITE(unitFile_check,*) '###################################################'
  CLOSE(14)
  ! Defining the metrics characteristics.
  list_stats_obs(1)='Mean'
  list_stats_obs(2)='SD'
  list_stats_obs(3)='Min'
  list_stats_obs(4)='Max'
  list_stats_obs(5)='N_obs'
  list_stats_obs(6)='Day_1'
  list_stats_obs(7)='Day_N'
  list_stats_obs(8)='N_diff_updown'
  list_stats_obs(9)='N_flights'
  list_stats_obs(10)='Pressure'
  list_stats_model(1)='Mean'
  list_stats_model(2)='SD'
  list_stats_model(3)='Min'
  list_stats_model(4)='Max'
  ! Definition for the istat indexes. See function which_one_str defined at the end of the module.
  istat_mean=         which_one_str(c_Mean         , list_stats_obs, nstats_obs)
  istat_SD=           which_one_str(c_SD           , list_stats_obs, nstats_obs)
  istat_min=          which_one_str(c_Min          , list_stats_obs, nstats_obs)
  istat_max=          which_one_str(c_Max          , list_stats_obs, nstats_obs)
  istat_N=            which_one_str(c_N_obs        , list_stats_obs, nstats_obs)
  istat_day_1=        which_one_str(c_Day_1        , list_stats_obs, nstats_obs)
  istat_day_N=        which_one_str(c_Day_N        , list_stats_obs, nstats_obs)
  istat_N_diff_updown=which_one_str(c_N_diff_updown, list_stats_obs, nstats_obs)
  istat_N_flights    =which_one_str(c_N_flights    , list_stats_obs, nstats_obs)
  istat_Pressure     =which_one_str(c_Pressure     , list_stats_obs, nstats_obs)
  ! Same for the model.
  istat_mean_model=which_one_str(c_Mean_model, list_stats_model, nstats_model)
  istat_SD_model=  which_one_str(c_SD_model  , list_stats_model, nstats_model)
  istat_min_model= which_one_str(c_Min_model , list_stats_model, nstats_model)
  istat_max_model= which_one_str(c_Max_model , list_stats_model, nstats_model) 
  WRITE(unitFile_check,*) 'istat_mean=' ,istat_mean                 , list_stats_obs(istat_mean)
  WRITE(unitFile_check,*) 'istat_SD='   ,istat_SD                   , list_stats_obs(istat_SD)
  WRITE(unitFile_check,*) 'istat_min='  ,istat_min                  , list_stats_obs(istat_min)
  WRITE(unitFile_check,*) 'istat_max='  ,istat_max                  , list_stats_obs(istat_max)
  WRITE(unitFile_check,*) 'istat_N='    ,istat_N                    , list_stats_obs(istat_N)
  WRITE(unitFile_check,*) 'istat_day_1=',istat_day_1                , list_stats_obs(istat_day_1)
  WRITE(unitFile_check,*) 'istat_day_N=',istat_day_N                , list_stats_obs(istat_day_N)
  WRITE(unitFile_check,*) 'istat_N_diff_updown=',istat_N_diff_updown, list_stats_obs(istat_N_diff_updown)
  WRITE(unitFile_check,*) 'istat_N_flights=',istat_N_flights        , list_stats_obs(istat_N_flights)
  WRITE(unitFile_check,*) 'istat_Pressure=',istat_Pressure          , list_stats_obs(istat_Pressure)
  WRITE(unitFile_check,*) 'istat_mean_model=' ,istat_mean_model     , list_stats_obs(istat_mean_model)
  WRITE(unitFile_check,*) 'istat_SD_model='   ,istat_SD_model       , list_stats_obs(istat_SD_model)
  WRITE(unitFile_check,*) 'istat_min_model='  ,istat_min_model      , list_stats_obs(istat_min_model)
  WRITE(unitFile_check,*) 'istat_max_model='  ,istat_max_model      , list_stats_obs(istat_max_model)
  ! This array below is the argument used in the tonetcdf subroutine for observations.
  ! It couples the names of the variables and the names of the metrics.
  ALLOCATE(list_variables_with_stats(nvariables_with_stats))
  DO ivar=1,nvariables
     DO istat=1,nstats_obs
        list_variables_with_stats(ivar+(istat-1)*nvariables)= &
             TRIM(list_variables_nc(ivar))//'_'//TRIM(list_stats_obs(istat))
     END DO
  END DO
  DEALLOCATE(list_variables_nc)
  WRITE(unitFile_check,*) 'list_variables_with_stats=', list_variables_with_stats
  ! Finds the column containing the variables full names.
  icol_model_full_name=which_one_str(c_Full_name,header_models,n_models_registered+2)
  ! Finds the right column to adopt the right names for variables and the right conversion factor.
  ! If it is in a scientific program registered in the tables, then we look for the column
  ! of this scientific program rather than the individual model.
  icol_program = which_one_str(sci_program,header_models,n_models_registered+2)
  IF (icol_program .ne. 0) THEN
     icol_model = icol_program
  ELSE
     icol_model=which_one_str(cmodel,header_models,n_models_registered+2)
  END IF
  WRITE(unitFile_check,*) 'icol_model=', icol_model, header_models(icol_model)
  ! The following index corresponds to the list of models present in the table files.
  ! It is easier to use this one to find the right columns in the units table. 
  ! -1 because the first column is for the variable names entered by the user, and is not proper to a model.
  imodel=icol_model-1
  ! Finds the right column for the variables output names.
  icol_model_output=which_one_str(c_Output_name, header_models, n_models_registered+2)
  WRITE(unitFile_check,*) 'icol_model_output=',icol_model_output, header_models(icol_model_output)
  ! Finds the right column for the variables output units.
  icol_unit_model_output=which_one_str(c_Output_unit, header_units_models, 2*n_models_registered+2)
  WRITE(unitFile_check,*) 'icol_unit_model_output=',icol_unit_model_output, header_units_models(icol_unit_model_output)
  ! Same for the observations.
  icol_IAGOS_output     =which_one_str(c_Output_name, header_packages   , npackages+2)
  icol_unit_IAGOS_output=which_one_str(c_Output_unit, header_units_IAGOS, 4)
  icol_unit_IAGOS_input =which_one_str(c_Input_unit , header_units_IAGOS, 4)
  ! The next two arrays are the dimensions/variables names as they
  ! should be found in the model output files.
  ALLOCATE(model_dimensions(ndim_model))
  ALLOCATE(model_variables(nvar_model))
  ! Units for the model variables, without accounting for the stats.
  ALLOCATE(model_units_dim_output(ndim_model))
  ALLOCATE(model_units_var_output(nvar_model))
  ! Units for the model variables, accounting for the stats. To be written in the NetCDF files.
  nvar_model_with_stats = nvar_model * nstats_model
  ALLOCATE(model_units_output(nvar_model_with_stats + ndim_model))
  ALLOCATE(model_conv_factor(nvar_model))
  model_dimensions=tab_model_names(1:ndim_model,icol_model)
  WRITE(unitFile_check,*) 'model_dimensions=', model_dimensions
  model_variables =tab_model_names((ndim_model+1):(ndim_model+nvar_model),icol_model)
  model_var_output=tab_model_names((ndim_model+1):(ndim_model+nvar_model),icol_model_output)
  ! Particular case: pressure. Last row from tab_model_names, but we do not want it as an output from this routine.
  irow_model_P_S=which_one_str(c_Surface_pressure,tab_model_names(:,icol_model_full_name),ndim_model+nvar_model+2)
  irow_model_P  =which_one_str(c_Pressure        ,tab_model_names(:,icol_model_full_name),ndim_model+nvar_model+2)
  WRITE(*,*) 'icol_model=',icol_model
  WRITE(*,*) 'irow_model_P_S=',irow_model_P_S
  WRITE(*,*) 'irow_model_P=',irow_model_P
  model_P_S=tab_model_names(irow_model_P_S,icol_model)
  model_P  =tab_model_names(irow_model_P  ,icol_model)
  WRITE(*,*) 'model_P_S=', model_P_S
  WRITE(*,*) 'model_P='  , model_P
  ALLOCATE(list_var_model_with_stats(nvar_model_with_stats))
  DO istat=1,nstats_model
     DO ivar=1,nvar_model
        list_var_model_with_stats(ivar+(istat-1)*nvar_model)= &
             TRIM(model_var_output(ivar))//'_'//TRIM(list_stats_model(istat))
     END DO
  END DO
  WRITE(unitFile_check,*) 'model_var_output=',model_var_output  
  WRITE(unitFile_check,*) 'list_var_model_with_stats=',list_var_model_with_stats
  coord_names_output=tab_model_names(1:ndim_model,icol_model_output)
  model_units_dim_output=tab_var_units_models(1:ndim_model,icol_unit_model_output)
  model_units_var_output= &
       tab_var_units_models((ndim_model+1) : (ndim_model+nvar_model), icol_unit_model_output)
  WRITE(unitFile_check,*) 'model_units_dim_output========', model_units_dim_output
  WRITE(unitFile_check,*) 'model_units_var_output========', model_units_var_output
  model_units_output(1:ndim_model) = model_units_dim_output
  DO istat=1,nstats_model
     model_units_output( ((istat-1)*nvar_model+1+ndim_model) : &
          (istat*nvar_model+ndim_model)) = model_units_var_output
  END DO
  WRITE(unitFile_check,*) 'model_units_output========', model_units_output
  ! Defines the name for the time record dimension inquiry.
  idimmod_time=which_one_str(c_Time, list_dim_model, ndim_model)
  ctime_model=model_dimensions(idimmod_time)
  WRITE(unitFile_check,*) 'ctime_model=', ctime_model
  ! Also finds the names for longitude and latitude in the model files.
  idimmod_lon=which_one_str(c_Longitude, list_dim_model, ndim_model)
  clon_model=model_dimensions(idimmod_lon)
  WRITE(unitFile_check,*) 'clon_model=', clon_model
  idimmod_lat=which_one_str(c_Latitude, list_dim_model, ndim_model)
  clat_model=model_dimensions(idimmod_lat)
  WRITE(unitFile_check,*) 'clat_model=', clat_model
  idimmod_lev=which_one_str(c_Vertical_grid_level, list_dim_model, ndim_model)
  clev_model=model_dimensions(idimmod_lev)
  WRITE(unitFile_check,*) 'clev_model=', clev_model
  ! Conversion string --> real.
  DO ivarmod=1,nvar_model
     READ(tab_var_units_models(ivarmod+ndim_model,imodel*2+2),*) model_conv_factor(ivarmod)
  END DO
  ALLOCATE(IAGOS_units_output(nvariables))
  ALLOCATE(IAGOS_units_input(nvariables))
  ALLOCATE(IAGOS_conv_factor(nvariables))
  IAGOS_units_output=tab_var_units_IAGOS((ncoord+1):nvar_tot,icol_unit_IAGOS_output)
  IAGOS_units_input =tab_var_units_IAGOS((ncoord+1):nvar_tot,icol_unit_IAGOS_input)
  WRITE(unitFile_check,*) 'IAGOS_units_output=',IAGOS_units_output
  WRITE(unitFile_check,*) 'IAGOS_units_input=' ,IAGOS_units_input
  ! Conversion string --> real.
  DO ivar=1,nvariables
     READ(tab_var_units_IAGOS(ivar+ncoord,4),*) IAGOS_conv_factor(ivar)
  END DO
  ALLOCATE(IAGOS_units_output_with_stats(nvariables_with_stats))
  ! For each stat, we fill all the variable boxes with a units array, depending on
  ! the kind of stat.
  DO istat=1,nstats_obs
     select case (TRIM(list_stats_obs(istat)))
        case ('Mean', 'SD', 'Min', 'Max')
           IAGOS_units_output_with_stats( ((istat-1)*nvariables+1) : (istat*nvariables) ) = &
                IAGOS_units_output
        case ('N_obs','N_diff_updown')
           IAGOS_units_output_with_stats( ((istat-1)*nvariables+1) : (istat*nvariables) ) = &
                'data'
        case ('N_flights')
           IAGOS_units_output_with_stats( ((istat-1)*nvariables+1) : (istat*nvariables) ) = &
                'flights'
        case('Day_1', 'Day_N')
           IAGOS_units_output_with_stats( ((istat-1)*nvariables+1) : (istat*nvariables) ) = &
                'th_day'
        case('Pressure')
           IAGOS_units_output_with_stats( ((istat-1)*nvariables+1) : (istat*nvariables) ) = &
                'hPa'
        end select
     END DO
  DEALLOCATE(IAGOS_units_output)
  WRITE(unitFile_check,*) 'IAGOS_units_output_with_stats=', IAGOS_units_output_with_stats
  ! ############## Generates the horizontal grids ###########
  ! The lon/lat_netcdf vectors are the coordinates as they appear in the output NetCDF file.
  ! They refer to the cells southwest limit if we use the MOCAGE model.
  WRITE(*,*) 'Defining the horizontal grid.'
  rlat(1)=lat1+shift_lat
  rlat_netcdf(1)=lat1
  DO k=2,nlat
     rlat(k)=rlat(k-1)+dlat
     rlat_netcdf(k)=rlat_netcdf(k-1)+dlat
  END DO
  WRITE(*,*) rlat(1), '°N --> ', rlat(nlat),'°N'
  rlon(1)=lon1+shift_lon
  rlon_netcdf(1)=lon1
  DO i=2,nlon
     rlon(i)=rlon(i-1)+dlon
     rlon_netcdf(i)=rlon_netcdf(i-1)+dlon
  END DO
  WRITE(*,*) rlon_netcdf(1), '°E --> ', rlon_netcdf(nlon),'°E'
  WRITE(*,*) nlev, tmpdir
  ! ############## Generates the vertical grid #############
  call INT2CHAR(nlev,tmpdir,cnlev)
  WRITE(*,*) 'Reading the vertical grid coefficients.'
  WRITE(*,*) cnlev,'vertical levels'
  IF (.not. logic_P .and. .not. logic_cst_P) THEN
     unitFile=20
     fileName=TRIM(ref_dir)//'zvamh_'//TRIM(cnlev)//'.txt'
     OPEN(unitFile,FILE=fileName,FORM='formatted',STATUS='old')
     DO l=1,nlev
        READ(unitFile,*) A_inter(l)
     END DO
     CLOSE(unitFile)
     fileName=TRIM(ref_dir)//'zvbmh_'//TRIM(cnlev)//'.txt'
     OPEN(unitFile,FILE=fileName,FORM='formatted',STATUS='old')
     DO l=1,nlev
        READ(unitFile,*) B_inter(l)
     END DO
     CLOSE(unitFile)
  END IF
  ! Just for testing the fields on a mid-latitude gridcell.
  ilon_check=int(nlon/2) ; klat_check=int(nlat*4/5) ; llev_check=INT(nlev*3/5)
  ! ilon_ok, klat_ok, llev_ok: 
  ! just computing the last sampled coordinates, in order to check the arrays on them later.
  ALLOCATE(ilon_ok(nlayers))
  ALLOCATE(klat_ok(nlayers))
  ALLOCATE(llev_ok(nlayers))
  !     #####################################
  !     The loop on months starts here.  ####
  DO imth=1,nmonths                     !####
     iyyyymm=list_months(imth)          !####
     !  #####################################
     WRITE(unitFile_check,*) 'imth, list_months(imth)=', imth, list_months(imth)
     WRITE(*,*) 'imth=',imth
     WRITE(*,*) 'iyyyymm=',iyyyymm
     ALLOCATE(mean_5D(nvariables,nlayers,nlon,nlat,nlev))
     ALLOCATE(  SD_5D(nvariables,nlayers,nlon,nlat,nlev))
     ALLOCATE( min_5D(nvariables,nlayers,nlon,nlat,nlev))
     ALLOCATE( max_5D(nvariables,nlayers,nlon,nlat,nlev))
     ALLOCATE( rcount(nvariables,nlayers,nlon,nlat,nlev))
     ALLOCATE(day1_5D(nvariables,nlayers,nlon,nlat,nlev))
     ALLOCATE(dayN_5D(nvariables,nlayers,nlon,nlat,nlev))
     ALLOCATE(n_diff_updown(nvariables,nlayers,nlon,nlat,nlev))
     ALLOCATE( N_flights_5D(nvariables,nlayers,nlon,nlat,nlev))
     ALLOCATE(   P_5D(nvariables,nlayers,nlon,nlat,nlev))
     
     call INT2CHAR(iyyyymm,tmpdir,cyyyymm) ! Converts the integer date into characters
     iyyyy=nint(REAL(iyyyymm)/100)
     call INT2CHAR(iyyyy,tmpdir,cyyyy) ! Converts the integer date into characters
     ! First, reads the address of the IAGOS files in a dedicated file.
     fileName=TRIM(tmpdir)//'list_iagos_files/'//TRIM(cyyyy)//'/'//TRIM(cyyyymm)
     ALLOCATE(list_data_files(list_nfiles(imth)))
     OPEN(unitFile_list,FILE=fileName,FORM='formatted')
     DO iFile=1,list_nfiles(imth)
        READ(unitFile_list,*) list_data_files(iFile)
     END DO
     CLOSE(unitFile_list)
     WRITE(*,*) '#####################################################'
     WRITE(*,*) 'Reads the list of IAGOS files for the current month: ', &
          TRIM(cyyyymm),' (imth=',imth,')'
     WRITE(*,*) TRIM(fileName)
     WRITE(*,*) 'Number of files for this month: ', list_nfiles(imth)
     ! First reading the model time dimension for the current month.
     ! If the pressure levels are not predefined, then we use the ones from the current model.
     IF (logic_P .EQV. .TRUE. .OR. logic_cst_P .EQV. .TRUE.) THEN
        filename_P=TRIM(init_model_output_dir)//TRIM(cyyyy)//'/P_'//TRIM(cyyyymm)//'.nc'
     ELSE ! if (logic_cst_P .EQV. .TRUE.) THEN ! Then, no pressure file to 
        filename_P=TRIM(init_model_output_dir)//TRIM(cyyyy)//'/P_S_'//TRIM(cyyyymm)//'.nc'
     END IF
     WRITE(*,*) TRIM(filename_P)
     il_err=NF90_OPEN (TRIM(filename_P),NF90_NOWRITE,il_ncfile)
     il_err=NF90_INQ_DIMID (il_ncfile, TRIM(ctime_model), il_ncdim_time)
     WRITE(*,*) TRIM(ctime_model)
     IF (il_err .NE. 0) THEN
        WRITE(*,*) NF90_STRERROR(il_err)
     END IF
     il_err=NF90_INQUIRE_DIMENSION (il_ncfile, il_ncdim_time, LEN = ntime)
     ! ntime=1 if the model outputs are monthly resolved.
     WRITE(unitFile_check,*) 'ntime= ', ntime
     ALLOCATE(P_mid(nlon,nlat,nlev,ntime))
!!!!!! Reading the pressure fields (at the middle of grid cells) in the model file.
     IF (logic_P .EQV. .TRUE. .AND. .NOT. logic_cst_P) THEN
        il_err=NF90_INQ_VARID(il_ncfile, model_P, il_ncvar)
        IF (il_err.eq.0) THEN
           pres_file_exists=.TRUE.
        ELSE
           pres_file_exists=.FALSE.
           WRITE(unitFile_check,*) 'variable P not found. NF90_INQ_VARID returned',il_err
           WRITE(*,*) 'variable P not found. NF90_INQ_VARID returned',il_err
           STOP
        END IF
        il_err=NF90_GET_VAR(il_ncfile,il_ncvar,P_mid)
        WRITE(unitFile_check,*) 'filename_P=', filename_P
        ! Checking the daily pressures.
        DO itime=1,ntime
           WRITE(unitFile_check,*) 'P_mid(nlon/2,4/5*nlat,nlev,itime)=', P_mid(ilon_check,klat_check,nlev,itime)
        END DO
     END IF
     ALLOCATE(P_surf(nlon,nlat,ntime))
     ALLOCATE(P_edge(nlon,nlat,nlev,ntime))
     ALLOCATE(P_full(nlon,nlat,nlev,ntime))
     P_edge(:,:,:,:)=rwrong_val
     P_full(:,:,:,:)=rwrong_val
     IF (.NOT. logic_cst_P .AND. .NOT. logic_P) THEN
        ! If the pressure file was the surface pressures, then we generate here the 3D pressure fields,
        ! both for half and full levels.
        il_err=NF90_INQ_VARID(il_ncfile,model_P_S,il_ncvar)
        IF (il_err.NE.0) THEN
           WRITE(*,*) 'variable P_S not found. NF90_INQ_VARID returned',il_err
           WRITE(unitFile_check,*) 'variable P_S not found. NF90_INQ_VARID returned',il_err
        END IF
        il_err=NF90_GET_VAR(il_ncfile,il_ncvar,P_surf)
        WRITE(unitFile_check,*) 'filename_Psurf=', filename_Psurf
        ! Checking the daily surface pressures.
        DO itime=1,ntime
           WRITE(unitFile_check,*) 'P_surf(nlon/2,4/5*nlat,itime)=', P_surf(ilon_check,klat_check,itime)
        END DO
        IF (il_err.EQ.0) THEN
           pres_file_exists=.TRUE.
        ELSE
           pres_file_exists=.FALSE.
        END IF
        IF (pres_file_exists .AND. .NOT. logic_P) THEN
           CALL COMPUTE_PRESSURES(nlon, nlat, nlev, ntime, P_surf(:,:,:), P_0, &
                A_inter(:), B_inter(:), P_edge(:,:,:,:), P_full(:,:,:,:))
        END IF
        IF (.NOT. logic_P) THEN
           P_mid=P_full
        END IF
        ! The following loop is just in order to check the program behaviour,
        ! and also the intra-monthly pressure variability in case of daily resolution.
        ! It allows us to verify that changes in surface pressure have very
        ! small influence on the pressure at the UTLS levels. In other words,
        ! choosing a daily or monthly resolution in surface pressure
        ! has no significant impact on the cruise data distribution.
        ! For vertical profiles however, it shows that this impact would be strong,
        ! at least in the lower troposphere.
        DO itime=1,ntime
           WRITE(unitFile_check,*) 'P_mid(nlon/2,4/5*nlat,l=24,imth,itime)=', &
                P_mid(ilon_check,klat_check,llev_check,itime),itime
        END DO
     ELSE IF (logic_cst_P) THEN
        P_edge(:,:,:,:)=rwrong_val
        P_full(:,:,:,:)=rwrong_val
        ALLOCATE(cst_P_levels(nlev))
        OPEN(15, FILE=TRIM(filename_cst_P_levels), FORM='formatted')
        READ(15,*) cst_P_levels(:)
        CLOSE(15)
        DO l=1,nlev
           P_mid(:,:,l,:) = cst_P_levels(l)
        END DO
        DEALLOCATE(cst_P_levels)
     END IF ! if (.not. logic_cst_P)
     il_err=NF90_CLOSE(il_ncfile)
!!!!!! Reading (or not) the potential vorticity fields.
     ! We hopp this step if we ask the program to treat the layer 'all_layers' or/and 'UTLS' only,
     ! for the climatology maps, for example. Or for comparing with gridded IASI measurements
     ! between 9 and 12 km.
     ! To be used later, for the separation between the different layers.
     ALLOCATE(PV(nlon,nlat,nlev,ntime))
     IF (nlayers.eq.1) THEN 
        logic_PV = TRIM(list_layers(1)).ne.'UTLS' .and. TRIM(list_layers(1)).ne.'all_layers'
     ELSE IF (nlayers.eq.2) THEN
        logic_PV = (TRIM(list_layers(1)).ne.'UTLS' .or. TRIM(list_layers(2)).ne.'all_layers') .and. &
             (TRIM(list_layers(2)).ne.'UTLS' .or. TRIM(list_layers(1)).ne.'all_layers')
     ELSE 
        logic_PV=.TRUE.
     END IF
     IF (logic_PV .EQV. .TRUE.) THEN
        WRITE(*,*) 'Reading the PV files.'
        filename_PV=TRIM(model_output_dir)//TRIM(cyyyy)//'/PV_'//TRIM(cyyyymm)//'.nc'
        il_err=NF90_OPEN(TRIM(filename_PV),NF90_NOWRITE,il_ncfile)
        il_err=NF90_INQ_VARID(il_ncfile,'PV',il_ncvar)
        IF (il_err.ne.0) THEN
           WRITE(*,*) 'variable PV not found. NF90_INQ_VARID returned',il_err
        END IF
        il_err=NF90_GET_VAR(il_ncfile,il_ncvar,PV)
        il_err=NF90_CLOSE(il_ncfile)
     ELSE
        WRITE(*,*) 'Only one layer, without PV-based distinction.'
        PV(:,:,:,:)=0.
     END IF
     ! The following loop checks the intramonthly PV variations.
     ! The results clearly show the importance of a daily resolution
     ! in the UT/LS separation!
     DO itime=1,ntime
        WRITE(*,*) 'PV(',ilon_check,klat_check,llev_check,itime,')=', &
             PV(ilon_check,klat_check,llev_check,itime), 'itime=', itime
        WRITE(unitFile_check,*) 'PV(ilon_check,klat_check,llev_check,itime)=', &
             PV(ilon_check,klat_check,llev_check,itime), 'itime=', itime
     END DO
     ! rcount_daily: data counting vector in order to filter out non-sampled days 
     ! during the model monthly averaging step.
     ALLOCATE(rcount_daily(nvariables,nlayers,nlon,nlat,nlev,ntime))
     ALLOCATE(var_model_daily(nvar_model,nlon,nlat,nlev,ntime))
     ALLOCATE(mean_5D_model(nvar_model,nlayers,nlon,nlat,nlev))
     ALLOCATE(SD_5D_model(nvar_model,nlayers,nlon,nlat,nlev))
     ALLOCATE(min_5D_model(nvar_model,nlayers,nlon,nlat,nlev))
     ALLOCATE(max_5D_model(nvar_model,nlayers,nlon,nlat,nlev))

     data_dir_date = TRIM(cyyyy)//'/'//TRIM(cyyyymm)//'/'
     data_dir = TRIM(data_dir_cst)//TRIM(data_dir_date)
     filename_model=TRIM(model_output_dir)//TRIM(cyyyy)//'/' &
          //TRIM(experiment)//'_'//TRIM(cyyyymm)//'.nc'
     IF (.NOT. pres_file_exists .AND. .NOT. logic_cst_P) THEN
        WRITE(*,*) filename_Psurf
        WRITE(*,*) '---> Surface pressure field missing: cannot distribute the observations on a grid for this month.'
        WRITE(unitFile_check,*) 'Surface pressure field missing: cannot distribute the observations on a grid for this month.'
        mean_5D(:,:,:,:,:)=rwrong_val
        SD_5D(:,:,:,:,:)  =rwrong_val
        min_5D(:,:,:,:,:) =rwrong_val
        max_5D(:,:,:,:,:) =rwrong_val
        rcount(:,:,:,:,:) =0.
        day1_5D(:,:,:,:,:)=rwrong_val
        dayN_5D(:,:,:,:,:)=rwrong_val
        n_diff_updown(:,:,:,:,:)=rwrong_val
        N_flights_5D(:,:,:,:,:)    =0.
        P_5D(:,:,:,:,:)         =rwrong_val
        var_model_daily(:,:,:,:,:)=rwrong_val
        mean_5D_model  (:,:,:,:,:)=rwrong_val
        SD_5D_model    (:,:,:,:,:)=rwrong_val
        min_5D_model   (:,:,:,:,:)=rwrong_val
        max_5D_model   (:,:,:,:,:)=rwrong_val
     ELSE
        WRITE(*,*) '#######################################'
        WRITE(*,*) 'Calling mean_values_interpol subroutine'
        WRITE(*,*) '#######################################'
        SD_5D      (:,:,:,:,:)=rwrong_val
        SD_5D_model(:,:,:,:,:)=rwrong_val
        CALL MEAN_VALUES_INTERPOL(                                                                   &
!!! 1/ Constant parameters
!!!!! Observations
             filename_prefix_IAGOS, filename_suffix_IAGOS,                                           &
             list_coord(:), list_variables(:), list_der_variables(:), list_der_variables_nc(:),      &
             thres_val(:), IAGOS_conv_factor(:), header_packages(:), tab_var_names(:,:),             &
             list_packages(:),                                                                       &
!!!!! Model
             list_dim_model(:), list_var_model(:), model_variables(:), header_models(:),             &
             model_units_output(:), model_conv_factor(:),                                            &
!!!!! Layer-selection parameters.
             list_layers(:), PV_down_layers(:), PV_up_layers(:), P_down_layers(:), P_up_layers(:),   &
             ozone_max_layers(:), ozone_min_layers(:), logic_PV, PV_tropopause,                      &
!!! 2/ Monthly parameters
             iyyyymm, ntime, list_nfiles(imth), list_data_files(:),                                  &
!!!!! Grid definition
             rlon(:), rlat(:), dlon, dlat, shift_lon, shift_lat,                                     &
!!!!! Input meteorological fields
             P_edge(:,:,:,:), P_mid(:,:,:,:), PV(:,:,:,:),                                           &
!!!!! Daily outputs, temporary arrays useful for SD derivation.
             rcount_daily(:,:,:,:,:,:), var_model_daily(:,:,:,:,:),                                  &
!!!!! 5D-outputs, to be written in the output files.
             mean_5D(:,:,:,:,:), rcount(:,:,:,:,:), n_diff_updown (:,:,:,:,:),                       &
             min_5D(:,:,:,:,:), max_5D(:,:,:,:,:), day1_5D(:,:,:,:,:), dayN_5D(:,:,:,:,:),           &
             N_flights_5D(:,:,:,:,:), P_5D(:,:,:,:,:),                                               &
             mean_5D_model(:,:,:,:,:), min_5D_model(:,:,:,:,:), max_5D_model(:,:,:,:,:),             &
             ilon_ok(:), klat_ok(:), llev_ok(:)                                                      &
             )
        DEALLOCATE(P_mid)
        DEALLOCATE(PV)
        DEALLOCATE(var_model_daily)
        DEALLOCATE(list_data_files)
        DEALLOCATE(rcount_daily)
        DEALLOCATE(P_surf, P_edge, P_full)
        WRITE(*,*) 'Checking all the metrics, in the case of ozone.'
        DO ilay=1,nlayers
           WRITE(*,*) '-------------'
           WRITE(*,*) TRIM(list_layers(ilay))//' : ilon_ok, klat_ok, llev_ok=', &
                ilon_ok(ilay), klat_ok(ilay), llev_ok(ilay)
           WRITE(*,*) '-------------'
           IF (ilon_ok(ilay)*klat_ok(ilay)*llev_ok(ilay) .ne. 0) THEN
              WRITE(*,*) 'N obs='     ,rcount       (ivarO3,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
              WRITE(*,*) 'Mean='      ,mean_5D      (ivarO3,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
              WRITE(*,*) 'SD='        ,SD_5D        (ivarO3,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
              WRITE(*,*) 'Min='       ,min_5D       (ivarO3,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
              WRITE(*,*) 'Max='       ,max_5D       (ivarO3,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
              WRITE(*,*) 'Day 1='     ,day1_5D      (ivarO3,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
              WRITE(*,*) 'Day N='     ,dayN_5D      (ivarO3,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
              WRITE(*,*) 'N_diff_updown=',N_diff_updown(ivarO3,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
              WRITE(*,*) 'N_flights_5D=',N_flights_5D(ivarO3,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
              WRITE(*,*) 'P_5D='      ,P_5D         (ivarO3,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
              WRITE(*,*) 'Mean model=',Mean_5D_model(ivarO3,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
              WRITE(*,*) 'Min model=' ,min_5D_model (ivarO3,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
              WRITE(*,*) 'Max model=' ,max_5D_model (ivarO3,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
           END IF
        END DO
!!!!!!!! Now gathering all the final data (coupling variables and metrics) into one single matrix.
        ALLOCATE(mean_5D_with_stats(nvariables_with_stats,nlayers,nlon,nlat,nlev))
        DO ivar=1,nvariables
           mean_5D_with_stats(ivar+(istat_mean -1)*nvariables,:,:,:,:)=mean_5D(ivar,:,:,:,:)
           mean_5D_with_stats(ivar+(istat_SD   -1)*nvariables,:,:,:,:)=SD_5D  (ivar,:,:,:,:)
           mean_5D_with_stats(ivar+(istat_min  -1)*nvariables,:,:,:,:)=min_5D (ivar,:,:,:,:)
           mean_5D_with_stats(ivar+(istat_max  -1)*nvariables,:,:,:,:)=max_5D (ivar,:,:,:,:)
           mean_5D_with_stats(ivar+(istat_N    -1)*nvariables,:,:,:,:)=rcount (ivar,:,:,:,:)
           mean_5D_with_stats(ivar+(istat_day_1-1)*nvariables,:,:,:,:)=day1_5D(ivar,:,:,:,:)
           mean_5D_with_stats(ivar+(istat_day_N-1)*nvariables,:,:,:,:)=dayN_5D(ivar,:,:,:,:)
           mean_5D_with_stats(ivar+(istat_N_diff_updown-1)*nvariables,:,:,:,:)= &
                n_diff_updown(ivar,:,:,:,:)
           mean_5D_with_stats(ivar    +(istat_N_flights-1)*nvariables,:,:,:,:)= &
                N_flights_5D(ivar,:,:,:,:)
           mean_5D_with_stats(ivar+(istat_Pressure-1)*nvariables,:,:,:,:)= &
                P_5D(ivar,:,:,:,:)
        END DO
        DEALLOCATE(mean_5D)
        DEALLOCATE(SD_5D)
        DEALLOCATE(min_5D)
        DEALLOCATE(max_5D)
        DEALLOCATE(rcount)
        DEALLOCATE(day1_5D, dayN_5D)
        DEALLOCATE(n_diff_updown, N_flights_5D)
        DEALLOCATE(P_5D)
        ALLOCATE(mean_5D_model_with_stats(nvar_model_with_stats,nlayers,nlon,nlat,nlev))
        DO ivar=1,nvar_model
           mean_5D_model_with_stats(ivar+(istat_mean_model-1)*nvar_model,:,:,:,:)= &
                mean_5D_model(ivar,:,:,:,:)
           mean_5D_model_with_stats(ivar+(istat_SD_model  -1)*nvar_model,:,:,:,:)= &
                SD_5D_model  (ivar,:,:,:,:)
           mean_5D_model_with_stats(ivar+(istat_min_model -1)*nvar_model,:,:,:,:)= &
                min_5D_model (ivar,:,:,:,:)
           mean_5D_model_with_stats(ivar+(istat_max_model -1)*nvar_model,:,:,:,:)= &
                max_5D_model (ivar,:,:,:,:)
        END DO
        DEALLOCATE(mean_5D_model)
        DEALLOCATE(SD_5D_model)
        DEALLOCATE(min_5D_model)
        DEALLOCATE(max_5D_model)
        IF (ilon_ok(ilay)*klat_ok(ilay)*llev_ok(ilay) .ne. 0) THEN
           WRITE(*,*) 'mean_5D_model'
           IF (SUM(mean_5D_with_stats(ivarO3+nvariables*(istat_N-1),:,:,:,:)).gt.1) THEN
              DO ivar=1,nvariables
                 DO ilay=1,nlayers
                    WRITE(*,*)'ilon_ok,klat_ok,llev_ok=',ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay)
                    WRITE(*,*) TRIM(list_variables(ivar)),' ',TRIM(list_layers(ilay)),'(ilon_ok,klat_ok,llev_ok)=', &
                         mean_5D_with_stats(ivar,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
                    WRITE(unitFile_check,*)'ilon_ok,klat_ok,llev_ok=',ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay)
                    WRITE(unitFile_check,*) TRIM(list_variables(ivar)),' ', &
                         TRIM(list_layers(ilay)),'(ilon_ok,klat_ok,llev_ok)=', &
                         mean_5D_with_stats(ivar,ilay,ilon_ok(ilay),klat_ok(ilay),llev_ok(ilay))
                 END DO
              END DO
           ELSE
              WRITE(*,*)'ilon_ok,k_ok,l_ok do not exist because no data this month'
              WRITE(unitFile_check,*)'ilon_ok,klat_ok,llev_ok do not exist because no data this month'
           END IF
        END IF
        DO ilay=1,nlayers
           ! Writes gridded IAGOS data in separate layers.
           outfile_NC=TRIM(outdir_NC)//TRIM(list_layers(ilay))//'/'// &
                TRIM(cyyyy)//'/IAGOS_'//TRIM(cyyyymm)//'.nc'
           CALL TONETCDF(mean_5D_with_stats(:,ilay,:,:,:), nlon, nlat, nlev, 1, coord_names_output(:), &
                nvariables_with_stats, list_variables_with_stats(:), IAGOS_units_output_with_stats(:), &
                rlon_netcdf(:), rlat_netcdf(:), cyyyymm, rwrong_val, outfile_NC)
           ! Writes model monthly fields in separate layers.
           outfile_NC=TRIM(outdir_NC)//TRIM(list_layers(ilay))//'/'// &
                TRIM(cyyyy)//'/'//TRIM(experiment)//'_'//TRIM(cyyyymm)//'.nc'
           CALL TONETCDF(mean_5D_model_with_stats(:,ilay,:,:,:), nlon, nlat, nlev, 1, coord_names_output(:),&
                nvar_model_with_stats, list_var_model_with_stats(:),                                            &
                model_units_output((ndim_model+1):(ndim_model+nvar_model_with_stats)),                          &
                rlon_netcdf(:), rlat_netcdf(:), cyyyymm, rwrong_val, outfile_NC)
           WRITE(*,*) list_layers(ilay), 'done'
        END DO
     END IF ! END IF .not. pres_file_exists
     DEALLOCATE(mean_5D_with_stats)
     DEALLOCATE(mean_5D_model_with_stats)
  END DO ! end imth=1,nmonths
  DEALLOCATE(ilon_ok, klat_ok, llev_ok)
  DEALLOCATE(A_inter, B_inter)
  DEALLOCATE(list_der_variables)
  DEALLOCATE(list_der_variables_nc)
  DEALLOCATE(list_layers)
  DEALLOCATE(PV_down_layers)
  DEALLOCATE(PV_up_layers)
  DEALLOCATE(P_down_layers)
  DEALLOCATE(P_up_layers)
  WRITE(*,*) '### Interpol_IAGOS: done. ###'
  WRITE(unitFile_check,*) '### Interpol_IAGOS: done. ###'
  CLOSE(unitFile_check)
END PROGRAM INTERPOL_IAGOS
