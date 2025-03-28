MODULE GLOBAL_VAR_INTERPOL_MOD
  IMPLICIT NONE
!!!!!! Parameters !!!!!!
  ! Current directory.
  CHARACTER*2,PARAMETER :: pdir='./'
  ! Tmp directory.
  ! Constants.
  REAL, PARAMETER :: Pi = ACOS(-1.)! 3.14159265
  
  REAL, PARAMETER :: M_O3 = 48. ! Ozone molar mass, in g.mol-1
  REAL, PARAMETER :: M_air = 29. ! Ambient air molar mass, in g.mol-1
  REAL, PARAMETER :: n_seconds_a_year = 31557600. ! 365.25 days into seconds
  REAL, PARAMETER :: N_avo = 6.022E23 ! Avogadro number, in mol-1
  ! character*100,parameter :: pdir='/home/scratch01/ycohen/Pack_code/routines/grid_IAGOS_model/'
  ! Integers referring to different files.
  INTEGER,PARAMETER :: unitFile_list =11 ! list of variables
  INTEGER,PARAMETER :: unitFile_check=12 ! debug file
  ! For the IAGOS data extrapolation. If the order of the extrapolation is n, then nweight=2*n+1.
  ! Actually, it is likely that it only works with nweight=3.
  INTEGER,PARAMETER :: nweight=3
  ! Constant integers.
  INTEGER,PARAMETER :: n_seconds_a_day=86400
  INTEGER,PARAMETER :: n4seasons=4 ! The four seasons.
  INTEGER,PARAMETER :: n5seasons=5 ! The four seasons with the whole year.
  INTEGER,PARAMETER :: iseaDJF=1
  INTEGER,PARAMETER :: iseaMAM=2
  INTEGER,PARAMETER :: iseaJJA=3
  INTEGER,PARAMETER :: iseaSON=4
  INTEGER,PARAMETER :: n12months=12 ! The 12 months in the year.
  INTEGER,PARAMETER :: nstats_obs=10
  INTEGER,PARAMETER :: nstats_obs_column=8 ! The statistics for observations in the column_average routine.
  INTEGER,PARAMETER :: nstats_model=4
  
  ! The following temp_conversion is a factor that depends on the initial measured temperature units.
  ! Now, it equals 0 because the measured temperature is already in Kelvin degrees. It was not the case
  ! in the previous IAGOS files version.
  REAL,PARAMETER :: temp_conversion=0
  ! Default missing value.
  REAL,PARAMETER :: rwrong_val=-999.
  ! real,parameter :: N_min_daily=0.00033
  ! Creating 30-char strings for consistency with the functions that takes them in arguments.
  ! Seasons.
  CHARACTER(LEN=30),PARAMETER :: c_DJF='DJF',                        &
       c_MAM='MAM', c_JJA='JJA', c_SON='SON', c_ANN='ANN'
  ! Variables (with PV and its unit PVU).
  CHARACTER(LEN=30),PARAMETER :: c_Time_UTC='Time_UTC',              &
       c_Latitude='Latitude'          , c_Longitude='Longitude',     &
       c_Baro_altitude='Baro_altitude', c_Pressure='Pressure', c_Pa='Pa'
  CHARACTER(LEN=30),PARAMETER :: c_Time='Time',                      &
       c_Surface_pressure='Surface_pressure',                        &
       c_Vertical_grid_level='Vertical_grid_level',                  &
       c_Potential_vorticity='Potential_vorticity',                  &
       c_PV='PV', c_PVU='PVU'
  
  CHARACTER(LEN=30),PARAMETER :: c_Ozone='Ozone',                    &
       c_Carbon_monoxide='Carbon_monoxide',                          &
       c_Water_vapour='Water_vapour',                                &
       c_Relative_humidity_liquid='Relative_humidity_liquid',        &
       c_Relative_humidity_ice='Relative_humidity_ice',              &
       c_Cloud_liquid_water_particles='Cloud_liquid_water_particles',&
       c_Nitrogen_oxides='Nitrogen_oxides',                          &
       c_Nitrogen_monoxide='Nitrogen_monoxide',                      &
       c_Nitrogen_dioxide='Nitrogen_dioxide',                        &
       c_Nitrogen_reactive_species='Nitrogen_reactive_species',      &
       c_NOy_ACACIA='NOy_ACACIA',                                    &
       c_NOy='NOy',                                                  &
       c_NO_mask_NOy='NO_mask_NOy',                                  &
       c_NO2_mask_NOy='NO2_mask_NOy',                                &
       c_NOx_mask_NOy='NOx_mask_NOy',                                &
       c_HNO3_mask_NOy='HNO3_mask_NOy',                              &
       c_PAN_mask_NOy='PAN_mask_NOy',                                &
       c_O3_mask_NOy='O3_mask_NOy',                                  &
       c_Temperature='Temperature',                                  &
       c_Potential_temperature='Potential_temperature',              &
       c_Ozone_to_CO_ratio='Ozone_to_CO_ratio',                      &
       c_Zonal_wind_speed='Zonal_wind_speed',                        &
       c_Meridional_wind_speed='Meridional_wind_speed'
       
  ! Statistics.
  CHARACTER(LEN=30),PARAMETER :: c_Mean='Mean',               &
       c_SD='SD'      , c_Min='Min'    , c_Max='Max',         &
       c_N_obs='N_obs', c_Day_1='Day_1', c_Day_N='Day_N',     &
       c_N_diff_updown='N_diff_updown', c_N_flights='N_flights'
  CHARACTER(LEN=30),PARAMETER :: c_Mean_model='Mean',         &
       c_SD_model='SD', c_Min_model='Min', c_Max_model='Max'
  CHARACTER(LEN=30),PARAMETER :: c_Output_name='Output_name', &
       c_Output_unit='Output_unit', c_Input_unit='Input_unit'
  ! Column names in the dimensions/variables tables.
  CHARACTER(LEN=30),PARAMETER :: c_Full_name='Full_name'
  
  ! List of seasons and, by extent, the whole year.
  CHARACTER(LEN=3) seasons(n5seasons)
!!!!!!!!!!!!!!!!!!!!!!!!
  ! Temporary directories.
  CHARACTER(LEN=200) tmpdir_there
  ! Paths used to read the IAGOS files.
  CHARACTER(LEN=200) data_dir, data_dir_cst, data_dir_date, column_names_dir_cst
  CHARACTER(LEN=200) column_names_dir, column_names_file
  CHARACTER(LEN=200) data_file
  CHARACTER(LEN=200) tmpdir
  ! Location of the tabular referring all the possible names for each variable.
  ! Indeed, they depend on the data file, more exactly on the IAGOS package
  ! that generated this data file.
  CHARACTER(LEN=200) var_names_file, var_units_file
  ! IAGOS files format. Can be ASCII or NetCDF.
  CHARACTER(LEN=6) iagos_format
  ! Time resolution for model outputs. Can be monthly or daily.
  CHARACTER(LEN=7) time_resol
  CHARACTER(LEN=20) cwrong_val ! Default missing value.
  ! Paths used to read the required information in the model files:
  ! surface pressure Psurf and Ertel's potential vorticity PV.
  ! They are located in model_output_dir.
  ! ref_dir contains several variables used as reference in the whole code,
  ! e.g. the vertical grid coefficients A_i and B_i.
  CHARACTER(LEN=200) filename_P, filename_Psurf, filename_PV, filename_model
  CHARACTER(LEN=300) filename_cst_P_levels ! This one is for the predefined vertical grid, if used.
  CHARACTER(LEN=200) ref_dir, model_output_dir, init_model_output_dir, cmodel, experiment, sci_program
  ! Paths for the output files.
  CHARACTER(LEN=300) outfile_NC, outdir_NC
  ! Some other file.
  CHARACTER(LEN=300) fileName
  ! Name for the IAGOS-package that generated the current data file.
  ! packName is the global attribute, but it is sometimes erroneous and leads to an "oblivion" of data.
  ! packName_correct is retrieved from the name of the ozone variable. More robust than packName.
  CHARACTER(LEN=30) packName, packName_correct
  ! ctime_model is to be the name of the time record dimension in the model NetCDF files.
  CHARACTER(LEN=30) ctime_model, clon_model, clat_model, clev_model
  ! Names for the coordinates, to be written in the output files.
  CHARACTER(LEN=30) coord_names_output(4)
  ! Statistics
  CHARACTER(LEN=30) list_stats_obs(nstats_obs), list_stats_model(nstats_model)

  INTEGER quicktest
  INTEGER unitFile ! Every other file.
  ! Logicals used to truncate the files which contains data from both current and next month.
  logical bin_day_passed, bin_allocate_P
  ! Observations
  ! Amount of packages processed in the IAGOS database.
  INTEGER npackages
  ! A whole set of variable classifications.
  INTEGER ncoord, nvariables, n_obs_variables, n_der_variables ! Thus there will be nvar_iagos IAGOS variables read in the data, and ncoord dimensions in the output NCDF file.
  INTEGER n_counted_variables, nvariables_with_stats, nvar_iagos, nvar_tot
  INTEGER nvar_IAGOS_with_stats ! For the average.f90 subroutine. Size of the input files, concerning the "variables" dimension.
  INTEGER nvar_IAGOS_with_supp_stats, nvar_IAGOS_with_all_stats ! For the same routine. Size of the output files, also concerning the var dimension.
  ! Model
  INTEGER ndim_model, nvar_model, n_models_registered, nvar_model_with_stats
  ! ... and more generally
  INTEGER ndim, nvar, idim, ivar
  ! Grid definition
  ! nlev, nlat, nlon is the spatial size of the domain, in matter of gridcells.
  INTEGER nlon, nlat, nlev
  ! Indexes for checking the diverse fields.
  ! Below in the code, you will find them defined at nlon/2,
  ! nlat*4/5 and 24 (for a 39-level or 47-level vertical grid) respectively.
  INTEGER ilon_check, klat_check, llev_check
  character*30 clon_check, clat_check, clev_check ! To be displayed in a character string.
  INTEGER nlayers
  ! Those indexes will be used to deal with every dimension or variable individually.
  ! They have to be defined once only. They depend on the list of variables required by the user.
  ! This one is not temporary and thus does not change between two IAGOS files.
  INTEGER idimUTC, idimlat, idimlon, idimalt, idimP ! dimensions
  INTEGER ivarO3, ivarCO, ivarP, ivarT, ivarNOx, ivarNO, ivarNO2 ! observed variables
  INTEGER ivarCloud,ivarH2O,ivarRHI,ivarRHL, ivarNOy, ivarU, ivarV ! observed variables again
  INTEGER ivarTheta, ivarO3toCO ! ancillary variables derived from previous ones
  INTEGER list_ivar_NOx(3) ! The vector gathering ivarNOx, ivarNO, ivarNO2.
  INTEGER ivarmod_H2O ! The water vapour model variable.
  ! Another index referring to the surface pressure amongst the list of model variables.
  INTEGER irow_model_P_S, irow_model_P
  INTEGER irow_model_PV
  ! Same for the different metrics.
  INTEGER istat_mean, istat_SD, istat_min, istat_max, istat_N, istat_N_flights
  INTEGER istat_day_1, istat_day_N, istat_N_diff_updown, istat_Pressure
  INTEGER istat_mean_model, istat_SD_model, istat_min_model, istat_max_model
  ! Indexes linked with the Netcdf files processing.
  INTEGER il_err, il_ncfile, il_ncdim_time, il_ncdim_lon, il_ncdim_lat, il_ncdim_lev
  INTEGER il_ncvar, il_ncdim
  ! Character string displaying the error messages from the NF90 instructions.
  CHARACTER*300 char_err_status

END MODULE GLOBAL_VAR_INTERPOL_MOD
