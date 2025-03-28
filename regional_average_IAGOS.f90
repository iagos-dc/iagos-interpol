! This routine is made to derive monthly time series from gridded IAGOS and model data,
! for each region required by the user, which must have pretreated files with only
! their local gridcells only. The time series are written in two dimensions: months from
! yyyymm1 until yyyymm2, and the "five" seasons (the fifth is the whole year,
! i.e. without any seasonal separation).
! Concerning the observations, the averaging process sums up the amounts of measurements
! made in each gridcell. It counts the duration (in days) separating the first and the last
! measurements at the regional scale in order to have a look on the temporal representativity
! of each monthly mean on its corresponding month. It also counts the amount of involved
! gridcells in each average, this time for spatial representativity.
! Concerning the model, the averages are calculated following a IAGOS mask, defined for each variable.
! The latter mask is based on the corresponding observed variable if it exists and if it has been
! required by the user.
! Else, observed ozone is the default mask variable.

PROGRAM REGIONAL_AVERAGE_IAGOS
  USE NETCDF
  USE GLOBAL_VAR_INTERPOL_MOD
  USE FUNCTIONS_MOD
  IMPLICIT NONE
!!!!!! Parameters !!!!!!
  ! Stats for counts.
  INTEGER,PARAMETER :: istat_count_Ncells=1
  ! INTEGER,PARAMETER :: istat_count_Nflights=2
  INTEGER,PARAMETER :: istat_count_Dt=2
  INTEGER,PARAMETER :: nstats_count=2
  ! Season for each month of the year.
  CHARACTER*3 seasons_monthly(n12months)
  ! Current season.
  CHARACTER*3 this_season
  ! List of NEW counting statistics (thus excluding the amount of measurements).
  CHARACTER*20 list_stats_count(nstats_count)
  ! Names for some variables. The 30-character length is important to match
  ! with the which_one_str function input argument declaration.
  CHARACTER*30 c_O3
!!!!!! Standard variables. !!!!!!
  ! Length of the character string read in a NetCDF file.
  INTEGER char_len
  ! Various indexes. ivar for the variable;
  ! ivar_mask for the observed variable to be used as a mask on the current model variable;
  ! istat for the metric; istat_count for a new counting metric;
  ! ivar_stat for the couple variable x metric;
  INTEGER ivar_mask, istat, istat_count, ivar_stat
  ! imth for the current month; iyyyymm, iyyyy and imm for the date; isea for the season;
  ! i, k, l for the longitude, latitude and vertical coordinate respectively.
  INTEGER imth, iyyyymm, iyyyy, imm, isea, i, k, l
  ! Amount of observations for the current variable.
  INTEGER ivar_stat_N
  ! Amount of flights contributing to the mean value, and of validated gridcells.
  INTEGER nflights, ncells_ok
  ! Variables and their units.
  CHARACTER*30 cvar, cunit
!!!!!!
  ! Variables extracted from the (regional) namelist file
  ! and their derived products.
  ! Period limits (included).
  INTEGER iyyyymm1, iyyyymm2
  ! Period duration, in months.
  INTEGER nmonths
  ! List of months, format yyyymm, as integers and as characters.
  INTEGER,ALLOCATABLE :: list_months(:)
  CHARACTER*20,ALLOCATABLE :: clist_months(:)
  ! Top and bottom vertical grid levels to be involved in the averaging process.
  INTEGER l_top, l_down
  ! Same, in case we include some layers from the troposphere.
  INTEGER l_top_UTLS, l_down_UTLS, l_top_FT, l_down_FT
  ! Dimension names in input files.
  CHARACTER*30 lon_dim_name, lat_dim_name, lev_dim_name
  ! Diverse lists of variables.
  ! listvar_IAGOS and listvar_model are the variables names, both for input and output files.
  ! Not that in both kinds of files, a whole "variable" name is made of one variable
  ! and one metric. These whole names are stored in the listvar_..._with_stats 
  ! or with_all_stats arrays. The _with_stats suffix means it gathers all the metrics
  ! present in the input files. The _with_all_stats suffix adds the new metrics
  ! calculated in this routine in the _with_stats arrays. They are the names of this routine outputs.
  ! The counting metrics are only meaningful for observations, and do not have to be
  ! copied into the model files. Post-treatment routines know it has to be read
  ! from IAGOS files only.
  ! listvar_count is the coupling between variable and new counting metrics.
  CHARACTER*30,ALLOCATABLE :: listvar_IAGOS(:), listvar_model(:), listvar_count(:,:)
  ! Same organisation for the units, since they depend on the variable and the metric.
  ! For example, "ppb" is for ozone if the metric is "Mean", "SD", "Min", etc. but
  ! "data" if the metric is "N_obs".
  CHARACTER*30,ALLOCATABLE :: listvar_units_model(:), listvar_units_count(:,:)
  CHARACTER*30,ALLOCATABLE :: listvar_units_IAGOS(:)
  CHARACTER*30,ALLOCATABLE :: listvar_IAGOS_with_stats(:), listvar_units_IAGOS_with_stats(:)
  CHARACTER*30,ALLOCATABLE :: listvar_IAGOS_with_all_stats(:), listvar_units_IAGOS_with_all_stats(:)
  CHARACTER*30,ALLOCATABLE :: listvar_model_with_stats(:), listvar_units_model_with_stats(:)
  ! Min equivalent amount of data for a given gridcell, on a given month.
  ! Just to ensure a minimum representativeness in each gridcell monthly mean.
  REAL rcount_min
  CHARACTER*30 region
  CHARACTER*30 layer
  ! Directories. dir_regions for the input files, outdir_time_series for the output files.
  CHARACTER*300 dir_regions, outdir_time_series
  ! Regional IAGOS and model files.
  CHARACTER*300 input_file_IAGOS, input_file_model
  CHARACTER*300 output_file_IAGOS, output_file_model
  ! Character strings for integers.
  CHARACTER*30 cyyyy
  CHARACTER*30 cyyyymm
  ! First and last days in the month with a measurement.
  INTEGER day_1, day_N
  ! Data.
  ! var_4D: the array containing monthly statistics for each variable, and for each gridcell
  ! in the current region.
  REAL, ALLOCATABLE :: var_4D(:,:,:,:), var_4D_model(:,:,:,:)
  ! real, ALLOCATABLE :: var_3D(:,:,:),rcount_3D(:,:,:)
  ! output_var_3D: the array containing monthly statistics for each variable, 
  ! for each month and for each season.
  REAL, ALLOCATABLE :: output_var_3D(:,:,:), output_var_3D_model(:,:,:)
  ! Same as the previous ones, including the new counted metrics for each variable,
  ! i.e. the following rsupp_stats_4D array.
  REAL, ALLOCATABLE :: output_var_3D_with_all_stats(:,:,:)
  REAL, ALLOCATABLE :: rsupp_stats_4D(:,:,:,:)
  ! istat_count=1 for the equivalent amount of data
  ! istat_count=2 for the amount of sampled gridcells

  tmpdir=TRIM(pdir)//'tmp/'
  CALL FILL_STR_CHAR(tmpdir,6,LEN(tmpdir),tmpdir)

  c_O3='O3'
  ! Seasons of the year and, by extent, the whole year.
  seasons(1)='DJF'
  seasons(2)='MAM'
  seasons(3)='JJA'
  seasons(4)='SON'
  seasons(5)='ANN'
  ! Corresponding season for each month of the year.
  seasons_monthly(12)='DJF'
  seasons_monthly(1:2)='DJF'
  seasons_monthly(3:5)='MAM'
  seasons_monthly(6:8)='JJA'
  seasons_monthly(9:11)='SON'
  ! Statistics
  list_stats_obs(1) ='Mean'
  list_stats_obs(2) ='SD'
  list_stats_obs(3) ='Min'
  list_stats_obs(4) ='Max'
  list_stats_obs(5) ='N_obs'
  list_stats_obs(6) ='N_flights'
  list_stats_obs(7) ='Day_1'
  list_stats_obs(8) ='Day_N'
  list_stats_obs(9) ='N_diff_updown'
  list_stats_obs(10)='Pressure'
  list_stats_obs(11)='Delta_P_wrt_TP'
  list_stats_obs(12)='Delta_z_wrt_TP'
  list_stats_model(1)='Mean'
  list_stats_model(2)='SD'
  list_stats_model(3)='Min'
  list_stats_model(4)='Max'
  ! Definition for the istat indexes. See function which_one_str defined in functions_mod.f90.
  istat_mean=         which_one_str(c_Mean         , list_stats_obs, nstats_obs)
  istat_SD=           which_one_str(c_SD           , list_stats_obs, nstats_obs)
  istat_min=          which_one_str(c_Min          , list_stats_obs, nstats_obs)
  istat_max=          which_one_str(c_Max          , list_stats_obs, nstats_obs)
  istat_N=            which_one_str(c_N_obs        , list_stats_obs, nstats_obs)
  istat_N_flights=    which_one_str(c_N_flights    , list_stats_obs, nstats_obs)
  istat_day_1=        which_one_str(c_Day_1        , list_stats_obs, nstats_obs)
  istat_day_N=        which_one_str(c_Day_N        , list_stats_obs, nstats_obs)
  istat_N_diff_updown=which_one_str(c_N_diff_updown, list_stats_obs, nstats_obs)
  istat_Pressure=     which_one_str(c_Pressure     , list_stats_obs, nstats_obs)
  istat_Delta_P=      which_one_str(c_Delta_P_wrt_TP,list_stats_obs, nstats_obs)
  istat_Delta_z=      which_one_str(c_Delta_z_wrt_TP,list_stats_obs, nstats_obs)
  ! Same for the model.
  istat_mean_model=which_one_str(c_Mean_model, list_stats_model, nstats_model)
  istat_SD_model=  which_one_str(c_SD_model  , list_stats_model, nstats_model)
  istat_min_model= which_one_str(c_Min_model , list_stats_model, nstats_model)
  istat_max_model= which_one_str(c_Max_model , list_stats_model, nstats_model) 
  ! The new counting metrics.
  list_stats_count(istat_count_Ncells)='N_gridcells'
  list_stats_count(istat_count_Dt)='Dt_days'
  ! list_stats_count(istat_count_Nflights)='N_flights'
  
  OPEN(14,FILE=TRIM(pdir)//'namelist_regional',FORM='formatted')
  WRITE(unitFile_check,*) 'Reading namelist_regional file.'
  READ(14,*) experiment
  READ(14,*) iyyyymm1, iyyyymm2
  READ(14,*) nmonths
  ALLOCATE(list_months(nmonths))
  READ(14,*) list_months(1:nmonths)
  READ(14,*) l_top_UTLS, l_down_UTLS
  READ(14,*) l_top_FT, l_down_FT
  READ(14,*) lon_dim_name, lat_dim_name, lev_dim_name ! The latter: name of the vertical dimension.
  READ(14,*) nvar_IAGOS
  ALLOCATE(listvar_IAGOS(nvar_IAGOS))
  READ(14,*) listvar_IAGOS(1:nvar_IAGOS)
  READ(14,*) nvar_model
  ALLOCATE(listvar_model(nvar_model))
  READ(14,*) listvar_model(1:nvar_model)
  READ(14,*) rcount_min
  READ(14,*) tmpdir_there
  READ(14,*) region
  READ(14,*) layer
  READ(14,*) dir_regions
  READ(14,*) outdir_time_series
  CLOSE(14)
  OPEN(unitFile_check, FILE=TRIM(tmpdir_there)//'check_average_iagos_'// &
       TRIM(layer)//'_'//TRIM(region), FORM='formatted')
  WRITE(unitFile_check,*) 'experiment=', TRIM(experiment)
  WRITE(unitFile_check,*) 'iyyyymm1, iyyyymm2=', iyyyymm1, iyyyymm2
  WRITE(unitFile_check,*) 'nmonths=', nmonths
  WRITE(unitFile_check,*) 'list_months(1:nmonths)=', list_months(1:nmonths)
  WRITE(unitFile_check,*) 'l_top_UTLS, l_down_UTLS=', l_top_UTLS, l_down_UTLS
  WRITE(unitFile_check,*) 'l_top_FT, l_down_FT=', l_top_FT, l_down_FT
  WRITE(unitFile_check,*) 'lon_dim_name, lat_dim_name, lev_dim_name=', lon_dim_name, lat_dim_name, lev_dim_name
  WRITE(unitFile_check,*) 'nvar_IAGOS=', nvar_IAGOS
  WRITE(unitFile_check,*) 'listvar_IAGOS(1:nvar_IAGOS)= ', listvar_IAGOS(1:nvar_IAGOS)
  WRITE(unitFile_check,*) 'nvar_model=', nvar_model
  WRITE(unitFile_check,*) 'listvar_model(1:nvar_model)=', listvar_model(1:nvar_model)
  WRITE(unitFile_check,*) 'rcount_min=', rcount_min
  WRITE(unitFile_check,*) 'tmpdir=', tmpdir
  WRITE(unitFile_check,*) 'tmpdir_there=', tmpdir_there
  WRITE(unitFile_check,*) 'region=', region
  WRITE(unitFile_check,*) 'layer=', layer
  WRITE(unitFile_check,*) 'dir_regions=', dir_regions
  WRITE(unitFile_check,*) 'outdir_time_series=', outdir_time_series
  WRITE(*,*) 'region=', region
  WRITE(*,*) 'layer=', layer
!!! Defining the layer boundaries.
  IF (TRIM(layer) .EQ.'UT' .OR. &
       TRIM(layer).EQ.'TPL'.OR. &
       TRIM(layer).EQ.'LS' .OR. &
       TRIM(layer).EQ.'UTLS') THEN
     l_top = l_top_UTLS
     l_down = l_down_UTLS
  ELSE IF (TRIM(layer).EQ.'LFT'        .OR. &
       TRIM(layer)    .EQ.'MFT'        .OR. &
       TRIM(layer)    .EQ.'UBL'        .OR. &
       TRIM(layer)    .EQ.'all_layers') THEN
     l_top = l_top_FT
     l_down = l_down_FT
  ELSE
     WRITE(*,*) 'Unrecognized layer. Stopping the execution.'
     STOP
  END IF
!!! 1/ Defining the variables and units arrays, with regard to the lists of metrics.
  ALLOCATE(listvar_units_IAGOS(nvar_IAGOS))
  ALLOCATE(listvar_units_count(nvar_IAGOS, nstats_count))
  ALLOCATE(listvar_units_model(nvar_model))
!!! Defining the count variables.
  ALLOCATE(listvar_count(nvar_IAGOS,nstats_count))
  DO ivar=1,nvar_IAGOS
     DO istat_count=1,nstats_count
        listvar_count(ivar, istat_count)=TRIM(listvar_IAGOS(ivar))//'_'// &
             TRIM(list_stats_count(istat_count))
     END DO
  END DO
  WRITE(unitFile_check,*) 'listvar_count=', listvar_count
  ! Defining the whole sets of couples variables x stats.
  nvar_IAGOS_with_stats=nvar_IAGOS*nstats_obs
  nvar_model_with_stats=nvar_model*nstats_model
  ! WRITE(*,*) 'nvar_IAGOS_with_stats=', nvar_IAGOS_with_stats, nvar_IAGOS, nstats_obs
  ALLOCATE(listvar_IAGOS_with_stats(nvar_IAGOS_with_stats))
  ALLOCATE(listvar_model_with_stats(nvar_model_with_stats))
  ALLOCATE(listvar_units_IAGOS_with_stats(nvar_IAGOS_with_stats))
  ALLOCATE(listvar_units_model_with_stats(nvar_model_with_stats))

  DO ivar=1,nvar_IAGOS
     DO istat=1,nstats_obs
        ivar_stat=ivar+(istat-1)*nvar_IAGOS
        listvar_IAGOS_with_stats(ivar_stat)=TRIM(listvar_IAGOS(ivar))//'_'//TRIM(list_stats_obs(istat))
     END DO
  END DO
  WRITE(unitFile_check,*) 'listvar_IAGOS_with_stats=', listvar_IAGOS_with_stats
  DO ivar=1,nvar_model
     DO istat=1,nstats_model
        ivar_stat=ivar+(istat-1)*nvar_model
        listvar_model_with_stats(ivar_stat)=TRIM(listvar_model(ivar))//'_'//TRIM(list_stats_model(istat))
     END DO
  END DO
  WRITE(unitFile_check,*) 'listvar_model_with_stats=', listvar_model_with_stats
!!! 2/ Filling the data arrays.
  ALLOCATE(output_var_3D(nvar_IAGOS_with_stats,n5seasons,nmonths))
  ALLOCATE(rsupp_stats_4D(nvar_IAGOS,nstats_count,n5seasons,nmonths))
  ALLOCATE(output_var_3D_model(nvar_model_with_stats,n5seasons,nmonths))
  ! Initialization.
  output_var_3D(:,:,:)      =rwrong_val
  output_var_3D_model(:,:,:)=rwrong_val
  rsupp_stats_4D(:,:,:,:)   =0
  DO imth=1,nmonths
!!! 2.1/ Reads the regional files, for each data set.
     iyyyymm=list_months(imth)
     imm=MODULO(iyyyymm,100)
     ! What is the season please?
     this_season=seasons_monthly(imm)
     isea=1
     DO WHILE (seasons(isea).NE.this_season)
        isea=isea+1
     END DO
     CALL INT2CHAR(iyyyymm,tmpdir,cyyyymm) ! Converts the integer date into characters
     iyyyy=NINT(REAL(iyyyymm)/100)
     CALL INT2CHAR(iyyyy,tmpdir,cyyyy) ! Converts the integer year into characters
     input_file_IAGOS=TRIM(dir_regions)//TRIM(region)//'/'//TRIM(layer)//'/'//TRIM(cyyyy)// &
          '/IAGOS_'//TRIM(cyyyymm)//'.nc'
     input_file_model=TRIM(dir_regions)//TRIM(region)//'/'//TRIM(layer)//'/'//TRIM(cyyyy)// &
          '/'//TRIM(experiment)//'_'//TRIM(cyyyymm)//'.nc'
!!! 2.1.1/ We start by reading the IAGOS data, especially the amounts of data. 
!!! It will be used as a reference to pick up the gridcells we account for in the average computation.
     il_err=NF90_OPEN(TRIM(input_file_IAGOS),NF90_NOWRITE,il_ncfile)
     ! Longitude
     il_err=NF90_INQ_DIMID(il_ncfile,TRIM(lon_dim_name),il_ncdim)
     il_err=NF90_INQUIRE_DIMENSION (il_ncfile, il_ncdim, LEN = nlon)
     ! Latitude
     il_err=NF90_INQ_DIMID(il_ncfile,TRIM(lat_dim_name),il_ncdim)
     il_err=NF90_INQUIRE_DIMENSION (il_ncfile, il_ncdim, LEN = nlat)
     ! Vertical grid level
     il_err=NF90_INQ_DIMID(il_ncfile,TRIM(lev_dim_name),il_ncdim)
     il_err=NF90_INQUIRE_DIMENSION (il_ncfile, il_ncdim, LEN = nlev)
     ALLOCATE(var_4D(nvar_IAGOS_with_stats,nlon,nlat,nlev))
     DO ivar=1,nvar_IAGOS
        ! 1st: we need to know how many gridcells are sampled enough.
        ivar_stat_N=ivar+(istat_N-1)*nvar_IAGOS ! Index for the amount of data in the coupled list.
        il_err=NF90_INQ_VARID(il_ncfile,TRIM(listvar_IAGOS_with_stats(ivar_stat_N)),il_ncvar)
        il_err=NF90_GET_VAR(il_ncfile,il_ncvar,var_4D(ivar_stat_N,:,:,:))
        ncells_ok=COUNT(var_4D(ivar_stat_N,:,:,l_top:l_down) .GE. rcount_min)
        ! 2nd: we calculate the time series.
        DO istat=1,nstats_obs
           ivar_stat=ivar+(istat-1)*nvar_IAGOS
           il_err=NF90_INQ_VARID(il_ncfile,TRIM(listvar_IAGOS_with_stats(ivar_stat)),il_ncvar)
           il_err=NF90_GET_VAR(il_ncfile,il_ncvar,var_4D(ivar_stat,:,:,:))
           ! Copying the units attribute.
           ! il_err = NF_INQ_ATTLEN (il_ncfile, il_ncvar, 'units', char_len)
           ! il_err = NF_GET_ATT_TEXT (il_ncfile, il_ncvar, 'units', cunit(1:char_len))
           il_err = NF90_GET_ATT_STRING (il_ncfile, il_ncvar, 'units', cunit)
           ! listvar_units_IAGOS(ivar,1:char_len)=cunit(1:char_len)
           listvar_units_IAGOS_with_stats(ivar_stat)=cunit(1:char_len)
           ! ... and also for the amounts of data.
           listvar_units_count(ivar,istat_count_Ncells)='avail_gridcells'
           listvar_units_count(ivar,istat_count_Dt)='days'
           ! listvar_units_count(ivar,istat_count_Nflights)='flights'
!!! 2.1.2/ Computing averages and other metrics.
           IF (ncells_ok.GE.1) THEN
              ! Several options for the last operation, depending on the current stat.
              SELECT CASE (TRIM(list_stats_obs(istat)))
              CASE ('Mean', 'SD', 'Min', 'Max', 'Pressure')
                 output_var_3D(ivar_stat,isea,imth)=SUM(var_4D(ivar_stat,:,:,l_top:l_down), &
                      var_4D(ivar_stat_N,:,:,l_top:l_down) .GE. rcount_min)/ncells_ok
                 ! Careful: works if rcount_min > 1 and not 0, because of SD which
                 ! denominator was n-1 instead of n.
              CASE ('N_obs')
                 output_var_3D(ivar_stat,isea,imth)=SUM(var_4D(ivar_stat,:,:,l_top:l_down), &
                      var_4D(ivar_stat_N,:,:,l_top:l_down) .GE. rcount_min)
              CASE ('N_flights') ! This simplified definition will underestimate N_flights.
                 ! Still better than the sum, which could overestimate by a factor N_cells.
                 output_var_3D(ivar_stat,isea,imth)=MAXVAL(var_4D(ivar_stat,:,:,l_top:l_down), &
                      var_4D(ivar_stat_N,:,:,l_top:l_down) .GE. rcount_min)
              CASE ('Day_1')
                 output_var_3D(ivar_stat,isea,imth)=MINVAL(var_4D(ivar_stat,:,:,l_top:l_down), &
                      var_4D(ivar_stat_N,:,:,l_top:l_down) .GE. rcount_min)
                 day_1=MINVAL(var_4D(ivar_stat,:,:,l_top:l_down), &
                      var_4D(ivar_stat_N,:,:,l_top:l_down) .GE. rcount_min)
              CASE ('Day_N')
                 output_var_3D(ivar_stat,isea,imth)=MAXVAL(var_4D(ivar_stat,:,:,l_top:l_down), &
                      var_4D(ivar_stat_N,:,:,l_top:l_down) .GE. rcount_min)
                 day_N=MAXVAL(var_4D(ivar_stat,:,:,l_top:l_down), &
                      var_4D(ivar_stat_N,:,:,l_top:l_down) .GE. rcount_min)
                 rsupp_stats_4D(ivar,istat_count_Ncells,isea,imth)=ncells_ok
                 rsupp_stats_4D(ivar,istat_count_Dt,isea,imth)= day_N-day_1
              END SELECT
           ELSE
              output_var_3D(ivar_stat,isea,imth)=rwrong_val
              rsupp_stats_4D(ivar,:,isea,imth)  =0
           END IF
           output_var_3D(ivar_stat,n5seasons,imth)=output_var_3D(ivar_stat,isea,imth)
           rsupp_stats_4D(ivar,:,n5seasons,imth)=rsupp_stats_4D(ivar,:,isea,imth)
        END DO
     END DO
     il_err=NF90_CLOSE(il_ncfile)
!!! 2.1.3/ Now we read the model data. We just care that we apply the IAGOS mask on the model cells.
     il_err=NF90_OPEN(TRIM(input_file_model),NF90_NOWRITE,il_ncfile)
     ! No need to read the dimension sizes. It has to be the same as in the IAGOS file.
     ALLOCATE(var_4D_model(nvar_model_with_stats,nlon,nlat,nlev))
     DO ivar=1,nvar_model
        ivar_mask=1
!!! Finds the IAGOS masking variable. If the current one is not amongst the IAGOS ones,
!!! or does not have a relevant observed peer to be filtered with, then we use ozone.
        IF (ANY(listvar_IAGOS(:) .EQ. listvar_model(ivar))) THEN
           ivar_mask=which_one_str(listvar_model(ivar), listvar_IAGOS, nvar_IAGOS)
        ELSE
           ivar_mask=which_one_str(c_O3, listvar_IAGOS, nvar_IAGOS)
        END IF
        ivar_stat_N=ivar_mask+(istat_N-1)*nvar_IAGOS ! Index for the amount of data in the (obs) coupled list.
        DO istat=1,nstats_model
           ivar_stat=ivar+(istat-1)*nvar_model
           il_err=NF90_INQ_VARID(il_ncfile,TRIM(listvar_model(ivar))//'_'//TRIM(list_stats_model(istat)),il_ncvar)
           il_err=NF90_GET_VAR(il_ncfile,il_ncvar,var_4D_model(ivar_stat,:,:,:))
           ! Copying the units attribute.
           ! il_err = NF_INQ_ATTLEN (il_ncfile, il_ncvar, 'units', char_len)
           ! il_err = NF_GET_ATT_TEXT (il_ncfile, il_ncvar, 'units', cunit(1:char_len))

           il_err = NF90_GET_ATT_STRING (il_ncfile, il_ncvar, 'units', cunit)
           
           listvar_units_model_with_stats(ivar_stat)=cunit(1:char_len)
!!! 2.2/ Computing averages and other metrics.
           ncells_ok=COUNT(var_4D(ivar_stat_N,:,:,l_top:l_down) .GE. rcount_min .AND. &
                var_4D_model(ivar_stat,:,:,l_top:l_down) .NE. rwrong_val)
           IF (ncells_ok .GE. 1) THEN
              output_var_3D_model(ivar_stat,isea,imth)=SUM(var_4D_model(ivar_stat,:,:,l_top:l_down), &
                   var_4D(ivar_stat_N,:,:,l_top:l_down) .GE. rcount_min .AND. &
                   var_4D_model(ivar_stat,:,:,l_top:l_down) .NE. rwrong_val)/ncells_ok
              ! WRITE(*,*) 'Counting anomalous missval: '
              ! WRITE(*,*) count(var_4D_model(ivar_stat,:,:,l_top:l_down) .lt. 0 .and. &
              !      var_4D(ivar_mask+(istat_SD-1)*nvar_IAGOS,:,:,l_top:l_down) .gt. 0)
              IF (output_var_3D_model(ivar_stat,isea,imth).LT.0 .AND. &
                   TRIM(listvar_model(ivar)).NE.'Net_O3') THEN ! It is ok to find negative values for Net_O3. 
                 WRITE(*,*) 'Negative model average: '//TRIM(listvar_model(ivar))
                 WRITE(*,*) layer, output_var_3D_model(ivar_stat,isea,imth)
                 WRITE(*,*) MINVAL(var_4D_model(ivar_stat,:,:,l_top:l_down), &
                      var_4D(ivar_stat_N,:,:,l_top:l_down) .GE. rcount_min)
                 WRITE(*,*) MINVAL(var_4D(ivar_mask+(istat_SD-1)*nvar_IAGOS,:,:,l_top:l_down), &
                      var_4D(ivar_stat_N,:,:,l_top:l_down) .GE. rcount_min)
                 ! WRITE(*,*) sum(var_4D_model(ivar_stat,:,:,l_top:l_down), &
                 !   var_4D(ivar_stat_N,:,:,l_top:l_down) .ge. rcount_min)/
              END IF
              ! Converting the zero-sums into missing values.
              IF (output_var_3D_model(ivar_stat,isea,imth).EQ.0) THEN
                 output_var_3D_model(ivar_stat,isea,imth)=rwrong_val
              END IF
           ELSE
              output_var_3D_model(ivar_stat,isea,imth)=rwrong_val
           END IF
           output_var_3D_model(ivar_stat,n5seasons,imth)=output_var_3D_model(ivar_stat,isea,imth)
        END DO ! istat
     END DO  ! ivar
     il_err=NF90_CLOSE(il_ncfile)
     DEALLOCATE(var_4D, var_4D_model)
  END DO ! imth=1,nmonths
  ! 3/ Writing the time series in the output file.
  output_file_IAGOS=TRIM(outdir_time_series)//'IAGOS_'//TRIM(region)//'.nc'
  nvar_IAGOS_with_supp_stats = nvar_IAGOS*nstats_count
  nvar_IAGOS_with_all_stats = nvar_IAGOS_with_stats + nvar_IAGOS_with_supp_stats
  ALLOCATE(listvar_IAGOS_with_all_stats(nvar_IAGOS_with_all_stats))
  ALLOCATE(listvar_units_IAGOS_with_all_stats(nvar_IAGOS_with_all_stats))
  ALLOCATE(output_var_3D_with_all_stats(nvar_IAGOS_with_all_stats,n5seasons,nmonths))
  ! Gathering all the calculated data into one single array.
  ! Careful while incrementing the ivar_stat index: 
  ! the variable loop has to be included into the stats loop.
  ivar_stat=0
  DO istat=1,nstats_obs
     DO ivar=1,nvar_IAGOS
        ivar_stat=ivar_stat+1
        listvar_IAGOS_with_all_stats(ivar_stat) = listvar_IAGOS_with_stats(ivar_stat)
        listvar_units_IAGOS_with_all_stats(ivar_stat) = listvar_units_IAGOS_with_stats(ivar_stat)
        output_var_3D_with_all_stats(ivar_stat,:,:) = output_var_3D(ivar_stat,:,:)
     END DO
  END DO
  DO istat_count=1,nstats_count
     DO ivar=1,nvar_IAGOS
        ivar_stat=ivar_stat+1
        ! Variable names.
        listvar_IAGOS_with_all_stats(ivar_stat)=listvar_count(ivar,istat_count)
        ! Variable units.
        listvar_units_IAGOS_with_all_stats(ivar_stat)=listvar_units_count(ivar,istat_count)
        ! Fields.
        output_var_3D_with_all_stats(ivar_stat,:,:)=rsupp_stats_4D(ivar,istat_count,:,:)
     END DO
  END DO
  DEALLOCATE(output_var_3D,rsupp_stats_4D)
  ! WRITE(*,*) '-----------------------------'
  ! WRITE(*,*) 'list_months =', list_months
  ! WRITE(*,*) 'output_var_3D_with_all_stats(1,iseaDJF,:)=', output_var_3D_with_all_stats(1,iseaDJF,:)
  ! WRITE(*,*) listvar_units_IAGOS_with_all_stats(:)
  ! WRITE(*,*) '-----------------------------'
  ! Last: calling the routines that write the arrays into NetCDF files.
  CALL TONETCDF_TIME_SERIES(output_var_3D_with_all_stats(:,:,:), n5seasons, nmonths,&
       nvar_IAGOS_with_all_stats, listvar_IAGOS_with_all_stats(:),                  &
       listvar_units_IAGOS_with_all_stats(:),                                       &
       list_months(:), rwrong_val, tmpdir, output_file_IAGOS)
  output_file_model=TRIM(outdir_time_series)//TRIM(experiment)//'_'//TRIM(region)//'.nc'
  CALL TONETCDF_TIME_SERIES(output_var_3D_model(:,:,:), n5seasons, nmonths,         &
       nvar_model_with_stats, listvar_model_with_stats(:),                          &
       listvar_units_model_with_stats(:),                                           &
       list_months(:), rwrong_val, tmpdir, output_file_model)
  DEALLOCATE(output_var_3D_with_all_stats)
  DEALLOCATE(output_var_3D_model)
END PROGRAM REGIONAL_AVERAGE_IAGOS
