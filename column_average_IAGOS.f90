! This routine is made to...

PROGRAM COLUMN_AVERAGE_IAGOS
  USE GLOBAL_VAR_INTERPOL_MOD
  USE FUNCTIONS_MOD
  USE NETCDF
  IMPLICIT NONE
  
!!!!!! Parameters !!!!!!
  ! Stats for counts.
  ! istat_count_Ncells is the amount of grid cells accounted for in a given climatological mean.
  INTEGER,PARAMETER :: istat_count_Ncells=1
  INTEGER,PARAMETER :: nstats_count=1
  ! List of NEW counting statistics (thus excluding the amount of measurements, notably).
  CHARACTER*30 list_stats_count(nstats_count)
  ! Amount of PREVIOUS counting statistics.
  INTEGER,PARAMETER :: nstats_N_sth=2
  ! Season for each month of the year.
  CHARACTER*3 seasons_monthly(n12months)
  ! Current season.
  CHARACTER*3 this_season
  ! Names for some variables. The 30-character length is important to match
  ! with the which_one_str function input argument declaration.
  CHARACTER*30 c_O3
!!!!!! Standard variables. !!!!!!
  ! Length of the character string read in a NetCDF file.
  INTEGER char_len
  ! Various indexes. ivar for the variable;
  ! ivar_mask for the observed variable to be used as a mask on the current model variable;
  ! istat for the metric;
  ! ivar_stat for the couple variable x metric;
  INTEGER ivar_mask, istat, ivar_stat, ivar_mean, ivar_mean_model, ilon, ilat, ilev
  ! imth for the current month; iyyyymm, iyyyy and imm for the date; isea for the season;
  ! i, k, l for the longitude, latitude and vertical coordinate respectively.
  INTEGER imth, iyyyymm, iyyyy, imm, isea, i, k, l
  ! Amount of observations for the current variable.
  INTEGER ivar_stat_N, ivar_stat_N_flights
  ! Amount of validated gridcells.
  INTEGER ncells_ok
  ! The amount of vertical grid levels to take into account.
  INTEGER nlev_less
  ! The amount of all output statistics (* obs variables), 
  ! combining the input and the ones added from this routine.
  INTEGER nvar_with_all_stats
  ! Variables and their units.
  CHARACTER*30 cvar, cunit
  ! output_dim_names(:) is the names of 4 dimensions, required by the tonetcdf.f90 writing routine.
  ! Here, the third and fourth dimensions (vertical axis and time, resp.) have one element only.
  CHARACTER*30 output_dim_names(2)
!!!!!!
  ! Variables extracted from the (regional) namelist file
  ! and their derived products.
  ! Period limits (included).
  INTEGER iyyyymm1, iyyyymm2
  ! Converted into characters.
  CHARACTER*6 cyyyymm1, cyyyymm2
  ! Top and bottom vertical grid levels to be involved in the averaging process.
  INTEGER l_top, l_down
  ! Same, in case we include some layers from the troposphere.
  INTEGER l_top_UTLS, l_down_UTLS, l_top_FT, l_down_FT
  ! Dimension names in input files.
  CHARACTER*30 lon_dim_name, lat_dim_name, lev_dim_name
  INTEGER nchar, ichar
  CHARACTER*30 tmp_char
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
  CHARACTER*30,ALLOCATABLE :: listvar_IAGOS(:), listvar_model(:)
  ! Same organisation for the units, since they depend on the variable and the metric.
  ! For example, "ppb" is for ozone if the metric is "Mean", "SD", "Min", etc. but
  ! "data" if the metric is "N_obs" and "flights" if "N_flights".
  CHARACTER*30,ALLOCATABLE :: listvar_units_model(:)
  CHARACTER*30,ALLOCATABLE :: listvar_units_IAGOS(:)
  CHARACTER*30,ALLOCATABLE :: listvar_IAGOS_with_stats(:), listvar_units_IAGOS_with_stats(:)
  CHARACTER*30,ALLOCATABLE :: listvar_IAGOS_with_all_stats(:), listvar_units_IAGOS_with_all_stats(:)
  ! nmonths is the duration (in months) of the period defined by the user.
  ! listvar_nmonths is the duration (in months) of the measurement period for each species.
  ! It allows to adjust the filter based on the amount of data.
  INTEGER                     nmonths
  INTEGER     ,ALLOCATABLE :: listvar_nmonths(:)
  CHARACTER*30,ALLOCATABLE :: listvar_model_with_stats(:), listvar_units_model_with_stats(:)
  ! Min equivalent amount of data for a given gridcell, on a given month.
  ! Just to ensure a minimum representativeness in each gridcell monthly mean.
  REAL rcount_min_ref
  ! time_ratio is the ratio between the whole time period length and the one for a given species.
  REAL time_ratio
  ! rcount_min_var is the rcount_min_ref reference threshold multiplied by a factor
  ! decreasing with latitude and with the measurement time period.
  REAL, ALLOCATABLE :: rcount_min_var(:,:)
  CHARACTER*30 season
  CHARACTER*30 layer
  ! Directories. dir_level_by_level for the input files, outdir_time_period for the output files.
  CHARACTER*300 dir_level_by_level, outdir_time_period
  ! Global IAGOS and model files.
  CHARACTER*300 input_file_IAGOS, input_file_model
  CHARACTER*300 output_file_IAGOS, output_file_model
  ! Character strings for integers.
  CHARACTER*4 cyyyy
  CHARACTER*6 cyyyymm
  ! Latitudes and subsequent latitude ratios, that account for the varying size of the grid cells
  ! with respect to latitude.
  REAL, ALLOCATABLE :: lon(:), lat(:), lat_ratio(:)
  ! Data.
  ! var_4D: the array containing monthly statistics for each variable, and for each gridcell
  ! in the current region.
  REAL, ALLOCATABLE :: var_3D_all_lev(:,:,:), var_3D_model_all_lev(:,:,:)
  REAL, ALLOCATABLE :: var_4D(:,:,:,:), var_4D_model(:,:,:,:)
  ! Output arrays.
  REAL, ALLOCATABLE :: outvar_3D(:,:,:), outvar_3D_model(:,:,:), outvar_3D_all_stats(:,:,:)
  INTEGER, ALLOCATABLE :: mask_3D(:,:,:)
  ! REAL, ALLOCATABLE :: var_3D(:,:,:),rcount_3D(:,:,:)
  ! output_var_3D: the array containing monthly statistics for each variable, 
  ! for each month and for each season.
  REAL, ALLOCATABLE :: output_var_3D(:,:,:), output_var_3D_model(:,:,:)
  ! Amount of gridcells involved in each column average.
  INTEGER, ALLOCATABLE :: count_N_gridcells(:,:,:)
  ! Same as the previous ones, including the new counted metrics for each variable,
  ! i.e. the following rsupp_stats_4D array.
  REAL, ALLOCATABLE :: output_var_3D_with_all_stats(:,:,:)
  
  tmpdir=TRIM(pdir)//'tmp/'
  CALL FILL_STR_CHAR(tmpdir,6,LEN(tmpdir),tmpdir)
  
  c_O3='O3' ! 30-char long string for O3. To be changed with the output names in the IAGOS-DM files.
  ! Statistics
  list_stats_obs(1)='Mean'
  list_stats_obs(2)='SD'
  list_stats_obs(3)='Min'
  list_stats_obs(4)='Max'
  list_stats_obs(5)='N_diff_updown'
  list_stats_obs(6)='Pressure'
  list_stats_obs(7)='N_flights'
  list_stats_obs(8)='N_obs'
  list_stats_model(1)='Mean'
  list_stats_model(2)='SD'
  list_stats_model(3)='Min'
  list_stats_model(4)='Max'
  ! Statistics added from this routine.
  list_stats_count(1)='N_gridcells'
  ! Definition for the istat indexes. See function which_one_str defined in functions_mod.f90.
  istat_mean  =which_one_str(c_Mean , list_stats_obs, nstats_obs_column)
  istat_SD    =which_one_str(c_SD   , list_stats_obs, nstats_obs_column)
  istat_min   =which_one_str(c_Min  , list_stats_obs, nstats_obs_column)
  istat_max   =which_one_str(c_Max  , list_stats_obs, nstats_obs_column)
  istat_N_diff_updown=which_one_str(c_N_diff_updown, list_stats_obs, nstats_obs_column)
  istat_Pressure=which_one_str(c_Pressure, list_stats_obs, nstats_obs_column)
  istat_N     =which_one_str(c_N_obs, list_stats_obs, nstats_obs_column)
  istat_N_flights=which_one_str(c_N_flights, list_stats_obs, nstats_obs_column)
  IF (istat_N .NE. nstats_obs_column) THEN
     WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     WRITE(*,*) '!!!!!! Careful: istat_N is not the last stat index,          !!!!!!'
     WRITE(*,*) '!!!!!! whereas a part of the program assumes istat_N=nstats. !!!!!!'
     WRITE(*,*) '!!!!!! A shift in the arrays might be expected.              !!!!!!'
     WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  END IF
  ! Same for the model.
  istat_mean_model=which_one_str(c_Mean_model, list_stats_model, nstats_model)
  istat_SD_model=  which_one_str(c_SD_model  , list_stats_model, nstats_model)
  istat_min_model= which_one_str(c_Min_model , list_stats_model, nstats_model)
  istat_max_model= which_one_str(c_Max_model , list_stats_model, nstats_model) 
  OPEN(14,FILE=TRIM(pdir)//'namelist_columns',FORM='formatted')
  WRITE(unitFile_check,*) 'Reading namelist_columns file.'
  WRITE(*,*) 'Reading namelist_columns file.'
  READ(14,*) experiment
  READ(14,*) cyyyymm1, cyyyymm2
  READ(14,*) l_top_UTLS, l_down_UTLS
  READ(14,*) l_top_FT, l_down_FT
  READ(14,*) lon_dim_name, lat_dim_name, lev_dim_name ! The latter: name of the vertical dimension.
  READ(14,*) nlat
  READ(14,*) nvar_IAGOS
  ALLOCATE(listvar_IAGOS(nvar_IAGOS))
  READ(14,*) listvar_IAGOS(1:nvar_IAGOS)
  ALLOCATE(listvar_nmonths(nvar_IAGOS))
  READ(14,*) listvar_nmonths(1:nvar_IAGOS)
  READ(14,*) nmonths
  READ(14,*) nvar_model
  ALLOCATE(listvar_model(nvar_model))
  READ(14,*) listvar_model(1:nvar_model)
  READ(14,*) rcount_min_ref ! Reference threshold for ozone.
  READ(14,*) tmpdir_there
  READ(14,*) season
  READ(14,*) layer
  READ(14,*) dir_level_by_level
  READ(14,*) outdir_time_period
  READ(14,*) input_file_IAGOS
  READ(14,*) input_file_model
  CLOSE(14)
  READ(cyyyymm1,'(i6)') iyyyymm1
  READ(cyyyymm2,'(i6)') iyyyymm2
  OPEN(unitFile_check, FILE=TRIM(tmpdir_there)//'check_average_iagos_'// &
       TRIM(layer)//'_'//TRIM(season), FORM='formatted')
  WRITE(*,*) 'experiment=', experiment
  WRITE(*,*) 'iyyyymm1, iyyyymm2=', iyyyymm1, iyyyymm2
  WRITE(*,*) 'l_top_UTLS, l_down_UTLS=', l_top_UTLS, l_down_UTLS
  WRITE(*,*) 'l_top_FT, l_down_FT=', l_top_FT, l_down_FT
  WRITE(*,*) 'lon_dim_name, lat_dim_name, lev_dim_name=', lon_dim_name, lat_dim_name, lev_dim_name
  WRITE(*,*) 'nlat=',nlat
  WRITE(*,*) 'nvar_IAGOS=', nvar_IAGOS
  WRITE(*,*) 'listvar_IAGOS(1:nvar_IAGOS)= ', listvar_IAGOS(1:nvar_IAGOS)
  WRITE(*,*) 'nvar_model=', nvar_model
  WRITE(*,*) 'listvar_model(1:nvar_model)=', listvar_model(1:nvar_model)
  WRITE(*,*) 'rcount_min_ref=', rcount_min_ref
  WRITE(*,*) 'tmpdir=', tmpdir
  WRITE(*,*) 'tmpdir_there=', tmpdir_there
  WRITE(*,*) 'season=', season
  WRITE(*,*) 'layer=', layer
  WRITE(*,*) 'dir_level_by_level=', dir_level_by_level
  WRITE(*,*) 'outdir_time_period=', outdir_time_period
  WRITE(*,*) 'input_file_IAGOS=', input_file_IAGOS
  WRITE(*,*) 'input_file_model=', input_file_model
!!! Defining the layer boundaries.
  IF (TRIM(layer) .EQ.'UT' .OR. &
       TRIM(layer).EQ.'TPL'.OR. &
       TRIM(layer).EQ.'LS' .OR. &
       TRIM(layer).EQ.'UTLS') THEN
     l_top = l_top_UTLS
     l_down = l_down_UTLS
  ELSE IF (TRIM(layer).EQ.'UBL'        .OR. &
       TRIM(layer)    .EQ.'LFT'        .OR. &
       TRIM(layer)    .EQ.'MFT'        .OR. &
       TRIM(layer)    .EQ.'all_layers') THEN
     l_top = l_top_FT
     l_down = l_down_FT
  ELSE
     WRITE(*,*) TRIM(layer),': Unrecognized layer. Stopping the execution.'
     STOP
  END IF
  nlev_less=l_down - l_top + 1
!!! 1/ Defining the variables and units arrays, with regard to the lists of metrics.
  ALLOCATE(listvar_units_IAGOS(nvar_IAGOS))
  ALLOCATE(listvar_units_model(nvar_model))
  ! Defining the whole sets of couples variable x stat.
  nvar_IAGOS_with_stats=nvar_IAGOS*nstats_obs_column
  nvar_model_with_stats=nvar_model*nstats_model
  nvar_with_all_stats=(nstats_obs_column + nstats_count)*nvar_IAGOS
  ALLOCATE(listvar_IAGOS_with_stats(nvar_with_all_stats))
  ALLOCATE(listvar_model_with_stats(nvar_model_with_stats))
  ALLOCATE(listvar_units_IAGOS_with_stats(nvar_with_all_stats))
  ALLOCATE(listvar_units_model_with_stats(nvar_model_with_stats))
  DO ivar=1,nvar_IAGOS
     DO istat=1,nstats_obs_column
        ivar_stat=ivar+(istat-1)*nvar_IAGOS
        listvar_IAGOS_with_stats(ivar_stat)=TRIM(listvar_IAGOS(ivar))//'_'//TRIM(list_stats_obs(istat))
     END DO
     DO istat=1,nstats_count
        ivar_stat=ivar+(nstats_obs_column + istat - 1)*nvar_IAGOS
        listvar_IAGOS_with_stats(ivar_stat)=TRIM(listvar_IAGOS(ivar))//'_'//TRIM(list_stats_count(istat))
     END DO
  END DO
  WRITE(*,*) 'listvar_IAGOS_with_stats=', listvar_IAGOS_with_stats
  DO ivar=1,nvar_model
     DO istat=1,nstats_model
        ivar_stat=ivar+(istat-1)*nvar_model
        listvar_model_with_stats(ivar_stat)=TRIM(listvar_model(ivar))//'_'//TRIM(list_stats_model(istat))
     END DO
  END DO
  WRITE(*,*) 'listvar_model_with_stats=', listvar_model_with_stats
!!! 2/ Filling the data arrays.
  ALLOCATE(output_var_3D(nvar_IAGOS_with_stats,nlon,nlat))
  ALLOCATE(output_var_3D_model(nvar_model_with_stats,nlon,nlat))
  ! Initialization.
  output_var_3D(:,:,:)      =rwrong_val
  output_var_3D_model(:,:,:)=rwrong_val
!!! 2.1/ Reads the regional files, for each data set.
!!! 2.1.1/ We start by reading the IAGOS data, especially the amounts of data. 
!!! It will be used as a reference to pick up the gridcells we account for in the average computation.
  il_err = NF90_OPEN (TRIM(input_file_IAGOS), NF90_NOWRITE, il_ncfile)
  il_err = NF90_INQ_DIMID (il_ncfile, TRIM(lon_dim_name), il_ncdim) ! longitude
  il_err = NF90_INQUIRE_DIMENSION (il_ncfile, il_ncdim, LEN = nlon)
  ALLOCATE(lon(nlon))
  il_err = NF90_GET_VAR (il_ncfile, il_ncdim, lon(:))
  il_err = NF90_INQ_DIMID (il_ncfile, TRIM(lat_dim_name), il_ncdim) ! latitude
  il_err = NF90_INQUIRE_DIMENSION (il_ncfile, il_ncdim, LEN = nlat)
  ALLOCATE(lat(nlat))
  il_err = NF90_GET_VAR (il_ncfile, il_ncdim, lat(:))
  ALLOCATE(lat_ratio(nlat))
  lat_ratio(:) = 0.
  CALL LATITUDE_SURFACE_RATIO(nlat, lat(:), lat_ratio(:))
  IF (TRIM(season) .EQ. 'ANN') THEN
     rcount_min_ref = rcount_min_ref*4
  ELSE
     ! In case there are not always 3 months per season.
     ! Just be careful and call a season by the months initials.
     rcount_min_ref = rcount_min_ref*REAL(LEN_TRIM(season))/3
  END IF
  ! Ozone (we expect the ozone column to be contained in every file)
  if (ANY(listvar_IAGOS=='O3')) THEN
     ivarO3= which_one_str(c_O3, listvar_IAGOS, nvar_IAGOS)
     WRITE(*,*) 'O3:::::::::', listvar_IAGOS(ivarO3)
  END IF
  ALLOCATE(rcount_min_var(nvar_IAGOS,nlat))
  DO ivar=1,nvar_IAGOS
     time_ratio = REAL(listvar_nmonths(ivar))/REAL(nmonths)
     rcount_min_var(ivar,:) = lat_ratio(:)*time_ratio*rcount_min_ref
     IF (TRIM(listvar_IAGOS(ivar)) .EQ. 'NOy') THEN
        rcount_min_var(ivar,:) = rcount_min_var(ivar,:)/6 ! 1/6 is the ratio used in Cohen et al. (2023) to account for the lower sampling frequency of NOy.
     END IF
  END DO
  il_err = NF90_INQ_DIMID (il_ncfile, TRIM(lev_dim_name), il_ncdim) ! vertical grid level
  il_err = NF90_INQUIRE_DIMENSION (il_ncfile, il_ncdim, LEN = nlev)
  ALLOCATE(var_4D(nvar_IAGOS_with_stats,nlon,nlat,nlev_less))
  ALLOCATE(count_N_gridcells(nvar_IAGOS,nlon,nlat))
  WRITE(*,*) 'Starting the averaging over the column'
  count_N_gridcells(:,:,:)=0
  DO ivar=1,nvar_IAGOS
     ! Indexes for the amount of data/flights in the coupled list.
     ivar_stat_N        =ivar+(istat_N-1)        *nvar_IAGOS
     ivar_stat_N_flights=ivar+(istat_N_flights-1)*nvar_IAGOS
     ALLOCATE(var_3D_all_lev(nlon,nlat,nlev))
     il_err = NF90_INQ_VARID (il_ncfile, TRIM(listvar_IAGOS_with_stats(ivar_stat_N)), il_ncvar)
     il_err = NF90_GET_VAR (il_ncfile, il_ncvar, var_3D_all_lev(:,:,:))
     var_4D(ivar_stat_N,:,:,:) = var_3D_all_lev(:,:,l_top:l_down)
     DEALLOCATE(var_3D_all_lev)
     ALLOCATE(mask_3D(nlon,nlat,nlev_less))
     mask_3D(:,:,:)=0
     DO ilat=1,nlat
        WHERE (var_4D(ivar_stat_N,:,ilat,:) >= rcount_min_var(ivar,ilat))
           mask_3D(:,ilat,:)=1
        END WHERE
     END DO
     ! This condition ensures that a given average is representative of 2 vertical grid levels at least.
     DO ilon=1,nlon
        DO ilat=1,nlat
           count_N_gridcells(ivar,ilon,ilat)=SUM(mask_3D(ilon,ilat,:))
        END DO
     END DO
     ! WRITE(*,*) listvar_IAGOS_with_stats(ivar_stat_N), meanok(var_4D(ivar_stat_N,:,1:20,:),nlev_less*nlon*20,0.)
     DO istat = 1,(istat_N-1) ! No need to read the N_obs stat again.
        ivar_stat=ivar+(istat-1)*nvar_IAGOS
        il_err = NF90_INQ_VARID (il_ncfile, TRIM(listvar_IAGOS_with_stats(ivar_stat)), il_ncvar)
        ALLOCATE(var_3D_all_lev(nlon,nlat,nlev))
        il_err = NF90_GET_VAR (il_ncfile, il_ncvar, var_3D_all_lev(:,:,:))
        var_4D(ivar_stat,:,:,:) = var_3D_all_lev(:,:,l_top:l_down)
        DEALLOCATE(var_3D_all_lev)
        ! WRITE(*,*) listvar_IAGOS_with_stats(ivar_stat), meanok(var_4D(ivar_stat,:,:,:),nlev_less*nlon*nlat,-200.)
        ! Copying the units attribute.
        il_err = NF90_GET_ATT_STRING (il_ncfile, il_ncvar, 'units', cunit)
        listvar_units_IAGOS_with_stats(ivar_stat) = TRIM(cunit)
        ! Filtering the grid cells with respect to the sampling level.
        WHERE(mask_3D(:,:,:) == 0)
           var_4D(ivar_stat,:,:,:)=rwrong_val
        END WHERE
     END DO
     WHERE(mask_3D(:,:,:) == 0)
        var_4D(ivar_stat_N,:,:,:)=0
     END WHERE
     DEALLOCATE(mask_3D)
  END DO
  
  ALLOCATE(var_4D_model(nvar_model_with_stats,nlon,nlat,nlev_less))
  il_err = NF90_OPEN (TRIM(input_file_model), NF90_NOWRITE, il_ncfile)
  DO ivar=1,nvar_model
     ivar_mask=1
!!! Finds the IAGOS masking variable. If the current one is not amongst the IAGOS ones,
!!! or does not have a relevant observed peer to be filtered with, then we use ozone.
     IF (ANY(listvar_IAGOS(:) .EQ. listvar_model(ivar))) THEN
        ivar_mask=which_one_str(listvar_model(ivar), listvar_IAGOS, nvar_IAGOS)
     ELSE
        ivar_mask=which_one_str(c_O3, listvar_IAGOS, nvar_IAGOS)
     END IF
     IF (ANY(listvar_IAGOS(:) .EQ. listvar_model(ivar))) THEN
        ivar_mask=which_one_str(listvar_model(ivar), listvar_IAGOS, nvar_IAGOS)
     ELSE IF ((TRIM(listvar_model(ivar)) .EQ. TRIM(c_NOy_ACACIA)) .OR. &
          (TRIM(listvar_model(ivar)) .EQ. TRIM(c_NO_mask_NOy)) .OR. &
          (TRIM(listvar_model(ivar)) .EQ. TRIM(c_NO2_mask_NOy)) .OR. &
          (TRIM(listvar_model(ivar)) .EQ. TRIM(c_NOx_mask_NOy)) .OR. &
          (TRIM(listvar_model(ivar)) .EQ. TRIM(c_HNO3_mask_NOy)) .OR. &
          (TRIM(listvar_model(ivar)) .EQ. TRIM(c_PAN_mask_NOy)) .OR. &
          (TRIM(listvar_model(ivar)) .EQ. TRIM(c_O3_mask_NOy)) &
          .AND. ANY(listvar_IAGOS==c_NOy)) THEN
        ivar_mask=which_one_str(c_NOy, listvar_IAGOS, nvar_IAGOS)
        WRITE(*,*) 'listvar_IAGOS(ivar_mask), listvar_model(ivar)=',listvar_IAGOS(ivar_mask), listvar_model(ivar)
     ELSE
        ivar_mask=which_one_str(c_O3, listvar_IAGOS, nvar_IAGOS)
     END IF
     ivar_stat_N=ivar_mask+(istat_N-1)*nvar_IAGOS ! Index for the amount of data in the (obs) coupled list.
     ALLOCATE(mask_3D(nlon,nlat,nlev_less))
     mask_3D(:,:,:)=0
     DO ilat=1,nlat
        WHERE (var_4D(ivar_stat_N,:,ilat,:) >= rcount_min_var(ivar_mask,ilat)+0.3)
           mask_3D(:,ilat,:)=1
        END WHERE
     END DO
     DO istat = 1,nstats_model
        ivar_stat=ivar+(istat-1)*nvar_model
        ALLOCATE(var_3D_model_all_lev(nlon,nlat,nlev))
        il_err = NF90_INQ_VARID (il_ncfile, TRIM(listvar_model_with_stats(ivar_stat)), il_ncvar)
        il_err = NF90_GET_VAR (il_ncfile, il_ncvar, var_3D_model_all_lev(:,:,:))
        var_4D_model(ivar_stat,:,:,:)=var_3D_model_all_lev(:,:,l_top:l_down)
        DEALLOCATE(var_3D_model_all_lev)
        ! Copying the units attribute.
        il_err = NF90_GET_ATT_STRING (il_ncfile, il_ncvar, 'units', cunit)
        listvar_units_model_with_stats(ivar_stat) = TRIM(cunit)
        ! Filtering the grid cells with respect to the sampling level.
        WHERE(mask_3D(:,:,:) == 0)
           var_4D_model(ivar_stat,:,:,:)=rwrong_val
        END WHERE
     END DO
     DEALLOCATE(mask_3D)
  END DO
  il_err=NF90_CLOSE(il_ncfile)
  ALLOCATE(outvar_3D(nvar_IAGOS_with_stats,nlon,nlat))
  DO ivar=1,nvar_IAGOS
     DO istat = 1,(nstats_obs_column - nstats_N_sth) ! Special treatment for the N_obs and N_flights stats (sum instead of mean value).
        ivar_stat=ivar+(istat-1)*nvar_IAGOS
        DO ilon=1,nlon
           DO ilat=1,nlat
              outvar_3D(ivar_stat,ilon,ilat)=meanok(var_4D(ivar_stat,ilon,ilat,:),nlev_less,-200.)
           END DO
        END DO
     END DO
     DO istat = (nstats_obs_column - nstats_N_sth + 1), nstats_obs_column
        ivar_stat=ivar+(istat-1)*nvar_IAGOS
        DO ilon=1,nlon
           DO ilat=1,nlat
              outvar_3D(ivar_stat,ilon,ilat) = sumok(var_4D(ivar_stat,ilon,ilat,:),nlev_less,0.)
           END DO
        END DO
     END DO
  END DO
  ALLOCATE(outvar_3D_all_stats(nvar_with_all_stats,nlon,nlat))
  outvar_3D_all_stats(1:nvar_IAGOS_with_stats,:,:)=outvar_3D
  DO ivar=1,nvar_IAGOS
     ivar_stat=ivar+(nstats_obs_column + istat_count_Ncells - 1)*nvar_IAGOS
     DO ilat=1,nlat
        DO ilon=1,nlon
           outvar_3D_all_stats(ivar_stat,ilon,ilat)=count_N_gridcells(ivar,ilon,ilat)
        END DO
     END DO
  END DO
  DEALLOCATE(outvar_3D)
  ALLOCATE(outvar_3D_model(nvar_model_with_stats,nlon,nlat))
  DO ivar=1,nvar_model
     ivar_stat_N=ivar_mask+(istat_N-1)*nvar_IAGOS ! Index for the amount of data in the (obs) coupled list.
     DO istat = 1,nstats_model
        ivar_stat=ivar+(istat-1)*nvar_model
        DO ilon=1,nlon
           DO ilat=1,nlat
              outvar_3D_model(ivar_stat,ilon,ilat)=MEANOK(var_4D_model(ivar_stat,ilon,ilat,:),nlev_less,-200.)
           END DO
        END DO
     END DO
  END DO
  DEALLOCATE(var_4D,var_4D_model)
  ! Adding the units for the counting variables.
  DO ivar=1,nvar_IAGOS
     DO istat=1,nstats_count
        ivar_stat=ivar+(nstats_obs_column + istat - 1)*nvar_IAGOS
        IF (istat .EQ. istat_count_Ncells) THEN
           listvar_units_IAGOS_with_stats(ivar_stat)(1:LEN_TRIM('grid cells'))='grid cells'
           listvar_units_IAGOS_with_stats(ivar_stat)((LEN_TRIM('grid cells')+1):30)=' '
        ELSE
           WRITE(*,*) 'Unit: unknown counting stat. Exit program.'
        END IF
     END DO
  END DO
  il_err=NF90_CLOSE(il_ncfile)
  output_dim_names(1)=lon_dim_name
  output_dim_names(2)=lat_dim_name
  output_file_IAGOS=TRIM(outdir_time_period)//'IAGOS_'//TRIM(layer)//'_'//TRIM(season)//&
       '_'//cyyyymm1//'_'//cyyyymm2//'.nc'
  output_file_model=TRIM(outdir_time_period)//TRIM(experiment)//'_'//TRIM(layer)//'_'//TRIM(season)//&
       '_'//cyyyymm1//'_'//cyyyymm2//'.nc'
  
  WRITE(*,*) '---------------- Call tonetcdf_climato ------------------'
  CALL TONETCDF_CLIMATO(outvar_3D_all_stats(:,:,:), nlon, nlat, output_dim_names(:),&
       nvar_with_all_stats, listvar_IAGOS_with_stats(:),                     &
       listvar_units_IAGOS_with_stats(:),                                    &
       lon(:), lat(:), rwrong_val, output_file_IAGOS)
  CALL TONETCDF_CLIMATO(outvar_3D_model(:,:,:), nlon, nlat, output_dim_names(:), &
       nvar_model_with_stats, listvar_model_with_stats(:),                   &
       listvar_units_model_with_stats(:),                                    &
       lon(:), lat(:), rwrong_val, output_file_model)
  DEALLOCATE(outvar_3D_all_stats, outvar_3D_model)
END PROGRAM COLUMN_AVERAGE_IAGOS
