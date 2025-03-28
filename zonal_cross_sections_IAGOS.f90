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

PROGRAM ZONAL_CROSS_SECTIONS_IAGOS
  USE GLOBAL_VAR_INTERPOL_MOD
  USE FUNCTIONS_MOD
  USE NETCDF
  IMPLICIT NONE
  
!!!!!! Parameters !!!!!!
  ! Stats for counts.
  ! istat_supp_Ncells refers to the amount of data cells involved in a zonal mean.
  ! istat_supp_D1, Q1, Q3 and D9 correspond to metrics related to the gaussian distribution (deciles and quartiles).
  INTEGER,PARAMETER :: istat_supp_D1=1
  INTEGER,PARAMETER :: istat_supp_Q1=2
  INTEGER,PARAMETER :: istat_supp_Q3=3
  INTEGER,PARAMETER :: istat_supp_D9=4
  INTEGER,PARAMETER :: istat_supp_Ncells=5
  INTEGER,PARAMETER :: nstats_supp=5
  ! List of NEW counting statistics (thus excluding the amount of measurements).
  CHARACTER*30 list_stats_supp(nstats_supp)
  ! Amount of PREVIOUS counting statistics.
  INTEGER,PARAMETER :: nstats_N_sth=2
  CHARACTER*30,PARAMETER :: c_Decile_1='Decile_1', c_Quartile_1='Quartile_1', &
       c_Quartile_3='Quartile_3', c_Decile_9='Decile_9'
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
  ! Amount of observations for the current variable, mean value and several quantiles.
  INTEGER ivar_stat_N, ivar_stat_mean, ivar_stat_D1, ivar_stat_Q1, ivar_stat_Q3, ivar_stat_D9
  ! The amount of vertical grid levels to take into account.
  INTEGER nlev_less
  ! The amount of grid cells in the current longitude band.
  INTEGER nlon_domain
  ! The amount of all output statistics (* obs variables), 
  ! combining the input and the ones added from this routine.
  INTEGER nvar_with_all_stats, nvar_model_with_all_stats
  ! Variables and their units.
  CHARACTER*30 cvar, cunit
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
  ! Boundaries for the zonal band
  REAL rlon1, rlon2
  CHARACTER*30 clon1, clon2
  INTEGER ilon1, ilon2, nlon1, nlon2
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
  ! listvar_supp is the coupling between variable and new counting metrics.
  CHARACTER*30,ALLOCATABLE :: listvar_IAGOS(:), listvar_model(:)
  ! Same organisation for the units, since they depend on the variable and the metric.
  ! For example, "ppb" is for ozone if the metric is "Mean", "SD", "Min", etc. but
  ! "data" if the metric is "N_obs".
  CHARACTER*30,ALLOCATABLE :: listvar_IAGOS_with_stats(:), listvar_units_IAGOS_with_stats(:)
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
  REAL, ALLOCATABLE :: outvar_2D(:,:), outvar_2D_model(:,:)
  REAL, ALLOCATABLE :: outvar_2D_all_stats(:,:), outvar_2D_model_all_stats(:,:)
  logical, ALLOCATABLE :: mask_3D(:,:,:)
  ! Amount of gridcells involved in each column average.
  INTEGER, ALLOCATABLE :: count_N_gridcells(:,:,:)
  REAL, ALLOCATABLE :: var_D1_2D(:,:), var_Q1_2D(:,:), var_Q3_2D(:,:), var_D9_2D(:,:)
  REAL, ALLOCATABLE :: var_D1_2D_model(:,:), var_Q1_2D_model(:,:)
  REAL, ALLOCATABLE :: var_Q3_2D_model(:,:), var_D9_2D_model(:,:)
  REAL :: a(1001) = (/(i, i=2,2002, 2)/)

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
  list_stats_supp(1)='Decile_1'
  list_stats_supp(2)='Quartile_1'
  list_stats_supp(3)='Quartile_3'
  list_stats_supp(4)='Decile_9'
  list_stats_supp(5)='N_gridcells'
  ! Definition for the istat indexes. See function which_one_str defined in functions_mod.f90.
  ! Only concerns the input stats, not the supplementary ones (already defined as parameters).
  istat_mean  =which_one_str(c_Mean , list_stats_obs, nstats_obs_column)
  istat_SD    =which_one_str(c_SD   , list_stats_obs, nstats_obs_column)
  istat_min   =which_one_str(c_Min  , list_stats_obs, nstats_obs_column)
  istat_max   =which_one_str(c_Max  , list_stats_obs, nstats_obs_column)
  istat_N_diff_updown=which_one_str(c_N_diff_updown, list_stats_obs, nstats_obs_column)
  istat_Pressure=which_one_str(c_Pressure, list_stats_obs, nstats_obs_column)
  istat_N        =which_one_str(c_N_obs, list_stats_obs, nstats_obs_column)
  istat_N_flights=which_one_str(c_N_flights, list_stats_obs, nstats_obs_column)
  IF (istat_N .NE. nstats_obs_column) THEN
     WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     WRITE(*,*) '!!!!!! Careful: istat_N is not the last stat index,          !!!!!!'
     WRITE(*,*) '!!!!!! whereas a part of the program assumes istat_N=nstats. !!!!!!'
     WRITE(*,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
  END IF
  ! Same for the model.
  istat_mean_model=which_one_str(c_Mean_model, list_stats_model, nstats_model)
  istat_SD_model=  which_one_str(c_SD_model  , list_stats_model, nstats_model)
  istat_min_model= which_one_str(c_Min_model , list_stats_model, nstats_model)
  istat_max_model= which_one_str(c_Max_model , list_stats_model, nstats_model) 
  OPEN(14,FILE=TRIM(pdir)//'namelist_zonal',FORM='formatted')
  READ(14,*) experiment, cmodel
  READ(14,*) cyyyymm1, cyyyymm2
  READ(14,*) l_top_UTLS, l_down_UTLS
  READ(14,*) l_top_FT, l_down_FT
  READ(14,*) rlon1, rlon2
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
  close(14)
  WRITE(*,*) '***************************************'
  WRITE(*,*) TRIM(experiment)
  WRITE(*,*) cyyyymm1//' '//cyyyymm2
  WRITE(*,*) l_top_UTLS, l_down_UTLS
  WRITE(*,*) l_top_FT, l_down_FT
  WRITE(*,*) rlon1, rlon2
  WRITE(*,*) lon_dim_name, lat_dim_name, lev_dim_name
  WRITE(*,*) nlat
  WRITE(*,*) nvar_IAGOS
  WRITE(*,*) listvar_IAGOS(1:nvar_IAGOS)
  WRITE(*,*) listvar_nmonths(1:nvar_IAGOS)
  WRITE(*,*) nmonths
  WRITE(*,*) nvar_model
  WRITE(*,*) listvar_model(1:nvar_model)
  WRITE(*,*) rcount_min_ref
  WRITE(*,*) tmpdir_there
  WRITE(*,*) season
  WRITE(*,*) layer
  WRITE(*,*) dir_level_by_level
  WRITE(*,*) outdir_time_period
  WRITE(*,*) input_file_IAGOS
  WRITE(*,*) input_file_model
  WRITE(*,*) '***************************************'
  READ(cyyyymm1,'(i6)') iyyyymm1
  READ(cyyyymm2,'(i6)') iyyyymm2
  open(unitFile_check, FILE=TRIM(tmpdir_there)//'check_average_iagos_'// &
       TRIM(layer)//'_'//TRIM(season), FORM='formatted')
  WRITE(unitFile_check,*) 'namelist_zonal file content:'
  WRITE(unitFile_check,*) 'experiment=', experiment
  WRITE(unitFile_check,*) 'iyyyymm1, iyyyymm2=', iyyyymm1, iyyyymm2
  WRITE(unitFile_check,*) 'l_top_UTLS, l_down_UTLS=', l_top_UTLS, l_down_UTLS
  WRITE(unitFile_check,*) 'l_top_FT, l_down_FT=', l_top_FT, l_down_FT
  WRITE(unitFile_check,*) 'rlon1, rlon2=', rlon1, rlon2
  WRITE(unitFile_check,*) 'lon_dim_name, lat_dim_name, lev_dim_name=', lon_dim_name, lat_dim_name, lev_dim_name
  WRITE(unitFile_check,*) 'nlat=',nlat
  WRITE(unitFile_check,*) 'nvar_IAGOS=', nvar_IAGOS
  WRITE(unitFile_check,*) 'listvar_IAGOS(1:nvar_IAGOS)= ', listvar_IAGOS(1:nvar_IAGOS)
  WRITE(unitFile_check,*) 'nvar_model=', nvar_model
  WRITE(unitFile_check,*) 'listvar_model(1:nvar_model)=', listvar_model(1:nvar_model)
  WRITE(unitFile_check,*) 'rcount_min_ref=', rcount_min_ref
  WRITE(unitFile_check,*) 'tmpdir=', tmpdir
  WRITE(unitFile_check,*) 'tmpdir_there=', tmpdir_there
  WRITE(unitFile_check,*) 'season=', season
  WRITE(unitFile_check,*) 'layer=', layer
  WRITE(unitFile_check,*) 'dir_level_by_level=', dir_level_by_level
  WRITE(unitFile_check,*) 'outdir_time_period=', outdir_time_period
  WRITE(unitFile_check,*) 'input_file_IAGOS=', input_file_IAGOS
  WRITE(unitFile_check,*) 'input_file_model=', input_file_model
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
  ! Defining the whole sets of couples variable x stat.
  nvar_IAGOS_with_stats=nvar_IAGOS*nstats_obs_column
  nvar_model_with_stats=nvar_model*nstats_model
  nvar_with_all_stats=(nstats_obs_column + nstats_supp)*nvar_IAGOS
  nvar_model_with_all_stats=(nstats_model + nstats_supp)*nvar_model
  ALLOCATE(listvar_IAGOS_with_stats(nvar_with_all_stats))
  ALLOCATE(listvar_model_with_stats(nvar_model_with_all_stats))
  ALLOCATE(listvar_units_IAGOS_with_stats(nvar_with_all_stats))
  ALLOCATE(listvar_units_model_with_stats(nvar_model_with_all_stats))
  DO ivar=1,nvar_IAGOS
     DO istat=1,nstats_obs_column
        ivar_stat=ivar+(istat-1)*nvar_IAGOS
        listvar_IAGOS_with_stats(ivar_stat)=TRIM(listvar_IAGOS(ivar))//'_'//TRIM(list_stats_obs(istat))
     END DO
     DO istat=1,nstats_supp
        ivar_stat=ivar+(nstats_obs_column + istat - 1)*nvar_IAGOS
        listvar_IAGOS_with_stats(ivar_stat)=TRIM(listvar_IAGOS(ivar))//'_'//TRIM(list_stats_supp(istat))
     END DO
  END DO
  WRITE(unitFile_check,*) 'listvar_IAGOS_with_stats=', listvar_IAGOS_with_stats
  DO ivar=1,nvar_model
     DO istat=1,nstats_model
        ivar_stat=ivar+(istat-1)*nvar_model
        listvar_model_with_stats(ivar_stat)=TRIM(listvar_model(ivar))//'_'//TRIM(list_stats_model(istat))
     END DO
     DO istat=1,nstats_supp ! Ok, here there is one useless stat for the models (N_cells), but it is easier this way.
        ivar_stat=ivar+(nstats_model + istat - 1)*nvar_model
        listvar_model_with_stats(ivar_stat)=TRIM(listvar_model(ivar))//'_'//TRIM(list_stats_supp(istat))
     END DO
  END DO
  WRITE(unitFile_check,*) 'listvar_model_with_stats=', listvar_model_with_stats
!!! 2/ Filling the data arrays.
!!! 2.1/ Reads the regional files, for each data set.
!!! 2.1.1/ We start by reading the IAGOS data, especially the amounts of data. 
!!! It will be used as a reference to pick up the gridcells we account for in the average computation.
  il_err = NF90_OPEN (TRIM(input_file_IAGOS), NF90_NOWRITE, il_ncfile)
  il_err = NF90_INQ_DIMID (il_ncfile, TRIM(lon_dim_name), il_ncdim) ! longitude
  WRITE(*,*) 'nlon before dimlen=', nlon, il_err
  il_err = NF90_INQUIRE_DIMENSION (il_ncfile, il_ncdim, LEN = nlon)
  WRITE(*,*) 'nlon after dimlen=', nlon, il_err
  ALLOCATE(lon(nlon))
  il_err = NF90_GET_VAR (il_ncfile, il_ncdim, lon(:))
  il_err = NF90_INQ_DIMID (il_ncfile, TRIM(lat_dim_name), il_ncdim) ! latitude
  WRITE(*,*) 'nlat before dimlen=', nlat, il_err
  il_err = NF90_INQUIRE_DIMENSION (il_ncfile, il_ncdim, LEN = nlat)
  WRITE(*,*) 'nlat after dimlen=', nlat
  WRITE(*,*) TRIM(input_file_IAGOS)
  ! Within the longitude axis, we locate the zonal band where we want to average the grid cells.
  ilon1=1 ; ilon2=nlon
  DO WHILE (lon(ilon1) <= rlon1)
     ilon1=ilon1+1
  END DO
  DO WHILE (lon(ilon2) >= rlon2)
     ilon2=ilon2-1
  END DO
  nlon_domain=ilon2-ilon1+1
  WRITE(unitFile_check,*) 'rlon1, rlon2=', rlon1, rlon2
  WRITE(unitFile_check,*) 'nlon_domain, ilon1, ilon2, lon(ilon1), lon(ilon2)=', &
       nlon_domain, ilon1, ilon2, lon(ilon1), lon(ilon2)
  ALLOCATE(lat(nlat))
  il_err = NF90_GET_VAR(il_ncfile, il_ncdim, lat(:))
  ALLOCATE(lat_ratio(nlat))
  CALL LATITUDE_SURFACE_RATIO(nlat, lat(:), lat_ratio(:))
  IF (TRIM(season) == 'ANN') then
     rcount_min_ref = rcount_min_ref*4
  ELSE
     ! In case there are not always 3 months per season.
     ! Just be careful and call a season by the months initials.
     rcount_min_ref = rcount_min_ref*REAL(LEN_TRIM(season))/3
  END IF
  ! Ozone (we expect the ozone column to be contained in every file)
  IF (ANY(listvar_IAGOS=='O3')) THEN
     ivarO3= which_one_str(c_O3, listvar_IAGOS, nvar_IAGOS)
     WRITE(unitFile_check,*) 'O3:::::::::', listvar_IAGOS(ivarO3)
  END IF
  ALLOCATE(rcount_min_var(nvar_IAGOS,nlat))
  DO ivar=1,nvar_IAGOS
     time_ratio = REAL(listvar_nmonths(ivar))/REAL(nmonths)
     rcount_min_var(ivar,:) = lat_ratio(:)*time_ratio*rcount_min_ref
     IF (TRIM(listvar_IAGOS(ivar)) == 'NOy') THEN
        rcount_min_var(ivar,:) = rcount_min_var(ivar,:)/6 ! 1/6 is the ratio used in Cohen et al. (2023) to account for the lower sampling frequency of NOy.
     END IF
  END DO
  il_err = NF90_INQ_DIMID (il_ncfile, TRIM(lev_dim_name), il_ncdim) ! vertical grid level
  il_err = NF90_INQUIRE_DIMENSION (il_ncfile, il_ncdim, LEN = nlev)
  ALLOCATE(var_4D(nvar_IAGOS_with_stats,nlon,nlat,nlev_less))
  ALLOCATE(count_N_gridcells(nvar_IAGOS,nlon,nlat))
  
  count_N_gridcells(:,:,:)=0
  DO ivar=1,nvar_IAGOS
     ivar_stat_N=ivar+(istat_N-1)*nvar_IAGOS ! Index for the amount of data in the coupled list.
     ALLOCATE(var_3D_all_lev(nlon,nlat,nlev))
     il_err = NF90_INQ_VARID (il_ncfile, TRIM(listvar_IAGOS_with_stats(ivar_stat_N)), il_ncvar)
     il_err = NF90_GET_VAR (il_ncfile, il_ncvar, var_3D_all_lev(:,:,:))
     var_4D(ivar_stat_N,:,:,:)=var_3D_all_lev(:,:,l_top:l_down)
     DEALLOCATE(var_3D_all_lev)
     ALLOCATE(mask_3D(nlon,nlat,nlev_less))
     mask_3D(:,:,:)=.TRUE.
     DO ilat=1,nlat
        WHERE (var_4D(ivar_stat_N,:,ilat,:) < rcount_min_var(ivar,ilat)+0.3) ! 0.3 is a margin for the gridcells at the poles.
           mask_3D(:,ilat,:)=.FALSE.
        END WHERE
     END DO
     ! How many grid cells represent this horizontal grid cell. If only one, then it will be considered insufficient.
     DO ilat=1,nlat
        DO ilon=1,nlon
           count_N_gridcells(ivar,ilon,ilat)=COUNT(mask_3D(ilon,ilat,:))
     ! The following condition ensures that a given average is representative of 2 vertical grid levels at least.
           IF (count_N_gridcells(ivar,ilon,ilat) < 2) THEN
              mask_3D(ilon,ilat,:)=.FALSE.
           END IF
        END DO
     END DO
     DO istat = 1,(istat_N-1) ! No need to read the N_obs stat again.
        ivar_stat=ivar+(istat-1)*nvar_IAGOS
        il_err = NF90_INQ_VARID (il_ncfile, TRIM(listvar_IAGOS_with_stats(ivar_stat)), il_ncvar)
        ALLOCATE(var_3D_all_lev(nlon,nlat,nlev))
        il_err = NF90_GET_VAR (il_ncfile, il_ncvar,var_3D_all_lev(:,:,:))
        var_4D(ivar_stat,:,:,:)=var_3D_all_lev(:,:,l_top:l_down)
        DEALLOCATE(var_3D_all_lev)
        WRITE(unitFile_check,*) listvar_IAGOS_with_stats(ivar_stat), meanok(var_4D(ivar_stat,:,:,:),nlev_less*nlon*nlat,-200.)
        ! Copying the units attribute.
        il_err = NF90_GET_ATT_STRING (il_ncfile, il_ncvar, 'units', cunit)
        listvar_units_IAGOS_with_stats(ivar_stat) = TRIM(cunit)
     END DO
     ! Filtering the grid cells with respect to the sampling level.
     DO istat = 1,(nstats_obs_column - nstats_N_sth)
        ivar_stat = ivar+(istat-1)*nvar_IAGOS
        WHERE(.NOT. mask_3D(:,:,:))
           var_4D(ivar_stat,:,:,:)=rwrong_val
        END WHERE
     END DO
     ! Last: looping the loop. The non-selected grid cells are set to a 0-sampling.
     DO istat = (nstats_obs_column - nstats_N_sth + 1), nstats_obs_column
        ivar_stat = ivar+(istat-1)*nvar_IAGOS
        WHERE (.NOT. mask_3D(:,:,:))
           var_4D(ivar_stat,:,:,:)=0.
        END WHERE
     END DO
     DEALLOCATE(mask_3D)
  END DO
  WRITE(*,*) 'Checking some mean values:', meanok(var_4D(1,ilon1:ilon2,:,:),nlon_domain*nlat*nlev_less,-200.)
  il_err = NF90_CLOSE(il_ncfile)
  
  ALLOCATE(var_4D_model(nvar_model_with_stats,nlon,nlat,nlev_less))
  il_err = NF90_OPEN (TRIM(input_file_model), NF90_NOWRITE, il_ncfile)
  DO ivar=1,nvar_model
     ivar_mask=1
!!! Finds the IAGOS masking variable. If the current one is not amongst the IAGOS ones,
!!! or does not have a relevant observed peer to be filtered with, then we use ozone.
     IF (ANY(listvar_IAGOS(:) == listvar_model(ivar))) THEN
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
     ! Repeating the mask definition for the modelled fields.
     mask_3D(:,:,:)=.FALSE.
     DO ilev=1,nlev_less
        DO ilat=1,nlat
           WHERE (var_4D(ivar_stat_N,:,ilat,ilev) >= rcount_min_var(ivar_mask,ilat)+0.3 &
                .AND. count_N_gridcells(ivar_mask,:,ilat) >= 2)
              mask_3D(:,ilat,ilev)=.TRUE.
           END WHERE
        END DO
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
        WHERE(mask_3D(:,:,:) .EQV. .FALSE.)
           var_4D_model(ivar_stat,:,:,:)=rwrong_val
        END WHERE
     END DO
     DEALLOCATE(mask_3D)
  END DO
! Now we can compute the output meridional profiles.
  ALLOCATE(outvar_2D(nvar_IAGOS_with_stats,nlat))
  DO ivar=1,nvar_IAGOS
     DO istat = 1,(nstats_obs_column - nstats_N_sth) ! Special treatment for the N_obs and N_flights stats (sum instead of mean value).
        ivar_stat=ivar+(istat-1)*nvar_IAGOS
        DO ilat=1,nlat
           outvar_2D(ivar_stat,ilat)=meanok(var_4D(ivar_stat,ilon1:ilon2,ilat,:),nlon_domain*nlev_less,-200.)
        END DO
     END DO
     DO istat = (nstats_obs_column - nstats_N_sth + 1), nstats_obs_column
        ivar_stat=ivar+(istat-1)*nvar_IAGOS
        DO ilat=1,nlat
           outvar_2D(ivar_stat,ilat) = sumok(var_4D(ivar_stat,ilon1:ilon2,ilat,:),nlon_domain*nlev_less,0.)
        END DO
     END DO
  END DO
  ! Integrating the supplementary statistics to the output arrays.
  ALLOCATE(var_D1_2D(nvar_IAGOS,nlat))
  ALLOCATE(var_Q1_2D(nvar_IAGOS,nlat))
  ALLOCATE(var_Q3_2D(nvar_IAGOS,nlat))
  ALLOCATE(var_D9_2D(nvar_IAGOS,nlat))
  ALLOCATE(outvar_2D_all_stats(nvar_with_all_stats,nlat))
  outvar_2D_all_stats(1:nvar_IAGOS_with_stats,:)=outvar_2D
  DO ivar=1,nvar_IAGOS
     ivar_stat=ivar+(nstats_obs_column + istat_supp_Ncells - 1)*nvar_IAGOS
     DO ilat=1,nlat
        outvar_2D_all_stats(ivar_stat,ilat) = SUM( count_N_gridcells(ivar,ilon1:ilon2,ilat), &
             count_N_gridcells(ivar,ilon1:ilon2,ilat) >= 2 )
     END DO
     ! Now addind the main quantiles.
     ivar_stat_mean=ivar+(istat_mean-1)*nvar_IAGOS
     DO ilat=1,nlat
        var_D1_2D(ivar,ilat)=quantileok(var_4D(ivar_stat_mean,ilon1:ilon2,ilat,:), &
             nlon_domain*nlev_less, -200., 0.1)
        var_Q1_2D(ivar,ilat)=quantileok(var_4D(ivar_stat_mean,ilon1:ilon2,ilat,:), &
             nlon_domain*nlev_less, -200., 0.25)
        var_Q3_2D(ivar,ilat)=quantileok(var_4D(ivar_stat_mean,ilon1:ilon2,ilat,:), &
             nlon_domain*nlev_less, -200., 0.75)
        var_D9_2D(ivar,ilat)=quantileok(var_4D(ivar_stat_mean,ilon1:ilon2,ilat,:), &
             nlon_domain*nlev_less, -200., 0.9)
     END DO
     ivar_stat_D1 = ivar + (istat_supp_D1 + nstats_obs_column - 1)*nvar_IAGOS
     ivar_stat_Q1 = ivar + (istat_supp_Q1 + nstats_obs_column - 1)*nvar_IAGOS
     ivar_stat_Q3 = ivar + (istat_supp_Q3 + nstats_obs_column - 1)*nvar_IAGOS
     ivar_stat_D9 = ivar + (istat_supp_D9 + nstats_obs_column - 1)*nvar_IAGOS
     outvar_2D_all_stats(ivar_stat_D1,:)=var_D1_2D(ivar,:)
     outvar_2D_all_stats(ivar_stat_Q1,:)=var_Q1_2D(ivar,:)
     outvar_2D_all_stats(ivar_stat_Q3,:)=var_Q3_2D(ivar,:)
     outvar_2D_all_stats(ivar_stat_D9,:)=var_D9_2D(ivar,:)
  END DO
  outvar_2D_all_stats(1:nvar_IAGOS_with_stats,:)=outvar_2D(:,:)
  DEALLOCATE(outvar_2D)
  ALLOCATE(outvar_2D_model(nvar_model_with_stats,nlat))
  ALLOCATE(outvar_2D_model_all_stats(nvar_model_with_all_stats,nlat))
  ALLOCATE(var_D1_2D_model(nvar_model,nlat))
  ALLOCATE(var_Q1_2D_model(nvar_model,nlat))
  ALLOCATE(var_Q3_2D_model(nvar_model,nlat))
  ALLOCATE(var_D9_2D_model(nvar_model,nlat))
  DO ivar=1,nvar_model
     ivar_mask=1
!!! Finds the IAGOS masking variable. If the current one is not amongst the IAGOS ones,
!!! or does not have a relevant observed peer to be filtered with, then we use ozone.
     IF (ANY(listvar_IAGOS(:) == listvar_model(ivar))) THEN
        ivar_mask=which_one_str(listvar_model(ivar), listvar_IAGOS, nvar_IAGOS)
     ELSE
        ivar_mask=which_one_str(c_O3, listvar_IAGOS, nvar_IAGOS)
     END IF
     ivar_stat_N=ivar_mask+(istat_N-1)*nvar_IAGOS ! Index for the amount of data in the (obs) coupled list.
     DO istat = 1,nstats_model
        ivar_stat=ivar+(istat-1)*nvar_model
        DO ilat=1,nlat
           outvar_2D_model(ivar_stat,ilat)=meanok(var_4D_model(ivar_stat,ilon1:ilon2,ilat,:), &
                nlon_domain*nlev_less,-200.)
        END DO
     END DO
     ! Now addind the main quantiles.
     ivar_stat_mean=ivar+(istat_mean_model-1)*nvar_model
     DO ilat=1,nlat
        var_D1_2D_model(ivar,ilat)=quantileok(var_4D_model(ivar_stat_mean,ilon1:ilon2,ilat,:), &
             nlon_domain*nlev_less, -200., 0.1)
        var_Q1_2D_model(ivar,ilat)=quantileok(var_4D_model(ivar_stat_mean,ilon1:ilon2,ilat,:), &
             nlon_domain*nlev_less, -200., 0.25)
        var_Q3_2D_model(ivar,ilat)=quantileok(var_4D_model(ivar_stat_mean,ilon1:ilon2,ilat,:), &
             nlon_domain*nlev_less, -200., 0.75)
        var_D9_2D_model(ivar,ilat)=quantileok(var_4D_model(ivar_stat_mean,ilon1:ilon2,ilat,:), &
             nlon_domain*nlev_less, -200., 0.9)
     END DO
     ivar_stat_D1 = ivar + (istat_supp_D1 + nstats_model - 1)*nvar_model
     ivar_stat_Q1 = ivar + (istat_supp_Q1 + nstats_model - 1)*nvar_model
     ivar_stat_Q3 = ivar + (istat_supp_Q3 + nstats_model - 1)*nvar_model
     ivar_stat_D9 = ivar + (istat_supp_D9 + nstats_model - 1)*nvar_model
     outvar_2D_model_all_stats(ivar_stat_D1,:)=var_D1_2D_model(ivar,:)
     outvar_2D_model_all_stats(ivar_stat_Q1,:)=var_Q1_2D_model(ivar,:)
     outvar_2D_model_all_stats(ivar_stat_Q3,:)=var_Q3_2D_model(ivar,:)
     outvar_2D_model_all_stats(ivar_stat_D9,:)=var_D9_2D_model(ivar,:)
  END DO
  outvar_2D_model_all_stats(1:nvar_model_with_stats,:)=outvar_2D_model(:,:)
  DEALLOCATE(var_4D,var_4D_model)
  ! Adding the units for the supplementary variables.
  DO ivar=1,nvar_IAGOS
     ivar_stat_mean=ivar+(istat_mean - 1)*nvar_IAGOS
     DO istat=1,nstats_supp
        ivar_stat=ivar+(nstats_obs_column + istat - 1)*nvar_IAGOS
        IF (istat == istat_supp_Ncells) THEN
           listvar_units_IAGOS_with_stats(ivar_stat)(1:LEN_TRIM('grid cells'))='grid cells'
           listvar_units_IAGOS_with_stats(ivar_stat)((LEN_TRIM('grid cells')+1):30)=' '
        ELSE
           listvar_units_IAGOS_with_stats(ivar_stat)=listvar_units_IAGOS_with_stats(ivar_stat_mean)
        END IF
     END DO
  END DO
  DO ivar=1,nvar_model
     ivar_stat_mean=ivar+(istat_mean_model - 1)*nvar_model
     DO istat=1,nstats_supp
        ivar_stat=ivar+(nstats_model + istat - 1)*nvar_model
        IF (istat == istat_supp_Ncells) THEN
           listvar_units_model_with_stats(ivar_stat)(1:LEN_TRIM('grid cells'))='grid cells'
           listvar_units_model_with_stats(ivar_stat)((LEN_TRIM('grid cells')+1):30)=' '
        ELSE
           listvar_units_model_with_stats(ivar_stat)=listvar_units_model_with_stats(ivar_stat_mean)
        END IF
     END DO
  END DO
  il_err = NF90_CLOSE(il_ncfile)
  nlon1=INT(rlon1)
  nlon2=INT(rlon2)
  CALL INT2CHAR(nlon1,tmpdir,clon1)
  CALL INT2CHAR(nlon2,tmpdir,clon2)
  ! 3/ Writing the columns in the output file.
  output_file_IAGOS=TRIM(outdir_time_period)//'IAGOS_'//TRIM(layer)//'_'//TRIM(season)//&
       '_'//TRIM(clon1)//'E'//'_'//TRIM(clon2)//'E'//'_'//cyyyymm1//'_'//cyyyymm2//'.nc'
  output_file_model=TRIM(outdir_time_period)//TRIM(experiment)//'_'//TRIM(layer)//'_'//TRIM(season)//&
       '_'//TRIM(clon1)//'E'//'_'//TRIM(clon2)//'E'//'_'//cyyyymm1//'_'//cyyyymm2//'.nc'
  
  WRITE(*,*) '---------------- Call tonetcdf_zonal ------------------'
  WRITE(*,*) nlat
  WRITE(*,*) lat
  CALL TONETCDF_ZONAL(outvar_2D_all_stats(:,:), nlat, lat_dim_name, &
       nvar_with_all_stats, listvar_IAGOS_with_stats(:),            & ! The listvar arrays contain all the stats, no need to specify "all", unlike nvar and outvar_2D.
       listvar_units_IAGOS_with_stats(:),                           &
       rlon1, rlon2, lat(:), rwrong_val, output_file_IAGOS)
  CALL TONETCDF_ZONAL(outvar_2D_model_all_stats(:,:), nlat, lat_dim_name, &
       nvar_model_with_all_stats, listvar_model_with_stats(:),            &
       listvar_units_model_with_stats(:),                                 &
       rlon1, rlon2, lat(:), rwrong_val, output_file_model)
  DEALLOCATE(outvar_2D_all_stats, outvar_2D_model_all_stats)
END PROGRAM ZONAL_CROSS_SECTIONS_IAGOS
