SUBROUTINE MEAN_VALUES_INTERPOL(                                           &
!!! 1/ Constant parameters
!!! Observations
     filename_prefix_IAGOS, filename_suffix_IAGOS,                         &
     list_coord, list_variables, list_der_variables, list_der_variables_nc,&
     thres_val, IAGOS_conv_factor, header_packages, tab_var_names,         &
     list_packages,                                                        &
!!! Model
     list_dim_model, list_var_model, model_variables, header_models,       &
     model_units_output, model_conv_factor,                                &
!!! Layer-selection parameters.
     list_layers, PV_down_layers, PV_up_layers, P_down_layers, P_up_layers,&
     ozone_max_layers, ozone_min_layers, logic_PV, PV_tropopause,          &
!!! 2/ Monthly parameters
     iyyyymm, ntime, nfiles, list_data_files,                              &
!!! Grid definition
     rlon, rlat, dlon, dlat, shift_lon, shift_lat,                         &
!!! Input meteorological fields
     P_edge, P_mid, PV,                                                    &
!!! Daily outputs, temporary arrays useful for SD calculation.
     rcount_daily, var_model_daily,                                        &
!!! 5D-outputs, to be written in the output files.
     mean_5D, rcount, n_diff_updown, min_5D, max_5D, day1_5D, dayN_5D,     &
     N_flights_5D, P_5D,                                                   &
     mean_5D_model, min_5D_model, max_5D_model,                            &
     ilon_ok, klat_ok, llev_ok                                             &
     )
  USE GLOBAL_VAR_INTERPOL_MOD
  USE FUNCTIONS_MOD
  USE NETCDF
  ! USE FUNCTION NF90_GET_ATT_STRING
  IMPLICIT NONE
  ! Monthly parameters
  INTEGER,INTENT(IN):: iyyyymm, ntime, nfiles
  REAL   ,INTENT(IN):: thres_val(nvariables),                                                  &
       IAGOS_conv_factor(nvariables), model_conv_factor(nvar_model)
  ! Observations
  CHARACTER(LEN=50) ,INTENT(IN):: filename_prefix_IAGOS, filename_suffix_IAGOS
  CHARACTER(LEN=30) ,INTENT(IN):: list_coord(ncoord), list_variables(nvariables),              &
       list_der_variables(n_der_variables), list_der_variables_nc(n_der_variables),            &
       header_packages(npackages+2), tab_var_names(nvar_tot,npackages+2),                      &
       list_packages(npackages)
  CHARACTER(LEN=200),INTENT(IN):: list_data_files(nfiles)
  ! Model
  CHARACTER(LEN=30) ,INTENT(IN):: list_dim_model(ndim_model),                                  &
       list_var_model(nvar_model), model_variables(nvar_model),                                &
       header_models(n_models_registered+2), model_units_output(nvar_model_with_stats)
  ! Layers definition.
  CHARACTER(LEN=30) ,INTENT(IN):: list_layers(nlayers)
  LOGICAL           ,INTENT(IN):: logic_PV 
  REAL, INTENT(IN)  :: PV_down_layers(nlayers), PV_up_layers(nlayers), P_down_layers(nlayers), &
       P_up_layers(nlayers), ozone_max_layers(nlayers), ozone_min_layers(nlayers), PV_tropopause
  ! Grid definition.
  REAL, INTENT(IN)  :: rlon(nlon), rlat(nlat), dlon, dlat, shift_lon, shift_lat
  ! Input meteorological fields.
  REAL, INTENT(IN)  :: P_edge(nlon,nlat,nlev,ntime),                                           &
       P_mid(nlon,nlat,nlev,ntime),                                                            &
       PV(nlon,nlat,nlev,ntime)
  ! Daily outputs, temporary arrays useful for SD derivation.
  REAL, INTENT(OUT):: rcount_daily(nvariables,nlayers,nlon,nlat,nlev,ntime),                   &
       var_model_daily(nvar_model,nlon,nlat,nlev,ntime)
  ! 5D-outputs, to be written in the output files.
  REAL, INTENT(OUT):: mean_5D(nvariables,nlayers,nlon,nlat,nlev),                              &
       rcount(nvariables,nlayers,nlon,nlat,nlev),                                              &
       n_diff_updown(nvariables,nlayers,nlon,nlat,nlev),                                       &       
       min_5D(nvariables,nlayers,nlon,nlat,nlev),                                              &
       max_5D(nvariables,nlayers,nlon,nlat,nlev),                                              &
       N_flights_5D(nvariables,nlayers,nlon,nlat,nlev),                                        &
       P_5D(nvariables,nlayers,nlon,nlat,nlev),                                                &
       mean_5D_model(nvar_model,nlayers,nlon,nlat,nlev),                                       &
       min_5D_model(nvar_model,nlayers,nlon,nlat,nlev),                                        &
       max_5D_model(nvar_model,nlayers,nlon,nlat,nlev)
  INTEGER, INTENT(OUT):: day1_5D(nvariables,nlayers,nlon,nlat,nlev),                           &
       dayN_5D(nvariables,nlayers,nlon,nlat,nlev)
  INTEGER, INTENT(OUT):: ilon_ok(nlayers), klat_ok(nlayers), llev_ok(nlayers)
  ! Local variables.
  ! The integer i_sign is a tool, just for determining the sign of a number x with the SIGN function,
  ! without creating a supplementary loop. It works with real numbers too. n_1_or_0 works similarly.
  INTEGER i_sign,n_1_or_0
  ! N_min_daily is the threshold required on the daily amount of data, for filtering out the days with a very low sampling.
  ! If the time resolution is monthly, then it will be a monthly threshold.
  REAL N_min_daily
  ! i_sign real equivalent.
  REAL r_sign  ! Length of the character string. Used while reading packName global attribute.
  REAL rlon_eq ! Meant to facilitate the distance calculation between 179E and 179W, by replacing the true longitude by an imaginary one, 360° away.
  INTEGER char_len
  INTEGER xtype
  INTEGER ichar ! Character position into a string.
  LOGICAL, ALLOCATABLE :: bin_sampled_NOx(:)
  ! Just to check a bug only once during the current month.
  LOGICAL bin_this_mth
  ! Index pointing to a given IAGOS file in each monthly list of files.
  INTEGER iFile
  ! Dates as integers, then as characters.
  ! Their lengths have to be the same as the variable in INT2CHAR subroutine (here: 20).
  INTEGER iyyyy
  CHARACTER(LEN=30) cyyyymm, cyyyy
  ! Integers relative to the IAGOS files dimensions.
  ! nlines and nlines_quarter are the amount of rows in the current IAGOS file, and its rounded quarter.
  ! Its use only consists in checking data. Unless the flight is very short,
  ! nlines_quarter generally points at measurements made in the UTLS.
  ! ncol is the total amount of columns in the current IAGOS file. ncol_to_read is the amount of
  ! columns to be read in the current IAGOS file, depending on the user request and the
  ! data availability in the current IAGOS file.
  INTEGER nlines,nlines_quarter,ncol,ncol_to_read,nvariables_this_flight,nvar_NOx_this_flight,nvar_NOx
  INTEGER nlines_this_flight, nlines_until_midnight, nlines_after_midnight, nlines_after_midnight_concat
  ! Indexes repairing in the IAGOS files, and in the lists of variables and dimensions (also for the model outputs). Also, ipack is the index for the current file IAGOS package. ilay is the index for the atmospheric layer.
  INTEGER iline, icol, jcol, ivarmod, ivar_masking, idimmod_time, ipack, ipack_correct, ilay
  INTEGER i_which_NOx, ivar_which_NOx
  ! Integers relative to the model grid dimension.
  ! nlev, nlat, nlon is the spatial size of the domain, in matter of gridcells.
  INTEGER nfiles_daily_sum, ntime_sampled
  ! l,k,i are used respectively for the height, latitude, and longitude dimensions.
  ! itime is for the day into the current month. It remains 1 if the time resolution is monthly.
  ! i_current_day is the same, but is read from the IAGOS file name, undependently on the time resolution.
  ! Only for a better clarity in the code.
  INTEGER l,k,i,itime,i_current_day
  CHARACTER*2 c_current_day
  ! Indexes repairing the gridcell closest to the current IAGOS measurement point.
  INTEGER ilon,klat,llev
  ! l_downward is defined as the first level where PV < PV_tropopause, from top to bottom.
  ! l_upward is the same, but from bottom to top.
  ! Using those two indexes allows to account for the double tropopause cases.
  ! Of course, they are the same in case of single tropopause.
  ! l_TP_upper and l_TP_lower are the indexes for the vertical gridcells levels that have to
  ! be filtered out from the upper troposphere. They correspond to the gridcells levels
  ! containing the tropopause pressure (the upper tropopause and the lower one respectively, 
  ! in case of double tropopause).
  INTEGER l_downward, l_upward, l_TP_upper, l_TP_lower
  ! l_DTP is the 2 element array that gathers the indexes of the gridcells to be filtered out
  ! according to this criterion: l_TP_upper and l_TP_lower. Most of the time, they are equal.
  INTEGER l_DTP(2)
  ! Test for water vapour outliers.
  INTEGER n_WV_300
!!! Weight distribution on the closest gridcells for each measurement point. !!!
  ! During this step, we generate a small matrix (nweight x nweight x nweight)
  ! containing the weights of the current measurement point.
  ! The integers with w refer to the position in this matrix. The ones with c define the central coordinate.
  ! The ones with the _end suffix refer to the relative coordinates compared to the central coordinate.
  INTEGER iw,kw,lw,iwc,kwc,lwc,iw_rel,kw_rel,lw_rel
  ! lw_ini and _end define the limits of the pressure interpolation.
  ! Useful in case we deal with the atmosphere boundaries as central vertical grid levels.
  INTEGER lw_ini,lw_end
  ! rotation and i_rot_lon are in case we have to turn the world to find the latitudinal neighbouring cell at the poles.
  LOGICAL rotation
  INTEGER i_rot_lon ! Added longitude cells, making sure the world is not flat.
  ! Values useful to move into the weight matrix from the center.
  INTEGER sign_x, sign_y, sign_P ! Values are either 1 or -1.
  REAL rsign_x, rsign_y ! real equivalents
  ! Useful to avoid writing twice the same code once for an addition, once for a subtraction.
  REAL rsign_add_NOx
  ! Distance (or pressure difference) between the obs coordinates and the gridcells
  ! in the weight matrix, divided by the gridcell size.
  REAL dx_i, dy_k, dz_l
  ! Same as above, this time with the adjacent gridcell pointed out by sign_x/y/P.
  REAL dx_adj, dy_adj, dz_adj
  ! diff_z is the vertical size of the current gridcell, in meters.
  REAL diff_z
  ! Pressures and potential vorticities at the levels l and l-1, l referring here to the
  ! tropopause level.
  ! They are temporary variables, used to derive the tropopause pressure.
  REAL P_l, P_lm1, PV_l, PV_lm1
  ! Tropopauses pressures, derived from an interpolation in the PV field
  ! between the l-1 and l levels.
  REAL P_TP_lower, P_TP_upper
  ! alpha_PV is just planned to be an exponent for the log-pressure interpolation that allows
  ! to derive the tropopause pressure.
  REAL alpha_PV
  REAL val_valid_or_not(2)
!!!!!!! Dynamic arrays !!!!!!!!
  ! Array representing the current IAGOS file.
  CHARACTER*20, ALLOCATABLE :: iagos_data(:,:)
  ! Dimensions read in the current IAGOS file.
  INTEGER, ALLOCATABLE :: dim_val_UTC(:)
  REAL, ALLOCATABLE :: dim_val_lat(:),dim_val_lon(:)
  REAL, ALLOCATABLE :: dim_val_alt(:),dim_val_P(:)
  REAL, ALLOCATABLE :: var_val(:,:) ! First dimension for the variable, second for the data point
  REAL, ALLOCATABLE :: var_val_NOx_1(:),var_val_NOx_2(:) ! The two observed NOx variables, if they exist.
  REAL, ALLOCATABLE :: var_val_NOx_derived(:) ! The last NOx column, to derive from the two others.
  ! Temporary arrays.
  REAL, ALLOCATABLE :: tmp_real(:)
  ! For the daily resolution, we also define daily arrays for each month.
  INTEGER,ALLOCATABLE :: list_nfiles_daily(:)
  ! Listing the variables to be read in the current IAGOS file. Its size varies between two files.
  CHARACTER*30, ALLOCATABLE :: list_var_to_read(:)
  ! Listing the variables to be read in the current IAGOS file or derived. Varying size too.
  CHARACTER*30, ALLOCATABLE :: list_var_this_flight(:)
  ! Listing the corresponding indexes for the columns in the current IAGOS file.
  INTEGER,ALLOCATABLE :: list_col_to_read(:)
  ! Variables that shorten the instructions lines.
  REAL PV_tmp, P_tmp, W_tmp, tmp_value, tmpcount, tmpvar
  ! Checking the model latitude and the model vertical axis.
  REAL, ALLOCATABLE :: check_lon(:), check_lat(:), check_lev(:)
  ! Weight matrix, and weighting coefficients depending on the distance
  ! between the current gridcell and the current measurement point.
  REAL,ALLOCATABLE :: weight(:,:,:), alpha_x(:), alpha_y(:), alpha_P(:)
  ! An array of binaries, for a given flight: 1 if sampled, 0 if not. In the end, it is summed up.
  REAL,ALLOCATABLE :: bin_sampled(:,:,:,:,:)
  ! Daily means derived from observations. In case of a daily resolution, they are used to filter out 
  ! the grid cells where observed ozone is in disagreement with the PV field.
  ! The P array is the mean pressure, there is one per variable.
  REAL,ALLOCATABLE :: sum_6D(:,:,:,:,:,:), min_6D(:,:,:,:,:,:), max_6D(:,:,:,:,:,:), P_6D(:,:,:,:,:,:)
  REAL,ALLOCATABLE :: daily_ozone(:,:,:,:,:), n_diff_updown_daily(:,:,:,:,:,:), N_flights_6D(:,:,:,:,:,:)
  REAL,ALLOCATABLE :: rcount_daily_dry(:,:,:,:,:)
  ! Initialization.
  bin_this_mth=.TRUE.
  ALLOCATE( sum_6D(nvariables,nlayers,nlon,nlat,nlev,ntime))
  ALLOCATE(   P_6D(nvariables,nlayers,nlon,nlat,nlev,ntime))
  ALLOCATE( min_6D(nvariables,nlayers,nlon,nlat,nlev,ntime))
  ALLOCATE( max_6D(nvariables,nlayers,nlon,nlat,nlev,ntime))
  ALLOCATE(n_diff_updown_daily(nvariables,nlayers,nlon,nlat,nlev,ntime))
  ALLOCATE(N_flights_6D(nvariables,nlayers,nlon,nlat,nlev,ntime))
  ALLOCATE(rcount_daily_dry(nlayers,nlon,nlat,nlev,ntime))
  sum_6D(:,:,:,:,:,:)=0.
  P_6D(:,:,:,:,:,:)=0.
  min_6D (:,:,:,:,:,:)=1E+20
  max_6D (:,:,:,:,:,:)=0.
  mean_5D(:,:,:,:,:)=0.
  min_5D (:,:,:,:,:)=1E+20
  max_5D (:,:,:,:,:)=0.
  rcount (:,:,:,:,:)=0.
  day1_5D(:,:,:,:,:)=32
  dayN_5D(:,:,:,:,:)=0
  rcount_daily(:,:,:,:,:,:)=0.
  n_diff_updown(:,:,:,:,:)=0.
  n_diff_updown_daily(:,:,:,:,:,:)=0.
  N_flights_5D(:,:,:,:,:)=0.
  N_flights_6D(:,:,:,:,:,:)=0.
  
  mean_5D_model(:,:,:,:,:)=0.
  min_5D_model (:,:,:,:,:)=1E+20
  max_5D_model (:,:,:,:,:)=0.

  rcount_daily_dry(:,:,:,:,:)=0.

  !     ####################################################
  !     1. Reads a IAGOS file and stocks it into a matrix ##
  !     ####################################################
  CALL INT2CHAR(iyyyymm,tmpdir,cyyyymm)
  cyyyy=cyyyymm(1:4)
  ! If we work on a daily resolution, we need to know
  ! how many files correspond to each day. The program will read it later.
  ! This code also works with the monthly resolution, with ntime=1 and list_nfiles_daily=nfiles.
  ALLOCATE(list_nfiles_daily(ntime))
  IF (TRIM(time_resol) .EQ. 'daily') THEN
     ! First, reads the amount of IAGOS files during each day in the current month.
     fileName=TRIM(tmpdir)//'list_iagos_files/'//TRIM(cyyyy)// &
          '/list_nfiles_daily_'//TRIM(cyyyymm)
     OPEN(14,FILE=fileName,FORM='formatted')
     IF (quicktest.EQ.1) THEN
        READ(14,*) list_nfiles_daily(1)
     ELSE
        READ(14,*) list_nfiles_daily(1:ntime)
     END IF
     CLOSE(14)
  ELSE IF (TRIM(time_resol) .EQ. 'monthly') THEN
     list_nfiles_daily(1)=nfiles
  END IF
  N_min_daily = 0.01/ntime
  ! The integer below is to be the amount of lines stored in a concatenated array,
  ! gathering all the measurements made in the current month from flights that
  ! started during the end of the previous month.
  nlines_after_midnight_concat=0
  DO iFile=1,nfiles
     ! Finds the day corresponding to the current file.
     ! It is necessary also because we need to adapt the code to the simulation calendar:
     ! if it is defined without leap, then we have to forget the IAGOS data taken
     ! during each 29th of February.
     ! The next lines work for both daily and monthly resolutions.
     ! In case of a daily resolution, they also work with data-empty days.
     ! In this way, first we do not need an if loop each time the matter of choosing
     ! between daily and monthly resolutions appears, and second we can keep on coding
     ! quietly, without writing huge conditions on the treatment for days without data.
     ! 1. Initialization.
     itime=1
     nfiles_daily_sum=list_nfiles_daily(itime)
     ! 2. Incrementation.
     DO WHILE(nfiles_daily_sum .LT. iFile) ! This never happens with the monthly resolution.
        itime=itime+1
        nfiles_daily_sum=nfiles_daily_sum+list_nfiles_daily(itime)
     END DO
     WRITE(*,*) 'itime=',itime
     WRITE(*,*) 'nfiles_daily_sum=',nfiles_daily_sum
     ! And it is over. Now, we know the daily coordinate for the current file.
     ! We will only have to be careful whether the UTC time passes through midnight during the flight.
     ! The data file analysis starts right now.
     nlines=0
     data_file=TRIM(data_dir)//TRIM(list_data_files(iFile))
     ! The current file is named ${filename_prefix_IAGOS}${yyyymm}ddhhmmssxx.nc.
     ! We first have to find the position of the dd characters.
     i=LEN_TRIM(filename_prefix_IAGOS) ; i=i+6 ! +6 for yyyymm
     ! Then we extract the day numbers.
     c_current_day=list_data_files(iFile)(i+1:i+2)
     READ(c_current_day,'(i2)') i_current_day
     il_err = NF90_OPEN(TRIM(data_file),NF90_NOWRITE,il_ncfile)
     il_err = NF90_GET_ATT_STRING (il_ncfile, NF90_GLOBAL, 'program', packName)
     IF (il_err .NE. 0) THEN
        char_err_status = NF90_STRERROR(il_err)
        WRITE(*,*) char_err_status
     END IF
     il_err = NF90_INQ_DIMID (il_ncfile, tab_var_names(idimUTC,2), il_ncdim)
     il_err = NF90_INQUIRE_DIMENSION (il_ncfile, il_ncdim, LEN = nlines)
     WRITE(*,*) 'Analysing file:'
     WRITE(*,*) TRIM(data_file)
     ! Finds out which package is the current file from, using the ozone variable.
     ipack_correct = 1
     il_err = NF90_INQ_VARID (il_ncfile, tab_var_names(ivarO3+ncoord,ipack_correct), il_ncvar)
     DO WHILE (il_err .NE. 0 .AND. ipack_correct .LE. npackages+1)
        ipack_correct = ipack_correct + 1
        il_err = NF90_INQ_VARID (il_ncfile, tab_var_names(ivarO3+ncoord,ipack_correct), il_ncvar)
     END DO
     packName_correct = header_packages(ipack_correct)
     WRITE(*,*) il_err, ipack_correct, packName_correct
     
     ! Now finds out the package indicated in the file global attribute.
     WRITE(*,*) 'packName:', packName_correct
     WRITE(unitFile_check,*) 'Analysing file:'
     WRITE(unitFile_check,*) TRIM(data_file)
     WRITE(unitFile_check,*) 'packName:', packName
     ipack=which_one_str(packName, list_packages, size(list_packages))
     ipack=ipack+1 ! The first column is dedicated to the full names denomination.
     WRITE(*,*) 'ipack=',ipack
     IF (ipack .NE. 1) THEN ! Make sure we do not treat files from unwished packages.
        IF (TRIM(header_packages(ipack)).EQ.'IAGOS-CARIBIC') THEN
           il_err = NF90_INQ_VARID (il_ncfile, tab_var_names(ncoord+ivarO3,ipack), il_ncvar)
           IF (il_err.NE.0) THEN
              ipack=ipack+1
           END IF
        END IF
     END IF
     
     ! Last: compares the two packName variables, in order to detect inconsistencies between
     ! the ozone variable name and the package name written in the file global attribute.
     IF (TRIM(packName) .NE. TRIM(packName_correct)) THEN
        WRITE(*,*) 'TRUE packName: ', TRIM(packName_correct), ' instead of ', TRIM(packName)
     END IF

     ! Now, in order to simplify the following, we replace ipack and packName by their robust version.
     packName = packName_correct ; ipack = ipack_correct
     ! There is one bugged file in the database, where the latitude probably has the values from another variable.
     IF (TRIM(list_data_files(iFile)) .eq. "IAGOS_timeseries_2017070501001116_L2_3.1.0.nc4" .OR. &
          TRIM(list_data_files(iFile)) .eq. "IAGOS_timeseries_2017070500573511_L2_3.1.0.nc4") THEN
        ipack=1
     END IF
     IF (ipack .NE. 1) THEN ! Make sure we do not treat files from unwished packages.
!!! Defines the dynamical array list_var_to_read(ncol_to_read).
!!! It is the list of available variables amongst the ones queried by the user.
!!! First: counts the amount of lines (record dimension range).
        il_err = NF90_INQ_DIMID (il_ncfile, tab_var_names(idimUTC,ipack), il_ncdim)
        il_err = NF90_INQUIRE_DIMENSION (il_ncfile, il_ncdim, LEN = nlines)
        WRITE(*,*) 'nlines=',nlines
!!! First: counts the amount of available variables ncol_to_read.
        ncol_to_read=0
        DO ivar=1,(n_obs_variables+ncoord)
           ! il_err = NF90_INQ_VAR
           il_err = NF90_INQ_VARID (il_ncfile, tab_var_names(ivar,ipack), il_ncvar)
           IF (il_err.EQ.0) THEN
              ncol_to_read = ncol_to_read + 1
           END IF
        END DO
!!! nvariables_this_flight gives the amount of observed (during the current flight) and derived variables.
        nvariables_this_flight = ncol_to_read + n_der_variables
!!! Second: defines the corresponding variables that are available in the current file.
!!! Stocks their NetCDF ID in a table.
        ALLOCATE(list_var_to_read(ncol_to_read))
        ALLOCATE(list_col_to_read(ncol_to_read))
        ALLOCATE(list_var_this_flight(nvariables_this_flight))
        icol=0
        DO ivar=1,ncoord
           il_err = NF90_INQ_VARID (il_ncfile, tab_var_names(ivar,ipack), il_ncvar)
           IF (il_err.EQ.0) THEN
              icol=icol+1
              list_var_this_flight(icol)=list_coord(ivar)
              list_var_to_read(icol)=list_coord(ivar)
              list_col_to_read(icol)=il_ncvar
              WRITE(unitFile_check,*) list_coord(ivar), 'found: il_ncvar=',il_ncvar
           ELSE
              WRITE(unitFile_check,*) list_coord(ivar), 'not found'
           END IF
        END DO
        DO ivar=1,n_obs_variables
           il_err = NF90_INQ_VARID (il_ncfile, tab_var_names(ivar+ncoord,ipack), il_ncvar)
           IF (il_err.EQ.0) THEN
              icol=icol+1
              list_var_this_flight(icol)=list_variables(ivar)
              list_var_to_read(icol)=list_variables(ivar)
              list_col_to_read(icol)=il_ncvar
              WRITE(unitFile_check,*) list_variables(ivar), 'found: il_ncvar=',il_ncvar
           ELSE
              WRITE(unitFile_check,*) list_variables(ivar), 'not found'
           END IF
        END DO
        DO ivar=1,n_der_variables
           icol=icol+1
           ! No need to search for the derived variables in the input NetCDF file, it does not exist.
           list_var_this_flight(icol)=list_der_variables(ivar)
           WRITE(unitFile_check,*) list_der_variables(ivar)
        END DO
     END IF
     IF (ipack .NE. 1) THEN ! Make sure we do not treat files from unwished packages.
        ALLOCATE(var_val(nvariables,nlines))
        ALLOCATE(bin_sampled(nvariables,nlayers,nlon,nlat,nlev))
        ivar=0 ! this index will start its incrementation as soon as the coordinates are done
        nlines_quarter=INT(nlines/4)
        bin_allocate_P=.TRUE.
        ALLOCATE(tmp_real(nlines))
        DO icol=1,ncol_to_read
           il_ncvar=list_col_to_read(icol)
           WRITE(unitFile_check,*) list_var_to_read(icol),'var. n.',icol
           IF (icol.eq.idimUTC) THEN
              !     Conversion character string --> real
              ALLOCATE(dim_val_UTC(nlines))
              il_err = NF90_GET_VAR (il_ncfile, il_ncvar, dim_val_UTC)
              WRITE(unitFile_check,*) dim_val_UTC(1), &
                   '-->',dim_val_UTC(nlines)
           ELSE IF (icol.eq.idimlat) THEN
              !     Conversion character string --> real
              ALLOCATE(dim_val_lat(nlines))
              il_err = NF90_GET_VAR (il_ncfile, il_ncvar, dim_val_lat)
              WRITE(unitFile_check,*) dim_val_lat(1), &
                   '-->',dim_val_lat(nlines)
           ELSE IF (icol.eq.idimlon) THEN
              !     Conversion character string --> real
              ALLOCATE(dim_val_lon(nlines))
              il_err = NF90_GET_VAR (il_ncfile, il_ncvar, dim_val_lon)
              WRITE(unitFile_check,*) dim_val_lon(1), &
                   '-->',dim_val_lon(nlines)
           ELSE IF (icol.eq.idimalt) THEN   
              ALLOCATE(dim_val_alt(nlines))
              il_err = NF90_GET_VAR (il_ncfile, il_ncvar, dim_val_alt)
              WRITE(unitFile_check,*) dim_val_alt(1), &
                   '-->',dim_val_alt(nlines)
              WRITE(unitFile_check,*) dim_val_alt(nlines_quarter), &
                   '-->',dim_val_alt(3*nlines_quarter)
           ELSE IF (icol.eq.idimP .and. bin_allocate_P.EQV..TRUE.) THEN
              ALLOCATE(dim_val_P(nlines))
              bin_allocate_P=.false.
              il_err = NF90_GET_VAR (il_ncfile, il_ncvar, dim_val_P)
              WRITE(unitFile_check,*) dim_val_P(1),'-->',dim_val_P(nlines)
              WRITE(unitFile_check,*) dim_val_P(nlines_quarter), &
                   '-->',dim_val_P(3*nlines_quarter)
           ELSE IF (icol.GT.ncoord) THEN
              ! Knowing the current variable and its corresponding column in the data file,
              ! now we need to know the index of the variable in the list_variables vector.
              ivar=1
              DO WHILE (tab_var_names(ivar+ncoord,1).NE.list_var_to_read(icol))
                 ivar=ivar+1
              END DO
              WRITE(unitFile_check,*) 'Checking variable consistency.'
              WRITE(unitFile_check,*) 'tab_var_names(ivar+ncoord,ipack), ivar, list_var_to_read(icol)='
              WRITE(unitFile_check,*) tab_var_names(ivar+ncoord,ipack), ivar, list_var_to_read(icol)
              IF (TRIM(tab_var_names(ivar+ncoord,1)).NE.'Pressure') THEN
                 !     Conversion character string --> real
                 il_err = NF90_GET_VAR (il_ncfile, il_ncvar, tmp_real)
                 var_val(ivar,1:nlines)=tmp_real(1:nlines)
                 DO iline=1,nlines
                    ! This loop is made for the fields filled with '********', found in 201411 notably.
                    ! We have to set these missing values at rwrong_val before converting it in real numbers.
                    IF (MIN(dim_val_lon(iline),dim_val_lat(iline)).LE.-180.1) THEN
                       var_val(ivar,iline)=rwrong_val
                    END IF
                 END DO
                 WRITE(unitFile_check,*) var_val(ivar,1), &
                      '-->',var_val(ivar,nlines)
                 WRITE(unitFile_check,*) var_val(ivar,nlines_quarter), &
                      '-->',var_val(ivar,3*nlines_quarter)
              END IF
           END IF
        END DO
        ! Derives the complementary variables (potential temperature and O3-to-CO ratio)
        WRITE(unitFile_check,*) 'ivarP=',ivarP, tab_var_names(ivarP+ncoord,1), thres_val(ivarP)
        DO iline=1,nlines
           IF (dim_val_P(iline).GT.thres_val(ivarP)) THEN
              var_val(ivarP,iline)=dim_val_P(iline)
           ELSE
              var_val(ivarP,iline)=rwrong_val
           END IF
           IF (dim_val_P(iline).GT.thres_val(ivarP) .and. &
                var_val(ivarT,iline).GT.thres_val(ivarT)) THEN
              var_val(ivarT,iline)=var_val(ivarT,iline)+temp_conversion ! Conversion into Kelvin degrees if not yet done.
              var_val(ivarTheta,iline)= &
                   var_val(ivarT,iline)* &
                   (1.E+05/dim_val_P(iline))**0.288
           ELSE
              var_val(ivarTheta,iline)=rwrong_val
           END IF
           IF (var_val(ivarCO,iline).GT.thres_val(ivarCO) .AND. &
                var_val(ivarO3,iline).GT.thres_val(ivarO3)) THEN
              var_val(ivarO3toCO,iline)= &
                   var_val(ivarO3,iline)/var_val(ivarCO,iline)
           ELSE
              var_val(ivarO3toCO,iline)=rwrong_val
           END IF
        END DO
        ! Derives the missing observed variable among these three: NOx, NO, NO2.
        ! Indeed, the ones that are measured change with the IAGOS package.
        nvar_NOx_this_flight=0
        list_ivar_NOx=(/ivarNOx,ivarNO,ivarNO2/)
        nvar_NOx=SIZE(list_ivar_NOx)
        ALLOCATE(bin_sampled_NOx(nvar_NOx))
        ! First: verifies that two of the three variables are sampled.
        ! If it is the case, identifies which one is missing.
        DO i_which_NOx=1,nvar_NOx
           ! For one of these variables, finds its corresponding column in the permanent IAGOS variables array.
           ivar_which_NOx=list_ivar_NOx(i_which_NOx)
           IF (ANY(list_var_to_read==list_variables(ivar_which_NOx))) THEN
              bin_sampled_NOx(i_which_NOx)=.TRUE.
              nvar_NOx_this_flight=nvar_NOx_this_flight+1
           END IF
        END DO
        IF (nvar_NOx_this_flight.EQ.2) THEN
           WRITE(*,*) 'Deriving the last NOx species:'
           ALLOCATE(var_val_NOx_1(nlines))
           ALLOCATE(var_val_NOx_2(nlines))
           ALLOCATE(var_val_NOx_derived(nlines))
           ! Which one is not sampled?
           i_which_NOx=1
           DO WHILE (bin_sampled_NOx(i_which_NOx).EQV. .TRUE.)
              i_which_NOx=i_which_NOx+1
           END DO
           IF (i_which_NOx .EQ. 1) THEN
              ! Then NOx has to be deduced from the addition between NO and NO2.
              rsign_add_NOx=1.
              var_val_NOx_1=var_val(ivarNO,:)
              var_val_NOx_2=var_val(ivarNO2,:)
           ELSE
              ! Then NO/NO2 (to determine) has to be deduced from the subtraction between NOx and NO2/NO.
              var_val_NOx_1=var_val(ivarNOx,:)
              rsign_add_NOx=-1.
              IF (i_which_NOx .EQ. 2) THEN
                 var_val_NOx_2=var_val(ivarNO2,:) ! NO to deduce.
              ELSE
                 var_val_NOx_2=var_val(ivarNO,:)  ! NO2 to deduce.
              END IF
           END IF
           WRITE(*,*) list_variables(list_ivar_NOx(:))
           WRITE(*,*) 'i_which_NOx, its variable=', &
                i_which_NOx,list_variables(list_ivar_NOx(i_which_NOx))
           DO iline=1,nlines
              IF (MIN(var_val_NOx_1(iline),var_val_NOx_2(iline)) .GT. thres_val(ivarNOx)) THEN
                 var_val_NOx_derived(iline)=var_val_NOx_1(iline)+ &
                      rsign_add_NOx*var_val_NOx_2(iline)
              ELSE
                 var_val_NOx_derived(iline)=rwrong_val
              END IF
           END DO
           ! Last: replaces the right column in the data array by the temporary vector var_val_NOx_derived.
           ! Let's find which one it is, knowing its index in the list_ivar_NOx array is i_which_NOx.
           ivar_which_NOx=list_ivar_NOx(i_which_NOx)
           ! Replacement.
           var_val(ivar_which_NOx,:)=var_val_NOx_derived(:)
           DEALLOCATE(var_val_NOx_1)
           DEALLOCATE(var_val_NOx_2)
           DEALLOCATE(var_val_NOx_derived)
        END IF
        DEALLOCATE(bin_sampled_NOx)
        WRITE(unitFile_check,*) list_variables(ivarTheta)
        WRITE(unitFile_check,*) var_val(ivarTheta,1), &
             '-->',var_val(ivarTheta,nlines)
        WRITE(unitFile_check,*) var_val(ivarTheta,nlines_quarter), &
             '-->',var_val(ivarTheta,3*nlines_quarter)
        WRITE(unitFile_check,*) list_variables(ivarO3toCO)
        WRITE(unitFile_check,*) var_val(ivarO3toCO,1), &
             '-->',var_val(ivarO3toCO,nlines)
        WRITE(unitFile_check,*) var_val(ivarO3toCO,nlines_quarter), &
             '-->',var_val(ivarO3toCO,3*nlines_quarter)
        DEALLOCATE(tmp_real)
        DEALLOCATE(list_col_to_read)
        il_err = NF90_CLOSE(il_ncfile)
        !     ###############################################################
        !     2. Repairing the adjacent grid levels, and assigning weights ##
        !     ###############################################################
        !     Temporarily: choses the inferior bound of the belonged grid cell
        !     First: removes the data taken at the end of the month, after midnight.
        tmpvar=0
        tmpcount=0
        IF (itime .EQ. ntime .AND. ntime .NE. 1) THEN
           IF (dim_val_UTC(nlines) .GE. n_seconds_a_day) THEN
              nlines_until_midnight=1
              DO WHILE(dim_val_UTC(nlines_until_midnight).LT.n_seconds_a_day)
                 nlines_until_midnight = nlines_until_midnight + 1
              END DO
              nlines_this_flight = nlines
              nlines = nlines_until_midnight
              nlines_after_midnight = nlines_this_flight - nlines_until_midnight
              nlines_after_midnight_concat = nlines_after_midnight_concat + nlines_after_midnight
              ! Stocks next month data into a temporary array.
              ! Finally not. Just removes them, for the moment.
              ! If later we try to code it, we shoud think about concatenating the past-midnight data
              ! from all the month-end flights into this temporary array.
              ! This would avoid us to allocate additionally next month fields.
              ! However, be VERY careful with the order between the variables.
              ! Remember it changes from one file to another.
           END IF
        END IF
        bin_day_passed=.FALSE.
        bin_sampled(:,:,:,:,:) = 0.
        IF (itime .LE. ntime) THEN ! In this way, we forget every 29th of February if it does not exist in the simulation.
           DO iline=1,nlines
              ! First: using UTC time, we account for the nighttime flights.
              ! Indeed, the current method assigns a date to the flight by using the first measurement only.
              ! We have to correct it by assigning the measurements made after midnight to the morrow.
              IF (dim_val_UTC(iline) .GE. n_seconds_a_day) THEN
                 IF ((itime .ne. ntime).AND.(bin_day_passed .EQV..FALSE.)) THEN
                    itime=itime+1 ! itime itself can be changed, since it is reinitialized for each flight.
                 ELSE
                    itime=ntime ! Avoids segmentation faults. It is an approximation that simplifies the code.
                 END IF
                 bin_day_passed=.TRUE.
              END IF
              ! Longitude and latitude: The while loops (lon and lat) look for
              ! the gridcell that contains the measurement point.
              ! But first: boundary condition on longitudes, in case the first gridcell
              ! has a part on each side of the 180°E meridian, which can be problematic
              ! if the measurement point is on the extreme East.
              IF (dim_val_lon(iline) > rlon(nlon) + dlon/2 ) THEN
                 ilon = 1
                 rlon_eq = rlon(ilon) + 360
              ELSE
                 ilon = 1
                 DO WHILE ((rlon(ilon) + dlon/2) < dim_val_lon(iline) .AND. ilon < nlon)
                    ilon = ilon + 1
                 END DO
                 rlon_eq = rlon(ilon)
              END IF
              klat = 1
              DO WHILE ((rlat(klat) + dlat/2) < dim_val_lat(iline) .AND. klat < nlat)
                 klat = klat + 1
              END DO
              r_sign=SIGN(1.,dim_val_lat(iline))
              i_sign=INT(r_sign)
              ALLOCATE(alpha_x(nweight)) ! The alpha coefficients are an intermediate step to implement the interpolation matrix "weight"
              ALLOCATE(alpha_y(nweight))
              ALLOCATE(alpha_P(nweight))
              ALLOCATE(weight(nweight,nweight,nweight))
              ! Vertical levels
              l=1
              iwc=(nweight+1)/2 ! the letter c stands for "central"
              kwc=(nweight+1)/2
              lwc=(nweight+1)/2
              DO iw=1,nweight
                 DO kw=1,nweight
                    DO lw=1,nweight
                       weight(iw,kw,lw)=0.
                    END DO
                 END DO
              END DO
              
              alpha_x=0 !(/0,0,0/) more precisely, for a linear interpolation with distance
              alpha_y=0 !(/0,0,0/)
              alpha_P=0 !(/0,0,0/)
              IF (dim_val_P(iline) .GT. 0) THEN
                 DO WHILE ((P_mid(ilon,klat,l,itime).LT.dim_val_P(iline)).AND.(l.LT.nlev))
                    l=l+1
                 END DO
                 ! The last thing we have to keep in mind during this retrieval:
                 ! the model levels can be given in monthly means. In this case,
                 ! it happens quite often that the ground model level is not
                 ! representative of the instantaneous lower pressure levels.
                 ! In order to avoid any bug, I chose to associate the orphan points
                 ! to the last level (ex.: normal ground pressure, during a month
                 ! dominated by a cyclone)
                 llev=l
                 ! Defining the relative positions of the observation between the cells ilon and ilon+/-1.
                 ! The distances (in degrees, it's ok) necessarily refer to the centers of the cells.
                 ! Thus, the adjacent cell can be the next one or the previous one depending on
                 ! the position of the observation compared to the center of the cell ilon.
                 ! ###############
                 dx_i=ABS(dim_val_lon(iline)-rlon_eq)/dlon
                 rsign_x=SIGN(1.,dim_val_lon(iline)-rlon_eq)
                 sign_x=INT(rsign_x)
                 dx_adj=1.-dx_i
                 alpha_x(iwc)=1.-dx_i
                 alpha_x(iwc+sign_x)=1.-dx_adj ! I know, it equals dx_i, but I want to keep the calculations clear. At least, it shows explicitely the meaning of the alpha coefficients.
                 ! ###############
                 ! Idem with the latitudes
                 dy_k=ABS(dim_val_lat(iline)-rlat(klat))/dlat
                 rsign_y=SIGN(1.,dim_val_lat(iline)-rlat(klat))
                 sign_y=INT(rsign_y)
                 dy_adj=1.-dy_k
                 alpha_y(kwc)=1.-dy_k
                 alpha_y(kwc+sign_y)=1.-dy_adj ! I know, it equals dy_k, but I want to keep the calculations clear. At least, it shows explicitely the meaning of the alpha coefficients.
                 !###############
                 ! Now, idem with the pressures. One has to keep in mind that :
                 ! - the vertical grid is not necessarily regular. I thus have to compute the height difference between the pressure levels.
                 ! - the distance (in Pa) refers to the full pressure levels, not the intermediate levels.
                 ! This is where the vectors P_mid appear to be useful.
                 ! 1/ First, we treat the atmosphere boundaries: the ground and the roof pressures
                 IF (llev .EQ. nlev) THEN
                    lw_ini=1 ; lw_end=nweight-lwc+1 ! This decides which part of the moving grid will be superimposed with the global one.
                    ! Indeed, if the moving grid is centered on one boundary cell, then a part of the moving grid
                    ! will not have any support. This part has to be shaded.
                    ! New sub-condition: if the measurement is located below the full vertical grid level,
                    ! the surface level swallows the weight from the removed part.
                    IF (dim_val_P(iline)>=P_mid(ilon,klat,llev,itime)) THEN
                       alpha_P(lwc)=1
                    ELSE ! Then: normal treatment. See the "else" condition below (commented "Normal condition").
                       diff_z=ABS(LOG(P_mid(ilon,klat,llev,itime)/P_mid(ilon,klat,llev-1,itime)))
                       dz_l=ABS(LOG(dim_val_P(iline)/P_mid(ilon,klat,llev,itime)))/diff_z
                       dz_adj=ABS(LOG(dim_val_P(iline)/P_mid(ilon,klat,llev-1,itime)))/diff_z
                       alpha_P(lwc)=1.-dz_l
                       alpha_P(lwc-1)=1.-dz_adj ! I know, it equals dz_l, but I want to keep the calculations clear. At least, it shows explicitely the meaning of the alpha coefficients.
                    END IF
                 ELSE IF (llev.EQ.1 .AND. dim_val_P(iline)<P_mid(ilon,klat,llev,itime)) THEN
                    alpha_P(lwc)=1
                    lw_ini=lwc ; lw_end=nweight
                 ELSE ! Normal condition.
                    lw_ini=1 ; lw_end=nweight
                    ! 2/ Now, we treat the main grid domain.
                    ! diff_P=ABS(P_mid(ilon,klat,llev)-P_mid(ilon,klat,llev-1,itime))
                    diff_z=ABS(LOG(P_mid(ilon,klat,llev,itime)/P_mid(ilon,klat,llev-1,itime)))
                    dz_l=ABS(LOG(dim_val_P(iline)/P_mid(ilon,klat,llev,itime)))/diff_z
                    dz_adj=ABS(LOG(dim_val_P(iline)/P_mid(ilon,klat,llev-1,itime)))/diff_z
                    alpha_P(lwc)=1.-dz_l
                    alpha_P(lwc-1)=1.-dz_adj ! I know, it equals dz_l, but I want to keep the calculations clear. At least, it shows explicitely the meaning of the alpha coefficients.
                 END IF
                 DO iw=1,nweight
                    DO kw=1,nweight
                       DO lw=1,nweight
                          ! Last: in order to avoid some numeric bugs due to too many decimals.
                          alpha_x(iw)=aINT(alpha_x(iw)*1000)/1000
                          alpha_y(kw)=aINT(alpha_y(kw)*1000)/1000
                          alpha_P(lw)=aINT(alpha_P(lw)*1000)/1000
                          weight(iw,kw,lw)=alpha_x(iw)*alpha_y(kw)*alpha_P(lw)
                          ! Each element of the matrix weight is a distribution factor. Respectively on each axis,
                          ! the central element (iwc,kwc,lwc) corresponds to the grid point (ilon,klat,llev).
                          ! For a given variable "ivar", the sum of all the elements in this matrix has to be 1.
                          IF (weight(iw,kw,lw) .LT. 0 .and. &
                               dim_val_lon(iline) .GT. -180.1 .and. &
                               dim_val_lat(iline) .GT. -90.1) THEN ! This second condition excludes the cases already taken into account in the following code (next paragraph).
                             WRITE(*,*) 'rsign_y=', rsign_y
                             WRITE(*,*) 'sign_y=', sign_y
                             WRITE(*,*) 'negative weight at the location:', &
                                  dim_val_lon(iline), dim_val_lat(iline), dim_val_P(iline)
                             WRITE(*,*) 'P_mid(llev-1), P_edge(llev), P_mid(llev)= ', &
                                  P_mid(ilon,klat,llev-1,itime), P_edge(ilon,klat,l,itime), &
                                  P_mid(ilon,klat,llev,itime)
                             WRITE(*,*) 'vertical weight vector: ', &
                                  weight(iw,kw,:)
                             WRITE(*,*) 'alpha_x(:)=',alpha_x(:)
                             WRITE(*,*) 'alpha_y(:)=',alpha_y(:)
                             WRITE(*,*) 'alpha_P(:)=',alpha_P(:)
                             STOP
                          END IF
                       END DO
                    END DO
                 END DO
              END IF
              ! The superimposition of the global grid and the moving weight grid is not done yet.
              ! iw_rel, kw_rel and lw_rel are the indexes iw, kw and lw expressed relatively to the central point (iwc,kwc,lwc)
              ! Those ones are used to locate the observation relatively to the current grid point.
              IF (MIN(dim_val_lon(iline),dim_val_lat(iline)).GT.-180.1) THEN
                 ! This condition on longitudes because of the apparison of missing values in this coordinate.
                 ! Last condition (on latitudes) because now, we never know what can happen.
                 ! Knowing the current variable and its corresponding column in the data file,
                 ! now we need to know the index of the variable in the list_variables vector.
                 ! The previous version's command ivar=ivar+1 does not work anymore,
                 ! because the order between the columns changes from one file to another,
                 ! depending on which variables are measured.
                 DO icol=(ncoord+1),nvariables_this_flight ! Variables only. No need for coordinates.
                    ivar=1
                    DO WHILE (tab_var_names(ivar+ncoord,1) .NE. list_var_this_flight(icol))
                       ivar=ivar+1
                    END DO
                    ! WRITE(*,*) icol,list_var_this_flight(icol), tab_var_names(ivar+ncoord,1),ivar
                    IF (var_val(ivar,iline) .GT. thres_val(ivar)) THEN
                       DO kw=1,nweight
                          ! Boundary conditions, at 90°S or 90°N for the global grid.
                          ! I know, IAGOS does not perform measurements at these latitudes,
                          ! but maybe one day...
                          IF ((klat.EQ.nlat .AND. sign_y.EQ.1) .OR. (klat.EQ.1 .AND. sign_y.EQ.-1)) THEN
                             kw_rel=0
                             rotation=.TRUE.
                          ELSE
                             kw_rel=kw-kwc
                             rotation=.FALSE.
                          END IF
                          DO iw=1,nweight
                             ! Boundary conditions: make sure Earth planet is not flat.
                             IF (ilon.EQ.nlon .AND. sign_x.EQ.1) THEN
                                iw_rel=-nlon+1 ! This will ensure that the cell after "nlon" is the number 1
                             ELSE IF (ilon.EQ.1 .AND. sign_x.EQ.-1) THEN
                                iw_rel=nlon-1 ! This will ensure that the cell before the number 1 is "nlon"
                             ELSE
                                iw_rel=iw-iwc ! Normal conditions
                             END IF
                             
                             IF (rotation.EQV..TRUE.) THEN
                                IF (ilon+iw_rel<=nlon/2) THEN
                                   i_rot_lon=nlon/2 ! This is the condition on the longitudes too: it targets the opposite meridian.
                                ELSE
                                   i_rot_lon=-nlon/2
                                END IF
                             ELSE
                                i_rot_lon=0
                             END IF
                             IF (logic_PV .EQV. .TRUE.) THEN ! Only if we ask for the discrimination.
                                ! Using the PV field to derive the tropopause(s) pressure.
                                IF (ABS(rlat(klat+kw_rel)) .GT. 25) THEN ! Preliminatory condition: extratropics.
                                   ! Finds the first vertical grid level below the dynamical tropopause, from the top of atmosphere.
                                   ! I start it at 2 because the PV fields that I computed from INCA are 0 at the top level.
                                   l_downward=2
                                   DO WHILE (ABS(PV(ilon+iw_rel+i_rot_lon,klat+kw_rel,l_downward,itime)).GT.PV_tropopause &
                                        .AND. l_downward < nlev-1)
                                      l_downward=l_downward+1
                                   END DO
                                   ! Shortening the line commands...
                                   PV_l  =ABS(PV(ilon+iw_rel+i_rot_lon,klat+kw_rel,l_downward  ,itime))
                                   PV_lm1=ABS(PV(ilon+iw_rel+i_rot_lon,klat+kw_rel,l_downward-1,itime))
                                   P_l   =P_mid(ilon+iw_rel+i_rot_lon,klat+kw_rel,l_downward  ,itime)
                                   P_lm1 =P_mid(ilon+iw_rel+i_rot_lon,klat+kw_rel,l_downward-1,itime)
                                   ! Defining the exponent for the log-pressure interpolation.
                                   alpha_PV  =(PV_lm1 - PV_tropopause) / (PV_lm1 - PV_l)
                                   ! Deriving the upper tropopause pressure.
                                   P_TP_upper=P_lm1 * (P_l / P_lm1)**alpha_PV
                                   ! Same step, this time going upward. Useful in case of a double tropopause event.
                                   ! I start it at nlev-2 because they are hardly realistic at the surface level.
                                   ! The -5 is just in case errors are also present in the other neighouring levels.
                                   l_upward=nlev-2
                                   DO WHILE (ABS(PV(ilon+iw_rel+i_rot_lon,klat+kw_rel,l_upward,itime)).LT.PV_tropopause .AND. l_upward > 1)
                                      l_upward=l_upward-1
                                   END DO
                                   l_upward=l_upward+1 ! The level just below, consistently with l_downward definition.
                                   ! Shortening the line commands...
                                   PV_l  =ABS(PV(ilon+iw_rel+i_rot_lon,klat+kw_rel,l_upward  ,itime))
                                   PV_lm1=ABS(PV(ilon+iw_rel+i_rot_lon,klat+kw_rel,l_upward-1,itime))
                                   P_l   =P_mid(ilon+iw_rel+i_rot_lon,klat+kw_rel,l_upward  ,itime)
                                   P_lm1 =P_mid(ilon+iw_rel+i_rot_lon,klat+kw_rel,l_upward-1,itime)
                                   ! Defining the exponent for the log-pressure interpolation.
                                   alpha_PV  =(PV_lm1 - PV_tropopause) / (PV_lm1 - PV_l)
                                   ! Deriving the lower tropopause pressure.
                                   P_TP_lower=P_lm1 * (P_l / P_lm1)**alpha_PV
                                   ! Last: we define the gridcell to be filtered out for the upper troposphere.
                                   ! It has to be the one(s) that contain(s) (one of) the tropopause pressure(s).
                                   ! And yes, this sentence full of s between parentheses is awful. Sorry for that. Anyway.
                                   ! In this way, we filter out what happens just below the tropopause,
                                   ! in a height spreading from a half grid cell to a full grid cell.
                                   ! Sufficient but not too much.
                                   ! Upper tropopause first.
                                   IF (P_edge(ilon+iw_rel+i_rot_lon,klat+kw_rel,l_downward-1,itime) .GT. P_TP_upper) THEN
                                      l_TP_upper=l_downward-1
                                   ELSE
                                      l_TP_upper=l_downward
                                   END IF
                                   ! Lower tropopause then.
                                   IF (P_edge(ilon+iw_rel+i_rot_lon,klat+kw_rel,l_upward-1,itime) .GT. P_TP_lower) THEN
                                      l_TP_lower=l_upward-1
                                   ELSE
                                      l_TP_lower=l_upward
                                   END IF
                                ELSE
                                   l_TP_lower=1
                                   l_TP_upper=1
                                END IF
                                l_DTP=(/ l_TP_upper, l_TP_lower/)
                             ELSE
                                l_DTP=(/ 0, 0/)
                             END IF ! IF (logic_PV .EQV. .TRUE.) THEN
                             DO lw=lw_ini,lw_end
                                lw_rel=lw-lwc
                                ! Now, we can start to distribute the observed values on the grid.
                                ! By the way, we account for the PV to distribute the data over
                                ! the atmospheric layers also (UT, TPL, LS,...).
                                ! This is the last moment we keep the information on the day.
                                ! It is time to use it.
                                PV_tmp=ABS(PV(ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime)) ! Current gridcell's PV.
                                P_tmp=P_mid(ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime) ! Current gridcell's full P.
                                DO ilay=1,nlayers
                                   W_tmp=0
                                   ! Filtering with PV and pressure.
                                   IF ((PV_tmp.LE.PV_up_layers(ilay)) .AND. (PV_tmp.GT.PV_down_layers(ilay)) .AND. &
                                        (P_tmp.LT.P_down_layers(ilay)) .AND. (P_tmp.GE.P_up_layers(ilay))) THEN
                                      ! New filtering: water vapour outliers. Here, we test their impact on the climatologies.
                                      ! Filtering with ozone volume mixing ratio.
                                      ! IF ((var_val(ivarO3,iline).lt.ozone_max_layers(ilay)).and. &
                                      !      (var_val(ivarO3,iline).gt.ozone_min_layers(ilay))) THEN
                                      ! For the upper troposphere, filtering out the gridcell
                                      ! which lower boundary is the first one below the/each tropopause.
                                      IF ((TRIM(list_layers(ilay)).NE.'UT') .OR. (all( l_DTP.NE.(llev+lw_rel) ))) THEN
                                         min_6D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime)= &
                                              MIN(min_6D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime), &
                                              var_val(ivar,iline))
                                         max_6D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime)= &
                                              MAX(max_6D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime), &
                                              var_val(ivar,iline))
                                         IF (ntime .EQ. 1) THEN ! Special case for the monthly resolution.
                                            day1_5D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel)= &
                                                 MIN(i_current_day,day1_5D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel))
                                            dayN_5D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel)= &
                                                 MAX(i_current_day,dayN_5D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel))
                                         END IF
                                         ! Last thing: assigning W_tmp the weight coefficient if the conditions are ok.
                                         W_tmp=weight(iw,kw,lw)
                                      END IF
                                   END IF
                                   rcount_daily(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime)= &
                                        rcount_daily(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime)+ &
                                        W_tmp
                                   sign_P=SIGN(1.,dim_val_P(iline)-P_tmp)
                                   n_diff_updown_daily(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime)= &
                                        n_diff_updown_daily(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime)+ &
                                        W_tmp*sign_P
                                   sum_6D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime)= &
                                        sum_6D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime)+ &
                                        var_val(ivar,iline)*W_tmp
                                   P_6D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime)= &
                                        P_6D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime)+ &
                                        var_val(ivarP,iline)*W_tmp
                                   bin_sampled(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel) = &
                                        bin_sampled(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel)+ &
                                        W_tmp
                                   IF (ivar.EQ.ivarO3 .AND. &
                                        (sum_6D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime).LT.0 .OR. &
                                        rcount_daily(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime).LT.0)) THEN
                                      WRITE(*,*) 'sum_6D=', sum_6D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime)
                                      WRITE(*,*) 'rcount_daily=', rcount_daily(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel,itime)
                                      WRITE(*,*) 'Negative ozone or count value :', &
                                           mean_5D(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel), &
                                           rcount(ivar,ilay,ilon+iw_rel+i_rot_lon,klat+kw_rel,llev+lw_rel)
                                      WRITE(*,*) 'Negative value reached for var_val and weight =', &
                                           var_val(ivar,iline), weight(iw,kw,lw), W_tmp
                                      WRITE(*,*) 'Coordinates lon, lat, P =', &
                                           dim_val_lon(iline), dim_val_lat(iline), dim_val_P(iline)
                                      WRITE(*,*) 'sign_P = ', sign_P
                                      WRITE(*,*) 'P_mid(llev+sign_P), P_edge(llev), P_mid(llev) = ', &
                                           P_mid(ilon,klat,llev+sign_P,itime), P_edge(ilon,klat,l,itime),  &
                                           P_mid(ilon,klat,llev,itime)
                                      WRITE(*,*) 'alpha P-1, alpha P and alpha P+1 = ', &
                                           alpha_P
                                      STOP
                                   END IF
                                END DO
                             END DO
                          END DO
                       END DO
                    END IF
                 END DO
              END IF
              DEALLOCATE(weight,alpha_x,alpha_y,alpha_P)
           END DO ! Ends the loop on lines
           WHERE (bin_sampled(:,:,:,:,:) .GT. 0.01)
              N_flights_6D(:,:,:,:,:,itime) = N_flights_6D(:,:,:,:,:,itime) + 1.
           END WHERE
        ELSE IF (itime .GT. ntime) THEN
           WRITE(*,*) '29th of February does not exist in the simulation, so we skip the corresponding IAGOS files.'
        END IF
        DEALLOCATE(var_val)
        DEALLOCATE(dim_val_UTC)
        DEALLOCATE(dim_val_lat)
        DEALLOCATE(dim_val_lon)
        DEALLOCATE(dim_val_alt)
        DEALLOCATE(dim_val_P)
        DEALLOCATE(list_var_to_read)
        DEALLOCATE(list_var_this_flight)
        DEALLOCATE(bin_sampled)
     ELSE
        WRITE(*,*) 'Skipping this non-asked '//packName(1:char_len)//' file.'
     END IF
  END DO ! END IFile=1,nfiles
  
  ! Now the sums are done while deriving the weight coefficients,
  ! we just have to get the mean values.
  ! Filtering step, day-by-day if the PV resolution is daily. 
  ! We derive ozone daily means from the IAGOS measurements, and we use it as a criterion
  ! to validate or to filter out the days following the consistency with the PV field.
  ! This filter assumes the daily means from observations
  ! to be representative of the whole day.
  WRITE(unitFile_check,*) '########## Computes IAGOS ozone daily mean values ########'
  WRITE(*,*)              '########## Computes IAGOS ozone daily mean values ########'
  ALLOCATE(daily_ozone(nlayers,nlon,nlat,nlev,ntime))
  DO ilay=1,nlayers
     WHERE (rcount_daily(ivarO3,ilay,:,:,:,:).GT.N_min_daily)
        daily_ozone(ilay,:,:,:,:)= &
             sum_6D(ivarO3,ilay,:,:,:,:)/rcount_daily(ivarO3,ilay,:,:,:,:)*IAGOS_conv_factor(ivarO3)
     END WHERE
     ! Filtering out the inconsistent grid cells.
     DO ivar=1,nvariables
        WHERE (rcount_daily(ivar,ilay,:,:,:,:) .LT. N_min_daily .OR. &
             (daily_ozone(ilay,:,:,:,:) .LT. ozone_min_layers(ilay) .AND. &
             daily_ozone(ilay,:,:,:,:) .GT. thres_val(ivarO3) ) .OR. &
             daily_ozone(ilay,:,:,:,:) .GT. ozone_max_layers(ilay))
           rcount_daily(ivar,ilay,:,:,:,:)=0.
        END WHERE
     END DO
  END DO
  WRITE(unitFile_check,*) '########## Computes IAGOS monthly mean values ########'
  WRITE(*,*)              '########## Computes IAGOS monthly mean values ########'
  DO ivar=1,nvariables
     ! Summing up the daily sums.
     DO ilay=1,nlayers
        DO i=1,nlon
           DO k=1,nlat
              DO l=1,nlev
                 mean_5D(ivar,ilay,i,k,l) = sum(sum_6D(ivar,ilay,i,k,l,:), &
                      rcount_daily(ivar,ilay,i,k,l,:) .GT. N_min_daily)
                 rcount(ivar,ilay,i,k,l) = sum(rcount_daily(ivar,ilay,i,k,l,:), &
                      rcount_daily(ivar,ilay,i,k,l,:) .GT. N_min_daily)
                 max_5D(ivar,ilay,i,k,l) = maxval(max_6D(ivar,ilay,i,k,l,:), &
                      rcount_daily(ivar,ilay,i,k,l,:) .GT. N_min_daily)
                 min_5D(ivar,ilay,i,k,l) = minval(min_6D(ivar,ilay,i,k,l,:), &
                      rcount_daily(ivar,ilay,i,k,l,:) .GT. N_min_daily)
                 n_diff_updown(ivar,ilay,i,k,l) = sum(n_diff_updown_daily(ivar,ilay,i,k,l,:), &
                      rcount_daily(ivar,ilay,i,k,l,:) .GT. N_min_daily)
                 N_flights_5D(ivar,ilay,i,k,l) = sum(N_flights_6D(ivar,ilay,i,k,l,:), &
                      rcount_daily(ivar,ilay,i,k,l,:) .GT. N_min_daily)
                 P_5D(ivar,ilay,i,k,l) = sum(P_6D(ivar,ilay,i,k,l,:), &
                      rcount_daily(ivar,ilay,i,k,l,:) .GT. N_min_daily)
              END DO
           END DO
        END DO
     END DO
     IF (ntime .NE. 1) THEN ! This would not work with the monthly resolution.
        DO ilay=1,nlayers
           DO i=1,nlon
              DO k=1,nlat
                 DO l=1,nlev
                    DO itime=1,ntime
                       IF (rcount_daily(ivar,ilay,i,k,l,itime) .GT. N_min_daily) THEN
                          dayN_5D(ivar,ilay,i,k,l) = itime
                          IF (day1_5D(ivar,ilay,i,k,l) .EQ. 32) THEN ! We increment this index only once.
                             day1_5D(ivar,ilay,i,k,l) = itime
                          END IF
                       END IF
                    END DO
                 END DO
              END DO
           END DO
        END DO
     END IF
     ! Deriving the last metrics: mean values, min, max...
     WHERE (rcount(ivar,:,:,:,:) .GT. 0.01)
        mean_5D(ivar,:,:,:,:)=mean_5D(ivar,:,:,:,:)/rcount(ivar,:,:,:,:)*IAGOS_conv_factor(ivar)
        min_5D (ivar,:,:,:,:)=min_5D (ivar,:,:,:,:)                     *IAGOS_conv_factor(ivar)
        max_5D (ivar,:,:,:,:)=max_5D (ivar,:,:,:,:)                     *IAGOS_conv_factor(ivar)
        P_5D   (ivar,:,:,:,:)=P_5D   (ivar,:,:,:,:)/rcount(ivar,:,:,:,:)*IAGOS_conv_factor(ivarP)
        ! Normalizing n_diff_updown. We want it as a real between -1 and 1.
        n_diff_updown(ivar,:,:,:,:)=n_diff_updown(ivar,:,:,:,:)/rcount(ivar,:,:,:,:)
     ELSEWHERE (rcount(ivar,:,:,:,:) .LE. 0.01)
        max_5D(ivar,:,:,:,:) =rwrong_val
        mean_5D(ivar,:,:,:,:)=rwrong_val
        min_5D(ivar,:,:,:,:) =rwrong_val
        rcount(ivar,:,:,:,:) =0.
        day1_5D(ivar,:,:,:,:)=rwrong_val
        dayN_5D(ivar,:,:,:,:)=rwrong_val
        n_diff_updown(ivar,:,:,:,:)=0.
        N_flights_5D(ivar,:,:,:,:)=0.
        P_5D(ivar,:,:,:,:)=rwrong_val
     END WHERE
  END DO ! end ivar
  
  ilon_ok(:)=0
  klat_ok(:)=0
  llev_ok(:)=0
  IF (ANY(list_variables=='Ozone')) THEN
     ivar=ivarO3
  ELSE
     ivar=1 ! Whatever the 1st variable is.
  END IF
  DO ilay=1,nlayers
     DO i=1,nlon
        DO k=1,nlat
           DO l=1,nlev
              IF (rcount(ivar,ilay,i,k,l) .GT. 0.5) THEN
                 ilon_ok(ilay)=i
                 klat_ok(ilay)=k
                 llev_ok(ilay)=l
              END IF
           END DO
        END DO
     END DO
  END DO
  WRITE(unitFile_check,*) '########## Computes model monthly mean values (using temporal IAGOS mask) ########'
  WRITE(*,*) '########## Computes model monthly mean values (using temporal IAGOS mask) ########'
!!!!!! Reading the model chemical fields.
  il_err = NF90_OPEN (TRIM(filename_model), NF90_NOWRITE, il_ncfile)
  il_err = NF90_INQ_DIMID (il_ncfile, TRIM(ctime_model), il_ncdim_time)
  WRITE(*,*) ctime_model, il_ncdim_time
  WRITE(*,*) 'Reading the model fields in ', filename_model
  WRITE(unitFile_check,*) 'Reading the model fields in ', filename_model
  IF (il_err .NE. 0) THEN
     WRITE(*,*) 'Model file not found.'
     WRITE(unitFile_check,*) 'Model file not found. NF90_INQ_VARID returned',il_err
  END IF
  ! Checking the coordinates in the model file.
  ALLOCATE(check_lon(nlon),check_lat(nlat))
  il_err = NF90_INQ_VARID (il_ncfile, TRIM(clon_model), il_ncdim_lon)
  il_err = NF90_GET_VAR   (il_ncfile,     il_ncdim_lon, check_lon(:))
  il_err = NF90_INQ_VARID (il_ncfile, TRIM(clat_model), il_ncdim_lat)
  il_err = NF90_GET_VAR   (il_ncfile,     il_ncdim_lat, check_lat(:))
  IF (check_lat(2)-check_lat(1) .LE. 0. .OR. & ! northward latitude axis
       P_mid(ilon_check,klat_check,2,1)-P_mid(ilon_check,klat_check,1,1) .LE. 0. .OR. & ! downward vertical axis
     ABS(check_lon(nlon))-ABS(check_lon(1)) .GE. 300) THEN ! longitude axis centered around 0 degree
     WRITE(*,*) "----------------------------------------------------"
     WRITE(*,*) "ERROR: reversed axis in the model file: ", TRIM(cyyyymm)
     WRITE(*,*) "----------------------------------------------------"
     WRITE(*,*) P_mid(ilon_check,klat_check,1,1), P_mid(ilon_check,klat_check,2,1)
     WRITE(*,*) check_lon(1), check_lon(nlon)
     WRITE(*,*) check_lat(1), check_lat(2)
     STOP
  END IF
  
  il_err = NF90_INQ_VARID (il_ncfile, TRIM(list_coord(idimlat)), il_ncdim_lev)
  DO ivarmod=1,nvar_model
     il_err = NF90_INQ_VARID (il_ncfile, TRIM(model_variables(ivarmod)), il_ncvar)
     WRITE(*,*) ivarmod, model_variables(ivarmod), il_ncvar, "err=", il_err
     IF (il_err .EQ. 0) THEN
        il_err = NF90_GET_VAR (il_ncfile, il_ncvar, var_model_daily(ivarmod,:,:,:,:))
     ELSE
        WRITE(*,*) "----------------------------------------------------"
        WRITE(*,*) 'Problem with the variables! ==> Exit'
        char_err_status = NF90_STRERROR(il_err)
        WRITE(*,*) char_err_status
        WRITE(*,*) "----------------------------------------------------"
        STOP
     END IF
     ! Just checking.
     WRITE(unitFile_check,*) TRIM(model_variables(ivarmod)),'(ilon_check,klat_check,llev_check,ntime)=', &
          var_model_daily(ivarmod,ilon_check,klat_check,llev_check,ntime)
  END DO
  il_err = NF90_CLOSE(il_ncfile)
  mean_5D_model(:,:,:,:,:)=0.
  val_valid_or_not(1)=rwrong_val
  DO ivarmod=1,nvar_model
     ! First: which IAGOS variable do we have to take as reference for the mask?
     ! If the current model variable is amongst the IAGOS variables too, then we use its observed twin.
     ! Else, we take IAGOS ozone as reference (unless an observed variable directly has to do with
     ! the current one, like the stratospheric ozone tracer O3S for ozone or the model's NOy_ACACIA for the Nitrogen_reactive_species).
     IF (ANY(list_variables==TRIM(list_var_model(ivarmod)))) THEN ! Find which one it is.
        ivar_masking=which_one_str(list_var_model(ivarmod), list_variables, nvariables)
     ELSE IF (TRIM(list_var_model(ivarmod)) .EQ. TRIM(c_NOy_ACACIA) .AND. ANY(list_variables==c_Nitrogen_reactive_species)) THEN
        ivar_masking=ivarNOy
     ELSE
        ivar_masking=ivarO3
     END IF
     ! Another case for the mask as Nitrogen_reactive_species: NOx, HNO3 and PAN (the input and output with the NOy mask
     ! then have another name, though it has the same name in the model output file).
     IF (ANY(list_variables==c_Nitrogen_reactive_species)) THEN
        IF (TRIM(list_var_model(ivarmod)) == "NO_mask_NOy" .OR. &
             TRIM(list_var_model(ivarmod)) == "NO2_mask_NOy" .OR. &
             TRIM(list_var_model(ivarmod)) == "NOx_mask_NOy" .OR. &
             TRIM(list_var_model(ivarmod)) == "HNO3_mask_NOy" .OR. &
             TRIM(list_var_model(ivarmod)) == "PAN_mask_NOy" .OR. &
             TRIM(list_var_model(ivarmod)) == "O3_mask_NOy") THEN
           ivar_masking=ivarNOy
        END IF
     END IF
     WRITE(unitFile_check,*) 'ivarmod=',ivarmod
     WRITE(unitFile_check,*) 'ivar_masking,list_variables(ivar_masking)=',ivar_masking,list_variables(ivar_masking)
     WRITE(*,*) 'ivarmod=',ivarmod, list_var_model(ivarmod), &
          model_units_output(ivarmod+ndim_model)
     WRITE(*,*) 'ivar_masking,list_variables(ivar_masking)=',ivar_masking,list_variables(ivar_masking)
     WRITE(*,*) 'model_conv_factor(ivarmod+ndim_model)=',model_conv_factor(ivarmod)
     DO ilay=1,nlayers
        WRITE(*,*) list_layers(ilay)
        DO i=1,nlon
           DO k=1,nlat
              DO l=1,nlev
                 ntime_sampled=0
                 ! First summing and counting.
                 DO itime=1,ntime
                    tmp_value=0.
                    ! The following integer allows us to apply a mask without any if loop.
                    n_1_or_0=INT(( SIGN(1., rcount_daily(ivar_masking,ilay,i,k,l,itime)-N_min_daily ) +1 )/2)
                    tmp_value=var_model_daily(ivarmod,i,k,l,itime)*n_1_or_0
                    mean_5D_model(ivarmod,ilay,i,k,l)=mean_5D_model(ivarmod,ilay,i,k,l)+ &
                         tmp_value
                    max_5D_model(ivarmod,ilay,i,k,l)=MAX(max_5D_model(ivarmod,ilay,i,k,l),tmp_value)
                    ntime_sampled=ntime_sampled+n_1_or_0
                    ! Need another definition for the minimum derivation while filtering.
                    ! n_1_or_0=1 if it's ok, and =1E10 if not. In this way, it's out of reach for being picked by the min function.
                    n_1_or_0=(n_1_or_0-1)*(-1000000)+1
                    tmp_value=var_model_daily(ivarmod,i,k,l,itime)*n_1_or_0
                    min_5D_model(ivarmod,ilay,i,k,l)=MIN(min_5D_model(ivarmod,ilay,i,k,l),tmp_value)
                 END DO
                 ! Averaging now.
                 n_1_or_0=INT((SIGN(1,ntime_sampled-1)+1)/2) ! Filtering from the amount of sampled days.
                 n_1_or_0 = n_1_or_0 * INT((SIGN(1., rcount(ivar_masking,ilay,i,k,l)-0.01 ) +1 )/2) ! And on the total amount of data during this month: making sure it is not tiny.
                 tmp_value=mean_5D_model(ivarmod,ilay,i,k,l)/REAL(MAX(1,ntime_sampled))
!!!!!!!! Conversion into the required output units.
                 val_valid_or_not(2)=tmp_value*model_conv_factor(ivarmod)
                 mean_5D_model(ivarmod,ilay,i,k,l)=val_valid_or_not(n_1_or_0 + 1)
                 ! Same for min.
                 tmp_value=min_5D_model(ivarmod,ilay,i,k,l)
                 val_valid_or_not(2)=tmp_value*model_conv_factor(ivarmod)
                 min_5D_model(ivarmod,ilay,i,k,l) =val_valid_or_not(n_1_or_0 + 1)
                 ! Same for max.
                 tmp_value=max_5D_model(ivarmod,ilay,i,k,l)
                 val_valid_or_not(2)=tmp_value*model_conv_factor(ivarmod)
                 max_5D_model(ivarmod,ilay,i,k,l) =val_valid_or_not(n_1_or_0 + 1)
              END DO
           END DO
        END DO
     END DO
  END DO ! END DO ivarmod=1,nvar_model
  DEALLOCATE(list_nfiles_daily)
  DEALLOCATE(rcount_daily_dry)
END SUBROUTINE MEAN_VALUES_INTERPOL
