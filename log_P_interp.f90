PROGRAM LOG_P_INTERP
! This program aims to change a model grid into another.
  USE GLOBAL_VAR_INTERPOL_MOD
  USE FUNCTIONS_MOD
  USE NETCDF
  IMPLICIT NONE
  ! Checks the existence of the pressure files. If they do not exist, then we cannot do anything.
  LOGICAL pres_file_exists
  ! Depending on the layers asked by the user, logic_PV decides whether or not the program reads the PV fields.
  ! logic_P decided whether the pressure is read from the model 3D field or calculated as sigma-P hybrid levels,
  ! thus only depending on the surface pressure field and the vertical coefficients (A_i) and (B_i).
  LOGICAL logic_PV, logic_P
  ! Output directory and output file.
  CHARACTER(LEN=300) outdir, outfile, outfile_P, outfile_PV, outfile_VMR
  ! Input files.
  CHARACTER(LEN=200) P_S_file, P_file, PV_file, VMR_file
  ! Naming the reference grid.
  CHARACTER(LEN=30) ref_grid
  ! Character strings for pressure units.
  CHARACTER(LEN=30) unit_model_P_S, unit_model_P, c_Input_unit_cmodel
  ! Index pointing to a given IAGOS file in each monthly list of files.
  INTEGER iFile
  ! Dates as integers, then as characters.
  ! Their lengths have to be the same as the variable in INT2CHAR subroutine (here: 20).
  INTEGER iyyyymm, iyyyy, iyyyymm1, iyyyymm2
  CHARACTER(LEN=30) cyyyymm, cyyyymm1, cyyyy, cnlev, cnlat, cnlon
  INTEGER ivarmod, idimmod_time, idimmod_lon, idimmod_lat, idimmod_lev
  INTEGER imodel, icol_model, icol_model_output, icol_unit_model_output, icol_model_full_name
  INTEGER istat, ilay
  ! Integers relative to the model grid dimension.
  ! Amount of days during the current month. Set at 1 if the time resolution is monthly.
  INTEGER ntime, nlon_out, nlat_out, nlev_out
  ! l,k,i and imth are used respectively for the height, latitude, longitude and time dimensions.
  ! itime is for the day into the current month.
  ! imth_this_data should be created to account for flights that take place at nighttime, at the end of the current month.
  INTEGER l, k, i, imth, itime
  ! Reals corresponding to the domain borders and to the grid horizontal resolution.
  REAL lat1, lat2, lon1, lon2, dlat, dlon
  ! Reals corresponding to the domain borders and to the ref grid horizontal resolution.
  REAL lon1_out, lon2_out, lat1_out, lat2_out, dlat_out, dlon_out
  ! Reference pressure scalar. Necessary for building the vertical grid.
  REAL P_0
  ! The following shift_lon and shift_lat are a correction for the grid,
  ! in case the gridcells are defined by their center coordinates (as for the INCA model) instead of their
  ! southwestern border (as for MOCAGE).
  REAL shift_lon, shift_lat
  ! Averaging step: indexes repairing the gridcells validated for being averaged.
  INTEGER ilon_ok, klat_ok, llev_ok
  ! Which version of the VERTICAL_INTERP subroutine we intend to use.
  ! Version 0 is straightforward, but only works from a coarse resolution to a much finer resolution.
  ! Version 1 was intended to be more general (and it works well for the fine to a much coarser resolution),
  ! but I have been spending weeks on it and losing all my time, with no progress on the project.
  ! It was awful, I was not far from burn-out, and I don't want to spend more time on it.
  INTEGER version_vertical_interp
!!!!!!! Dynamic arrays !!!!!!!!
  ! Lists of dimensions and variables. The suffix _nc is for the variable names as outputs from this routine.
  CHARACTER*30,ALLOCATABLE :: list_coord(:), list_variables(:), list_variables_nc(:)
!!!! Permanent lists.
  ! Lists of variables and dimensions as queried by the user (standard name).
  CHARACTER*30,ALLOCATABLE :: list_dim_model(:), list_var_model(:)
  ! Lists of variables and dimensions as named in the model output files (i.e. as inputs for this routine).
  CHARACTER*30,ALLOCATABLE :: model_dimensions(:), model_variables(:), model_var_output(:)
  ! Pressure, excluded from the model lists cause it is not a variable we want in output for this routine.
  ! Also the Netcdf file containing the vertical hybrid coefficients, if needed. The latter is in the ./ repertory.
  CHARACTER*30                model_P, model_P_S, filename_hybrid_coef
  ! Corresponding units in the input and output model files. Below, corresponding input-output conversion factors.
  CHARACTER*30,ALLOCATABLE :: model_units_output(:), model_units_dim_output(:), model_units_var_output(:)
  ! Permanent tabular gathering all the possible names for each required variable.
  CHARACTER*30,ALLOCATABLE :: tab_var_names(:,:)
  ! Names of the different models which variables have been registered into a reference tabular.
  ! The size of model_names is n_models_registered. +2 for header models.
  CHARACTER*30,ALLOCATABLE :: model_names(:), header_models(:)
  ! Headers for the two other tables, gathering the names of the units in input (depending on the data set)
  ! and the units that are to be given to the outputs. Consequently, they also contain the appropriate
  ! conversion factors.
  CHARACTER*30,ALLOCATABLE :: header_units_models(:)
  ! Permanent tabular gathering all the possible names for each required variable.
  CHARACTER*30,ALLOCATABLE :: tab_model_names(:,:)
  ! Permanent tabulars gathering the units in input and output files in IAGOS and models, 
  ! with the corresponding conversion factors.
  CHARACTER*30,ALLOCATABLE ::tab_var_units_models(:,:)
  ! Threshold values for each variable. Below it, the data point is filtered out.
  REAL, ALLOCATABLE :: thres_val(:)
!!!!!!!!! Gridded data. !!!!!!!!!!
  ! Gridded dimensions. The _netcdf suffix is for the dimension as it will be written in this routine output.
  REAL, ALLOCATABLE :: rlat(:), rlon(:), rlat_netcdf(:), rlon_netcdf(:)
  REAL, ALLOCATABLE :: check_lon(:), check_lat(:)
  ! Output horizontal dimensions.
  REAL, ALLOCATABLE :: rlat_out(:), rlon_out(:)
  ! Arrays needed to build the pressure grid. A and B coefficients, and surface pressure for the current month.
  ! The dimensions for P_surf are longitude, latitude and the day in the current month.
  REAL, ALLOCATABLE :: A_inter(:), B_inter(:), P_surf(:,:,:)
  ! Half and full vertical grid levels.
  REAL, ALLOCATABLE :: P_edge(:,:,:,:), P_full(:,:,:,:)
  REAL, ALLOCATABLE :: P_mid(:,:,:,:)
  ! Arrays stocking the model fields, except pressure.
  REAL, ALLOCATABLE :: PV(:,:,:,:), VMR(:,:,:,:,:)
  ! Arrays interpolated onto the new grid.
  REAL, ALLOCATABLE :: PV_out(:,:,:,:), VMR_out(:,:,:,:,:)
  REAL, ALLOCATABLE :: P_out_test(:,:,:,:)
  ! Temporary arrays interpolated onto the new vertical grid, but not the horizontal yet.
  REAL, ALLOCATABLE :: PV_P_out(:,:,:,:), VMR_P_out(:,:,:,:,:)
  REAL, ALLOCATABLE :: P_mid_P_out(:,:,:,:)
  
  ! Pressure levels that we want in output.
  REAL, ALLOCATABLE :: P_out(:)
  ! Reads all the input parameters asked in make_interpol_iagos.sh
  OPEN(14,FILE=TRIM(pdir)//'namelist_interp',FORM='formatted')
  READ(14,*) iyyyymm, iyyyy
  WRITE(*,*) 'iyyyymm, iyyyy =', iyyyymm, iyyyy
  READ(14,*) nlon, nlat, nlev
  WRITE(*,*) 'nlon, nlat, nlev =', nlon, nlat, nlev
  READ(14,*) dlon, dlat
  WRITE(*,*) 'dlon, dlat =', dlon, dlat
  READ(14,*) lon1, lon2, lat1, lat2
  WRITE(*,*) 'lon1, lon2, lat1, lat2 =', lon1, lon2, lat1, lat2
  READ(14,*) nlon_out, nlat_out, nlev_out
  WRITE(*,*) 'nlon_out, nlat_out, nlev_out =', nlon_out, nlat_out, nlev_out
  READ(14,*) dlon_out, dlat_out
  WRITE(*,*) 'dlon_out, dlat_out =', dlon_out, dlat_out
  READ(14,*) lon1_out, lon2_out, lat1_out, lat2_out
  WRITE(*,*) 'lon1_out, lon2_out, lat1_out, lat2_out=', lon1_out, lon2_out, lat1_out, lat2_out
  allocate(P_out(nlev_out))
  READ(14,*) P_out(:)
  WRITE(*,*) 'P_out(:) =', P_out(:)
  READ(14,*) logic_PV
  WRITE(*,*) 'logic_PV =',logic_PV
  READ(14,*) logic_P
  WRITE(*,*) 'logic_P =',logic_P
  READ(14,*) ndim_model
  allocate(list_dim_model(ndim_model))
  READ(14,*) list_dim_model(1:ndim_model)
  READ(14,*) nvar_model
  allocate(list_var_model(nvar_model))
  READ(14,*) list_var_model(1:nvar_model)
  READ(14,*) n_models_registered
  READ(14,*) version_vertical_interp
  READ(14,*) ref_grid
  WRITE(*,*) 'ref_grid =', TRIM(ref_grid)
  READ(14,*) cmodel
  WRITE(*,*) 'cmodel =',TRIM(cmodel)
  READ(14,*) experiment
  WRITE(*,*) 'experiment =',TRIM(experiment)
  READ(14,*) tmpdir
  WRITE(*,*) 'tmpdir =',TRIM(tmpdir)
  READ(14,*) tmpdir_there
  WRITE(*,*) 'tmpdir_there =',TRIM(tmpdir_there)
  READ(14,*) ref_dir
  WRITE(*,*) 'ref_dir =',TRIM(ref_dir)
  READ(14,*) model_output_dir
  WRITE(*,*) 'model_output_dir =', TRIM(model_output_dir)
  READ(14,*) VMR_file
  WRITE(*,*) 'VMR_file =', TRIM(VMR_file)
  READ(14,*) PV_file
  WRITE(*,*) 'PV_file =', TRIM(PV_file)
  READ(14,*) P_file
  WRITE(*,*) 'P_file =', TRIM(P_file)
  READ(14,*) P_S_file
  WRITE(*,*) 'P_S_file =', TRIM(P_S_file)
  READ(14,*) filename_hybrid_coef
  WRITE(*,*) 'filename_hybrid_coef =', TRIM(filename_hybrid_coef)
  open(unitFile_check, FILE=TRIM(tmpdir_there)//'check_log_P_interp',   &
       FORM='formatted')
  WRITE(*,*) 'ndim_model =',ndim_model
  WRITE(*,*) 'list_dim_model(1:ndim_model) =',list_dim_model(1:ndim_model)
  WRITE(*,*) 'nvar_model =',nvar_model
  WRITE(*,*) 'list_var_model(1:nvar_model) =',list_var_model(1:nvar_model)
  WRITE(*,*) 'n_models_registered =',n_models_registered
  READ(14,*) outdir
  WRITE(*,*) 'outdir =',TRIM(outdir)
  close(14)
  ! Converting the dimension sizes into characters.
  call INT2CHAR(nlev,tmpdir,cnlev)
  call INT2CHAR(nlat,tmpdir,cnlat)
  call INT2CHAR(nlon,tmpdir,cnlon)
  allocate(A_inter(nlev),B_inter(nlev))
  allocate(rlon(nlon), rlat(nlat), rlon_netcdf(nlon), rlat_netcdf(nlat))
  allocate(rlon_out(nlon_out), rlat_out(nlat_out))
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!! Reading the 2 reference tables (model variable names and their respective units),
!!!!!!!!!!! with their corresponding conversion factors.
  WRITE(unitFile_check,*) '###################################################'
  var_names_file=TRIM(tmpdir)//'table_variable_names_in_models.txt'
  open(14,FILE=TRIM(var_names_file),FORM='formatted')
  allocate(header_models(n_models_registered+2))
  allocate(model_names(n_models_registered))
  READ(14,*) header_models
  WRITE(unitFile_check,*) 'n_models_registered=', n_models_registered
  WRITE(unitFile_check,*) 'ndim_model=',ndim_model
  WRITE(unitFile_check,*) 'nvar_model=',nvar_model
  allocate(tab_model_names(nvar_model+ndim_model+2,n_models_registered+2))
  WRITE(unitFile_check,*) '############# Models variable names ###############'
  WRITE(unitFile_check,*) header_models(:)
  do ivarmod=1,(nvar_model+ndim_model+2) ! +2 stands for surface pressure and pressure. Neither accounted for in variables nor dimensions.
     READ(14,*) tab_model_names(ivarmod,:)
     WRITE(*,*) tab_model_names(ivarmod,:)
     WRITE(unitFile_check,*) tab_model_names(ivarmod,:)
  end do
  WRITE(unitFile_check,*) '###################################################'
  model_names=header_models(2:n_models_registered+1)
  WRITE(unitFile_check,*) model_names
  close(14)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Model variable units.
  var_units_file=TRIM(tmpdir)//'table_variable_units_in_models.txt'
  open(14,FILE=TRIM(var_units_file),FORM='formatted')
  allocate(header_units_models(2*n_models_registered+2)) ! 2 columns per model, and 2 additional columns.
  READ(14,*) header_units_models(:)
  allocate(tab_var_units_models(nvar_model+ndim_model+2,2*n_models_registered+2))
  WRITE(unitFile_check,*) '############# Models variable units ###############'
  WRITE(unitFile_check,*) header_units_models(:)
  do ivarmod=1,(nvar_model+ndim_model+2) ! +2 stands for surface pressure and pressure. Neither accounted for in variables nor dimensions.
     READ(14,*) tab_var_units_models(ivarmod,:)
     WRITE(*,*) tab_var_units_models(ivarmod,:)
     WRITE(unitFile_check,*) tab_var_units_models(ivarmod,:)
  end do
  WRITE(unitFile_check,*) '###################################################'
  close(14)
  ! Finds the column containing the variables full names.
  icol_model_full_name=which_one_str(c_Full_name,header_models,n_models_registered+2)
  ! Finds the right column to adopt the right names for variables and the right conversion factor.
  WRITE(*,*) experiment,header_models
  if (which_one_str(experiment,header_models,n_models_registered+2) .ne. 0) then
     icol_model=which_one_str(experiment,header_models,n_models_registered+2)
  else
     icol_model=which_one_str(cmodel,header_models,n_models_registered+2)
  end if
  WRITE(*,*) 'icol_model=', icol_model, header_models(icol_model)
  ! The following index corresponds to the list of models present in the table files.
  ! It is easier to use this one to find the right columns in the units table. 
  ! -1 because the first column is for the variable names entered by the user, and is not proper to a model.
  imodel=icol_model-1
  ! Finds the right column for the variables output names. Same as input, for this routine.
  icol_model_output=icol_model
  ! Finds the right column for the variables output units.
  call fill_str_char(TRIM(c_Input_unit)//'_'//TRIM(cmodel), len(TRIM(c_Input_unit)//'_'//TRIM(cmodel)), &
       30, c_Input_unit_cmodel)
  icol_unit_model_output=which_one_str(c_Input_unit_cmodel, header_units_models, 2*n_models_registered+2)
  WRITE(unitFile_check,*) 'icol_unit_model_output=',icol_unit_model_output, header_units_models(icol_unit_model_output)
  WRITE(*,*) c_Input_unit_cmodel
  WRITE(*,*) header_units_models
  WRITE(*,*) 'icol_unit_model_output=',icol_unit_model_output, header_units_models(icol_unit_model_output)
  ! The next two arrays are the dimensions/variables names as they
  ! should be found in the model output files.
  allocate(model_dimensions(ndim_model))
  allocate(model_variables(nvar_model))
  ! Units for the model variables, without accounting for the stats.
  allocate(model_units_dim_output(ndim_model))
  allocate(model_units_var_output(nvar_model))
  allocate(model_units_output(nvar_model + ndim_model))
  model_dimensions(:) = tab_model_names(1:ndim_model,icol_model_output)
  WRITE(unitFile_check,*) 'model_dimensions=', model_dimensions
  model_variables(:)  = tab_model_names((ndim_model+1):(ndim_model+nvar_model),icol_model_output)
  ! This program does not have to modify the variable and coordinates names.
  model_var_output=model_variables
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
  WRITE(unitFile_check,*) 'model_var_output=',model_var_output
  ! Here, model_var_output is the same as the input (the homogenization step 
  ! through the variables and coordinates names comes later, in the interpol_IAGOS program).
  ! Same for coordcoord_names_output.
  coord_names_output=tab_model_names(1:ndim_model,icol_model)
  model_units_dim_output=tab_var_units_models(1:ndim_model,icol_unit_model_output)
  model_units_var_output= &
       tab_var_units_models((ndim_model+1) : (ndim_model+nvar_model), icol_unit_model_output)
  WRITE(unitFile_check,*) 'model_units_dim_output========', model_units_dim_output
  WRITE(unitFile_check,*) 'model_units_var_output========', model_units_var_output
  model_units_output(1:ndim_model) = model_units_dim_output
  model_units_output( (1+ndim_model) : (nvar_model+ndim_model)) = model_units_var_output
  WRITE(unitFile_check,*) 'model_units_output========', model_units_output
  ! Units for pressure.
  unit_model_P_S = model_units_dim_output(irow_model_P_S)
  unit_model_P   = model_units_dim_output(irow_model_P  )
  ! Defines the name for the time record dimension inquiry.
  idimmod_time=which_one_str(c_Time, list_dim_model, ndim_model)
  ctime_model=model_dimensions(idimmod_time)
  IF (cmodel .NE. "UKESM1.1") THEN
     ! Also finds the names for longitude and latitude in the model files.
     idimmod_lon=which_one_str(c_Longitude, list_dim_model, ndim_model)
     clon_model=model_dimensions(idimmod_lon)
     idimmod_lat=which_one_str(c_Latitude, list_dim_model, ndim_model)
     clat_model=model_dimensions(idimmod_lat)
     idimmod_lev=which_one_str(c_Vertical_grid_level, list_dim_model, ndim_model)
     clev_model=model_dimensions(idimmod_lev)
  ELSE
     clon_model='longitude'
     clat_model='latitude'
     clev_model='model_level_number'
  END IF
  WRITE(unitFile_check,*) 'ctime_model, clon_model, clat_model, clev_model=', ctime_model, clon_model, clat_model, clev_model
  ! ############## Generates the horizontal grids ###########
  ! The lon/lat_netcdf vectors are the coordinates as they appear in the output NetCDF file.
  ! They refer to the cells southwest limit if we use the MOCAGE model.
  WRITE(*,*) 'Defining the horizontal grid.'
  shift_lat=0. ; shift_lon=0.
  rlat(1)=lat1+shift_lat
  rlat_netcdf(1)=lat1
  do k=2,nlat
     rlat(k)=rlat(k-1)+dlat
     rlat_netcdf(k)=rlat_netcdf(k-1)+dlat
  end do
  WRITE(*,*) rlat(1), '°N --> ', rlat(nlat),'°N'
  rlon(1)=lon1+shift_lon
  rlon_netcdf(1)=lon1
  do i=2,nlon
     rlon(i)=rlon(i-1)+dlon
     rlon_netcdf(i)=rlon_netcdf(i-1)+dlon
  end do
  WRITE(*,*) rlon_netcdf(1), '°E --> ', rlon_netcdf(nlon),'°E'
  ! Now: the output grid.
  WRITE(*,*) 'Defining the output horizontal grid.'
  rlat_out(1)=lat1_out
  do k=2,nlat_out
     rlat_out(k)=rlat_out(k-1)+dlat_out
  end do
  WRITE(*,*) rlat_out(1), '°N --> ', rlat_out(nlat_out),'°N'
  rlon_out(1)=lon1_out
  do i=2,nlon_out
     rlon_out(i)=rlon_out(i-1)+dlon_out
  end do
  WRITE(*,*) rlon_out(1), '°E --> ', rlon_out(nlon_out),'°E'
  ! ############## Generates the vertical grid #############
  WRITE(*,*) nlev, TRIM(tmpdir)
  call INT2CHAR(nlev,tmpdir,cnlev)
  WRITE(*,*) 'Reading the vertical grid coefficients.'
  WRITE(*,*) cnlev,'vertical levels'
  if (.not. logic_P) then
     il_err = NF90_OPEN(TRIM(filename_hybrid_coef), NF90_NOWRITE, il_ncfile)
     il_err = NF90_INQ_VARID (il_ncfile, 'hyam', il_ncvar)
     il_err = NF90_GET_VAR   (il_ncfile, il_ncvar, A_inter(:))
     il_err = NF90_INQ_VARID (il_ncfile, 'hybm', il_ncvar)
     il_err = NF90_GET_VAR   (il_ncfile, il_ncvar, B_inter(:))
     il_err = NF90_CLOSE (il_ncfile)
  end if
  ! Just for testing the fields on a mid-latitude gridcell.
  ilon_check=int(nlon/2) ; klat_check=int(nlat*4/5) ; llev_check=int(3*nlev/5)
  CALL INT2CHAR(ilon_check, tmpdir, clon_check)
  CALL INT2CHAR(klat_check, tmpdir, clat_check)
  CALL INT2CHAR(llev_check, tmpdir, clev_check)
  call INT2CHAR(iyyyymm,tmpdir,cyyyymm) ! Converts the integer date into characters
  iyyyy=nint(real(iyyyymm)/100)
  call INT2CHAR(iyyyy,tmpdir,cyyyy) ! Converts the integer date into characters

  data_dir_date = TRIM(cyyyy)//'/'//TRIM(cyyyymm)//'/'
  data_dir = TRIM(data_dir_cst)//TRIM(data_dir_date)
!!!!!! Reading the surface pressure fields.
  if (logic_P) then
     filename_P = P_file
  else
     filename_P = P_S_file
  end if
  ! filename_Psurf=TRIM(model_output_dir)//TRIM(cyyyy)//'/P_S_'//TRIM(cyyyymm)//'.nc'
  WRITE(*,*) 'filename_P=', TRIM(filename_P)
  il_err = NF90_OPEN(TRIM(filename_P), NF90_NOWRITE, il_ncfile)
  IF (il_err == 0) THEN
     pres_file_exists=.TRUE.
  ELSE
     pres_file_exists=.FALSE.
  END IF
  il_err = NF90_INQ_DIMID (il_ncfile, TRIM(ctime_model), il_ncdim_time)
  il_err = NF90_INQUIRE_DIMENSION (il_ncfile, il_ncdim_time, LEN = ntime)
  WRITE(*,*) 'nlon, nlat, nlev, ntime =', nlon, nlat, nlev, ntime
!!!!!! Generating the 3D pressure fields, both for half and full levels.
  ALLOCATE( P_mid(nlon, nlat, nlev, ntime) )
  IF (.NOT. logic_P) THEN
     ALLOCATE( P_surf(nlon,nlat,ntime) )
     ALLOCATE( P_edge(nlon,nlat,nlev,ntime) )
     il_err = NF90_INQ_VARID (il_ncfile, TRIM(model_P_S), il_ncvar)
     IF (il_err .NE. 0) THEN
        char_err_status = NF90_STRERROR(il_err)
        WRITE(*,*) char_err_status
     END IF
     il_err = NF90_GET_VAR (il_ncfile, il_ncvar, P_surf)
     il_err = NF90_CLOSE (il_ncfile)
     WRITE(*,*) 'filename_Psurf=', filename_P
     ! Checking the daily surface pressures.
     DO itime=1,ntime
        WRITE(unitFile_check,*) 'P_surf(nlon/2,4/5*nlat,itime)=', &
             P_surf(ilon_check, klat_check, itime)
     END DO
     IF (pres_file_exists) THEN
        !! Fixing the reference pressure P_0. 1E5 Pa for MOZART3, 1 Pa elsewise.
        IF (TRIM(cmodel) .EQ. 'MOZART3') THEN
           P_0 = 100000.
        ELSE
           P_0 = 1.
        END IF
        CALL COMPUTE_PRESSURES(nlon, nlat, nlev, ntime, P_surf(:,:,:), P_0, &
             A_inter(:), B_inter(:), P_edge(:,:,:,:), P_mid(:,:,:,:))
        IF (TRIM(cmodel) .EQ. 'MOZART3') THEN ! In this case, the full pressure levels are directly derived linearly from the surface pressure.
           P_mid(:,:,:,:) = P_edge(:,:,:,:)
        END IF
     ELSE
        P_edge(:,:,:,:) = rwrong_val
        P_mid (:,:,:,:) = rwrong_val
     END IF
     DEALLOCATE(P_surf, P_edge)
  END IF ! .NOT. logic_P
  IF (.NOT. pres_file_exists) THEN
     WRITE(*,*) filename_P
     WRITE(*,*) '---> Pressure field missing: cannot interpolate the output onto another grid for this month.'
     WRITE(unitFile_check,*) filename_P
     WRITE(unitFile_check,*) '---> Pressure field missing: cannot interpolate the output onto another grid for this month.'
  ELSE IF (logic_P) THEN
!!!!!! Reading the pressure fields (at the middle of grid cells) in the model file.
     il_err = NF90_INQ_VARID(il_ncfile, TRIM(model_P), il_ncvar)
     IF (il_err .NE. 0) THEN
        WRITE(*,*) NF90_STRERROR(il_err), "model_P =", model_P
     END IF
     il_err = NF90_GET_VAR(il_ncfile, il_ncvar, P_mid)
     il_err = NF90_CLOSE  (il_ncfile)
     WRITE(unitFile_check,*) 'filename_P=', filename_P
     ! Checking the daily pressures.
     DO itime=1, ntime
        WRITE(unitFile_check,*) 'P_mid(nlon/2,4/5*nlat,nlev,itime)=', P_mid(ilon_check,klat_check,nlev,itime)
     END DO
  END IF
!!!!!! Reading the model chemical fields.
  il_err = NF90_OPEN(TRIM(VMR_file), NF90_NOWRITE, il_ncfile)
  WRITE(*,*) 'ntime = ', ntime, il_err
  IF (il_err == 0) THEN
     WRITE(*,*) 'Reading the model fields in ', TRIM(VMR_file)
     WRITE(unitFile_check,*) 'Reading the model fields in ', TRIM(VMR_file)
     ! Checking the coordinates in the model file.
     ALLOCATE( check_lon(nlon), check_lat(nlat) )
     il_err = NF90_INQ_VARID(il_ncfile, TRIM(clon_model), il_ncdim_lon)
     il_err = NF90_GET_VAR  (il_ncfile,     il_ncdim_lon, check_lon(:))
     il_err = NF90_INQ_VARID(il_ncfile, TRIM(clat_model), il_ncdim_lat)
     il_err = NF90_GET_VAR  (il_ncfile,     il_ncdim_lat, check_lat(:))
     IF (check_lat(2)-check_lat(1) <= 0. .OR. & ! northward latitude axis
       P_mid(ilon_check,klat_check,2,1)-P_mid(ilon_check,klat_check,1,1) <= 0. .OR. &! downward vertical axis
        ABS(check_lon(nlon))-ABS(check_lon(1)) >= 300) THEN ! longitude axis centered around 0 degree
        WRITE(*,*) "ERROR: reversed axis in the model file: ", TRIM(cyyyymm)
        WRITE(*,*) "P_mid(1), P_mid(2)=", P_mid(ilon_check,klat_check,1,1), P_mid(ilon_check,klat_check,2,1)
        WRITE(*,*) "check_lon(1), check_lon(nlon)=", check_lon(1), check_lon(nlon)
        WRITE(*,*) "check_lat(1), check_lat(2)   =", check_lat(1), check_lat(2)
        STOP
     END IF
     il_err = NF90_INQ_VARID(il_ncfile,TRIM(clev_model),il_ncdim_lev)
     ALLOCATE(VMR(nvar_model,nlon,nlat,nlev,ntime))
     DO ivarmod=1,nvar_model
        il_err=NF90_INQ_VARID(il_ncfile,TRIM(model_variables(ivarmod)),il_ncvar)
        WRITE(*,*) ivarmod, model_variables(ivarmod), il_ncvar, "err=", il_err
        IF (il_err.EQ.0) THEN
           il_err=NF90_GET_VAR(il_ncfile,il_ncvar,VMR(ivarmod,:,:,:,:))
        ELSE
           WRITE(*,*) 'OOOOOOOOOOOOOOOOOOOOOOOPS var problem! ==> Exit'
           STOP
        END IF
        ! Just checking.
        WRITE(unitFile_check,*) TRIM(model_variables(ivarmod)),'(ilon_check,klat_check,llev_check,ntime)=', &
             VMR(ivarmod,ilon_check,klat_check,llev_check,ntime)
     END DO
     il_err=NF90_CLOSE(il_ncfile)
  ELSE
     WRITE(*,*) TRIM(VMR_file)
     WRITE(*,*) '--> Model file not found.'
     WRITE(unitFile_check,*) TRIM(VMR_file)
     WRITE(unitFile_check,*) '--> Model file not found. NF90_INQ_VARID returned',il_err
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
!!!!!! Reading (or not) the potential vorticity fields.
  ! We hopp this step if we ask the program to treat the layer 'all_layers' or/and 'UTLS' only,
  ! for the climatology maps, for example. Or for comparing with gridded IASI measurements
  ! between 9 and 12 km.
  ! To be used later, for the separation between the different layers.
  IF (logic_PV) THEN
     ALLOCATE(PV(nlon,nlat,nlev,ntime))
     WRITE(*,*) 'Reading the PV files.'
     il_err=NF90_OPEN(TRIM(PV_file),NF90_NOWRITE,il_ncfile)
     il_err=NF90_INQ_VARID(il_ncfile,TRIM(c_PV),il_ncvar)
     IF (il_err.NE.0) THEN
        WRITE(unitFile_check,*) 'variable PV not found. NF90_INQ_VARID returned',il_err
     END IF
     il_err=NF90_GET_VAR(il_ncfile,il_ncvar,PV)
     il_err=NF90_CLOSE(il_ncfile)
     ! The following loop checks the intramonthly PV variations.
     DO itime=1,ntime
        WRITE(*,*) 'PV('//TRIM(clon_check)//','//TRIM(clat_check)//','//TRIM(clev_check)//',','itime)=', &
             PV(ilon_check,klat_check,llev_check,itime), 'itime=', itime
        WRITE(unitFile_check,*) 'PV(ilon_check,klat_check,llev_check,itime)=', &
             PV(ilon_check,klat_check,llev_check,itime), 'itime=', itime
     END DO
  END IF
  WRITE(*,*) '##################################'
  WRITE(*,*) 'Calling VERTICAL_INTERP subroutine'
  WRITE(*,*) '##################################'
  IF (logic_PV) THEN
     ALLOCATE(PV_P_out(nlon,nlat,nlev_out,ntime))
  END IF
  ALLOCATE(P_mid_P_out(nlon,nlat,nlev_out,ntime))
  ALLOCATE(VMR_P_out(nvar_model,nlon,nlat,nlev_out,ntime))
  DO ivarmod=1,nvar_model
     DO i=1,nlon
        DO k=1,nlat
           DO itime=1,ntime
                 CALL VERTICAL_INTERP(nlev,nlev_out,                &
                      P_mid(i,k,:,itime), VMR(ivarmod,i,k,:,itime), &
                      P_out(:), VMR_P_out(ivarmod,i,k,:,itime), version_vertical_interp)
                 IF (MINVAL(VMR_P_out(ivarmod,i,k,:,itime), MASK=VMR_P_out(ivarmod,i,k,:,itime) .NE. rwrong_val) < 0.) THEN
                    WRITE(*,*) "WARNING VERTICAL INTERP "//TRIM(cyyyymm)//" "//TRIM(list_var_model(ivarmod))
                    WRITE(*,*) minval(VMR_P_out(ivarmod,i,k,:,itime), MASK=VMR_P_out(ivarmod,i,k,:,itime) .NE. rwrong_val)
                    STOP
                 END IF
                 IF (logic_PV) THEN
                    CALL VERTICAL_INTERP(nlev,nlev_out,             &
                         P_mid(i,k,:,itime), PV(i,k,:,itime),       &
                         P_out(:), PV_P_out(i,k,:,itime), version_vertical_interp )
                 END IF
                 CALL VERTICAL_INTERP(nlev,nlev_out,                &
                      P_mid(i,k,:,itime), P_mid(i,k,:,itime),       &
                      P_out(:), P_mid_P_out(i,k,:,itime), version_vertical_interp )
                 IF (MINVAL(P_mid_P_out(i,k,:,itime), MASK=P_mid_P_out(i,k,:,itime) .NE. rwrong_val) < 0.) THEN
                    WRITE(*,*) "WARNING VERTICAL INTERP "//TRIM(cyyyymm)//" "//TRIM(list_var_model(ivarmod))
                    WRITE(*,*) minval(P_mid_P_out(i,k,:,itime), MASK=P_mid_P_out(i,k,:,itime) .NE. rwrong_val)
                    STOP
                 END IF
           END DO
        END DO
     END DO
  END DO
  DEALLOCATE(P_mid, VMR)
  IF (logic_PV) THEN
     DEALLOCATE(PV)
  END IF
  ! Writes the new model fields.
  outfile_VMR=TRIM(outdir)//'tmp_'//TRIM(experiment)//'_'//TRIM(cyyyymm)//'.nc'
  CALL TONETCDF(VMR_P_out(:,:,:,:,:), nlon, nlat, nlev_out, ntime, &
       coord_names_output(:), nvar_model, model_var_output(:),     &
       model_units_var_output(:),                                  &
       rlon(:), rlat(:), cyyyymm, rwrong_val, outfile_VMR)
  IF (logic_PV) THEN
     outfile_PV=TRIM(outdir)//'tmp_PV_'//TRIM(cyyyymm)//'.nc'
     CALL TONETCDF(PV_P_out(:,:,:,:), nlon, nlat, nlev_out, ntime, &
          coord_names_output(:), 1, c_PV, c_PVU,                   &
          rlon(:), rlat(:), cyyyymm, rwrong_val, outfile_PV)
  END IF
  outfile_P=TRIM(outdir)//'tmp_P_'//TRIM(cyyyymm)//'.nc'
  CALL TONETCDF(P_mid_P_out(:,:,:,:), nlon, nlat, nlev_out, ntime, &
       coord_names_output(:), 1, model_P, unit_model_P,            &
       rlon(:), rlat(:), cyyyymm, rwrong_val, outfile_P)
  WRITE(*,*) '### log_P_interpol: done. ###'
  WRITE(unitFile_check,*) '### log_P_interpol: done. ###'
  CLOSE(unitFile_check)

END PROGRAM LOG_P_INTERP
