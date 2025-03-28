MODULE FUNCTIONS_MOD
  IMPLICIT NONE

CONTAINS
  ! which_one_str locates, if it exists, a given
  ! string character in a vector of string characters.
  ! If it does not exist, then it yields 0.
  INTEGER FUNCTION WHICH_ONE_STR(str,vect,n_size_vect)
    IMPLICIT NONE
    CHARACTER*30, INTENT(IN) :: str
    INTEGER,      INTENT(IN) :: n_size_vect
    CHARACTER*30, INTENT(IN) :: vect(n_size_vect)
    INTEGER i_search
    i_search=1
    DO WHILE ( (TRIM(vect(i_search)) .NE. TRIM(str)) .AND. &
         (i_search .LE. n_size_vect) )
       i_search=i_search+1
    END DO
    ! Case str not found: i_search = n_size_vect + 1.
    ! We bring it back to 0.
    i_search = MODULO(i_search, (n_size_vect+1) )
    which_one_str = i_search
  END FUNCTION WHICH_ONE_STR

  ! All the ...ok functions below are classical operations
  ! that exclude the values below a given threshold.
  ! Useful if the missing values are like -999, for variables
  ! like mixing ratios or temperatures, etc.
  REAL FUNCTION MEANOK(vect, n_size_vect, min_val)
    USE GLOBAL_VAR_INTERPOL_MOD, ONLY : rwrong_val
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_size_vect
    REAL   , INTENT(IN) :: vect(n_size_vect)
    REAL   , INTENT(IN) :: min_val
    REAL mu
    INTEGER n_valid
    mu = SUM(vect, vect(:).GT.min_val)
    n_valid = COUNT(vect(:).GT.min_val)
    IF (n_valid .GT. 0) THEN
       mu = mu / n_valid
    ELSE
       mu = rwrong_val
    END IF
    meanok = mu
  END FUNCTION MEANOK
  

  REAL FUNCTION SUMOK(vect, n_size_vect, min_val)
    USE GLOBAL_VAR_INTERPOL_MOD, ONLY : rwrong_val
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_size_vect
    REAL   , INTENT(IN) :: vect(n_size_vect)
    REAL   , INTENT(IN) :: min_val
    REAL Sigma
    INTEGER n_valid
    Sigma = SUM(vect, vect(:).GT.min_val)
    n_valid = COUNT(vect(:).GT.min_val)
    IF (n_valid .EQ. 0) THEN
       Sigma = rwrong_val
    END IF
    sumok = Sigma
  END FUNCTION SUMOK

  
  ! The 2 following functions are the same as before,
  ! but with several dimensions.
  REAL FUNCTION MEANOK_4D(array, n1, n2, n3, n4, min_val)
    USE GLOBAL_VAR_INTERPOL_MOD, ONLY : rwrong_val
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n1, n2, n3, n4
    REAL   , INTENT(IN) :: array(n1, n2, n3, n4)
    REAL   , INTENT(IN) :: min_val
    REAL mu
    INTEGER n_valid
    mu = SUM(array, array(:,:,:,:) .GT. min_val)
    n_valid = COUNT(array(:,:,:,:) .GT. min_val)
    IF (n_valid .GT. 0) THEN
       mu = mu / n_valid
    ELSE
       mu = rwrong_val
    END IF
    meanok_4d = mu
  END FUNCTION MEANOK_4D
  

  REAL FUNCTION SUMOK_4D(array, n1, n2, n3, n4, min_val)
    USE GLOBAL_VAR_INTERPOL_MOD, ONLY : rwrong_val
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n1, n2, n3, n4
    REAL   , INTENT(IN) :: array(n1, n2, n3, n4)
    REAL   , INTENT(IN) :: min_val
    REAL Sigma
    INTEGER n_valid
    Sigma = SUM(array, array(:,:,:,:) .GT. min_val)
    n_valid = COUNT(array(:,:,:,:) .GT. min_val)
    IF (n_valid .EQ. 0) THEN
       Sigma = rwrong_val
    END IF
    sumok_4d = Sigma
  END FUNCTION SUMOK_4D
  

  SUBROUTINE CLASSIFY(vect, ntot, x)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: ntot
    REAL, INTENT(IN) :: vect(ntot)
    REAL, INTENT(OUT) :: x(ntot)
    ! INTEGER, INTENT(OUT) :: indice(ntot)
    INTEGER :: i, itemp, itest
    REAL :: temp
    ! do i=1,ntot
    !    indice(i) = i
    ! END DO
    x(:) = vect(:)
    itest = 1
    DO WHILE (itest .NE. 0)
       itest=0
       DO i = 1, (ntot-1)
          IF (x(i) > x(i+1)) THEN

             itest = itest+1

             temp = x(i)
             x(i) = x(i+1)
             x(i+1) = temp
             ! itemp = indice(i)
             ! indice(i) = indice(i+1)
             ! indice(i+1) = itemp
          END IF
       END DO
    END DO
    ! IF (itest == 0) then 
    !    exit
    ! END IF
  END SUBROUTINE CLASSIFY

  REAL FUNCTION QUANTILEOK(vect, n_size_vect, min_val, quant_val)
    USE GLOBAL_VAR_INTERPOL_MOD, ONLY : rwrong_val
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: n_size_vect
    REAL   , INTENT(IN) :: vect(n_size_vect)
    REAL   , INTENT(IN) :: min_val
    REAL   , INTENT(IN) :: quant_val
    REAL, ALLOCATABLE :: vect_incr(:), vect_incr_ok(:)
    INTEGER i, n_ok, i_ok, i_quantile
    REAL r_quantile, delta_x, delta_y

    ALLOCATE(vect_incr(n_size_vect))
    CALL CLASSIFY(vect, n_size_vect, vect_incr)
    ! Excluding the missing values.
    n_ok=COUNT(vect_incr .GT. min_val)
    ! Then we build an array of these valid values only.
    ALLOCATE(vect_incr_ok(n_ok))
    i_ok=0
    DO i=1,n_size_vect
       IF (vect_incr(i) .GT. min_val) THEN
          i_ok = i_ok + 1
          vect_incr_ok(i_ok) = vect_incr(i)
       END IF
    END DO
    DEALLOCATE(vect_incr)
    IF (real(n_ok) .LT. MIN(quant_val,1.- quant_val)**(-1)) THEN
       ! Sample too small. Returns the mean value.
       quantileok=meanok(vect,n_size_vect,min_val)
    ELSE
       ! Next step: we calculate the required quantile.
       ! First, we define its location in the array.
       r_quantile = quant_val*REAL(n_ok)
       ! Last: we interpolate if the sample size is not multiple of the quantile.
       IF (MODULO(INT(1/quant_val),n_ok) .EQ. 0) THEN
          i_quantile = ANINT(r_quantile) ! r_quantile is INTEGER, so it equals i_quantile.
          quantileok= vect_incr_ok(i_quantile)
       ELSE
          delta_x = r_quantile - AINT(r_quantile)
          ! Temporary i_quantile, as the INTEGER part of r_quantile.
          i_quantile = AINT(r_quantile)
          ! Slope from i_quantile to i_quantile + 1.
          delta_y = vect_incr_ok(i_quantile+1) - vect_incr_ok(i_quantile)
          ! Now we use this slope to interpolate between these two values.
          quantileok = vect_incr_ok(i_quantile) + delta_y*delta_x
       END IF
    END IF
  END FUNCTION QUANTILEOK

  ! This subroutine fills a character string with empty spaces 
  ! with the required string length.
  SUBROUTINE FILL_STR_CHAR(str_in, nchar_in, nchar_out, str_out)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nchar_in, nchar_out
    CHARACTER(LEN=nchar_in), INTENT(IN)  :: str_in
    INTEGER ichar
    CHARACTER(LEN=nchar_out), INTENT(OUT) :: str_out
    str_out(1:nchar_in) = str_in(1:nchar_in)
    DO ichar=(nchar_in+1),nchar_out
       str_out(ichar:ichar)=' '
    END DO
  END SUBROUTINE FILL_STR_CHAR

  ! This subroutine is written to automatize the conversion of
  ! an INTEGER into a string of characters.
  ! The input is a n-sized vector.
  SUBROUTINE INT2CHAR(ivar,pdir,cvar)
    IMPLICIT NONE
    INTEGER i,iabs,ilog10,nchar
    INTEGER unitFileConv
    REAL rlog10,rabs,rilog10
    CHARACTER(LEN=200) pdir,convert_file
    INTEGER ivar
    CHARACTER(LEN=30) cvar   ! final result
    CHARACTER(LEN=30) FMTvar ! An adaptable format, depending on the size of the variable
    PARAMETER(UNITFILECONV=14)
    convert_file = TRIM(pdir)//'tmp_convert_int2char.txt'
    OPEN(UNIT=unitFileConv,FILE=convert_file,FORM='formatted')
    ! Computes the length 'nchar' of the character string,
    ! using the 10-logarithm function
    IF (ivar .NE. 0) THEN
       iabs=ABS(ivar)
       rabs=REAL(iabs)
       rlog10=LOG10(rabs)
       rilog10=AINT(rlog10)
       ilog10=INT(rilog10)
    ELSE
       ilog10=0
    END IF
    IF (ivar .GE. 0) THEN
       nchar=ilog10+1
    ELSE
       nchar=ilog10+2  !! for the addition of the minus sign if negative
    END IF
    ! Prints the INTEGER in a temporary file, then reads it by filling
    ! another variable, this time string of characters, with an adaptable length.
    WRITE(FMTvar,'("(i", i0,".", i0,")")') nchar, (nchar-1)
    WRITE(unitFileConv,FMTvar) ivar
    CLOSE(unitFileConv)
    OPEN(UNIT=unitFileConv,FILE=convert_file, &
         FORM='formatted',ACTION='read')
    READ(UNIT=unitFileConv,FMT='(A)') cvar
    CLOSE(unitFileConv)
  END SUBROUTINE INT2CHAR

  ! This subroutine is written to automatize the conversion of
  ! a REAL into a string of characters.
  ! Careful: with more than 8 significant numbers,
  ! some rounds are made and the exact REAL value cannot be restored.

  ! I added a truncation at the third decimal, in order to avoid problems that can happen with too long numbers.

  SUBROUTINE REAL2CHAR(rvar,pdir,cvar)
    IMPLICIT NONE
    REAL rabs,rlog10,rilog10,rvar_tmp,rvar_3dec
    INTEGER i,ndec,nchar,ilog10,nchartot,unitFileConv
    ! ndec=nb of decimals, and nchar=nb of char before the "." (including the - sign)
    ! nchartot designs the total length, including - sign and the "."
    CHARACTER(LEN=100) pdir, convert_file
    REAL rvar
    CHARACTER(LEN=30) cvar ! final result
    CHARACTER(LEN=30) FMTvar  ! An adaptable format, depending on the size of the variable
    PARAMETER(UNITFILECONV=11)

    convert_file = TRIM(pdir)//'tmp_convert_real2char.txt'
    ! Computes the length of the character string nchar,
    ! using the 10-logarithm function      
    IF (rvar .NE. 0) THEN
       rabs=ABS(rvar)
       rlog10=LOG10(rabs)
       rilog10=AINT(rlog10)
       ilog10=INT(rilog10)
    ELSE
       ilog10=0
    END IF
    IF (rvar .GT. 0) THEN
       nchar=ilog10+1
    ELSE
       nchar=ilog10+2 !! for the addition of the minus sign if negative
    END IF
    ! Computes the length of the character string ndec
    ndec=0
    rvar_tmp=rvar
    IF (ABS(rvar).GT.0.001) THEN
       DO WHILE (rvar_tmp.NE.AINT(rvar_tmp))
          rvar_tmp=rvar_tmp*10
          ndec=ndec+1
          IF (ndec.EQ.3) THEN
             rvar_tmp=AINT(rvar_tmp)
          END IF
       END DO
       IF (ndec.EQ.3) THEN
          rvar_3dec=rvar_tmp/(10**ndec) ! truncates after the third decimal
       ELSE IF (ndec.LT.3) THEN
          rvar_3dec=rvar
       END IF
       IF (nchar.LE.1) THEN
          nchar=2
       END IF
       nchartot=nchar+ndec+1 ! +1 for the '.'
    END IF
    ! Prints the REAL in a temporary file, then reads it by filling
    ! another variable, this time string of characters, with an adaptable length.
    IF (ABS(rvar).GT.0.001) THEN
       OPEN(UNIT=unitFileConv,FILE=convert_file,FORM='formatted')
       WRITE(FMTvar,'("(f", i0,".", i0,")")') nchartot, ndec
       WRITE(unitFileConv,FMTvar) rvar_3dec
       CLOSE(unitFileConv)
       OPEN(UNIT=unitFileConv,FILE=convert_file, &
            FORM='formatted',ACTION='read')
       READ(UNIT=unitFileConv,FMT='(A)') cvar
       CLOSE(unitFileConv)
       IF (cvar.EQ.'***') THEN
          WRITE(*,*) cvar
          WRITE(*,*) i
          WRITE(*,*) 'rvar, rabs=',rvar,',',rabs
          WRITE(*,*) 'rlog10=',rlog10

          WRITE(*,*) 'rilog10=',rilog10
          WRITE(*,*) 'ilog10=',ilog10
          WRITE(*,*) 'nchar=',nchar
          WRITE(*,*) 'ndec=',ndec
          STOP
       END IF
    ELSE
       cvar='0.'
    END IF
  END SUBROUTINE REAL2CHAR

  ! This routine prints an input data field into a NetCDF file.
  ! It is called for each month. Thus, this one it does not make time series.

  SUBROUTINE TONETCDF(field, nlon, nlat, nlev, ntime, &
        dim_name, nvar, var_name, units_name,         &
        rlon, rlat, cyyyymm, rwrong_val, ncfile)

    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
!!!!!!!! Input parameters !!!!!!!!
    ! 3-dimension sizes and amount of variables to write.
    INTEGER, INTENT(IN) :: nlon, nlat, nlev, ntime, nvar
    ! Data array to write.
    REAL, INTENT(IN) :: field(nvar,nlon,nlat,nlev,ntime)
    ! Current month.
    CHARACTER*6, INTENT(IN) :: cyyyymm
    ! Default missing value.
    REAL, INTENT(IN) :: rwrong_val
    ! Output file name.
    CHARACTER*300, INTENT(IN) :: ncfile
    ! dim_name has 4 dimensions: space and time, though time dimension has
    ! 1 element only. 
    CHARACTER*30, INTENT(IN) :: var_name(nvar), units_name(nvar), dim_name(4)
    ! Horizontal dimensions values. No need to import the vertical levels, 
    ! their values are from 1 to nlev, simply.
    REAL, INTENT(IN) :: rlon(nlon), rlat(nlat)
    CHARACTER*48 cltori, cldat
!!!!!!!!!! Indexes for the NetCDF file manipulation.
    INTEGER il_err
    ! File and dimensions ID.
    INTEGER il_ncfile, il_dlev, il_dlat, il_dlon, il_dtime
    INTEGER il_flev, il_flat, il_flon, il_ftime
    ! One index for each variable...
    INTEGER il_nc_field(nvar)
    INTEGER ndim
    PARAMETER(NDIM=4)
    ! ... and for each dimension.
    INTEGER ila_ncdims(ndim) 
    INTEGER ila_ncstart(ndim), ila_nccount(ndim)
!!!!!!!!!!
    ! k is conventionally the index running the vertical grid levels.
    INTEGER k, ivar, itime
    ! Vertical grid levels.
    REAL rlev(nlev)
    ! Time axis.
    REAL, ALLOCATABLE :: rl_daynum(:)
    !
    !** OPEN NETCDF FILE
    !
    il_err = NF_CREATE(TRIM(ncfile), NF_WRITE, il_ncfile)
    !
    !** DEFINE DIMENSIONS
    !
    il_err = NF_DEF_DIM(il_ncfile, TRIM(dim_name(4)), nlev, il_dlev)
    il_err = NF_DEF_DIM(il_ncfile, TRIM(dim_name(3)), nlat, il_dlat)
    il_err = NF_DEF_DIM(il_ncfile, TRIM(dim_name(2)), nlon, il_dlon)
    il_err = NF_DEF_DIM(il_ncfile, TRIM(dim_name(1)), NF_UNLIMITED, il_dtime)

    !** DEFINE DIMENSION VARIABLES
    ! 
    il_err = NF_DEF_VAR(il_ncfile, TRIM(dim_name(1)), NF_REAL, &
         1, il_dtime, il_ftime)
    il_err = NF_DEF_VAR(il_ncfile, TRIM(dim_name(2)), NF_REAL, &
         1, il_dlon , il_flon)
    il_err = NF_DEF_VAR(il_ncfile, TRIM(dim_name(3)), NF_REAL, &
         1, il_dlat , il_flat)
    il_err = NF_DEF_VAR(il_ncfile, TRIM(dim_name(4)), NF_REAL, &
         1, il_dlev , il_flev)

    !
    !**  DEFINE ATTRIBUTES
    !

    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_flon, 'standard_name', &
         LEN('Longitude'), 'Longitude')
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_flon, 'units',         &
         LEN('degrees_east'), 'degrees_east')
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_flat, 'standard_name', &
         LEN('Latitude'), 'Latitude')
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_flat, 'units',         &
         LEN('degrees_north'), 'degrees_north')
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_flev, 'units',         &
         LEN('level'), 'level')
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_ftime, 'standard_name',&
         LEN('Time_axis'), 'Time_axis')
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_ftime,'axis',          &
         1, 't')
    !** DEFINE FIELD VARIABLE
    !
    ! First setting the array of dimensions indexes. 
    ila_ncdims(1) = il_dlon
    ila_ncdims(2) = il_dlat
    ila_ncdims(3) = il_dlev
    ila_ncdims(4) = il_dtime
    ! Then we define the fields for each variable.
    DO ivar=1,nvar
       write(*,*) var_name(ivar), units_name(ivar)
       il_err = NF_DEF_VAR(il_ncfile, TRIM(var_name(ivar)),                  &
            NF_REAL, ndim, ila_ncdims, il_nc_field(ivar))
       il_err = NF_PUT_ATT_TEXT(il_ncfile, il_nc_field(ivar), 'units',       &
            LEN_TRIM(units_name(ivar)), TRIM(units_name(ivar)))
       il_err = NF_PUT_ATT_real(il_ncfile, il_nc_field(ivar),'missing_value',&
            NF_REAL, 1, rwrong_val)
       il_err = NF_PUT_ATT_real(il_ncfile, il_nc_field(ivar),'_FillValue',   &
            NF_REAL, 1, rwrong_val)
    END DO
    !** END OF DEFINE MODE
    !
    il_err = NF_ENDDEF(il_ncfile)
    ! write(*,*) NF_STRERROR(il_err)

    !
    !** WRITE DATA TO FILE
    !
    ! First writing the dimension arrays.
    ! Time.
    
    ! Horizontal grid.
    il_err = NF_PUT_VAR(il_ncfile,il_flon,rlon)
    il_err = NF_PUT_VAR(il_ncfile,il_flat,rlat)
    ! Vertical grid.
    DO k = 1, nlev
       rlev(k) = REAL(k)
    END DO
    il_err = NF_PUT_VAR(il_ncfile,il_flev,rlev)
    
    ! Setting the first coordinate number for each dimension.
    ila_ncstart(1) = 1
    ila_ncstart(2) = 1
    ila_ncstart(3) = 1
    ila_ncstart(4) = 1
    ! Setting the length for each dimension axis.
    ila_nccount(1) = nlon
    ila_nccount(2) = nlat
    ila_nccount(3) = nlev
    ila_nccount(4) = ntime
    ! Time value.
    ALLOCATE(rl_daynum(ntime))
    DO itime = 1, ntime
       rl_daynum(itime) = REAL(itime)
    END DO
    il_err = NF_PUT_VAR(il_ncfile,il_ftime,ila_ncstart, &
         ila_nccount,rl_daynum)
    il_err = NF_PUT_VAR(il_ncfile,il_ftime,rl_daynum)
    DO ivar=1,nvar
       il_err = NF_PUT_VARA(il_ncfile,il_nc_field(ivar), &
            ila_ncstart, ila_nccount, FIELD(ivar,:,:,:,:))
    END DO
    !
    !** CLOSE THE NETCDF FILE
    !
    il_err = NF_CLOSE(il_ncfile)
    !
    WRITE(*,*) 'Tonetcdf: file available as ', TRIM(ncfile)
    RETURN
  END SUBROUTINE TONETCDF

  ! This routine prints an input data field into a NetCDF file.
  ! It is called for the whole time period. Thus, this one it does not make time series.

  SUBROUTINE TONETCDF_CLIMATO(field, nlon, nlat, dim_name, &
       nvar, var_name, units_name,                         &
       rlon, rlat, rwrong_val, ncfile)
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'

!!!!!!!! Input parameters !!!!!!!!
    ! 3-dimension sizes and amount of variables to write.
    INTEGER, INTENT(IN) :: nlon, nlat, nvar
    ! Data array to write.
    REAL, INTENT(IN) :: field(nvar,nlon,nlat)
    ! Default missing value.
    REAL, INTENT(IN) :: rwrong_val
    ! Output file name.
    CHARACTER*300, INTENT(IN) :: ncfile
    ! dim_name has 2 dimensions: lon and lat.
    CHARACTER*30, INTENT(IN) :: var_name(nvar), units_name(nvar), dim_name(2)
    ! Horizontal dimensions values.
    REAL, INTENT(IN) :: rlon(nlon), rlat(nlat)
!!!!!!!!!! Indexes for the NetCDF file manipulation.
    INTEGER il_err
    ! File and dimensions ID.
    INTEGER il_ncfile, il_dlat, il_dlon
    INTEGER il_flat, il_flon
    ! One index for each variable...
    INTEGER il_nc_field(nvar)
    INTEGER ndim
    PARAMETER(NDIM=2)
    ! ... and for each dimension.
    INTEGER ila_ncdims(ndim) 
    INTEGER ila_ncstart(ndim), ila_nccount(ndim)
!!!!!!!!!!
    INTEGER ivar
    !** OPEN NETCDF FILE
    !
    il_err = NF_CREATE(TRIM(ncfile), NF_WRITE, il_ncfile)
    !
    !** DEFINE DIMENSIONS
    !
    il_err = NF_DEF_DIM(il_ncfile, TRIM(dim_name(2)), nlat, il_dlat)
    il_err = NF_DEF_DIM(il_ncfile, TRIM(dim_name(1)), nlon, il_dlon)
    !
    !** DEFINE DIMENSION VARIABLES
    ! 
    il_err = NF_DEF_VAR(il_ncfile, TRIM(dim_name(2)), NF_REAL, &
         1, il_dlat , il_flat)
    il_err = NF_DEF_VAR(il_ncfile, TRIM(dim_name(1)), NF_REAL, &
         1, il_dlon , il_flon)
    !
    !**  DEFINE ATTRIBUTES
    !
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_flon, 'standard_name', &
         LEN('Longitude'), 'Longitude')
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_flon, 'units',         &
         LEN('degrees_east'), 'degrees_east')
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_flat, 'standard_name', &
         LEN('Latitude'), 'Latitude')
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_flat, 'units',         &
         LEN('degrees_north'), 'degrees_north')
    !
    !** DEFINE FIELD VARIABLE
    !
    ! First setting the array of dimensions indexes. 
    ila_ncdims(1) = il_dlon
    ila_ncdims(2) = il_dlat
    ! Then we define the fields for each variable.
    DO ivar=1,nvar
       il_err = NF_DEF_VAR(il_ncfile, TRIM(var_name(ivar)),                  &
            NF_REAL, ndim, ila_ncdims, il_nc_field(ivar))
       il_err = NF_PUT_ATT_TEXT(il_ncfile, il_nc_field(ivar), 'units',       &
            LEN_TRIM(units_name(ivar)), TRIM(units_name(ivar)))
       il_err = NF_PUT_ATT_real(il_ncfile, il_nc_field(ivar),'missing_value',&
            NF_REAL, 1, rwrong_val)
       il_err = NF_PUT_ATT_real(il_ncfile, il_nc_field(ivar),'_FillValue',   &
            NF_REAL, 1, rwrong_val)
    END DO

    !** END OF DEFINE MODE
    !
    il_err = NF_ENDDEF(il_ncfile)
    !
    !** WRITE DATA TO FILE
    !
    ! First writing the dimension arrays.
    ! Horizontal grid.
    il_err = NF_PUT_VAR(il_ncfile,il_flon,rlon)
    il_err = NF_PUT_VAR(il_ncfile,il_flat,rlat)
    ! Setting the first coordinate number for each dimension.
    ila_ncstart(1) = 1
    ila_ncstart(2) = 1
    ! Setting the size of each dimension axis.
    ila_nccount(1) = nlon
    ila_nccount(2) = nlat
    ! Writing the field.
    DO ivar=1,nvar
       il_err = NF_PUT_VARA(il_ncfile,il_nc_field(ivar), &
            ila_ncstart, ila_nccount, FIELD(ivar,:,:))
    END DO
    !
    !** CLOSE THE NETCDF FILE
    !
    il_err = NF_CLOSE(il_ncfile)
    !
    WRITE(*,*) 'Tonetcdf_climato: file available as ', TRIM(ncfile)
    RETURN
  END SUBROUTINE TONETCDF_CLIMATO

  ! This routine prints an input data field into a NetCDF file.
  ! It is called for the whole time period. Thus, this one it does not make time series.

  SUBROUTINE TONETCDF_ZONAL(field, nlat, dim_name, &
       nvar, var_name, units_name,                 &
       rlon1, rlon2, rlat, rwrong_val, ncfile)
    USE GLOBAL_VAR_INTERPOL_MOD, ONLY : tmpdir
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'

!!!!!!!! Input parameters !!!!!!!!
    ! 3-dimension sizes and amount of variables to write.
    INTEGER, INTENT(IN) :: nlat, nvar
    ! Data array to write.
    REAL, INTENT(IN) :: field(nvar,nlat)
    ! Default missing value.
    REAL, INTENT(IN) :: rwrong_val
    ! Output file name.
    CHARACTER*300, INTENT(IN) :: ncfile
    ! dim_name has 1 dimension: lat.
    CHARACTER*30, INTENT(IN) :: var_name(nvar), units_name(nvar), dim_name
    ! Horizontal dimensions values.
    REAL, INTENT(IN) :: rlat(nlat)
    ! Western and eastern boundaries for the current zonal band.
    REAL, INTENT(IN) :: rlon1, rlon2
    CHARACTER*30 clon1, clon2, c_boundaries
!!!!!!!!!! Indexes for the NetCDF file manipulation.
    INTEGER il_err
    ! File and dimensions ID.
    INTEGER il_ncfile, il_dlat
    INTEGER il_flat
    ! One index for each variable...
    INTEGER il_nc_field(nvar)
    INTEGER ndim
    PARAMETER(NDIM=1)
    ! ... and for each dimension.
    INTEGER ila_ncdims(ndim) 
    INTEGER ila_ncstart(ndim), ila_nccount(ndim)
!!!!!!!!!!
    INTEGER ivar, ichar
    ! Preparing the global-attribute longitudes interval, as a character string.
    CALL REAL2CHAR(rlon1,tmpdir,clon1)
    CALL REAL2CHAR(rlon2,tmpdir,clon2)
    c_boundaries=TRIM(clon1)//' '//TRIM(clon2)//' degrees east'
    DO ichar=LEN_TRIM(c_boundaries),LEN(c_boundaries)
       c_boundaries=c_boundaries//' '
    END DO
    !
    !** OPEN NETCDF FILE
    !
    il_err = NF_CREATE(TRIM(ncfile), NF_WRITE, il_ncfile)
    !
    !** DEFINE DIMENSIONS
    !
    il_err = NF_DEF_DIM(il_ncfile, TRIM(dim_name), nlat, il_dlat)
    !
    !** DEFINE DIMENSION VARIABLES
    ! 
    il_err = NF_DEF_VAR(il_ncfile, TRIM(dim_name), NF_REAL, &
         1, il_dlat , il_flat)
    !
    !**  DEFINE ATTRIBUTES
    !
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_flat, 'standard_name', &
         LEN('Latitude'), 'Latitude')
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_flat, 'units',         &
         LEN('degrees_north'), 'degrees_north')
    !
    !** DEFINE FIELD VARIABLE
    !
    ! First setting the array of dimensions indexes. 
    ila_ncdims = il_dlat
    ! Then we define the fields for each variable.
    DO ivar=1,nvar
       il_err = NF_DEF_VAR(il_ncfile, TRIM(var_name(ivar)),                  &
            NF_REAL, ndim, ila_ncdims, il_nc_field(ivar))
       il_err = NF_PUT_ATT_TEXT(il_ncfile, il_nc_field(ivar), 'units',       &
            LEN_TRIM(units_name(ivar)), TRIM(units_name(ivar)))
       il_err = NF_PUT_ATT_TEXT(il_ncfile, il_nc_field(ivar), 'longitude_domain', &
            LEN_TRIM(c_boundaries), TRIM(c_boundaries))
       il_err = NF_PUT_ATT_real(il_ncfile, il_nc_field(ivar),'missing_value',&
            NF_REAL, 1, rwrong_val)
       il_err = NF_PUT_ATT_real(il_ncfile, il_nc_field(ivar),'_FillValue',   &
            NF_REAL, 1, rwrong_val)
    END DO
    !** END OF DEFINE MODE
    !
    il_err = NF_ENDDEF(il_ncfile)
    !
    !** WRITE DATA TO FILE
    !
    ! First writing the dimension arrays.
    ! Meridional grid.
    il_err = NF_PUT_VAR(il_ncfile,il_flat,rlat)
    ! Setting the first coordinate number for each dimension.
    ila_ncstart = 1
    ! Setting the size of each dimension axis.
    ila_nccount = nlat
    ! Writing the field.
    DO ivar=1,nvar
       il_err = NF_PUT_VARA(il_ncfile,il_nc_field(ivar), &
            ila_ncstart, ila_nccount, FIELD(ivar,:))
    END DO
    !
    !** CLOSE THE NETCDF FILE
    !
    il_err = NF_CLOSE(il_ncfile)
    !
    WRITE(*,*) 'Tonetcdf_zonal: file available as ', TRIM(ncfile)
    RETURN
  END SUBROUTINE TONETCDF_ZONAL

  ! This subroutine writes monthly time series into NetCDF files.

  SUBROUTINE TONETCDF_TIME_SERIES(FIELD, nseasons, nmonths, nvar, &
       var_name, units_name, list_months, rwrong_val, tmpdir, ncfile)
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
!!!!! Declaring the subroutine input arguments.
    ! Dimension sizes: the variable, the season, and the month.
    INTEGER, INTENT(IN) :: nvar, nseasons, nmonths
    REAL,    INTENT(IN) :: FIELD(nvar, nseasons, nmonths)
    ! List of variables, and their corresponding units.
    CHARACTER*30, INTENT(IN) :: var_name(nvar), units_name(nvar)
    ! List of months in the time series to build.
    INTEGER, INTENT(IN) :: list_months(nmonths)
    ! Default missing value.
    REAL,    INTENT(IN) :: rwrong_val
    ! Temporary directory.
    CHARACTER*200,INTENT(IN) :: tmpdir
    ! Output file.
    CHARACTER*300,INTENT(IN) :: ncfile
!!!!! Other declarations.
    ! Year as character.
    CHARACTER*30 cyyyy(nmonths)
    ! Month as character.
    CHARACTER*30 cmm(nmonths)
    ! Year and month as integers, one value for each month composing the time series.
    INTEGER iyyyy(nmonths), imm(nmonths)
    ! Index running the monthly time series.
    INTEGER imth
    ! Value for time origin.
    CHARACTER*48 cltori
    ! Values for the time series axis.
    CHARACTER*48 cldat(nmonths)
!!!!!!!!!! Indexes for the NetCDF file manipulation.
    INTEGER il_err
    ! File and dimensions ID.
    INTEGER il_ncfile, il_dim_mth, il_dim_seas
    INTEGER il_fdim_mth, il_fdim_seas
    ! One index for each variable...
    INTEGER il_nc_field(nvar)
    INTEGER,PARAMETER :: ndim=2
    ! ... and for each dimension.
    INTEGER ila_ncdims(ndim) 
    INTEGER ila_ncstart(ndim), ila_nccount(ndim)
!!!!!!!!!!
    INTEGER ivar
    ! Indexes for seasons during one year.
    INTEGER iseasons(nseasons)
    iseasons=(/1,2,3,4,5/)
    !
    !** OPEN NETCDF FILE
    !
    il_err = NF_CREATE(TRIM(ncfile), NF_WRITE, il_ncfile)
    !
    !** DEFINE DIMENSIONS
    !
    il_err = NF_DEF_DIM(il_ncfile, 'season', nseasons, il_dim_seas)
    il_err = NF_DEF_DIM(il_ncfile, 'time', NF_UNLIMITED, il_dim_mth) ! Record dimension.
    ! Reading the year and the month, for each month composing the time series.
    DO imth=1,nmonths
       iyyyy(imth)=NINT(REAL(list_months(imth))/100)
       imm(imth)=MODULO(list_months(imth),100)
    END DO
    ! Conversion into characters for their first element.
    CALL INT2CHAR(iyyyy(1),tmpdir,cyyyy(1))
    CALL INT2CHAR(imm(1),tmpdir,cmm(1))
    ! Writing the time origin.
    WRITE(cltori,FMT= &
         '(A4,"-",A2,"-",A2," ",A2,":",A2,":"A2)') &
         TRIM(cyyyy(1)), TRIM(cmm(1)), '01', '00', '00', '00'
    !
    !** DEFINE DIMENSION VARIABLES
    ! 
    ! The 1 stands for the amount of dimension defining the current variable.
    ! It is evident it has to be 1 if the variable is a dimension.
    il_err = NF_DEF_VAR(il_ncfile, 'time', NF_INT, 1, il_dim_mth , &
         il_fdim_mth) 
    IF (nseasons .NE. 1) THEN
       il_err = NF_DEF_VAR(il_ncfile, 'season', NF_INT, 1, il_dim_seas, &
            il_fdim_seas) 
    END IF
    !
    !**  DEFINE ATTRIBUTES
    !
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_fdim_mth, 'standard_name', &
         LEN('Time_axis'), 'Time_axis')
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_fdim_mth, 'axis', 1, 't')
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_fdim_mth, 'time_origin', &
         LEN_TRIM(cltori), cltori)
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_fdim_mth, 'calendar', &
         LEN('noleap'), 'noleap')
    ! Explanation for the seasons dimension.
    IF (nseasons .NE. 1) THEN
       il_err = NF_PUT_ATT_TEXT(il_ncfile, il_fdim_seas, 'order', &
            LEN('DJF MAM JJA SON ANN'), 'DJF MAM JJA SON ANN')
    END IF
    !
    !** DEFINE FIELD VARIABLE
    !
    ! First setting the array of dimensions indexes.
    ila_ncdims(1) = il_dim_seas ! seasons
    ila_ncdims(2) = il_dim_mth  ! months
    ! Then we define the fields for each variable.
    DO ivar=1,nvar
       il_err = NF_DEF_VAR(il_ncfile, TRIM(var_name(ivar)),                   &
            NF_REAL, ndim, ila_ncdims(:), il_nc_field(ivar))
       il_err = NF_PUT_ATT_TEXT(il_ncfile, il_nc_field(ivar), 'units',        &
            LEN_TRIM(units_name(ivar)), TRIM(units_name(ivar)))
       ! The 1 stands for the amount of values to introduce as attributes.
       il_err = NF_PUT_ATT_REAL(il_ncfile, il_nc_field(ivar), 'missing_value',&
            NF_REAL, 1, rwrong_val)
       il_err = NF_PUT_ATT_REAL(il_ncfile, il_nc_field(ivar), '_FillValue',   &
            NF_REAL, 1, rwrong_val)
    END DO
    !** END OF DEFINE MODE
    !
    il_err= NF_ENDDEF(il_ncfile)
    !
    !** WRITE DATA TO FILE
    !
    ! First writing the dimension arrays.
    ! Setting the first coordinate number for each dimension.
    ila_ncstart(1) = 1
    ila_ncstart(2) = 1
    ! Setting the length for each dimension axis.
    ila_nccount(1) = nseasons
    ila_nccount(2) = nmonths
    ! Writing the dimensions.
    il_err = NF_PUT_VARA(il_ncfile, il_fdim_seas, ila_ncstart(1), ila_nccount(1), iseasons(:))
    il_err = NF_PUT_VARA(il_ncfile, il_fdim_mth , ila_ncstart(2), ila_nccount(2), list_months(:))
    ! Writing the field.
    DO ivar=1,nvar
       il_err = NF_PUT_VARA(il_ncfile, il_nc_field(ivar), &
            ila_ncstart(:), ila_nccount(:), FIELD(ivar,:,:))
    END DO
    !
    !** CLOSE THE NETCDF FILE
    !
    il_err = NF_CLOSE(il_ncfile)
    !
    WRITE(*,*) 'tonetcdf_time_series: file available as ', TRIM(ncfile)
    RETURN
  END SUBROUTINE TONETCDF_TIME_SERIES

  SUBROUTINE TONETCDF_1D(FIELD, ntime, nvar, var_name, units_name, &
       packName, rtime, ncfile)
    IMPLICIT NONE
    INCLUDE 'netcdf.inc'
!!!!! Declaring the subroutine input arguments.
    ! Dimension sizes: the variable, and time.
    INTEGER, INTENT(IN) :: nvar, ntime
    REAL,    INTENT(IN) :: FIELD(nvar, ntime)
    ! List of variables, and their corresponding units.
    CHARACTER*30, INTENT(IN) :: var_name(nvar), units_name(nvar)
    CHARACTER*30, INTENT(IN) :: packName
    INTEGER, INTENT(IN) :: rtime(ntime)
    ! Output file.
    CHARACTER*300,INTENT(IN) :: ncfile
!!!!!!!!!! Indexes for the NetCDF file manipulation.
    INTEGER il_err
    ! File and dimension ID.
    INTEGER il_ncfile, il_dim_time
    INTEGER il_fdim_time
    ! One index for each variable...
    INTEGER il_nc_field(nvar)
    INTEGER,PARAMETER :: ndim=1
    INTEGER ila_ncdims(ndim)
    INTEGER ila_ncstart(ndim),ila_nccount(ndim)
    INTEGER ivar
    REAL rl_daynum
    !
    !** OPEN NETCDF FILE
    !
    il_err = NF_CREATE(TRIM(ncfile), NF_WRITE, il_ncfile)
    !
    !** DEFINE DIMENSIONS
    !
    il_err = NF_DEF_DIM(il_ncfile, 'UTC_time', ntime, il_dim_time)
    !
    !** DEFINE DIMENSION VARIABLES
    ! 
    il_err = NF_DEF_VAR(il_ncfile, 'UTC_time', NF_INT, 1, il_dim_time, &
         il_fdim_time)
    !
    !**  DEFINE ATTRIBUTES
    !
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_fdim_time, 'standard_name', &
         8, 'UTC_time')
    il_err = NF_PUT_ATT_TEXT(il_ncfile, il_fdim_time, 'units', 24, &
         'UTC seconds into the day')
    il_err = NF_PUT_ATT_TEXT(il_ncfile, NF_GLOBAL, 'source', LEN_TRIM(packName), &
         TRIM(packName))
    !
    !** DEFINE FIELD VARIABLE
    !
    ila_ncdims(1) = il_dim_time
    DO ivar=1,nvar
       il_err = NF_DEF_VAR(il_ncfile, TRIM(var_name(ivar)), NF_REAL, ndim, &
            ila_ncdims, il_nc_field(ivar))
       il_err = NF_PUT_ATT_TEXT(il_ncfile, il_nc_field(ivar), 'units', &
            LEN_TRIM(units_name(ivar)), TRIM(units_name(ivar)))
    END DO
    !** END OF DEFINE MODE
    !
    il_err= NF_ENDDEF(il_ncfile)
    !
    !** WRITE DATA TO FILE
    !
    il_err = NF_PUT_VAR(il_ncfile,il_fdim_time,rtime)
    ila_ncstart(1) = 1
    ila_nccount(1) = 1
    rl_daynum=0.
    ila_ncstart(1) = 1
    ila_nccount(1) = ntime
    DO ivar=1,nvar
       il_err = NF_PUT_VARA(il_ncfile, il_nc_field(ivar), &
            ila_ncstart, ila_nccount, FIELD(ivar,:))
    END DO
    !
    !** CLOSE THE NETCDF FILE
    !
    il_err = NF_CLOSE(il_ncfile)
    WRITE(*,*) 'New IAGOS file available as ',ncfile
    RETURN
  END SUBROUTINE TONETCDF_1D
  
  !     This subroutine aims at writing the data read from IAGOS T-files
  !     into a CDL file. The input has to be an array of strings of characters.

  !     Unfortunately, it doesn't write when the input data are issued from
  !     the subroutine REAL2CHAR. Everything happens normally,
  !     except the write statement only.
  !
  !     ---> Resolved. I just have to take care not to use the same unitFile
  !     than an open file.

  SUBROUTINE WRITE_DATA_CDL(vec_data,n_data,unitFile,nreturn,nchar)
    IMPLICIT NONE

    INTEGER n_data,unitFile,ncol,nchar,nmar
    INTEGER nreturn,n_last_data
    INTEGER idata
    CHARACTER*100 outfile,FMTvar
    CHARACTER*20 vec_data(n_data)
    PARAMETER(nmar=2)
    !     Writing time with 10 elements per line
    !     For simplicity, returns line before the last element
    idata=1
    DO WHILE (idata.LE.(n_data-nreturn))
       WRITE(FMTvar,'( "(", i0, "a", i0,")" )') nreturn, nchar+nmar
       !     nchar+nmar just to create a marge (size nmar) between the data
       WRITE(unitFile,FMTvar) &
            vec_data(idata:(idata+nreturn-1))//','
       idata=idata+nreturn
    END DO
    IF (idata.NE.n_data) THEN
       WRITE(FMTvar,'( "(", i0, "a", i0,")" )') &
            n_data-idata,nchar+nmar
       WRITE(unitFile,FMTvar) &
            vec_data(idata:(n_data-1))//','
    END IF
    WRITE(FMTvar,'( "(a", i0,")" )') &
         nchar+nmar
    WRITE(unitFile,FMTvar) &
         vec_data(n_data)//';'
    WRITE(*,*) '-------------------------'   
    WRITE(*,*) 'vec_data=',vec_data(n_data)
    IF (n_data .LT. 200) THEN ! Check option
       WRITE(*,*) 'WRITE_DATA_CDL says:',vec_data(1:n_data)
    ELSE
       WRITE(*,*) 'WRITE_DATA_CDL says:',vec_data(1:200)
    END IF
    WRITE(*,*) 'length(vec_data)=',n_data
  END SUBROUTINE WRITE_DATA_CDL
  
  ! This routine aims at deriving the 3D hybrid sigma-pressure fields,
  ! based on a surface pressure field and two A and B coefficient arrays.
  SUBROUTINE COMPUTE_PRESSURES(nlon, nlat, nlev, ndays, P_surf, P_0, &
       A_inter, B_inter, P_half, P_full)
    IMPLICIT NONE
!!!!! Arguments declaration.
    ! Dimensions size.
    INTEGER, INTENT(IN) :: nlon, nlat, nlev, ndays
    ! Surface pressure field, and reference pressure scalar.
    REAL, INTENT(IN) :: P_surf(nlon, nlat, ndays), P_0
    ! Vertical grid coefficients.
    REAL, INTENT(IN) :: A_inter(nlev), B_inter(nlev)
    ! 3D half-level pressures.
    REAL, INTENT(OUT) :: P_half(nlon, nlat, nlev, ndays)
    ! 3D full-level pressures.
    REAL, INTENT(OUT) :: P_full(nlon, nlat, nlev, ndays)
!!!!! Other variables.
    ! Indexes running the model grid.
    INTEGER i, k, l, iday
    ! Pressures at half-levels l and l-1 respectively.
    REAL P_l, P_lm1

    DO iday=1,ndays
       DO l=1,nlev
          DO k=1,nlat
             DO i=1,nlon
                P_half(i,k,l,iday)= &
                     A_inter(l)*P_0 + &
                     B_inter(l)*P_surf(i,k,iday)
             END DO
          END DO
       END DO
       ! Definition of the full vertical grid levels.
       ! First: definition of the first level
       DO k=1,nlat
          DO i=1,nlon     
             P_full(i,k,1,iday)= &
                  EXP(LOG(P_half(i,k,1,iday))-1)
          END DO
       END DO
       ! Second: definition of the lower levels (depending on l-1)
       DO l=2,nlev
          DO k=1,nlat
             DO i=1,nlon
                P_l = P_half(i,k,l,iday)
                P_lm1 = P_half(i,k,l-1,iday)
                P_full(i,k,l,iday)= EXP( &
                     (P_l*LOG(P_l)-P_lm1*LOG(P_lm1)) &
                     /(P_l-P_lm1) -1.)
             END DO
          END DO
       END DO
    END DO ! end iday=1,ndays
  END SUBROUTINE COMPUTE_PRESSURES
  
  ! lat_ratio computes an array of factors
  ! that depend on latitude, because of the
  ! size of the grid cells that decreases with latitude.
  SUBROUTINE LATITUDE_SURFACE_RATIO(nlat, vec_lat, vec_ratio)
    USE GLOBAL_VAR_INTERPOL_MOD, ONLY : Pi
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nlat
    REAL,    INTENT(IN) :: vec_lat(nlat)
    REAL,    INTENT(OUT):: vec_ratio(nlat)
    REAL cos_lat, rad_lat
    REAL mean_cos_lat
    INTEGER ilat
    ! 1/ Computing the denominator, as the meridional average of the latitude cosine.
    mean_cos_lat = 0
    DO ilat=1,nlat
       rad_lat      = vec_lat(ilat)*Pi/180
       mean_cos_lat = mean_cos_lat+COS(rad_lat)
    END DO
    mean_cos_lat = mean_cos_lat/REAL(nlat)
    ! 2/ Computing the final vector, made of one factor for each latitude.
    DO ilat=1,nlat
       rad_lat         = vec_lat(ilat)*Pi/180
       cos_lat         = COS(rad_lat)
       vec_ratio(ilat) = cos_lat/mean_cos_lat
    END DO
  END SUBROUTINE LATITUDE_SURFACE_RATIO

!****s* Attributes/nf90_get_att_string
!
! NAME
!    nf90_get_att_string - Read a string attribute from a netCDF data file.
!
! SYNOPSIS
!    use ncdf
!      ...
!    call nf90_get_att_string(ncid, varid, attname, value) 
!
! DESCRIPTION
!    This subroutine reads string attribute data from the current netCDF file.
!    It extends the functionality of nf90_get_att function by adding a 
!    possibility to read single string attributes (NF90_STRING) which at the
!    moment netcdf-fortran library version 4.4.5 can not do.
!
! INPUTS
!
! OUTPUT
!    value
!
! AUTHOR
!    Leonid Butenko, Darmstadt, Germany            <leonid.butenko@eumetsat.int>
!
! COPYRIGHT
!
!    Copyright (c) 2019 EUMETSAT
!
!    All rights reserved.
!
!    Permission is hereby granted, free of charge, to any person obtaining
!    a copy of this software and associated documentation files (the
!    "Software"), to deal in the Software without restriction, including
!    without limitation the rights to use, copy, modify, merge, publish,
!    distribute, sublicense, and/or sell copies of the Software, and to
!    permit persons to whom the Software is furnished to do so, subject to
!    the following conditions:
!
!    The above copyright notice and this permission notice shall be
!    included in all copies or substantial portions of the Software as well
!    as in supporting documentation.
!
!    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
!    EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
!    MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
!    NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
!    LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
!    OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
!    WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
!    Available at https://trac.romsaf.org/ropp/
!****


FUNCTION NF90_GET_ATT_STRING(ncid, varid, attname, value) RESULT(status)

  USE typeSizes
  USE NCDF, NOT_THIS => NF90_GET_ATT_STRING
  USE NETCDF, ONLY: NF90_INQUIRE_ATTRIBUTE
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_size_t, c_f_pointer, c_int
  IMPLICIT NONE

  INTERFACE
    FUNCTION NC_GET_ATT_STRING(ncid, varid, name, ip) BIND(c)
      USE ISO_C_BINDING, ONLY: c_int, c_char, c_ptr

      INTEGER(c_int),         VALUE       :: ncid, varid
      CHARACTER(KIND=c_char), INTENT(IN)  :: name
      TYPE(c_ptr), INTENT(OUT)            :: ip

      INTEGER(c_int     )                 :: nc_get_att_string

    END FUNCTION NC_GET_ATT_STRING
  END INTERFACE

  INTERFACE
    FUNCTION STRLEN(s) BIND(c, NAME='strlen')
       USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_ptr, c_size_t
       IMPLICIT NONE
       !----
       TYPE(c_ptr), INTENT(IN), VALUE :: s
       INTEGER(c_size_t) :: strlen
    END FUNCTION STRLEN
  END INTERFACE

  INTEGER,                  INTENT(IN ) :: ncid
  INTEGER,                  INTENT(IN ) :: varid
  CHARACTER(LEN = *),       INTENT(IN ) :: attname
  CHARACTER(LEN = *),       INTENT(OUT) :: value

  INTEGER                                :: status, xtype, nlen, attid, i

  INTEGER(c_int)                         :: c_ncid, c_varid, c_status, c_nlen
  TYPE(c_ptr)                            :: c_str
  CHARACTER(LEN=:), ALLOCATABLE          :: c_aname
  CHARACTER, POINTER                     :: f_str(:)

  status = nf90_inquire_attribute( ncid, varid, attname, xtype, nlen, attid )
  IF (status /= nf90_noerr) WRITE(*,*) NF90_STRERROR(status) !call ncdf_error_handler(status)

  IF (xtype == NF90_STRING .AND. nlen == 1) THEN

     c_ncid = ncid
     c_varid = varid-1 ! C-library counts variables starting with 0
     ALLOCATE( CHARACTER(LEN_TRIM(attname)+1) :: c_aname )
     c_aname = TRIM(attname)//CHAR(0)
     value = ADJUSTL("")

     c_status = nc_get_att_string( c_ncid, c_varid, c_aname, c_str )
     status = c_status

     CALL C_F_POINTER( c_str, f_str,  [STRLEN(c_str)] )
     
     ! convert char array to string
     DO i=1,SIZE(f_str)
        value(i:i) = f_str(i)
     END DO
  ELSE !! all others
    status = NF90_GET_ATT( ncid, varid, attname, value )
  ENDIF

END FUNCTION NF90_GET_ATT_STRING

END MODULE FUNCTIONS_MOD
