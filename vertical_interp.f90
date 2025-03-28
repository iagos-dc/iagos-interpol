SUBROUTINE VERTICAL_INTERP(nlev_in, nlev_out, P_in, v_in, P_out, v_out, version_vertical_interp)
  USE GLOBAL_VAR_INTERPOL_MOD
  USE FUNCTIONS_MOD
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: nlev_in, nlev_out, version_vertical_interp
  REAL,    INTENT(IN) :: P_in(nlev_in), P_out(nlev_out)
  REAL,    INTENT(IN) :: v_in(nlev_in)
  REAL,    INTENT(OUT):: v_out(nlev_out)
  REAL,    ALLOCATABLE:: P_out_up(:), P_in_up(:) ! Pressure at the upper limit of a gridcell.
  INTEGER, ALLOCATABLE:: sampled_lev(:)
  INTEGER k, j, jp, llev_out_inf, l_up, l_down, i_count
  
  ! P_edge is the pressure at the edge of an output grid cell.
  ! It is calculated as the middle between two P levels, in altitude.
  REAL alpha_j_k, log_P_in_ratio, log_P_ratio, wght1, wght2, wght3, wght_sum, wght_check
  
  ! First: defining the boundary pressure for the input gridcells and output gridcells too,
  ! assuming it is equidistant between two full levels. We use the full pressure then.
  ALLOCATE( P_out_up(nlev_out + 1), P_in_up(nlev_in + 1) )
  DO k = 2, nlev_out
     P_out_up(k) = (P_out(k) * P_out(k-1))**0.5
  END DO
  ! Now, we derive the 1 and nlev_out gridcell limits.
  P_out_up(1) = P_out(1)**2 / P_out_up(2)
  P_out_up(nlev_out + 1) = P_out(nlev_out)**2 / P_out_up(nlev_out) ! supposed to be the surface pressure.

  DO j = 2, nlev_in
     P_in_up(j) = (P_in(j) * P_in(j-1))**0.5
     ! P_in_down(j) = (P_in(j_up+1) * P_in(j_up))**0.5
  END DO
  P_in_up(1) = P_in(1)**2 / P_in_up(2)
  P_in_up(nlev_in + 1) = P_in(nlev_in)**2 / P_in_up(nlev_in) ! supposed to be the surface pressure.

  ! First: calculates the output lowest level corresponding to an input level.
  ! The point is that the output levels are constant pressure levels, contrary to the input levels.
  ! This way, we apply a mask on the output levels that cut across mountains chains.
  llev_out_inf = 0
  DO WHILE ( P_out_up(llev_out_inf + 1) < P_in_up(nlev_in + 1) .AND. llev_out_inf < nlev_out )
     llev_out_inf = llev_out_inf + 1
  END DO
  ! llev_out_inf is thus the output level that corresponds both to the atmosphere and to the ground.
  ! Below, there is no atmosphere:
  v_out( (llev_out_inf + 1) : nlev_out) = rwrong_val
  ! And at this level specifically, we have to project the input gridcells that share some altitudes
  ! with this output level. In this case exceptionally, we do not apply the mass conservation,
  ! because the volume of the output gridcell being overestimated, it would result in an underestimation
  ! of the new mixing ratio (because too much dilution).
  ! And the IAGOS instruments measure the mixing ratio, not the mass.

  ! sampled_lev is a vector that is used in case the output grid is more resolved than the input grid.
  ! Its role is to say which output level has remained empty after the interpolation, so we can fill
  ! the output grid with an interpolation between its own sampled levels.
  ALLOCATE(sampled_lev(llev_out_inf))
  sampled_lev(:) = 0
  v_out( 1 : llev_out_inf) = 0.
  IF (version_vertical_interp == 0) THEN
     DO k = 1, llev_out_inf
        ! First, we locate the boundary pressure of the output grid cell.
        ! Supposing negligible changes in temperature between the 2 grid cells.
        ! The max and min commands are only here to avoid if loops for boundary conditions.
        ! P_edge_up = (P_out(k) * P_out(max(1,k-1)))**0.5
        ! P_edge_down = (P_out(min(nlev_out,k+1)) * P_out(k))**0.5
        ! For the current output level, we memorize the first input level
        ! with a lower pressure.
        j = nlev_in
        ! Boundary condition: first we treat the case where there is no model grid level below.
        IF ( P_out(k) .GT. P_in(nlev_in) ) THEN
           v_out(k) = v_in(nlev_in)
        ELSE
           DO WHILE (P_in(j) .GE. P_out(k) .AND. j .GE. 1)
              j = j - 1
           END DO
           ! i_index(k) = j
           log_P_in_ratio = LOG(P_in(j+1)/P_in(j))
           log_P_ratio = LOG(P_out(k)/P_in(j))
           ! write(*,*) 'v_in up, v_in down =', v_in(j), v_in(j+1)
           ! Assigning weights and averaging, depending on the distance with the output level.
           alpha_j_k = ABS( log_P_ratio / log_P_in_ratio )
           v_out(k) = (1 - alpha_j_k)*v_in(j) + alpha_j_k*v_in(j+1)
        END IF
     END DO
  ELSE IF (version_vertical_interp == 1) THEN
     DO k = 1, llev_out_inf
        wght1 = 0.
        wght2 = 0.
        wght3 = 0.
        wght_sum = 0. ! This variable returns 1 if the mass is conserved.
        wght_check = 0. ! This variable returns 0 in absence of negative weights.
        ! First, we locate the boundary pressure of the output grid cell.
        ! Supposing negligible changes in temperature between the 2 grid cells.
        ! The max and min commands are only here to avoid if loops for boundary conditions.

        ! For the current output level, we memorize the first input level
        ! that superimposes the output grid cell.
        j = 1
        DO WHILE (P_in_up(j + 1) < P_out_up (k) .AND. j .NE. nlev_in)
           j = j + 1
        END DO

        ! This IF loop stems for the case when the j input gridcell
        ! is totally included into the k output gridcell. Then,
        ! we interpolate linearly using the distance between the input and output gridcell centers.
        ! If it is not the case (as we generally expect it), then we use weighting factors that
        ! depend on the intersection between the input and output gridcells,
        ! normalized by the size of the output gridcell.
        IF ( P_in_up(j + 1) > P_out_up(k + 1) ) THEN ! i.e. if the j input gridcell totally includes the k output gridcell. This case is treated after.
           IF (k == llev_out_inf .AND. P_in(j) > P_out(k)) THEN
              ! In this case, the k pressure level is higher in altitude than the j input level, so we involve the input gridcell above.
              log_P_in_ratio = LOG(P_in(j)/P_in(j-1))
              v_out(k) = v_in(j) * ( 1 - ABS(LOG(P_out(k)/P_in(j)))/ log_P_in_ratio ) + v_in(j-1) * ( 1 - ABS(LOG(P_out(k)/P_in(j-1)))/ log_P_in_ratio )
           ELSE IF (k == llev_out_inf .AND. P_in(j) <= P_out(k)) THEN
              ! In this case, the k pressure level is lower in altitude than the j input level, so we involve the input gridcell below.
              log_P_in_ratio = LOG(P_in(j+1)/P_in(j))
              v_out(k) = v_in(j) * ( 1 - ABS(LOG(P_out(k)/P_in(j)))/log_P_in_ratio ) + v_in(j+1) * ( 1 - ABS(LOG(P_out(k)/P_in(j+1)))/log_P_in_ratio )
           ELSE
              v_out(k) = 0.
           END IF
        ELSE
           ! This MIN option should take into account the cases where the input gridcell is bigger or smaller than the output gridcell.
           wght1 = LOG(P_out_up(k)/P_in_up(j + 1)) / LOG(P_out_up(k)/P_out_up(k + 1))
           v_out(k) = v_out(k) + v_in(j)*wght1
           wght_sum = wght_sum + wght1
           wght_check = wght_check + MIN(wght1, 0.)
           ! Secondary incrementation, in case the output gridcell is larger than the input gridcell (and thus needs several input gridcells to be calculated).
           jp = 1
           DO WHILE ( P_in_up(j + jp + 1) < P_out_up (k + 1) .AND. (j + jp) < nlev_in )
              wght2 = LOG( P_in_up(j+jp)/P_in_up(j+jp+1) ) / LOG( P_out_up(k)/P_out_up(k+1) )
              wght_sum = wght_sum + wght2
              wght_check = wght_check + MIN(wght2, 0.)
              v_out(k) = v_out(k) + v_in(j+jp)*wght2
              jp = jp + 1
           END DO
           ! Third: the lowest level j+jp, which lowest part contributes to the k gridcell in the output.
           ! There is also the possibility that the lowest input level remains totally included
           ! in the k output level, when the orography is elevated. Hence the MIN option below:
           ! in the normal case, the input volume is the interval between the pressures
           ! P_in_up(j + jp) and P_out_up(k + 1).
           ! In the high orography case, the input volume is the input gridcell's volume, thus:
           ! the interval between the pressures P_in_up(j + jp) and P_in_up(j + jp + 1).
           IF (j < nlev_in) THEN
              wght3 = LOG( P_in_up(j + jp)/ MIN(P_in_up(j + jp + 1), P_out_up(k + 1)) ) / &
                   LOG( P_out_up(k)/P_out_up(k+1) )
              wght_sum = wght_sum + wght3
              wght_check = wght_check + MIN(wght3, 0.)
              v_out(k) = v_out(k) + v_in(j+jp)*wght3
           END IF
           IF (j + jp <= nlev_in) THEN
              IF (P_in_up(j + jp + 1) < P_out_up(k + 1)) THEN
                 v_out(k) = v_out(k) * LOG( P_out_up(k)/P_out_up(k+1) ) / &
                      LOG( P_out_up(k)/P_in_up(j + jp + 1) )
              END IF
           ELSE IF (j + jp == nlev_in + 1) THEN
              IF (P_in_up(j + 1) < P_out_up(k + 1)) THEN
                 v_out(k) = v_out(k) * LOG( P_out_up(k)/P_out_up(k+1) ) / &
                      LOG( P_out_up(k)/P_in_up(j + 1) )
              END IF
           END IF
           sampled_lev(k) = 1
        END IF
        ! Last: if the output (constant) pressure level is never sampled by the input pressure levels.
        IF (P_in_up(j + jp + 1) < P_out_up(k)) THEN
           v_out(k) = rwrong_val
        END IF
        IF (wght_check .NE. 0) THEN
           WRITE(*,*) 'VERTICAL_INTERP: Negative weight somewhere at k, j, jp, wght_check=', k, j, jp, wght_check
           WRITE(*,*) 'P_in_up = ', P_in_up(:)
           WRITE(*,*) 'P_out_up = ', P_out_up(:)
           STOP
        END IF
     END DO ! k = 1, llev_out_inf

     ! Final step: the self completion, in case some output gridcells have been left empty.
     ! For each empty interval, the indexes l_up and l_down represent the top and the lowest empty output gridcells.

     l_up = 1 ; l_down = 1
     i_count = 0
     DO WHILE ( SUM(sampled_lev(l_down:llev_out_inf)) .NE. (llev_out_inf - l_down + 1) .AND. l_up < llev_out_inf .AND. i_count < 200) ! We launch the loop as long as there remain zeros in the sampled_lev array below l_down.
        i_count = i_count + 1

        DO WHILE (sampled_lev(l_up) == 1 .AND. l_up < llev_out_inf)
           l_up = l_up + 1
        END DO
        l_down = l_up
        DO WHILE (sampled_lev(l_down) == 0 .AND. l_down < llev_out_inf)
           l_down = l_down + 1
        END DO
        l_down = l_down - 1
        DO k = l_up, l_down
           log_P_ratio = LOG( P_out(l_down + 1)/P_out(l_up - 1) )
           v_out(k) = v_out(l_up - 1)* ABS(LOG(P_out(k)/P_out(l_down+1))) + v_out(l_down + 1)* ABS(LOG(P_out(k)/P_out(l_up-1)))
           v_out(k) = v_out(k) / log_P_ratio
        END DO
        l_down = l_down + 1
        l_up = l_down
     END DO
     IF (i_count >= 200) THEN
        WRITE(*,*) "VERTICAL_INTERP: the last while loop seems to have no end."
        WRITE(*,*) "i_count=", i_count
        STOP
     END IF
  END IF ! (version_vertical_interp == 0)
END SUBROUTINE VERTICAL_INTERP
