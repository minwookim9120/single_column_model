MODULE Mod_dyn_driver

  USE Mod_global
  IMPLICIT NONE

  CONTAINS

  SUBROUTINE Sub_dyn_driver
    IMPLICIT NONE
    INTEGER                     :: idbin
    IF (dyn_option .eq. 1) THEN
      CALL Sub_Finite_diff             &
           (                           &
            temp%dz, temp%sfc_dt(it),  &
            temp%top_dt(it),           &
            dz%dz, nz,                 &
            dt,                        &
            w%stag_dz,                 &
            temp%next_dz               &
           )

      CALL Sub_Finite_diff             &
           (                           &
            q%dz, q%sfc_dt(it),        &
            q%top_dt(it),              &
            dz%dz, nz,                 &
            dt,                        &
            w%stag_dz,                 &
            q%next_dz                  &
           )

    ELSE IF (dyn_option .eq. 2) THEN
      CALL Sub_Finite_volume           &
           (                           &
            temp%dz, temp%sfc_dt(it),  &
            temp%top_dt(it),           &
            dz%dz, nz, CFL%dz,         &
            dt,                        &
            w%stag_dz,                 &
            temp%next_dz               &
           )
  
      CALL Sub_Finite_volume           &
           (                           &
            q%dz, q%sfc_dt(it),        &
            q%top_dt(it),              &
            dz%dz, nz, CFL%dz,         &
            dt,                        &
            w%stag_dz,                 &
            q%next_dz                  &
           )
    ELSE IF (dyn_option .eq. 3) THEN
      CALL Sub_Finite_volume_PPM      &
           (                          &
            temp%dz, temp%sfc_dt(it), &                     
            temp%top_dt(it),          &                     
            dz%dz, nz, CFL%dz,        &                 
            dt,                       &
            w%stag_dz,                &                     
            temp%next_dz              &                     
           )                      

      CALL Sub_Finite_volume_PPM      &  
           (                          & 
            q%dz, q%sfc_dt(it),       &    
            q%top_dt(it),             &       
            dz%dz, nz, CFL%dz,        &                      
            dt,                       &                      
            w%stag_dz,                &                      
            q%next_dz                 &                      
           )                       

      DO idbin = 1, drop_column_num
        drop%num(:,idbin) = drop%num(:,idbin)*drop%dr(idbin)
        CALL Sub_Finite_volume_PPM      &
             (                          &
              drop%num(:,idbin), 0.,    &
              0.,                       &
              dz%dz, nz, CFL%dz,        &
              dt,                       &
              w%stag_dz,                &
              drop%next_num(:,idbin)    &
             )
         drop%next_num(:,idbin) = drop%next_num(:,idbin)/drop%dr(idbin)
      ENDDO

    ENDIF
  END SUBROUTINE Sub_dyn_driver

  SUBROUTINE Sub_Finite_diff      &
             (                    &
               var, sfc_var,      &
               top_var,           &
               dz, nz,            &
               dt,                &
               w,                 &
               next_var           &
             )

    IMPLICIT NONE
    ! IN
    REAL,                       INTENT(IN)    :: dt
    INTEGER,                    INTENT(IN)    :: nz
    REAL,                       INTENT(IN)    :: sfc_var
    REAL,                       INTENT(IN)    :: top_var
    REAL,    DIMENSION(nz),     INTENT(IN)    :: var
    REAL,    DIMENSION(0:nz),   INTENT(IN)    :: w
    REAL,    DIMENSION(nz),     INTENT(IN)    :: dz
    ! Local
    INTEGER                                   :: i
    REAL,    DIMENSION(nz)                    :: c 
    REAL,    DIMENSION(nz)                    :: interp_w
    ! OUT
    REAL,    DIMENSION(nz),     INTENT(OUT)   :: next_var

    Do i = 1, nz
    interp_w(i) = (w(i-1)+w(i))/2
    ENDDO
    IF ( w(nz) .eq. 0.0 ) interp_W(nz) = 0.0

    DO i = 1, nz

      c(i) = interp_w(i)*dt/dz(i)
      IF ( c(i) .le. 1) THEN !! CFL filter
        IF ( i .eq. 1 ) THEN
          next_var(i) = var(i) - interp_w(i)*dt*(var(i+1)-sfc_var)/(2*dz(i))
        ELSE IF ( i .eq. nz ) THEN
          next_var(i) = var(i) - interp_w(i)*dt*(top_var-var(i-1))/(2*dz(i))
        ELSE
          next_var(i) = var(i) - interp_w(i)*dt*(var(i+1)-var(i-1))/(dz(i)+dz(i+1))
        ENDIF
      ELSE
        write(*,*) " i    = ", i,"CFL  =  ", c(i)
        CALL FAIL_MSG("ERROR : You need to check CFL condition.")
      ENDIF

    ENDDO
 
  END SUBROUTINE Sub_Finite_diff

  !----------------------------------------------!
  SUBROUTINE Sub_Finite_volume    &
             (                    &
               var, sfc_var,      &
               top_var,           &
               dz, nz, CFL,       &
               dt,                &
               w,                 &
               next_var           &
             )

    IMPLICIT NONE
    ! IN
    REAL,                       INTENT(IN)    :: dt
    INTEGER,                    INTENT(IN)    :: nz
    REAL,                       INTENT(IN)    :: sfc_var
    REAL,                       INTENT(IN)    :: top_var
    REAL,    DIMENSION(nz),     INTENT(IN)    :: var,     &
                                                 dz,      &
                                                 CFL
    REAL,    DIMENSION(1:nz+1), INTENT(IN)    :: w
    ! Local
    INTEGER                                   :: i
    REAL                                      :: c, rst,  &
                                                 dvar  
    REAL,    DIMENSION(nz)                    :: slp            
    REAL,    DIMENSION(nz+1)                  :: flux
    ! OUT
    REAL,    DIMENSION(nz),     INTENT(OUT)   :: next_var

    flux=0.0
    ! Do outflow B.C
    CALL Sub_cal_slope ( var, dz, nz, slp )

    ! flux(1)    = 0.
    ! flux(nz+1) = 0.

    c   = dt*w(1)/dz(1)
    rst = (sfc_var + 0.5*slp(1)*(1.-c)) 
    flux(1) = w(1)*rst

    c   = dt*w(nz+1)/dz(nz)
    rst = (top_var + 0.5*slp(nz)*(1.-c)) 
    flux(nz+1) = w(nz+1)*rst

    DO i = 1, nz+1
      IF (w(i) >= 0.) THEN
        IF (i == 1) CYCLE          ! inflow
        c   = dt*w(i)/dz(i-1)  
        rst = (var(i-1) + 0.5*slp(i-1)*(1.-c))
      ELSE
        IF (i == nz+1) CYCLE       ! inflow
        c   = -dt*w(i)/dz(i)
        rst = (var(i) - 0.5*slp(i)*(1.-c)) 
      END IF
      flux(i) = w(i) * rst
      IF (c .gt. 1.)  CALL FAIL_MSG("Courant number > 1")
    ENDDO
 
    DO i = 1, nz
      dvar        = - (flux(i+1) - flux(i)) / dz(i) 
      next_var(i) = var(i) + dvar * dt
      IF ( next_var(i) .lt. 0. ) THEN !! mass conservation filter
        CALL FAIL_MSG("ERROR :: dynamics, Physical quantity cannot be negative")
      ENDIF
    END DO

  END SUBROUTINE Sub_Finite_volume

  !----------------------------------------------!
  SUBROUTINE Sub_Finite_volume_PPM    &
             (                        &
               var, sfc_var,          &
               top_var,               &
               dz, nz, CFL,           &
               dt,                    &
               w,                     &
               next_var               &
             ) 


    INTEGER,                    INTENT(IN) :: nz
    REAL,                       INTENT(IN) :: dt
    REAL,    DIMENSION(:),      INTENT(IN) :: dz, var
    REAL,    DIMENSION(1:nz+1), INTENT(IN) :: w
    REAL,    DIMENSION(nz),     INTENT(IN) :: CFL                                   
    REAL,    DIMENSION(nz)                 :: next_var
    REAL,    DIMENSION(nz)                 :: slp,                               & 
                                              var_left, var_right
    REAL,    DIMENSION(nz+1)               :: flux
    REAL,    DIMENSION(0:3,nz)             :: zwt
    REAL                                   :: xx, a, b, rm, r6, rst, wt
    REAL                                   :: tt, c1, c2
    REAL                                   :: sfc_var, top_var 
    LOGICAL                                :: test_1
    INTEGER                                :: i, j, k, ks, ke
    REAL                                   :: cn, rsum, dzsum, dtw
    REAL                                   :: cflerr, cflmaxx, cflmax, cflmaxcc
    REAL                                   :: dvar, c
    INTEGER                                :: kk

    ke=nz
    ks=1

    Call Sub_cal_weights ( dz, zwt )
    CALL Sub_cal_slope ( var, dz, nz, slp )

    DO k = 3, nz-1
      var_left(k) = var(k-1) + zwt(1,k)*(var(k)-var(k-1)) &
                       - zwt(2,k)*slp(k) + zwt(3,k)*slp(k-1)
      var_right(k-1) = var_left(k)
    ENDDO
   
   ! boundary values  
    var_left (1) = var(1) - 0.5*slp(1)
    var_right(1) = var(1) + 0.5*slp(1)
    var_left (nz) = var(nz) - 0.5*slp(nz)
    var_right(nz) = var(nz) + 0.5*slp(nz)

    ! make linear assumption near boundary
    ! var_left (2) = var(2) - 0.5*slp(2)
    ! var_right(nz-1) = var(nz-1) + 0.5*slp(nz-1)
    var_left (2) = var_right(1) 
    var_right(nz-1) = var_left(nz) 

    IF (.true.) THEN
      ! limiters from Lin (2003), Equation 6 (relaxed constraint)
      DO k = 1, nz
        var_left (k) = var(k) - sign(min(abs(slp(k)),abs(var_left(k)-var(k))), slp(k) )
        var_right(k) = var(k) + sign(min(abs(slp(k)),abs(var_right(k)-var(k))), slp(k) )
      ENDDO
    ELSE
      ! limiters from Colella and Woodward (1984), Equation 1.10
      DO k = ks, ke
        test_1 = (var_right(k)-var(k))*(var(k)-var_left(k)) <= 0.0
        IF (test_1) THEN
          var_left(k)  = var(k)
          var_right(k) = var(k)
        ENDIF
        IF (k == ks .or. k == ke) CYCLE
          rm = var_right(k) - var_left(k)
          a = rm*(var(k) - 0.5*(var_right(k) + var_left(k)))
          b = rm*rm/6.
        IF (a >  b) var_left (k) = 3.0*var(k) - 2.0*var_right(k)
        IF (a < -b) var_right(k) = 3.0*var(k) - 2.0*var_left (k)
      ENDDO
    ENDIF

    ! compute fluxes at interfaces
    ! flux(ks)   = w(ks)  *var(ks)/dz(1)
    ! flux(ke+1) = w(ke+1)*var(ke)/dz(nz+1)

    ! B.C = 0. 
    flux(ks)   = 0. 
    flux(ke+1) = 0. 

     c   = dt*w(1)/dz(1)
     rst = (sfc_var + 0.5*slp(1)*(1.-c)) !/ dz(1)
     flux(1) = w(ks)*rst 
     rst = (top_var + 0.5*slp(nz)*(1.-c)) !/ dz(nz)
     flux(ke+1) = w(ke+1)*rst

    ! Cal. flux
    tt = 2./3.
    DO k = 1, nz+1
      IF (w(k) >= 0.) THEN
        IF (k == ks) CYCLE ! inflow
        cn = dt*w(k)/dz(k-1)
        kk = k-1
        ! extension for Courant numbers > 1
        IF (cn > 1.) THEN
          rsum = 0.
          dzsum = 0.
          dtw = dt*w(k)
          DO WHILE (dzsum+dz(kk) < dtw)
            IF (kk == 1) THEN
              exit
            ENDIF
            dzsum = dzsum + dz(kk)
             rsum =  rsum +  var(kk)
            kk = kk-1
          ENDDO
          xx = (dtw-dzsum)/dz(kk)
        ELSE
          xx = cn
        ENDIF
        rm = var_right(kk) - var_left(kk)
        r6 = 6.0*(var(kk) - 0.5*(var_right(kk) + var_left(kk)))
        IF (kk == ks) r6 = 0.
        rst = ( var_right(kk) - 0.5*xx*(rm - (1.0 - tt*xx)*r6) ) !/ dz(k-1)
         
        ! extension for Courant numbers > 1
        IF (cn > 1.) rst = (xx*rst + rsum)/cn
      ELSE
        IF (k == ke+1) CYCLE ! inflow
        cn = - dt*w(k)/dz(k)
        kk = k
        ! extension for Courant numbers > 1
        IF (cn > 1.) THEN
          rsum = 0.
          dzsum = 0.
          dtw = -dt*w(k)
          DO WHILE (dzsum+dz(kk) < dtw)
            IF (kk == ks) THEN
              EXIT
            ENDIF
            dzsum = dzsum + dz(kk)
             rsum =  rsum + var(kk)
            kk = kk+1
          ENDDO
          xx = (dtw-dzsum)/dz(kk)
        ELSE
          xx = cn
        ENDIF
        rm = var_right(kk) - var_left(kk)
        r6 = 6.0*(var(kk) - 0.5*(var_right(kk) + var_left(kk)))
        IF (kk == ke) r6 = 0.
        rst = ( var_left(kk) + 0.5*xx*(rm + (1.0 - tt*xx)*r6) ) !/ dz(k)
        ! extension for Courant numbers > 1
        IF (cn > 1.) rst = (xx*rst + rsum)/cn
      ENDIF
      flux(k) = w(k)*rst
    ENDDO

    ! Cal. FV
    DO i = 1, nz
      dvar        = - (flux(i+1) - flux(i))/dz(i)
      ! dvar        = - (flux(i+1)/dz(i+1) - flux(i)/dz(i))
      next_var(i) = var(i) + dvar * dt
      ! IF ( next_var(i) .lt. 0. ) THEN !! mass conservation filter
      !   print*, next_var(i), var(i), dvar, dt
      !   CALL FAIL_MSG("ERROR :: dynamics, Physical quantity cannot be negative")
      ! ENDIF
    END DO

  END SUBROUTINE Sub_Finite_volume_PPM


  SUBROUTINE Sub_cal_slope ( var, dz, nz, slope)

    INTEGER,               INTENT(IN)   :: nz 
    REAL,    DIMENSION(:), INTENT(IN)   :: var, dz
    REAL,    DIMENSION(:), INTENT(OUT)  :: slope

    REAL    :: grad(2:nz)
    REAL    :: cmin, cmax
    INTEGER :: k
    LOGICAL :: limiters, dolinear

    limiters = .true.
    dolinear = .false.
     
    ! compute slope (weighted for unequal levels)
    DO k = 2, nz
      grad(k) = (var(k)-var(k-1))/(dz(k)+dz(k-1))
    ENDDO
    IF (dolinear) THEN
      DO k = 2, nz-1
        slope(k) = (grad(k+1)+grad(k))*dz(k)
      ENDDO
    ELSE
      DO k = 2, nz-1
        slope(k) = (grad(k+1)*(2.*dz(k-1)+dz(k)) + &
                    grad(k  )*(2.*dz(k+1)+dz(k)))  &
                     *dz(k)/(dz(k-1)+dz(k)+dz(k+1))
      ENDDO
    ENDIF
        slope(1) = (grad(2)*(2.*dz(1)+dz(1)) + &
                    grad(1)*(2.*dz(2)+dz(1)))  &
                     *dz(1)/(dz(1)+dz(1)+dz(2))
    ! slope(1) = 2.*grad(2)*dz(1)
    slope(nz) = 2.*grad(nz)*dz(nz)
  
    if (limiters) then
      do k = 1, nz
        if (k >= 2 .and. k <= nz-1) then
          Cmin = min(var(k-1), var(k), var(k+1))
          Cmax = max(var(k-1), var(k), var(k+1))
          slope(k) = sign(1.,slope(k)) *  &
                      min( abs(slope(k)), &
                        2.*(var(k)-Cmin), &
                        2.*(Cmax-var(k))  )   ! Equation 1.8
        else
          slope(k) = 0.  ! always slope=0
        endif
      enddo
    endif 

   END SUBROUTINE Sub_cal_slope

   SUBROUTINE Sub_cal_weights ( dz, zwt )

    IMPLICIT NONE
    REAL, DIMENSION(:),    INTENT(IN)   :: dz
    REAL, DIMENSION(0:,:), INTENT(OUT)  :: zwt
    REAL    :: denom1, denom2, denom3, denom4, num3, num4, x, y
    INTEGER :: k

    DO k = 3, size(dz,1)-1
      denom1 = 1.0/(dz(k-1) + dz(k))
      denom2 = 1.0/(dz(k-2) + dz(k-1) + dz(k) + dz(k+1))
      denom3 = 1.0/(2*dz(k-1) +   dz(k))
      denom4 = 1.0/(  dz(k-1) + 2*dz(k))
      num3   = dz(k-2) + dz(k-1)
      num4   = dz(k)   + dz(k+1)
      x      = num3*denom3 - num4*denom4
      y      = 2.0*dz(k-1)*dz(k)               ! everything up to this point is just
                                               ! needed to compute x1,x1,x3                      
      zwt(0,k) = dz(k-1)*denom1                ! = 1/2 in equally spaced case
      zwt(1,k) = zwt(0,k) + x*y*denom1*denom2  ! = 1/2 in equally spaced case
      zwt(2,k) = dz(k-1)*num3*denom3*denom2    ! = 1/6 ''
      zwt(3,k) = dz(k)*num4*denom4*denom2      ! = 1/6 ''
    ENDDO

   END SUBROUTINE Sub_cal_weights 

END MODULE Mod_dyn_driver
