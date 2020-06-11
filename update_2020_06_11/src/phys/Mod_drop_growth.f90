MODULE Mod_drop_growth

  USE Mod_global
  USE Mod_const
  USE Mod_dyn_driver

  IMPLICIT NONE
  
  CONTAINS

    !== !== !==
    SUBROUTINE compute_dmb_dt         &!{{{
               (                      &
                temp, qv, pres,       &
                r, rb,                &
                dm_dt, dmb_dt         &
               )

      IMPLICIT NONE
      ! In
      REAL,                               INTENT(IN)  :: temp, qv, pres
      REAL, DIMENSION(drop_column_num),   INTENT(IN)  :: r
      REAL, DIMENSION(drop_column_num+1), INTENT(IN)  :: rb
      ! Local
      REAL                               :: e, es, RH, S, Fk, Fd
      REAL, DIMENSION(drop_column_num)   :: Vf
      REAL, DIMENSION(drop_column_num+1) :: Vfb
      ! Out
      REAL, DIMENSION(drop_column_num),   INTENT(OUT) :: dm_dt
      REAL, DIMENSION(drop_column_num+1), INTENT(OUT) :: dmb_dt
      INTEGER :: iq

      CALL es_Fk_Fd(temp,pres,es,Fk,Fd) 
      
      e     = pres * qv/0.622      ! vapor pressure       [hPa]
      RH    = (e/es)               ! Relative humidity    

      S     = RH - 1               ! For test

      Vf = 1.; Vfb = 1.
      IF (ventilation_effect) THEN
        call ventilation(temp, Pres, r, Vf)
        call ventilation(temp, Pres, rb, Vfb)
      END IF

      dm_dt  = 4*Pi*r*(1./(Fd+Fk))*S*Vf
      dmb_dt = 4*Pi*rb*(1./(Fd+Fk))*S*Vfb

    END SUBROUTINE compute_dmb_dt         !}}}

    !== !== !==
    SUBROUTINE es_Fk_Fd(temp, pres, es, Fk, Fd) !{{{

      IMPLICIT NONE
      ! In
      REAL, INTENT(IN)  :: temp, pres
      ! Local
      REAL              :: temp_test
      ! Out
      REAL, INTENT(OUT) :: es, Fk, Fd

      es = 6.112 * exp(( 17.67*(temp-273.15) )/( (temp-273.15)+243.5 ))

      Fk = ((L/(Rv*temp))-1.)*(L/(Ka*temp))
      Fd = (Rv*temp) / ((Dv*es*100))

    END SUBROUTINE es_Fk_Fd!}}}

    !== !== !==
    SUBROUTINE terminal_velocity_beard_1976( T, p , r, Vt )!{{{

      IMPLICIT NONE

      ! In
      REAL, INTENT(IN)    :: T, p ,r
      ! Local
      REAL                :: T0, l0, p0, eta0,  &
                             rhov, drho, eta,   & 
                             C, C_sc, N_Da,     &
                             X_rg, Y_rg, N_Re,  &
                             sigma, d0_max, Bo, &
                             Np, l, d

      REAL, DIMENSION(7)  :: bn_rg2 
      REAL, DIMENSION(6)  :: bn_rg3
      ! Out
      REAL, INTENT(OUT)   :: Vt

      !---------------------------- parameter list -----------------------------
      T0      = 293.15                 ! [ K ] , 20C
      l0      = 6.62e-6    ! [ cm ]
      p0      = 1013.25                ! [ hPa = mb ] , standard sea level pressure   
      rhov    = p * 100. / ( Rd * T )  ! [ kg * m^-3  ] , density of air
      drho    = rho - rhov             ! [ kg * m^-3  ] , 
                                       ! difference between fluid density of drop and air

      ! eta0    = 1.818e-4     ! [  g * cm^-1 * s^-1 ] , at 20C 
      eta0 = 1.72 * (1./10.)**(5.) * ( 393. / ( T0 + 120. ) ) * ( ( T0 / 273. )**(3./2.) )

      ! dynamic viscosity, A short course in cloudphysics 3rd, 102 page
      eta     = 1.72e-5 * ( 393 / ( T + 120) ) * (( T / 273 )**(3./2.))

      ! [ cm ] mean free path of air molcules
      l        = l0 * ( eta / eta0 ) * ( p0 / p ) * (( T / T0 )**(1/2))
      ! l        = l * 0.01    ! [ m ]                                         

      d = 2.*r*1e+6

      ! slip correction factor 
      C_sc = 1. + 2.51 * ( l / d )  


      !------------------------------ Regime list ------------------------------
      ! IF ( d <= 5e-7 ) THEN !------------------------------------------ d <=5e-7
      IF ( d <= 0.5 ) THEN !------------------------------------------ d <=5e-7
        Vt = 0.
      ! ELSE IF ( 5e-7 <= d .and. d < 1.9e-5 ) THEN !-------------------- Regime 1 
      ELSE IF ( 0.5 <= d .and. d < 19 ) THEN !-------------------- Regime 1 
        ! [m^-1 * s^-1] = [kg * m^-3] * [m * s^-2] * [kg^-1 * m * s]
        C = ( drho * g ) / ( 18. * eta )

        ! [m * s^-1] = [m^-1 * s^-1] * [dimensionless] * [m^2]
        Vt = C * C_sc * ((d*1e-6)**2.)  

      ! ELSE IF ( 1.9e-5 <= d .and. d < 1.07e-3 ) THEN !------------------ Regime 2 
      ELSE IF ( 19 <= d .and. d < 1.07e+3 ) THEN !------------------ Regime 2 
        ! [m^-3] = [kg * m^-3] * [kg * m^-3]
        !            * [m * s^-2] * [kg^-2 * m^2 * s^2]
        C = 4. * rhov * ( drho ) * g / ( 3. * (eta**2.) )

        ! Davies number
        ! [ dimensionless ] = [ m^-3 ] * [ m^3 ] 
        N_Da = C * ((d*1e-6)**3.)  

        X_rg = log( N_Da )

        bn_rg2 = (/ -0.318657e+1 , 0.992696 , -0.153193e-2, &
                    -0.987059e-3, -0.578878e-3, 0.855176e-4, -0.327815e-5 /)

        Y_rg  = bn_rg2(1) +                &
                bn_rg2(2) * ( X_rg )**1 +  &
                bn_rg2(3) * ( X_rg )**2 +  &
                bn_rg2(4) * ( X_rg )**3 +  &
                bn_rg2(5) * ( X_rg )**4 +  &
                bn_rg2(6) * ( X_rg )**5 +  &
                bn_rg2(7) * ( X_rg )**6     

        ! slip correction factor
        C_sc = 1. + 2.51 * l / d           

        ! reynolds number
        N_Re = C_sc * exp(Y_rg)              

        ! terminal velocity at regime 2
        ! [ m * s^-1 ]    = [ kg * m^-1 * s^-1 ] 
        !                   * [ dimensionless ] * [ kg^-1 * m^3 ] * [ m^-1 ] 
        Vt = eta * N_Re / ( rhov * d*1e-6 )  
      ! ELSE IF ( 1.07e-3 <= d .and. d < 0.007 ) THEN !------------------- Regime 3 

      ELSE IF ( 1.07e+3 <= d .and. d < 7e+3 ) THEN !------------------- Regime 3 
        ! approx. surface tension , A short course in cloudphysics 3rd. 85 page
        ! [ N * m^-1 ] = [ kg * m * s^-2 * m^-1 ] = [ kg * s^-2 ]
        ! sigma = 7.5 * (1e-2) 
        ! Yau (1996) - 6.9 problem                                  
        sigma = (-1.55e-4)*T + 0.118   ! Resonable at -20 ~ 20 [K] temperature

        ! [ m^-2 ]     = [ kg * m^-3 ] * [ m * s^-2 ] * [ kg^-1 * s^2 ] 
        C = 4. *  drho  * g  / ( 3. * sigma )           

        ! modified bond umber 
        ! [ dimensionless ] = [ m^-2 ] * [ m^2 ]
        Bo = C * ( d*1e-6 )**2

        ! physical property number
        ! [ dimensionless ] = [ kg^3 * s^-6 ] * [ kg^2 * m^-6 ]
        !                     * [ kg^-4 * m^4 * s^4 ] * [ kg^-1 * m^3 ] 
        !                     * [ m^-1 * s^2 ] 
        Np = (sigma**3) * (rhov**2) / ( (eta**4) * drho * g )     

        ! [ dimensionless ] = [ dimensionless ] * [ dimensionless ]                                                               
        X_rg     = log( Bo * (Np**(1./6.)) )                      

        bn_rg3    = (/ -0.500015e+1 , 0.523778e+1, -0.204914e+1, &
                          0.475294, -0.542819e-1, 0.238449e-2/)

        Y_rg     = bn_rg3(1) +               &  
                   bn_rg3(2) * ( X_rg )**1 +     &
                   bn_rg3(3) * ( X_rg )**2 +     &
                   bn_rg3(4) * ( X_rg )**3 +     &
                   bn_rg3(5) * ( X_rg )**4 +     &
                   bn_rg3(6) * ( X_rg )**5     

        ! reynolds number
        N_Re = (Np**(1./6.)) * exp(Y_rg) 

        ! terminal velocity at regime 3
        ! [ m * s^-1 ]    = [ kg * m^-1 * s^-1 ] 
        !                   * [ dimensionless ] * [ kg^-1 * m^3 ] * [ m^-1 ]
        Vt  = ( eta * N_Re ) / ( rhov * d*1e-6 ) 
        ! write(*,*) bo, np, x_rg, y_rg, N_Re 
        ! write(*,*) C, c_sc, eta, vt
        ! stop
      ! ELSE IF ( d >= 0.007 ) THEN !------------------------------------ d >=0.007
      ELSE IF ( d >= 7e+3 ) THEN !------------------------------------ d >=0.007
        d = 7e+3
        ! sigma  = 7.5 * (1e-2)                                   
        sigma = (-1.55*1e-4)*T + 0.118   ! Resonable at -20 ~ 20 [K] temperature
        C      = 4. *  drho  * g  / ( 3. * sigma )           
        Bo     = C * ( d*1e-6 )**2
        Np     = (sigma**3) * (rhov**2) / ( (eta**4) * drho * g )     

        bn_rg3    = (/ -0.500015e+1 , 0.523778e+1, -0.204914e+1, &
                          0.475294, -0.542819e-1, 0.238449e-2/)

        X_rg = log( Bo * (Np**(1./6.)) )
        Y_rg = bn_rg3(1) +               & 
                   bn_rg3(2) * ( X_rg )**1 + &
                   bn_rg3(3) * ( X_rg )**2 + &
                   bn_rg3(4) * ( X_rg )**3 + &
                   bn_rg3(5) * ( X_rg )**4 + &
                   bn_rg3(6) * ( X_rg )**5  

        N_Re = (Np**(1./6.)) * exp(Y_rg)

        Vt  = ( eta * N_Re ) / ( rhov * d*1e-6 )     
      ENDIF

      END SUBROUTINE terminal_velocity_beard_1976 !}}}

    !== !== !==
    SUBROUTINE terminal_velocity_beard_1977( T, p, r, Vt )!{{{

      IMPLICIT NONE
      ! In
      REAL,                 INTENT(IN)  :: T, P
      REAL,                 INTENT(IN)  :: r   
      ! Local
      INTEGER                           :: irr
      REAL                              :: l0, T0, p0, rho0, eta0, &
                                           rhov, eta, l,           &
                                           eps_c, eps_s, f
      REAL                              :: V0, D
      ! Local
      REAL,                 INTENT(OUT) :: Vt
  
      l0          = 6.62e-6     ! [ cm ]
      rhov        = p  * 100. / ( Rd * T  )   ! density of air 
      eta         = 1.72e-5 * ( 393. / ( T + 120. ) ) * ( ( T / 273. )**(3./2.) )

      T0          = 293.15                    ! 20 'C ┐ Basic
      p0          = 1013.25                   ! 1 atm ┘ statement

      rho0        = p0 * 100. / ( Rd * T0 )
      eta0        = 1.72e-5 * ( 393. / ( T0 + 120. ) ) * ( ( T0 / 273. )**(3./2.) )

      ! V0
      ! beard 1976 : T = 20C, P = 1013.25
      CALL terminal_velocity_beard_1976( T0 , p0, r, V0 )

      !--------------------- beard 1977 --------------------------------------

      D = 2. * r * 1e+6        ! radius -> diameter & [ meter ] -> [ micrometer ]

      !------------------------- Regime 1 ------------------------------------

      IF ( D < 1. ) THEN
        Vt = 0.             ! [ m/s ]
      ELSE IF ( D >= 1. .and. D < 40. ) THEN
        ! [ m ]
        l = l0 * 0.01 * ( eta / eta0 ) * ( ( p0 * rho0 ) / ( p * rhov ) )**(1./2.)
        ! [ dimensionless ] velocity adjustment factor
        f = ( eta0 / eta ) * ( 1 + 2.51 * ( l / ( D * 1e-6 ) ) ) / &
                                   ( 1 + 2.51 * ( ( l0 * 0.01 ) / ( D * 1e-6 ) ) )

        ! [ m * s^-1 ] = [ dimensionless ] * [ m * s^-1 ] velocity at 1
        ! micrometer ~ 40 micrometer   
        Vt = V0 * f

      !------------------------- Regime 2 ------------------------------------

      !              The unit of input diameter must be [ cm ] !!!
      
      !-----------------------------------------------------------------------

      ELSE IF ( D >= 40  .and. D < 6e+3 ) THEN
        ! [ dimensionless ] = [ kg * m^-3 ] * [ kg^-1 * m^3 ]
        eps_c   = ( ( rho0 / rhov )**(1./2.) ) - 1.
        ! [ dimensionless ] = [ kg * m^-1 * s^-1 ] * [ kg^-1 * m * s ]
        eps_s   = ( eta0 / eta ) - 1.
        ! [ dimensionless ]
        f       = 1.104 * eps_s + ( ( ( 1.058 * eps_c ) - ( 1.104 * eps_s) ) &
                                                * ( 5.52 + log( D*1e-4 ) ) / 5.01 )   + 1.

        Vt      = V0 * f

      ELSE IF ( D >= 6e+3 )  THEN
        ! [ dimensionless ] = [ kg * m^-3 ] * [ kg^-1 * m^3 ]
        eps_c   = ( ( rho0 / rhov )**(1./2.) ) - 1.
        ! [ dimensionless ] = [ kg * m^-1 * s^-1 ] * [ kg^-1 * m * s ]
        eps_s   = ( eta0 / eta ) - 1.
        ! [ dimensionless ]
        f     = 1.104 * eps_s + ( ( ( 1.058 * eps_c ) - ( 1.104 * eps_s) ) &
                                             * ( 5.52 + log( 6e+3*1e-4 ) ) / 5.01 )   + 1.

        Vt    = V0 * f

      ENDIF 

    END SUBROUTINE terminal_velocity_beard_1977!}}}

    !== !== !==
   SUBROUTINE ventilation(T, P, r, Vf) !{{{

      IMPLICIT NONE
      ! In
      REAL,               INTENT(IN )   :: T,   & ! Temperature [K]
                                           P      ! Pressure    [hPa]
      ! Local
      REAL, DIMENSION(:), INTENT(IN )   :: r      ! radius      [m]
      REAL, DIMENSION(SIZE(r))  :: Vt, Re
      REAL                              :: rhov
      REAL                              :: eta 
      INTEGER                           :: irr
      REAL                              :: tt, pp 
      ! Out
      REAL, DIMENSION(:), INTENT(OUT) :: Vf     ! ventilation effect

      rhov    = p * 100. / ( Rd * T )  ! [kg m^-3] , density of air

      Vt = 0.
      ! tt = 293.15 ; pp = 1013.25
      ! tt = 263.15 ; pp = 500
      DO irr = 1, size(r)
        ! call terminal_velocity_beard_1976(tt, pp, r(irr), Vt(irr))   ! Beard (1977)
        call terminal_velocity_beard_1977(T, P, r(irr), Vt(irr))   ! Beard (1977)
      END DO

      ! write(*,*) shape(vt)
      ! write(*,*) "==========" 
      ! write(*,*) r 
      ! open(unit = 111, file = '1976.bin', status = "unknown", &
      ! open(unit = 111, file = 'r.bin', status = "unknown", &
      !       access="direct", form="unformatted", recl=100*4) 
      !
      ! write(unit = 111, rec=1) r
      ! write(unit = 111, rec=2) r
      !
      ! stop
      ! dynamic viscosity of air (See Yau (1996) 102-103p)
      eta = 1.72e-5 * ( 393./(T+120.) ) * ( T/273. )**(3./2.)    ! approximate formula
      ! mu = 1.717e-5          ! [kg m-1 s-1] (at 273 [K])

      ! Reynolds number
      Re = 2*rhov*r*Vt/eta         ! Yau (1996) 116p

      IF ( any(0 <= Re .and. Re < 2.5) ) THEN
        Vf = 1.0 + 0.09*Re
      ELSE !IF (Re > 2.5) then
        Vf = 0.78 + 0.28*sqrt(Re)
      ! ELSE
      !   CALL FAIL_MSG("check compute reynols number")
      END IF
      
    END SUBROUTINE ventilation !}}}

    !== !== !==
    SUBROUTINE compute_redist     &!{{{
               (                  &
                ref_m, next_m,    &
                ref_mb, next_mb,  & 
                dm_dt, dmb_dt,    & 
                dt, dr,           &
                nr, next_nr   &  
               )

      IMPLICIT NONE
      ! In
      REAL,                               INTENT(IN)  :: dt                    
      REAL, DIMENSION(drop_column_num),   INTENT(IN)  :: ref_m,        & ! const. mass (1st mass) 
                                                         next_m,       & ! after added dm/dt 
                                                         dm_dt,        &
                                                         nr,           &          
                                                         dr 
      REAL, DIMENSION(drop_column_num+1), INTENT(IN)  :: ref_mb,       & ! const. b. mass (1st b. mass) 
                                                         next_mb,      & ! after added dmb/dt 
                                                         dmb_dt
      ! Local
      REAL, DIMENSION(drop_column_num)                :: CFL_substep,  &
                                                         local_nr,     &
                                                         dm
      REAL, DIMENSION(drop_column_num+1)              :: ddmb_dt
      INTEGER                                         :: izz
      ! Out
      REAL, DIMENSION(drop_column_num),   INTENT(OUT) :: next_nr      
      
      DO izz = 1, nz
        dm(izz) = ref_mb(izz+1)-ref_mb(izz)
      ENDDO

      local_nr = nr!/(dr*2)  !! local_nr = nr -> nr/dr by han
!print*, '=========================='
!print*, local_nr 
!print*, '=========================='
!print*, nr 
!print*, '=========================='
!print*, dr 

      SELECT CASE (redist_option)
        CASE (1) ! reassign
          CALL reassign(local_nr, ref_m, next_m, next_nr)
        CASE (2) ! PPM
          WHERE (dmb_dt /= 0.)
            CFL_substep = courant_number*(dm/abs(dm_dt))
          ELSEWHERE
            CFL_substep = MAXVAL(CFL_substep)
          END WHERE
          substep_dt = MINVAL(CFL_substep)      ; IF (substep_dt == 0) substep_dt = 1
          substep_dt = 0.01                     ; IF (substep_dt == 0) substep_dt = 0.01
          substep_nt = INT(dt/substep_dt)       ; IF (substep_nt == 0) substep_nt = 1
        
          ddmb_dt = dmb_dt
          ! ddmb_dt(drop_column_num+1) = 0.
          IF ( dt >= substep_dt ) THEN
            DO itt = 1, substep_nt
              CALL Sub_Finite_volume_PPM_forphy ( local_nr, 0.,                      &
                                           0.,                                &
                                           dm, drop_column_num, CFL_substep,  &
                                           substep_dt,                        &
                                           dmb_dt,                            &
                                           next_nr                            &
                                                                             ) 
              !     write(*,*) "dm =", dm
              !     write(*,*) "dt =", dt
              !     write(*,*) "drop_column_num =", drop_column_num
              !     write(*,*) "CFL_substep =", substep_dt*dm_dt/dm
              !     write(*,*) "dmb_dt =", dmb_dt
              !     write(*,*) "nr =", nr
              !     write(*,*) "nr =", sum(nr)
              ! !     write(*,*) "next_nr =", next_nr
              !     write(*,*) "next_nr =", sum(next_nr)
              ! stop
              local_nr=next_nr
            ENDDO
          ELSE
            CALL Sub_Finite_volume_PPM_forphy ( local_nr, 0.,                      & !! nr -> local_nr by han
                                         0.,                                &
                                         dm, drop_column_num, CFL_substep,  &
                                         dt,                                &
                                         dmb_dt,                            &
                                         next_nr                            &
                                                                           )
                  ! write(*,*) "dm =", dm
                  ! write(*,*) "dt =", dt
                  ! write(*,*) "drop_column_num =", drop_column_num
                  ! write(*,*) "CFL_substep =", CFL_substep*dm_dt/dm
                  ! write(*,*) "dmb_dt =", dmb_dt
                  ! write(*,*) "nr =", nr
                  ! write(*,*) "next_nr =", next_nr
          ENDIF
        CASE DEFAULT
          CALL FAIL_MSG("phys redistribution error")
      END SELECT
!      next_nr = next_nr*2*dr !!2020.06.03 by han
                  ! write(*,*) "next_nr =", sum(next_nr)
                  ! stop

    END SUBROUTINE compute_redist!}}}

    !== !== !==
    SUBROUTINE reassign(nr, ref_m, next_m, next_nr)!{{{

      IMPLICIT NONE
      ! In
      REAL, DIMENSION(drop_column_num), INTENT(IN)  :: nr, ref_m, next_m
      ! Local
      INTEGER                                       :: im, inext_m
      REAL                                          :: x, y,        &
                                                       N, NM
      REAL                                          :: int_x, int_y
      REAL                                          :: e_x, e_y     
      ! Out
      REAL, DIMENSION(drop_column_num), INTENT(OUT) :: next_nr

      next_nr = 0.

       DO inext_m = 1, drop_column_num
         DO im = 1, drop_column_num-1

          IF (next_m(inext_m) > ref_m(im) .AND. &
              next_m(inext_m) < ref_m(im+1)) THEN

            N  = nr(inext_m)
            NM = nr(inext_m) * next_m(inext_m)

            x  = (NM - (ref_m(im+1) * N)) / (ref_m(im)-ref_m(im+1))
            y  = N - x 

            int_x = AINT(x)
            int_y = AINT(y)
            e_x   = x - int_x
            e_y   = y - int_y

            IF (e_x >= e_y) THEN
              x = x + e_y
              y = int_y
            ELSE IF (e_x < e_y) THEN
              y = y + e_x
              x = int_x
            ENDIF

            next_nr(im)   = next_nr(im)   + x 
            next_nr(im+1) = next_nr(im+1) + y



          ENDIF

         ENDDO
       ENDDO



        !print*, next_nr
!         write(*,*) sum(nr*ref_m)
!         write(*,*) sum(next_nr*ref_m)

    END SUBROUTINE reassign!}}}

    !== !== !==
  SUBROUTINE Sub_Finite_volume_PPM_forphy    &!{{{
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

    IF (.false.) THEN
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
    !
    ! B.C = 0. 
    ! flux(ks)   = 0. 
    ! flux(ke+1) = 0. 
    !
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
        rst = ( var_right(kk) - 0.5*xx*(rm - (1.0 - tt*xx)*r6) ) / dz(k-1)
         
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
        rst = ( var_left(kk) + 0.5*xx*(rm + (1.0 - tt*xx)*r6) ) / dz(k)
        ! extension for Courant numbers > 1
        IF (cn > 1.) rst = (xx*rst + rsum)/cn
      ENDIF
      flux(k) = w(k)*rst
    ENDDO
    flux(ks)   = 0. 
    flux(ke+1) = 0.
    ! Cal. FV
    DO i = 1, nz
      dvar        = - (flux(i+1) - flux(i)) !/ dz(i)
      ! dvar        = - (flux(i+1)/dz(i+1) - flux(i)/dz(i))
      next_var(i) = var(i) + dvar * dt
      IF ( next_var(i) .lt. 0. ) THEN !! mass conservation filter

      ENDIF
    END DO

  END SUBROUTINE Sub_Finite_volume_PPM_forphy


  ! SUBROUTINE Sub_cal_slope ( var, dz, nz, slope)
  !
  !   INTEGER,               INTENT(IN)   :: nz 
  !   REAL,    DIMENSION(:), INTENT(IN)   :: var, dz
  !   REAL,    DIMENSION(:), INTENT(OUT)  :: slope
  !
  !   REAL    :: grad(2:nz)
  !   REAL    :: cmin, cmax
  !   INTEGER :: k
  !   LOGICAL :: limiters, dolinear
  !
  !   limiters = .true.
  !   dolinear = .false.
  !    
  !   ! compute slope (weighted for unequal levels)
  !   DO k = 2, nz
  !     grad(k) = (var(k)-var(k-1))/(dz(k)+dz(k-1))
  !   ENDDO
  !   IF (dolinear) THEN
  !     DO k = 2, nz-1
  !       slope(k) = (grad(k+1)+grad(k))*dz(k)
  !     ENDDO
  !   ELSE
  !     DO k = 2, nz-1
  !       slope(k) = (grad(k+1)*(2.*dz(k-1)+dz(k)) + &
  !                   grad(k  )*(2.*dz(k+1)+dz(k)))  &
  !                    *dz(k)/(dz(k-1)+dz(k)+dz(k+1))
  !     ENDDO
  !   ENDIF
  !       slope(1) = (grad(2)*(2.*dz(1)+dz(1)) + &
  !                   grad(1)*(2.*dz(2)+dz(1)))  &
  !                    *dz(1)/(dz(1)+dz(1)+dz(2))
  !   ! slope(1) = 2.*grad(2)*dz(1)
  !   slope(nz) = 2.*grad(nz)*dz(nz)
  !
  !   if (limiters) then
  !     do k = 1, nz
  !       if (k >= 2 .and. k <= nz-1) then
  !         Cmin = min(var(k-1), var(k), var(k+1))
  !         Cmax = max(var(k-1), var(k), var(k+1))
  !         slope(k) = sign(1.,slope(k)) *  &
  !                     min( abs(slope(k)), &
  !                       2.*(var(k)-Cmin), &
  !                       2.*(Cmax-var(k))  )   ! Equation 1.8
  !       else
  !         slope(k) = 0.  ! always slope=0
  !       endif
  !     enddo
  !   endif 
  !
  !  END SUBROUTINE Sub_cal_slope
  !
  !  SUBROUTINE Sub_cal_weights ( dz, zwt )
  !
  !   IMPLICIT NONE
  !   REAL, DIMENSION(:),    INTENT(IN)   :: dz
  !   REAL, DIMENSION(0:,:), INTENT(OUT)  :: zwt
  !   REAL    :: denom1, denom2, denom3, denom4, num3, num4, x, y
  !   INTEGER :: k
  !
  !   DO k = 3, size(dz,1)-1
  !     denom1 = 1.0/(dz(k-1) + dz(k))
  !     denom2 = 1.0/(dz(k-2) + dz(k-1) + dz(k) + dz(k+1))
  !     denom3 = 1.0/(2*dz(k-1) +   dz(k))
  !     denom4 = 1.0/(  dz(k-1) + 2*dz(k))
  !     num3   = dz(k-2) + dz(k-1)
  !     num4   = dz(k)   + dz(k+1)
  !     x      = num3*denom3 - num4*denom4
  !     y      = 2.0*dz(k-1)*dz(k)               ! everything up to this point is just
  !                                              ! needed to compute x1,x1,x3                      
  !     zwt(0,k) = dz(k-1)*denom1                ! = 1/2 in equally spaced case
  !     zwt(1,k) = zwt(0,k) + x*y*denom1*denom2  ! = 1/2 in equally spaced case
  !     zwt(2,k) = dz(k-1)*num3*denom3*denom2    ! = 1/6 ''
  !     zwt(3,k) = dz(k)*num4*denom4*denom2      ! = 1/6 ''
  !   ENDDO
  !
  !  END SUBROUTINE Sub_cal_weights !}}}
  !


END MODULE Mod_drop_growth 

