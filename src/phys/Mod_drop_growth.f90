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
      REAL, DIMENSION(drop_column_num), INTENT(OUT) :: dm_dt
      REAL, DIMENSION(drop_column_num), INTENT(OUT) :: dmb_dt

      CALL es_Fk_Fd(temp,pres,es,Fk,Fd) 
      
      e     = pres * qv/0.622      ! vapor pressure       [hPa]
      RH    = (e/es)*100.          ! Relative humidity    [%]

      S     = 0.01                 ! For test

      Vf = 1.; Vfb = 1.
      IF (ventilation_effect) THEN
        call ventilation(temp, Pres, r, Vf)
        call ventilation(temp, Pres, rb, Vfb)
      END IF

      ! write(*,*) "*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***"
      ! write(*,*) "*** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** *** ***"
      !  write(*,*) temp, pres, r 
      ! write(*,*) shape(Vf)
      ! write(*,*) Vf
      ! write(*,*) "xxx " 
      ! write(*,*) shape(Vfb)
      ! write(*,*) Vfb
      ! stop

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
    SUBROUTINE terminal_velocity( T, p , r, Vt )!{{{

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

      eta0    = 1.818e-4     ! [  g * cm^-1 * s^-1 ] , at 20C 
      ! eta0    = eta0*0.1               ! [ kg * m^-1 * s^-1  ]

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

      END SUBROUTINE terminal_velocity!}}}

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
      ! Out
      REAL, DIMENSION(:), INTENT(OUT) :: Vf     ! ventilation effect

      rhov    = p * 100. / ( Rd * T )  ! [kg m^-3] , density of air

      Vt = 0.
      DO irr = 1, size(r)
        call terminal_velocity(T, P, r(irr), Vt(irr))   ! Beard (1976)
      END DO

      ! write(*,*) shape(vt)
      ! write(*,*) "==========" 
      ! write(*,*) r 
      ! open(unit = 111, file = 'test.txt', status = "unknown", &
      !       access="direct", form="unformatted", recl=100*4) 
      !
      ! write(unit = 111, rec=1) vt
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
                dt, nr, next_nr   &  
               )

      IMPLICIT NONE
      ! In
      REAL,                               INTENT(IN)  :: dt                    
      REAL, DIMENSION(drop_column_num),   INTENT(IN)  :: ref_m,        & ! const. mass (1st mass) 
                                                         next_m,       & ! after added dm/dt 
                                                         dm_dt,         &
                                                         nr            
      REAL, DIMENSION(drop_column_num+1), INTENT(IN)  :: ref_mb,       & ! const. b. mass (1st b. mass) 
                                                         next_mb,      & ! after added dmb/dt 
                                                         dmb_dt
      ! Local
      REAL, DIMENSION(drop_column_num)                :: CFL_substep,  &
                                                         local_nr,     &
                                                         dm
      INTEGER                                         :: izz
      ! Out
      REAL, DIMENSION(drop_column_num),   INTENT(OUT) :: next_nr      
      
      DO izz = 1, nz
        dm(izz) = ref_mb(izz+1)-ref_mb(izz)
      ENDDO

      SELECT CASE (redist_option)
        CASE (1) ! reassign
          CALL reassign(nr, ref_m, next_m, next_nr)
        CASE (2) ! PPM
          WHERE (dmb_dt /= 0.)
            CFL_substep = courant_number*(dm/abs(dm_dt))
          ELSEWHERE
            CFL_substep = MAXVAL(CFL_substep)
          END WHERE
          ! substep_dt = MINVAL(CFL_substep)      ; IF (substep_dt == 0) substep_dt = 1
          substep_dt = 0.02                      ; IF (substep_dt == 0) substep_dt = 0.01
          substep_nt = INT(dt/substep_dt)       ; IF (substep_nt == 0) substep_nt = 1

        
          local_nr = nr
          IF ( dt >= substep_dt ) THEN
            DO itt = 1, substep_nt
              CALL Sub_Finite_volume_PPM ( local_nr, q%sfc_dt(1),             &
                                           q%top_dt(1),                       &
                                           dm, drop_column_num, CFL_substep,  &
                                           substep_dt,                        &
                                           dmb_dt,                            &
                                           next_nr                            &
                                                                             ) 
                  ! write(*,*) "dm =", dm
                  ! write(*,*) "dt =", dt
                  ! write(*,*) "drop_column_num =", drop_column_num
                  ! write(*,*) "CFL_substep =", dt*dm_dt/dm
                  ! write(*,*) "dmb_dt =", dmb_dt
                  ! write(*,*) "nr =", nr
                  ! write(*,*) "next_nr =", next_nr
              local_nr=next_nr
            ENDDO
          ELSE
            CALL Sub_Finite_volume_PPM ( nr, q%sfc_dt(1),                   &
                                         q%top_dt(1),                       &
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
      ! Out
      REAL, DIMENSION(drop_column_num), INTENT(OUT) :: next_nr

       DO inext_m = 1, drop_column_num
         DO im = 1, drop_column_num-1

          IF (next_m(inext_m) > ref_m(im) .AND. &
              next_m(inext_m) < ref_m(im+1)) THEN

            N  = nr(inext_m)
            NM = nr(inext_m) * next_m(inext_m)

            x  = (NM - (ref_m(im+1) * N)) / (ref_m(im)-ref_m(im+1))
            y  = N - x 

            ! next_nr(im) = nr(im) + x
            ! next_nr(im+1) = nr(im+1) + y
            next_nr(im) = x
            next_nr(im+1) = y
          ENDIF

         ENDDO
       ENDDO
        ! write(*,*) sum(nr)
        ! write(*,*) sum(next_nr)
        ! stop

    END SUBROUTINE reassign!}}}
END MODULE Mod_drop_growth 
