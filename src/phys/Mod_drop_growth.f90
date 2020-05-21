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
    subroutine terminal_velocity(T, P, radius, Vt) !{{{
    ! Input
    ! - T      : Temperature [K]
    ! - P      : Pressure [hPa]
    ! - radius : [m]
    !
    ! Output
    ! - Vt     : Terminal velocity [m s-1]
    ! 
    ! Reference 
    ! - Beard (1976)
        implicit none
        real, intent(in)  :: T, P
        real, intent(in)  :: radius
        real, intent(out) :: Vt

        real :: d0  ! diameter [um]
        real :: R, g, T0, P0, l0, mu0
        real :: rho_liquid, rho_air, drho, mu
        real :: l, C1, C2, C3, Csc, Cl
        real :: b0, b1, b2, b3, b4, b5, b6
        real :: Da, X, Y, Re, sigma, Bo, Np

        d0 = radius * 2 * 1.e6  ! radius [m] -> diameter [um]

        ! constant
        R   = 287.       ! [J kg-1 K-1]
        g   = 9.8        ! [m s-2]
        T0  = 293.15     ! [K]
        P0  = 1013.25    ! [hPa]
        l0  = 6.62e-6    ! [cm]
        mu0 = 0.0001818  ! [g cm-1 s-1]

        rho_liquid = 1000.      ! [kg m-3] water density
        rho_air    = (P*100.)/(R*T)    ! [kg m-3] air density  ; p = rho R T -> rho = P/RT
        drho = rho_liquid - rho_air    ! drop - air

        ! dynamic viscosity ( Approximate formula, See Yau (1996) - 102p )
        mu = 1.72e-5 * ( 393/(T+120.) ) * ( (T/273)**(3./2.) )   ! [kg m-1 s-1]

        l   = l0 * (mu/mu0) * (P0/P) * sqrt(T/T0)  ! [cm] mean free path of air molecules
        C1  = drho * g / (18*mu)        ! [m-1 s-1] = [kg m-3] * [m s-2] / [kg m-1 s-1]
        Csc = 1 + 2.51*l/d0             ! [dimensionless]

        ! dynamic viscosity ( Approximate formula, See Yau (1996) - 102p )
        mu = 1.72e-5 * ( 393/(T+120.) ) * ( (T/273)**(3./2.) )   ! [kg m-1 s-1]

        l   = l0 * (mu/mu0) * (P0/P) * sqrt(T/T0)  ! [cm] mean free path of air molecules
        C1  = drho * g / (18*mu)        ! [m-1 s-1] = [kg m-3] * [m s-2] / [kg m-1 s-1]
        Csc = 1 + 2.51*l/d0             ! [dimensionless]

        ! Calculate terminal velocity in each regime
        if (d0 .lt. 0.5 ) then
            Vt = 0.     ! ignore
        else if (d0 .le. 19) then
            ! Regime 1
            Vt = C1 * Csc * (d0*1.e-6)**2   ! [m s-1] = [m-1 s-1] [dimensionless] [um^2]
        else if (d0 .le. 1.07e3) then
            ! Regime 2
            b0 = -0.318657e+1; b1 =  0.992696;    b2 = -0.153193e-2
            b3 = -0.987059e-3; b4 = -0.578878e-3; b5 =  0.855176e-4
            b6 = -0.327815e-5
            C2  = 4 * rho_air * drho* g / ( 3 * mu**2 )  ! [] = [kg2 m-6] [m s-2] / [um2]
            Da  = C2 * (d0*1e-6)**3  ! Davies number [kg2 s-2] = [kg2 m-3 s-2] [m-3]
            X   = log( Da )
            Y   = b0 + b1*X + b2*X**2 + b3*X**3 + b4*X**4 + b5*X**5 + b6*X**6
            Re  = Csc * exp(Y)              ! Reynolds number
            Vt  = mu * Re / (rho_air * d0*1.e-6)
        else if (d0 .le. 7.e3) then
            ! Regime 3
            b0 = -0.500015e+1; b1 =  0.523778e+1; b2 = -0.204914e+1
            b3 =  0.475294;    b4 = -0.542819e-1; b5 =  0.238449e-2

            ! surface tension [N m-1]
            ! sigma = 7.5 * 1e-2  ; Yau (1996) - 85p
            ! See Yau (1996) - 6.9 problem
            Cl = -1.55 * 1e-4   ! [N m-1 K-1]
            C2 = 0.118          ! [N m-1]
            sigma = Cl*T + C2   ! Resonable at -20 ~ 20 [K] temperature

            C3  = 4 * drho * g / (3*sigma)
            Bo  = C3 * (d0*1.e-6)**2.       ! modified Bond number
            Np  = sigma**3. * rho_air**2. / (mu**4. * drho * g)
            X   = log( Bo * Np**(1./6.) )
            Y   = b0 + b1*X + b2*X**2 + b3*X**3 + b4*X**4 + b5*X**5

            Re  = Np**(1./6.) * exp(Y)
            Vt  = mu * Re / (rho_air*(d0*1.e-6))
        else
            ! terminal velocity [m s-1] at d0 = 7.e3 [um] (using regime 3)
            d0 = 7.e3
            b0 = -0.500015e+1; b1 =  0.523778e+1; b2 = -0.204914e+1
            b3 =  0.475294;    b4 = -0.542819e-1; b5 =  0.238449e-2

            Cl = -1.55 * 1e-4   ! [N m-1 K-1]
            C2 = 0.118          ! [N m-1]
            sigma = Cl*T + C2   ! Resonable at -20 ~ 20 [K] temperature

            C3  = 4 * drho * g / (3*sigma)
            Bo  = C3 * (d0*1.e-6)**2.       ! modified Bond number
            Np  = sigma**3. * rho_air**2. / (mu**4. * drho * g)
            X   = log( Bo * Np**(1./6.) )
            Y   = b0 + b1*X + b2*X**2 + b3*X**3 + b4*X**4 + b5*X**5

            Re  = Np**(1./6.) * exp(Y)
            Vt  = mu * Re / (rho_air*(d0*1.e-6))
        end if

    end subroutine terminal_velocity!}}}

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
        write(*,*) sum(nr)
        write(*,*) sum(next_nr)
        stop

    END SUBROUTINE reassign!}}}
END MODULE Mod_drop_growth 
