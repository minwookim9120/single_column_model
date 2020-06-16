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
                dm_dt, dmb_dt, S      &
               )

      IMPLICIT NONE
      ! In
      REAL,                               INTENT(IN)  :: temp, qv, pres
      REAL, DIMENSION(drop_column_num),   INTENT(IN)  :: r
      REAL, DIMENSION(drop_column_num+1), INTENT(IN)  :: rb
      ! Local
      REAL                               :: e, es, RH, Fk, Fd,q ,qs
      REAL, DIMENSION(drop_column_num)   :: Vf
      REAL, DIMENSION(drop_column_num+1) :: Vfb
      ! Out
      REAL, DIMENSION(drop_column_num),   INTENT(OUT) :: dm_dt
      REAL, DIMENSION(drop_column_num+1), INTENT(OUT) :: dmb_dt
      REAL,                               INTENT(OUT) :: S      
      INTEGER :: iq

      CALL es_Fk_Fd(temp,pres,es,Fk,Fd) 
 
      e     = pres * qv/0.622      ! vapor pressure       [hPa]
      RH    = (e/es)               ! Relative humidity    

      S     = RH - 1               ! For test
      ! S     = 0.01               ! For test

      ! write(*,*) es, e, S 

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
      Fd = ( Rv*temp ) / ( (Dv*(1000./Pres)) * (es*100.) )

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

      WHERE (Vf .gt. 5.0)
       VF = 5.0
      END WHERE

      
    END SUBROUTINE ventilation !}}}

    !== !== !==
    SUBROUTINE compute_redist     &!{{{za
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
      
      DO izz = 1, drop_column_num
        dm(izz) = ref_mb(izz+1)-ref_mb(izz)
      ENDDO

      local_nr = nr!/(dr*2)  !! local_nr = nr -> nr/dr by han
      SELECT CASE (redist_option)
        CASE (1) ! reassign
          local_nr=local_nr*drop%dr
          CALL reassign(local_nr, ref_m, next_m, next_nr)
          ! CALL redistribution( ref_m, dm_dt, local_nr, next_nr ) 

          next_nr=next_nr/drop%dr
        CASE (2) ! PPM
          WHERE (dmb_dt /= 0.)
            CFL_substep = courant_number*(dm/abs(dm_dt))
          ELSEWHERE
            CFL_substep = MAXVAL(CFL_substep)
          END WHERE
          substep_dt = MINVAL(CFL_substep)      ; IF (substep_dt == 0) substep_dt = 1
          substep_dt = 0.01                      ; IF (substep_dt == 0) substep_dt = 0.01
          substep_nt = INT(dt/substep_dt)       ; IF (substep_nt == 0) substep_nt = 1

!           write(*,*) dmb_dt
!           write(*,*) dm
! stop
          ddmb_dt = dmb_dt
          ! ddmb_dt(drop_column_num+1) = 0.
          IF ( dt >= substep_dt ) THEN
            DO itt = 1, substep_nt
              CALL Sub_Finite_volume_PPM_forphy ( substep_dt, dmb_dt, dm, local_nr, next_nr )
                  ! write(*,*) "dm =", dm
                  ! write(*,*) "dt =", dt
                  ! write(*,*) "drop_column_num =", drop_column_num
                  ! write(*,*) "CFL_substep =", substep_dt*dm_dt/dm
                  ! write(*,*) "dmb_dt =", dmb_dt
                  ! write(*,*) "nr =", nr
              !     write(*,*) "nr =", sum(nr)
              ! !     write(*,*) "next_nr =", next_nr
              !     write(*,*) "next_nr =", sum(next_nr)
              ! stop
              local_nr=next_nr
            ENDDO
          ELSE
            CALL Sub_Finite_volume_PPM_forphy ( dt, dmb_dt, dm, nr, next_nr )
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

    subroutine redistribution( mass, dm_dt, Nr, next_n ) !{{{
        implicit none
        real, dimension(:) ,intent(in)    :: mass
        real, dimension(:) ,intent(in)    :: dm_dt
        real, dimension(:) ,intent(in) :: Nr

        integer                   :: i, j, nbin
        real, dimension(size(Nr)),intent(out) :: next_N
        real, dimension(size(Nr)) :: x, y
        real, dimension(size(mass)) :: next_m

        next_m = mass + dm_dt*dt

        nbin = size(Nr)
        y    = 0.
        ! next_N = 0.

        do i = 1, nbin-1
            ! print*,
            ! "i=",i,"mass(i)=",mass(i),"next_m(i)",next_m(i),"mass(i+1)=",mass(i+1)
            if ( next_m(i) <= 0. ) then
                next_N(i) = 0.
            else
                do j = 1, nbin-1
                    if ( mass(j) < next_m(i) ) then
                        if ( mass(j+1) > next_m(i) ) then
                            x(i)      = ( Nr(i)*(next_m(i)-mass(j+1)) ) &
                                      / (          mass(j)-mass(j+1) )
                            y(i+1)    = Nr(i) - x(i)
                            next_N(i) = x(i) + y(i)
                            exit
                        end if
                    end if
                end do
            end if
            ! print*,"i=",i,"Nr(i)=",Nr(i),"x=",x(i),"y=",y(i),"next_N(i)=",next_N(i)
            ! print*, mass(i),mass(i+1),next_N(i)
        end do
        next_N(nbin) = Nr(nbin)+y(nbin)

        ! print*, "TODO! redistribution"
        ! print*, "     Porting from ncl code... (+
        ! doc/redistribution/reassign.ncl)"
        ! stop

        ! Nr = next_N

    end subroutine redistribution   !}}}

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


    subroutine Sub_Finite_volume_PPM_forphy ( dt, w_half, dz, C, next_C )  ! {{{
    !-- Input
    ! dt     = length of time step
    ! w_half = vertical velocity at half coordinate(nz+1)
    ! nt     = size of time for iteration
    ! dz     = depth of model layers
    !
    ! from namelist
    ! dyn_adv_scheme = differencing scheme, use one of these values:
    !  1: finite_difference = second-order centered finite difference
    !  2: finite_volume     = finite volume method
    !  3: PPM               = piecewise parabolic method
    !                        See Colella and Woodward (1984)
    !  4: PPM               = PPM but Lin (2003) limeter used
    !
    !-- Output
    ! C = advected quantity
    !
    ! Note! Here, flux form is used for advection term
    ! FLUX_FORM      = solves for -d(wr)/dt
    ! ADVECTIVE_FORM = solves for -w*d(r)/dt
    !
    ! Here, we use Lorenz configuration
    ! See Figure 1 in Holdaway et al., (2012) 
    ! https://rmets.onlinelibrary.wiley.com/doi/epdf/10.1002/qj.2016

    implicit none
    real,               intent(in   ) :: dt
    real, dimension(:), intent(in   ) :: w_half
    real, dimension(:), intent(in   ) :: dz
    real, dimension(:), intent(in   ) :: C
    real, dimension(:), intent(out)   :: next_C

    integer :: k, kk
    integer :: ks, ke, kstart, kend
    real    :: wgt   ! weight dz
    real    :: Cwgt  ! weight variable
    real    :: dC_dt, zbottom = 0.
    real    :: tt, cn, Csum, dzsum, dtw
    real    :: xx, a, b, Cm, C6, Cst, Cdt
    real, dimension(0:3,size(C)) :: zwt
    real, dimension(size(C))     :: slp, C_left, C_right
    real, dimension(size(C))     :: slope
    real, dimension(size(C)+1)   :: flux
    character(len=20)       :: eqn_form = "FLUX_FORM"
    ! character(len=20)       :: eqn_form = "ADVECTIVE_FORM"
    logical :: linear, test_1
    logical :: do_outflow_bnd = .true.

    ! set default values for optional arguments
    ! if (phy_adv_scheme == 1) do_outflow_bnd = .false.

    ! vertical indexing
    ks = 1;  ke = size(C)

    ! start and end indexing for finite volume fluxex
    kstart = ks+1; kend = ke
    if (do_outflow_bnd) then
        kstart = ks;  kend = ke+1
    end if

    ! Make stagged grid for advection
    if ( size(w_half) /= size(C)+1 )then
        call fail_msg("vertical dimension of input arrays inconsistent")
    end if

    ! Set Boundary Condition (Homogeneous Dirichlet BC)
    ! most likely w = 0 at these points
    if (do_outflow_bnd) then
        flux(ks)   = 0.
        flux(ke+1) = 0.
    else
        flux(ks)   = w_half(ks)*C(ks)
        flux(ke+1) = w_half(ke+1)*C(ke)
    end if


        ! 3) Piecewise Parabolic Method, Colella and Woodward (1984) {{{
            zwt = 0 ! Note! Avoid for Nan value occur
            call compute_weights(dz, zwt)
            call slope_z(C, dz, slp, linear=.false.)        ! Equation 1.7
            do k = ks+2, ke-1
                C_left(k) = C(k-1) + zwt(1,k)*(C(k)-C(k-1)) &
                                   - zwt(2,k)*slp(k)        &
                                   + zwt(3,k)*slp(k-1)      ! Equation 1.6 
                C_right(k-1) = C_left(k)
                ! Or, we can use Equation 1.9 (Need condition)
                ! C_rihgt(k) = (7./12.)*(a(k)+a(k+1)) - (1./12.)*(a(k+2)+a(k-1))
                ! coming out of this loop, all we need is r_left and r_right
            enddo

            ! boundary values  ! masks ???????
            C_left (ks+1) = C(ks+1) - 0.5*slp(ks+1)
            C_right(ke-1) = C(ke-1) + 0.5*slp(ke-1)

            ! pure upstream advection near boundary
            ! r_left (ks) = r(ks)
            ! r_right(ks) = r(ks)
            ! r_left (ke) = r(ke)
            ! r_right(ke) = r(ke)

            ! make linear assumption near boundary
            ! NOTE: slope is zero at ks and ks therefore
            !       this reduces to upstream advection near boundary
            C_left (ks) = C(ks) - 0.5*slp(ks)
            C_right(ks) = C(ks) + 0.5*slp(ks)
            C_left (ke) = C(ke) - 0.5*slp(ke)
            C_right(ke) = C(ke) + 0.5*slp(ke)

                ! limiters from Lin (2003), Equation 6 (relaxed constraint)
                do k = ks, ke
                C_left (k) = C(k) - sign( min(abs(slp(k)),       &
                                          abs(C_left (k)-C(k))), &
                                          slp(k) )  ! (B3)
                C_right(k) = C(k) + sign( min(abs(slp(k)),       &
                                          abs(C_right(k)-C(k))), &
                                          slp(k) )  ! (B4)
                enddo

            ! compute fluxes at interfaces {{{
            tt = 2./3.
            do k = kstart, kend ! ks+1, nz
                if (w_half(k) >= 0.) then ! w = positive
                    if (k == ks) cycle    ! inflow
                        cn = dt*w_half(k)/dz(k-1)   ! Courant number
                        kk = k-1
                    ! extension for Courant numbers > 1
                    if (cn > 1.) then
                        Csum = 0.; dzsum = 0.
                        dtw  = dt*w_half(k)
                        do while (dzsum+dz(kk) < dtw)
                            if (kk == 1) exit
                            dzsum = dzsum + dz(kk)
                            Csum  =  Csum +  C(kk)
                            kk    =    kk -1
                        enddo
                        xx = (dtw-dzsum)/dz(kk)
                    else
                        xx = cn     ! y = u*dt (1.13)
                    endif
                    Cm = C_right(kk) - C_left(kk)
                    C6 = 6.0*(C(kk) - 0.5*(C_right(kk) + C_left(kk)))   ! (1.5)
                    if (kk == ks) C6 = 0.
                    Cst = C_right(kk) - 0.5*xx*(Cm - (1.0 - tt*xx)*C6)  ! (1.12)
                    ! extension for Courant numbers > 1
                    if (cn > 1.) Cst = (xx*Cst + Csum)/cn
                else                      ! w = negative
                    if (k == ke+1) cycle  ! inflow
                    cn = - dt*w_half(k)/dz(k)
                    kk = k
                    ! extension for Courant numbers > 1
                    if (cn > 1.) then
                        Csum = 0.; dzsum = 0.
                        dtw  = -dt*w_half(k)
                        do while (dzsum+dz(kk) < dtw)
                            if (kk == ks) exit
                            dzsum = dzsum + dz(kk)
                            Csum  =  Csum +  C(kk)
                            kk    =    kk + 1
                        enddo
                        xx = (dtw-dzsum)/dz(kk)
                    else
                        xx = cn
                    endif
                    Cm = C_right(kk) - C_left(kk)
                    C6 = 6.0*(C(kk) - 0.5*(C_right(kk) + C_left(kk)))
                    if (kk == ke) C6 = 0.
                    Cst = C_left(kk) + 0.5*xx*(Cm + (1.0 - tt*xx)*C6)
                    ! extension for Courant numbers > 1
                    if (cn > 1.) Cst = (xx*Cst + Csum)/cn
                endif
                flux(k) = w_half(k)*Cst
                ! if (xx > 1.) cflerr = cflerr+1
                ! cflmaxx = max(cflmaxx,xx)
                ! cflmaxc = max(cflmaxc,cn)
                ! }}}
            enddo ! }}}
    ! vertical advective tendency
    select case (eqn_form)
        case ("FLUX_FORM")
            do k = ks, ke
                 dC_dt     = - ( flux(k+1) - flux(k) ) / dz(k)
                !dC_dt     = - ( flux(k+1)/dz(k+1) - flux(k)/dz(k) )
                next_C(k) = C(k) + dC_dt * dt
            end do
        case ("ADVECTIVE_FORM")
            do k = ks, ke
                dC_dt     = - ( flux(k+1) - flux(k) - &
                                C(k)*(w_half(k+1)-w_half(k)) ) / dz(k)
                next_C(k) = C(k) + dC_dt * dt
            end do
        case default
            call fail_msg("No setup physics equation form.")
    end select

    end subroutine Sub_Finite_volume_PPM_forphy  ! }}}

    subroutine slope_z(C, dz, slope, limit, linear) ! {{{
    real, dimension(:), intent(in)  :: C, dz
    real, dimension(:), intent(out) :: slope
    logical,  optional, intent(in)  :: limit, linear

    integer :: k, n
    real    :: grad(2:size(C))
    real    :: Cmin, Cmax
    logical :: limiters = .true.
    logical :: dolinear = .true.

    if (present( limit)) limiters = limit
    if (present(linear)) dolinear = linear

    n = size(C)

    ! compute slope (weighted for unequal levels)
    do k = 2, n
        grad(k) = (C(k)-C(k-1))/(dz(k)+dz(k-1))
    enddo
    if (dolinear) then
        do k = 2, n-1
            slope(k) = (grad(k+1)+grad(k))*dz(k)
        enddo
    else
        do k = 2, n-1
            slope(k) = ( grad(k+1)*(2.*dz(k-1)+dz(k)) + &   ! Equation 1.7
                         grad(k  )*(2.*dz(k+1)+dz(k)) ) * dz(k) &
                     / (   dz(k-1) + dz(k) + dz(k+1)  )
        enddo
    endif
    slope(1) = 2.*grad(2)*dz(1)
    slope(n) = 2.*grad(n)*dz(n)

    ! apply limiters to slope
    if (limiters) then
        do k = 1, n
            if (k >= 2 .and. k <= n-1) then
                Cmin = min(C(k-1), C(k), C(k+1))
                Cmax = max(C(k-1), C(k), C(k+1))
                slope(k) = sign(1.,slope(k)) *  &
                            min( abs(slope(k)), &
                                2.*(C(k)-Cmin), &
                                2.*(Cmax-C(k))  )   ! Equation 1.8
            else
                slope(k) = 0.  ! always slope=0
            endif
        enddo
    endif

    end subroutine slope_z  ! }}}

    subroutine compute_weights(dz, zwt) ! {{{
    real, intent(in),  dimension(:)            :: dz
    real, intent(out), dimension(0:3,size(dz)) :: zwt
    real    :: denom1, denom2, denom3, denom4, num3, num4, x, y
    integer :: k

    do k = 3, size(dz)-1
        denom1 = 1.0/(  dz(k-1) +   dz(k))
        denom2 = 1.0/(  dz(k-2) +   dz(k-1) + dz(k) + dz(k+1))
        denom3 = 1.0/(2*dz(k-1) +   dz(k))  
        denom4 = 1.0/(  dz(k-1) + 2*dz(k))  
        num3   = dz(k-2) + dz(k-1)          
        num4   = dz(k)   + dz(k+1)        
        x      = num3*denom3 - num4*denom4        
        y      = 2.0*dz(k-1)*dz(k)  ! everything up to this point is just
                                    ! needed to compute x1,x1,x3                      
        zwt(0,k) = dz(k-1)*denom1               ! = 1/2 in equally spaced case
        zwt(1,k) = zwt(0,k) + x*y*denom1*denom2 ! = 1/2 in equally spaced case
        zwt(2,k) = dz(k-1)*num3*denom3*denom2   ! = 1/6 ''
        zwt(3,k) = dz(k)*num4*denom4*denom2     ! = 1/6 ''
    enddo

    end subroutine compute_weights  ! }}}

!     !== !== !==
! !   SUBROUTINE Sub_Finite_volume_PPM_forphy    &!{{{
! !              (                        &
! !                var, sfc_var,          &
! !                top_var,               &
! !                dz, nz, CFL,           &
! !                dt,                    &
! !                w,                     &
! !                next_var               &
! !              ) 
! !
! !
! !
! !     INTEGER,                    INTENT(IN) :: nz
! !     REAL,                       INTENT(IN) :: dt
! !     REAL,    DIMENSION(:),      INTENT(IN) :: dz, var
! !     REAL,    DIMENSION(1:nz+1), INTENT(IN) :: w
! !     REAL,    DIMENSION(nz),     INTENT(IN) :: CFL                                   
! !     REAL,    DIMENSION(nz)                 :: next_var
! !     REAL,    DIMENSION(nz)                 :: slp,                               & 
! !                                               var_left, var_right
! !     REAL,    DIMENSION(nz+1)               :: flux
! !     REAL,    DIMENSION(0:3,nz)             :: zwt
! !     REAL                                   :: xx, a, b, rm, r6, rst, wt
! !     REAL                                   :: tt, c1, c2
! !     REAL                                   :: sfc_var, top_var 
! !     LOGICAL                                :: test_1
! !     INTEGER                                :: i, j, k, ks, ke
! !     REAL                                   :: cn, rsum, dzsum, dtw
! !     REAL                                   :: cflerr, cflmaxx, cflmax, cflmaxcc
! !     REAL                                   :: dvar, c
! !     INTEGER                                :: kk
! !
! !     ke=nz
! !     ks=1
! !
! !     Call Sub_cal_weights ( dz, zwt )
! !     CALL Sub_cal_slope ( var, dz, nz, slp )
! !
! !     DO k = 3, nz-1
! !       var_left(k) = var(k-1) + zwt(1,k)*(var(k)-var(k-1)) &
! !                        - zwt(2,k)*slp(k) + zwt(3,k)*slp(k-1)
! !       var_right(k-1) = var_left(k)
! !     ENDDO
! !    
! !    ! boundary values  
! !     var_left (1) = var(1) - 0.5*slp(1)
! !     var_right(1) = var(1) + 0.5*slp(1)
! !     var_left (nz) = var(nz) - 0.5*slp(nz)
! !     var_right(nz) = var(nz) + 0.5*slp(nz)
! !
! !     ! make linear assumption near boundary
! !     ! var_left (2) = var(2) - 0.5*slp(2)
! !     ! var_right(nz-1) = var(nz-1) + 0.5*slp(nz-1)
! !     var_left (2) = var_right(1) 
! !     var_right(nz-1) = var_left(nz) 
! !
! !     IF (.false.) THEN
! !       ! limiters from Lin (2003), Equation 6 (relaxed constraint)
! !       DO k = 1, nz
! !         var_left (k) = var(k) - sign(min(abs(slp(k)),abs(var_left(k)-var(k))), slp(k) )
! !         var_right(k) = var(k) + sign(min(abs(slp(k)),abs(var_right(k)-var(k))), slp(k) )
! !       ENDDO
! !     ELSE
! !       ! limiters from Colella and Woodward (1984), Equation 1.10
! !       DO k = ks, ke
! !         test_1 = (var_right(k)-var(k))*(var(k)-var_left(k)) <= 0.0
! !         IF (test_1) THEN
! !           var_left(k)  = var(k)
! !           var_right(k) = var(k)
! !         ENDIF
! !         IF (k == ks .or. k == ke) CYCLE
! !           rm = var_right(k) - var_left(k)
! !           a = rm*(var(k) - 0.5*(var_right(k) + var_left(k)))
! !           b = rm*rm/6.
! !         IF (a >  b) var_left (k) = 3.0*var(k) - 2.0*var_right(k)
! !         IF (a < -b) var_right(k) = 3.0*var(k) - 2.0*var_left (k)
! !       ENDDO
! !     ENDIF
! !
! !     ! compute fluxes at interfaces
! !     ! flux(ks)   = w(ks)  *var(ks)/dz(1)
! !     ! flux(ke+1) = w(ke+1)*var(ke)/dz(nz+1)
! !     !
! !     ! B.C = 0. 
! !     ! flux(ks)   = 0. 
! !     ! flux(ke+1) = 0. 
! !     !
! !      c   = dt*w(1)/dz(1)
! !      rst = (sfc_var + 0.5*slp(1)*(1.-c)) !/ dz(1)
! !      flux(1) = w(ks)*rst 
! !      rst = (top_var + 0.5*slp(nz)*(1.-c)) !/ dz(nz)
! !      flux(ke+1) = w(ke+1)*rst
! !
! !     ! Cal. flux
! !     tt = 2./3.
! !     DO k = 1, nz+1
! !       IF (w(k) >= 0.) THEN
! !         IF (k == ks) CYCLE ! inflow
! !         cn = dt*w(k)/dz(k-1)
! !         kk = k-1
! !         ! extension for Courant numbers > 1
! !         IF (cn > 1.) THEN
! !           rsum = 0.
! !           dzsum = 0.
! !           dtw = dt*w(k)
! !           DO WHILE (dzsum+dz(kk) < dtw)
! !             IF (kk == 1) THEN
! !               exit
! !             ENDIF
! !             dzsum = dzsum + dz(kk)
! !              rsum =  rsum +  var(kk)
! !             kk = kk-1
! !             ! write(*,*) KK
! !             ! stop
! !           ENDDO
! !           xx = (dtw-dzsum)/dz(kk)
! !         ELSE
! !           xx = cn
! !         ENDIF
! !         rm = var_right(kk) - var_left(kk)
! !         r6 = 6.0*(var(kk) - 0.5*(var_right(kk) + var_left(kk)))
! !         IF (kk == ks) r6 = 0.
! !         rst = ( var_right(kk) - 0.5*xx*(rm - (1.0 - tt*xx)*r6) ) !/ dz(k)
! !          
! !         ! extension for Courant numbers > 1
! !         IF (cn > 1.) rst = (xx*rst + rsum)/cn
! !         ! write(*,*) kk
! !         ! write(*,*) rst
! !         ! write(*,*) var_right(kk)
! !         
! !       ELSE
! !         IF (k == ke+1) CYCLE ! inflow
! !         cn = - dt*w(k)/dz(k)
! !         kk = k
! !         ! extension for Courant numbers > 1
! !         IF (cn > 1.) THEN
! !           rsum = 0.
! !           dzsum = 0.
! !           dtw = -dt*w(k)
! !           DO WHILE (dzsum+dz(kk) < dtw)
! !             IF (kk == ks) THEN
! !               EXIT
! !             ENDIF
! !             dzsum = dzsum + dz(kk)
! !              rsum =  rsum + var(kk)
! !             kk = kk+1
! !           ENDDO
! !           xx = (dtw-dzsum)/dz(kk)
! !         ELSE
! !           xx = cn
! !         ENDIF
! !         rm = var_right(kk) - var_left(kk)
! !         r6 = 6.0*(var(kk) - 0.5*(var_right(kk) + var_left(kk)))
! !         IF (kk == ke) r6 = 0.
! !         rst = ( var_left(kk) + 0.5*xx*(rm + (1.0 - tt*xx)*r6) ) !/ dz(k)
! !         ! extension for Courant numbers > 1
! !         IF (cn > 1.) rst = (xx*rst + rsum)/cn
! !       ENDIF
! !       flux(k) = w(k)*rst !/ dz(k)
! !       write(*,*) k
! !       write(*,*) flux(k)
! !     ENDDO
! !     ! write(*,*) var_right
! !     ! write(*,*) w
! !
! !     ! stop   
! !    ! flux(ks)   = 0. 
! !    ! flux(ke+1) = 0.
! !    ! Cal. FV
! !     DO i = 1, nz
! !        dvar        = - (flux(i+1) - flux(i)) / dz(i)
! !       ! dvar        = - ((flux(i+1) - flux(i)) - &
! !       !                 var(i)*(w(i+1)-w(i))) / dz(i)
! !       ! dvar        = - (flux(i+1)/dz(i+1) - flux(i)/dz(i))
! !       next_var(i) = var(i) + dvar * dt
! !       IF ( next_var(i) .lt. 0. ) THEN !! mass conservation filter
! !         CALL FAIL_MSG("problem : phys ppm")
! !       ENDIF
! !     END DO
! !
! !   END SUBROUTINE Sub_Finite_volume_PPM_forphy
! !
END MODULE Mod_drop_growth 
