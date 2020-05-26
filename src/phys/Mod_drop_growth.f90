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
      REAL, DIMENSION(drop_column_num+1), INTENT(OUT) :: dmb_dt !! 2020.05.26 "add +1" by han

      CALL es_Fk_Fd(temp,pres,es,Fk,Fd) 
      
      e     = pres * qv/0.622      ! vapor pressure       [hPa]
      RH    = (e/es)*100.          ! Relative humidity    [%]

      !S     = 0.01                 ! For test
      S     = -0.01                 ! For test

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

      next_nr = 0.

       DO inext_m = 1, drop_column_num
         DO im = 1, drop_column_num-1

          IF (next_m(inext_m) > ref_m(im) .AND. &
              next_m(inext_m) < ref_m(im+1)) THEN

            N  = nr(inext_m)
            NM = nr(inext_m) * next_m(inext_m)

            x  = (NM - (ref_m(im+1) * N)) / (ref_m(im)-ref_m(im+1))
            y  = N - x 

             next_nr(im) = next_nr(im) + x !! modified by han
             next_nr(im+1) = next_nr(im+1) + y !! modified by han

          ENDIF

         ENDDO
       ENDDO

        !print*, next_nr
        write(*,*) sum(nr)
        write(*,*) sum(next_nr)


    END SUBROUTINE reassign!}}}
END MODULE Mod_drop_growth 
