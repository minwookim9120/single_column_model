MODULE Mod_init_driver

  USE Mod_global
  USE Mod_const
  USE Mod_read 
  USE Mod_idealcase
  USE Mod_realcase
  USE Mod_distribution 

  IMPLICIT NONE

    CONTAINS
    
    SUBROUTINE Sub_init_vars
    
    IMPLICIT NONE

      !! prescribed vertical-wind

      IF ( incase == 2 ) THEN
        CALL Sub_ideal_init         &
            (                       &
             nz, dz%dz, gamma_dry,  &
             283.15,                &
             temp%dz, q%dz          &
            )  
      ELSE IF ( incase == 1 ) THEN
        CALL Sub_read_NC_file( in_path, z_in_name, &
                                            z%din, &
                                  slat,      elat, &
                                  slon,      elon   )
        CALL Sub_read_NC_file( in_path, t_in_name, &
                                         temp%din, &
                                  slat,      elat, &
                                  slon,      elon   )
        CALL Sub_read_NC_file( in_path, q_in_name, &
                                            q%din, &
                                  slat,      elat, &
                                  slon,      elon   )
        CALL Sub_read_NC_file( in_path, w_in_name, &
                                            w%din, &
                                  slat,      elat, &
                                  slon,      elon   )


        CALL Sub_real_init                        &
             (                                    &
              read_pres, read_psfc, read_temp,    &
              z%dz, input_nz, nz, Ps, sfc_t,      &
              q%din, temp%din, z%din,    &
              q%dz,  temp%dz,  p%dz     &
             )
      ELSE IF ( incase == 3 ) THEN
        CALL Sub_read_txt_file                &
             (                                &
              in_path, input_name,             &
              input_nz,                        &
              sfc_p, sfc_t, sfc_q,             &
              z%din, temp%din, q%din           &
             )

        IF (ALLOCATED(w%din)) DEALLOCATE(w%din)
        ALLOCATE(w%din(input_nz))
        w%din = 0.


        CALL Sub_real_init                        &
             (                                    &
              read_pres, read_psfc, read_temp,    &
              z%dz, input_nz, nz, sfc_p, sfc_t,   &
              q%din, temp%din, z%din,      &
              q%dz,  temp%dz,  p%dz       &
             )
      ENDIF
      w%dz(:) = 0.5

      ! open(unit = 111, file = 'r.bin', status = "unknown", &
      !       access="direct", form="unformatted", recl=100*4) 
      ! write(unit = 111, rec=1) temp%dz
      ! write(*,*) "qwerqwerqwerqwerqwerqwer"
      !
      ! stop

      CALL Sub_set_W ( nz , dz%dz , w%dz , w%stag_dz )

      IF ( vertical_advection ) THEN
        w%stag_dz(:) = 0.5
      ELSE
        w%stag_dz(:) = 0.0
      ENDIF

      !! computed dt considering CFL condition 
      CALL Sub_set_dt
      CALL Sub_allocate_dt
      IF (incase == 1 .OR. incase == 2) THEN
        Temp%sfc_dt(:) = 10. + 273.5
           q%sfc_dt(:) = 0.
      ELSE IF (incase == 3) THEN
        Temp%sfc_dt(:) = sfc_t
           q%sfc_dt(:) = sfc_q
      ENDIF

      DO it = 1, nt   !! top condition
        Temp%top_dt(it) = Temp%dz(nz)
           q%top_dt(it) = 0.
      ENDDO 

      DO iz = 1, nz
        CALL Sub_drop_distributions        &
             (                             &
              dist_option,                 &
              drop_column_num,             &
              drop_min_diameter,           &
              drop_max_diameter,           &
              drop%num(iz,:),              &
              drop%r(iz,:), drop%rb(iz,:), &
              drop%dr,                     & 
              drop%m(iz,:), drop%mb(iz,:)  &
             ) 
      ENDDO
      IF (vertical_advection == .TRUE.) drop%num(2:nz,:) = 0.

      write(*,*) "!== checking initial distribution"
      SELECT CASE(dist_option)
        CASE(1)
          ! Log-normal distribution
          write(*,*) "gamma dist. in model (qc) =  ", sum(drop%m(1,:)*drop%num(1,:)*drop%dr(:))
          write(*,*) "gamma dist. const.   (qc) =  ", qc
          write(*,*) "log-normal dist. in model (nc) =  ", sum(drop%num(1,:)*drop%dr(:))
          write(*,*) "log-normal dist. const.   (nc) =  ", nc
          write(*,*) "accuracy  = ", (sum(drop%num(1,:)*drop%dr(:))/nc) * 100., "%"
        CASE(2)
          ! Gamma distribution
          write(*,*) "gamma dist. in model (qc) =  ", sum(drop%m(1,:)*drop%num(1,:)*drop%dr(:))
          write(*,*) "gamma dist. const.   (qc) =  ", qc
          write(*,*) "gamma dist. in model (nc) =  ", sum(drop%num(1,:)*drop%dr(:))
          write(*,*) "gamma dist. const.   (nc) =  ", nc
          write(*,*) "accuracy  = ", (sum(drop%num(1,:)*drop%dr(:))/nc) * 100., "%"
      END SELECT 
      write(*,*) "!==!==!==!=="
    END SUBROUTINE Sub_init_vars

    !!---------------------------------------------!!
    !!  Cal. dt regared CFL condition              !!
    !!---------------------------------------------!!
    SUBROUTINE Sub_set_dt

    IMPLICIT NONE

      ! Courant-Friedrichs-Lewy condition
      WHERE (w%dz /= 0.)
        CFL%dz = courant_number*(dz%dz/abs(w%dz))
      ELSEWHERE
        CFL%dz = MAXVAL(CFL%dz)
      END WHERE
      dt = REAL(INT(MINVAL(CFL%dz)))
      dt = 1. 
      nt = INT(integrated_time/dt)

      IF ( nt*dt .ne. integrated_time ) then
        write(*,*) "  "
        write(*,*) "********WARMING"
        write(*,*) "Calculated Total integrated time is different from the namelist integrated_time"
        write(*,*) "Total integrated time     =  ", nt*dt
        write(*,*) "Namelist integrated_time  =  ", integrated_time
        write(*,*) "********"
        write(*,*) "  "
      ENDIF

      write(*,*) "dt =  ", dt
      write(*,*) "nt =  ", nt
      write(*,*) "  "

    END SUBROUTINE Sub_set_dt

    !!---------------------------------------------!!
    !!  Cal. vertical coordinate                   !!
    !!---------------------------------------------!!
    SUBROUTINE Sub_set_grid

    IMPLICIT NONE

      ! Cal. dz      
      IF ( dz_option .eq. 1) THEN
        dz%dz(:) = z_top/nz
      ELSE IF ( dz_option .eq. 2) THEN
        dz%dz(1) = ((dzr-1)*z_top)/((dzr**nz)-1) 
        DO iz = 2, nz
          dz%dz(iz) = dz%dz(iz-1)*dzr
        ENDDO !! z
      ENDIF

      ! Cal. height     
      z%dz(1)= dz%dz(1)
      DO iz = 2, nz
        z%dz(iz)= z%dz(iz-1) + dz%dz(iz)
      ENDDO

    END SUBROUTINE Sub_set_grid
   
    !!---------------------------------------------!!
    !!  Cal. vertical wind                         !!
    !!---------------------------------------------!!
    SUBROUTINE Sub_set_W ( nz , dz , grid_w , stag_w )

      IMPLICIT NONE

      !In
      INTEGER,               INTENT(IN)   :: nz
      REAL, DIMENSION(nz),   INTENT(IN)   :: dz 
      REAL, DIMENSION(nz),   INTENT(IN)   :: grid_w
      !Local
      INTEGER                             :: iz
      REAL                                :: wgt 
      !Out
      REAL, DIMENSION(nz+1), INTENT(OUT)  :: stag_w

      ! Make stagged grid for advection
      DO iz = 2, nz
          wgt = dz(iz-1) / (dz(iz-1)+dz(iz))
          stag_w(iz) = grid_w(iz-1) + wgt*(grid_w(iz)-grid_w(iz-1))
      end do

      ! Set Boundary Condition
      ! most likely w = 0 at these points
      stag_w(1) = 0.; stag_w(nz+1) = 0.     ! Homogeneous Dirichlet BC

    END SUBROUTINE Sub_set_W  


    SUBROUTINE Sub_set_boundary
    END SUBROUTINE Sub_set_boundary  
  
END MODULE Mod_init_driver
