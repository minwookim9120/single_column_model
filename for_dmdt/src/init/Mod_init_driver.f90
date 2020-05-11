MODULE Mod_init_driver

  USE Mod_global
  USE Mod_const
  USE Mod_idealcase
  USE Mod_realcase
  USE Mod_distribution 

  IMPLICIT NONE

    CONTAINS
    
    SUBROUTINE Sub_init_vars
    
    IMPLICIT NONE

      !! prescribed vertical-wind

      !!  
      IF ( incase == 2 ) THEN
        CALL Sub_ideal_init         &
            (                       &
             nz, dz%dz, gamma_dry,  &
             283.15,                &
             temp%dz, q%dz          &
            )  
      ELSE IF ( incase == 1 ) THEN
        CALL Sub_real_init                        &
             (                                    &
              read_pres, read_psfc, read_temp,    &
              z%dz, input_nz, nz,                 &
              q%din, temp%din, z%din,  w%din,     &
              q%dz,  temp%dz,  p%dz,    w%dz     &
             )
      ENDIF
      ! w%dz(:) = 0.5 

      CALL Sub_set_W ( nz , dz%dz , w%dz , w%stag_dz )

        write(*,*) w%dz
      ! w%stag_dz(0) = 0.9 

      !! computed dt considering CFL condition 
      CALL Sub_set_dt

      nt = 50
      dt = 72
      CALL Sub_allocate_dt

   
      !! surface condition
      DO it = 1, nt
        Temp%sfc_dt(it) = 10. + 273.5
      ENDDO !! time do
      !! top condition
      DO it = 1, nt
        Temp%top_dt(it) = Temp%dz(nz)
      ENDDO !! time do

      q%sfc_dt(:)=0.
      q%top_dt(:)=0.

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
      dt = INT(MINVAL(CFL%dz))
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
