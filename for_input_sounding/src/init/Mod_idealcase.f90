MODULE Mod_idealcase

  USE Mod_global

  IMPLICIT NONE

    CONTAINS
    
    !!---------------------------------------------!!
    !!  Cal. Temperature                           !!
    !!---------------------------------------------!!
    SUBROUTINE Sub_ideal_init      &
              (                    &
               nz, dz, gamma_dry,  &
               temp_sfc,           &
               temp, q             & 
              )

      IMPLICIT NONE
      ! In
      INTEGER,             INTENT(IN)  :: nz
      REAL,                INTENT(IN)  :: temp_sfc, gamma_dry
      REAL, DIMENSION(nz), INTENT(IN)  :: dz
      ! Out
      REAL, DIMENSION(nz), INTENT(OUT) :: temp, q
      ! Local
      INTEGER                          :: iz

      !! initial  condition
      DO iz = 1, nz
        IF (iz .EQ. 1) THEN
          Temp(iz) = Temp_sfc - gamma_dry * dz(iz)
        ELSE
          Temp(iz) = Temp(iz-1) - gamma_dry * dz(iz)
        ENDIF 
        q(iz) = 100*sin(real(iz)) + 100.
      ENDDO 

    END SUBROUTINE Sub_ideal_init

END MODULE Mod_idealcase
