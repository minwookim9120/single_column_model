MODULE Mod_phys_driver 

  USE Mod_global
  USE Mod_const
  USE Mod_drop_growth

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE Sub_phys_driver

  IMPLICIT NONE
  INTEGER :: i_m
  INTEGER :: iz
  REAL    :: dtotal_m, dqv, dtemp  

    ! write(*,*) "==============================================================="
    ! write(*,*) "for checking conservation after condensation and redistribution"
    ! write(*,*) " before "
    ! write(*,*) " num of bot layer at 1st time  =", sum(drop%num(1,:))
    ! write(*,*) " num of top layer at 1st time  =", sum(drop%num(100,:))
    DO iz = 1, nz 

      ! condensation and evaporation
      CALL compute_dmb_dt                        &
           (                                     &
            temp%dz(iz), q%dz(iz), p%dz(iz),     &
            drop%r(iz,:), drop%rb(iz,:), &
            drop%dm_dt(iz,:), drop%dmb_dt(iz,:)  &
           ) 

      drop%m(iz,:)  = drop%ref_m(:)  + drop%dm_dt(iz,:)*dt
      drop%mb(iz,:) = drop%ref_mb(:) + drop%dmb_dt(iz,:)*dt

      DO i_m = 1, SIZE(drop%mb(iz,:)) 
        IF (drop%mb(iz,i_m) < 0.) drop%mb(iz,i_m) = 0
        IF (i_m == SIZE(drop%mb(iz,:))) EXIT
        IF (drop%m(iz,i_m) < 0.) drop%m(iz,i_m) = 0
      ENDDO


      CALL compute_redist                            &
           (                                         &
            drop%ref_m(:),  drop%m(iz,:),            &
            drop%ref_mb(:), drop%mb(iz,:),           &
            drop%dm_dt(iz,:), drop%dmb_dt(iz,:),     &
            dt, drop%dr(:),                          &
            drop%num(iz,:), drop%next_num(iz,:)      &
           )


      IF ( phys_feedback ) THEN  ! feedback by phys
        dtotal_m = SUM(drop%next_num(iz,:)*drop%ref_m(:)*drop%dr(:)) - &
                   SUM(drop%num(iz,:)*drop%ref_m(:)*drop%dr(:))
        !! 2020.06.03 by han ====== start
        dqv    = -1 * dtotal_m            ! kg kg-1
        !! 2020.06.03 by han ====== end   
        dTemp  = -1 * (L*dqv) / Cp 
      ELSE ! non - feedback by phys
        dtemp = 0.
        dqv   = 0.
      ENDIF

      drop%num(iz,:) = drop%next_num(iz,:)

      temp%next_dz(iz)=temp%dz(iz) + dtemp
      q%next_dz(iz)=q%dz(iz) + dqv

      ! collision
      ! CALL compute_collision
    ENDDO


    !print*, 'qv [kg/kg]            qc[kg/kg]         qv+qc [kg/kg]'
    !print*, q%next_dz(10), SUM(drop%next_num(10,:)*drop%ref_m(:)*drop%dr(:)), &
    !        q%next_dz(10) + SUM(drop%next_num(10,:)*drop%ref_m(:)*drop%dr(:))
    ! write(*,*) " after "
    ! ! write(*,*) drop%dm_dt(1,:)
    ! write(*,*) maxval(drop%next_num(10,:))
    !
    !write(*,*) " num of bot layer at 1st time  =", sum(drop%next_num(10,:))
    ! write(*,*) " num of top layer at 1st time  =", sum(drop%next_num(100,:))
    ! write(*,*) "==============================================================="

  END SUBROUTINE Sub_phys_driver

END MODULE Mod_phys_driver  
