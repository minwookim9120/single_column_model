MODULE Mod_phys_driver 

  USE Mod_global
  USE Mod_const
  USE Mod_drop_growth

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE Sub_phys_driver

  IMPLICIT NONE
  INTEGER :: i_m
  REAL    :: dtotal_m, dqv, dtemp  

    ! write(*,*) "==============================================================="
    ! write(*,*) "for checking conservation after condensation and redistribution"
    ! write(*,*) " before "
    ! write(*,*) " num of bot layer at 1st time  =", sum(drop%num(1,:))
    ! write(*,*) " num of top layer at 1st time  =", sum(drop%num(100,:))
    DO iz = 1, nz 

      ! condensation and evaporation
      CALL compute_dmb_dt                       &
           (                                    &
            temp%dz(iz), q%dz(iz), p%dz(iz),    &
            drop%r(iz,:), drop%rb(iz,:),        &
            drop%dm_dt(iz,:), drop%dmb_dt(iz,:) &
           )

      drop%m(iz,:)  = drop%m(iz,:)  + drop%dm_dt(iz,:)*dt
      drop%mb(iz,:) = drop%mb(iz,:) + drop%dmb_dt(iz,:)*dt
      drop%r(iz,:)  = (3 * drop%m(iz,:) * 1. / (rho * 4 * pi))**(1./3.)
      drop%rb(iz,:) = (3 * drop%mb(iz,:) * 1. / (rho * 4 * pi))**(1./3.)

      DO i_m = 1, SIZE(drop%mb(iz,:)) 
        IF (drop%mb(iz,i_m) < 0.) drop%mb(iz,i_m) = 0
        IF (i_m == SIZE(drop%mb(iz,:))) EXIT
        IF (drop%m(iz,i_m) < 0.) drop%m(iz,i_m) = 0
      ENDDO
      !
      ! write(*,*) iz 
       !write(*,*) temp%dz(iz), q%dz(iz)
       !write(*,*) drop%dm_dt(iz,:)   

      CALL compute_redist                            &
           (                                         &
            drop%ref_m(:),  drop%m(iz,:),            &
            drop%ref_mb(:), drop%mb(iz,:),           &
            drop%dm_dt(iz,:), drop%dmb_dt(iz,:),     &
            dt, drop%ref_num(:), drop%next_num(iz,:) &
           )

      ! feedback by phys
      !dtemp = 0.
      !dqv   = 0.

!print*, '--------' 
!print*, drop%ref_m(:)
!print*, '--------' 
!print*, drop%m(iz,:)
!print*, '--------' 
!print*, iz 
!print*, '--------' 
!!print*, drop%ref_num(:)
!print*, drop%next_num(iz,:)

       dtotal_m = SUM(drop%next_num(iz,:)*drop%ref_m(:)) - &
                  SUM(drop%num(iz,:)*drop%ref_m(:))
      
       drop%num(iz,:) = drop%next_num(iz,:)
      
       dqv    = -1 * dtotal_m / rho      ! kg kg-1
       dTemp  = -1 * (L*dqv) / (rho*Cp) 

       temp%next_dz(iz)=temp%dz(iz) + dtemp
       q%next_dz(iz)=q%dz(iz) + dqv
!       temp%next_dz(iz)=temp%dz(iz)
!       q%next_dz(iz)=q%dz(iz) 

      ! collision
      ! CALL compute_collision
    ENDDO

    !drop%next_num(1,:) = drop%ref_num(:)
    !drop%drop%m(1,:) = drop%ref_m(:)

    ! write(*,*) " after "
    ! ! write(*,*) drop%dm_dt(1,:)
    ! write(*,*) drop%next_num(1,:)
    !
    ! write(*,*) " num of bot layer at 1st time  =", sum(drop%next_num(1,:))
    ! write(*,*) " num of top layer at 1st time  =", sum(drop%next_num(100,:))
    ! write(*,*) "==============================================================="

  END SUBROUTINE Sub_phys_driver

END MODULE Mod_phys_driver  
