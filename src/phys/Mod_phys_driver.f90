MODULE Mod_phys_driver 

  USE Mod_global
  USE Mod_const
  USE Mod_drop_growth

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE Sub_phys_driver

  IMPLICIT NONE
  INTEGER :: i_m

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

      DO i_m = 1, SIZE(drop%mb(iz,:)) 
        IF (drop%mb(iz,i_m) < 0.) drop%mb(iz,i_m) = 0
        IF (i_m == SIZE(drop%mb(iz,:))) EXIT
        IF (drop%m(iz,i_m) < 0.) drop%m(iz,i_m) = 0
      ENDDO

      CALL compute_redist                           &
           (                                        &
            drop%ref_m(:),  drop%m(iz,:),           &
            drop%ref_mb(:), drop%mb(iz,:),          &
            drop%dm_dt(iz,:), drop%dmb_dt(iz,:),    &
            dt, drop%num(iz,:), drop%next_num(iz,:) &
           )

      drop%num(iz,:) = drop%next_num(iz,:)

      ! total_m(it) = SUM(drop%next_num(iz,:)*drop%m(iz,:))
      ! drop%dqv    = -1 * (total_m(it) - total_m(it-1))/rho
      ! drop%dTemp  = -1 * (L* drop%dqv)/(rho*Cp) 
      !
      ! temp%next_dz=temp%dz !- dl /dt
      ! q%next_dz=q%dz + drop%dqv
      ! print*, 'total mass :', SUM(drop%next_num(iz,:)*drop%m(iz,:))

      ! collision./com
      ! CALL compute_collision
    ENDDO
    ! write(*,*) " after "
    ! write(*,*) " num of bot layer at 1st time  =", sum(drop%num(1,:))
    ! write(*,*) " num of top layer at 1st time  =", sum(drop%num(100,:))
    ! write(*,*) "==============================================================="

  END SUBROUTINE Sub_phys_driver

END MODULE Mod_phys_driver  
