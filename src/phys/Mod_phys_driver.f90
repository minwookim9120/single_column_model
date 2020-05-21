MODULE Mod_phys_driver 

  USE Mod_global
  USE Mod_const
  USE Mod_drop_growth

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE Sub_phys_driver

  IMPLICIT NONE

    write(*,*) "==============================================================="
    write(*,*) "for checking conservation after condensation and redistribution"
    write(*,*) " before "
    write(*,*) " num of bot layer at 1st time  =", sum(drop%num(1,:))
    write(*,*) " num of top layer at 1st time  =", sum(drop%num(100,:))
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

      CALL compute_redist                           &
           (                                        &
            drop%ref_m(:),  drop%m(iz,:),           &
            drop%ref_mb(:), drop%mb(iz,:),          &
            drop%dm_dt(iz,:), drop%dmb_dt(iz,:),    &
            dt, drop%num(iz,:), drop%next_num(iz,:) &
           )

      drop%num(iz,:) = drop%next_num(iz,:)

      ! drop%dqv    = -1 * (total_m(it) - total_m(it-1))
      ! drop%dTemp  = -1 * (L* drop%dqv)/(rho*Cp) 

      ! temp%next_dz=temp%dz !- dl/dt
      ! q%next_dz=q%dz !- dRH/dt
      ! collision./com
      ! CALL compute_collision
    ENDDO
    ! write(*,*) drop%num(5,:)
    write(*,*) " after "
    write(*,*) " num of bot layer at 1st time  =", sum(drop%num(1,:))
    write(*,*) " num of top layer at 1st time  =", sum(drop%num(100,:))
    write(*,*) "==============================================================="

  END SUBROUTINE Sub_phys_driver

END MODULE Mod_phys_driver  
