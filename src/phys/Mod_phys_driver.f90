MODULE Mod_phys_driver 

  USE Mod_global
  USE Mod_const
  USE Mod_drop_growth
  USE Mod_collision

  IMPLICIT NONE

  CONTAINS

  SUBROUTINE Sub_phys_driver

  IMPLICIT NONE
  INTEGER :: i_m, ittt
  INTEGER :: iz
  REAL    :: dtotal_m, dqv, dtemp, sub_dt, sub_nt  
  REAL,dimension(nz)    :: qc 
    ! sub_dt=0.01
    !
    ! sub_nt=dt/sub_dt
    DO iz = 1, nz 
      ! do ittt = 1, sub_nt     
      ! condensation and evaporation
      CALL compute_dmb_dt                        &
           (                                     &
            temp%dz(iz), q%dz(iz), p%dz(iz),     &
            drop%r(iz,:), drop%rb(iz,:), &
            drop%dm_dt(iz,:), drop%dmb_dt(iz,:), drop%S(iz)  &
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
        dqv    = -1 * dtotal_m            ! kg kg-1
        dTemp  = -1 * (L*dqv) / Cp 
      ELSE ! non - feedback by phys
        dtemp = 0.
        dqv   = 0.
      ENDIF

      temp%next_dz(iz)=temp%dz(iz) + dtemp
      q%next_dz(iz)=q%dz(iz) + dqv

      ! temp%dz(iz)=temp%next_dz(iz)
      ! q%dz(iz)=q%next_dz(iz)
      drop%num(iz,:) = drop%next_num(iz,:)

      ! collision
      ! Change unit for collision
      qc(iz) = sum(drop%num(:,iz)*drop%m(:,iz)) * rho ! [kg kg-1] -> [kg m-3]

      ! Compute Stochastic Collision Equation
      if ( collision_effect )  then
          call coad1d( dt, drop_column_num, r0, qc(iz), drop%num(iz,:) )
      end if
    ! CALL compute_collision
    ENDDO
     ! write(*,*) drop%dmb_dt(2,:)

  END SUBROUTINE Sub_phys_driver

END MODULE Mod_phys_driver  
