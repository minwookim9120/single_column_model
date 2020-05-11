MODULE Mod_integration

  USE Mod_global
  USE Mod_dyn_driver

  IMPLICIT NONE

    CONTAINS

    !!---------------------------------------------!!
    !!---------------------------------------------!!
    !!---------------------------------------------!!
    SUBROUTINE Sub_Integration_FV

      IMPLICIT NONE

      ! temp%dout     = 0. 
      ! q%dout        = 0.

      temp%dout(:,1)=temp%dz(:)
      q%dout(:,1)=q%dz(:)

      DO it = 1, nt

        CALL Sub_Finite_volume ( temp%dz, temp%sfc_dt(it),  &
                                 temp%top_dt(it),           &
                                 dz%dz, nz, CFL%dz,         &
                                 dt,                        &
                                 w%stag_dz,                 &
                                 temp%next_dz               &
                                                           )

        CALL Sub_Finite_volume ( q%dz, q%sfc_dt(it),        &
                                 q%top_dt(it),              &
                                 dz%dz, nz, CFL%dz,         &
                                 dt,                        &
                                 w%stag_dz,                 &
                                 q%next_dz                  & 
                                                           )

        IF (ALLOCATED(temp%dz)) DEALLOCATE(temp%dz)
        IF (ALLOCATED(q%dz   )) DEALLOCATE(q%dz   )
        IF (.NOT. ALLOCATED(temp%dz)) ALLOCATE(temp%dz(nz))
        IF (.NOT. ALLOCATED(q%dz   )) ALLOCATE(q%dz   (nz))

        temp%dz(:)=temp%next_dz(:)
        q%dz(:)=q%next_dz(:)
 
        temp%dout(:,it+1)=temp%next_dz(:)
        q%dout(:,it+1)=q%next_dz(:)

        ! write(*,*) "variable"
        ! write(*,*) temp%dz
        ! IF ( it == 2 )  stop
        !!CALL cloud_pysics
      ENDDO !! time do

    END SUBROUTINE Sub_Integration_FV

    !!---------------------------------------------!!
    !!---------------------------------------------!!
    !!---------------------------------------------!!
    SUBROUTINE Sub_Integration_FD

      IMPLICIT NONE

      temp%dout(:,1)=temp%dz(:)
      q%dout(:,1)=q%dz(:)

      DO it = 1, nt

        CALL Sub_Finite_diff ( temp%dz, temp%sfc_dt(it),    &
                                 temp%top_dt(it),           &
                                 dz%dz, nz,                 &
                                 dt,                        &
                                 w%stag_dz,                 &
                                 temp%next_dz               &
                                                           )
  
        CALL Sub_Finite_diff ( q%dz, q%sfc_dt(it),          &
                                 q%top_dt(it),              &
                                 dz%dz, nz,                 &
                                 dt,                        &
                                 w%stag_dz,                 &
                                 q%next_dz                  & 
                                                           )
        IF (ALLOCATED(temp%dz)) DEALLOCATE(temp%dz)
        IF (ALLOCATED(q%dz   )) DEALLOCATE(q%dz   )
        IF (.NOT. ALLOCATED(temp%dz)) ALLOCATE(temp%dz(nz))
        IF (.NOT. ALLOCATED(q%dz   )) ALLOCATE(q%dz   (nz))

        temp%dz(:)=temp%next_dz(:)
        q%dz(:)=q%next_dz(:)

        temp%dout(:,it+1)=temp%next_dz(:)
        q%dout(:,it+1)=q%next_dz(:)

        ! CALL Sub_Cal_P
        ! CALL Sub_Cal_W
        !!CALL cloud_pysics
      ENDDO !! time do

    END SUBROUTINE Sub_Integration_FD

    !!---------------------------------------------!!
    !!---------------------------------------------!!
    !!---------------------------------------------!!
    SUBROUTINE Sub_Integration_PPM

      IMPLICIT NONE

      ! temp%dout     = 0. 
      ! q%dout        = 0.

      temp%dout(:,1)=temp%dz(:)
      q%dout(:,1)=q%dz(:)

      DO it = 1, nt

        CALL Sub_Finite_volume_PPM ( temp%dz, temp%sfc_dt(it),  &
                                     temp%top_dt(it),           &
                                     dz%dz, nz, CFL%dz,         &
                                     dt,                        &
                                     w%stag_dz,                 &
                                     temp%next_dz               &
                                                               )
        CALL Sub_Finite_volume_PPM ( q%dz, q%sfc_dt(it),        &
                                     q%top_dt(it),              &
                                     dz%dz, nz, CFL%dz,         &
                                     dt,                        &
                                     w%stag_dz,                 &
                                     q%next_dz                  &
                                                               )
        IF (ALLOCATED(temp%dz)) DEALLOCATE(temp%dz)
        IF (ALLOCATED(q%dz   )) DEALLOCATE(q%dz   )
        IF (.NOT. ALLOCATED(temp%dz)) ALLOCATE(temp%dz(nz))
        IF (.NOT. ALLOCATED(q%dz   )) ALLOCATE(q%dz   (nz))
        temp%dz(:)=temp%next_dz(:)
        q%dz(:)=q%next_dz(:)

        temp%dout(:,it+1)=temp%next_dz(:)
        q%dout(:,it+1)=q%next_dz(:)

        ! CALL Sub_Cal_P
        ! CALL Sub_Cal_W
        !!CALL cloud_pysics
      ENDDO !! time do

    END SUBROUTINE Sub_Integration_PPM

    SUBROUTINE Sub_Integration_PPM_dist

      IMPLICIT NONE

      ! temp%dout     = 0. 
      ! q%dout        = 0.

      drop%dout(:,1)=drop%num(:)
      drop%dout(:,1)=drop%num(:)

      do i = 1, drop_column_num
      drop%dmb = drop%mb(i+1)-drop%mb(i)
      enddo

      do it  = 0, nt-1
        drop%dmb_dt(it) = 4.*pi*rb(:)*1. / (((Rv*T)/(Dv*es)) + (((L/(Rv*T))-1)*(L/(Ka*T)))) *S*dt  !  ; dm@units = "m-2 s-1"
      end do 

      DO it = 1, nt

        CALL Sub_Finite_volume_PPM ( drop%num, q%sfc_dt(it),  &
                                     q%top_dt(it),           &
                                     drop%dmb, drop_column_num, CFL%dz,         &
                                     1,                        &
                                     w%stag_dz,                 &
                                     drop%next_num               &
                                                               )
        temp%dz(:)=temp%next_dz(:)
        q%dz(:)=q%next_dz(:)

        temp%dout(:,it+1)=temp%next_dz(:)
        q%dout(:,it+1)=q%next_dz(:)

        ! CALL Sub_Cal_P
        ! CALL Sub_Cal_W
        !!CALL cloud_pysics
      ENDDO !! time do

    END SUBROUTINE Sub_Integration_PPM_dist

ENDMODULE
