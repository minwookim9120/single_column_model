MODULE Mod_integration

  USE Mod_global
  USE Mod_const
  USE Mod_dyn_driver
  USE Mod_drop_growth
  USE Mod_phys_driver

  IMPLICIT NONE

    CONTAINS

    !!---------------------------------------------!!
    !!---------------------------------------------!!
    !!---------------------------------------------!!
    SUBROUTINE Sub_Integration_time

      IMPLICIT NONE
      INTEGER :: iar, iaz
      ! temp%dout     = 0. 
      ! q%dout        = 0.

      ! defirne output 1st time variables.
      temp%dout(:,1)        = temp%dz(:)
      q%dout(:,1)           = q%dz(:)
      DO iar = 1, drop_column_num
      DO iaz = 1, nz
        drop%drop_dout(iar,iaz,1)  = drop%num(iaz,iar)*drop%dr(iar)
      ENDDO
      ENDDO
      ! ref_m
      drop%ref_m     = drop%m(1,:) 
      drop%ref_mb    = drop%mb(1,:)
      drop%ref_num   = drop%num(1,:)

      DO it = 1, nt
        print*, "============================" 
        print*, it 
        !print*, temp%dz
        print*, "============================" 
        ! write(*,*)  it
        ! Compute dyn
        CALL Sub_dyn_driver          ! in  : temp%dz, q%dz
                                     ! out : temp%next_dz, q%next_dz

        ! update temp. and q by dyn.
        temp%dz     = temp%next_dz
        q%dz        = q%next_dz
        drop%num    = drop%next_num
  
        drop%num(1,:) = drop%ref_num(:)
        drop%m  (1,:) = drop%ref_m(:)

        ! Compute phys
        CALL Sub_phys_driver         ! in  : temp%dz, q%dz, drop%num, RH
                                     ! out : temp%next_dz, q%next_dz, drop%next_num, RH 
        drop%num(1,:) = drop%ref_num(:)
        drop%m  (1,:) = drop%ref_m(:)

        ! update temp. and q by phys.
        temp%dz = temp%next_dz
        q%dz    = q%next_dz

        ! save output vars
        temp%dout(:,it+1)         = temp%dz(:)
        q%dout(:,it+1)            = q%dz(:) 
        DO iar = 1, drop_column_num
        DO iaz = 1, nz
          drop%drop_dout(iar,iaz,it+1)  = drop%num(iaz,iar)*drop%dr(iar)
        ENDDO
        ENDDO
      ENDDO !! time do

    END SUBROUTINE Sub_Integration_time

ENDMODULE
