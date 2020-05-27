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

      ! temp%dout     = 0. 
      ! q%dout        = 0.

      ! defirne output 1st time variables.
      ! temp%dout(:,1)=temp%dz(:)
      temp%dout(:,1) = drop%num(50,:)
         q%dout(:,1) = q%dz(:)
      drop%ref_m     = drop%m(1,:) 
      drop%ref_mb    = drop%mb(1,:)

      DO it = 1, nt

        write(*,*) 'it times = ', it
        ! Compute dyn
        CALL Sub_dyn_driver          ! in  : temp%dz, q%dz
                                     ! out : temp%next_dz, q%next_dz

        ! update temp. and q by dyn.
        temp%dz = temp%next_dz
        q%dz    = q%next_dz
  
        write(*,*) 'it times = ', it
        ! CALL SUCCESS_MSG("dyn success") 
        ! write(*,*) 'check tmep. =', temp%dz

        ! Compute phys
        CALL Sub_phys_driver         ! in  : temp%dz, q%dz, drop%num, RH
                                     ! out : temp%next_dz, q%next_dz, drop%next_num, RH 

        write(*,*) 'it times = ', it
        ! CALL SUCCESS_MSG("phys success") 
        ! write(*,*) 'check drop. =', drop%num(50,:) 

        ! write(*,*) it
        ! write(*,*) drop%num(1,:)
        ! update temp. and q by dyn.
        temp%dz = temp%next_dz
        q%dz    = q%next_dz

        ! save output vars
        ! test temp.
        ! temp%dout(:,it+1)=temp%next_dz(:)
        temp%dout(:,it+1)=drop%next_num(50,:)
        q%dout(:,it+1)=q%next_dz(:)

      ENDDO !! time do

    END SUBROUTINE Sub_Integration_time

ENDMODULE
