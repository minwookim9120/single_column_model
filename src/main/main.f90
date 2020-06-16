PROGRAM main_prog

  USE NETCDF
  USE Mod_global
  USE Mod_const
  USE Mod_read
  USE Mod_init_driver
  USE Mod_dyn_driver    
  USE Mod_phys_driver    
  USE Mod_integration
  USE Mod_write
  
  IMPLICIT NONE
  


!== initialization
  CALL Sub_read_namelist
  CALL Sub_allocate_dz
  CALL Sub_set_grid

  CALL Sub_init_vars

!== do time
  CALL Sub_Integration_time

!== write output

   write(*,*) "it =            2", " total Q =",  sum(q%dout(:,2)*dz%dz(:))
   write(*,*) "it = ",         nt,  "total Q =", sum(q%dout(:,nt)*dz%dz(:))
   write(*,*) "it =            2", " total t =",  sum(temp%dout(:,2)*dz%dz(:))
   write(*,*) "it = ",         nt,  "total t =", sum(temp%dout(:,nt)*dz%dz(:))

   ! ddmm = drop%mb(1,2:drop_column_num+1)-drop%mb(1,1:drop_column_num)
   !
   ! write(*,*) "it =            2", " total t =",  sum(drop%drop_dout(:,2,1)*ddmm(:))
   ! write(*,*) "it = ",         nt,  "total t =", sum(drop%drop_dout(:,2,nt)*ddmm(:))
  CALL Sub_write_netcdf ( nz, nt, dz%dz, z%dz,      &
                          temp%dout, q%dout,        &
                          w%dz(1:nz),               &
                          output_path, output_name )     

  CALL Sub_write_netcdf_drop ( nz, nt, z%dz,                 &
                               drop_column_num, drop%r(1,:), drop%dr, dz%dz, &
                               drop%drop_dout, drop%ref_m, drop%S_out,   &
                               output_path, drop_name )    


  ! CALL Sub_deallocate 

END PROGRAM main_prog
