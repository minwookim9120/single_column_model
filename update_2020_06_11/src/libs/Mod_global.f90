MODULE Mod_global

  USE NETCDF
  IMPLICIT NONE
   
  INTEGER :: it, iz, itt   !! do parameter

    !! for namelist val
  INTEGER            :: integrated_time,    &
                        nz,                 &
                        input_nz,           &
                        ionum,              &
                        output_interval

  INTEGER            :: read_pres,     &
                        read_temp,     &        
                        read_psfc,     &        
                        incase

  INTEGER            :: slon,     &
                        elon,     &
                        slat,     &
                        elat

  INTEGER            :: dyn_option,         &
                        dz_option,          &
                        dist_option,        &
                        drop_column_num,    &
                        redist_option

  REAL               :: z_top,              &
                        drop_ratio,         &
                        gamma_dry,          &
                        dzr,                &
                        drop_max_diameter,  &
                        drop_min_diameter,  &
                        nc,                 &
                        qc,                 &
                        r0

  INTEGER            :: aaa,bbb,ccc,ddd

  CHARACTER(LEN=256) ::     in_path, & 
                          t_in_name, &
                          q_in_name, &
                          w_in_name, &
                          z_in_name, &
                         input_name, &
                        output_path, &
                          drop_name, &
                        output_name, &
                        phys_unit

 LOGICAL             :: ventilation_effect, &
                             phys_feedback, & 
                        vertical_advection, &
                         collision_effect

    ! Declare variables 
  REAL                                  :: dt, substep_dt
  INTEGER                               :: nt     
  INTEGER                               :: substep_nt
  INTEGER                               :: varid
  REAL                                  :: sfc_p,          &
                                           sfc_t,          &
                                           sfc_q


  TYPE varinfo
    INTEGER                             :: varid
    REAL, DIMENSION(:),     ALLOCATABLE :: dz,            & ! vars in delta z                   (nz)  
                                           next_dz,       & ! next vars in delta z              (nz)
                                           stag_dz,       & ! stagged vars                      (nz + 1)
                                           dt,            & ! vars in time array                (nt)
                                           sfc_dt,        & ! boundary value in time array      (nt)
                                           top_dt,        & ! boundary value in time array      (nt)
                                           ref_m,         & ! ref. mass range for bin           (drop_column_num)
                                           ref_mb,        & ! ref. boundary mass range for bin  (drop_column_num)
                                           ref_num,       & ! ref. number conc. for bin         (drop_column_num)
                                           dr,            & ! dr = rb(i+1) - rb(i)              (drop_column_num)
                                           din              ! var from input data or prescribed (input_nz)
    REAL, DIMENSION(:,:),   ALLOCATABLE :: dout,          & ! output                            (nz, nt)
                                           r,             & ! drop radius                       (nz, drop_column_num)
                                           m,             & ! drop mass                         (nz, drop_column_num)
                                           rb,            & ! boundary drop radius              (nz, drop_column_num + 1)
                                           mb,            & ! boundary drop mass                (nz, drop_column_num + 1)
                                           num,           & ! bin array                         (nz, drop_column_num)
                                           next_num,      & ! next var in bin array             (nz, drop_column_num)
                                           dm_dt,         & ! dm/dt                             (nz, drop_column_num)
                                           dmb,           & ! dmb                               (nz, drop_column_num)
                                           dmb_dt           ! dmb/dt                            (nz, drop_column_num + 1)
    REAL, DIMENSION(:,:,:), ALLOCATABLE :: drop_dout,     & ! output                            (nz, drop_column_num, nt) 
                                           mass_dout,     &
                                              r_dout
    CHARACTER(LEN=256)                  :: vname, axis,   &
                                           desc, units
  END TYPE varinfo

  TYPE(varinfo) ::    Time,      & !! for nc write
                      Temp,      & !! Temperature [K] 
                      th,        & !! Theta 
                      q,         & !! Specific humidity
                      w,         & !! Vertical velocity
                      dz,        & !! Difference z
                      drop,      & !! droplet 
                      p,         & !! pressure
                      CFL,       & !! CFL
                      z            !! Height 

    ! For nc file 
  INTEGER                       :: ncid,                       &
                                   rec_dimid, lev_dimid,       &
                                   lat_dimid, lon_dimid
  

  INTEGER,            PARAMETER :: dim1     = 1,               &
                                   dim4     = 4
  
  INTEGER, DIMENSION(dim1)      :: dimid1
  INTEGER, DIMENSION(dim4)      :: dimid4
  INTEGER, DIMENSION(dim1)      :: dim1_start, dim1_count
  INTEGER, DIMENSION(dim4)      :: dim4_start, dim4_count
  
  CHARACTER(LEN=256), PARAMETER :: des      = "description"
  CHARACTER(LEN=256), PARAMETER :: un       = "units"
  CHARACTER(LEN=256), PARAMETER :: ax       = "axis"

  CONTAINS

  !!-----------------------------!!
  SUBROUTINE Sub_allocate_dz

    IF (.NOT. ALLOCATED(Temp%din      )) ALLOCATE(Temp%din    (input_nz))
    IF (.NOT. ALLOCATED(Temp%dz       )) ALLOCATE(Temp%dz           (nz))
    IF (.NOT. ALLOCATED(Temp%next_dz  )) ALLOCATE(Temp%next_dz      (nz))

    IF (.NOT. ALLOCATED(q%din         )) ALLOCATE(q%din       (input_nz))
    IF (.NOT. ALLOCATED(q%dz          )) ALLOCATE(q%dz              (nz))
    IF (.NOT. ALLOCATED(q%next_dz     )) ALLOCATE(q%next_dz         (nz))

    IF (.NOT. ALLOCATED(w%din         )) ALLOCATE(w%din       (input_nz))
    IF (.NOT. ALLOCATED(w%dz          )) ALLOCATE(w%dz              (nz))
    IF (.NOT. ALLOCATED(w%stag_dz     )) ALLOCATE(w%stag_dz       (0:nz))

    IF (.NOT. ALLOCATED(z%din         )) ALLOCATE(z%din       (input_nz))
    IF (.NOT. ALLOCATED(z%dz          )) ALLOCATE(z%dz              (nz))

    IF (.NOT. ALLOCATED(dz%dz         )) ALLOCATE(dz%dz             (nz))
    IF (.NOT. ALLOCATED(p%dz          )) ALLOCATE(p%dz              (nz))

    IF (.NOT. ALLOCATED(CFL%dz        )) ALLOCATE(CFL%dz            (nz))
    
    IF (.NOT. ALLOCATED(drop%num      )) ALLOCATE(drop%num      (nz, drop_column_num))
    IF (.NOT. ALLOCATED(drop%next_num )) ALLOCATE(drop%next_num (nz, drop_column_num))
    IF (.NOT. ALLOCATED(drop%ref_num  )) ALLOCATE(drop%ref_num  (    drop_column_num))

    IF (.NOT. ALLOCATED(drop%r        )) ALLOCATE(drop%r        (nz, drop_column_num))
    IF (.NOT. ALLOCATED(drop%m        )) ALLOCATE(drop%m        (nz, drop_column_num))
    IF (.NOT. ALLOCATED(drop%dr       )) ALLOCATE(drop%dr       (    drop_column_num))
    IF (.NOT. ALLOCATED(drop%ref_m    )) ALLOCATE(drop%ref_m    (    drop_column_num))
    IF (.NOT. ALLOCATED(drop%ref_mb   )) ALLOCATE(drop%ref_mb   (    drop_column_num))
    IF (.NOT. ALLOCATED(drop%dm_dt    )) ALLOCATE(drop%dm_dt    (nz, drop_column_num))

    IF (.NOT. ALLOCATED(drop%rb       )) ALLOCATE(drop%rb      (nz, drop_column_num+1))
    IF (.NOT. ALLOCATED(drop%mb       )) ALLOCATE(drop%mb      (nz, drop_column_num+1))
    IF (.NOT. ALLOCATED(drop%dmb      )) ALLOCATE(drop%dmb     (nz, drop_column_num+1))
    IF (.NOT. ALLOCATED(drop%dmb_dt   )) ALLOCATE(drop%dmb_dt  (nz, drop_column_num+1))

  END SUBROUTINE Sub_allocate_dz
 
  !!-----------------------------!!
  SUBROUTINE Sub_allocate_dt

    IF (.NOT. ALLOCATED(Temp%sfc_dt  )) ALLOCATE(Temp%sfc_dt       (nt))
    IF (.NOT. ALLOCATED(Temp%top_dt  )) ALLOCATE(Temp%top_dt       (nt))
    IF (.NOT. ALLOCATED(Temp%dout    )) ALLOCATE(Temp%dout    (nz,nt+1))

    IF (.NOT. ALLOCATED(q%sfc_dt     )) ALLOCATE(q%sfc_dt          (nt))
    IF (.NOT. ALLOCATED(q%top_dt     )) ALLOCATE(q%top_dt          (nt))
    IF (.NOT. ALLOCATED(q%dout       )) ALLOCATE(q%dout       (nz,nt+1))

    IF (.NOT. ALLOCATED(drop%drop_dout )) ALLOCATE(drop%drop_dout (drop_column_num,nz,nt+1))
    IF (.NOT. ALLOCATED(drop%mass_dout )) ALLOCATE(drop%mass_dout (drop_column_num,nz,nt+1))
    IF (.NOT. ALLOCATED(drop%r_dout    )) ALLOCATE(drop%r_dout    (drop_column_num,nz,nt+1))

  END SUBROUTINE Sub_allocate_dt

  !!-----------------------------!!
  SUBROUTINE Sub_deallocate

    DEALLOCATE(Temp%din      )
    DEALLOCATE(Temp%dz       )
    DEALLOCATE(Temp%next_dz  )

    DEALLOCATE(q%din         )
    DEALLOCATE(q%dz          )
    DEALLOCATE(q%next_dz     )

    DEALLOCATE(w%din         )
    DEALLOCATE(w%dz          )
    DEALLOCATE(w%stag_dz     )

    DEALLOCATE(z%din         )
    DEALLOCATE(z%dz          )

    DEALLOCATE(dz%dz         )

    DEALLOCATE(Temp%sfc_dt   )
    DEALLOCATE(Temp%top_dt   )
    DEALLOCATE(Temp%dout     )
    DEALLOCATE(q%sfc_dt      )
    DEALLOCATE(q%top_dt      )
    DEALLOCATE(q%dout        )

  END SUBROUTINE Sub_deallocate

  !!-----------------------------!!
  SUBROUTINE Sub_nc_attri

    ! Set name of variables.
    Temp%vname     = "T"
    q%vname        = "Q"
    w%vname        = "W"
    z%vname        = "Lev"
    time%vname     = "Time"

    ! Set "Description" attributes.
    Temp%desc      = "Temperature"
    q%desc         = "mass of water droplets"
    w%desc         = "Vertical velocity"
    z%desc         = "Height"

    ! Set "units" attributes.
    Temp%units     = "K"
    q%units        = "kg kg-1"
    w%units        = "m s-1"
    z%units        = "m"
    time%units     = "minutes since 2000-01-01 00:00:00"

    ! Set "axis" attributes.
    z%axis         = "Z"
    time%axis      = "T"

  END SUBROUTINE Sub_nc_attri

  !!-----------------------------!!
  SUBROUTINE CHECK(status)

    IMPLICIT NONE

    INTEGER, INTENT(IN) :: status

    !Check errors.
    IF (status .ne. nf90_noerr) THEN
     PRINT *, trim(nf90_strerror(status))
     PRINT *, "    ERROR :: CHECK NC CODE       "
     STOP "##### ERROR: PROGRAM ABORTED. #####"
    END IF

  END SUBROUTINE CHECK

  !!-----------------------------!!
  SUBROUTINE FAIL_MSG(f_msg)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: f_msg

    WRITE (*,'("FAIL: ", a)') f_msg
    STOP "##### ERROR: PROGRAM ABORTED. #####"

  END SUBROUTINE FAIL_MSG

  !!-----------------------------!!
  SUBROUTINE SUCCESS_MSG(s_msg)

    IMPLICIT NONE
    CHARACTER(LEN=*), INTENT(IN) :: s_msg
    
    WRITE (*,'("SUCCESS: ", a)') s_msg

  END SUBROUTINE SUCCESS_MSG 


END MODULE Mod_global
