MODULE Mod_global

  USE NETCDF
  IMPLICIT NONE
   
  INTEGER :: it, iz   !! do parameter

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
                        drop_column_num

  REAL               :: z_top,              &
                        drop_ratio,         &
                        gamma_dry,          &
                        dzr,                &
                        drop_max_diameter,  &
                        drop_min_diameter

  INTEGER            :: aaa,bbb,ccc,ddd

  CHARACTER(LEN=256) ::     in_path, & 
                          t_in_name, &
                          q_in_name, &
                          w_in_name, &
                          z_in_name, &
                        output_path, &
                        output_name

    ! Declare variables 
  INTEGER                               :: nt, dt
  INTEGER                               :: varid
  TYPE varinfo
    INTEGER                             :: varid
    REAL, DIMENSION(:),     ALLOCATABLE :: dz, next_dz,      &
                                           stag_dz, dt,      &
                                           sfc_dt, top_dt,   &
                                           din, r, m, rb, mb, dmb, &
                                           num, next_num, dmb_dt
    REAL, DIMENSION(:,:),   ALLOCATABLE :: dout       
    CHARACTER(LEN=256)                  :: vname, axis,      &
                                           desc, units
  END TYPE varinfo

  TYPE(varinfo) ::    Temp,      & !! Temperature [K] 
                      th,        & !! Theta 
                      q,         & !! Number of water droplets
                      w,         & !! Vertical velocity
                      dz,        & !! Difference z
                      time,      & !! Difference z
                      drop,      & !! Difference z
                      p,         & !! Difference z
                      CFL,       & !! Difference z
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

    IF (.NOT. ALLOCATED(Temp%din     )) ALLOCATE(Temp%din    (input_nz))
    IF (.NOT. ALLOCATED(Temp%dz      )) ALLOCATE(Temp%dz           (nz))
    IF (.NOT. ALLOCATED(Temp%next_dz )) ALLOCATE(Temp%next_dz      (nz))

    IF (.NOT. ALLOCATED(q%din        )) ALLOCATE(q%din       (input_nz))
    IF (.NOT. ALLOCATED(q%dz         )) ALLOCATE(q%dz              (nz))
    IF (.NOT. ALLOCATED(q%next_dz    )) ALLOCATE(q%next_dz         (nz))

    IF (.NOT. ALLOCATED(w%din        )) ALLOCATE(w%din       (input_nz))
    IF (.NOT. ALLOCATED(w%dz         )) ALLOCATE(w%dz              (nz))
    IF (.NOT. ALLOCATED(w%stag_dz    )) ALLOCATE(w%stag_dz       (0:nz))

    IF (.NOT. ALLOCATED(z%din        )) ALLOCATE(z%din       (input_nz))
    IF (.NOT. ALLOCATED(z%dz         )) ALLOCATE(z%dz              (nz))

    IF (.NOT. ALLOCATED(dz%dz        )) ALLOCATE(dz%dz             (nz))
    IF (.NOT. ALLOCATED(p%dz         )) ALLOCATE(p%dz              (nz))

    IF (.NOT. ALLOCATED(CFL%dz       )) ALLOCATE(CFL%dz            (nz))

    IF (.NOT. ALLOCATED(drop%num     )) ALLOCATE(drop%num (drop_column_num))
    IF (.NOT. ALLOCATED(drop%next_num     )) ALLOCATE(drop%next_num (drop_column_num))
    IF (.NOT. ALLOCATED(drop%r       )) ALLOCATE(drop%r   (drop_column_num))
    IF (.NOT. ALLOCATED(drop%m       )) ALLOCATE(drop%m   (drop_column_num))
    IF (.NOT. ALLOCATED(drop%dmb       )) ALLOCATE(drop%dmb   (drop_column_num))

    IF (.NOT. ALLOCATED(drop%rb      )) ALLOCATE(drop%rb  (drop_column_num+1))
    IF (.NOT. ALLOCATED(drop%mb      )) ALLOCATE(drop%mb  (drop_column_num+1))
    IF (.NOT. ALLOCATED(drop%dmb_dt      )) ALLOCATE(drop%dmb_dt  (drop_column_num+1))

  END SUBROUTINE Sub_allocate_dz
 
  !!-----------------------------!!
  SUBROUTINE Sub_allocate_dt

    IF (.NOT. ALLOCATED(Temp%sfc_dt  )) ALLOCATE(Temp%sfc_dt       (nt))
    IF (.NOT. ALLOCATED(Temp%top_dt  )) ALLOCATE(Temp%top_dt       (nt))
    IF (.NOT. ALLOCATED(Temp%dout    )) ALLOCATE(Temp%dout    (nz,nt+1))

    IF (.NOT. ALLOCATED(q%sfc_dt     )) ALLOCATE(q%sfc_dt          (nt))
    IF (.NOT. ALLOCATED(q%top_dt     )) ALLOCATE(q%top_dt          (nt))
    IF (.NOT. ALLOCATED(q%dout       )) ALLOCATE(q%dout       (nz,nt+1))

    IF (.NOT. ALLOCATED(drop%dout    )) ALLOCATE(drop%dout (drop_column_num,nt))

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
