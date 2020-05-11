MODULE Mod_read  


  USE NETCDF
  USE Mod_global

  IMPLICIT NONE

  NAMELIST /Time_control/ integrated_time, &
                          output_interval 

  NAMELIST /var_info    / incase,    &
                          read_pres, &
                          read_temp, &
                          read_psfc

  NAMELIST /Domain      / input_nz,  &
                          nz,        &
                          slon,      &
                          elon,      &
                          slat,      &
                          elat,      &
                          read_pres, &
                          read_psfc, &
                          z_top

  NAMELIST /Dyn_options/ gamma_dry,     &

                         dyn_option,    &
                         dz_option,     &
                         dzr

  NAMELIST /Phys_options/ dist_option,       &
                          drop_column_num,   &
                          drop_min_diameter, &
                          drop_max_diameter

  NAMELIST /file_info/ in_path,     &
                       t_in_name,   &
                       q_in_name,   &
                       w_in_name,   &
                       z_in_name,   &
                       output_path, &
                       output_name


    CONTAINS

      !!---------------------------------------------!!
      !!---Sub_name : Sub_read_namelist      --------!!
      !!---------------------------------------------!!
    SUBROUTINE Sub_read_namelist

      IMPLICIT NONE

      OPEN(10,FILE='./namelist.info', iostat=ionum)

      IF ( ionum .ne. 0 ) CALL FAIL_MSG("error namelist")
      
      READ(10, Time_control)
      READ(10, var_info    )
      READ(10, Domain      )
      READ(10, Dyn_options )
      READ(10, Phys_options)
      READ(10, file_info   )

    END SUBROUTINE Sub_read_namelist

    !!---------------------------------------------!!
    !!---Sub_name :                        --------!!
    !!---------------------------------------------!!
    SUBROUTINE Sub_read_NC_file( inpath, inname, &
                                        out_var, &
                                  slat,    elat, &
                                  slon,    elon   )


      IMPLICIT NONE
      CHARACTER(LEN=*),                INTENT(IN)    :: inpath, inname
      INTEGER,                         INTENT(IN)    :: slat, elat, &
                                                        slon, elon     
                                                        
      REAL, ALLOCATABLE, DIMENSION(:), INTENT(INOUT) :: out_var

      !! for reading nc file
      INTEGER               :: ncid, i
      CHARACTER(LEN=50)     :: zname, xname, yname, tname, vname
      INTEGER               :: nz, nx, ny, nt
      INTEGER, DIMENSION(4) :: dim_count, dim_start 
      INTEGER, DIMENSION(4) :: dimids
      INTEGER               :: xtype,ndims,   &
                               ynum, xnum  

      REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: var_area
     
      IF (ALLOCATED(out_var)) DEALLOCATE(out_var)

      CALL check(nf90_open(TRIM(inpath)//TRIM(inname), nf90_nowrite, ncid))
      CALL SUCCESS_MSG("Open input")
      CALL check(nf90_inquire_dimension(ncid, 1, tname, nt))
      CALL check(nf90_inquire_dimension(ncid, 2, zname, nz))
      CALL check(nf90_inquire_dimension(ncid, 3, xname, nx))
      CALL check(nf90_inquire_dimension(ncid, 4, yname, ny))

      !!!!-----------------------------  var(ny, nx, nz, nt)
      IF ( elat .lt. slat )  ynum = slat-elat ; dim_start =(/ elat  , slon    , 1 , 8 /)
      IF ( elat .gt. slat )  ynum = elat-slat ; dim_start =(/ slat  , slon    , 1 , 8 /)
      xnum = elon - slon

      dim_count =(/ ynum  , xnum    , nz , 1 /)

      IF (.NOT. ALLOCATED(var_area)) ALLOCATE(var_area (ynum , xnum , nz, 1 ))
      IF (.NOT. ALLOCATED(out_var))  ALLOCATE(out_var (nz))

      ! dim_count =(/  1    ,   1    , nz, 1 /)
      ! dim_start =(/lat_iy , lon_ix , 1 , 1 /)

      CALL check(nf90_inquire_variable(ncid,1,vname,xtype,ndims,dimids))
      CALL check(nf90_inq_varid(ncid, vname ,varid))
      CALL check(nf90_get_var(ncid,varid,var_area,start=dim_start, count=dim_count))

      DO i = 1, nz
        out_var(i) = sum(var_area(:,:,nz+1-i,1))/(ynum*xnum)
      ENDDO
      ! out_var(:) = var(1,1,:,1)

    END SUBROUTINE Sub_read_NC_file 

END MODULE Mod_read  
