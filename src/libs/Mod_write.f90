  MODULE Mod_write

    USE NETCDF
    USE Mod_global

    IMPLICIT NONE
    
    CONTAINS

      !!---------------------------------------------!!
      SUBROUTINE Sub_write_netcdf ( nz, nt, dz, z_out,   &!{{{
                                    temp_out, q_out,     &
                                    w_out,               &
                                    out_path, out_name  )

        USE NETCDF
        IMPLICIT NONE
        ! In
        INTEGER,                INTENT(IN)   ::    nz, nt
        REAL, DIMENSION(nz),    INTENT(IN)   ::        dz,  &
                                                    z_out,  &
                                                    w_out
        REAL, DIMENSION(nz,nt+1), INTENT(IN) ::  temp_out,  &
                                                    q_out
        CHARACTER(LEN=*),       INTENT(IN)   ::  out_path,  &
                                                 out_name
        ! Local
        INTEGER                              ::      irec
        INTEGER                              ::       nnt,  &
                                                        x,  &
                                                        y
        INTEGER, DIMENSION(nt+1)             ::      date

        nnt=nt+1
        x=1;y=1

        DO it = 1, nnt
          date(it) = it
        ENDDO

        CALL Sub_nc_attri 

        ! Create the file.
        CALL CHECK( nf90_create(trim(out_path)//"/"//trim(out_name),nf90_clobber,ncid) )
        CALL SUCCESS_MSG("Creation of new output data")

        ! Define the dimensions. The time dimension has no limit.
        CALL check( nf90_def_dim(ncid,      "Lon",              x, lon_dimid) )
        CALL check( nf90_def_dim(ncid,      "Lat",              y, lat_dimid) )
        CALL check( nf90_def_dim(ncid,    z%vname,             nz, lev_dimid) )
        CALL check( nf90_def_dim(ncid, time%vname, nf90_unlimited, rec_dimid) )
  
        ! The dimids array is used to pass the dimids of the dimension of
        ! the netCDF variables. In Fortran, the unlimited dimension
        ! must come last on the list of dimids.
        dimid4 = (/ lon_dimid, lat_dimid, lev_dimid, rec_dimid /)
        dimid1 = (/ lev_dimid /)

        ! Define "x","y" and its attributes 
        CALL check( nf90_def_var(ncid, "Lon", nf90_double, lon_dimid, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 1,  ax,  "X") )
        CALL check( nf90_def_var(ncid, "Lat", nf90_double, lat_dimid, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 2,  ax,  "Y") )

        ! Define "lev" and its attributes 
        CALL check( nf90_def_var(ncid, z%vname, nf90_double, lev_dimid, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 3, des,  z%desc) )
        CALL CHECK( nf90_put_att(ncid, 3,  un, z%units) )
        CALL CHECK( nf90_put_att(ncid, 3,  ax,  z%axis) )

        ! Define "time" and its attributes 
        CALL check( nf90_def_var(NCID, time%vname, nf90_double, rec_dimid, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 4,  un, time%units) )
        CALL CHECK( nf90_put_att(ncid, 4,  ax,  time%axis) )

        ! Define "W" and its attributes 
        CALL CHECK( nf90_def_var(ncid, w%vname, nf90_float, dimid1, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 5, des,  w%desc) )
        CALL CHECK( nf90_put_att(ncid, 5,  un, w%units) )

        ! Define "T" and its attributes 
        CALL CHECK( nf90_def_var(ncid, temp%vname, nf90_double, dimid4, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 6, des,   temp%desc) )
        CALL CHECK( nf90_put_att(ncid, 6,  un,  temp%units) )

        ! Define "Q" and its attributes 
        CALL CHECK( nf90_def_var(ncid, q%vname, nf90_float, dimid4, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 7, des,  q%desc) )
        CALL CHECK( nf90_put_att(ncid, 7,  un, q%units) )

                ! End define mode.
        CALL CHECK( nf90_enddef(ncid) )

        ! These settings tell netCDF to write one timestep of data.
        dim1_count = (/ nz /)
        dim1_start = (/ 1  /)

        ! Write the pretend data.
        CALL CHECK( nf90_put_var(ncid, 1,     x) )
        CALL CHECK( nf90_put_var(ncid, 2,     y) )
        CALL CHECK( nf90_put_var(ncid, 3, z_out) )
        CALL CHECK( nf90_put_var(ncid, 4,  date) )

        CALL CHECK( nf90_put_var(ncid, 5, w_out, start=dim1_start, count=dim1_count) )
        CALL SUCCESS_MSG("w")

        dim4_count = (/ 1,  1, nz, nnt /)
        dim4_start = (/ 1,  1,  1,  1  /)
        ! do irec = 1, nnt
          ! dim4_start(4) = 1
          ! dim4_count(4) = irec

          CALL CHECK( nf90_put_var(ncid, 6, temp_out(:,:), start=dim4_start, count=dim4_count) )
          CALL CHECK( nf90_put_var(ncid, 7,    q_out(:,:), start=dim4_start, count=dim4_count) )
        ! end do 
        CALL SUCCESS_MSG("T")
        CALL SUCCESS_MSG("Q")

        ! Close the file. This causes netCDF to flush all buffers.
        CALL CHECK( nf90_close(ncid) ) 

      END SUBROUTINE Sub_write_netcdf!}}}

      SUBROUTINE Sub_write_netcdf_drop ( nz, nt, z_out,       &!{{{
                                         nbin, ref_r, dr, dz,     &
                                         num, m, S,            &
                                         out_path, out_name  )

        USE NETCDF
        IMPLICIT NONE
        ! In
        INTEGER,                  INTENT(IN)  ::    nz, nt, nbin
        REAL, DIMENSION(nz),      INTENT(IN)  ::     z_out, dz
        REAL, DIMENSION(nbin),    INTENT(IN)  ::     ref_r, &
                                                         m, &
                                                        dr

        REAL, DIMENSION(nz,nt+1),      INTENT(IN) ::    S
        REAL, DIMENSION(nbin,nz,nt+1), INTENT(IN) ::  num
        CHARACTER(LEN=*),       INTENT(IN)   ::  out_path,  &
                                                 out_name
        ! Local
        INTEGER                              ::      irec
        INTEGER                              ::       nnt,  &
                                                        y
        INTEGER, DIMENSION(nt+1)             ::      date

        nnt=nt+1
        y=1

        DO it = 1, nnt
          date(it) = it
        ENDDO

        CALL Sub_nc_attri 

        ! Create the file.
        CALL CHECK( nf90_create(trim(out_path)//"/"//trim(out_name),nf90_clobber,ncid) )
        CALL SUCCESS_MSG("Creation of new output data")

        ! Define the dimensions. The time dimension has no limit.
        CALL check( nf90_def_dim(ncid,      "Lon",           nbin, lon_dimid) )
        CALL check( nf90_def_dim(ncid,      "Lat",              y, lat_dimid) )
        CALL check( nf90_def_dim(ncid,    z%vname,             nz, lev_dimid) )
        CALL check( nf90_def_dim(ncid, time%vname, nf90_unlimited, rec_dimid) )

        ! The dimids array is used to pass the dimids of the dimension of
        ! the netCDF variables. In Fortran, the unlimited dimension
        ! must come last on the list of dimids.
        dimid4 = (/ lon_dimid, lat_dimid, lev_dimid, rec_dimid /)
        dimid2 = (/ lev_dimid, rec_dimid /)
        dimid1 = (/ lon_dimid /)
        dimid11 = (/ lev_dimid /)


        ! Define "x","y" and its attributes 
        CALL check( nf90_def_var(ncid, "Lon", nf90_double, lon_dimid, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 1,  ax,  "X") )
        CALL check( nf90_def_var(ncid, "Lat", nf90_double, lat_dimid, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 2,  ax,  "Y") )

        ! Define "lev" and its attributes 
        CALL check( nf90_def_var(ncid, z%vname, nf90_double, lev_dimid, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 3, des,  z%desc) )
        CALL CHECK( nf90_put_att(ncid, 3,  un, z%units) )
        CALL CHECK( nf90_put_att(ncid, 3,  ax,  z%axis) )

        ! Define "time" and its attributes 
        CALL check( nf90_def_var(NCID, time%vname, nf90_double, rec_dimid, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 4,  un, time%units) )
        CALL CHECK( nf90_put_att(ncid, 4,  ax,  time%axis) )

        ! Defin "num" and its attributes 
        CALL CHECK( nf90_def_var(ncid, "num", nf90_double, dimid4, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 5, des,   "number conc.") )
        CALL CHECK( nf90_put_att(ncid, 5,  un,          "# m-3") )

        ! Defin "r" and its attributes 
        CALL CHECK( nf90_def_var(ncid, "r", nf90_double, dimid1, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 6, des,   "radius.") )
        CALL CHECK( nf90_put_att(ncid, 6,  un,          "m ") )

        ! Defin "m" and its attributes 
        CALL CHECK( nf90_def_var(ncid, "m", nf90_double, dimid1, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 7, des,   "mass.") )
        CALL CHECK( nf90_put_att(ncid, 7,  un,          "kg ") )

        ! Defin "S" and its attributes 
        CALL CHECK( nf90_def_var(ncid, "S", nf90_double, dimid2, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 8, des,   "RH") )
        CALL CHECK( nf90_put_att(ncid, 8,  un,          "% ") )

        ! Defin "dr" and its attributes 
        CALL CHECK( nf90_def_var(ncid, "dr", nf90_double, dimid1, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 9, des,   "delta raidus") )
        CALL CHECK( nf90_put_att(ncid, 9,  un,          "m") )

        ! Defin "dr" and its attributes 
        CALL CHECK( nf90_def_var(ncid, "dz", nf90_double, dimid11, varid=varid) )
        CALL CHECK( nf90_put_att(ncid, 10, des,   "delta z") )
        CALL CHECK( nf90_put_att(ncid, 10,  un,          "m") )
                ! End define mode.
        CALL CHECK( nf90_enddef(ncid) )

        ! These settings tell netCDF to write one timestep of data.
        dim1_count = (/ nbin /)
        dim1_start = (/ 1  /)
        dim2_count = (/ nz, nnt /)
        dim2_start = (/ 1, 1  /)

        ! Write the pretend data.
        CALL CHECK( nf90_put_var(ncid, 1, ref_r) )
        CALL CHECK( nf90_put_var(ncid, 2,     y) )
        CALL CHECK( nf90_put_var(ncid, 3, z_out) )
        CALL CHECK( nf90_put_var(ncid, 4,  date) )

        dim4_count = (/ nbin,  1, nz, nnt /)
        dim4_start = (/ 1,  1,  1,  1  /)
        ! do irec = 1, nnt
          ! dim4_start(4) = 1
          ! dim4_count(4) = irec
  
          CALL CHECK( nf90_put_var(ncid, 5, num(:,:,:), start=dim4_start, count=dim4_count) )
          CALL CHECK( nf90_put_var(ncid, 6, ref_r, start=dim1_start, count=dim1_count) )
          CALL CHECK( nf90_put_var(ncid, 7, m, start=dim1_start, count=dim1_count) )
          CALL CHECK( nf90_put_var(ncid, 8, S, start=dim2_start, count=dim2_count) )
          CALL CHECK( nf90_put_var(ncid, 9, dr, start=dim1_start, count=dim1_count) )
        dim1_count = (/ nz /)
        dim1_start = (/ 1  /)
          CALL CHECK( nf90_put_var(ncid, 10, dr, start=dim1_start, count=dim1_count) )
        ! end do 
        CALL SUCCESS_MSG("m")
        CALL SUCCESS_MSG("num")

        ! Close the file. This causes netCDF to flush all buffers.
        CALL CHECK( nf90_close(ncid) ) 

      END SUBROUTINE Sub_write_netcdf_drop!}}}
  END MODULE Mod_write
