MODULE Mod_realcase

  USE Mod_global

  IMPLICIT NONE

    CONTAINS

    !!---------------------------------------------!!
    !!  Cal. vertical coordinate                   !!
    !!---------------------------------------------!!
    SUBROUTINE Sub_real_init                        & 
               (                                    &
                read_pres, read_psfc, read_temp,    &
                z, input_nz, nz,                    &
                qv_in,  temp_in, vert_in, w_in,     &
                qv_out, temp_out,  p_out, w_out       &
               )

    USE Mod_const, only: Ps, Rd, g, Cp
    IMPLICIT NONE
    ! In
    INTEGER,                   INTENT(IN)  :: read_pres, read_psfc,  & !!
                                              read_temp, input_nz,   &
                                              nz
    REAL, DIMENSION(input_nz), INTENT(IN)  ::     qv_in,          &  !! [kg/kg]
                                                temp_in,          &  !! [K]
                                                vert_in              !! P [hPa] or Z [m]
    REAL, DIMENSION(input_nz), INTENT(IN)  ::      w_in
    REAL, DIMENSION(nz),       INTENT(IN)  ::         z              !! [kg/kg]
    ! Out
    REAL, DIMENSION(nz),        INTENT(OUT) :: p_out, temp_out, qv_out
    REAL, DIMENSION(nz),        INTENT(OUT) :: w_out
    ! Local
    REAL, DIMENSION(input_nz) :: z_in, p_in, t_in, th_in,        &
                                   Tv,                           &  !! [K], virtual temperature
                                    H                               !! [m] scale height
    INTEGER :: i, k
    REAL    :: t 

    if (read_pres == 1 ) then
        if (read_temp == 1 ) then        ! input data : P[hPa] & theta[K]
            th_in = temp_in
            p_in  = vert_in
            t_in  = th_in*((P_in/Ps)**(Rd/Cp))
            Tv    = t_in*(1+(0.61*qv_in))
            H     = (Rd*Tv)/g
            z_in  = -H*(log(P_in/Ps))
        else                            ! input data : P[hPa] & T[K]
            P_in = vert_in
            T_in = temp_in
            Tv   = T_in*(1+(0.61*qv_in))
            H    = (Rd*Tv)/g
            z_in = -H*(log(P_in/Ps))
            Th_in = T_in*((Ps/P_in)**(Rd/cp))
        endif
    else
        if (read_temp == 1 ) then        ! input data : z[m] & theta[K]
            print*, ":: INPUT DATA VARIABLE ERROR ::"
            print*, ":: 'Pres' is calculated using 'Temp', and"
            print*, ":: 'Temp' is calculated using 'Pres'."
            print*, ":: Please check the input data variable."
        else                            ! input data : z[m] & T[K]
            z_in  = vert_in/g
            t_in  = temp_in
            Tv    = t_in*(1+(0.61*qv_in))
            H     = (Rd*Tv)/g
            p_in  = Ps*exp(-(z_in/H))
            th_in = t_in*((Ps/p_in)**(Rd/cp))
        endif
    endif

!!! :: Interpolate to fit nz

    do i = 1, nz

      if ( z(i) <= z_in(1) ) then
        t = (z(i) - z_in(1)) / (z_in(2) - z(1))
        temp_out(i) = (1.0D+00 - t)*t_in(1)  + t*t_in(2)
           p_out(i) = (1.0D+00 - t)*p_in(1)  + t*p_in(2)
          qv_out(i) = (1.0D+00 - t)*qv_in(1) + t*qv_in(2)
           w_out(i) = (1.0D+00 - t)*w_in(1)  + t*w_in(2)
      else if ( z_in(input_nz) <= z(i) ) then
        t = ( z(i) - z_in(input_nz-1) ) / (z_in(input_nz) - z_in(input_nz-1))
        temp_out(i) = (1.0D+00 - t)*t_in(input_nz-1)  + t*t_in(input_nz)
           p_out(i) = (1.0D+00 - t)*p_in(input_nz-1)  + t*p_in(input_nz)
          qv_out(i) = (1.0D+00 - t)*qv_in(input_nz-1) + t*qv_in(input_nz)
           w_out(i) = (1.0D+00 - t)*w_in(input_nz-1)  + t*w_in(input_nz)
      else
        do k = 2, input_nz

          if ( z_in(k-1) <= z(i) .and. z(i) <= z_in(k) ) then
            t = (z(i) - z_in(k-1)) / (z_in(k) - z_in(k-1))
            temp_out(i) = (1.0D+00 - t)*t_in(k-1)  + t*t_in(k)
               p_out(i) = (1.0D+00 - t)*p_in(k-1)  + t*p_in(k)
              qv_out(i) = (1.0D+00 - t)*qv_in(k-1) + t*qv_in(k)
               w_out(i) = (1.0D+00 - t)*w_in(k-1) + t*w_in(k)
            exit
          end if

        end do
      end if

    end do 


    END SUBROUTINE Sub_real_init 
   
END MODULE Mod_realcase
