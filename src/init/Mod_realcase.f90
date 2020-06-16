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
                z, input_nz, nz, input_ps, input_t, &
                qv_in,  temp_in, vert_in,     &
                qv_out, temp_out,  p_out     &
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
    REAL,                      INTENT(IN)  :: input_Ps, input_t 
    ! REAL, DIMENSION(input_nz), optional, INTENT(IN)  ::      w_in
    REAL, DIMENSION(nz),       INTENT(IN)  ::         z              !! [kg/kg]
    ! Out
    REAL, DIMENSION(nz),        INTENT(OUT) :: p_out, temp_out, qv_out
    ! REAL, DIMENSION(nz),   optional,     INTENT(OUT) :: w_out
    ! Local
    REAL, DIMENSION(input_nz) :: z_in, p_in, t_in, th_in,        &
                                   Tv, rho_in, dp_in, dz_in,     &  !! [K], virtual temperature
                                    H                               !! [m] scale height
    INTEGER :: i, k
    REAL    :: t, ka, R, p0, pss

    if (read_pres == 1 ) then
      if (read_temp == 1 ) then        ! input data : P[hPa] & theta[K]
      !   th_in = temp_in
      !   p_in  = vert_in
      !   t_in  = th_in*((P_in/Ps)**(Rd/Cp))
      !   Tv    = t_in*(1+(0.61*qv_in))
      !   H     = (Rd*Tv)/g
      !   z_in  = -H*(log(P_in/Ps))
      ! else                            ! input data : P[hPa] & T[K]
      !   P_in = vert_in
      !   T_in = temp_in
      !   Tv   = T_in*(1+(0.61*qv_in))
      !   H    = (Rd*Tv)/g
      !   z_in = -H*(log(P_in/Ps))
      !   Th_in = T_in*((Ps/P_in)**(Rd/cp))
      endif
    else
      if (read_temp == 1 ) then        ! input data : z[m] & theta[K]
ka=0.286   !(J/kg k)
p0=100000. !(Pa)
R=287.     !(J/kg k)
! g=9.8      !(m/s^2)
pss=input_ps*100.  !convert unit hPa->pa

!!convert z->p
      p_in(1) = pss
      rho_in(1)=pss/(R*input_t)
      dz_in(1) = vert_in(1)-0
      dp_in(1)=-rho_in(1)*g*dz_in(1)
      p_in(1)=pss+dp_in(1)

  do k=1, input_nz-1
      t_in(k)=temp_in(k)*(p_in(k)/p0)**(ka)
      rho_in(k)=p_in(k)/(R*t_in(k))
      if ( k == 1 ) then
        dz_in(k)=vert_in(k)-0.
      else
        dz_in(k)=vert_in(k)-vert_in(k-1)
      endif
      dp_in(k)=-rho_in(k)*g*dz_in(k)
      if ( k == 1 ) then
        p_in(k+1)=pss+dp_in(k)
      else
        p_in(k+1)=p_in(k)+dp_in(k)
      endif

    ! if ( k == 70 ) then
    !   write(*,*) t_in(k)
    !   write(*,*) rho_in(k)
    !   write(*,*) dz_in(k)
    !   write(*,*) dp_in(k)
    !   write(*,*) p_in(k)
    !   stop
    ! endif
  enddo
      t_in(input_nz)=temp_in(input_nz)*(p_in(input_nz)/p0)**(ka)
      z_in  = vert_in

    !   write(*,*) t_in, shape(t_in)
    !   write(*,*) p_in/100
    !    write(*,*) vert_in 
    ! stop

 p_in=p_in/100.
 
        ! z_in  = vert_in
        ! H     = 8000.
        ! p_in  = input_Ps*exp(-(z_in/H))
        ! Th_in = temp_in
        ! T_in = Th_in*((P_in/Ps)**(Rd/Cp))
      else                            ! input data : z[m] & T[K]
        ! z_in  = vert_in/g
        ! t_in  = temp_in
        ! Tv    = t_in*(1+(0.61*qv_in))
        ! H     = (Rd*Tv)/g
        ! p_in  = Ps*exp(-(z_in/H))
        ! th_in = t_in*((Ps/p_in)**(Rd/cp))
      endif
    endif
!!! :: Interpolate to fit nz

    do i = 1, nz

      if ( z(i) <= z_in(1) ) then
        t = (z(i) - z_in(1)) / (z_in(2) - z(1))
        temp_out(i) = (1.0D+00 - t)*t_in(1)  + t*t_in(2)
           p_out(i) = (1.0D+00 - t)*p_in(1)  + t*p_in(2)
          qv_out(i) = (1.0D+00 - t)*qv_in(1) + t*qv_in(2)
           ! w_out(i) = (1.0D+00 - t)*w_in(1)  + t*w_in(2)
      else if ( z_in(input_nz) <= z(i) ) then
        t = ( z(i) - z_in(input_nz-1) ) / (z_in(input_nz) - z_in(input_nz-1))
        temp_out(i) = (1.0D+00 - t)*t_in(input_nz-1)  + t*t_in(input_nz)
           p_out(i) = (1.0D+00 - t)*p_in(input_nz-1)  + t*p_in(input_nz)
          qv_out(i) = (1.0D+00 - t)*qv_in(input_nz-1) + t*qv_in(input_nz)
           ! w_out(i) = (1.0D+00 - t)*w_in(input_nz-1)  + t*w_in(input_nz)
      else
        do k = 2, input_nz

          if ( z_in(k-1) <= z(i) .and. z(i) <= z_in(k) ) then
            t = (z(i) - z_in(k-1)) / (z_in(k) - z_in(k-1))
            temp_out(i) = (1.0D+00 - t)*t_in(k-1)  + t*t_in(k)
               p_out(i) = (1.0D+00 - t)*p_in(k-1)  + t*p_in(k)
              qv_out(i) = (1.0D+00 - t)*qv_in(k-1) + t*qv_in(k)
               ! w_out(i) = (1.0D+00 - t)*w_in(k-1) + t*w_in(k)
            exit
          end if

        end do
      end if

    end do 

    END SUBROUTINE Sub_real_init 
   
END MODULE Mod_realcase
