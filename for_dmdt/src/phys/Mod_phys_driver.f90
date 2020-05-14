MODULE Mod_phys_driver 

  USE Mod_global
  USE Mod_distribution
  USE Mod_const

  IMPLICIT NONE

CONTAINS

   SUBROUTINE Sub_cal_dm_dt_reasign
     IMPLICIT NONE
     INTEGER :: im, inext_m, it
     REAL    :: T, es
     REAL    :: N, NM
     REAL    :: re_num_m, re_num_m_1
     REAL, DIMENSION(drop_column_num) :: before_m, before_r

     drop%next_num = 0 
     drop%next_num(1,:) = drop%num(:)
     drop%next_m(1,:) = drop%m
     before_r = drop%r
!     drop%dTemp = 0
!print*, 'mass', drop%m
!print*, 'num ', drop%num
     DO it = 2, nt
!     DO it = 2, 3 

!       T = temp%dout(30,it)
!       T = SIN(REAL(it))* 10 !+ 280.
!       T = 100. + drop%dTemp 
       T = 100. 

       es = 611.2*EXP((17.67*T)/(T+243.5))
       drop%dm_dt = 4.*pi*before_r*1. / (((Rv*T)/(Dv*es)) + &
                         (((L/(Rv*T))-1)*(L/(Ka*T)))) *S

       drop%next_m(it,:) = drop%next_m(it-1,:)+ drop%dm_dt * dt
!print*, 'mass'
print*, 'time', it
!print*, drop%m
!print*, 'next mass', drop%next_m(it,:)
!print*, drop%next_m  
!print*, drop%num
       DO inext_m = 1, drop_column_num
       DO im = 1, drop_column_num-1

         IF (drop%next_m(it, inext_m) > drop%m(im) .AND. &
             drop%next_m(it, inext_m) < drop%m(im+1)) THEN

           N  = drop%next_num(it-1,inext_m)
           NM = drop%next_num(it-1,inext_m) * drop%next_m(it,inext_m)

           re_num_m   = (NM - (drop%m(im+1) * N)) / &
                        (drop%m(im)-drop%m(im+1))
           re_num_m_1 = N - re_num_m 

           drop%next_num(it, im) = drop%next_num(it, im) + &
                                        re_num_m
           drop%next_num(it, im+1) = drop%next_num(it, im+1) + &
                                          re_num_m_1
!print*, 'r    : ', inext_m
!print*, 'next mass', drop%next_m(it,inext_m)
!print*, 'mass', drop%m(im), drop%m(im+1)
!print*, "N    : ",N
!print*, "NM   : ", NM
!print*, "Ni, Ni+1  : ", re_num_m, re_num_m_1 
!print*, "Ni+ Ni+1  : ", re_num_m+ re_num_m_1 
!print*, "Ni, Ni+1  : ", drop%next_num(it, im), drop%next_num(it,im+1)
!print*, "Ni + Ni+1 :  ", drop%next_num(it, im) + drop%next_num(it,im+1)
!print*, "miNi + mi+1Ni+1 :  ", re_num_m*drop%m(im)+re_num_m_1*drop%m(im+1) 
         ENDIF

       ENDDO
       ENDDO

       before_r = ( 3 * drop%next_m(it,:) * 1. / ( rho * 4 * pi ) )**( 1./3. ) 
!       drop%dqv   = -1 * SUM(drop%num(:)*drop%next_m(:)) 
!       drop%dTemp  = -1 * (L* drop%dqv)/(rho*Cp)

!print*, 'next_num   : ', drop%next_num(it,:)
print*, 'Total_n    : ', SUM(drop%next_num(it,:))
print*, 'Total_mass(mass const) : ',SUM(drop%next_num(it,:)*drop%m(:))
print*, 'Total_mass(n const) : ', SUM(drop%num(:)*drop%next_m(it, :))
print*, 'Total_mass(test   ) : ', SUM(drop%next_num(it-1,:)*drop%next_m(it, :))

     ENDDO

   END SUBROUTINE Sub_cal_dm_dt_reasign

END MODULE Mod_phys_driver 
