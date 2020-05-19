MODULE Mod_phys_driver 

  USE Mod_global
  USE Mod_distribution
  USE Mod_const

  IMPLICIT NONE

CONTAINS

   SUBROUTINE Sub_cal_dm_dt_reasign
     IMPLICIT NONE
     INTEGER :: im, inext_m, it
     REAL    :: T, es, e, p, z, H, Tv, qv, RH, S
     REAL    :: N, NM
     REAL    :: re_num_m, re_num_m_1
     REAL, DIMENSION(drop_column_num) ::before_r
     REAL, DIMENSION(nt) :: total_m

     drop%next_num = 0 
     drop%dTemp = 0
   
     !!------to input data
     drop%next_num(1,:) = drop%num(:)
     drop%next_m(1,:) = drop%m
     before_r = drop%r
     total_m(1) = SUM(drop%next_num(1,:)*drop%m(:))

     T = 100.
     qv = 0.2 !!! qv ???? 
     z = 300
     !!---------------------------------!!

     DO it = 2, nt

       Tv = T*(1+(0.61*qv))
       H = (Rd*Tv)/g
       p = Ps*exp(-z/H)
       e = p*qv/0.622 !!! ??? e =?????
       es = 611.2*EXP((17.67*T)/(T+243.5)) !! ok

       RH = e/es 
!       S = RH -1
        S = 0.01

       drop%dm_dt = 4.*pi*before_r*1. / (((Rv*T)/(Dv*es)) + &
                         (((L/(Rv*T))-1)*(L/(Ka*T)))) *S

       drop%next_m(it,:) = drop%next_m(it-1,:)+ drop%dm_dt * dt
print*, 'time', it


       DO inext_m = 1, drop_column_num
       DO im = 1, drop_column_num-1

         IF (drop%next_m(it, inext_m) > drop%m(im) .AND. &
             drop%next_m(it, inext_m) < drop%m(im+1)) THEN

       !    N  = drop%next_num(it-1,inext_m)
       !    NM = drop%next_num(it-1,inext_m) * drop%next_m(it,inext_m)
           N  = drop%next_num(1,inext_m)
           NM = drop%next_num(1,inext_m) * drop%next_m(it,inext_m)

           re_num_m   = (NM - (drop%m(im+1) * N)) / &
                        (drop%m(im)-drop%m(im+1))
           re_num_m_1 = N - re_num_m 

           drop%next_num(it, im) = drop%next_num(it, im) + &
                                        re_num_m
           drop%next_num(it, im+1) = drop%next_num(it, im+1) + &
                                          re_num_m_1

         ENDIF

       ENDDO
       ENDDO

       before_r    = ( 3 * drop%next_m(it,:) * 1. / ( rho * 4 * pi ) )**( 1./3. ) 
       total_m(it) = SUM(drop%next_num(it,:)*drop%m(:))
       drop%dqv    = -1 * (total_m(it) - total_m(it-1))
       drop%dTemp  = -1 * (L* drop%dqv)/(rho*Cp)

       qv = qv + drop%dqv
       T = T + drop%dTemp 
print*, 'RH      :', RH
print*, 'S       :', S
print*, 'total_m :', total_m(it)
print*, 'dqv     :',drop%dqv
print*, 'qv      :',qv      
print*, 'T       :',T 

     ENDDO

   END SUBROUTINE Sub_cal_dm_dt_reasign

END MODULE Mod_phys_driver 
