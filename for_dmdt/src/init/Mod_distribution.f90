MODULE Mod_distribution

  USE Mod_global
  USE Mod_const
  USE Mod_read, only: Sub_read_NC_file

  IMPLICIT NONE

    CONTAINS
    
    !!---------------------------------------------!!
    !!  Set. size distributions                    !!
    !!---------------------------------------------!!
    SUBROUTINE Sub_drop_distributions     &
               (                          &
                distribution_option,      &
                drop_column_num,          & !!  drop_column_num
                rbmin,                    & !!  drop_1st_diameter
                rbmax,                    & !!  drop_1st_diameter
                nr,r,rb,m,mb              & !!  drop_num
               )

    USE Mod_const, only: r0,rho,nc,qc,pi
    IMPLICIT NONE

    ! IN
    INTEGER, INTENT(IN)                      :: distribution_option, &
                                                drop_column_num              
    REAL,    INTENT(IN)                      :: rbmin,               &
                                                rbmax
    ! Local
    REAL, DIMENSION(drop_column_num)         :: D,  &
                                                dr, &
                                                dD, &
                                                dlnr

    INTEGER                                  :: ir, nbin
    REAL                                     :: ratio, &
                                                mu, sigma, lambda
    ! OUT
    ! REAL, DIMENSION(nz,column_num), INTENT(OUT) :: conc
    REAL, DIMENSION(drop_column_num),    INTENT(OUT) :: nr, r, m
    REAL, DIMENSION(drop_column_num+1),  INTENT(OUT) :: rb, mb
    
    nbin  = drop_column_num

    rb(1) = rbmin ; rb(nbin+1) = rbmax       ! boundary of drop diameter
    ratio = (rb(nbin+1)/rb(1))**(1./nbin)
    DO ir = 2, nbin
      rb(ir) = ratio*rb(ir-1)  
    ENDDO

    do ir   = 1, nbin+1
      mb(ir) = rho * 4. / 3. * pi * ( rb(ir) )**3
    end do

    do ir   = 1, nbin
        m(ir) = ( mb(ir+1) + mb(ir) ) / 2.
        r(ir) = ( 3 * m(ir) * 1. / ( rho * 4 * pi ) )**( 1./3. )
       dr(ir) = ( rb(ir+1) - rb(ir) )
     dlnr(ir) = ( dr(ir)/r(ir) )
    end do 

    SELECT CASE(distribution_option)
      CASE(1)

        ! Log-normal distribution
        mu = LOG(r0)
        sigma = SQRT((2./9.)*LOG(qc/(nc*rho*(4./3.)*pi*(r0)**3)))
        DO ir = 1, nbin
          nr(ir) = (nc/((SQRT(2.*pi))*sigma*r(ir))) * EXP((-1.*(LOG(r(ir))-mu)**2)/(2.*(sigma**2))) &
                    * dr(ir)
        ENDDO
        write(*,*) "gamma dist. in model (qc) =  ", sum(m*nr)
        write(*,*) "gamma dist. const.   (qc) =  ", qc
        write(*,*) "log-normal dist. in model (nc) =  ", sum(nr)
        write(*,*) "log-normal dist. const.   (nc) =  ", nc
        write(*,*) "accuracy  = ", (sum(nr)/nc) * 100., "%"

      CASE(2)

        ! Gamma distribution
        D = 2*r
        dD = 2*dr
        mu = (1.0e+9/nc) + 2.
        lambda = ((nc/qc)*(GAMMA(mu+4.)/GAMMA(mu+1.))*pi*(1./6.)*rho) ** (1./3.)
        !n0 = (nc*lamb**(mu+1.))/GAMMA(mu+1.)

        DO ir = 1, nbin
          !nr(ir) = n0*r(ir)^mu*EXP(-1*lamb*r(ir))
          nr(ir) = (nc/gamma(mu+1.)) * lambda*((lambda*D(ir))**mu) * exp(-1*lambda*D(ir))*dD(ir)
        ENDDO
  
        write(*,*) "gamma dist. in model (qc) =  ", sum(m*nr)
        write(*,*) "gamma dist. const.   (qc) =  ", qc
        write(*,*) "gamma dist. in model (nc) =  ", sum(nr)
        write(*,*) "gamma dist. const.   (nc) =  ", nc
        write(*,*) "accuracy  = ", (sum(nr)/nc) * 100., "%"

    END SELECT

    ! DO ir = 1, nbin
    !   rr(ir) = dr(ir)*nr(ir)
    ! ENDDO
    !
    ! write(*,*) sum(rr)
    ! write(*,*) maxval(nr)

    ! conc(:,:)=0.
    END SUBROUTINE Sub_drop_distributions 

END MODULE Mod_distribution
