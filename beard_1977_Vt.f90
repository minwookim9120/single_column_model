MODULE beard_Vt

    USE Mod_const

    SUBROUTINE beard_1977( T, p, r, V )

    IMPLICIT NONE

    !----- Mod_const -----

    REAL, PARAMETER                 ::      g            , & ! 9.8 [m*s^-2] , acceleration of gravity
                                            pi           , & ! pi = 3.141591
                                            Rd           , & ! 287 [J/(kg*K)] ; J = N*m =(kg*m/s^2)*m
                                            rho              ! density of water [kg * m^-3]
    !---------------------

    
    !----- Requires definition or calculation -----

    REAL                            ::      T0           , & ! T0 =  293.15 K   = 20'C  ┐ Basic
                                            p0           , & ! p0 = 1013.25 hPa = 1 atm ┘ statement
                                            l0           , & ! [ cm ]
                                             l           , & ! Mean free path of air molcules [m]
                                          rho0           , & ! density of air at basic statment
                                                           & ! [ kg * m^-3 ] at 20'C and 1 atm 
                                          rhov           , & ! density of air [ kg * m^-3 ]
                                           eta           , & ! [ kg * m^-1 * s^-1 ] , dynamic viscosity
                                                           & ! A short course in cloudphysics 3rd, 102 page
                                          eta0           , & ! dynamic viscosity at 20'C and 1 atm
                                             D           , & ! diameter of droplet [ micrometer ]
                                             f           , & ! adjustment factor
                                         eps_c           , & ! velocity deviation for a 
                                                           & ! constant drag coefficient
                                         eps_s               ! velocity deviation for
                                                             ! Stokes drag

    !----------------------------------------------

    !---- beard 1977 ----

    REAL,               INTENT(IN)  ::       T           , & ! Temperature          [ K ]
                                             p               ! Atmospheric pressure [ hPa ]


    REAL, DIMENSION(:), INTENT(IN)  ::       r               ! radius [ m ]

    REAL, DIMENSION(:), INTENT(OUT) ::       V               ! terminal velocity

    !--------------------


    !-------- V0 --------

    REAL,               INTENT(IN)  ::      T0           , & ! 293.15 K    ┐ Basic
                                            p0               ! 1013.25 hPa ┘ statement 

    REAL, DIMENSION(:), INTENT(OUT) ::      V0               ! terminal velocity at basic statement

    !--------------------


    l0          = 6.62 * (1./10.)**(6.)     ! [ cm ]

    rhov        = p  * 100. / ( Rd * T  )   ! density of air 

    eta         = 1.72 * (1./10.)**(5.) * ( 393. / ( T + 120. ) ) * ( ( T / 273. )**(3./2.) )


    ! density of air and dynamic viscosity at basic statement ( 1 atm, 20C )

    T0          = 293.15                    ! 20 'C ┐ Basic
    p0          = 1013.25                   ! 1 atm ┘ statement

    rho0        = p0 * 100. / ( Rd * T0 )  
    eta0        = 1.72 * (1./10.)**(5.) * ( 393. / ( T0 + 120. ) ) * ( ( T0 / 273. )**(3./2.) )


    ! V0
    DO irr  = 1, size(r)        
        CALL terminal_velocity( T0 , p0, r(irr), V0(irr) ) ! beard 1976 : T = 20C, P = 1013.25
    END DO



    !---------------------------------------- beard 1977 ------------------------------------------

    D           = 2. * r * 1e+6             ! radius -> diameter & [ meter ] -> [ micrometer ]

    !------------------------------------- Two diameter range -------------------------------------

    !--------------------------------------- Regime 1 ---------------------------------------------

    IF ( D < 1. ) THEN


        V = 0.             ! [ m/s ]


    ELSE IF ( ( D >= 1. ) .and. ( D < 40. ) ) THEN

        ! [ m ]
        l       = l0 * 0.01 * ( eta / eta0 ) * ( ( p0 * rho0 ) / ( p * rhov ) )**(1./2.) 

        ! [ dimensionless ] velocity adjustment factor
        f       = ( eta0 / eta ) * ( 1 + 2.51 * ( l / ( D * 1e-6 ) ) ) / &
                                   ( 1 + 2.51 * ( ( l0 * 0.01 ) / ( D * 1e-6 ) ) ) 

        ! [ m * s^-1 ] = [ dimensionless ] * [ m * s^-1 ] velocity at 1 micrometer ~ 40 micrometer   
        V       = V0 * f         


    !--------------------------------------- Regime 2 --------------------------------------------!
    !                                                                                             !
    !                      The unit of input diameter must be [ cm ] !!!                          !
    !                                                                                             !
    !---------------------------------------------------------------------------------------------!


    ELSE IF ( ( D >= 40 ) .and. ( D < 6e+3 ) ) THEN

        ! [ dimensionless ] = [ kg * m^-3 ] * [ kg^-1 * m^3 ]
        eps_c   = ( ( rho0 / rhov )**(1./2.) ) - 1.        

        ! [ dimensionless ] = [ kg * m^-1 * s^-1 ] * [ kg^-1 * m * s ]
        eps_s   = ( eta0 / eta ) - 1.                    

        ! [ dimensionless ]
        f       = 1.104 * eps_s + ( ( ( 1.058 * eps_c ) - ( 1.104 * eps_s) ) &
                                                * ( 5.52 + log( D*1e-4 ) ) / 5.01 )   + 1. 

        V       = V0 * f

    ELSE IF ( ( D >= 6e+3) )  THEN

        ! [ dimensionless ] = [ kg * m^-3 ] * [ kg^-1 * m^3 ]
        eps_c   = ( ( rho0 / rhov )**(1./2.) ) - 1.        

        ! [ dimensionless ] = [ kg * m^-1 * s^-1 ] * [ kg^-1 * m * s ]
        eps_s   = ( eta0 / eta ) - 1.                     

        ! [ dimensionless ]
        f     = 1.104 * eps_s + ( ( ( 1.058 * eps_c ) - ( 1.104 * eps_s) ) & 
                                             * ( 5.52 + log( 6e+3*1e-4 ) ) / 5.01 )   + 1. 

        V     = V0 * f


    ENDIF


    end do


    END SUBROUTINE


END MODULE 
