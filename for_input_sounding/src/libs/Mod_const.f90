MODULE Mod_const

  IMPLICIT NONE

    REAL, PARAMETER     :: Rd  = 287.        !!! Rd = 287 J/(kg*K); J = N*m =(kg*m/s^2)*m
    REAL, PARAMETER     :: Cp  = 1005.       !!! Cp = 1005.J/(kg*K)
    REAL, PARAMETER     :: g   = 9.8         !!! unit = m/s
    REAL, PARAMETER     :: pi  = 4*ATAN(1.)  !!! pi = 3.141591   
    REAL, PARAMETER     :: Ps  = 1013        !!! hPa 
    REAL, PARAMETER     :: rho = 1000        !!! kg/m^3 
    REAL, PARAMETER     :: nc  = 1.0e+8      !!! #/m^3 
    REAL, PARAMETER     :: qc  = 0.002       !!! [kg / kg]
    REAL, PARAMETER     :: r0  = 1.0e-5      !!! [kg / kg]
         !== Refer to Rogers & Yau (1996), 103p - Table 7.1 (T=273K)
    REAL, PARAMETER     :: Dv  = 0.0000221   !!! vapor's diffusion coefficient [m^2 * s^-1]
    REAL, PARAMETER     :: Ka  = 0.024       !!! thermal conductivity          [J m^-1 s^-1 K^-1]
         !== 
    REAL, PARAMETER     :: Rv  = 461.5       !!! gas constant of water vapor   [J * Kg^-1 * K^-1]
    REAL, PARAMETER     :: S   = 0.01        !!! supersaturation             
    REAL, PARAMETER     :: L   = 2501000.    !!! latent heat at 0C             [J kg^-1]
    REAL, PARAMETER     :: courant_number = 1.   !!!



  CONTAINS

END MODULE Mod_const
