;------------------------------------ parameter list ---------------------------------------


    nbin    = 300               ; drop_column_num

    g       = 9.8               ; [ m * s^-2 ] , acceleration of gravity 

    pi      = 4*atan(1.0)       ; pi = 3.141591

    T       = 293.15            ; [ K ] , Thermodynamic temperature
    Tc      = T - 273.15        ; [ C ]
    T0      = 293.15            ; [ K ] , 20C

    l0      = 6.62 * 10^(-6)    ; [ cm ]

     p      = 1013.25           ; [ hPa = mb ]
    p0      = 1013.25           ; [ hPa = mb ] , standard sea level pressure   

    R       = 287               ; [ J kg^-1 K^-1 ]

    rho_i   = 1000                ; [ kg * m^-3  ] , density of water
    rho     = p * 100 / ( R * T ) ; [ kg * m^-3  ] , density of air
    del_rho = rho_i - rho         ; [ kg * m^-3  ] , difference between fluid density of drop and air
    rho0    = 1.204               ; [ kg * m^-3  ] , beard 1977

    eta0    = 1.818 * 10^(-4)   ; [  g * cm^-1 * s^-1 ] , at 20C 
    eta0   := eta0*0.1          ; [ kg * m^-1 * s^-1  ]

    eta     = 1.72 * 10^(-5) * ( 393 / ( T + 120) ) * ( T / 273 )^(3./2.)   ; [ kg * m^-1 * s^-1 ] , dynamic viscosity, A short course in cloudphysics 3rd, 102 page


;-------------------------------- Calculate radius of droplet -------------------------------


    rbmin   = 1e-6       ; min boundary r [m]
    rbmax   = 1e-2       ; max boundary r [m]

    rb   = new( nbin+1  , "double" )
    mb   = new( nbin+1  , "double" )


    m        = new( nbin  , "double" )
    r        = new( nbin  , "double" )

    rb(0)     = rbmin  ; boundary of drop radius
    rb(nbin)  = rbmax


    ratio   = ( rbmax / rbmin )^(1. / (nbin-1) )

        do ir   = 1, nbin-1

            rb(ir)  = rb(0) * ( ratio^(ir) )

        end do

        do ir   = 0, nbin

            mb(ir) = rho * 4. / 3. * pi * ( rb(ir) )^3

        end do


        do ir   = 0, nbin-1

            m(ir) = ( mb(ir+1) + mb(ir) ) / 2.
            r(ir) = ( 3 * m(ir) * 1. / ( rho * 4 * pi ) )^( 1./3. )

        end do


     D = 2*r   ; [ m ] radius -> diameter



;----------------------------------- Regime list ----------------------------------------


;------------------- Regime 1 ) 5e-7 <= d0 < 1.9e-5 [ m ] ----------------------------------


    d0_rg1_ind  = ind( ( D .ge. 5e-7 ) .and. ( D .lt. 1.9e-5 ) )  ; [ m ]


    d0_rg1      = new( dimsizes(d0_rg1_ind) , "double" )
    d0_rg1      = D(d0_rg1_ind)   


    l       = l0 * ( eta / eta0 ) * ( p0 / p ) * ( T / T0 )^(1/2)   ; [ cm ] mean free path of air molcules
    l      := l * 0.01                                              ; [ m  ] mean free path of air molcules


    C_sc_rg1= 1 + 2.51 * ( l / d0_rg1 )          ; slip correction factor
                                                 ; [ dimensionless ] = [ m ] * [ m^-1 ] 

    C_1     = ( del_rho * g ) / ( 18 * eta )     ; [ m^-1 * s^-1 ] = [ kg * m^-3 ] * [ m * s^-2 ] * [ kg^-1 * m * s ]

    V_inf1  = C_1 * C_sc_rg1 * d0_rg1^2          ; terminal Velocity at regime 1 
                                                 ; [ m * s^-1 ] = [ m^-1 * s^-1 ] * [ dimensionless ] * [ m^2 ]


;------------------ Regime 2 ) 1.9e-5 <= d0 < 0.00107 [ m ] ----------------------------------


    d0_rg2_ind  = ind( ( D .ge. 1.9e-5 ) .and. ( D .lt. 0.00107 ) ) ; [ m ]


    d0_rg2      = new( dimsizes(d0_rg2_ind) , "double" )
    d0_rg2      = D(d0_rg2_ind)
       
    C_2         = 4 * rho * ( del_rho ) * g / ( 3 * eta^2 ) ; [ m^-3 ] = [ kg * m^-3 ] * [ kg * m^-3 ] * [ m * s^-2 ] * [ kg^-2 * m^2 * s^2 ]

    N_Da        = C_2 * (d0_rg2)^3                          ; Davies number
                                                            ; [ dimensionless ] = [ m^-3 ] * [ m^3 ] 

    X_rg2       = log( N_Da )                               ; [ dimensionless ] 


    bn_rg2 = (/ -0.318657e+1 , 0.992696 , -0.153193e-2, -0.987059e-3, -0.578878e-3, 0.855176e-4, -0.327815e-5 /)


            Y_rg2 = bn_rg2(0) + bn_rg2(1) * X_rg2^1 + bn_rg2(2) * X_rg2^2  + bn_rg2(3) * X_rg2^3  + bn_rg2(4) * X_rg2^4 + bn_rg2(5) * X_rg2^5 + bn_rg2(6) * X_rg2^6 



    C_sc_rg2   = 1 + 2.51 * l / (d0_rg2)            ; slip correction factor
                                                    ; [ dimensionless ] = [ m ] * [ m^-1 ] 

    N_Re_rg2   = C_sc_rg2 * exp(Y_rg2)              ; reynolds number
                                                    ; [ dimensionless ] = [ dimensionless ] * [ dimensionless ]

    V_inf2     = eta * N_Re_rg2 / ( rho * d0_rg2 )  ; terminal velocity at regime 2
                                                    ; [ m * s^-1 ]    = [ kg * m^-1 * s^-1 ] * [ dimensionless ] * [ kg^-1 * m^3 ] * [ m^-1 ] 



;------------------ regime 3 ) 0.00107 <= d0 < 0.007 + alpha [m] ----------------------------------


    d0_rg3_ind  = ind( ( D .ge. 0.00107 ) .and. ( D .lt. 0.007 ) )


    d0_rg3      = new( dimsizes(d0_rg3_ind) , "double" )
    d0_rg3      = D(d0_rg3_ind)

    d0_max      = 0.007

    sigma       = 7.5 * 10^-2                                   ; approx. surface tension , A short course in cloudphysics 3rd. 85 page
                                                                ; [ N * m^-1 ] = [ kg * m * s^-2 * m^-1 ] = [ kg * s^-2 ]

    C_3         = 4 *  del_rho  * g  / ( 3 * sigma )            ; [ m^-2 ]     = [ kg * m^-3 ] * [ m * s^-2 ] * [ kg^-1 * s^2 ] 


    Bo          = C_3 * ( d0_rg3 )^2                            ; modified bond number 
                                                                ; [ dimensionless ] = [ m^-2 ] * [ m^2 ]

    Bo_max      = C_3 * ( d0_max )^2


    Np          = sigma^3 * rho^2 / ( eta^4 * del_rho * g )     ; physical property number
                                                                ; [ dimensionless ] = [ kg^3 * s^-6 ] * [ kg^2 * m^-6 ] * [ kg^-4 * m^4 * s^4 ] * [ kg^-1 * m^3 ] * [ m^-1 * s^2 ] 



    X_rg3       = log( Bo * Np^(1./6.) )                        ; [ dimensionless ] = [ dimensionless ] * [ dimensionless ]


    X_rg3_max   = log( Bo_max * Np^(1./6.) )


    bn_rg3      = (/ -0.500015e+1 , 0.523778e+1, -0.204914e+1, 0.475294, -0.542819e-1, 0.238449e-2/)


            Y_rg3 = bn_rg3(0) + bn_rg3(1) * X_rg3^1 + bn_rg3(2) * X_rg3^2  + bn_rg3(3) * X_rg3^3 + bn_rg3(4) * X_rg3^4 + bn_rg3(5) * X_rg3^5                     ; [ dimensionless ]
        Y_rg3_max = bn_rg3(0) + bn_rg3(1) * X_rg3_max^1 + bn_rg3(2) * X_rg3_max^2  + bn_rg3(3) * X_rg3_max^3 + bn_rg3(4) * X_rg3_max^4 + bn_rg3(5) * X_rg3_max^5 ; [ dimensionless ]


    N_Re_rg3    = Np^(1./6.) * exp(Y_rg3)                       ; reynolds number
                                                                ; [ dimensionless ] = [ dimensionless ] * [ dimensionless ]

    N_Re_rg3_max= Np^(1./6.) * exp(Y_rg3_max)


    V_inf3      = ( eta * N_Re_rg3 ) / ( rho * d0_rg3 )         ; terminal velocity at regime 3
                                                                ; [ m * s^-1 ]    = [ kg * m^-1 * s^-1 ] * [ dimensionless ] * [ kg^-1 * m^3 ] * [ m^-1 ]

    V_inf3_max  = ( eta * N_Re_rg3_max ) / ( rho * d0_max )     


;------------------------------ Total terminal velocity -------------------------------------


        V_inf_arr = new( dimsizes(V_inf1) + dimsizes(V_inf2) + dimsizes(V_inf3) , "double" )

        V_inf_arr(0:dimsizes(V_inf1)-1) = V_inf1
        V_inf_arr(dimsizes(V_inf1):dimsizes(V_inf1)+dimsizes(V_inf2)-1) = V_inf2
        V_inf_arr(dimsizes(V_inf1)+dimsizes(V_inf2):dimsizes(V_inf1)+dimsizes(V_inf2)+dimsizes(V_inf3)-1) = V_inf3


;--------------------------------- Various diameter -----------------------------------------


    V_inf  = new( nbin , "double" )


    if ( any( D .lt. 5e-7 ) .and. any( D .ge. 0.007 ) ) then


        D_min = max(ind( D .lt. 5e-7  ) )
        D_max = min(ind( D .ge. 0.007 ) )

        V_inf(0:D_min)           = 0
        V_inf(D_max:nbin-1)      = V_inf3_max
        V_inf(D_min+1:D_max-1)   = V_inf_arr


    else if ( any( D .lt. 5e-7 ) .and. all( D .lt. 0.007 ) ) then


        D_min = max(ind( D .lt. 5e-7  ) )
        D_max = nbin-1
        
        V_inf(0:D_min)           = 0
        V_inf(D_min+1:D_max)     = V_inf_arr


    else if ( all( D .ge. 5e-7 ) .and. any( D .ge. 0.007 ) ) then


        D_min = 0
        D_max = min(ind( D .ge. 0.007 ) )

        V_inf(D_min:D_max-1) = V_inf_arr
        V_inf(D_max:nbin-1)  = V_inf3_max


    else


        V_inf = V_inf_arr


    end if
    end if
    end if


;---------------------------------------- beard 1977 ------------------------------------------


;------------------------------------- Two diameter range -------------------------------------

;--------------------------------------- Regime 1 ---------------------------------------------



    brd_ind_1    = ind( ( D .ge. 1e-6 ) .and. ( D .lt. 40e-6 ) ) ; [ m ] 
    brd_val_1    = new( dimsizes(brd_ind_1) , "double" )
    brd_val_1    = D(brd_ind_1)                                  ; 1 micrometer ~ 40 micrometer


    l_brd   = l0 * ( eta / eta0 ) * ( ( p0 * rho0 ) / ( p * rho ) )^(1./2.) ; [ m ]

    brd_f_1 = ( eta0 / eta ) * ( 1 + 2.51 * ( l_brd / brd_val_1 ) ) / ( 1 + 2.51 * ( l0 / brd_val_1 ) ) ; [ dimensionless ] velocity adjustment factor

    brd_1_inf = brd_f_1 * V_inf(0:dimsizes(brd_ind_1)-1)         ; [ m * s^-1 ] = [ dimensionless ] * [ m * s^-1 ] velocity at 1 micrometer ~ 40 micrometer



;--------------------------------------- Regime 2 ---------------------------------------------

    brd_ind_2   = ind( ( D .ge. 40e-6 ) .and. ( D .lt. 0.006 ) )
    brd_val_2   = new( dimsizes(brd_ind_2) , "double" )
    brd_val_2   = D(brd_ind_2)                                   ; 40 micrometer ~ 6 millimeter


;    brd_rho     = 0.348 * p / T                                       ; [ kg * m^-3 ]        approx.
;
;    brd_eta     = 1.832 * 10^(-5) * ( 1 + 0.00266 * ( T - 296 ) )     ; [ kg * m^-1 * s^-1 ] approx.
;        
;    brd_eps_c   = ( rho / brd_rho )^(1./2.) - 1                       ; [ dimensionless ] = [ kg * m^-3 ] * [ kg^-1 * m^3 ]
;
;    brd_eps_s   = ( eta0 / brd_eta ) - 1                              ; [ dimensionless ] = [ kg * m^-1 * s^-1 ] * [ kg^-1 * m * s ]


    brd_eps_c   = ( ( rho0 / rho )^(1./2.) ) - 1.         ; [ dimensionless ] = [ kg * m^-3 ] * [ kg^-1 * m^3 ]

    brd_eps_s   = ( eta0 / eta ) - 1.                 ; [ dimensionless ] = [ kg * m^-1 * s^-1 ] * [ kg^-1 * m * s ]




    brd_f_2     = 1.104 * brd_eps_s + ( ( ( 1.058 * brd_eps_c ) - ( 1.104 * brd_eps_s) ) * ( 5.52 + log( brd_val_2 ) ) / 5.01 )  -1. 
    brd_f_2_max = 1.104 * brd_eps_s + ( ( ( 1.058 * brd_eps_c ) - ( 1.104 * brd_eps_s) ) * ( 5.52 + log( 0.006 ) )     / 5.01 )  -1. 



    brd_2_inf     = -1. * brd_f_2 * V_inf(dimsizes(brd_ind_1):dimsizes(brd_ind_1)+dimsizes(brd_ind_2)-1)  ; [ m * s^-1 ] 
    brd_2_inf_max = -1. * brd_f_2_max * V_inf(dimsizes(brd_ind_1)+dimsizes(brd_ind_2))                    ; [ m * s^-1 ] 





;--------------------------------------- Total velocity ---------------------------------------



    d0_brd_ind  = ind( ( D .ge. 1e-6 ) .and. ( D .lt. 0.006 ) )
    d0_brd      = new( dimsizes( brd_val_1 ) + dimsizes( brd_val_2 ) , "double" )
    d0_brd      = D(d0_brd_ind)


    brd_inf     = new( dimsizes( brd_1_inf ) + dimsizes( brd_2_inf ) , "double" )
    brd_inf(0:dimsizes(brd_1_inf)-1)                                        = brd_1_inf
    brd_inf(dimsizes(brd_1_inf):dimsizes(brd_1_inf)+dimsizes(brd_2_inf)-1 ) = brd_2_inf



;--------------------------------- Various diameter -----------------------------------------


    V_inf_77  = new( nbin , "double" )


    if ( any( D .lt. 1e-6 ) .and. any( D .ge. 0.006 ) ) then


        D_min = max(ind( D .lt. 1e-6  ) )
        D_max = min(ind( D .ge. 0.006 ) )

        V_inf_77(0:D_min)           = 0
        V_inf_77(D_max:nbin-1)      = brd_2_inf_max
        V_inf_77(D_min+1:D_max-1)   = brd_inf


    else if ( any( D .lt. 1e-6 ) .and. all( D .lt. 0.006 ) ) then


        D_min = max(ind( D .lt. 1e-6  ) )
        D_max = nbin-1

        V_inf_77(0:D_min)           = 0
        V_inf_77(D_min+1:D_max)     = brd_inf


    else if ( all( D .ge. 1e-6 ) .and. any( D .ge. 0.006 ) ) then


        D_min = 0
        D_max = min(ind( D .ge. 0.006 ) )

        V_inf_77(D_min:D_max-1) = brd_inf
        V_inf_77(D_max:nbin-1)  = brd_2_inf_max


    else


        V_inf_77 = brd_2_inf_max


    end if
    end if
    end if


    wks = gsn_open_wks("x11" , "brd" )

    res                 = True
    res@gsnLeftString   = "Pressure : "+p+" hPa"
    res@gsnRightString  = "Teamperature : "+(T-273.15)+"~S~o~N~C"
    res@tmXTOn          = False
    res@tmYROn          = False
    res@trXMinF         = 0
    res@trXMaxF         = 0.009
    res@tiYAxisString   = "Terminal Velocity [ m/s ] "
    res@tiXAxisString   = "Diameter [ mm ] "
    res@tmXBMode        = "Explicit"
    res@tmXBValues      = (/ 0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009/)
    res@tmXBLabels      = (/ 0.001,0.002,0.003,0.004,0.005,0.006,0.007,0.008,0.009/)*1000


    plot = gsn_csm_xy( wks, D, V_inf_77, res)





