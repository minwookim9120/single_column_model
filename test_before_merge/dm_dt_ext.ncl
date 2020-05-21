;================================ Define Radius of droplet ==================================


;-------------------------------- parameter list --------------------------------------------


        nbin    = 100           ; drop_column_num

        Rd      = 287.          ; Rd = 287 J/(kg*K); J = N*m =(kg*m/s^2)*m
        Cp      = 1005.         ; Cp = 1005.J/(kg*K)
        g       = 9.8           ; unit = kg * m * s^-2
        pi      = 4*atan(1.0)   ; pi = 3.141591
        Ps      = 1013          ; hPa
        rho     = 1000          ; kg/m^3
        nc      = 1.0e+8        ; 1.0e+8
        qc      = 0.002         ; #/m^3
        r0      = 1e-5          ; meter

        nz      = 100           ; unit = number of vertical layer
     z_top      = 10000.0       ; unit = m
       dzr      = 1.05          ; You need this optionsm when you choose geometric sequence

courant_number  = 1.0

integrated_time = 600

       Tc       = 0.            ; temperature [C]
       T        = Tc + 273.15   ; temperature [K]

       Dv       = 0.0000221     ; vapor's diffusion coefficient [m^2 * s^-1]
       Rv       = 461.5         ; gas constant of water vapor   [J * Kg^-1 * K^-1]
       Ka       = 0.024         ; thermal conductivity          [J m^-1 s^-1 K^-1]
       S        = 0.01          ; supersaturation             
       L        = 2501000.      ; latent heat at 0C             [J kg^-1]



;------------------------------------------ option ------------------------------------------


    dist_opt = 2    ; 1) Log-normal distribution    2) Gamma distribution 

    dz_opt   = 2    ; 1) arithmetic progression     2) geometric sequence


    draw_opt = 1    ; 1) X Axis : N     2) X Axis : N       3) X Axis : N        4) X Axis : r 
                    ;    Y Axis : m        Y Axis : dm/dr      Y Axis : dm/dlnr     Y Axis : N 
                                       ; Do not work!!      ; Do not work!!

                    ; 5) X Axis : r     6) X Axis : r       7) X Axis : r        8) X Axis : r
                    ;    Y Axis : n * dr   Y Axis : dN/dlnr    Y Axis : dN/dr       Y Axis : dm/dlnr

                    
                    ; 9) X Axis : t
                    ;    Y Axis : r
                    




;-------------------------------- Calculate radius of droplet -------------------------------


n_nr    = nbin-1   ; drop_column_num -1

rbmin   = 1e-6       ; min boundary r [m]
rbmax   = 1e-2       ; max boundary r [m]

rb   = new( nbin  , "double" ) 
r    = new( n_nr  , "double" ) 
dr   = new( n_nr  , "double" )
dlnr = new( n_nr  , "double" )

rb(0)      = rbmin  ; boundary of drop diameter
rb(nbin-1) = rbmax

ratio   = ( rbmax / rbmin )^(1. / (nbin-1) )

    do ir   = 1, dimsizes(rb)-2 

        rb(ir)  = rb(0) * ( ratio^(ir) )

    end do

    do ir   = 0, dimsizes(r)-1
 
        r(ir)   =      ( rb(ir+1) + rb(ir) ) / 2.
       dr(ir)   =      ( rb(ir+1) - rb(ir) ) 
     dlnr(ir)   = ( log(rb(ir+1)) - log(rb(ir)) )

    end do






;---------------------------------------- Distribution --------------------------------------


    nr      = new( dimsizes(r) , "double" )

if ( dist_opt .eq. 1 ) then

    print("Log-normal distribution")
    
    mu      = log(r0)

    sigma   = sqrt( ( 2. / 9. ) * log( qc / ( nc * rho * ( 4. / 3. ) * pi * ( (r0)^3.) ) ) )


    do  ir = 0, dimsizes(r)-1


        nr(ir)   = ( nc / ( ( sqrt( 2. * pi ) ) * sigma * r(ir) ) ) * exp( ( -1. * ( log(r(ir)) - mu )^2. ) / ( 2. * ( sigma )^2  ) ) * dr(ir) 
                                                                        
;    if (any(isnan_ieee(nr))) then               ; Avoid NaN values.
;        replace_ieeenan (nr, nr@_FillValue, 0)
;    end if

    end do



    print( "log-normal dist. in model (nc) =  "+sum(nr) )
    print( "log-normal dist. const.   (nc) =  "+nc )
    print( "                     accuracy  =  "+( (sum(nr)/nc) * 100.)+"%" )


else if ( dist_opt .eq. 2 ) then


    print("Gamma distribution")

    D       = 2*r

    dD      = 2*dr

    mu      = ( 1.0e+9 / nc ) + 2.

    lambda  = ( ( nc / qc ) * ( gamma( mu + 4. ) / gamma( mu + 1. ) ) * pi * ( 4. / 3. ) * rho )^( 1. / 3. )


    do ir   = 0, dimsizes(r)-1

        nr(ir)  = ( nc / gamma( mu + 1. ) ) * lambda * ( ( lambda * D(ir) )^(mu) ) * exp( -1. * lambda * D(ir) ) * dD(ir)
                

;    if (any(isnan_ieee(nr))) then               ; Avoid NaN values.
;        replace_ieeenan (nr, nr@_FillValue, 0)
;    end if


    end do


    print("gamma dist. in model (nc) =  "+sum(nr) )
    print("gamma dist. const.   (nc) =  "+nc)
    print("                accuracy  =  "+( (sum(nr)/nc) * 100.)+"%")


end if
end if


;---------------------------- Cal. vertical coordinate ----------------------------------------


dz  = new( nz, "double" )


if ( dz_opt .eq. 1 ) then

    dz(:)   = z_top / nz

else if ( dz_opt .eq. 2 ) then

    dz(0)   = ( ( dzr - 1) * z_top ) / ( ( dzr^nz ) -1 )
    
    do i = 1, nz-1

        dz(i)   = dz(i-1) * dzr

    end do

end if
end if


;--------------------------- Cal. dt regarded CFL condition -----------------------------------


CFL     = new( nz, "double" )

w       = new( nz, "double" )

w(:)    = 1.0

    do i    = 0, nz-1

        CFL(i) = where( w(i) .ne. 0. , courant_number * ( dz(i) / abs(w(i))  ) , max(CFL) )

    end do

dt  = toint(min(CFL))
nt  = toint(integrated_time/dt)

    if ( nt*dt .ne. integrated_time ) then

        print("  ")
        print("********WARMING")
        print("Calculated Total integrated time is different from the namelist integrated_time")
        print( "Total integrated time     =  "+nt*dt)
        print( "Namelist integrated_time  =  "+integrated_time)
        print("********")
        print("  ")

    end if


print("dt = "+dt)
print("nt = "+nt)

;----------------------------------------- Graph ----------------------------------------------


es      = new( (/nt, n_nr/) , "double" )

dm      = new( (/nt, n_nr/) , "double" )

mass    = new( n_nr, "double" ) 
mass    = rho * 4. / 3. * pi * r^3

m_dm      = new( (/nt+1, n_nr/) , "double" )
m_dm(0,:) = mass


n_dr    = new( n_nr, "double" )

dm_dr   = new( (/nt, n_nr/) , "double" )

dm_dlnr = new( (/nt, n_nr/) , "double" )



        res               = True        ; Basic setting
        res@vpWidthF      = 1.2
        res@vpHeightF     = 0.6
        res@gsnMaximize   = True
       

;-------------------------------------- m N ---------------------------------------------

        if ( draw_opt .eq. 1 ) then

                wks   = gsn_open_wks("png", "m_dm_0")

                res@gsnLeftString  = "time : 0"
                res@trYMaxF        = 1e-9
                res@trYMinF        = 1e-14
;                res@trXMaxF        = 1.4e7
                res@trYLog         = True
;                res@trXLog         = True
                res@tiYAxisString  = "m ( logscale ) "
                res@tiXAxisString  = "N "
;                res@tiXAxisString  = "N ( logscale ) "
;                res@xyComputeXMin  = True
;                res@xyComputeYMin  = True

                p   = gsn_csm_xy( wks, nr, m_dm(0,:), res )


        do tim   =   0, nt-1

        
                do ra  = 0, n_nr-1

                   es(tim, ra)    = 611.2*exp( (17.67 * Tc ) / (  Tc  + 243.5 ) )  ; saturated water vapor pressure  
                                                                                    ; es [Pa] , Tc [C] 

                   dm(tim, ra)    = 4. * pi * r(ra) * 1.  / ( ( ( Rv * T )  / ( Dv * es(tim,ra) ) ) + ( ( L / ( Rv * T ) -1 ) * ( L / ( Ka * T ) ) ) ) * S * dt      ; dm@units = "m-2 s-1"\

        
                   m_dm(tim+1,ra)     = m_dm(tim,ra) + dm(tim, ra)


                end do

                      wks   = gsn_open_wks("png", "m_dm_"+(tim+1) )

                      res@gsnLeftString  = "time : "+(tim+1)
                      res@trYMaxF        = 1e-9
                      res@trYMinF        = 1e-14
;                      res@trXMaxF        = 1.4e7
                      res@trYLog         = True
;                     res@trXLog         = True
                      res@tiYAxisString  = "m ( logscale ) "
                      res@tiXAxisString  = "N "
;                      res@tiXAxisString  = "N ( logscale ) "
;                      res@xyComputeXMin  = True
;                      res@xyComputeYMin  = True

   
                p   = gsn_csm_xy( wks, nr, m_dm(tim+1,:), res )
;                p   = gsn_csm_xy( wks, dm(0,:), nr, res )



         end do  


         end if

;---------------------------------------- m ---------------------------------------------------

        do tim   =   0, nt-1


                do ra  = 0, n_nr-1

                es(tim, ra)    = 611.2*exp( (17.67 * Tc ) / (  Tc  + 243.5 ) )  ; saturated water vapor pressure  
                                                                            ; es [Pa] , Tc [C] 

                   dm(tim, ra)    = 4. * pi * r(ra) * 1.  / ( ( ( Rv * T )  / ( Dv * es(tim,ra) ) ) + ( ( L / ( Rv * T ) -1 ) * ( L / ( Ka * T ) ) ) ) * S * dt      ; dm@units = "m-2 s-1"\


                m_dm(tim+1,ra)     = m_dm(tim,ra) + dm(tim, ra)


                end do

         end do   
