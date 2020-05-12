;================================ Define Radius of droplet ==================================


;-------------------------------- parameter list --------------------------------------------


        nbin    = 300           ; drop_column_num

        Rd      = 287.          ; Rd = 287 J/(kg*K); J = N*m =(kg*m/s^2)*m
        Cp      = 1005.         ; Cp = 1005.J/(kg*K)
        g       = 9.8           ; unit = kg * m * s^-2
        pi      = 4*atan(1.0)   ; pi = 3.141591
        Ps      = 1013          ; hPa
        rho     = 1000          ; kg/m^3


;------------------------------- nc opt.-----------------------------------------------------


        nc1      = 1.0e+8        ; 1.0e+8
        nc2      = 1.0e+8        ; 1.0e+8
        nc3      = 1.0e+8        ; 1.0e+8
        
        qc1      = 0.002        ; #/m^3
        qc2      = 0.001         ; #/m^3
        qc3      = 0.0005        ; #/m^3


;--------------------------------------------------------------------------------------------

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


    dist_opt = 1    ; 1) Log-normal distribution    2) Gamma distribution 

    dz_opt   = 2    ; 1) arithmetic progression     2) geometric sequence


    draw_opt = 6   ; 1) X Axis : N     2) X Axis : N       3) X Axis : N        4) X Axis : r 
                    ;    Y Axis : m        Y Axis : dm/dr      Y Axis : dm/dlnr     Y Axis : N 
                                       ; Do not work!!      ; Do not work!!

                    ; 5) X Axis : r     6) X Axis : r       7) X Axis : r        8) X Axis : r
                    ;    Y Axis : n * dr   Y Axis : dN/dlnr    Y Axis : dN/dr       Y Axis : dm/dlnr

                    
                    ; 9) X Axis : t    10) X Axis : m      11) X Axis : r       12) X Axis : r (growth ver.)
                    ;    Y Axis : r        Y Axis : dN/dlnr    Y Axis : dm/dr       Y Axis : dN/dlnr
                                                                                ; Do not work!!                    

                    ; 13) X Axis : r
                    ;     Y Axis : N*dm




;-------------------------------- Calculate radius of droplet -------------------------------

;n_nr    = nbin-1   ; drop_column_num -1

rbmin   = 1e-6       ; min boundary r [m]
rbmax   = 1e-2       ; max boundary r [m]

rb   = new( nbin+1  , "double" ) 
r    = new( nbin  , "double" ) 
dr   = new( nbin  , "double" )
dlnr = new( nbin  , "double" )
mb = new( nbin+1  , "double" )
m = new( nbin  , "double" )
dm = new( nbin  , "double" )

rb(0)      = rbmin  ; boundary of drop diameter
rb(nbin-1) = rbmax

ratio   = ( rbmax / rbmin )^(1. / nbin )

    do ir   = 1, (nbin-1)

        rb(ir)  = rb(0) * ( ratio^(ir) )

    end do

   do ir = 0, nbin

      mb(ir) = rho * (4. / 3.) * pi * (rb(ir))^3
   end do

    do ir   = 0, (nbin-1)
      m(ir) = (mb(ir) + mb(ir+1)) / 2.
        r(ir)   = ((3. * m(ir)) / (4 * rho * pi))^(1. / 3.)
       dr(ir)   =      ( rb(ir+1) - rb(ir) ) 
      dm(ir) = mb(ir+1) - mb(ir)
     dlnr(ir)   = ( log(rb(ir+1)) - log(rb(ir)) )

    end do



;---------------------------------------- Distribution --------------------------------------


    nr1      = new( nbin , "double" )
    nr2      = new( nbin , "double" )
    nr3      = new( nbin , "double" )

    N_lnd1      = new( nbin , "double" )
    N_lnd2      = new( nbin , "double" )
    N_lnd3      = new( nbin , "double" )

if ( dist_opt .eq. 1 ) then

    print("Log-normal distribution")
    
    mu      = log(r0)

    sigma1   = sqrt( ( 2. / 9. ) * log( qc1 / ( nc1 * rho * ( 4. / 3. ) * pi * ( (r0)^3.) ) ) )
    sigma2   = sqrt( ( 2. / 9. ) * log( qc2 / ( nc2 * rho * ( 4. / 3. ) * pi * ( (r0)^3.) ) ) )
    sigma3   = sqrt( ( 2. / 9. ) * log( qc3 / ( nc3 * rho * ( 4. / 3. ) * pi * ( (r0)^3.) ) ) )


    do  ir = 0, (nbin-1)


        N_lnd1(ir)   = ( nc1 / ( ( sqrt( 2. * pi ) ) * sigma1 * r(ir) ) ) * exp( ( -1. * ( log(r(ir)) - mu )^2. ) / ( 2. * ( sigma1 )^2  ) ) * dr(ir) 
        N_lnd2(ir)   = ( nc2 / ( ( sqrt( 2. * pi ) ) * sigma2 * r(ir) ) ) * exp( ( -1. * ( log(r(ir)) - mu )^2. ) / ( 2. * ( sigma2 )^2  ) ) * dr(ir) 
        N_lnd3(ir)   = ( nc3 / ( ( sqrt( 2. * pi ) ) * sigma3 * r(ir) ) ) * exp( ( -1. * ( log(r(ir)) - mu )^2. ) / ( 2. * ( sigma3 )^2  ) ) * dr(ir) 
                                                                        
;    if (any(isnan_ieee(nr))) then               ; Avoid NaN values.
;        replace_ieeenan (nr, nr@_FillValue, 0)
;    end if

    end do



    print( "log-normal dist. in model (nc1) =  "+sum(N_lnd1) )
    print( "log-normal dist. const.   (nc1) =  "+nc1 )
    print( "                     accuracy  =  "+( (sum(N_lnd1)/nc1) * 100.)+"%" )

    print( "log-normal dist. in model (nc2) =  "+sum(N_lnd2) )
    print( "log-normal dist. const.   (nc2) =  "+nc2 )
    print( "                     accuracy  =  "+( (sum(N_lnd2)/nc2) * 100.)+"%" )

    print( "log-normal dist. in model (nc3) =  "+sum(N_lnd3) )
    print( "log-normal dist. const.   (nc3) =  "+nc3 )
    print( "                     accuracy  =  "+( (sum(N_lnd3)/nc3) * 100.)+"%" )

    print( "log-normal dist. in model (qc1) =  "+sum(m*N_lnd1) )
    print( "log-normal dist. const.   (qc1) =  "+qc1 )
    print( "                     accuracy  =  "+( (sum(m*N_lnd1)/qc1) * 100.)+"%" )

    print( "log-normal dist. in model (qc2) =  "+sum(m*N_lnd2) )
    print( "log-normal dist. const.   (qc2) =  "+qc2 )
    print( "                     accuracy  =  "+( (sum(m*N_lnd2)/qc2) * 100.)+"%" )

    print( "log-normal dist. in model (qc3) =  "+sum(m*N_lnd3) )
    print( "log-normal dist. const.   (qc3) =  "+qc3 )
    print( "                     accuracy  =  "+( (sum(m*N_lnd3)/qc3) * 100.)+"%" )

else if ( dist_opt .eq. 2 ) then


    print("Gamma distribution")

    D       = 2*r

    dD      = 2*dr

    mu1      = ( 1.0e+9 / nc1 ) + 2.
    mu2      = ( 1.0e+9 / nc2 ) + 2.
    mu3      = ( 1.0e+9 / nc3 ) + 2.

    lambda1  = ( ( nc1 / qc1 ) * ( gamma( mu1 + 4. ) / gamma( mu1 + 1. ) ) * pi * ( 4. / 3. ) * rho )^( 1. / 3. )
    lambda2  = ( ( nc2 / qc2 ) * ( gamma( mu2 + 4. ) / gamma( mu2 + 1. ) ) * pi * ( 4. / 3. ) * rho )^( 1. / 3. )
    lambda3  = ( ( nc3 / qc3 ) * ( gamma( mu3 + 4. ) / gamma( mu3 + 1. ) ) * pi * ( 4. / 3. ) * rho )^( 1. / 3. )


    do ir   = 0, nbin-1

        nr1(ir)  = ( nc1 / gamma( mu1 + 1. ) ) * lambda1 * ( ( lambda1 * D(ir) )^(mu1) ) * exp( -1. * lambda1 * D(ir) ) * dD(ir)
        nr2(ir)  = ( nc2 / gamma( mu2 + 1. ) ) * lambda2 * ( ( lambda2 * D(ir) )^(mu2) ) * exp( -1. * lambda2 * D(ir) ) * dD(ir)
        nr3(ir)  = ( nc3 / gamma( mu3 + 1. ) ) * lambda3 * ( ( lambda3 * D(ir) )^(mu3) ) * exp( -1. * lambda3 * D(ir) ) * dD(ir)
                

;    if (any(isnan_ieee(nr))) then               ; Avoid NaN values.
;        replace_ieeenan (nr, nr@_FillValue, 0)
;    end if


    end do
    print(" ")
    print("gamma dist. in model (nc1) =  "+sum(nr1) )
    print("gamma dist. const.   (nc1) =  "+nc1)
    print("                accuracy  =  "+( (sum(nr1)/nc1) * 100.)+"%")

    print("gamma dist. in model (nc2) =  "+sum(nr2) )
    print("gamma dist. const.   (nc2) =  "+nc2)
    print("                accuracy  =  "+( (sum(nr2)/nc2) * 100.)+"%")

    print("gamma dist. in model (nc3) =  "+sum(nr3) )
    print("gamma dist. const.   (nc3) =  "+nc3)
    print("                accuracy  =  "+( (sum(nr3)/nc3) * 100.)+"%")

    print("gamma dist. in model (qc1) =  "+sum(m*nr1) )
    print("gamma dist. const.   (qc1) =  "+qc1)
    print("                accuracy  =  "+( (sum(m*nr1)/qc1) * 100.)+"%")

    print("gamma dist. in model (qc2) =  "+sum(m*nr2) )
    print("gamma dist. const.   (qc2) =  "+qc2)
    print("                accuracy  =  "+( (sum(m*nr2)/qc2) * 100.)+"%")
    
    print("gamma dist. in model (qc3) =  "+sum(m*nr3) )
    print("gamma dist. const.   (qc3) =  "+qc3)
    print("                accuracy  =  "+( (sum(m*nr3)/qc3) * 100.)+"%")
    print(" ")

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


es      = new( (/nt, nbin/) , "double" )

dm2      = new( (/nt, nbin/) , "double" )

mass    = new( nbin, "double" ) 
mass    = rho * 4. / 3. * pi * r^3

m_dm      = new( (/nt+1, nbin/) , "double" )
m_dm(0,:) = mass


n_dr    = new( nbin, "double" )

dm_dr   = new( (/nt, nbin/) , "double" )

dm_dlnr = new( (/nt+1, nbin/) , "double" )



        res               = True        ; Basic setting
        res@vpWidthF      = 1.2
        res@vpHeightF     = 0.6
        res@gsnMaximize   = True
       

;-------------------------------------- m N ---------------------------------------------

        if ( draw_opt .eq. 1 ) then

                wks   = gsn_open_wks("x11", "m_dm_0")


                res@gsnLeftString  = "time : 0"
                res@trXMaxF        = 1e-9
                res@trXMinF        = 1e-14
;                res@trXMaxF        = 1.4e7
                res@trXLog         = True
;                res@trXLog         = True
                res@tiXAxisString  = "m ( logscale ) "
                res@tiYAxisString  = "N "
;                res@tiXAxisString  = "N ( logscale ) "
;                res@xyComputeXMin  = True
;                res@xyComputeYMin  = True





                p   = gsn_csm_xy( wks, m_dm(0,:), nr, res )


        do tim   =   0, nt-1

        
                do ra  = 0, n_nr-1

                   es(tim, ra)    = 611.2*exp( (17.67 * Tc ) / (  Tc  + 243.5 ) )  ; saturated water vapor pressure  
                                                                                    ; es [Pa] , Tc [C] 

                   dm2(tim, ra)    = 4. * pi * r(ra) * 1.  / ( ( ( Rv * T )  / ( Dv * es(tim,ra) ) ) + ( ( ( L / ( Rv * T ) ) -1 ) * ( L / ( Ka * T ) ) ) ) * S * dt      ; dm@units = "m-2 s-1"\

        
                   m_dm(tim+1,ra)     = m_dm(tim,ra) + dm2(tim, ra)


                end do

                      wks   = gsn_open_wks("x11", "m_dm_"+(tim+1) )

   
                res@gsnLeftString  = "time : "+(tim+1)
                res@trXMaxF        = 1e-9
                res@trXMinF        = 1e-14
;                res@trXMaxF        = 1.4e7
                res@trXLog         = True
;                res@trXLog         = True
                res@tiXAxisString  = "m ( logscale ) "
                res@tiYAxisString  = "N "
;                res@tiXAxisString  = "N ( logscale ) "
;                res@xyComputeXMin  = True
;                res@xyComputeYMin  = True

                p   = gsn_csm_xy( wks, m_dm(tim+1,:), nr, res )




         end do  


         end if

;---------------------------------------- m ---------------------------------------------------

        do tim   =   0, nt-1


                do ra  = 0, (nbin-1)

                es(tim, ra)    = 611.2*exp( (17.67 * Tc ) / (  Tc  + 243.5 ) )  ; saturated water vapor pressure  
                                                                            ; es [Pa] , Tc [C] 

                   dm2(tim, ra)    = 4. * pi * r(ra) * 1.  / ( ( ( Rv * T )  / ( Dv * es(tim,ra) ) ) + ( ( ( L / ( Rv * T ) ) -1 ) * ( L / ( Ka * T ) ) ) ) * S * dt      ; dm@units = "m-2 s-1"\


                m_dm(tim+1,ra)     = m_dm(tim,ra) + dm2(tim, ra)


                end do

         end do   

         
;------------------------------------------ dm/dr N ----------------------------------------------


             if ( draw_opt .eq. 2 ) then

                do tim  = 0, nt-1

                    dm_dr(tim,:)  = dm(tim,:) / dr(:)

                    wks   = gsn_open_wks("x11", "dm_dr")

                    res@gsnLeftString = ""
                    res@tiYAxisString = "dm/dr"
                    res@tiXAxisString = "N"

                    p     = gsn_csm_xy( wks, nr, dm_dr(tim,:), res )

                end do

             end if

;------------------------------------------ dm/dlnr nr ----------------------------------------------

             if ( draw_opt .eq. 3 ) then


               do tim = 0, nt-1

                    dm_dlnr(tim,:) = m_dm(tim,:) / dlnr(:)

;                    wks   = gsn_open_wks("x11" , "dm_dlnr_"+tim)
;
;                    res@gsnLeftString = ""
;                    res@tiYAxisString = "dm/dlnr"
;                    res@tiXAxisString = "N"
;                    res@trYMaxF       = 1e-6
;                    res@trYMinF       = 1e-17
;                    res@trXMaxF       = 1.4e7
;                    res@trYLog        = True
;                    res@trXLog        = True
;
;
;                p   = gsn_csm_xy( wks, nr, dm_dlnr(tim, :), res)
                
               end do

            end if


;------------------------------------------ dN/dr r----------------------------------------------


            if ( draw_opt .eq. 4 ) then


                  wks   = gsn_open_wks("x11", "N_r" )

                  res@gsnLeftString = ""
                  res@tiYAxisString = "N"
                  res@tiXAxisString = "r"
                  res@trXLog        = True
                  res@xyComputeXMin = True
                  res@xyComputeYMin = True

                  p = gsn_csm_xy( wks, r, nr , res)

            end if


;------------------------------------------ N*dr r----------------------------------------------


            if ( draw_opt .eq. 5 ) then


                  n_dr = nr * dr

                  wks   = gsn_open_wks("x11", "N*dr")

                  res@gsnLeftString = ""
                  res@tiYAxisString = "N*dr"
                  res@tiXAxisString = "r"
                  res@trXLog        = True
                  res@xyComputeXMin = True
                  res@xyComputeYMin = True

                  p = gsn_csm_xy( wks, r, n_dr, res)

           end if 

;------------------------------------------ dN/dlnr r----------------------------------------------

fn = (/"gd_nc3", "gd_qc3", "nld_nc3", "nld_qc3"/)

            if ( draw_opt .eq. 6 ) then


                label_n = (/"N~B~c~N~ = " + nc1 + " / m~S~3~N~", "N~B~c~N~ = " + nc2 + " / m~S~3~N~", "N~B~c~N~ = " + nc3 + " / m~S~3~N~"/)
                label_q = (/"q~B~c~N~ = " + qc1 + " kg / kg", "q~B~c~N~ = " + qc2 + " kg / kg", "q~B~c~N~ = " + qc3 + " kg / kg"/)
                    ;dN_dlnr1 = nr1 / dlnr
                    ;dN_dlnr2 = nr2 / dlnr
                    ;dN_dlnr3 = nr3 / dlnr

                    dN_dlnr1 = N_lnd1 / dlnr
                    dN_dlnr2 = N_lnd2 / dlnr
                    dN_dlnr3 = N_lnd3 / dlnr

                    wks   = gsn_open_wks("png" , fn(3))

                    res@gsnMaximize    = True
                    ;res@gsnRightString = "q~B~c~N~ = " + qc1 + " kg / kg"
                    res@gsnRightString = "N~B~c~N~ = " + nc1 + " m~S~3~N~"
                    ;res@tiMainString   = "Gamma Dist."
                    res@tiMainString   = "Log-Normal Dist."
                    res@tiYAxisString  = "dN/dlnr "
                    res@tiXAxisString  = "r (logscale) [m]"
                    res@trXLog         = True
                    res@tmXTOn         = False
                    res@tmYROn         = False
                    res@xyLineColors   = (/"red", "blue", "green"/)
                    res@xyLineThicknessF = 2.0
                    res@xyMonoDashPattern = 1
                    res@xyExplicitLegendLabels = label_q
                    res@pmLegendDisplayMode    = "Always"
                    res@pmLegendOrthogonalPosF = -1.3
                    res@pmLegendParallelPosF   = 0.8
                    res@pmLegendWidthF  = 0.4
                    res@pmLegendHeightF = 0.2
                    res@lgPerimOn       = False
                    res@lgLabelFontHeightF = 0.02




                ;p   = gsn_csm_xy( wks, r, (/dN_dlnr1, dN_dlnr2, dN_dlnr3/), res)
                p   = gsn_csm_xy( wks, r, (/dN_dlnr1, dN_dlnr2, dN_dlnr3/), res)

            system("convert -alpha off -background white -density 300 -trim " + fn(1) + ".png " + fn(1) + ".png")
                r2 = True




                

            end if

                
;----------------------------------------- N/dr r ---------------------------------------------------


            if ( draw_opt .eq. 7 ) then


                    dN_dr  = nr / dr

                    wks   = gsn_open_wks("x11" , "dN_dr")

                    res@gsnLeftString = ""
                    res@tiYAxisString = "dN/dr"
                  ; res@tiXAxisString = "r (logscale)"
                    res@trXLog        = True                                 
                    res@trYLog        = True
                    res@trXMinF       = 0.
                    res@trYMinF       = 0.
        
                  ;  res@xyComputeXMin = True
                  ;  res@xyComputeYMin = True

                       
                p   = gsn_csm_xy( wks, r, dN_dr, res)


            end if


;----------------------------------------- dm/dlnr r -------------------------------------------------


            if ( draw_opt .eq. 8 ) then

                do tim = 0, nt-1


   ;                dm_dlnr(tim,:) =   m_dm(tim,:) / dlnr
                   dm_dlnr(tim,:) =   dm(tim,:) / dlnr


                    wks = gsn_open_wks("x11", "dm_dlnr"+tim )

                    res@gsnLeftString = "time : "+tim
                    res@tiYAxisString = "dm/dlnr"
                    res@tiXAxisString = "r"
                    res@trXLog        = True
          
                    p = gsn_csm_xy(wks, r, dm_dlnr, res)

                end do

             end if

;------------------------------------------ r t ----------------------------------------------------


            if ( draw_opt .eq. 9 ) then

                dr_1 = new( nt, "double")
                dr_2 = new( nt, "double")
                dr_3 = new( nt, "double")

                 r_1 = new( nt+1, "double" )
                 r_2 = new( nt+1, "double" )
                 r_3 = new( nt+1, "double" )

                r_1(0) = r(7)
                r_2(0) = r(22)
                r_3(0) = r(30)

                es := 611.2
            
                T  := 273.15

                do i = 0, nt-1


                    dr_1(i) = 1. / ( rho * ( ( ( Rv * T ) / ( Dv * es ) )  + ( ( L / ( Rv * T ) ) -1 ) * ( L / ( Ka * T ) ) ) ) * ( 1 / ( r_1(i) ) ) * S * dt
                    dr_2(i) = 1. / ( rho * ( ( ( Rv * T ) / ( Dv * es ) )  + ( ( L / ( Rv * T ) ) -1 ) * ( L / ( Ka * T ) ) ) ) * ( 1 / ( r_2(i) ) ) * S * dt
                    dr_3(i) = 1. / ( rho * ( ( ( Rv * T ) / ( Dv * es ) )  + ( ( L / ( Rv * T ) ) -1 ) * ( L / ( Ka * T ) ) ) ) * ( 1 / ( r_3(i) ) ) * S * dt

                    r_1(i+1) = r_1(i) + dr_1(i)
                    r_2(i+1) = r_2(i) + dr_2(i)
                    r_3(i+1) = r_3(i) + dr_3(i)

                   
                end do


                    wks = gsn_open_wks("x11", "r_t")

                    res@trYLog          = True
                    res@trYMaxF         = 1e-4
                    res@tiXAxisString   = "time"
                    res@tiYAxisString   = "r (m)"


                    p = gsn_csm_xy( wks, ispan(0, nt, 1), (/r_1, r_2, r_3/), res)


                end if



;---------------------------------------- dN/dlnr m --------------------------------------


                if ( draw_opt .eq. 10 ) then


                    do tim = 0, nt


                        dN_dlnr = nr / dlnr

                        wks = gsn_open_wks("x11" , "dN_dlnr_m"+tim)

                        res@gsnLeftString = "time : "+tim
                        res@trXLog        = True
                        res@trXMinF       = 1e-14
                        res@trXMaxF       = 1e-9
                        res@tiYAxisString = "dN/dlnr"
                        res@tiXAxisString = "m"
    
                        p   = gsn_csm_xy(wks, m_dm(tim, :), dN_dlnr , res)


                    end do
                        

                end if


;---------------------------------------- dm/dr r ----------------------------------------


                if ( draw_opt .eq. 11 ) then

                m_dr = new( (/nt+1, n_nr /), "double" )


                    do tim = 0, nt

                        m_dr(tim,:) = m_dm(tim,:) / dr

                        wks = gsn_open_wks("x11" , "dm_dr_r"+tim)

                        res@gsnLeftString = "time : "+tim
                        res@trXLog        = True
                  ;      res@trXMinF       = 1e-14
                  ;      res@trXMaxF       = 1e-9
                        res@tiYAxisString = "dm/dr"
                        res@tiXAxisString = "r"
    
                        p   = gsn_csm_xy(wks, r, m_dr , res)

                        

                    end do

                end if


;----------------------------------------- dN/dlnr r(growth ver.)-------------------------


                if ( draw_opt .eq. 12 ) then


                dr_tot = new( (/ nt, n_nr /) , "double")
                 r_tot = new( (/nt+1, n_nr/) , "double" )


                r_tot(0,:) = r(:)


                es := 611.2
            
                T  := 273.15



                do i = 0, nt-1


                    dr_tot(i,:) = 1. / ( rho * ( ( ( Rv * T ) / ( Dv * es ) )  + ( ( L / ( Rv * T ) ) -1 ) * ( L / ( Ka * T ) ) ) ) * ( 1 / ( r_tot(i,:) ) ) * S * dt


                    r_tot(i+1,:) = r_tot(i,:) + dr_tot(i,:)


                   
                end do


                    dlnr_tot = log(r_tot(i+1,:)) - log(r_tot(i+1,:))


                end if


;----------------------------------------- N*dm r --------------------------------------------


                if ( draw_opt .eq. 13 ) then

                
                Ndm_r = new( (/nt, n_nr /), "double" )


                    do tim = 0, nt-1

                        Ndm_r(tim,:) = dm(tim,:) * nr(:)

                        wks = gsn_open_wks("png" , "Ndm_r"+tim)

                        res@gsnLeftString = "time : "+tim
                        res@trXLog        = True
                        res@tiYAxisString = "N*dm"
                        res@tiXAxisString = "r"
    
                        p   = gsn_csm_xy(wks, r, Ndm_r(tim,:) , res)

                        

                    end do

                end if
