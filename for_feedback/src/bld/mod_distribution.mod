	  ·!  ?   k820309              14.0        ôoÏ^                                                                                                           
       ../src/init/Mod_distribution.f90 MOD_DISTRIBUTION                                                    
                                                          
                                                           
       SUB_READ_NC_FILE                   @               @                '                    #VARID    #DZ    #NEXT_DZ    #STAG_DZ    #DT 	   #SFC_DT 
   #TOP_DT    #REF_M    #REF_MB    #REF_NUM    #DIN    #DOUT    #R    #M    #RB    #MB    #NUM    #NEXT_NUM    #DM_DT    #DMB    #DMB_DT    #DROP_DOUT    #MASS_DOUT    #R_DOUT    #VNAME    #AXIS    #DESC    #UNITS                                                                                                                                           
            &                                                                                                 P                 
            &                                                                                                                  
            &                                                                                     	            à                 
            &                                                                                     
            (                
            &                                                                                                 p                
            &                                                                                                 ¸                
            &                                                                                                               	   
            &                                                                                                 H             
   
            &                                                                                                                 
            &                                                                                                 Ø                
            &                   &                                                                                                 8                
            &                   &                                                                                                                 
            &                   &                                                                                                 ø                
            &                   &                                                                                                 X                
            &                   &                                                                                                 ¸                
            &                   &                                                                                                                 
            &                   &                                                                                                 x                
            &                   &                                                                                                 Ø                
            &                   &                                                                                                 8                
            &                   &                                                                                                                 
            &                   &                   &                                                                                                                 
            &                   &                   &                                                                                                                 
            &                   &                   &                                                                                                                                                                        	                                                                     
                                                                                   #         @                                   !                   #SUB_READ_NC_FILE%ALLOCATED "   #SUB_READ_NC_FILE%TRIM #   #SUB_READ_NC_FILE%SUM $   #INPATH %   #INNAME &   #OUT_VAR '   #SLAT (   #ELAT )   #SLON *   #ELON +                                              "     ALLOCATED                                            #     TRIM                                            $     SUM           
                                 %                    1           
                                 &                    1         
                                 '                   
               &                                                     
                                  (                     
                                  )                     
                                  *                     
                                  +                                                        ,     
                 
                 ñhãµøä>        1.0E-5                                             -     
                 
                      @@        1000                                             .     
                 
                     ×A        5.0E+7                                             /     
                 
                 íµ ÷ÆÀ>        2.0E-6                                             0     
                   
                  -DTû!	@        #         @                                   1                	   #SUB_DROP_DISTRIBUTIONS%GAMMA 2   #SUB_DROP_DISTRIBUTIONS%EXP 3   #SUB_DROP_DISTRIBUTIONS%SQRT 4   #SUB_DROP_DISTRIBUTIONS%LOG 5   #DISTRIBUTION_OPTION 6   #DROP_COLUMN_NUM 7   #RBMIN 8   #RBMAX 9   #NR :   #R ;   #RB <   #M =   #MB >                                                                                                                                                                                                          2     GAMMA                                            3     EXP                                            4     SQRT                                            5     LOG           
                                  6                     
                                  7                     
                                  8     
                
                                  9     
               D                                 :                    
     p          5  p        r 7       5  p        r 7                              D                                 ;                    
     p          5  p        r 7       5  p        r 7                              D                                 <                    
     p           5  p        r 7   n                                       1     5  p        r 7   n                                      1                                    D                                 =                    
     p          5  p        r 7       5  p        r 7                              D                                 >                    
 	    p           5  p        r 7   n                                       1     5  p        r 7   n                                      1                                  :      fn#fn    Ú   @   J   MOD_GLOBAL      @   J   MOD_CONST    Z  Q   J  MOD_READ #   «  }      VARINFO+MOD_GLOBAL )   (  H   a   VARINFO%VARID+MOD_GLOBAL &   p     a   VARINFO%DZ+MOD_GLOBAL +        a   VARINFO%NEXT_DZ+MOD_GLOBAL +        a   VARINFO%STAG_DZ+MOD_GLOBAL &   ,     a   VARINFO%DT+MOD_GLOBAL *   À     a   VARINFO%SFC_DT+MOD_GLOBAL *   T     a   VARINFO%TOP_DT+MOD_GLOBAL )   è     a   VARINFO%REF_M+MOD_GLOBAL *   |     a   VARINFO%REF_MB+MOD_GLOBAL +        a   VARINFO%REF_NUM+MOD_GLOBAL '   ¤     a   VARINFO%DIN+MOD_GLOBAL (   8	  ¬   a   VARINFO%DOUT+MOD_GLOBAL %   ä	  ¬   a   VARINFO%R+MOD_GLOBAL %   
  ¬   a   VARINFO%M+MOD_GLOBAL &   <  ¬   a   VARINFO%RB+MOD_GLOBAL &   è  ¬   a   VARINFO%MB+MOD_GLOBAL '     ¬   a   VARINFO%NUM+MOD_GLOBAL ,   @  ¬   a   VARINFO%NEXT_NUM+MOD_GLOBAL )   ì  ¬   a   VARINFO%DM_DT+MOD_GLOBAL '     ¬   a   VARINFO%DMB+MOD_GLOBAL *   D  ¬   a   VARINFO%DMB_DT+MOD_GLOBAL -   ð  Ä   a   VARINFO%DROP_DOUT+MOD_GLOBAL -   ´  Ä   a   VARINFO%MASS_DOUT+MOD_GLOBAL *   x  Ä   a   VARINFO%R_DOUT+MOD_GLOBAL )   <  P   a   VARINFO%VNAME+MOD_GLOBAL (     P   a   VARINFO%AXIS+MOD_GLOBAL (   Ü  P   a   VARINFO%DESC+MOD_GLOBAL )   ,  P   a   VARINFO%UNITS+MOD_GLOBAL *   |  ê       SUB_READ_NC_FILE+MOD_READ >   f  B      SUB_READ_NC_FILE%ALLOCATED+MOD_READ=ALLOCATED 4   ¨  =      SUB_READ_NC_FILE%TRIM+MOD_READ=TRIM 2   å  <      SUB_READ_NC_FILE%SUM+MOD_READ=SUM 1   !  L   a   SUB_READ_NC_FILE%INPATH+MOD_READ 1   m  L   a   SUB_READ_NC_FILE%INNAME+MOD_READ 2   ¹     a   SUB_READ_NC_FILE%OUT_VAR+MOD_READ /   E  @   a   SUB_READ_NC_FILE%SLAT+MOD_READ /     @   a   SUB_READ_NC_FILE%ELAT+MOD_READ /   Å  @   a   SUB_READ_NC_FILE%SLON+MOD_READ /     @   a   SUB_READ_NC_FILE%ELON+MOD_READ    E  v       R0+MOD_CONST    »  t       RHO+MOD_CONST    /  v       NC+MOD_CONST    ¥  v       QC+MOD_CONST      p       PI+MOD_CONST '     Ñ      SUB_DROP_DISTRIBUTIONS -   \  >      SUB_DROP_DISTRIBUTIONS%GAMMA +     <      SUB_DROP_DISTRIBUTIONS%EXP ,   Ö  =      SUB_DROP_DISTRIBUTIONS%SQRT +     <      SUB_DROP_DISTRIBUTIONS%LOG ;   O  @   a   SUB_DROP_DISTRIBUTIONS%DISTRIBUTION_OPTION 7     @   a   SUB_DROP_DISTRIBUTIONS%DROP_COLUMN_NUM -   Ï  @   a   SUB_DROP_DISTRIBUTIONS%RBMIN -     @   a   SUB_DROP_DISTRIBUTIONS%RBMAX *   O  ´   a   SUB_DROP_DISTRIBUTIONS%NR )     ´   a   SUB_DROP_DISTRIBUTIONS%R *   ·  &  a   SUB_DROP_DISTRIBUTIONS%RB )   Ý  ´   a   SUB_DROP_DISTRIBUTIONS%M *      &  a   SUB_DROP_DISTRIBUTIONS%MB 