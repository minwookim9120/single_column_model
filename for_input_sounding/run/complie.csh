#!/bin/csh -f 

# [ DO NOT CHANGE ] ==================================================
set rundir=`pwd`
set srcdir=../src
set blddir=../src/bld

rm -f ${blddir}/*.o ${blddir}/*.mod 
rm -f ${srcdir}/output/output.nc main.exe 

ifort -r8 ${srcdir}/libs/Mod_global.f90          \
          ${srcdir}/libs/Mod_const.f90           \
          ${srcdir}/libs/Mod_read.f90            \
          ${srcdir}/init/Mod_idealcase.f90       \
          ${srcdir}/init/Mod_realcase.f90        \
          ${srcdir}/init/Mod_distribution.f90    \
          ${srcdir}/init/Mod_init_driver.f90     \
          ${srcdir}/dyn/Mod_dyn_driver.f90       \
          ${srcdir}/phys/Mod_drop_growth.f90     \
          ${srcdir}/phys/Mod_phys_driver.f90     \
          ${srcdir}/main/Mod_integration.f90     \
          ${srcdir}/libs/Mod_write.f90           \
          ${srcdir}/main/main.f90                \
      -o  main.exe                               \
      -L/$NETCDF/lib -I/$NETCDF/include -lnetcdf -lnetcdff

mv *.o *.mod ${blddir}

./main.exe 
#2>/dev/null

exit 0
#=====================================================================

