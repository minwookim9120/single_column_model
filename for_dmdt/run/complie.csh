#!/bin/csh -f 

# [ DO NOT CHANGE ] ==================================================
set rundir=`pwd`
set srcdir=../src
set blddir=../src/bld

rm -f ${blddir}/*.o ${blddir}/*.mod 
# rm -f ${srcdir}/mod_global.mod  ${srcdir}/Mod_global.o  ${srcdir}/mod_read.mod  M${srcdir}/od_read.o 
# rm -f ${srcdir}/Mod_init.o ${srcdir}/mod_init.mod ${srcdir}/mod_integration.mod ${srcdir}/Mod_integration.o 
# rm -f ${srcdir}/mod_dyn_driver.mod ${srcdir}/Mod_dyn_driver.o ${srcdir}/mod_write.mod ${srcdir}/Mod_write.o
# rm -f ${srcdir}/mod_phys_driver.mod ${srcdir}/Mod_dyn_driver.o  ${srcdir}/mod_const.mod 
rm -f ${srcdir}/output.nc main.exe 
# rm -f ${srcdir}/Mod_init_driver.o ${srcdir}/Mod_init_driver.mod ${srcdir}/Mod_idealcase.o
# rm -f ${srcdir}/Mod_idealcase.mod ${srcdir}/Mod_realcase.o ${srcdir}/Mod_realcase.mod 
# rm -f ${srcdir}/Mod_distribution.mod  ${srcdir}/Mod_distribution.o   

# ifort -r8 ${srcdir}/libs/Mod_global.f90          \
ifort -r8 ${srcdir}/libs/Mod_global.f90          \
          ${srcdir}/libs/Mod_const.f90           \
          ${srcdir}/libs/Mod_read.f90            \
          ${srcdir}/init/Mod_idealcase.f90       \
          ${srcdir}/init/Mod_realcase.f90        \
          ${srcdir}/init/Mod_distribution.f90    \
          ${srcdir}/init/Mod_init_driver.f90     \
          ${srcdir}/dyn/Mod_dyn_driver.f90       \
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
