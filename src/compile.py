import os
#os.system('rm -rf terrain_tools.pyf')
os.system('rm -rf *.dSYM')
#subroutines
subroutines = 'calculate_d8_acc \
              calculate_mfd_acc \
              calculate_channels \
              delineate_basins \
              delineate_hillslopes \
              calculate_hillslope_properties \
              calculate_depth2channel \
              calculate_hillslopesd8 \
              assign_clusters_to_hillslopes \
              calculate_hru_properties \
              retrieve_basin_properties'

#Create library
cmd = 'f2py -c only: %s : -m terrain_tools_fortran terrain_tools.f90 -lgomp --fcompiler=gnu95 --f90flags="-w -fopenmp -O3"' % subroutines
print cmd
os.system(cmd)

#Move to the previos directory
os.system('mv terrain_tools_fortran.so ../libraries/.')
