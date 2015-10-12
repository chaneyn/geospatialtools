import os
#os.system('rm -rf terrain_tools.pyf')
os.system('rm -rf terrain_tools.so.dSYM')
#subroutines
subroutines = 'calculate_d8_acc \
              calculate_mfd_acc \
              calculate_channels'

#Create library
cmd = 'f2py -c only: %s : -m terrain_tools terrain_tools.f90 -lgomp --fcompiler=gnu95 --f90flags="-w -fopenmp -O3"' % subroutines
print cmd
os.system(cmd)

#Move to the previos directory
os.system('mv terrain_tools.so ../libraries/.')
