import os
#Terrain tools
#os.system('rm -rf terrain_tools.pyf')
os.system('rm -rf *.dSYM')
#compile modules
#os.system('gfortran -fPIC -c planchon_2001.f90')
os.system('gfortran -fPIC -c planchon_2001.f90')
#subroutines
subroutines = 'calculate_d8_acc \
              calculate_d8_acc_pfdir \
              calculate_d8_acc_wipoints \
	      calculate_d8_acc_neighbors \
              calculate_mfd_acc \
              neighbr_check_d8 \
              cleanup_hillslopes \
              calculate_channels \
              calculate_channels_wocean \
              calculate_channels_wocean_wprop \
              delineate_basins \
              delineate_hillslopes \
              assign_properties_to_hillslopes \
              calculate_basin_properties \
              calculate_depth2channel \
              calculate_depth2channel_mfd \
              assign_clusters_to_hillslopes \
              calculate_hru_properties \
              retrieve_basin_properties \
              calculate_slope_and_aspect \
              remove_pits \
              remove_pits_planchon \
              polygonize_raster \
              compute_polygon_info \
              gap_fill_hrus'

#Create library
cmd = 'f2py -c only: %s : *.o -m terrain_tools_fortran terrain_tools.f90 -lgomp --fcompiler=gnu95 --f90flags="-Wall -pedantic -fopenmp -O3"' % subroutines
#print cmd
os.system(cmd)

#Move to the previos directory
os.system('mv terrain_tools_fortran.*.so terrain_tools_fortran_v2.so')
#os.system('mv terrain_tools_fortran_v2.so ../libraries/.')

#Upscaling tools
#os.system('rm -rf terrain_tools.pyf')
os.system('rm -rf *.dSYM')
#subroutines
subroutines = 'time_average'

#Create library
cmd = 'f2py -c only: %s : -m upscaling_tools_fortran upscaling_tools.f90 -lgomp --fcompiler=gnu95 --f90flags="-w -fopenmp -O3"' % subroutines
#print cmd
os.system(cmd)

#Move to the previos directory
os.system('mv upscaling_tools_fortran.*.so upscaling_tools_fortran.so')
os.system('mv upscaling_tools_fortran.so ../libraries/.')
