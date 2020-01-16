import setuptools

def configuration(parent_package='', top_path=None):

    from numpy.distutils.misc_util import Configuration
    from numpy.distutils.core import Extension

    config = Configuration('geospatialtools', parent_package, top_path)

    config.add_extension('terrain_tools_fortran',
                         sources=['src/planchon_2001.f90','src/terrain_tools.f90'],
                         extra_f90_compile_args = ['-fPIC','-lgomp','-Wall','-pedantic','-fopenmp','-O3']
                        ),

    config.add_extension('upscaling_tools_fortran',
                         sources=['src/upscaling_tools.f90'],
                         extra_f90_compile_args = ['-fPIC','-lgomp','-Wall','-pedantic','-fopenmp','-O3']
                        ),

    config.add_subpackage('',subpackage_path='libraries')
                        

    return config

if __name__ == '__main__':
    from numpy.distutils.core import setup
    setup(**configuration(top_path='').todict())
