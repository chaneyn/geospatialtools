#First compile the source code
import os
cwd = os.getcwd()
os.chdir('src')
os.system('python compile.py')
os.chdir(cwd)

#Write the libraries
from distutils.core import setup
setup(name='geospatialtools',
      version='1.0',
      package_dir={'geospatialtools': 'libraries'},
      packages=['geospatialtools'],
      package_data={'geospatialtools': ['*.py','*.so']}
      )
