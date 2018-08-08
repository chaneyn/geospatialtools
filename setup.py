import os
from setuptools import setup, find_packages

#First compile the code
cwd = os.getcwd()
os.chdir('src')
os.system('python compile.py')
os.chdir(cwd)

#base_dir = os.path.dirname(os.path.abspath(__file__))
#requirements_file = open(os.path.join(base_dir, 'requirements.txt'))
#requirements = requirements_file.read().splitlines()
setup(name='geospatialtools',
      version='1.0',
      package_dir={'geospatialtools': 'libraries'},
      packages=['geospatialtools'],
      package_data={'geospatialtools': ['*.py','*.so']}
      )

'''setup(
    name='geospatialtools',
    version='0.1',
    #cmdclass=versioneer.get_cmdclass(),
    author='Nathaniel Chaney',
    author_email='chaneyna@gmail.com',
    packages=find_packages(),
    install_requires=requirements,
    ##package_dir={'geospatialtools': 'libraries'},
    #packages=['geospatialtools'],
    #package_data={'geospatialtools': ['*.py','*.so']}
    #entry_points={
    #    'console_scripts': [
    #        'anaconda = binstar_client.scripts.cli:main',
    #        'binstar = binstar_client.scripts.cli:main',
    #        'conda-server = binstar_client.scripts.cli:main'
    #    ]
    #},
    #license='BSD License',
    classifiers=[
        'Development Status :: 5 - Production/Stable',
        'Environment :: Console',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python',
    ]
)'''
