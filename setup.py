#!/usr/bin/env python
from setuptools import setup
from setuptools import find_packages
import numpy as np

setup(name='sph_models',
      version='0.1',
      description='Package for handling global tomography models in sph (spherical harmonic)\
      parameterization of Ritsema et al. (2011)',
      author='Ross Maguire',
      author_email='rmaguire@umd.edu',
      url='www.github.com/romaguir/sph_tools',
      packages=find_packages(),
      install_requires=['pyshtools','fortranformat','cartopy'],
      license='GNU')
