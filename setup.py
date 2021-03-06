#!/usr/bin/env python
"""
Setup script for Python DDD library.
Started by M. Stricker, 2015

install as
python setup.py install --u
"""

from setuptools import setup, find_packages
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))



setup(name='DDDutils',
      version='1.0',
      
      description='DDD postprocessing routines in python',

      # URL
      url = "git@github.com:mastricker/DDDutils.git",
      
      # Author
      author='Markus Stricker',
      author_email='mail@markusstricker.de',

      license='MIT',

      classifiers=[
          # How mature is this project? Common values are
          # 3 - Alpha
          # 4 - Beta
          # 5 - Production/Stable
          'Development Status :: 3 - Alpha'

          'License :: OSI Approved :: MIT License'

      ],

      keywords='DDD, KIT, postprocessing',

      packages = find_packages(),
     )


