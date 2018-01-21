#!/usr/bin/env python


from __future__ import division, absolute_import, print_function

__author__ = 'andrea tramacere'




from setuptools import setup, find_packages,Extension
from distutils.extension import Extension
import glob

version='1.2.0'

src_files=['jetset/jetkernel/jetkernel.i']
src_files.extend(glob.glob ('jetkernel_src/src/*.c'))
_module=Extension('_jetkernel',
                  sources=src_files,
                  extra_compile_options='-fPIC  -v  -c -m64 -I',
                  extra_link_options='-suppress',
                  swig_opts=['-v',],
                  include_dirs=['jetkernel_src/include'])

setup(name='jetset',
      version=version,
      author='Andrea Tramacere',
      author_email='andrea.tramacere@gmail.com',
      packages=['jetset','leastsqbound','jetset.jetkernel'],
      package_data={'jetset':['Spectral_Templates_Repo/*.dat','test_data/SEDs_data/*dat','jetkernel/mathkernel/*dat']},
      scripts=['bin/test_interactive.py'],
      requires=['scipy','numpy','astropy'],
      ext_modules = [_module],
      py_modules=['jetkernel'], )
