#!/usr/bin/env python


from __future__ import division, absolute_import, print_function

__author__ = 'andrea tramacere'




from setuptools import setup, find_packages,Extension
from setuptools.command.install import install
#from distutils.extension import Extension
import distutils.command.install as orig
from distutils.command.build import build
from setuptools.command.install import install
from distutils.sysconfig import get_python_lib
import os
import glob
import shutil
import fnmatch
import json
import sys


class CustomBuild(build):
    def run(self):
        print('----> custom build')

        self.run_command('build_ext')
        build.run(self)


class CustomInstall(install):
    def run(self):
        print('----> custom install',self.get_command_name())


        self.run_command('build_ext')

        #if 'pip' in __file__:
        #    install.run(self)
        #else:
        #    self.do_egg_install()

        try:
            self.do_egg_install()
        except:
            install.run(self)




class CustomClean(install):
    def run(self):



        try:
            shutil.rmtree('dist')
        except:
            pass
        try:
            shutil.rmtree('build')
        except:
            pass
        try:
            shutil.rmtree(glob.glob('*.egg-info')[0])
        except:
            pass
        try:
            os.remove('jetkernel/jetkernel.py')
        except:
            pass
        try:
            os.remove('jetkernel/jetkernel_wrap.c')
        except:
            pass
        try:
            shutil.rmtree('jetkernel/__pycache__')
        except:
            pass

        #to remove files installed by old versions
        site_p=get_python_lib()
        #print('path',site_p)
        #print(site_p, glob.glob(site_p + '/*_jetkernel*'))
        for f in glob.glob(site_p+'/*_jetkernel*'):
            print ('found .so object:', f)
            print ('removing it')
            try:
                os.remove(f)
            except:
                pass



custom_cmdclass = {'build': CustomBuild,
                   'install': CustomInstall,
                   'clean':CustomClean}





with open('jetset/pkg_info.json') as fp:
    _info = json.load(fp)

__version__ = _info['version']


f = open("./requirements.txt",'r')
req=f.readlines()
f.close()
req=[n.strip() for n in req]

if  os.getenv('USE_PIP')=='TRUE':
    install_req=req
else:
    install_req=None

src_files=['jetkernel/jetkernel.i']
src_files.extend(glob.glob ('jetkernel_src/src/*.c'))
_module=Extension('jetkernel/_jetkernel',
                  sources=src_files,
                  #extra_compile_options='-fPIC  -v  -c -m64 -I',
                  #extra_link_options='-suppress',
                  swig_opts=['-v',],
                  include_dirs=['jetkernel_src/include'])

#'jetkernel/mathkernel/*dat'

with open("proj_descr.md", "r") as f:
    long_description = f.read()

print(__version__)

setup(name='jetset',
      version=__version__,
      author='Andrea Tramacere',
      url='https://github.com/andreatramacere/jetset',
      long_description=long_description,
      long_description_content_type='text/markdown',
      description="A framework for self-consistent modeling and fitting of  astrophysical relativistic jets SEDs",
      author_email='andrea.tramacere@gmail.com',
      packages=['jetset', 'leastsqbound', 'jetkernel','jetset.tests'],
      package_data={'jetset':['Spectral_Templates_Repo/*.dat','test_data/SEDs_data/*ecsv','./requirements.txt'],'jetkernel':['mathkernel/*dat']},
      include_package_data = True,
      cmdclass=custom_cmdclass,
      #not to use with setup tools
      #requires=req,
      ext_modules = [_module],
      install_requires=install_req,
      py_modules=['jetkernel/jetkernel'],
      python_requires='>=3.5',
      zip_safe=False)
