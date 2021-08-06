#!/usr/bin/env python


from __future__ import division, absolute_import, print_function

__author__ = 'andrea tramacere'


from setuptools import setup, find_packages,Extension
# from setuptools.command.install import install
# from distutils.extension import Extension
# import distutils.command.install as orig
from distutils.command.build import build
from setuptools.command.install import install
from distutils.sysconfig import get_python_lib
import os
import glob
import shutil
import fnmatch
import json
import sys


def check_swig():
    command = 'swig'
    if shutil.which(command) is None:
        _mess = """ 
         ***********************************************************************************************************
         ***  you need to install swig v>=3.0.0 to install from source                                           ***
         ***                                                                                                     ***
         ***  - on linux Ubuntu: sudo apt-get install swig                                                       ***
         ***                                                                                                     ***
         ***  - on linux Debian: sudo aptitude install swig                                                      ***
         ***                                                                                                     ***
         ***  - on linux Fedora: sudo yum install swig                                                           ***
         ***                                                                                                     ***
         ***  - on mac:                                                                                          ***
         ***   /bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"   ***
         ***   brew install swig                                                                                 ***
         ***                                                                                                     ***
         ***   visit: http://www.swig.org/  for more info                                                        ***
         ***********************************************************************************************************
        """

        raise RuntimeError(_mess)

class CustomBuild(build):
    def run(self):
        print('-> custom build')
        check_swig()
        self.run_command('build_ext')
        build.run(self)


class CustomInstall(install):
    def run(self):
        print('-> custom install',self.get_command_name())
        check_swig()
        self.run_command('build_ext')
        install.run(self)
        print ('JETSETBESSELBUILD',os.getenv('JETSETBESSELBUILD') == 'TRUE')
        if os.getenv('JETSETBESSELBUILD') == 'TRUE':
            self.run_command('test')
        else:
            pass




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
            os.remove('jetset/jetkernel/jetkernel.py')
        except:
            pass
        try:
            os.remove('jetset/jetkernel/jetkernel_wrap.c')
        except:
            pass
        try:
            shutil.rmtree('jetset/jetkernel/__pycache__')
        except:
            pass

        # to remove files installed by old versions
        site_p=get_python_lib()

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
req=[n.strip() for n in req  if n.startswith('#') is False]


src_files=['jetset/jetkernel/jetkernel.i']
src_files.extend(glob.glob ('jetkernel_src/src/*.c'))
_module=Extension('jetset.jetkernel/_jetkernel',
                  sources=src_files,
                  #extra_compile_options='-fPIC  -v  -c -m64 -I',
                  #extra_link_options='-suppress',
                  swig_opts=['-v','-threads'],
                  include_dirs=['jetkernel_src/include'])


if os.getenv('JETSETBESSELBUILD') == 'TRUE':
    _test_suite = 'jetset.tests.test_build_functions'
else:
    _test_suite = None

with open("proj_descr.md", "r") as f:
    long_description = f.read()


# this to skip that pip install packages when installing src from conda
is_conda = os.path.exists(os.path.join(sys.prefix, 'conda-meta'))

if is_conda:
    install_req=None
else:
    install_req=req

print('-> version', __version__, install_req)

setup(name='jetset',
      version=__version__,
      author='Andrea Tramacere',
      url='https://github.com/andreatramacere/jetset',
      long_description=long_description,
      long_description_content_type='text/markdown',
      description="A framework for self-consistent modeling and fitting of  astrophysical relativistic jets SEDs",
      author_email='andrea.tramacere@gmail.com',
      packages=['jetset', 'jetset.leastsqbound', 'jetset.jetkernel','jetset.tests'],
      package_data={'jetset':['Spectral_Templates_Repo/*.dat','test_data/SEDs_data/*ecsv','./requirements.txt','ebl_data/*','mathkernel/*dat'],'jetkernel':['mathkernel/*dat']},
      include_package_data = True,
      cmdclass=custom_cmdclass,
      ext_modules = [_module],
      install_requires=install_req,
      py_modules=['jetset.jetkernel/jetkernel'],
      python_requires='>=3.7',
      test_suite =_test_suite,
      zip_safe=True)
