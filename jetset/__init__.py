"""
JetSeT package
"""


import pkgutil
import os
#import json
from .utils import get_info

from . import jetkernel

if jetkernel._refactored_ is True:
    pass
else:
    raise RuntimeError('please, uninstall the current version of jetekernel')

__author__ = "Andrea Tramacere"


pkg_dir = os.path.abspath(os.path.dirname(__file__))
pkg_name = os.path.basename(pkg_dir)
__all__ = []

_info = get_info()
__version__ = _info['version']

if 'label'  in _info.keys():
    __label__= _info['label']
else:
    __label__= None

for importer, modname, ispkg in pkgutil.walk_packages(path=[pkg_dir],
                                                      prefix=pkg_name+'.',
                                                      onerror=lambda x: None):
    if ispkg is True:
        __all__.append(modname)
    else:
        pass


data_dir = os.path.dirname(__file__)+'/data'

