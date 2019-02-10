"""
JetSeT package
"""


import pkgutil
import os
import json

__author__ = "Andrea Tramacere"


pkg_dir = os.path.abspath(os.path.dirname(__file__))
pkg_name = os.path.basename(pkg_dir)
__all__=[]

with open(os.path.dirname(__file__) + '/pkg_info.json') as fp:
    _info = json.load(fp)

__version__ = _info['version']

if 'label'  in _info.keys():
    __label__= _info['label']
else:
    __label__= None

for importer, modname, ispkg in pkgutil.walk_packages(path=[pkg_dir],
                                                      prefix=pkg_name+'.',
                                                      onerror=lambda x: None):

    if ispkg == True:
        __all__.append(modname)
    else:
        pass


    data_dir=os.path.dirname(__file__)+'/data'