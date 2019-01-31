
from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

__author__ = "Andrea Tramacere"

__all__=['commands']

def commands(obj):
    cmd_list= dir(obj)
    for item in cmd_list:
        if item[0]!='_':
            print (item)