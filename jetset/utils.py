
from __future__ import absolute_import, division, print_function

from builtins import (bytes, str, open, super, range,
                      zip, round, input, int, pow, object, map, zip)

__author__ = "Andrea Tramacere"

__all__=['check_frame','unexpetced_behaviour']
import re




def check_frame(frame):
    allowed=['obs','src','blob']
    if frame not in allowed:
        raise RuntimeError('rest frame', frame, 'not allowed',allowed)


def unexpetced_behaviour():
    raise RuntimeError('the code reached a condition that should never happen!')




def clean_var_name(s):

   s.replace('-','_')
   s.replace(' ', '_')

   # Remove invalid characters
   s = re.sub('[^0-9a-zA-Z_]', '', s)

   # Remove leading characters until we find a letter or underscore
   s = re.sub('^[^a-zA-Z_]+', '', s)

   return s