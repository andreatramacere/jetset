#!/usr/bin/env python

from __future__ import absolute_import, division, print_function

from builtins import (bytes, open, super, range,
                      zip, round, input, pow, object, map, zip)

try:
    # Python 2
    from __builtin__ import str as builtin_str
except ImportError:
    # Python 3
    from builtins import str as builtin_str

import yaml,json
import sys
import argparse
import subprocess
from collections import OrderedDict

def ordered_dict_representer(self, value):  # can be a lambda if that's what you prefer
    return self.represent_mapping('tag:yaml.org,2002:map', value.items())


def get_tag_from_version():

    with open('jetset/pkg_info.json') as fp:
        _info = json.load(fp)
      
        version=_info['version']
        
    print(version)
    return version

def do_tag(tag,remote=False):
    l=[]
    l.append(["git tag  -d %s"%tag])
    l.append(['git push --delete origin %s'%tag])
    #l.append(['git push --delete origin-github %s'%tag])
    l.append(['''git tag -a %s  -m '%s' '''%(tag,tag)])
    l.append(['git push origin %s'%tag])
    #l.append(['git push origin-github %s' % tag])

    res = subprocess.call(l[0], shell=True)
    for  i in l[1::]:
        try:
            res = subprocess.call(i, shell=True)
        except:
            pass
    print()

def main(argv=None):
    parser = argparse.ArgumentParser()

    parser.add_argument('-tag', type=str, )
    parser.add_argument('-do_remote_tag', action='store_true')
    parser.add_argument('-dry', action='store_true')

    args = parser.parse_args()

    tag=args.tag

    tag=get_tag_from_version()
    print("tag is ",tag)
    if args.dry is False:
        do_tag(tag,remote=args.do_remote_tag)
    
       


if __name__ == "__main__":
    main(sys.argv)
