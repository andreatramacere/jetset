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

def update_version(version):

    with open('jetset/pkg_info.json') as fp:
        _info = json.load(fp)
        _info[version]

    with open('jetset/pkg_info.json', 'w') as json_file:
        json.dumps(_info, json_file, default_flow_style=False)


    with open('jetset/pkg_info.json') as fp:
        _info = json.load(fp)


    with open("conda-pipeline/meta.yaml", 'r') as stream:
        data_loaded = yaml.load(stream)

        data_loaded['package']['version']=version

    with open('conda-pipeline/meta.yaml', 'w') as yaml_file:
        yaml.dump(data_loaded, yaml_file, default_flow_style=False)

def do_remote_tag(version):
    l=[]
    l.append(["git", "tag", "-d %s"%version])
    l.append(['git', 'push', '--delete origin %s'%version])
    l.append(['git', 'push', '--delete origin-github %s'%version])
    l.append(['git', 'tag' '-a %s'%version ,'''-m '%s' '''%version])
    l.append(['git', 'push', 'origin %s'%version])
    l.append(['git', 'push', 'origin-github %s' % version])



    for c in l:
        #print (l)
        res = subprocess.check_output(c)
        for line in res.splitlines():
            print (c)


def main(argv=None):
    parser = argparse.ArgumentParser()

    parser.add_argument('version', type=str, )
    parser.add_argument('-do_remote_tag', action='store_true')

    args = parser.parse_args()


    update_version(args.verions)
    if args.do_remote_tag==True:
        do_remote_tag(args.version)

if __name__ == "__main__":
    main(sys.argv)
