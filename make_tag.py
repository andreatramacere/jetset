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


def update_version(version):

    with open('jetset/pkg_info.json') as fp:
        _info = json.load(fp)
        #_info['version']=version
        label=subprocess.check_output(['git', 'rev-parse', 'HEAD']).decode('ascii').strip()
        _info['label']= '%s'%label
    print(_info)
    with open('jetset/pkg_info.json', 'w') as json_file:
        json.dump(_info, json_file)




    #with open("conda-pipeline/_template_meta.yaml", 'r') as stream:
    #    data_loaded = yaml.load(stream)

    #    data_loaded['package']['version']=version

    #yaml.add_representer(OrderedDict, ordered_dict_representer)
    #with open('conda-pipeline/meta.yaml', 'w') as yaml_file:
    #    yaml.dump(data_loaded, yaml_file, default_flow_style=False)

def do_remote_tag(version):
    l=[]
    l.append(["git tag  -d %s"%version])
    l.append(['git push --delete origin %s'%version])
    l.append(['git push --delete origin-github %s'%version])
    l.append(['''git tag -a %s  -m '%s' '''%(version,version)])
    l.append(['git push origin %s'%version])
    l.append(['git push origin-github %s' % version])

    res = subprocess.call(l[0], shell=True)
    for  i in l[1::]:
        try:
            res = subprocess.call(i, shell=True)
        except:
            pass

def do_remote_push(message):
    l=[]
    l.append(["git add . "])
    l.append(['''git commit -m '%s' ''' % (message)])
    l.append(['git push origin py23'])
    l.append(['git push origin-github py23 %s'])

    res = subprocess.call(l[0], shell=True)
    for  i in l[1::]:
        try:
            res = subprocess.call(i, shell=True)
        except:
            pass


def main(argv=None):
    parser = argparse.ArgumentParser()

    parser.add_argument('version', type=str, )
    parser.add_argument('-do_remote_tag', action='store_true')
    parser.add_argument('-push_message', type=str, default=None)

    args = parser.parse_args()


    update_version(args.version)
    if args.push_message is not None:
        do_remote_push(args.push_message)
    if args.do_remote_tag==True:
        do_remote_tag(args.version)

if __name__ == "__main__":
    main(sys.argv)
