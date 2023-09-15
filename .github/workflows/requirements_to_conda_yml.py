
_str_end="""
about:
  home: https://github.com/andreatramacere/jetset
  license: BSD-3
  summary: ''
  license_family: BSD

extra:
  recipe-maintainers:
    - andreatramacere
"""




_skip_list=['pyqt','swig']
f = open("./requirements.txt",'r')
req=f.readlines()
f.close()

req=[n.strip() for n in req]
for r in req[:]:
    for s in _skip_list[:]:

        if s in r:
            req.remove(r)


np_str=''
pkg_str_list=[]
for r in req:
    if r.startswith('#') is False:
      if 'numpy' in r:
          np_str=r
      pkg_str_list.append(r)

_str_start="""
{% set data = load_setup_py_data(setup_file='../../../setup.py', from_recipe_dir=True) %}
{% set version = data.get('version')  %}

package:
  name: jetset
  version:  {{ version }}

source:
  path: ../../../

build:
  preserve_egg_dir: True
  script_env:
    - JETSETBESSELBUILD

requirements:

  build:
    - swig>3.0.0
    - python {{ python }}
    - setuptools"""

f = open(".github/conda-pipeline/github/meta.yaml",'w')
print(_str_start,file=f)
print( '    - %s'%np_str, file=f)

print('',file=f)
print('  run:',file=f)
print('    - python>=3.8', file=f)
for pkg_str in pkg_str_list:
    print('    - %s'%pkg_str, file=f)
print(_str_end,file=f)
f.close()