import pathlib

_str_end = """
about:
  home: https://github.com/andreatramacere/jetset
  license: BSD-3
  summary: ''
  license_family: BSD

extra:
  recipe-maintainers:
    - andreatramacere
"""

_skip_list = ["pyqt", "swig"]

# Read requirements.txt
req_file = pathlib.Path("./requirements.txt")
req = [line.strip() for line in req_file.read_text().splitlines()]
req = [r for r in req if r and not any(skip in r for skip in _skip_list) and not r.startswith("#")]

np_str = ""
pkg_str_list = []
for r in req:
    if "numpy" in r:
        np_str = r
    pkg_str_list.append(r)

_str_start = f"""
{{% set data = load_setup_py_data(setup_file='../../../setup.py', from_recipe_dir=True) %}}
{{% set version = data.get('version')  %}}

package:
  name: jetset
  version: {{ version }}

source:
  path: ../../../

build:
  preserve_egg_dir: True
  script_env:
    - JETSETBESSELBUILD
  script: python -m pip install . --no-deps -vv

requirements:
  host:
    - python {{ '{{ python }}' }}
    - setuptools
    - {np_str}
"""

meta_path = pathlib.Path(".github/conda-pipeline/github/meta.yaml")
with meta_path.open("w") as f:
    print(_str_start, file=f)
    print("\n  run:", file=f)
    print("    - python", file=f)
    for pkg_str in pkg_str_list:
        print(f"    - {pkg_str}", file=f)
    print(_str_end, file=f)

print(f"Generated meta.yaml at {meta_path}")
