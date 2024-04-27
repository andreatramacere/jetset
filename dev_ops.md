## local operations

-  update version: update the tag in  `jetset/pkg_info.json`

-  clean `__pycache__`: find . -type d -name '__pycache__' | xargs rm -r 

-  If you want to make a tag you can use the script, the  `-do_remote_tag` option will create the tag remotely,
   deleting remotely the one with the same name, e.g.: 
    - `./make_tag.py 1.1.2 -do_remote_tag`
   
   To use as tag the current version
    - `./make_tag.py  -do_remote_tag`
   
   To just print the tag without creating it:
   - `./make_tag.py  -dry`
   
## operations on the action workflow
- If you want  to create a release from the same branch of the workflow
  1) set the release tag in the 'tag to create a release' field

- If you want to create a release from a specific tag or branch (different from the branch of the workflow)
  1) set the git tag in the `checkout this tag` field of the action wf
  2) set the release tag in the `tag to create a release` field of action wf

 - when publishing a new version create both the tag for the version and the stable

- Anaconda `meta.yaml` is built from `requirements.txt` using: `python .github/workflows/requirements_to_conda_yml.py` directly in the action wf  

- if you do not pass a `tag to create a release` value, the Release will not be created
- to run test on action leave the `skip test` form empty, otherwise type `yes` to skip the tests
## testing locally
<!-- - python -c"import iminuit; print('iminuit',iminuit.__version__); import jetset; print('jetset',jetset.#__version__)" -->

- user test:
  - `pytest  --pyargs  -vvv jetset.tests.test_users::TestUser`

- specific module test (eg: test_jet_model):
  - `pytest  --pyargs -vvv -s jetset.tests.test_jet_model`

- integration test :
  - `pytest  --pyargs -vvv jetset.tests.test_integration::TestIntegration`

- test a specificmethod:
 - `pytest  --pyargs -vvv jetset.tests.test_integration::TestIntegration::tets_method`

## Installation 
- installing from source without setup.py install
   - `rm -rf ./build`
   - `pip install --verbose .`