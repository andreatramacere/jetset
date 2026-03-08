# Instruction to build sphinx doc



Notebooks are built in the `documentation_notebooks/user_guide`, update them before building doc!

To avoid to import the jetkernel module otherwise the build will crash: `export READTHEDOCS='True'`


## Prerequisites 
1) build the requirements file: `python build_rtd_requirements.py`

2) and install the requirements

3) to build proper api strucutre: `python make_apidoc_and_uml_graphs.py`
 
   
## Building doc in separate steps: 

### Cleaning

to clean everything under documentation_notebooks/notebooks: `./scripts/clean_rst_and_images.sh` 

to clean everything under documentation_notebooks/notebooks/dir-name: `./scripts/clean_rst_and_images.sh dir-name`


### Building notebooks rst
3) run step a or b
 
    a) run all notebooks and build rst file with script:
      - run all notebooks found in documentation_notebooks/notebooks  :`./scripts/build_rst_files.sh -e`

        OR     

      - dir-name, runs only documentation_notebooks/notebooks/dir-name: `./scripts/build_rst_files.sh -e`  



    b) only build rst file with script
      - `./scripts/build_rst_files.sh` 

        OR
      
      - `./scripts/build_rst_files.sh dir-name`



### Update the rst in user_guide from the notebooks in documentation_notebooks
  
  - `./scripts/clean_rst_and_images.sh`
  - `./scripts/update_rts_images.sh`

### builds the sphinx docs
 - `sphinx-build -b html ./ build`
 
 OR

 - `sphinx-build -j 10 -b html ./ build`   to use 10 parallel jobs


## Building doc in separate steps: 

argument to add:
 - `-c` will clean rst/png prods
 - `-e` the notebook will be executed
 - `-d` will build the doc

to build without executing:
- `./scripts/build.sh  -cb`

OR

- ./scripts/build.sh  -cb dir-name

to build executing:

- `./scripts/build.sh -cbe` 

OR

- `./scripts/build.sh -cbe dir-name`
