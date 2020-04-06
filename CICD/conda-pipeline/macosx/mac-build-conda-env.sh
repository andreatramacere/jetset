echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> BUILD  jetset-cidc env <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<',$PWD
conda create --yes --name jetset-cidc python=3.7 ipython anaconda-client conda-build ipython>conda_env_build.log
