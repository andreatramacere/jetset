#!/bin/bash
source $CONDA_PREFIX/etc/profile.d/conda.sh


echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> prepoc <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
export JETSETBESSELBUILD='FALSE'
export USE_PIP='FALSE'

cd integration/jetset

conda activate root
echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> BUILD  jetset-cidc env <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<',$PWD
conda create --yes --name jetset-cidc python=3.7 ipython anaconda-client conda-build ipython>conda_env_build.log
conda activate jetset-cidc
conda install --yes   -c conda-forge emcee">=3.0.0"
conda install --yes   -c astropy --file requirements.txt

python setup.py clean
python setup.py install > install.log 2>install.err

cd CICD/conda-pipeline/macosx

echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> SET VERSION <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

export PKG_VERSION=$(cd ../../../ && python -c "import jetset;print(jetset.__version__)")
rm -rf ../../../jetset/__pycache__/
echo  $PKG_VERSION

#now using env var
#set the proper branch/tag in: mata.yaml-> git_rev:

echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> CONDA BUILD  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<',$PKG_VERSION

echo '---> purge'
conda build purge
echo '---> build'
conda build .  -c defaults -c astropy -c conda-forge > build.log 2>build.err
export CONDABUILDJETSET=$(conda-build . --output)
echo  $CONDABUILDJETSET
echo '---> done'

#testing
echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> TESTING  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<',$CONDABUILDJETSET

conda install --yes --offline $CONDABUILDJETSET
cd ../../../../../test
python -c 'import os;os.environ["MPLBACKEND"]="Agg"; from jetset.tests import test_functions; test_functions.test_short()'
echo 'export CONDABUILDJETSET='$CONDABUILDJETSET>../deploy/CONDABUILD.sh
conda deactivate



