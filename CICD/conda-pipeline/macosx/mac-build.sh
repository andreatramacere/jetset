

echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> prepoc <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
export JETSETBESSELBUILD='FALSE'
export USE_PIP='FALSE'

cd integration/jetset
python setup.py clean

cd CICD/conda-pipeline/macosx

echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> BUILD  jetset-cidc env <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<',$PWD
conda env remove --name jetset-cidc
conda create --yes --name jetset-cidc python=3.7 ipython anaconda-client conda-build ipython>conda_env_build.log
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate jetset-cidc
conda install --yes   -c astropy -c conda-forge --file ../../../requirements.txt


echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> SET VERSION <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

export PKG_VERSION=$(cd ../../../ && python -c "import jetset;print(jetset.__version__)")
rm -rf ../../../jetset/__pycache__/
echo  $PKG_VERSION

#now using env var
#set the proper branch/tag in: mata.yaml-> git_rev:

echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> CONDA BUILD  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<',$PKG_VERSION

conda build purge
conda-build .  > build.log 2>build.err
export CONDABUILDJETSET=$(conda-build . --output)
echo  $CONDABUILDJETSET

#testing
echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> TESTING  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<',$CONDABUILDJETSET

conda install --yes --offline $CONDABUILDJETSET
cd ../../../../../test
python -c 'import os;os.environ["MPLBACKEND"]="Agg"; from jetset.tests import test_functions; test_functions.test_short()'
echo 'export CONDABUILDJETSET='$CONDABUILDJETSET>../deploy/CONDABUILD.sh
conda deactivate



