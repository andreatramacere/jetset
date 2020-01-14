cd /Users/orion/astro/Programmi/Projects/Active/JetSeT/JetSeT_CICD/MAC_OS/CONDA

#building
cd integration
echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> git  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

git clone https://github.com/andreatramacere/jetset.git
cd jetset
git checkout develop
git reset --hard HEAD
git pull origin develop

echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> prepoc <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

export JETSETBESSELBUILD='FALSE'
export USE_PIP='FALSE'

python setup.py clean

cd CICD/conda-pipeline/macosx

echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> BUILD  jetset-cidc env <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

conda create --yes --name jetset-cidc python=3.7 ipython anaconda-client conda-build ipython
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate jetset-cidc
conda install --yes  -c astropy --file ../../../requirements.txt


echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> SET VERSION <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

export PKG_VERSION=$(cd ../../../ && python -c "import jetset;print(jetset.__version__)")
rm -rf ../../../jetset/__pycache__/
echo  $PKG_VERSION

#now using env var
#set the proper branch/tag in: mata.yaml-> git_rev:

echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> CONDA BUILD <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

conda build purge
conda-build .  -c defaults -c astropy
export CONDABUILDJETSET=$(conda-build . --output)
echo  $CONDABUILDJETSET

#testing
echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> TESTING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

conda install --yes --offline $CONDABUILDJETSET
cd ../../../../test/
python -c 'import os;os.environ["MPLBACKEND"]="Agg"; from jetset.tests import test_functions; test_functions.test_short()'
conda deactivate



