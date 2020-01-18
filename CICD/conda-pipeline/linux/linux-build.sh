#building


cd jetset


echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> prepoc <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

export USE_PIP='FALSE'
export JETSETBESSELBUILD='TRUE'



echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> BUILD BESSESL <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
conda install --yes   -c astropy --file requirements.txt
rm jetkernel/mathkernel/F_Sync.dat
python setup.py clean
python setup.py install > install.log 2>install.err
python setup.py clean

#cd ..
#python -c 'import jetkernel; import os;p=os.path.join(jetkernel.__path__[0],"mathkernel"); os.system("cp jetset/jetkernel/mathkernel/F_Sync.dat %s"%p)'

export JETSETBESSELBUILD='FALSE'
cd CICD/conda-pipeline/linux

echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> SET VERSION <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

export PKG_VERSION=$(cd ../../ && python -c "import jetset;print(jetset.__version__)")
rm -rf ../../../jetset/__pycache__/
echo  $PKG_VERSION


echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> CONDA BUILD  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<',$PKG_VERSION
conda clean -a
conda build purge
conda build .  -c defaults -c astropy > build.log 2>build.err #for linux
export CONDABUILDJETSET=$(conda-build . --output)
echo  $CONDABUILDJETSET



echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> BUILD  jetset-cidc env <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<',$CONDABUILDJETSET
conda create --yes --name jetset-cidc python=3.7 ipython anaconda-client conda-build>conda_env_build.log
source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate jetset-cidc



echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> TESTING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<',$CONDABUILDJETSET
#testing
conda install --yes   -c astropy --file ../../../requirements.txt
conda install  --yes --offline $CONDABUILDJETSET
cd /workdir/test
python -c 'import os;os.environ["MPLBACKEND"]="Agg"; from jetset.tests import test_functions; test_functions.test_short()'
conda deactivate
