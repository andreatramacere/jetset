cd /Users/orion/astro/Programmi/Projects/Active/JetSeT/JetSeT_CICD/MAC_OS/LINUX


#building
cd integration
git clone https://github.com/andreatramacere/jetset.git
cd jetset
git checkout develop
git pull origin develop

#to build bessel fucntion locally


export USE_PIP='FALSE'
export JETSETSKIPBESSELBUILD='FALSE'
conda install --yes   -c astropy --file requirements.txt

python setup.py clean
python setup.py install
python setup.py clean
cd ..
python -c 'import jetkernel; import os;p=os.path.join(jetkernel.__path__[0],"mathkernel"); os.system("cp jetset/jetkernel/mathkernel/F_Sync.dat %s"%p)'

export JETSETSKIPBESSELBUILD='TRUE'
cd jetset/CICD/conda-pipeline/

conda create --yes --name jetset-cidc python=3.7 ipython anaconda-client conda-build ipython
conda activate jetset-cidc
conda install --yes   -c astropy --file requirements.txt
export PKG_VERSION=$(cd ../../ && python -c "import jetset;print(jetset.__version__)")
rm -rf ../../jetset/__pycache__/


#now using env var
#set the proper branch/tag in: mata.yaml-> git_rev:
 #for linux
conda build purge
export CONDABUILDJETSET=$(conda-build . --output)
conda build .  -c defaults -c astropy  #for linux




#testing
conda install  --yes --offline $CONDABUILDJETSET
cd /workdir/test
python -c 'import os;os.environ["MPLBACKEND"]="Agg"; from jetset.tests import test_functions; test_functions.test_short()'
