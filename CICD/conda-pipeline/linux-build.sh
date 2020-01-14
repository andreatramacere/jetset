cd /Users/orion/astro/Programmi/Projects/Active/JetSeT/JetSeT_CICD/MAC_OS/LINUX

export USE_PIP='FALSE'
export CONDABUILDSETUP='TRUE'
#building
cd integration
git clone https://github.com/andreatramacere/jetset.git
cd jetset
git checkout develop
git pull origin develop

#to build bessel fucntion locally
conda install --yes   -c astropy --file ../../requirements.txt
python setup.py clean
python setup.py insatll
python setup.py clean

cd ..
python -c 'import jetkernel; import os;p=os.path.join(jetkernel.__path__[0],'mathkernel','F_Sync.dat'); os.system('cp %s jetset/jetkernel/mathkernel'%p)'


cd jetset/CICD/conda-pipeline/

conda create --yes --name jetset-cidc python=3.7 ipython anaconda-client conda-build ipython
conda activate jetset-cidc
export PKG_VERSION=$(cd ../../ && python -c "import jetset;print(jetset.__version__)")
rm -rf ../../jetset/__pycache__/

python -c 'import os;os.environ["MPLBACKEND"]="Agg"; from jetset.tests.test_build_functions import test_build_bessel; test_build_bessel()'

#now using env var
#set the proper branch/tag in: mata.yaml-> git_rev:
 #for linux
conda build purge
conda build .  -c defaults -c astropy  #for linux
export CONDABUILDJETSET=$(conda-build . --output)



#testing
conda install  --yes --offline $CONDABUILDJETSET
python -c 'import os;os.environ["MPLBACKEND"]="Agg"; from jetset.tests import test_functions; test_functions.test_short()'
