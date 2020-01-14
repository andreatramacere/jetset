#deploy
cd /workdir


source $CONDA_PREFIX/etc/profile.d/conda.sh
conda activate jetset-cidc

anaconda login
anaconda upload --force $CONDABUILDJETSET

#tests after deploy
conda uninstall jetset
conda install --yes -c andreatramacere -c astropy jetset #for linux
cd test
python -c 'import os;os.environ["MPLBACKEND"]="Agg"; from jetset.tests import test_functions; test_functions.test_short()'

conda deactivate
conda env remove --name jetset-cidc
