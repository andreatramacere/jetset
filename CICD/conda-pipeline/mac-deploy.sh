#deploy
anaconda login
anaconda upload /path/produced/by/conda-build

#tests after deploy
conda uninstall jetset
conda install -c andreatramacere jetset
python -c 'import os;os.environ["MPLBACKEND"]="Agg"; from jetset.tests import test_functions; test_functions.test_short()'

conda deactivate
conda env remove --name jetset-cidc
