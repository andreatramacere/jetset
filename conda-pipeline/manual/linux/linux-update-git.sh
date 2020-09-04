#path on miniconda docker
cd /workdir


#building
echo  '>>>>>>>>>>>>>>>>>>>>>>>>>>> git  <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

cd integration
git clone https://github.com/andreatramacere/jetset.git
cd jetset
git checkout develop
git reset --hard HEAD
git pull origin develop
cd /workdir
