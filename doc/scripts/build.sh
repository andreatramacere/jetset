execute=0
clean=0
build=0
while getopts "h?ecb" opt; do
    case "$opt" in
    h|\?)
        echo ''
        exit 0
        ;;
    e)  execute=1
        ;;
    c)  clean=1
        ;;
    b)  build=1
    esac
done

shift $((OPTIND-1))

if [ "$#" -eq  "0" ]
then
	echo "No arguments supplied, please provide the directory with documentation_notebooks to convert"
fi


echo '--------------------------------'
echo 'cleam rst/png files'
if [ $clean -eq 1 ]
then 
   ./scripts/clean_rst_and_images.sh $1
fi
echo '--------------------------------'
echo
echo

echo '--------------------------------'
echo 'generating rst/png files'
if [ $execute -eq 1 ]
then 
    echo 'execute'
    ./scripts/build_rst_files.sh -e $1
else
    echo 'non execute'
    ./scripts/build_rst_files.sh $1

fi
echo '--------------------------------'
echo
echo

echo '--------------------------------'
echo 'copying rs/png files and images'
./scripts/update_rts_images.sh $1

echo '--------------------------------'
echo
echo



if [ $build -eq 1 ]
then 
    echo '--------------------------------'
    echo 'running sphinx-build'
    sphinx-build -j 10 -b html ./ build
    echo '--------------------------------'
else
    echo 'not running sphinx-build '
fi

