
execute=0
while getopts "h?e" opt; do
    case "$opt" in
    h|\?)
        echo ''
        exit 0
        ;;
    e)  execute=1
        ;;
    esac
done

shift $((OPTIND-1))

if [ "$#" -eq  "0" ]
then
	echo "No arguments supplied, please provide the directory with documentation_notebooks to convert"
fi

echo '--------------------------------'
echo 'generating rst files'
if [ $execute -eq 1 ]
then 
    echo 'execute'
    ./build_rst_files.sh -e $1
else
    echo 'non execute'
    ./build_rst_files.sh $1

fi
echo '--------------------------------'
echo
echo

echo '--------------------------------'
echo 'copying rst files'
find ./documentation_notebooks -name '*.rst' | cpio -pdm ./user_guide
find ./documentation_notebooks -name '*.png' | cpio -pdm ./user_guide
echo '--------------------------------'
echo
echo

echo '--------------------------------'
echo 'running sphinx-build'
sphinx-build -b html ./ build
echo '--------------------------------'