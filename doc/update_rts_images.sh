execute=0
while getopts "h?e" opt; do
    case "$opt" in
    h|\?)
        echo ''
        exit 0
        ;;
    esac
done

shift $((OPTIND-1))



if [ "$#" -eq  "0" ]
then
	echo "No arguments supplied"
    find documentation_notebooks/notebooks -name '*.rst' # | cpio -pdm ./user_guide
    find documentation_notebooks/notebooks -name '*.png' #| cpio -pdm ./user_guide

    find documentation_notebooks/slides -name '*.png' #| cpio -pdm ./user_guide
    find documentation_notebooks/images -name '*.png' #| cpio -pdm ./user_guide



     	
else
    find documentation_notebooks/notebooks/$1 -name '*.rst' | cpio -pdm ./user_guide
    find documentation_notebooks/notebooks/$1 -name '*.png' | cpio -pdm ./user_guide

    find documentation_notebooks/slides/ -name '*.png' | cpio -pdm ./user_guide
    find documentation_notebooks/images/ -name '*.png' | cpio -pdm ./user_guide
fi	





