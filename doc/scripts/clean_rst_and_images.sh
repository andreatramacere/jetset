
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
    find documentation_notebooks/notebooks -name '*.png' -delete
    find user_guide/documentation_notebooks/notebooks -name '*.png' -delete
    find user_guide/documentation_notebooks/notebooks -name '*.rst' -delete


     	
else
    find documentation_notebooks/notebooks/$1 -name  '*.png' -delete
    find user_guide/documentation_notebooks/notebooks/$1 -name '*.png' -delete
    find user_guide/documentation_notebooks/notebooks/$1 -name '*.rst' -delete
fi	

