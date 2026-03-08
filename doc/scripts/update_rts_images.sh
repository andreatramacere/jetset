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


mkdir -p ./user_guide/documentation_notebooks_rst


if [ "$#" -eq  "0" ]
then
	echo "No arguments supplied"
    (cd documentation_notebooks && find notebooks -name '*.rst' | cpio -pdm ../user_guide/documentation_notebooks_rst)
    (cd documentation_notebooks && find notebooks -name '*.png' | cpio -pdm ../user_guide/documentation_notebooks_rst)

    (cd documentation_notebooks && find slides -name '*.png' | cpio -pdm ../user_guide/documentation_notebooks_rst)
    (cd documentation_notebooks && find images -name '*.png' | cpio -pdm ../user_guide/documentation_notebooks_rst)



     	
else
    (cd documentation_notebooks && find notebooks/"$1" -name '*.rst' | cpio -pdm ../user_guide/documentation_notebooks_rst)
    (cd documentation_notebooks && find notebooks/"$1" -name '*.png' | cpio -pdm ../user_guide/documentation_notebooks_rst)

    (cd documentation_notebooks && find slides -name '*.png' | cpio -pdm ../user_guide/documentation_notebooks_rst)
    (cd documentation_notebooks && find images -name '*.png' | cpio -pdm ../user_guide/documentation_notebooks_rst)
fi	



