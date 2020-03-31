
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
	echo "No arguments supplied"
     	files=($(find . -name   '*.ipynb' -not -path '*/\.*'))
else
	files=($(find $1 -name  '*.ipynb' -not -path '*/\.*'))
fi	

echo '---->' $execute
for file in "${files[@]}"; do
  echo $file
  if [ $execute -eq 1 ]
  then 
	  echo 'execute'
	  jupyter nbconvert --execute $file --to rst
  else
	  echo 'non execute'
      jupyter nbconvert $file --to rst 
 
  fi
done
