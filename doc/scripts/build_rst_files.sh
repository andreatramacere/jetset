#!/bin/sh

execute=0
while getopts "h?e" opt; do
    case "$opt" in
    h|\?)
        echo "Usage: ./scripts/build_rst_files.sh [-e] [dir-name]"
        exit 0
        ;;
    e)
        execute=1
        ;;
    esac
done

shift $((OPTIND-1))

if [ "$#" -eq "0" ]; then
    echo "No arguments supplied"
    search_dir="documentation_notebooks/notebooks"
else
    search_dir="documentation_notebooks/notebooks/$1"
fi

if [ ! -d "$search_dir" ]; then
    echo "Directory not found: $search_dir" >&2
    exit 1
fi

if ! command -v jupyter >/dev/null 2>&1; then
    echo "jupyter command not found" >&2
    exit 1
fi

echo "----> $execute"
find "$search_dir" -name '*.ipynb' -not -path '*/\.*' | while IFS= read -r file; do
    echo "$file"
    if [ "$execute" -eq 1 ]; then
        echo 'execute'
        jupyter nbconvert --execute "$file" --to rst
    else
        echo 'non execute'
        jupyter nbconvert "$file" --to rst
    fi
done
