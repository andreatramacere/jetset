#!/bin/sh

while getopts "h?" opt; do
    case "$opt" in
    h|\?)
        echo "Usage: ./scripts/clean_rst_and_images.sh [dir-name]"
        exit 0
        ;;
    esac
done

shift $((OPTIND-1))

delete_generated() {
    dir="$1"
    if [ -d "$dir" ]; then
        find "$dir" -type f \( -name '*.png' -o -name '*.rst' \) -delete
    fi
}

if [ "$#" -eq "0" ]; then
    echo "No arguments supplied, deleting .rst and .png files in all notebooks dirs"
    delete_generated "documentation_notebooks/notebooks"
    delete_generated "user_guide/documentation_notebooks_rst/notebooks"
else
    notebook_dir="documentation_notebooks/notebooks/$1"
    if [ ! -d "$notebook_dir" ]; then
        echo "Directory not found: $notebook_dir" >&2
        exit 1
    fi
    delete_generated "$notebook_dir"
    delete_generated "user_guide/documentation_notebooks_rst/notebooks/$1"
fi
