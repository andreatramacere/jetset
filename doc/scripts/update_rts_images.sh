#!/bin/sh

while getopts "h?" opt; do
    case "$opt" in
    h|\?)
        echo "Usage: ./scripts/update_rts_images.sh [dir-name]"
        exit 0
        ;;
    esac
done

shift $((OPTIND-1))

if ! command -v cpio >/dev/null 2>&1; then
    echo "cpio command not found" >&2
    exit 1
fi

src_root="documentation_notebooks"
dest_root="../user_guide/documentation_notebooks_rst"

mkdir -p ./user_guide/documentation_notebooks_rst

copy_by_pattern() {
    rel_dir="$1"
    pattern="$2"
    if [ -d "$src_root/$rel_dir" ]; then
        (
            cd "$src_root" &&
            find "$rel_dir" -type f -name "$pattern" | cpio -pdm "$dest_root"
        )
    fi
}

if [ "$#" -eq "0" ]; then
    echo "No arguments supplied"
    copy_by_pattern "notebooks" "*.rst"
    copy_by_pattern "notebooks" "*.png"
else
    scoped_dir="notebooks/$1"
    if [ ! -d "$src_root/$scoped_dir" ]; then
        echo "Directory not found: $src_root/$scoped_dir" >&2
        exit 1
    fi
    copy_by_pattern "$scoped_dir" "*.rst"
    copy_by_pattern "$scoped_dir" "*.png"
fi

copy_by_pattern "slides" "*.png"
copy_by_pattern "images" "*.png"
