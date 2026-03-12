#!/bin/sh

execute=0
clean=0
build=0

while getopts "h?ecb" opt; do
    case "$opt" in
    h|\?)
        echo "Usage: ./scripts/build.sh [-c] [-e] [-b] [dir-name]"
        exit 0
        ;;
    e)
        execute=1
        ;;
    c)
        clean=1
        ;;
    b)
        build=1
        ;;
    esac
done

shift $((OPTIND-1))

target_dir=""
if [ "$#" -ge 1 ]; then
    target_dir="$1"
fi

if [ -n "$target_dir" ] && [ ! -d "documentation_notebooks/notebooks/$target_dir" ]; then
    echo "Directory not found: documentation_notebooks/notebooks/$target_dir" >&2
    exit 1
fi

if [ -z "$target_dir" ]; then
    echo 'No directory supplied, operating on all notebooks'
fi

echo '--------------------------------'
echo 'clean rst/png files'
if [ "$clean" -eq 1 ]; then
    if [ -n "$target_dir" ]; then
        ./scripts/clean_rst_and_images.sh "$target_dir"
    else
        ./scripts/clean_rst_and_images.sh
    fi
fi
echo '--------------------------------'
echo
echo

echo '--------------------------------'
echo 'generating rst/png files'
if [ "$execute" -eq 1 ]; then
    echo 'execute'
    if [ -n "$target_dir" ]; then
        ./scripts/build_rst_files.sh -e "$target_dir"
    else
        ./scripts/build_rst_files.sh -e
    fi
else
    echo 'non execute'
    if [ -n "$target_dir" ]; then
        ./scripts/build_rst_files.sh "$target_dir"
    else
        ./scripts/build_rst_files.sh
    fi
fi
echo '--------------------------------'
echo
echo

echo '--------------------------------'
echo 'copying rst/png files and images'
if [ -n "$target_dir" ]; then
    ./scripts/update_rts_images.sh "$target_dir"
else
    ./scripts/update_rts_images.sh
fi
echo '--------------------------------'
echo
echo

if [ "$build" -eq 1 ]; then
    if ! command -v sphinx-build >/dev/null 2>&1; then
        echo "sphinx-build command not found" >&2
        exit 1
    fi
    echo '--------------------------------'
    echo 'running sphinx-build'
    sphinx-build -j 10 -b html ./ build
    echo '--------------------------------'
else
    echo 'not running sphinx-build'
fi
