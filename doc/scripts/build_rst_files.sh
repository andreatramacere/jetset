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

    rst_file="${file%.ipynb}.rst"
    if [ -f "$rst_file" ]; then
        # Avoid docutils strong-markup warnings caused by UltraNest progress
        # lines containing many '*' inside parsed-literal output blocks.
        perl -0777 -i -pe "s/\\.\\. parsed-literal::\\n\\n(\\s+====== ultranest script ========)/.. code-block:: text\\n\\n\$1/s" "$rst_file"

        # Keep OpenMPI UltraNest output compact in rendered docs by showing
        # only the first 20 lines of the huge stream block.
        awk '
        BEGIN {
            limit = 20
            seen_marker = 0
            cut_state = 0
            kept = 0
        }

        # Start truncation at the first UltraNest script output marker.
        /^    ====== ultranest script ========$/ && seen_marker == 0 {
            seen_marker = 1
            cut_state = 1
            kept = 1
            print
            next
        }

        # Keep only the first `limit` lines from the marked output block.
        cut_state == 1 {
            if (kept < limit) {
                print
                kept++
                next
            }
            print "    ... [output truncated: first " limit " lines shown] ..."
            print ""
            cut_state = 2
            next
        }

        # Skip the remaining indented output block lines.
        cut_state == 2 {
            if ($0 ~ /^    / || $0 ~ /^$/) {
                next
            }
            cut_state = 0
        }

        { print }
        ' "$rst_file" > "${rst_file}.tmp" && mv "${rst_file}.tmp" "$rst_file"
    fi
done
