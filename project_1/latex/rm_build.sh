#!/bin/bash

for f in $(ls); do
    f=$(basename -- "$f")
    if [ "$f" = "tikz" ] || [ "$f" = "tags" ] || [ "$f" = $(basename -- "$0") ]; then
        continue
    fi
    e="${f##*.}"
    if [ "$e" != "tex" ]; then
        echo "rm -f $f"
        rm -f $f
    fi
done

# vim: set ff=unix tw=79 sw=4 ts=4 et ic ai :
