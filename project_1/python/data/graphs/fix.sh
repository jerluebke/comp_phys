#!/bin/bash

# to get max extent, look picture with largest resolution, i.e.
#     for $f in *.png; do
#         file $f
#     done

rm extent*
for i in {001..063}; do
    f=out-$i.png
    s=$(sed "${i}q;d" key.file)
    convert $f -gravity center -extent 3617x827 tmp_$f
    montage -label "$s" tmp_$f -pointsize 64 -geometry +0+0 extent_$f
    rm tmp_$f
done

# vim: set ff=unix tw=79 sw=4 ts=4 et ic ai : 
