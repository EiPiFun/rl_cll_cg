#!/usr/bin/sh

for i in 2 3 4 5;do
for j in 1 3;do
rm -r ../v"$i$j"_physics/coarse_grained_cellulose_post_processing.py
cp ../coarse_grained_cellulose_post_processing.py ../v"$i$j"_physics/
done
done
