#!/usr/bin/sh

echo
echo 'timeout limit'
echo 'run this in vXY_physics'
echo

for i in $(ls -d ./physics_computation_*/inputs/);do
for j in coarse_grained_cellulose.in in_file_constant_head.txt in_file_constant_tail.txt;do
echo $i/$j

cat $i/$j | grep timeout

done
done

cat ./coarse_grained_cellulose_computation.py | grep timeout


