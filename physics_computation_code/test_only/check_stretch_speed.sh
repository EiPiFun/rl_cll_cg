#!/usr/bin/sh

echo
echo 'stretch speed'
echo 'run this in vXY_physics'
echo

for i in $(ls -d ./physics_computation_*/inputs/);do
for j in coarse_grained_cellulose.in in_file_constant_head.txt in_file_constant_tail.txt;do
echo $i/$j

echo 'preset stretch speed is 0.0001'
echo 'please check timeout settings'

cat $i/$j | grep 'run '
cat $i/$j | grep deform
cat $i/$j | grep cl_stretch
cat $i/$j | grep boxsizepressurestretch

done
done


