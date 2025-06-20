#!/usr/bin/sh

echo
echo "turn on suffix"
echo "run this in vXY_physics"
echo

for i in $(ls -d ./physics_computation_*/inputs/);do
for j in coarse_grained_cellulose.in in_file_constant_head.txt in_file_constant_tail.txt;do
echo $i/$j

sed -e "/suffix/ c suffix omp" -i $i/$j
cat $i/$j | grep suffix

done
done


