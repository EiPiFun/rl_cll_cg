#!/usr/bin/sh

echo
echo 'timeout limit'
echo 'run this in vXY_physics'
echo

for i in $(ls -d ./physics_computation_*/inputs/);do
for j in coarse_grained_cellulose.in in_file_constant_head.txt in_file_constant_tail.txt;do
echo $i/$j

sed -e '/timer timeout/ c timer timeout 0:29:00 every 1000' -i  $i/$j

cat $i/$j | grep timeout

done
done

sed -e '/coarse_grained_cellulose_computation_timeout =/ c coarse_grained_cellulose_computation_timeout = 1800.0' -i  ./coarse_grained_cellulose_computation.py
cat ./coarse_grained_cellulose_computation.py | grep timeout


