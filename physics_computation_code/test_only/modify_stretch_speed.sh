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

sed -e '/run 5000/ c run 50000' -i  $i/$j
sed -e '/run 10000/ c run 100000' -i  $i/$j
sed -e '/fix stretch all deform 1 y vel/ c fix stretch all deform 1 y vel 0.00001' -i  $i/$j
sed -e '/run 20000/ c run 200000' -i  $i/$j
sed -e '/run 60000/ c run 600000' -i  $i/$j
sed -e '/dump 2 cl custom 1000/ c dump 2 cl custom 1000 cl_stretch id mol type q x y z fx fy fz' -i $i/$j
sed -e '/fix boxsizepressurestretch all ave\/time 1000/ c fix boxsizepressurestretch all ave/time 1000 1 1000 v_lx v_ly v_lz v_px v_py v_pz file box_size-pressure_stretch mode scalar' -i $i/$j

cat $i/$j | grep 'run '
cat $i/$j | grep deform
cat $i/$j | grep cl_stretch
cat $i/$j | grep boxsizepressurestretch

done
done


