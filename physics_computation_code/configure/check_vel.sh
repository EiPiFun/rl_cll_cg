#!/usr/bin/sh

echo
echo "vX1 vX2 vX3 vX4"
echo

for i in 2 3 4 5;do
for j in 1 3;do
for k in 1 2 3 4;do
for l in coarse_grained_cellulose.in in_file_constant_head.txt in_file_constant_tail.txt;do

cat ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l | grep ' vel ';echo

done
done
done
done

echo
echo "vX5 vX6 vX7"
echo

for i in 2 3 4 5;do
for j in 3;do
for k in 5 6 7;do
for l in coarse_grained_cellulose.in in_file_constant_head.txt in_file_constant_tail.txt;do

cat ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l | grep ' vel ';echo

done
done
done
done


