#!/usr/bin/sh

echo
echo "vX1 vX2 vX3 vX4"
echo

for i in 2 3 4 5;do
for j in 1 3;do
for k in 1 2 3 4;do
for l in coarse_grained_cellulose.in in_file_constant_head.txt in_file_constant_tail.txt;do

sed -e "s|dump 1 cl custom 100 cl_relaxation|dump 1 cl custom 1000 cl_relaxation|g" -i ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l
sed -e "s|dump 2 cl custom 100 cl_stretch|dump 2 cl custom 1000 cl_stretch|g" -i ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l
sed -e "s|fix clmsd cl ave/time 100 1 100|fix clmsd cl ave/time 1000 1 1000|g" -i ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l
sed -e "s|fix boxsizepressurerelaxation all ave/time 100 1 100|fix boxsizepressurerelaxation all ave/time 1000 1 1000|g" -i ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l
sed -e "s|fix boxsizepressurestretch all ave/time 100 1 100|fix boxsizepressurestretch all ave/time 1000 1 1000|g" -i ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l
cat ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l | grep dump | grep custom
cat ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l | grep "fix clmsd "
cat ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l | grep "fix boxsizepressurerelaxation "
cat ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l | grep "fix boxsizepressurestretch "

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

sed -e "s|dump 1 cl custom 100 cl_relaxation|dump 1 cl custom 1000 cl_relaxation|g" -i ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l
sed -e "s|dump 2 cl custom 100 cl_stretch|dump 2 cl custom 1000 cl_stretch|g" -i ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l
sed -e "s|fix clmsd cl ave/time 100 1 100|fix clmsd cl ave/time 1000 1 1000|g" -i ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l
sed -e "s|fix boxsizepressurerelaxation all ave/time 100 1 100|fix boxsizepressurerelaxation all ave/time 1000 1 1000|g" -i ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l
sed -e "s|fix boxsizepressurestretch all ave/time 100 1 100|fix boxsizepressurestretch all ave/time 1000 1 1000|g" -i ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l
cat ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l | grep dump | grep custom
cat ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l | grep "fix clmsd "
cat ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l | grep "fix boxsizepressurerelaxation "
cat ../v"$i$j"_physics/physics_computation_"$k"/inputs/$l | grep "fix boxsizepressurestretch "

done
done
done
done


