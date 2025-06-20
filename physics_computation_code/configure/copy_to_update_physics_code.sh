#!/usr/bin/sh

for i in 2 3 4 5;do

for j in 1 2 3 4;do
rm ../v"$i"3_physics/physics_computation_$j/inputs/*.data
done

#cp ../v"$i"3_physics/physics_computation_7/inputs/*.data ../v"$i"3_physics/physics_computation_1/inputs/
cp ../v"$i"1_physics/physics_computation_1/inputs/*.data ../v"$i"3_physics/physics_computation_1/inputs/

for j in 2 3 4;do
cp ../v"$i"1_physics/physics_computation_$j/inputs/*.data ../v"$i"3_physics/physics_computation_$j/inputs/
done

done


