#!/usr/bin/sh

cc_name='gcc -lm'
if [ $# -gt 0 ];then
cc_name=$1
fi

$cc_name ../compute_end_to_end_distance_from_lammps_dump.c -o ../compute_end_to_end_distance_from_lammps_dump
chmod a+x ../compute_end_to_end_distance_from_lammps_dump
$cc_name ../compute_persistence_length_from_lammps_dump.c -o ../compute_persistence_length_from_lammps_dump
chmod a+x ../compute_persistence_length_from_lammps_dump

for i in 2 3 4 5;do
for j in 1 3;do
for k in 2 3 4;do

rm ../v"$i$j"_physics/physics_computation_"$k"/inputs/compute_end_to_end_distance_from_lammps_dump
cp ../compute_end_to_end_distance_from_lammps_dump ../v"$i$j"_physics/physics_computation_"$k"/inputs/compute_end_to_end_distance_from_lammps_dump
chmod a+x ../v"$i$j"_physics/physics_computation_"$k"/inputs/compute_end_to_end_distance_from_lammps_dump

rm ../v"$i$j"_physics/physics_computation_"$k"/inputs/compute_persistence_length_from_lammps_dump
cp ../compute_persistence_length_from_lammps_dump ../v"$i$j"_physics/physics_computation_"$k"/inputs/compute_persistence_length_from_lammps_dump
chmod a+x ../v"$i$j"_physics/physics_computation_"$k"/inputs/compute_persistence_length_from_lammps_dump

done
done
done

$cc_name ../compute_average_direction_angle_from_lammps_dump.c -o ../compute_average_direction_angle_from_lammps_dump
chmod a+x ../compute_average_direction_angle_from_lammps_dump

for i in 2 3 4 5;do
for j in 3;do
for k in 5 6 7;do

rm ../v"$i$j"_physics/physics_computation_"$k"/inputs/compute_average_direction_angle_from_lammps_dump
cp ../compute_average_direction_angle_from_lammps_dump ../v"$i$j"_physics/physics_computation_"$k"/inputs/compute_average_direction_angle_from_lammps_dump
chmod a+x ../v"$i$j"_physics/physics_computation_"$k"/inputs/compute_average_direction_angle_from_lammps_dump

done
done
done

rm ../compute_end_to_end_distance_from_lammps_dump
rm ../compute_persistence_length_from_lammps_dump
rm ../compute_average_direction_angle_from_lammps_dump


