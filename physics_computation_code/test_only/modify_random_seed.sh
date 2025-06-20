#!/usr/bin/bash

echo
echo "velocity create random seed"
echo "run this in vXY_physics"
echo

random_seed_index=0
if [ $# -gt 0 ];then
random_seed_index=$1
fi

random_seed_array=(114514 \
29189  1044 29337 21963  4225  5698 22285 28083 15229  9367 \
17119  9757 21552 12127 20845 28491 10627 25266 12789  3344 \
 3660 21026 12337  5195 19017 26809 24820 20112 26795 23751 \
14129 13353 20475  9001  7565  4223 32275 24015 14961  3742 \
24953  1678 25866 12310  7895 23610 10385 22951 29497  8021 \
17407 30926  8522 26622  3313 13239 15454 15924 32134 20642 \
 1766  2519 12910  5949  9227  2671 17634 11130 16093  3360 \
15549 25495  5944 31530 29113  2160 29635 15602  6334 32428 \
10230  6596  1431  2684  5864  4414 11177 30091  2466 13740 \
11205  6434  2827  3851  6301  2117  3273 27132 20392 24578)

random_seed=${random_seed_array[$random_seed_index]}

for i in $(ls -d ./physics_computation_*/inputs/);do
for j in coarse_grained_cellulose.in in_file_constant_head.txt in_file_constant_tail.txt;do
echo $i/$j

sed -e "/^velocity all create 300.0/ c velocity all create 300.0 $random_seed" -i $i/$j
cat $i/$j | grep velocity

done
done


