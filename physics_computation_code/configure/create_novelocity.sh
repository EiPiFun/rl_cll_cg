#!/usr/bin/sh

cd ../../../
echo $(pwd)

rm -rf ./coarse_grained_cellulose_coefficients_code_novelocity*
cp -r ./coarse_grained_cellulose_coefficients_code/ ./coarse_grained_cellulose_coefficients_code_novelocity/

cd ./coarse_grained_cellulose_coefficients_code_novelocity/physics_computation_code/configure/
sh ./modify_velocity_create_off.sh
cd ../../../
echo $(pwd)
zip -r -q ./coarse_grained_cellulose_coefficients_code_novelocity.zip ./coarse_grained_cellulose_coefficients_code_novelocity/

