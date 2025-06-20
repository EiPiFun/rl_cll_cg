#!/usr/bin/sh

cd ../../../
echo $(pwd)

rm -rf ./coarse_grained_cellulose_coefficients_code_nolangevin*
cp -r ./coarse_grained_cellulose_coefficients_code/ ./coarse_grained_cellulose_coefficients_code_nolangevin/

cd ./coarse_grained_cellulose_coefficients_code_nolangevin/physics_computation_code/configure/
sh ./turnoff_langevin.sh
cd ../../../
echo $(pwd)
zip -r -q ./coarse_grained_cellulose_coefficients_code_nolangevin.zip ./coarse_grained_cellulose_coefficients_code_nolangevin/


