#!/usr/bin/sh

rm ./v31_coarse_grained_cellulose_coefficients_env.py
cp ./v34_coarse_grained_cellulose_coefficients_env.py ./v31_coarse_grained_cellulose_coefficients_env.py
sed -e "/physics_version = / c physics_version = \'v31\'" -i ./v31_coarse_grained_cellulose_coefficients_env.py

rm ./v33_coarse_grained_cellulose_coefficients_env.py
cp ./v34_coarse_grained_cellulose_coefficients_env.py ./v33_coarse_grained_cellulose_coefficients_env.py
sed -e "/physics_version = / c physics_version = \'v33\'" -i ./v33_coarse_grained_cellulose_coefficients_env.py


