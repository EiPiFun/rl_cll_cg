#!/usr/bin/sh

rm ./v21_coarse_grained_cellulose_coefficients_env.py
cp ./v24_coarse_grained_cellulose_coefficients_env.py ./v21_coarse_grained_cellulose_coefficients_env.py
sed -e "/physics_version = / c physics_version = \'v21\'" -i ./v21_coarse_grained_cellulose_coefficients_env.py

rm ./v23_coarse_grained_cellulose_coefficients_env.py
cp ./v24_coarse_grained_cellulose_coefficients_env.py ./v23_coarse_grained_cellulose_coefficients_env.py
sed -e "/physics_version = / c physics_version = \'v23\'" -i ./v23_coarse_grained_cellulose_coefficients_env.py


