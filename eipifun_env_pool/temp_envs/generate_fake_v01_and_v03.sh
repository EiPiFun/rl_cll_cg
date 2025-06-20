#!/usr/bin/sh

rm ./v01_coarse_grained_cellulose_coefficients_env.py
cp ./v04_coarse_grained_cellulose_coefficients_env.py ./v01_coarse_grained_cellulose_coefficients_env.py
sed -e "/physics_version = / c physics_version = \'v01\'" -i ./v01_coarse_grained_cellulose_coefficients_env.py

rm ./v03_coarse_grained_cellulose_coefficients_env.py
cp ./v04_coarse_grained_cellulose_coefficients_env.py ./v03_coarse_grained_cellulose_coefficients_env.py
sed -e "/physics_version = / c physics_version = \'v03\'" -i ./v03_coarse_grained_cellulose_coefficients_env.py


