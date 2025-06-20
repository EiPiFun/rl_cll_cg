#!/usr/bin/sh

rm ./v41_coarse_grained_cellulose_coefficients_env.py
cp ./v44_coarse_grained_cellulose_coefficients_env.py ./v41_coarse_grained_cellulose_coefficients_env.py
sed -e "/physics_version = / c physics_version = \'v41\'" -i ./v41_coarse_grained_cellulose_coefficients_env.py

rm ./v43_coarse_grained_cellulose_coefficients_env.py
cp ./v44_coarse_grained_cellulose_coefficients_env.py ./v43_coarse_grained_cellulose_coefficients_env.py
sed -e "/physics_version = / c physics_version = \'v43\'" -i ./v43_coarse_grained_cellulose_coefficients_env.py


