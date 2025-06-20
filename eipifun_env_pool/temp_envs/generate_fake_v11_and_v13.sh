#!/usr/bin/sh

rm ./v11_coarse_grained_cellulose_coefficients_env.py
cp ./v14_coarse_grained_cellulose_coefficients_env.py ./v11_coarse_grained_cellulose_coefficients_env.py
sed -e "/physics_version = / c physics_version = \'v11\'" -i ./v11_coarse_grained_cellulose_coefficients_env.py

rm ./v13_coarse_grained_cellulose_coefficients_env.py
cp ./v14_coarse_grained_cellulose_coefficients_env.py ./v13_coarse_grained_cellulose_coefficients_env.py
sed -e "/physics_version = / c physics_version = \'v13\'" -i ./v13_coarse_grained_cellulose_coefficients_env.py


