#!/usr/bin/sh

rm ./v51_coarse_grained_cellulose_coefficients_env.py
cp ./v54_coarse_grained_cellulose_coefficients_env.py ./v51_coarse_grained_cellulose_coefficients_env.py
sed -e "/physics_version = / c physics_version = \'v51\'" -i ./v51_coarse_grained_cellulose_coefficients_env.py

rm ./v53_coarse_grained_cellulose_coefficients_env.py
cp ./v54_coarse_grained_cellulose_coefficients_env.py ./v53_coarse_grained_cellulose_coefficients_env.py
sed -e "/physics_version = / c physics_version = \'v53\'" -i ./v53_coarse_grained_cellulose_coefficients_env.py


