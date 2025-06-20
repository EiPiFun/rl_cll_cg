#!/usr/bin/sh

# comment

for i in $(ls ../*/*/inputs/coarse_grained_cellulose.in);do
sed -e "s|^fix 0 all langevin|#fix 0 all langevin|g" -i $i
done

for i in $(ls ../*/*/inputs/in_file_constant_tail.txt);do
sed -e "s|^fix 0 all langevin|#fix 0 all langevin|g" -i $i
done

# check

for i in $(ls ../*/*/inputs/coarse_grained_cellulose.in);do
cat $i | grep langevin
done

for i in $(ls ../*/*/inputs/in_file_constant_tail.txt);do
cat $i | grep langevin
done

