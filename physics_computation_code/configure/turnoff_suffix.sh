#!/usr/bin/sh

# comment

for i in $(ls ../*/*/inputs/coarse_grained_cellulose.in);do
sed -e '/suffix/ c suffix off' -i $i
sed -e 's|^package|#package|g' -i $i
done

for i in $(ls ../*/*/inputs/in_file_constant_head.txt);do
sed -e '/suffix/ c suffix off' -i $i
sed -e 's|^package|#package|g' -i $i
done

# check

for i in $(ls ../*/*/inputs/coarse_grained_cellulose.in);do
cat $i | grep suffix
cat $i | grep package
done

for i in $(ls ../*/*/inputs/in_file_constant_head.txt);do
cat $i | grep suffix
cat $i | grep package
done

