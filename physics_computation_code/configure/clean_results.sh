#!/usr/bin/sh

for i in $(ls -d ../*/debug/);do
rm -r $i
done

for i in $(ls -d ../*/__pycache__/);do
rm -r $i
done

for i in $(ls -d ../*/*/results/);do
rm $i/*
done


