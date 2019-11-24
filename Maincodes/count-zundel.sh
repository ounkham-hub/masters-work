#!/bin/bash

#Script created to count the number of OC-O in a single snapshot and 
#take the average.
#
#Made by Lelee Ounkham
#October 17, 2016

for((i=1; i<=34004; i++))
do
wc -l zundel$i-2.5m.txt > count-zundel$i.txt
done

cat $(find ./ -name "count-zundel*.txt" | sort -V) > total-zundel-count-2.5m.txt
#finds all files containing the word "count-zundel.txt" and merges them into a single 
#file in NUMERICAL order

awk '{h[$1]++} END { for(k in h) print k, h[k] }' total-zundel-count-2.5m.txt | sort -n > zundel-ED-2.5m.txt
