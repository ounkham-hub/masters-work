#!/bin/bash

#Script created to count the number of OC-O in a single snapshot and 
#take the average.
#
#Made by Lelee Ounkham
#October 17, 2016

#for((i=1; i<=34004; i++))
#do
#wc -l OC-O-list$i.txt > count-O$i-intra.txt
#done

cat $(find ./ -name "count*.txt" | sort -V) > total-count-O-intra.txt
#finds all files containing the word "count.txt" and merges them into a single 
#file in NUMERICAL order

awk '{h[$1]++} END { for(k in h) print k, h[k] }' total-count-O-intra.txt | sort -n > OC-ED.txt
