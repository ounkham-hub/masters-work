#!/bin/bash

#Script modified to determine the edge distribution using GraphGeod files 
#instead of Graph files. GraphGeod files have oxygen atoms listed in one 
#column and hydrogen/pair in the second.
#Made by Lelee Ounkham
#August 5, 2016

for((i=1; i<14; i++))
do

#Go to posi-# directory
cd /scratch/lelee/posi-$i

for((j=1; j<=61555; j++))
do

awk -F " " '{print $1}' R-HCl.input.solB$j.xyz.solC$j.xyz.GraphGeod | uniq -c > count-Hedges$j.txt
#counts the frequency of edges associated with each O indice

done

find . -name "count-Hedges*" -o -name "my*.txt" |xargs cat > total-edges-perOxy.txt
#finds all files containing the word "count" and merges them into a single 
#file

awk '{h[$1]++} END { for(k in h) print k, h[k] }' total-edges-perOxy.txt | sort -n > CB-edge-distribution.txt

for((j=1; j<=61555; j++))
do

awk -F " " '{print $2}' R-HCl.input.solB$j.xyz.solC$j.xyz.GraphGeod | uniq -c > count-Oedges$j.txt
#counts the frequency of edges associated with each O indice

done
#end of j-loop
find . -name "count-Oedges*" -o -name "my*.txt" |xargs cat > total-edges-perHyd.txt
#finds all files containing the word "count" and merges them into a single 
#file

awk '{h[$1]++} END { for(k in h) print k, h[k] }' total-edges-perHyd.txt | sort -n > CB-edge-distribution2.txt

done
#end of i-loop
