#!/bin/sh


for((i=1; i<=31136; i++))
do

#cp R-HCl.input.O$i.xyz.O$i.xyz.GraphGeod-updated set-10k-170-180

./post-select-angle.exe HCl.input.O$i.xyz.O$i.xyz.GraphGeod O$i.xyz H$i.xyz 150 180
done

#find . -name "*updated" -o -name "my*.txt" |xargs cat > total-2.6-3.0-GG.txt

#awk '{for(i=1;i<=NF;i++)sum[i]+=$i;};END{for(i in sum) print "for column "i" is "sum[i];}' 
