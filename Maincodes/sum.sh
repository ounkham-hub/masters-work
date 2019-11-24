#!/bin/bash

##############################################################
#Script created to calculate the sum of each column in a file.
#This should be used after the sort-delta-OO-dist-v4.c program
#that outputs 8 files (pertaining to varying O-O distance
#ranges, 2.3 to 2.7 and delta, -0.6 to 0.6).
#
#Created by Lelee Ounkham
#March 10, 2017
#
##############################################################

#Command that sums every column in a file. Note! The output is 
#not in numerical order. The sum of columns 1,2,3 are outputted 
#at the end of the file.

for((i=1; i<=8; i++))
do

awk '{for(i=1;i<=NF;i++)sum[i]+=$i;};END{for(i in sum) print "for column "i" is "sum[i];}' map$i.txt > sumofmaps$i.txt

sort -k3n sumofmaps$i.txt -o sumofmaps$i.txt
done

#This entire output is printed to screen for convenience and can be 
#easily pasted onto Excel. From the outputs above, the sum of all
#columns are combined.

printf 'c.4  c.5  c.6  c.7  c.8  c.9  c.10  c.11  c.12  c.1  c.2  c.3 \n'
#prints a header to remind user that the columns are not in order

for((i=1; i<=8; i++))
do

awk 'BEGIN{ORS=" "} {if(NR==1) {print $5} else {print $5}}' sumofmaps$i.txt

printf '\n'
done

for((i=1; i<=4; i++))
do

awk '{for(i=1;i<=NF;i++)sum[i]+=$i;};END{for(i in sum) print "for column "i" is "sum[i];}' angle$i.txt > sumofangles$i.txt

sort -k3n sumofangles$i.txt -o sumofangles$i.txt

done


printf '\n'
#prints blank line between both results - easier readability

for((i=1; i<=4; i++))
do

awk 'BEGIN{ORS=" "} {if(NR==1) {print $5} else {print $5}}' sumofangles$i.txt

printf '\n'

done
