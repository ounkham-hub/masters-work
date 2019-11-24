#!/bin/bash

for((i=7; i<=7; i++))
do

#Rename all O indices to X and H indices to Y
sed -i -e "s/O/X/g" orig-position-$i.xyz
sed -i -e "s/H/Y/g" orig-position-$i.xyz

#Change total atom count to account for added water molecule
sed -i 's/322/325/g' orig-position-$i.xyz

#Add a water molecule to end of every snapshot (every 324th line)
# in a trajectory file
sed 'n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;n;G;' orig-position-$i.xyz > output1.xyz
#By inserting a blank line every 324th line then writing an xyz coordinate

sed 's/^$/\       \O  0.00000e+00\ \ 0.00000e+00\ \ 0.00000e+00\n/' output1.xyz > output2.xyz

sed 's/^$/\       \H  0.00000e+00\ \ 0.00000e+00\ \ 0.00000e+00\n/' output2.xyz > output3.xyz

sed 's/^$/\       \H  0.00000e+00\ \ 0.00000e+00\ \ 0.00000e+00\n/' output3.xyz > output1.xyz

sed '/^$/d' output1.xyz > added-position-$i.xyz
#deletes any remaining/unused blank lines in the file

#This output corresponds to the executable wrapping code with defined
#input parameters
./a.out < input$i.txt

##################################################################
#Run ChemNetworks to create GraphGeod files to do analyses

for((j=1; j<=61556; j++))
do

#Rename all X indices to O and Y indices to H
sed -i -e "s/X/O/g" solB$j.xyz
sed -i -e "s/Y/H/g" solC$j.xyz

#HCl.input is set to print GraphGeods of O-O from 2.3 - 2.7A
./ChemNetworks-2.2.exe HCl.input solB$j.xyz solC$j.xyz solD$j.xyz

done
#####################################################################

#Start of analyses program
#1) O-O vs. delta
#2) O-O vs. O-H vs. angle

#Analysis #1 O-O vs. delta
#####################################################################
for((j=1; j<=61556; j++))
do

#Tiecheng's program that outputs the O-H distances and O-H-O angles from 0-180
./post-select-angle.exe HCl.input.solB$j.xyz.solB$j.xyz.GraphGeod solB$j.xyz solC$j.xyz 0 180

#Command used to take the difference between O1-H and O2-H, also known as delta
awk '{ $8 = $5 - $6 } 1' HCl.input.solB$j.xyz.solB$j.xyz.GraphGeod-updated > OH-diff$j.GG

done

#Compile sort program
gcc -lm sort-delta-OOdist-GG-v4.c -o sort-delta.exe

#Execute program, that will produce 8 outputs corresponding 8 O-O range
./sort-delta.exe

#Analysis #2 O-O vs. O-H vs. angle
############################################################

for((k=1; k<=61556; k++))
do

#Slight modification to get outputs containing O-H-O angles from 150-180
./select-150-180.exe HCl.input.solB$k.xyz.solB$k.xyz.GraphGeod solB$k.xyz solC$k.xyz 150 180

done

#compile sort program 
gcc ext-sort-OO-OH-dist.c -o sort-OO-angle.exe

#Execute program which will produce 4 outputs map2.3-2.4,map2.4-2.5,map2.5-2.6,map2.6-2.7
./sort-OO-angle.exe

#Makes a directory named posi-$i to store all corresponding filess
mkdir posi-$i

#Moves all xyz files into a single folder
mv sol* posi-$i
mv *GraphGeod posi-$i
mv map* posi-$i
mv angle* posi-$i
mv *angle posi-$i
mv *updated posi-$i
mv OH* posi-$i

done

rm *Graph
rm output*
rm water*
