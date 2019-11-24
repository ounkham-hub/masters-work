#!/bin/bash
#
#Script created to find the associated xyz coordinates
#for over-coordinated oxygen, Zundel and Eigen outputs to 
#visualize in VMD.
#
#Made by Lelee Ounkham
#October 27, 2016

for((i=2633; i<=5033;i++))
do

#Delete the first # LINES of header in the O.xyz and H.xyz files
#sed -i '1d' 1-eigen$i.txt 
#sed -i '1d' 3-eigen$i.txt
#sed -i 1,2d 1-eigen$i.txt 
#sed -i 1,2d Cl$i.xyz

#Concatenate O.xyz and H.xyz together
#cat O$i.xyz H$i.xyz Cl$i.xyz >> cOH$i.xyz

#An alternative to deleting lines, one that preserves the original file
#tail -n +1 1-eigen$i.txt > eigen1$i.txt

#Add number of rows to file, used to find pattern in Zunden/Eigen output
#nl cOH$i.xyz > O-H$i.xyz 

#Compares files zundel-2m.txt n-O-H-.xyz and prints the matching 
#line in O-H.xyz

#awk 'FNR==NR{a[$1];next}($1 in a){print}' 1-eigen$i.txt O-H$i.xyz > 1-traj$i.xyz

#Removes the first column of zundel-2m.txt, otherwise the above awk command 
#doesn't work. The commands area repeated until the end of the file. 
awk '!($1="")' 1HB-traj$i.xyz > 1HB$i.txt

#awk 'FNR==NR{a[$1];next}($1 in a){print}' 1z$i.txt O-H$i.xyz >> 1-traj$i.xyz

#awk '!($1="")' 1HB$i.txt > 2HB$i.txt

#awk 'FNR==NR{a[$1];next}($1 in a){print}' 2HB$i.txt CB-HCl.input.O$i.xyz.H$i.xyz.GraphGeod >> 1HB-traj$i.xyz

#awk '!($1="")' 2z$i.txt > 3z$i.txt

#awk 'FNR==NR{a[$1];next}($1 in a){print}' 3z$i.txt O-H$i.xyz >> 1-traj$i.xyz

#awk '!($1="")' 3z$i.txt > 4z$i.txt

#awk 'FNR==NR{a[$1];next}($1 in a){print}' 4z$i.txt O-H$i.xyz >> 1-traj$i.xyz

#Removes the first column of traj (the number of lines)
#awk '!($1="")' 1-traj$i.xyz > 1eigen-8m-traj$i.txt

done
#rm 1z*.txt
#rm 2z*.txt
#rm 3z*.txt
#rm 4z*.txt
#rm traj*.xyz
