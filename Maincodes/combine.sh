#!/bin/bash

#for((i=1; i<=34004;i++))
#do
#sed -i -e "1d" H$i.xyz
#sed -i -e "1d" Cl$i.xyz
#removes the first line of a file

#sed -i 's/111/343/g' O$i.xyz 
#changes the atom number from "111" to "343"

#sed -i '/^$/d' H$i.xyz
#sed -i '/^$/d' Cl$i.xyz
#deletes all blank lines in a file

#cat O$i.xyz H$i.xyz Cl$i.xyz > HCl-snapshot$i.xyz
#merges all xyz files into a single snapshot

#done

cat $(find ./ -name "HCl-snapshot*.xyz" | sort -V) > wrapped-34k-HCl.xyz
#finds all files containing the word "HCl-snapshot" and merges them into a single file #in NUMERICAL order

