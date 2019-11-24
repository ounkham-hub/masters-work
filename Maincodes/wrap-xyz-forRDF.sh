#/bin/bash

var="308"

for((i=1; i<55001; i++))
do
tail -n +3 H$i.xyz > r-H$i.xyz
done

for((i=1; i<55001; i++))
do
tail -n +3 Cl$i.xyz > r-Cl$i.xyz
done

for((i=1; i<55001; i++))
do
cat O$i.xyz r-H$i.xyz r-Cl$i.xyz > wrapped-$i.xyz
done

for((i=1; i<55001; i++))
do
sed -i "1s/.*/$var/" wrapped-$i.xyz
#the number 1 denotes the first line and the rest is substitution for everything (.*)
#in the first line by $var, the defined variable
done

cat $(find ./ -name "wrapped*.xyz" | sort -V) > 2m-wrapped-RDF.xyz
#finds all files containing the word "HCl-snapshot" and merges them into a single file #in NUMERICAL order

