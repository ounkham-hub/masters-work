#!/bin/bash

#made by L.Edens
#Modified by Lelee Ounkham
#October 11, 2017

#snapshot
t=20173

#Lines:
#OCO1  = 97
#H1CB1 = 161
#H1CB2 = 200
#OCO2  = 78
#H2CB1 = 140
#H2CB2 = 154
#shared H = 191

for ((i=1; i<=17; i++))
do
   
   cd /scratch/lelee/posi-$i
 
   awk 'NR==97 {printf "%s   %.3f %.3f %.3f\n", (($1)),(($2)),(($3)),(($4))}' solB$t.xyz >> ../32_step$t.xyz
   awk 'NR==78 {printf "%s   %.3f %.3f %.3f\n", (($1)),(($2)),(($3)),(($4))}' solB$t.xyz >> ../32_step$t.xyz
   awk 'NR==161 {printf "%s   %.3f %.3f %.3f\n", (($1)),(($2)),(($3)),(($4))}' solC$t.xyz >> ../32_step$t.xyz
   awk 'NR==200 {printf "%s   %.3f %.3f %.3f\n", (($1)),(($2)),(($3)),(($4))}' solC$t.xyz >> ../32_step$t.xyz
   awk 'NR==140 {printf "%s   %.3f %.3f %.3f\n", (($1)),(($2)),(($3)),(($4))}' solC$t.xyz >> ../32_step$t.xyz
   awk 'NR==154 {printf "%s   %.3f %.3f %.3f\n", (($1)),(($2)),(($3)),(($4))}' solC$t.xyz >> ../32_step$t.xyz
   awk 'NR==191 {printf "%s   %.3f %.3f %.3f\n", (($1)),(($2)),(($3)),(($4))}' solC$t.xyz >> ../32_step$t.xyz

done
