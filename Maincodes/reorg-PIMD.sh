#!/bin/bash

##################################################################
# Script created to create trajectories from PIMD data such that
# each atom has information on it's corresponding 32 neighbors 
# in a given snapshot.
#
# Made by: Lelee Ounkham
# Last Modified: November 3, 2017
#
################################################################## 


#for((i=1; i<=18; i++))
#do
#   cd /scratch/lelee/posi-$i

   for((j=1; j<=1; j++))
   do

      awk '{print NR " " $0}' solB$j.xyz > labeled-O-$j.xyz #O-atom
#      awk '{print NR " " $0}' solC$j.xyz > labeled-H-$j.xyz #H-atom
      #Add line number to add # index to atom for every single snapshot
   done
   #end j-loop
#done
#end i-loop


#for((i=1; i<=18; i++))
#do
#   cd /scratch/lelee/posi-$i

####For O atoms###
   for((j=1; j<=1; j++))
   do
      for((k=3; k<=3; k++)) 
      do
#         awk 'NR==$k {printf "%s   %.3f %.3f %.3f\n", (($1)),(($2)),(($3)),(($4))}' solB$t.xyz >> ../.xyz
	 awk -vk=$k 'NR==k {printf "%d %s     %f %f %f\n", (($1)),(($2)),(($3)),(($4)),(($5))}' labeled-O-$j.xyz >> 32rep-O$k-snap$j.xyz
#for every snapshot, print out associated coordinates per atom into unique files - the end result is a trajectory
#containing spatial and time information for a single atom
      done
      # end k-loop
   done
   #end j-loop

###For H atoms####
#   for((j=1; j<=61555; j++))
#   do
#      for((k=1; k<=212; k++)) 
#      do
#         awk 'NR==$k {printf "%s   %.3f %.3f %.3f\n", (($1)),(($2)),(($3)),(($4))}' solB$t.xyz >> ../.xyz
#	 awk -vk=$k 'NR==k {printf "%d %s   %.3f %.3f %.3f\n", (($1)),(($2)),(($3)),(($4)),(($5))}' labeled-H-$k.xyz >> ../32rep-H$k-snap$j.xyz
#for every snapshot, print out associated coordinates per atom into unique files - the end result is a trajectory
#containing spatial and time information for a single atom
#      done
      # end k-loop
#   done
   #end j-loop

#done
#end i-loop

