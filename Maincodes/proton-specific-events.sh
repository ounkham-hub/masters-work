#!/bin/bash

#######################################################
# Script created to organize PT events based on 
# the shared H index.
#
# Made by Lelee Ounkham
# Last Modified March 24, 2018
########################################################

   for((i=1; i<=1; i++))
   do
      
      for((j=103; j<=314; j++))
      do

	awk -v var1="$j" '$1==var1' rep-$i-zundels.txt > $j-Zundel-switches.txt
        
      done
   done

