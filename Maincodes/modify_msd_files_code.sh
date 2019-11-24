#!/bin/bash

#Script to change a pattern in specific lines of a file
#Made by E_Martinez-Baez
#Novemeber, 16, 2016

#original file name
input="tbp-com.txt"
output="tbp-com-modified.txt"
temp_backup="$input.bak"

sed -i.bak '/nan/d' ./$input	#this eliminates all lines containing 'nan' in them 
				#and saves a backup in .bak so the new file name will
 			 	#be the same as $input 
mv $input $output		#changes the resulting file name of previous step to
				#the name assiged to $output
mv $temp_backup $input 		#changes the resulting backup file name of previous 
				#step to the name of the original file $input

#original file name
input2="tbp-com-modified.txt"
output2="tbp-com-modified1.txt"
line="4"		#line where the sed command is gonna start acting  
step="43"		#number of lines ingnored before excuting the command again
var1="3540"		#pattern I want to change in the line selected
var2="42"		#pattern I want to be written instead of var1 


sed -e "${line}~${step}s/${var1}/${var2}/g" $input2 > $output2


