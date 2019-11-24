#!/bin/bash

#Script created to organize GraphGeod files in ascending order of oxygen 
#indices (first column).
#
#Made by Lelee Ounkham
#October 13, 2016

for((i=1; i<=31136;i++))
do

sort -n -k1 HCl.input.O$i.xyz.H$i.xyz.GraphGeod > R-HCl.input.O$i.xyz.H$i.xyz.GraphGeod

done
