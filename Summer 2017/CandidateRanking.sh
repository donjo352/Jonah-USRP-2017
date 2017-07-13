#!/bin/bash

files=("HTR111-006" "HTR250-006")
for ((i=0; i < ${#files[@]}; i++))
do  
	htr ${files[$i]} > "$i" + ".txt"
done