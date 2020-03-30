#! /bin/bash

target_dir_path="./www"

for file in images/*/*/*.jpg; do
        l1="${file%%/*}"
        l2="${file#*/}"
        l2="${l2%%/*}"
        l3="${file#*/*/}"
        l3="${l3%/*}" 
        filename="${file##*/}"
        target_file_name="${l1}_${l2}_${l3}_${filename}"
        cp "$file" "${target_dir_path}/${target_file_name}"
done
for file in images/*/*/*/*.jpg; do
        l1="${file%%/*}"
        l2="${file#*/}"
        l2="${l2%%/*}"
        l3="${file#*/*/}"
        l3="${l3%/*/*}"
        l4="${file#*/*/*/}"
        l4="${l4%/*}"        
        filename="${file##*/}"
        target_file_name="${l1}_${l2}_${l3}_${l4}_${filename}"
        cp "$file" "${target_dir_path}/${target_file_name}"
done

for file in images/*/*/*/*.png; do
        l1="${file%%/*}"
        l2="${file#*/}"
        l2="${l2%%/*}"
        l3="${file#*/*/}"
        l3="${l3%/*/*}"
        l4="${file#*/*/*/}"
        l4="${l4%/*}"        
        filename="${file##*/}"
        target_file_name="${l1}_${l2}_${l3}_${l4}_${filename}"
        cp "$file" "${target_dir_path}/${target_file_name}"
done