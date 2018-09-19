#!/bin/bash
set -e

declare -a arr=("conf0" "conf1"
                "conf2")

for i in "${arr[@]}"
do
   echo "./abel.sh $i"
done

exit 1

