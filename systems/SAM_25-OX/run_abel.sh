#!/bin/bash
set -e

declare -a arr=("conf1" "conf10" "conf20" "conf30" "conf40" "conf50"
                "conf60" "conf70" "conf90" "conf100")

for i in "${arr[@]}"
do
    sed "s/CONF/$i/g" abel.sh > run.tmp
    sbatch run.tmp
    rm run.tmp
done

exit 1

