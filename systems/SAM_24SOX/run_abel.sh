#!/bin/bash
set -e

mkdir -p log/
declare -a arr=("conf1" "conf10" "conf20" "conf30" "conf40" "conf50"
                "conf60" "conf70" "conf80" "conf90" "conf100"
		"conf110" "conf120" "conf130" "conf140" "conf150")

for i in "${arr[@]}"
do
    sed "s/CONF/$i/g" abel.sh > run.tmp
    sbatch run.tmp
    rm run.tmp
done

exit 1

