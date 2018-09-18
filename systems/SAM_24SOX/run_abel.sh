declare -a arr=("conf0" "conf1")

for i in "${arr[@]}"
do
   echo "./abel.sh $i"
done

