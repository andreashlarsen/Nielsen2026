
for p in b2 d1 d2 d22 kor2 hec23 hecate hecate23 atr
do
for rep in 1 2 3 4 5 6 7 8 9 10
do
    dir=${p}/rep${rep}/umb_${p}_bE_rep${rep}
    echo $dir
    if [ -d $dir ]
    then
        cd $dir
        taskset -c 8 gmx wham -it tpr_files_umbrella.dat -if pullf_files_umbrella.dat -hist -temp 320 -b 500 -bins 80 -min 1.3 -max 6.5 -quiet
        cd ../../..
    fi

    dir_idp=idp_${p}/CG_rep${rep}/umb_idp_${p}_bE_rep${rep}
    echo ${dir_idp}
    if [ -d $dir_idp ]
    then
        cd ${dir_idp}
        taskset -c 8 gmx wham -it tpr_files_umbrella.dat -if pullf_files_umbrella.dat -hist -temp 320 -b 500 -bins 80 -min 1.3 -max 6.5 -quiet
        cd ../../..
    fi

done
done

