# echo $0 $1

cat $1 | while read i || [ -n "$i" ]; do
    first_str=$(echo $i | cut -d" " -f3)
    secnd_str=$(echo $i | cut -d" " -f5)
    fatcat_resfile="./results_fatcat/"$(echo "$first_str"_"$secnd_str".chain.txt | sed "s/\.pdb//g")

    if [[ -f $fatcat_resfile ]]
    then
        # echo "###"
        continue
    fi

    ##########
#     echo $i
#     $i | sed "s/  */ /g" > $fatcat_resfile
    $i > /dev/null

done