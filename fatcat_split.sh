split ./fatcat_commands.txt "split_tmp_"
nb_cmd=$(cat fatcat_commands.txt | wc -l)

indexof()
{
    local word
    local item
    local idx

    word=$1
    shift
    item=$(printf '%s\n' "$@" | fgrep -nx "$word")
    let idx=${item%%:*}-1
    echo $idx
}

find ./results_fatcat -size 0 -delete

for f in ./"split_tmp_"*
do
    ##########
    # This block holds the background process count in check
    while [ $(jobs | wc -l) -ge 30 ]; do 
        sleep 0.1;

        files_count=$(find results_fatcat/ -type f -name "*txt" | wc -l)
        echo -ne $files_count "/" $nb_cmd" "
        echo -ne $(echo "scale=4; $files_count "/" $nb_cmd "*100 | bc)" "
        echo -ne $(echo "("$SECONDS "/3600)%24" | bc) $(echo "("$SECONDS "/60)%60" | bc)" "
        echo -ne $(echo $SECONDS "%60" | bc)" $(ls "./split_tmp_"* | grep -n $f )         \r"
    done; 

    ./fatcat_from_file.sh $f &
done
wait

rm ./"split_tmp_"*
sleep 1

# cd "./results_fatcat/"
# find . | xargs cat > ../fatcat_results.txt
# cd ..