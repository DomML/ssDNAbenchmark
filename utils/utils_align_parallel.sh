# utils_align_parallel.sh bound_list unbound_list boundchains_list unboundchains_list prefix suffix output nParal
### !python utils_align.py --bound $(cat "b_list.txt") --unbound $(cat "ub_list.txt") --prefix "./pdb/" --suffix "pdb" > align_analysis.txt

rm tmp_bound.split* tmp_unbound.split*

cat $1 | tr " " "\n" > tmp_bound
cat $2 | tr " " "\n" > tmp_unbound
cat $3 | tr " " "\n" > tmp_bound_ch
cat $4 | tr " " "\n" > tmp_unbound_ch

split -l$((`wc -l < tmp_bound`/$8))      tmp_bound      tmp_bound.split      -da $8
split -l$((`wc -l < tmp_unbound`/$8))    tmp_unbound    tmp_unbound.split    -da $8
split -l$((`wc -l < tmp_bound_ch`/$8))   tmp_bound_ch   tmp_bound_ch.split   -da $8
split -l$((`wc -l < tmp_unbound_ch`/$8)) tmp_unbound_ch tmp_unbound_ch.split -da $8

for i in tmp_*.split*; do
    cat $i | tr "\n" " " > tmp
    mv tmp $i
done

cur_dir=$(pwd)

for (( i=0; i<=$8; i++ )); do
    (
    mkdir tmp_"$i"
    cp tmp_bound.split*0$i      tmp_"$i"/
    cp tmp_unbound.split*0$i    tmp_"$i"/
    cp tmp_bound_ch.split*0$i   tmp_"$i"/
    cp tmp_unbound_ch.split*0$i tmp_"$i"/

    cd tmp_"$i"
    for j in tmp_bound.split*;      do mv $j bound.txt;      done
    for j in tmp_unbound.split*;    do mv $j unbound.txt;    done
    for j in tmp_bound_ch.split*;   do mv $j bound_ch.txt;   done
    for j in tmp_unbound_ch.split*; do mv $j unbound_ch.txt; done
    python ../utils_align.py --bound $(cat "bound.txt") --unbound $(cat "unbound.txt") --bound_chains $(cat "bound_ch.txt") --unbound_chains $(cat "unbound_ch.txt") --prefix "../pdb/" --suffix "pdb" > align_analysis.txt
    cd $cur_dir
    )&
done
wait


for (( i=0; i<=$8; i++ )); do
    cat tmp_"$i"/align_analysis.txt
    rm -rf tmp_"$i"
done | sort -n | uniq | grep -v "###" > $7

cat $7 | grep "error\|contain" > $7.err
cat $7 | grep -v "error\|contain" > tmp; mv tmp $7
cat $7 | sed "s/\#//g" > tmp; mv tmp $7

rm tmp_bound* tmp_unbound*
rm tmp_bound.split* tmp_unbound.split*