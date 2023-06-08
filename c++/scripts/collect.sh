echo 'The format is ../scripts/collect.sh [pattern] [pattern] [pattern]'
echo 'If the parameter pattern is present, file list will be grepped with this pattern'

if [ -d /eos/user/a/aluchins/SWAN_projects/bc_new/results/tmp/ ]; then
    out_dir=/eos/user/a/aluchins/SWAN_projects/bc_7pi/results/
else
    out_dir=../results_scripts/
fi
ls -l ${out_dir}/tmp/ | grep root | grep out | sed "s/[^_]*_\(.*\)_seed.*/\1/" | sort | uniq > list.txt
if [ "$#" -ge 1 ]; then
    echo 'Filtering the list1'
    cat list.txt | grep $1 > A.txt
    mv A.txt list.txt
fi
if [ "$#" -ge 2 ]; then
    echo 'Filtering the list2'
    cat list.txt | grep $2 > A.txt
    mv A.txt list.txt
fi
if [ "$#" -ge 3 ]; then
    echo 'Filtering the list3'
    cat list.txt | grep $3 > A.txt
    mv A.txt list.txt
fi
cat list.txt
for f in `cat list.txt`; do
    echo $f
    hadd -f ${out_dir}/${f}.root ${out_dir}/tmp/out_${f}_*.root
done
