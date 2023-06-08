echo " You should run  ../scripts/run.sh <in> <out> <mode> [ff=1] [seed = 0] [nEv = 100000] from the build directory"
echo " We have $# arguments"
if [ "$#" -lt 3 ]; then
    echo "Insufficient number of arguments"
    exit
fi
in=$1
out=$2
mode=$3
if [ ! -f ../templates/${mode}.dec ]; then
    echo "DECAY file ../templates/${mode}.dec does not exist, mode is not realized"
    exit
fi
ff=1
if [ "$#" -ge 4 ]; then
    ff=$4
fi
seed=0
if [ "$#" -ge 5 ]; then
    seed=$5
fi
nev=100000
if [ "$#" -ge 6 ]; then
    nev=$6
fi
out_name=${in}_${out}_${mode}_ff${ff}_seed${seed}

if [ -d /eos/user/a/aluchins/SWAN_projects/ ]; then
    out_dir=/eos/user/a/aluchins/SWAN_projects/bc_7pi/results/tmp/
else
    out_dir=../results_scripts/tmp/
fi
printf "============== Saving file to ${out_dir} ============ \n"
../scripts/_run.sh $1 $2 $3 $4 $5 $6 | tee log_${out_name}.txt
mv out_${out_name}.root ${out_dir}
mv log_${out_name}.txt ${out_dir}
