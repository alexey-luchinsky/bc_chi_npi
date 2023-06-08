echo " You should run  ../scripts/_run.sh <in> <out> <mode> [ff=1] [seed = 0] [nEv = 100000] from the build directory"
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
out_dir=../results/${in}_${out}_${mode}_ff${ff}

echo "**********************"
echo "*  mode = ${mode}"
echo "*  in = ${in}"
echo "*  out = ${out}"
echo "*  ff=${ff}"
echo "*  seed=${seed}"
echo "*  nev = ${nev}"
echo "*  out_name = ${out_name}"
echo "*  out_dir = ${out_dir}"
echo "**********************"

sed "s/IN/${in}/ ; s/OUT/${out}/; s/FF/${ff}/" ../templates/${mode}.dec | tee tmp.dec

./bc_chiJ.exe  ${in} tmp.dec ${nev} $seed
~/Work/rrF/c++/build/rrF.exe -v "[`cat ../templates/${mode}_vars.txt`]" -o out_${out_name}.root -s
mkdir -p ${out_dir}
mv m*.txt ${out_dir}


exit

