echo " You should run  ../scripts/batch_run.sh <in> <out> <mode> [ff=1] [nEv = 1000] [n=10] from the build directory" 
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
nev=1000
if [ "$#" -ge 5 ]; then
    nev=$5
fi
n=10
if [ "$#" -ge 6 ]; then
    n=$6
fi
echo "n=$n"
echo "====================================="
echo "====================================="
echo "====================================="
for i in $(seq $n); do 
    echo "============= batch  numner ${i} =========="
    ../scripts/run.sh ${in} ${out} ${mode} ${ff} $RANDOM ${nev}
done
mail -s "batch ${in} ${out} ${mode} ${ff} ${nev} ${n} is done" alexey.luchinsky@gmail.com < /dev/null

