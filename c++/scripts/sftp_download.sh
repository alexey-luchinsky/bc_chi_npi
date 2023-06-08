echo " You should run  ../scripts/stfp_download <mode> from the build directory"
echo " We have $# arguments"
if [ "$#" -lt 1 ]; then
    echo "Insufficient number of arguments"
    exit
fi
mode=$1
sed "s/OUT/${mode}/" ../templates/sftp.txt | tee sftp.txt
sftp aluchins@lxplus.cern.ch < sftp.txt


