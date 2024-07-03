
#!/bin/bash 
. ./projectProperties.conf
#script gets three arguments for mode, xseed, yseed
mode=$1
xSeed=$2
ySeed=$3
size=$4
if [ ! $# -eq 4 ]; then
	printf "\n--: There should be 4 arguments: Mode, xSeed, ySeed, Size \n"
	printf "*********************************************\n"
	exit
fi

ext=""
sizeH=$size
if [ $sizeH -ge 1024 ];then
	sizeH=$(($sizeH/1024))
	ext="k"
fi

if [ $sizeH -ge 1024 ];then
	sizeH=$(($sizeH/1024))
	ext="m"
fi

directoryName="data_"$mode"_"$xSeed"_"$ySeed"_"$sizeH$ext
printf "\nGenerating dataset: $directoryName\n"
# PROJECT_NAME="BigPrimeFieldFFT_3"
# PROJECT=$CUMODP_HOME/new_tests/$PROJECT_NAME


mkdir $PROJECT/test/data/$directoryName 2>/dev/null
# printf "\n++MAKING DIRECTORY "$directoryName"\n"
# printf "***********************\n"
cd $PROJECT/bin
rm -f vector
# printf "\n++GENERATING vector x\n"
./generateData.bin $mode $xSeed $size > /dev/null
mv -f vector $PROJECT/test/data/$directoryName/xData

# printf "\n++GENERATING vector y\n"
./generateData.bin $mode $ySeed $size > /dev/null

mv -f vector $PROJECT/test/data/$directoryName/yData
printf "\nGenerated vectors --> /test/data/$directoryName \n\n"
printf "*********************************************\n"
