
#!/bin/bash

#four argument, operation, data mode+xSeed+ySeed

# PROJECT_NAME="BigPrimeFieldFFT_4"
# PROJECt=$CUMODP_HOME/new_tests/$PROJECT_NAME

#########################################
source projectProperties.conf
#########################################
operation=$1
iteration=$2
shuffle=$3
padding=$4
dataSet=$5

if [ ! $# -eq 5 ]; then
	printf "\n--: There should be 5 arguments: operation, iteration, shuffle, padding, dataSet\n"
	printf "\n--: 0: Addition,\t 1: Multiplication,\t 2: Subtraction,\t 3: CyclicShift"
	printf "\n--: 4: FFT_16,\t\t 5: FFT_256,\t\t 6: FFT_4096,\t\t 7:FFT_16K \n"
	printf "********************************\n"
	exit
fi

dataPath=$PROJECt/test/data/$dataSet

echo $dataPath
if [ ! -d $dataPath ];then
	printf "********************************\n"
	printf "\n--: $dataSet doesn't exist!\n"
	printf "********************************\n"
	exit
fi

cd $dataPath
if [ ! -e "uData_cuda_"$operation"_"$iteration"_"$shuffle"_"$padding ];then
	printf "********************************\n"
	printf "\n--: "uData_cuda_"$operation"_"$iteration"_"$shuffle"_"$padding doesn't exist!\n"
	printf "********************************\n"
	exit
fi

if [ ! -e "uData_maple_"$operation ];then
	printf "********************************"
	printf "\n--: "uData_maple_"$operation doesn't exist!\n"
	printf "********************************\n"
	exit
fi

# printf "********************************\n"
printf "\n\t==:Result of Verification:\n\t"
diff -q "uData_cuda_"$operation"_"$iteration"_"$shuffle"_"$padding "uData_maple_"$operation
diff -s "uData_cuda_"$operation"_"$iteration"_"$shuffle"_"$padding "uData_maple_"$operation
printf "\n****************************************************************\n"
