#!/bin/bash

#four argument, operation, dataSet mode+xSeed+ySeed
# printf "+++++++++++++++++++       RUN CUDA        +++++++++++++++++++\n"
PROJECT_NAME="BigPrimeFieldFFT_3"
PROJECT="$CUMODP_HOME/new_tests/$PROJECT_NAME"

if [ ! $# -eq 6 ]; then
	printf "\n--: There should be 6 arguments: operation, iteration, shuffle, padding, dataSet, inputVectorSize\n"
	printf "\n--: 0: Addition,\t 1: Multiplication,\t 2: Subtraction,\t 3: CyclicShift"
	printf "\n--: 4: FFT_16,\t\t 5: FFT_256,\t\t 6: FFT_4096,\t\t 7:FFT_16K \n"
	printf "*******************\n"
	exit
fi

operation=$1
iteration=$2
shuffle=$3
padding=$4
dataSet=$5
inputVectorSize=$6


dataPath=$PROJECT/test/data/$dataSet

if [ ! -d $dataPath ];then
	printf "*******************\n"
	printf "\n--: Dataset = $dataSet doesn't exist!\n"
	printf "*******************\n"
	exit
fi

printf "\nRunning CUDA for \t operation = $operation,\t iteration = $iteration,\t shuffle = $shuffle,\t padding = $padding\n\n"

cp $dataPath/xData $PROJECT/test/tmp/xData
cp $dataPath/yData $PROJECT/test/tmp/yData

# printf "++COPY: xData ---> $PROJECT_NAME/test/tmp/ \n\n"
# printf "++COPY: yData ---> $PROJECT_NAME/test/tmp/ \n\n"
# printf "********************************\n\n"

cd $PROJECT/bin/

./bigArithmetic_test.bin $operation $iteration $shuffle $padding $dataSet $inputVectorSize
# printf "*******************\n"
# mv $PROJECT/test/tmp/uData_cuda $dataPath/"uData_cuda_"$operation"_"$iteration"_"$shuffle"_"$padding
# mv $PROJECT/test/tmp/uData_cuda $dataPath/"uData_cuda_"$operation"_"$iteration"_"$shuffle"_"$padding

mv $PROJECT/bin/uData_cuda $dataPath/"uData_cuda_"$operation"_"$iteration"_"$shuffle"_"$padding
mv $PROJECT/bin/xData_cuda $dataPath/"xData_cuda_"$operation"_"$iteration"_"$shuffle"_"$padding

rm $PROJECT/test/tmp/*

printf "\nResult:$dataSet/"uData_cuda_"$operation"_"$iteration"_"$shuffle"_"$padding\n\n"
printf "*********************************************************\n"
# printf "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n"
# printf "+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++\n\n"
