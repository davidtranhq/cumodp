#!/bin/bash

#four argument, operation, dataSet mode+xSeed+ySeed
PROJECT_NAME="BigPrimeFieldFFT_3"
PROJECT="$CUMODP_HOME/new_tests/$PROJECT_NAME"
INPUT_DATA_DIR=$PROJECT/test/tmp

#==================================================
if [ ! $# -eq 7 ]; then
	printf "\n--: There should be 6 arguments:"
	printf "\n--: [0]operation,\t[1]iteration,\t\t[2]shuffle,\t\t[3]padding"
	printf "\n--: [4]dataSet,\t\t[5]inputVectorSize,\t[6]outputFile\n"

	printf "\n--: operation should be one of following:"
	printf "\n--: 0: Addition,\t 1: Multiplication,\t 2: Subtraction,\t 3: CyclicShift"
	printf "\n--: 4: FFT_16,\t\t 5: FFT_256,\t\t 6: FFT_4096,\t\t 7:FFT_16K \n"
	printf "*******************\n"
	exit
fi
#==================================================
operation=$1
iteration=$2
shuffle=$3
padding=$4
dataSet=$5
inputVectorSize=$6
outputFile=$7
#==================================================
outputFile=$PROJECT/test/timingReports/$outputFile
dataPath=$PROJECT/test/data/$dataSet
#==================================================
if [ ! -d $dataPath ];then
	printf "*******************\n"
	printf "\n--: Dataset = $dataSet doesn't exist!\n"
	printf "*******************\n"
	exit
fi
#==================================================
printf "\nRunning CUDA for \t operation = $operation,\t iteration = $iteration,\t shuffle = $shuffle,\t padding = $padding\n" >tmp
printf "=============================\n" >> tmp
#==================================================
cp $dataPath/xData $INPUT_DATA_DIR
cp $dataPath/yData $INPUT_DATA_DIR

# printf "++COPY: xData ---> $PROJECT_NAME/test/tmp/ \n\n"
# printf "++COPY: yData ---> $PROJECT_NAME/test/tmp/ \n\n"
# printf "********************************\n\n"

cd $PROJECT/bin/
./bigArithmetic_test.bin $operation $iteration $shuffle $padding $dataSet $inputVectorSize >> tmp
# printf "$operation $iteration $shuffle $padding $dataSet $inputVectorSize\n" >> tmp
printf "=========================================================\n" >> tmp
# printf "*******************\n"
# mv $PROJECT/test/tmp/uData_cuda $dataPath/"uData_cuda_"$operation"_"$iteration"_"$shuffle"_"$padding
# mv $PROJECT/test/tmp/uData_cuda $dataPath/"uData_cuda_"$operation"_"$iteration"_"$shuffle"_"$padding

# mv $PROJECT/bin/uData_cuda $dataPath/"uData_cuda_"$operation"_"$iteration"_"$shuffle"_"$padding
# mv $PROJECT/bin/xData_cuda $dataPath/"xData_cuda_"$operation"_"$iteration"_"$shuffle"_"$padding

# rm $PROJECT/test/tmp/*
# printf "\nResult:$dataSet/"uData_cuda_"$operation"_"$iteration"_"$shuffle"_"$padding\n\n"
cat tmp >> $outputFile
cat tmp
cat /dev/null > tmp
