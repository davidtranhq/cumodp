#!/bin/bash

#four argument, operation, data mode+xSeed+ySeed
# printf "+++++++++++++++++++       RUN MAPLE        +++++++++++++++++++\n"
# PROJECT_NAME="BigPrimeFieldFFT_3"
# PROJECT="$CUMODP_HOME/new_tests/$PROJECT_NAME"

. ./projectProperties.conf
operation=$1
dataSet=$2
verbosity=$3

if [ ! $# -eq 2 ]; then
	printf "\n--: There should be 2 arguments: operation, dataSet"
	printf "\n--: 0: Addition,\t 1: Multiplication,\t 2: Subtraction,\t 3: CyclicShift"
	printf "\n--: 4: FFT_16,\t\t 5: FFT_256,\t\t 6: FFT_4096,\t\t 7:FFT_16K \n"
	printf "*******************\n"
	exit
fi

if [ $operation -gt 7 ]; then
	printf "operation $operation is not defined!\n"
	printf "\n*******************\n"
	exit
fi

dataPath=$PROJECT/test/data/$dataSet

if [ ! -d $dataPath ];then
	printf "*******************\n"
	printf "\n--: Dataset = $dataSet doesn't exist!\n"
	printf "*******************\n"
	exit
fi

cp $dataPath/xData $PROJECT/test/tmp/xData
cp $dataPath/yData $PROJECT/test/tmp/yData

cd $PROJECT/scripts/mapleScripts/

printf "\n++ Running Maple for \t operation = $operation, \t dataSet = $dataSet\n"

maple maple_operation_$operation.mpl
 # > /dev/null

mv $PROJECT/test/tmp/uData_maple_$operation $dataPath/

rm $PROJECT/test/tmp/*
printf "++ result: $dataSet/uData_maple_$operation \n\n"
printf "*********************************************************\n"