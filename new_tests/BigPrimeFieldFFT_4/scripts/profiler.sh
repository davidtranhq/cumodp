#!/bin/bash

#four argument, operation, dataSet mode+xSeed+ySeed
# printf "+++++++++++++++++++       PROFILING        +++++++++++++++++++\n"
# PROJECT_NAME="BigPrimeFieldFFT_3"
# PROJECT="$CUMODP_HOME/new_tests/$PROJECT_NAME"

####################
source projectProperties
####################
if [ ! $# -eq 6 ]; then
	printf "\n--: There should be 6 arguments: operation, iteration, shuffle, padding, dataSet, inputVectorSize\n"
	printf "\n--: 0: Addition,\t 1: Multiplication,\t 2: Subtraction,\t 3: CyclicShift"
	printf "\n--: 4: FFT_16,\t\t 5: FFT_256,\t\t 6: FFT_4096,\t\t 7:FFT_16K \n"
	printf "*********************************************************\n"
	exit
fi

operation=$1
iteration=$2
shuffle=$3
padding=$4
dataSet=$5
inputVectorSize=$6

outputFile=$PROJECT"/test/profilingReports/report_"$operation"_"$iteration"_"$shuffle"_"$padding"_"$dataSet".log"

dataPath=$PROJECT/test/data/$dataSet

if [ ! -d $dataPath ];then
	# printf "*******************\n"
	printf "\nDataset = $dataSet doesn't exist!\n\n"
	printf "**************************************\n"
	exit
fi

cd $PROJECT/bin/

printf "\nNVprof for operation = $operation, iteration = $iteration, shuffle = $shuffle, padding = $padding, dataSet = $dataSet, inputVectorSize = $inputVectorSize \n\n"

printf "\nNVprof for operation = $operation, iteration = $iteration, shuffle = $shuffle, padding = $padding, dataSet = $dataSet, inputVectorSize = $inputVectorSize \n\n" > header

# metrics="achieved_occupancy,shared_replay_overhead,local_replay_overhead,dram_read_throughput,dram_write_throughput,shared_efficiency,global_replay_overhead,global_cache_replay_overhead,dram_utilization,branch_efficiency,l1_shared_utilization,l1_cache_local_hit_rate,l1_cache_global_hit_rate,ipc,inst_replay_overhead"

metrics="achieved_occupancy,shared_replay_overhead,local_replay_overhead,dram_read_throughput,dram_write_throughput,shared_efficiency,global_cache_replay_overhead,dram_utilization,branch_efficiency,l1_shared_utilization,l1_cache_local_hit_rate,l1_cache_global_hit_rate,ipc,inst_replay_overhead"

events="elapsed_cycles_sm,l1_local_load_hit,l1_local_load_miss,l1_local_store_hit,l1_local_store_miss,local_load,local_store,shared_load,shared_store"

# nvprof --metrics $metrics --events $events --log-file $outputFile ./bigArithmetic_test.bin $operation $iteration $shuffle $padding $dataSet $inputVectorSize

# nvprof --metrics $metrics --log-file $outputFile ./bigArithmetic_test.bin $operation $iteration $shuffle $padding $dataSet $inputVectorSize > /dev/null

# nvprof  --log-file $outputFile ./bigArithmetic_test.bin $operation $iteration $shuffle $padding $dataSet $inputVectorSize > /dev/null

nvprof --metrics $metrics --log-file withMetrics ./bigArithmetic_test.bin $operation $iteration $shuffle $padding $dataSet $inputVectorSize > /dev/null

nvprof  --log-file withTiming ./bigArithmetic_test.bin $operation $iteration $shuffle $padding $dataSet $inputVectorSize > /dev/null

cat withMetrics withTiming > $outputFile
# cat header $outputFile > tmp
# mv tmp $outputFile && rm header
# outputFile="${outputFile/$PROJECT/}"
printf "*********************************************************\n"
