PROJECT_NAME:="BigPrimeFieldFFT_4"
PROJECT:=$(CUMODP_HOME)/new_tests/$(PROJECT_NAME)
# PROJECT:=$(CUMODP_HOME)/tests/$(PROJECT_NAME)

# include scripts/projectProperties.conf
# include scripts/operationList.conf

xSeed:=11
ySeed:=22

#operation list
opAdd:=0
opMultiplication:=1
opSubtraction:=2
opCyclicShift:=3
opFFT_16:=4
opFFT_256:=5
opFFT_4k:=6
opFFT_16k:=7
opFFT_64k:=14
opFFT_2:=9
opNewMult:=10
opFFT_general:=14

#number of iterations
nIterations:=1

#shuffling input data or not
noShuffling:=0
doShuffling:=1

#padding
padding_0:=0
padding_1:=1
padding_2:=2
padding_3:=3
padding_4:=4

#vectorSize = number of coefficients = 16, 256, 1k, 4k, 16k, 64k, 1M
vector_16:=$$((16))
vector_256:=$$((256))
vector_1k:=$$((1*1024))
vector_4k:=$$((4*1024))
vector_16k:=$$((16*1024))
vector_64k:=$$((64*1024))
vector_128k:=$$((128*1024))
vector_256k:=$$((256*1024))
vector_512k:=$$((512*1024))
vector_1m:=$$((1024*1024))
vector_4m:=$$((4*1024*1024))
vector_8m:=$$((8*1024*1024))
vector_16m:=$$((16*1024*1024))	

currentVectorSize:=$(vector_1m)
currentVectorSuffix:="1m"
dataSet_0:="data_0_$(xSeed)_$(ySeed)_$(currentVectorSuffix)"
dataSet_1:="data_1_$(xSeed)_$(ySeed)_$(currentVectorSuffix)"
dataSet_2:="data_2_$(xSeed)_$(ySeed)_$(currentVectorSuffix)"
dataSet_3:="data_3_$(xSeed)_$(ySeed)_$(currentVectorSuffix)"
GPU_NAME:=$(shell scripts/getGpuName.sh)
HOST_NAME:=$(GPU_NAME)"_"$(shell hostname)
############################
all: compile generateData measureTimeAddition
	@ printf "making all bigArithmeticFFT_3\n"
############################
test:
	@ printf "testing bigArithmeticFFT_3\n"
############################
check:
	@ printf "checking bigArithmeticFFT_3\n"
############################
getGpuName:
	@# @printf $(GPU_NAME)
	@printf "$(HOST_NAME)\n"
############################
switch:
	@ cd $(PROJECT)/scripts && sudo ./switch-cuda-7.0.sh on && printf "Please run: 'source ~/.bashrc'\n\n"
	@ cd $(nvcc -V)
############################
enablePM:
	@ cd $(PROJECT)/scripts && ./enablePM.sh
############################
compile: 
	@ cd $(PROJECT)/scripts && ./makeLinks.sh
	@ cd $(PROJECT)/scripts && ./compile.sh
	# @ cd $(PROJECT)/src/FFT_2 && make
	@ cd $(PROJECT)/scripts && ./precompute_pow_omega.sh
############################
benchmarkSmallFFT:
	@ cd $(PROJECT)/src/FFT/ && make benchmark
############################
generateData:
	@ #generate 1k datasets for verification purposes
	@ #generate random data (mode=0)
	@ cd scripts && ./generateData.sh 0 $(xSeed) $(ySeed) $(vector_1k) 

	@ #generate data for kernel case (mode=1)
	@ cd scripts && ./generateData.sh 1 $(xSeed) $(ySeed) $(vector_1k) 

	@ #generate all elements = seed (mode=2)
	@ cd scripts && ./generateData.sh 2 $(xSeed) $(ySeed) $(vector_1k)

	@ #generate all elements = seed (mode=3)
	@ cd scripts && ./generateData.sh 3 $(xSeed) $(ySeed) $(vector_1k)

	@ #generate random data (mode=0)
	@ cd scripts && ./generateData.sh 0 $(xSeed) $(ySeed) $(currentVectorSize)

	@ #generate data for kernel case (mode=1)
	@ cd scripts && ./generateData.sh 1 $(xSeed) $(ySeed) $(currentVectorSize)

	@ #generate all elements = seed (mode=2)
	@ cd scripts && ./generateData.sh 2 $(xSeed) $(ySeed) $(currentVectorSize)

	@ #generate all elements = seed (mode=3)
	@ cd scripts && ./generateData.sh 3 $(xSeed) $(ySeed) $(currentVectorSize)
########################################################
# Maple Verifying Scripts
########################################################
verifyAdditionMaple:
	@ #operation, dataSet ("data_Mode_xSeed_ySeed")
	@ cd scripts && ./runMaple.sh $(opAdd) "data_0_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./runMaple.sh $(opAdd) "data_1_$(xSeed)_$(ySeed)_1k" 
	@ cd scripts && ./runMaple.sh $(opAdd) "data_2_$(xSeed)_$(ySeed)_1k" 
############################
verifyAdditionCuda:
	@ #operation, iteration, shuffle, padding, dataSet ("data_Mode_xSeed_ySeed")
	@ cd scripts && ./runCuda.sh $(opAdd) $(nIterations) $(noShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./runCuda.sh $(opAdd) $(nIterations) $(noShuffling) $(padding_0) "data_1_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./runCuda.sh $(opAdd) $(nIterations) $(noShuffling) $(padding_0) "data_2_$(xSeed)_$(ySeed)_1k" $(vector_1k)

	@ cd scripts && ./compare.sh $(opAdd) $(nIterations) $(noShuffling) $(padding_0)  "data_0_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./compare.sh $(opAdd) $(nIterations) $(noShuffling) $(padding_0)  "data_1_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./compare.sh $(opAdd) $(nIterations) $(noShuffling) $(padding_0)  "data_2_$(xSeed)_$(ySeed)_1k"	
############################
verifyAdditionShuffling:
	@ #operation, iteration, shuffle, padding, dataSet ("data_Mode_xSeed_ySeed")
	@ cd scripts && ./runCuda.sh $(opAdd) $(nIterations) $(doShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k) > /dev/null
	@ cd scripts && ./runCuda.sh $(opAdd) $(nIterations) $(doShuffling) $(padding_0) "data_1_$(xSeed)_$(ySeed)_1k" $(vector_1k) > /dev/null
	@ cd scripts && ./runCuda.sh $(opAdd) $(nIterations) $(doShuffling) $(padding_0) "data_2_$(xSeed)_$(ySeed)_1k" $(vector_1k) > /dev/null

	@ cd scripts && ./compare.sh $(opAdd) $(nIterations) $(doShuffling) $(padding_0)  "data_0_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./compare.sh $(opAdd) $(nIterations) $(doShuffling) $(padding_0)  "data_1_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./compare.sh $(opAdd) $(nIterations) $(doShuffling) $(padding_0)  "data_2_$(xSeed)_$(ySeed)_1k"	
############################
verifyMultiplicationMaple:
	@ #operation, dataSet ("data_Mode_xSeed_ySeed")
	@ cd scripts && ./runMaple.sh $(opMultiplication) "data_0_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./runMaple.sh $(opMultiplication) "data_1_$(xSeed)_$(ySeed)_1k" 
	@ cd scripts && ./runMaple.sh $(opMultiplication) "data_2_$(xSeed)_$(ySeed)_1k" 
############################
verifyMultiplicationCuda:
	@ #operation, iteration, shuffle, padding, dataSet ("data_Mode_xSeed_ySeed")
	@ cd scripts && ./runCuda.sh $(opMultiplication) $(nIterations) $(noShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./runCuda.sh $(opMultiplication) $(nIterations) $(noShuffling) $(padding_0) "data_1_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./runCuda.sh $(opMultiplication) $(nIterations) $(noShuffling) $(padding_0) "data_2_$(xSeed)_$(ySeed)_1k" $(vector_1k)

	@ cd scripts && ./compare.sh $(opMultiplication) $(nIterations) $(noShuffling) $(padding_0)  "data_0_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./compare.sh $(opMultiplication) $(nIterations) $(noShuffling) $(padding_0)  "data_1_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./compare.sh $(opMultiplication) $(nIterations) $(noShuffling) $(padding_0)  "data_2_$(xSeed)_$(ySeed)_1k"
############################
verifySubtractionMaple:
		@ #operation, dataSet ("data_Mode_xSeed_ySeed")
	@ cd scripts && ./runMaple.sh $(opSubtraction) "data_0_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./runMaple.sh $(opSubtraction) "data_1_$(xSeed)_$(ySeed)_1k" 
	@ cd scripts && ./runMaple.sh $(opSubtraction) "data_2_$(xSeed)_$(ySeed)_1k" 

############################
verifySubtractionCuda:
	@ #operation, iteration, shuffle, padding, dataSet ("data_Mode_xSeed_ySeed")
	@ cd scripts && ./runCuda.sh $(opSubtraction) $(nIterations) $(noShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./runCuda.sh $(opSubtraction) $(nIterations) $(noShuffling) $(padding_0) "data_1_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./runCuda.sh $(opSubtraction) $(nIterations) $(noShuffling) $(padding_0) "data_2_$(xSeed)_$(ySeed)_1k" $(vector_1k)

	@ cd scripts && ./compare.sh $(opSubtraction) $(nIterations) $(noShuffling) $(padding_0)  "data_0_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./compare.sh $(opSubtraction) $(nIterations) $(noShuffling) $(padding_0)  "data_1_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./compare.sh $(opSubtraction) $(nIterations) $(noShuffling) $(padding_0)  "data_2_$(xSeed)_$(ySeed)_1k"

########################################################
# CUDA Verifying Scripts
########################################################
verifyCyclicShiftMaple:
	@ #operation, dataSet ("data_Mode_xSeed_ySeed")
	@ cd scripts && ./runMaple.sh $(opCyclicShift) "data_0_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./runMaple.sh $(opCyclicShift) "data_1_$(xSeed)_$(ySeed)_1k" 
	@ cd scripts && ./runMaple.sh $(opCyclicShift) "data_2_$(xSeed)_$(ySeed)_1k" 

############################
verifyCyclicShiftCuda:
	@ #operation, iteration, shuffle, padding, dataSet ("data_Mode_xSeed_ySeed")
	@ cd scripts && ./runCuda.sh $(opCyclicShift) $(nIterations) $(noShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./runCuda.sh $(opCyclicShift) $(nIterations) $(noShuffling) $(padding_0) "data_1_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./runCuda.sh $(opCyclicShift) $(nIterations) $(noShuffling) $(padding_0) "data_2_$(xSeed)_$(ySeed)_1k" $(vector_1k)

	@ cd scripts && ./compare.sh $(opCyclicShift) $(nIterations) $(noShuffling) $(padding_0)  "data_0_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./compare.sh $(opCyclicShift) $(nIterations) $(noShuffling) $(padding_0)  "data_1_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./compare.sh $(opCyclicShift) $(nIterations) $(noShuffling) $(padding_0)  "data_2_$(xSeed)_$(ySeed)_1k"

############################
verifyFFT_16Maple:
	@ #operation, dataSet ("data_Mode_xSeed_ySeed")
	@ cd scripts && ./runMaple.sh $(opFFT_16) "data_0_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./runMaple.sh $(opFFT_16) "data_1_$(xSeed)_$(ySeed)_1k" 
	@ cd scripts && ./runMaple.sh $(opFFT_16) "data_2_$(xSeed)_$(ySeed)_1k" 

############################
verifyFFT_16Cuda:
	@ #operation, iteration, shuffle, padding, dataSet ("data_Mode_xSeed_ySeed")
	@ cd scripts && ./runCuda.sh $(opFFT_16) $(nIterations) $(noShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./runCuda.sh $(opFFT_16) $(nIterations) $(noShuffling) $(padding_0) "data_1_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./runCuda.sh $(opFFT_16) $(nIterations) $(noShuffling) $(padding_0) "data_2_$(xSeed)_$(ySeed)_1k" $(vector_1k)

	@ cd scripts && ./runCuda.sh $(opFFT_16) $(nIterations) $(doShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./runCuda.sh $(opFFT_16) $(nIterations) $(doShuffling) $(padding_0) "data_1_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./runCuda.sh $(opFFT_16) $(nIterations) $(doShuffling) $(padding_0) "data_2_$(xSeed)_$(ySeed)_1k" $(vector_1k)

	@cd scripts && ./compareToPlain.sh $(opFFT_16) $(nIterations) $(doShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./compareToPlain.sh $(opFFT_16) $(nIterations) $(doShuffling) $(padding_0) "data_1_$(xSeed)_$(ySeed)_1k"
	@cd scripts && ./compareToPlain.sh $(opFFT_16) $(nIterations) $(doShuffling) $(padding_0) "data_2_$(xSeed)_$(ySeed)_1k"

############################
verifyFFT_256Maple:
	@ #operation, dataSet ("data_Mode_xSeed_ySeed")
	@ cd scripts && ./runMaple.sh $(opFFT_256) "data_0_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./runMaple.sh $(opFFT_256) "data_1_$(xSeed)_$(ySeed)_1k" 
	@ cd scripts && ./runMaple.sh $(opFFT_256) "data_2_$(xSeed)_$(ySeed)_1k" 

############################
verifyFFT_256Cuda:
	@ #operation, iteration, shuffle, padding, dataSet ("data_Mode_xSeed_ySeed")
	@ cd scripts && ./runCuda.sh $(opFFT_256) $(nIterations) $(noShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./runCuda.sh $(opFFT_256) $(nIterations) $(noShuffling) $(padding_0) "data_1_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./runCuda.sh $(opFFT_256) $(nIterations) $(noShuffling) $(padding_0) "data_2_$(xSeed)_$(ySeed)_1k" $(vector_1k)

	@ cd scripts && ./compare.sh $(opFFT_256) $(nIterations) $(noShuffling) $(padding_0)  "data_0_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./compare.sh $(opFFT_256) $(nIterations) $(noShuffling) $(padding_0)  "data_1_$(xSeed)_$(ySeed)_1k"
	@ cd scripts && ./compare.sh $(opFFT_256) $(nIterations) $(noShuffling) $(padding_0)  "data_2_$(xSeed)_$(ySeed)_1k"

############################
verifyFFT_4kMaple:
	@ #generate random data (mode=0)
	@ cd scripts && ./generateData.sh 0 $(xSeed) $(ySeed) $(vector_4k)
	@ #generate data for kernel case (mode=1)
	@ cd scripts && ./generateData.sh 1 $(xSeed) $(ySeed) $(vector_4k)
	@ #generate all elements = seed (mode=2)
	@ cd scripts && ./generateData.sh 2 $(xSeed) $(ySeed) $(vector_4k)
	@ #operation, dataSet ("data_Mode_xSeed_ySeed")
	# @ cd scripts && ./runMaple.sh $(opFFT_4k) "data_0_$(xSeed)_$(ySeed)_4k"
	@ cd scripts && ./runMaple.sh $(opFFT_4k) "data_1_$(xSeed)_$(ySeed)_4k" 
	@ cd scripts && ./runMaple.sh $(opFFT_4k) "data_2_$(xSeed)_$(ySeed)_4k" 

############################
verifyFFT_4kCuda:	
	@ #operation, iteration, shuffle, padding, dataSet ("data_Mode_xSeed_ySeed")
	@ cd scripts && ./runCuda.sh $(opFFT_4k) $(nIterations) $(noShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_4k" $(vector_4k)
	@ cd scripts && ./runCuda.sh $(opFFT_4k) $(nIterations) $(noShuffling) $(padding_0) "data_1_$(xSeed)_$(ySeed)_4k" $(vector_4k)
	@ cd scripts && ./runCuda.sh $(opFFT_4k) $(nIterations) $(noShuffling) $(padding_0) "data_2_$(xSeed)_$(ySeed)_4k" $(vector_4k)

	@ cd scripts && ./compare.sh $(opFFT_4k) $(nIterations) $(noShuffling) $(padding_0)  "data_0_$(xSeed)_$(ySeed)_4k"
	@ cd scripts && ./compare.sh $(opFFT_4k) $(nIterations) $(noShuffling) $(padding_0)  "data_1_$(xSeed)_$(ySeed)_4k"
	@ cd scripts && ./compare.sh $(opFFT_4k) $(nIterations) $(noShuffling) $(padding_0)  "data_2_$(xSeed)_$(ySeed)_4k"

########################################################
# Profiling Scripts
########################################################
measureTimeAddition:
	@ printf "\n"
	@ cd bin && ./bigArithmetic_test.bin $(opAdd) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	@ printf "\n=======================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opAdd) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	
	@ printf "\n======================================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opAdd) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	@ printf "\n=======================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opAdd) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)

	@ printf "\n======================================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opAdd) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)
	@ printf "\n=======================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opAdd) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	@ printf "\n======================================================\n"
########################################################
measureTimeMultiplication:
	@ printf "\n"
	@ cd bin && ./bigArithmetic_test.bin $(opNewMult) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	@ printf "\n=======================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opMultiplication) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	
	@ printf "\n======================================================\n"
		@ cd bin && ./bigArithmetic_test.bin $(opNewMult) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	@ printf "\n=======================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opMultiplication) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	
	@ printf "\n======================================================\n"

		@ cd bin && ./bigArithmetic_test.bin $(opNewMult) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)
	@ printf "\n=======================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opMultiplication) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)
	
	@ printf "\n======================================================\n"
########################################################
measureTimeCyclicShift:
	@ printf "\n"
	@ cd bin && ./bigArithmetic_test.bin $(opCyclicShift) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	@ printf "\n=======================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opCyclicShift) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	
	@ printf "\n======================================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opCyclicShift) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	@ printf "\n=======================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opCyclicShift) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)

	@ printf "\n======================================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opCyclicShift) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)
	@ printf "\n=======================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opCyclicShift) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	@ printf "\n======================================================\n"

########################################################
measureTimeFFT_16:
	@ rm -f test/timingReports/timing_FFT_16_$(HOST_NAME)

	@ cd scripts && ./measureTime.sh $(opFFT_16) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize) timing_FFT_16_$(HOST_NAME)
	@ cd scripts && ./measureTime.sh $(opFFT_16) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize) timing_FFT_16_$(HOST_NAME)

	@ cd scripts && ./measureTime.sh $(opFFT_16) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize) timing_FFT_16_$(HOST_NAME)
	@ cd scripts && ./measureTime.sh $(opFFT_16) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize) timing_FFT_16_$(HOST_NAME)

	@ cd scripts && ./measureTime.sh $(opFFT_16) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize) timing_FFT_16_$(HOST_NAME)
	@ cd scripts && ./measureTime.sh $(opFFT_16) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize) timing_FFT_16_$(HOST_NAME)

	@ less test/timingReports/timing_FFT_16_$(HOST_NAME)
########################################################
# measureTimeFFT64:
# 	@ printf "\n"
# 	@ cd bin && ./bigArithmetic_test.bin $(opFFT_64) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
# 	@ printf "\n=======================================\n"
# 	@ cd bin && ./bigArithmetic_test.bin $(opFFT_64) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	
# 	@ printf "\n======================================================\n"
# 	@ cd bin && ./bigArithmetic_test.bin $(opFFT_64) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
# 	@ printf "\n=======================================\n"
# 	@ cd bin && ./bigArithmetic_test.bin $(opFFT_64) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	
# 	@ printf "\n======================================================\n"
# 	@ cd bin && ./bigArithmetic_test.bin $(opFFT_64) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)
# 	@ printf "\n=======================================\n"
# 	@ cd bin && ./bigArithmetic_test.bin $(opFFT_64) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)
	
# 	@ printf "\n======================================================\n"
########################################################
measureTimeFFT_256:
	@ rm -f test/timingReports/timing_FFT_256_$(HOST_NAME)

	@ cd scripts && ./measureTime.sh $(opFFT_256) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize) timing_FFT_256_$(HOST_NAME)
	@ cd scripts && ./measureTime.sh $(opFFT_256) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize) timing_FFT_256_$(HOST_NAME)

	@ cd scripts && ./measureTime.sh $(opFFT_256) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize) timing_FFT_256_$(HOST_NAME)
	@ cd scripts && ./measureTime.sh $(opFFT_256) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize) timing_FFT_256_$(HOST_NAME)

	@ cd scripts && ./measureTime.sh $(opFFT_256) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize) timing_FFT_256_$(HOST_NAME)
	@ cd scripts && ./measureTime.sh $(opFFT_256) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize) timing_FFT_256_$(HOST_NAME)

	@ less test/timingReports/timing_FFT_256_$(HOST_NAME)

	########################################################
measureTimeFFT_64k:
	@ rm -f test/timingReports/timing_FFT_64k_$(HOST_NAME)

	@ cd scripts && ./measureTime.sh $(opFFT_64k) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize) timing_FFT_64k_$(HOST_NAME)
	@ cd scripts && ./measureTime.sh $(opFFT_64k) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize) timing_FFT_64k_$(HOST_NAME)

	@ cd scripts && ./measureTime.sh $(opFFT_64k) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize) timing_FFT_64k_$(HOST_NAME)
	@ cd scripts && ./measureTime.sh $(opFFT_64k) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize) timing_FFT_64k_$(HOST_NAME)

	@ cd scripts && ./measureTime.sh $(opFFT_64k) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize) timing_FFT_64k_$(HOST_NAME)
	@ cd scripts && ./measureTime.sh $(opFFT_64k) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize) timing_FFT_64k_$(HOST_NAME)

	@ less test/timingReports/timing_FFT_64k_$(HOST_NAME)
########################################################
measureTimeFFT_4k:
	@ printf "\n"
	@ cd bin && ./bigArithmetic_test.bin $(opFFT_4k) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	@ printf "\n=======================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opFFT_4k) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	
	@ printf "\n======================================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opFFT_4k) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	@ printf "\n=======================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opFFT_4k) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	
	@ printf "\n======================================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opFFT_4k) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)
	@ printf "\n=======================================\n"
	@ cd bin && ./bigArithmetic_test.bin $(opFFT_4k) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)
	
	@ printf "\n======================================================\n"
########################################################

profileAddition:
	@ cd ./	test/profilingReports/ && rm -f *.log
	@ cd scripts && ./profiler.sh $(opAdd) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opAdd) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opAdd) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	@ # @ cd scripts && ./profiler.sh $(opAdd) $(nIterations) $(doShuffling) $(padding_1) $(dataSet_0) $(currentVectorSize)
	@ # @ cd scripts && ./profiler.sh $(opAdd) $(nIterations) $(doShuffling) $(padding_2) $(dataSet_1) $(currentVectorSize)
	@ # @ cd scripts && ./profiler.sh $(opAdd) $(nIterations) $(doShuffling) $(padding_3) $(dataSet_2) $(currentVectorSize)

	@ cd scripts && ./profiler.sh $(opAdd) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opAdd) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opAdd) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	# @ cd scripts && ./profiler.sh $(opAdd) $(nIterations) $(noShuffling) $(padding_1) $(dataSet_0) $(currentVectorSize)
	# @ cd scripts && ./profiler.sh $(opAdd) $(nIterations) $(noShuffling) $(padding_1) $(dataSet_1) $(currentVectorSize)
	# @ cd scripts && ./profiler.sh $(opAdd) $(nIterations) $(noShuffling) $(padding_1) $(dataSet_2) $(currentVectorSize)

	@ cd scripts && ./generateReport.sh "profiling_report_addition"
	@ less test/profilingReports/profiling_report_addition
############################
profileMultiplication:
	@ cd ./	test/profilingReports/ && rm -f *.log
	#@ cd scripts && ./profiler.sh $(opNewMult) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	#@ cd scripts && ./profiler.sh $(opNewMult) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	#@ cd scripts && ./profiler.sh $(opNewMult) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	@ cd scripts && ./profiler.sh $(opMultiplication) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opMultiplication) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opMultiplication) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	@ cd scripts && ./profiler.sh $(opMultiplication) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opMultiplication) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opMultiplication) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	@ cd scripts && ./generateReport.sh "profiling_report_multiplication_$(HOST_NAME)"
	@ less test/profilingReports/profiling_report_multiplication_$(HOST_NAME)
############################
profileCyclicShift:
	@ cd ./	test/profilingReports/ && rm -f *.log
	@ cd scripts && ./profiler.sh $(opCyclicShift) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opCyclicShift) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opCyclicShift) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	@ cd scripts && ./profiler.sh $(opCyclicShift) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opCyclicShift) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opCyclicShift) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	@ cd scripts && ./generateReport.sh "profiling_report_cyclicShift"
	@ less test/profilingReports/profiling_report_cyclicShift
############################
profileFFT_16:
	@ cd ./	test/profilingReports/ && rm -f *.log
	@ cd scripts && ./profiler.sh $(opFFT_16) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opFFT_16) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opFFT_16) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	@ cd scripts && ./profiler.sh $(opFFT_16) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opFFT_16) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opFFT_16) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	@ cd scripts && ./generateReport.sh "profiling_report_FFT16_$(HOST_NAME)"
	@ less test/profilingReports/profiling_report_FFT16_$(HOST_NAME)
############################
# profileMultiplication:
# 	@ cd scripts && ./profiler.sh $(opMultiplication) 100 $(noShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
# 	# @ cd scripts && ./profiler.sh 1 100 0 1 "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
# 	# @ cd scripts && ./profiler.sh 1 100 0 2 "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
# 	@ cd scripts && ./generateReport.sh "profiling_report_multiplication"

############################
profileSubtraction:
	@ cd scripts && ./profiler.sh $(opSubtraction) 100 $(noShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	# @ cd scripts && ./profiler.sh 2 100 0 1 "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	# @ cd scripts && ./profiler.sh 2 100 0 2 "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	@ cd scripts && ./generateReport.sh "profiling_report_subtraction"

############################
# profileCyclicShift:
	# @ cd scripts && ./profiler.sh 3
	# @ cd scripts && ./profiler.sh $(opCyclicShift) 100 $(noShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
	# @ cd scripts && ./generateReport.sh "profiling_report_cyclicShift.log"

############################
# profileFFT_16:
# 	@ cd scripts && ./profiler.sh 2 100 0 0 "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
# 	@ cd scripts && ./profiler.sh $(opFFT_16) 100 $(noShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_1k" $(vector_1k)
# 	@ cd scripts && ./generateReport.sh "profiling_report_FFT_16"

# ############################
profileFFT_256:
	@ cd ./	test/profilingReports/ && rm -f *.log
	@ cd scripts && ./profiler.sh $(opFFT_256) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opFFT_256) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opFFT_256) $(nIterations) $(doShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	@ cd scripts && ./profiler.sh $(opFFT_256) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opFFT_256) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opFFT_256) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	@ cd scripts && ./generateReport.sh "profiling_report_FFT256_$(HOST_NAME)"
	@ less test/profilingReports/profiling_report_FFT256_$(HOST_NAME)
# ############################
# profileFFT_4k:
# 	@ cd scripts && ./profiler.sh $(opFFT_4k) 100 $(noShuffling) $(padding_0) "data_0_$(xSeed)_$(ySeed)_5k" $(vector_4k)
# 	@ cd scripts && ./generateReport.sh "profiling_report_FFT_4k	"
########################################################
profileFFT_64k:
	@ cd ./	test/profilingReports/ && rm -f *.log
	@ cd scripts && ./profiler.sh $(opFFT_64k) $(nIterations) $(doShuffling) 4 $(dataSet_0) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opFFT_64k) $(nIterations) $(doShuffling) 4 $(dataSet_1) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opFFT_64k) $(nIterations) $(doShuffling) 4 $(dataSet_2) $(currentVectorSize)

	#@ cd scripts && ./profiler.sh $(opFFT_64k) $(nIterations) $(noShuffling)  $(dataSet_0) $(currentVectorSize)
	#@ cd scripts && ./profiler.sh $(opFFT_64k) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	#@ cd scripts && ./profiler.sh $(opFFT_64k) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	@ cd scripts && ./generateReport.sh "profiling_report_FFT64k_$(HOST_NAME)"
	@ less test/profilingReports/profiling_report_FFT64k_$(HOST_NAME)
# ############################
profileFFT_1M:
	@ cd ./	test/profilingReports/ && rm -f *.log
	@ cd scripts && ./profiler.sh $(opFFT_general) $(nIterations) $(doShuffling) 5 $(dataSet_0) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opFFT_general) $(nIterations) $(doShuffling) 5 $(dataSet_1) $(currentVectorSize)
	@ cd scripts && ./profiler.sh $(opFFT_general) $(nIterations) $(doShuffling) 5 $(dataSet_2) $(currentVectorSize)

	#@ cd scripts && ./profiler.sh $(opFFT_64k) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_0) $(currentVectorSize)
	#@ cd scripts && ./profiler.sh $(opFFT_64k) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_1) $(currentVectorSize)
	#@ cd scripts && ./profiler.sh $(opFFT_64k) $(nIterations) $(noShuffling) $(padding_0) $(dataSet_2) $(currentVectorSize)

	@ cd scripts && ./generateReport.sh "profiling_report_FFT1m_$(HOST_NAME)"
	@ less test/profilingReports/profiling_report_FFT1m_$(HOST_NAME)
# ############################
# computeRootsOfUnity:
	# @cd scripts/roots_of_unity && maple compute_powers_omega.mpl
	# @cp scripts/roots_of_unity/powers_of_omega_K* src/
# ############################
# Clean
########################################################
clean:
	@ rm -rf $(PROJECT)/bin/*
	@ rm -f	 $(PROJECT)/test/compileMessage
	@ rm -rf $(PROJECT)/test/data/*
	@ rm -rf $(PROJECT)/test/tmp/*
	@ rm -rf $(PROJECT)/test/profilingReports/*
	@ rm -f $(PROJECT)/src/powers_of_omega_K*
############################
