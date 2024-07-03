#!/bin/bash

##################################################
# Check if MODPN is installed properly
##################################################
tabs 2
export MOPDN_IS_INSTALLED=0
source ./check_modpn_installation.sh
if [ $MOPDN_IS_INSTALLED -eq 0 ]; then
	printf "==================================================\n";
	printf "\tDependency error: MOPDN library is not installed!\n"
	printf "\tPlease retry after installing MODPN;\n"
	printf "\tYou can download the MODPN from the following link:\n"
	printf "\t\t http://cumodp.org/downloads.html\n"
	printf "==================================================\n\n";
	exit;
else
	printf "MODPN is installed and properly configured ... \n"
fi;

##################################################################
# Check if GPU card(s) and compiler are available and  configured
##################################################################
export GPU_PROPERLY_CONFIGURED=0
## the following script will set to GPU_PROPERLY_CONFIGURED=1 if 
## GPU and compiler is properly configured.
# source ./check_nvidia_gpu.sh
export GPU_PROPERLY_CONFIGURED=1

if [ $GPU_PROPERLY_CONFIGURED -eq 0 ]; then
	printf "\t Abort due to improper GPU configuration!\n";
	printf "====================================================\n\n"
	exit;
fi;

##################################################################
# Set path for CUMODP and MODPN
##################################################################
install_dir=$(pwd)
cumodp_variables_path=$install_dir"/cumodp_variables"

##################################################################
# Ask for installation path
##################################################################

printf "Current installation path=\n"
printf "\t $install_dir \n"

read -p "Do you like to install CUMODP in a different path? [y/n] (default=n): " choice
if [ "$choice" == "yes" ] || [ "$choice" == "y" ]; then
	is_inst_dir_set=0
	# while [ $is_inst_dir_set -eq 0 ]; do
		read -ep "Enter the new path: " new_install_dir
		install_dir=$new_install_dir
	# 	is_inst_dir_set=1;
	# done;
fi;

##################################################################
# Set up the installation path
##################################################################

mkdir -p $install_dir 
# && rm -r $install_dir/*
cd $install_dir

echo $install_dir
##################################################################
# Set environmental variables
##################################################################

cat /dev/null > $cumodp_variables_path

printf "export CUMODP_HOME=$install_dir \n" >> $cumodp_variables_path
printf "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$CUMODP_HOME/src \n" >> $cumodp_variables_path

# check if environmental variables are already set in .bashrc
if grep -Fxq "    source $cumodp_variables_path " ~/.bashrc 
then
	# printf "Path is already set! \n"
	printf ""
else
	printf "if [ -e $cumodp_variables_path ]; then \n" >> ~/.bashrc
	printf "    source $cumodp_variables_path \n" >> ~/.bashrc
	printf "fi;\n">> ~/.bashrc
fi

printf "=================================================\n"
printf "\nPath for CUMODP_HOME is set to the following:\n"
printf "\t $install_dir\n"
printf "\n"
printf "You can modify the environmental variables in the following:\n"
printf "\t$cumodp_variables_path \n"
printf "=================================================\n"

##################################################################
# Install the CUMODP
##################################################################
source ~/.bashrc
source $cumodp_variables_path
echo $CUMODP_HOME
read -p "Final step, press any key to compile CUMODP library ... "
cd $CUMODP_HOME
make

##################################################################
