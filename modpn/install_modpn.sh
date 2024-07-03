#!/bin/bash

##################################################################
# Set path for MODPN
##################################################################
install_dir=$(pwd)
modpn_variables_path=$install_dir"/modpn_variables"

##################################################################
# Ask for installation path
##################################################################

printf "Current installation path=\n"
printf "\t $install_dir \n"

read -p "Do you like to install modpn/MODPN in a different path? [y/n] (default=n): " choice
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

cat /dev/null > $modpn_variables_path

printf "export MODPN_HOME=$install_dir \n" >> $modpn_variables_path
printf "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:\$MODPN_HOME/src \n" >> $modpn_variables_path

# check if environmental variables are already set in .bashrc
if grep -Fxq "    source $modpn_variables_path " ~/.bashrc 
then
	# printf "Path is already set! \n"
	printf ""
else
	printf "if [ -e $modpn_variables_path ]; then \n" >> ~/.bashrc
	printf "    source $modpn_variables_path \n" >> ~/.bashrc
	printf "fi;\n">> ~/.bashrc
fi

printf "=================================================\n"
printf "\nPath for MODPN_HOME is set to the following:\n"
printf "\t $install_dir\n"
printf "\n"
printf "You can modify the environmental variables in the following:\n"
printf "\t$modpn_variables_path \n"
printf "=================================================\n"

##################################################################
# Install the modpn
##################################################################
source ~/.bashrc
source $modpn_variables_path
echo $modpn_HOME
read -p "Last step, press any key to compile modpn library ... "
cd $modpn_HOME
make

##################################################################
