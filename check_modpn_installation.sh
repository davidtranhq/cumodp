
printf "\n"
check=${MODPN_HOME:-default}
# printf "$check \n"
mv modpn ..
(cd ../modpn; ./install_modpn.sh;)

if [ "$check" == "default" ]; then
	printf "The variable for MODPN_HOME is NOT defined !\n"
	export MOPDN_IS_INSTALLED=1
else
	printf "The variable for MODPN_HOME is set as follows:\n"
	printf "\t $MODPN_HOME\n"
	export MOPDN_IS_INSTALLED=1
fi;
printf "=========================================\n"