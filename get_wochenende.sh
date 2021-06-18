# Dr Colin Davenport
# June 2021
# A helper script to copy Wochenende scripts/directory to the current path.
# Uses settings from the config yaml file

echo "INFO: setup.sh should have created  Wochenende env variable on this server and user account"
echo "INFO: If you get errors here check setup.sh was setup and run properly"

#path_we=/mnt/ngsnfs/tools/dev/Wochenende
#path_we=/mnt/ngsnfs/tools/dev/we_config2/wochenende
path_we=$WOCHENENDE_DIR

cp $path_we/README* .
cp $path_we/*.sh .
cp $path_we/*.py .
cp -R $path_we/extract/ .
cp -R $path_we/plots/ .
cp -R $path_we/reporting/ .
