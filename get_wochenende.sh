# Dr Colin Davenport
# April 2. 2020
# A helper script to copy Wochenende scripts/directory to the current path.
# Uses settings from the config yaml file

path_we=/mnt/ngsnfs/tools/dev/Wochenende
source scripts/parse_yaml.sh
eval $(parse_yaml config_cln.yml)

path_we=$wochenende_dir

cp $path_we/README* .
cp $path_we/*.sh .
cp $path_we/*.py .
cp -R $path_we/extract/ .
cp -R $path_we/plots/ .
cp -R $path_we/reporting/ .
