# Dr Colin Davenport
# April 2. 2020
# A helper script to copy Wochenende scripts/directory to the current path.

# change as appropriate
path_we=/mnt/ngsnfs/tools/dev/Wochenende

cp $path_we/README* .
cp $path_we/*.sh .
cp $path_we/*.py .
cp -R $path_we/extract/ .
cp -R $path_we/plots/ .
cp -R $path_we/reporting/ .
