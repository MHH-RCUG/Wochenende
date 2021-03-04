# Wochenende config file. Set environment variables for Wochenende https://github.com/MHH-RCUG/Wochenende
# v0.1 2021_03
# Only uncommented lines without a # are read

# Users: change this as appropriate
# Install_directory for Wochenende on your system, full path
wochenende_install_dir=/mnt/ngsnfs/tools/dev/Wochenende
# Install_directory for Haybaler on your system, full path
haybaler_install_dir=/mnt/ngsnfs/tools/dev/haybaler



# Users: please do not change this.
# add to user .bashrc so they are available in users environment after new login
echo "export WOCHENENDE_DIR=$wochenende_install_dir" >> ~/.bashrc
echo "export HAYBALER_DIR=$haybaler_install_dir" >> ~/.bashrc

