# Colin Oct 2021
# Run raspir by Marie Pust

echo "INFO: link BAM files in"
bash batch_create_links.sh

echo "INFO: Start preparing the files for raspir"
bash run_SLURM_file_prep.sh

echo "INFO: Run raspir"
bash run_raspir_SLURM.sh

echo "INFO: Remove soft linked BAM files"
batch_remove_links.sh

echo "INFO: Raspir module completed"
