# Test Wochenende pipeline locally without SLURM
# First edit the adapter locations to relevant adapter FASTA files in run_Wochenende.py

conda activate wochenende

python3 run_Wochenende.py --metagenome testdb --testWochenende testdb/reads_R1.fastq --force_restart
