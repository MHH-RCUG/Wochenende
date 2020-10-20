#!/bin/bash
## Supply the Wochenende bam.txt input as arg1, bash run_Wochenende_reporting_SLURM.sh in.bam.txt

# set partition
#SBATCH -p normal

# set run on x MB node only
#SBATCH --mem 20000

# set run on bigmem node only
#SBATCH --cpus-per-task 4

# share node
#SBATCH --share

# set max wallclock time
# SBATCH --time=4:00:00

# set name of job
#SBATCH --job-name=Wo_report

# mail alert at start, end and abortion of execution
#SBATCH --mail-type=ALL
#SBATCH --mail-user=<mailaddress>

echo "Input bam: " $1
bamtxt=$1

# Source miniconda3
. /mnt/ngsnfs/tools/miniconda3/etc/profile.d/conda.sh

# Activate env on cluster node
conda activate wochenende

#2020
#ref="/lager2/rcug/seqres/metagenref/bwa/2020_09_massiveref.fa"
ref="/lager2/rcug/seqres/metagenref/bwa/refSeqs_allKingdoms_2020_03.fa"
#2019
#ref_jan19="/ngsssd1/tuem_mp/refseqs/refseq2019_Jan/e_all_refseqs_Jan2019/all_kingdoms_refseq_2019_Jan_final.fasta"
#ref_test19="/ngsssd1/tuem_mp/refseqs/refseq2019_Jan/e_all_refseqs_Jan2019/reference_file_without_PseudomonasAeruginosa/all_kingdoms_refseq_2019_Jan_final_withoutPA.fasta"
#ref="/lager2/rcug/seqres/metagenref/bwa/all_kingdoms_refseq_2019_Jan_final.fasta"
#ref="/lager2/rcug/seqres/metagenref/bwa/refSeqs_allKingdoms_201910_3.fasta"
#ref="/lager2/rcug/seqres/metagenref/bwa/all_kingdoms_refseq_2019_Jan_final_mm10_no_human.fasta"
#ref="/lager2/rcug/seqres/metagenref/bwa/mm10_plus_ASF_OMM.fasta"
#ref="/lager2/rcug/seqres/metagenref/bwa/mm10_plus_ASF.fasta"
#ref="/lager2/rcug/seqres/metagenref/bwa/mm10_plus_OMM.fasta"
#pre2019
#ref2016="/lager2/rcug/seqres/metagenref/bwa/2016_06_PPKC_metagenome_test_1p_spec_change_cln.fa"
#ref_nov18="/ngsssd1/tuem_mp/refseq2018/f_all_refgenomes_2018/all_refgenomes_2018.fasta"
#ref="/working2/tuem/metagen/refs/2016/bwa/2016_06_PPKC_metagenome_test_1p_genus.fa"
#ref="/lager2/rcug/seqres/metagenref/bwa/2016_06_PPKC_metagenome_test_1p_spec_change_cln.fa"
#ref="/working2/tuem/metagen/refs/2016/bwa/2016_06_PPKC_metagenome_test_1p_spec_change.fa"

#genomes
#ref="/lager2/rcug/seqres/HS/bwa/hg19.fa"
#ref="/lager2/rcug/seqres/HS/bwa/GRCh37.fa"
#ref="/lager2/rcug/seqres/HS/bwa/Homo_sapiens.GRCh38.dna.toplevel.fa"
#ref="/lager2/rcug/seqres/HS/bwa/GRCh38_no_alt.fa"
#ref="/lager2/rcug/seqres/HS/bwa/Homo_sapiens.GRCh38.dna.chromosome.MT.fa"
#ref="/lager2/rcug/seqres/MM/bwa/mm10.fa"
#ref="/lager2/rcug/seqres/RN/bwa/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa"
#ref="/lager2/rcug/seqres/DR/bwa/GRCz10.fa"
#ref="/lager2/rcug/seqres/SS/bwa/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"
#ref="/lager2/rcug/seqres/PA/bwa/NC_008463.fna"


# Run script
srun python3 basic_reporting.py --input_file $bamtxt --reference $ref --sequencer illumina --output_name $bamtxt

#conda activate wochenende
#python3 basic_reporting.py --input_file tmp_R1.ndp.lc.trm.s.mq30.01mm.dup.bam.txt --refseq_file /lager2/rcug/seqres/metagenref/2016_06_PPKC_metagenome_test_1p_spec_change_cln.fa --sequencer illumina --sample_name test

