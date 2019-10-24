# Tobias Scheithauer, August 2018

# This script can be used for reporting the results of the Wochenende pipeline

import os, sys, time
from Bio import SeqIO, SeqUtils
import pysam
import numpy as np
import pandas as pd
import click

def solid_normalization(gc):
	# TODO: get normalization model
	return 1

# Command Line Argument Parsing
@click.command()
# Slow mode uses the bam file while standard mode uses the bam.text created by wochenende pipeline. Slow mode gives an advanced output.
@click.option('--slow', default=False, help='Use this flag if you want to process the bam file instead of *.bam.txt')
@click.option('--input_file', help='The output file of Wochenende. Either *.bam.txt or .bam (use with --slow flag only)')
@click.option('--refseq_file', help='The refseq file used by Wochenende.')
# While illumina does not, SOLID sequencing data requires special normalization. Therefore an extra step is required. The model has to be defined in the function above.
@click.option('--sequencer', help='Sequencer technology used (solid or illumina)')
@click.option('--sample_name', help='Name of the sample. Used for output file naming.')
def reporting(slow, input_file, refseq_file, sequencer, sample_name):
	if slow:
		click.echo('Started slow mode.') 
		click.echo(f'Using {input_file} as alignment file')
		click.echo(f'Using {refseq_file} as refseq file')
		click.echo()
		if sequencer == 'illumina':
			# slow illumina reporting
			click.echo('starting illumina reporting')
			# creating lists for dataframe creation
			species_list, chr_length_list, read_count_list, basecount_list, gc_ref_list, gc_reads_list = [], [], [], [], [], []
			for seq_record in SeqIO.parse(refseq_file, 'fasta'):
				species_list.append(seq_record.name)
				chr_length_list.append(len(seq_record.seq))
				read_count_list.append(pysam.AlignmentFile(input_file, 'rb').count(contig=seq_record.name))
				# joining all reads to get number of bases in experiment and gc content of reads
				joined_reads = ''.join([read.query_sequence for read in pysam.AlignmentFile(input_file, 'rb').fetch(contig=seq_record.name)])
				basecount_list.append(sum([len(joined_reads)]))
				gc_ref_list.append(SeqUtils.GC(seq_record.seq))
				gc_reads_list.append(SeqUtils.GC(joined_reads))
			res_df = pd.DataFrame(data={
				'species': species_list, 
				'chr_length': chr_length_list, 
				'gc_ref': gc_ref_list, 
				'gc_reads': gc_reads_list, 
				'read_count': read_count_list, 
				'basecount': basecount_list,
			})
			res_df['reads_per_million_ref_bases'] = res_df['read_count']/(res_df['chr_length']/1000000)
			res_df['reads_per_million_reads_in_experiment'] = res_df['read_count'] / (res_df['read_count'].sum()/1000000)
			# calculating bacteria per human cell
			human_refs = ['1','2','3','4','5','6','7','8','9','10', '11','12','13','14','15','16','17','18','19','20','21','22','x','y','mt']
			human_cov = res_df[res_df['species'].isin(human_refs)]['basecount'].sum()/res_df[res_df['species'].isin(human_refs)]['chr_length'].sum()
			print(human_cov)
			#bam.txtres_df['bacteria_per_human_cell'] = (res_df['basecount']/res_df['chr_length']) / human_cov
			# total normalization RPMM
			res_df['RPMM'] = res_df['read_count'] / (res_df['chr_length']/1000000 * res_df['read_count'].sum()/1000000)
			res_df.to_csv(f'{sample_name}.reporting.unsorted.csv', sep='\t', float_format='%.1f', index=False)
			res_df_filtered_and_sorted = res_df.loc[res_df['read_count'] >= 20].sort_values(by='RPMM', ascending=False)
			res_df_filtered_and_sorted.to_csv(f'{sample_name}.reporting.sorted.csv', sep='\t', float_format='%.1f', index=False)
		elif sequencer == 'solid':
			# slow solid reporting (works like illumina reporting but width normalization)
			click.echo('starting solid reporting')	
			species_list, chr_length_list, read_count_list, basecount_list, gc_ref_list, gc_reads_list = [], [], [], [], [], []
			for seq_record in SeqIO.parse(refseq_file, 'fasta'):
				species_list.append(seq_record.name)
				chr_length_list.append(len(seq_record.seq))
				read_count_list.append(pysam.AlignmentFile(input_file, 'rb').count(contig=seq_record.name))
				joined_reads = ''.join([read.query_sequence for read in pysam.AlignmentFile(input_file, 'rb').fetch(contig=seq_record.name)])
				basecount_list.append(sum([len(joined_reads)]))
				gc_ref_list.append(SeqUtils.GC(seq_record.seq))
				gc_reads_list.append(SeqUtils.GC(joined_reads))
			res_df = pd.DataFrame(data={
				'species': species_list, 
				'chr_length': chr_length_list, 
				'gc_ref': gc_ref_list, 
				'gc_reads': gc_reads_list, 
				'read_count': read_count_list, 
				'basecount': basecount_list,
			})
			res_df['reads_per_million_ref_bases'] = res_df['read_count']/(res_df['chr_length']/1000000)
			res_df['reads_per_million_reads_in_experiment'] = res_df['read_count'] / (res_df['read_count'].sum()/1000000)
			# special SOLID normalization steps
			res_df['norm_factor'] = [solid_normalization(gc) for gc in res_df['gc_ref']]
			res_df['read_count'] = res_df['read_count'] * res_df['norm_factor']
			res_df['basecount'] = res_df['basecount'] * res_df['norm_factor']
			res_df['reads_per_million_ref_bases'] = res_df['reads_per_million_ref_bases'] * res_df['norm_factor']
			res_df['reads_per_million_reads_in_experiment'] = res_df['reads_per_million_reads_in_experiment'] * res_df['norm_factor']
			human_refs = ['1','2','3','4','5','6','7','8','9','10', '11','12','13','14','15','16','17','18','19','20','21','22','x','y','mt']
			human_cov = res_df[res_df['species'].isin(human_refs)]['basecount'].sum()/res_df[res_df['species'].isin(human_refs)]['chr_length'].sum()
			print(human_cov)
			res_df['bacteria_per_human_cell'] =  (res_df['ibasecount']/res_df['chr_length']) / human_cov
			res_df['norm_factor'] = None
			# total normalization RPMM
			res_df['RPMM'] = res_df['read_count'] / (res_df['chr_length']/1000000 * res_df['read_count'].sum()/1000000)
			res_df.to_csv(f'{sample_name}.reporting.unsorted.csv', sep='\t', float_format='%.1f', index=False)
			res_df_filtered_and_sorted = res_df.loc[res_df['read_count'] >= 20].sort_values(by='RPMM', ascending=False)
			res_df_filtered_and_sorted.to_csv(f'{sample_name}.reporting.sorted.csv', sep='\t', float_format='%.1f', index=False)
		else: 
			click.echo('please specify sequencing technology')
			sys.exit(1)
	else:
		click.echo(f'Using {input_file} as alignment file')
		click.echo(f'Using {refseq_file} as refseq file')
		click.echo()
		if sequencer == 'illumina':
			# standard illumina reporting
			click.echo('starting illumina reporting')
			# reading in wochenende output file without last line (* as species name)
			res_df = pd.read_csv(input_file, sep='\t', header=None, names=['species', 'chr_length', 'read_count'], usecols=[0,1,2])[:-1]
			# get gc content of ref sequences
			gc_ref_dict = {}
			for seq_record in SeqIO.parse(refseq_file, 'fasta'):
				gc_ref_dict[seq_record.name] = SeqUtils.GC(str(seq_record.seq).replace('N', ''))
			res_df['gc_ref'] = [gc_ref_dict[s] for s in res_df['species']]
			res_df['reads_per_million_ref_bases'] = res_df['read_count']/(res_df['chr_length']/1000000)
			res_df['reads_per_million_reads_in_experiment'] = res_df['read_count'] / (res_df['read_count'].sum()/1000000)
			# total normalization RPMM
			res_df['RPMM'] = res_df['read_count'] / (res_df['chr_length']/1000000 * res_df['read_count'].sum()/1000000)

			#calculating bacteria per human cell
			#check for human_refs to be correct!
			human_refs = ['1_1_1_1','1_1_1_2','1_1_1_3','1_1_1_4','1_1_1_5','1_1_1_6','1_1_1_7','1_1_1_8','1_1_1_9','1_1_1_10', '1_1_1_11','1_1_1_12','1_1_1_13','1_1_1_14','1_1_1_15',\
					'1_1_1_16','1_1_1_17','1_1_1_18','1_1_1_19','1_1_1_20','1_1_1_21','1_1_1_22','1_1_1_X','1_1_1_Y']
			human_cov = res_df[res_df['species'].isin(human_refs)]['read_count'].sum()
			res_df['bacteria_per_human_cell'] = (6191.39 * res_df['reads_per_million_ref_bases']) / human_cov

			#rounding to 2 decimals, except for bacteria_per_human_cell, which gets 4 decimals
			cols = ['gc_ref', 'reads_per_million_ref_bases', 'reads_per_million_reads_in_experiment', 'RPMM']
			res_df[cols] = res_df[cols].round(2)
			res_df['bacteria_per_human_cell'] = res_df['bacteria_per_human_cell'].round(4)

			res_df.to_csv(f'{sample_name}.reporting.unsorted.csv', sep='\t', index=False)
			res_df_filtered_and_sorted = res_df.loc[res_df['read_count'] >= 20].sort_values(by='RPMM', ascending=False)
			res_df_filtered_and_sorted.to_csv(f'{sample_name}.reporting.sorted.csv', sep='\t', index=False)
		elif sequencer == 'solid':
			# standard solid reporting (works like illumina reporting but width normalization)
			click.echo('starting solid reporting')	
			res_df = pd.read_csv(input_file, sep='\t', header=None, names=['species', 'chr_length', 'read_count'], usecols=[0,1,2])[:-1]
			gc_ref_list = []
			for seq_record in SeqIO.parse(refseq_file, 'fasta'):
				gc_ref_dict[seq_record.name] = SeqUtils.GC(str(seq_record.seq).replace('N', ''))
			res_df['gc_ref'] = gc_ref_list
			res_df['reads_per_million_ref_bases'] = res_df['read_count']/(res_df['chr_length']/1000000)
			res_df['reads_per_million_reads_in_experiment'] = res_df['read_count'] / (res_df['read_count'].sum()/1000000)
			# special SOLID normalization
			res_df['norm_factor'] = [solid_normalization(gc) for gc in res_df['gc_ref']]
			res_df['read_count'] = res_df['read_count'] * res_df['norm_factor']
			res_df['basecount'] = res_df['basecount'] * res_df['norm_factor']
			res_df['reads_per_million_ref_bases'] = res_df['reads_per_million_ref_bases'] * res_df['norm_factor']
			res_df['reads_per_million_reads_in_experiment'] = res_df['reads_per_million_reads_in_experiment'] * res_df['norm_factor']
			res_df['norm_factor'] = None
			# total normalization RPMM
			res_df['RPMM'] = res_df['read_count'] / (res_df['chr_length']/1000000 * res_df['read_count'].sum()/1000000)
			
			#calculating bacteria per human cell
			#check for human_refs to be correct!
			human_refs = ['1_1_1_1','1_1_1_2','1_1_1_3','1_1_1_4','1_1_1_5','1_1_1_6','1_1_1_7','1_1_1_8','1_1_1_9','1_1_1_10', '1_1_1_11','1_1_1_12','1_1_1_13','1_1_1_14','1_1_1_15',\
					'1_1_1_16','1_1_1_17','1_1_1_18','1_1_1_19','1_1_1_20','1_1_1_21','1_1_1_22','1_1_1_X','1_1_1_Y']
			human_cov = res_df[res_df['species'].isin(human_refs)]['read_count'].sum()
			res_df['bacteria_per_human_cell'] = (6191.39 * res_df['reads_per_million_ref_bases']) / human_cov

			#rounding to 2 decimals, except for bacteria_per_human_cell, which gets 4 decimals
			cols = ['gc_ref', 'reads_per_million_ref_bases', 'reads_per_million_reads_in_experiment', 'RPMM']
			res_df[cols] = res_df[cols].round(2)
			res_df['bacteria_per_human_cell'] = res_df['bacteria_per_human_cell'].round(4)

			res_df.to_csv(f'{sample_name}.reporting.unsorted.csv', sep='\t', index=False)
			res_df_filtered_and_sorted = res_df.loc[res_df['read_count'] >= 20].sort_values(by='RPMM', ascending=False)
			res_df_filtered_and_sorted.to_csv(f'{sample_name}.reporting.sorted.csv', sep='\t', index=False)
		else: 
			click.echo('please specify sequencing technology')
			sys.exit(1)
 
	
if __name__ == '__main__':
	reporting()
