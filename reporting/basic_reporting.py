# Tobias Scheithauer, August 2018

# This script can be used for reporting the results of the Wochenende pipeline

import os, sys, time
from Bio import SeqIO, SeqUtils
import pysam
import numpy as np
import pandas as pd
import click

def solid_normalization(gc):
	# TODO: get model from Lutz
	return 1

@click.command()
@click.option('--alignment_file', help='The bam output file of Wochenende.')
@click.option('--refseq_file', help='The refseq file used by Wochenende.')
@click.option('--sequencer', help='Sequencer technology used (solid or illumina)')
@click.option('--sample_name', help='Name of the sample. Used for output file naming.')
def reporting(alignment_file, refseq_file, sequencer, sample_name):
	click.echo(f'Using {alignment_file} as alignment file')
	click.echo(f'Using {refseq_file} as refseq file')
	click.echo()
	
	if sequencer == 'illumina':
		# illumina reporting
		click.echo('starting illumina reporting')	
		species_list, size_list, readcount_list, basecount_list, gc_ref_list, gc_reads_list = [], [], [], [], [], []
		for seq_record in seqio.parse(refseq_file, 'fasta'):
			species_list.append(seq_record.name)
			size_list.append(len(seq_record.seq))
			readcount_list.append(pysam.alignmentfile(alignment_file, 'rb').count(contig=seq_record.name))
			joined_reads = ''.join([read.query_sequence for read in pysam.alignmentfile(alignment_file, 'rb').fetch(contig=seq_record.name)])
			basecount_list.append(sum([len(joined_reads)]))
			gc_ref_list.append(sequtils.gc(seq_record.seq))
			gc_reads_list.append(sequtils.gc(joined_reads))
		res_df = pd.dataframe(data={
			'species': species_list, 
			'size': size_list, 
			'gc_ref': gc_ref_list, 
			'gc_reads': gc_reads_list, 
			'readcount': readcount_list, 
			'basecount': basecount_list,
		})
		res_df['reads_per_million_ref_bases'] = res_df['readcount']/(res_df['size']/1000000)
		res_df['reads_per_million_reads_in_experiment'] = res_df['readcount'] / (res_df['readcount'].sum()/1000000)
		human_refs = ['1','2','3','4','5','6','7','8','9','10', '11','12','13','14','15','16','17','18','19','20','21','22','x','y','mt']
		human_cov = res_df[res_df['species'].isin(human_refs)]['basecount'].sum()/res_df[res_df['species'].isin(human_refs)]['size'].sum()
		print(human_cov)
		res_df['bacteria_per_human_cell'] =  (res_df['basecount']/res_df['size']) / human_cov
		res_df.to_csv(f'{sample_name}.basic_reporting.csv', sep='\t')

	elif sequencer == 'solid':
		# solid reporting
		click.echo('starting solid reporting')	
		species_list, size_list, readcount_list, basecount_list, gc_ref_list, gc_reads_list = [], [], [], [], [], []
		for seq_record in seqio.parse(refseq_file, 'fasta'):
			species_list.append(seq_record.name)
			size_list.append(len(seq_record.seq))
			readcount_list.append(pysam.alignmentfile(alignment_file, 'rb').count(contig=seq_record.name))
			joined_reads = ''.join([read.query_sequence for read in pysam.alignmentfile(alignment_file, 'rb').fetch(contig=seq_record.name)])
			basecount_list.append(sum([len(joined_reads)]))
			gc_ref_list.append(sequtils.gc(seq_record.seq))
			gc_reads_list.append(sequtils.gc(joined_reads))
		res_df = pd.dataframe(data={
			'species': species_list, 
			'size': size_list, 
			'gc_ref': gc_ref_list, 
			'gc_reads': gc_reads_list, 
			'readcount': readcount_list, 
			'basecount': basecount_list,
		})
		res_df['reads_per_million_ref_bases'] = res_df['readcount']/(res_df['size']/1000000)
		res_df['reads_per_million_reads_in_experiment'] = res_df['readcount'] / (res_df['readcount'].sum()/1000000)
		res_df['norm_factor'] = [solid_normalization(gc) for gc in res_df['gc_ref']]
		res_df['readcount'] = res_df['readcount'] * res_df['norm_factor']
		res_df['basecount'] = res_df['basecount'] * res_df['norm_factor']
		res_df['reads_per_million_ref_bases'] = res_df['reads_per_million_ref_bases'] * res_df['norm_factor']
		res_df['reads_per_million_reads_in_experiment'] = res_df['reads_per_million_reads_in_experiment'] * res_df['norm_factor']
		human_refs = ['1','2','3','4','5','6','7','8','9','10', '11','12','13','14','15','16','17','18','19','20','21','22','x','y','mt']
		human_cov = res_df[res_df['species'].isin(human_refs)]['basecount'].sum()/res_df[res_df['species'].isin(human_refs)]['size'].sum()
		print(human_cov)
		res_df['bacteria_per_human_cell'] =  (res_df['ibasecount']/res_df['size']) / human_cov
		res_df['norm_factor'] = None
		res_df.to_csv(f'{sample_name}.basic_reporting.csv', sep='\t')

	else: 
		click.echo('please specify sequencing technology')
		sys.exit(1)
	
if __name__ == '__main__':
	reporting()
