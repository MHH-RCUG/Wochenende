#!/usr/bin/env python
# coding: utf-8

# raspir
# Contact details: pust.marie-madlen@mh-hannover.de
# Last updated: 06 August 2021

# import
from __future__ import print_function

import os
import sys
import argparse
import logging
import itertools
import random
import warnings
from itertools import count, takewhile

warnings.simplefilter(action='ignore', category=Warning)

import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from scipy import fftpack, stats
from scipy.stats import linregress

# matplotlib init
matplotlib.use('Agg')

# logging
logging.basicConfig(format='%(asctime)s - %(message)s', level=logging.DEBUG)
logging.getLogger('matplotlib.font_manager').disabled = True

# argparse#
class CustomFormatter(argparse.ArgumentDefaultsHelpFormatter,
                      argparse.RawDescriptionHelpFormatter):
    pass

desc = 'raspir: rare species identifier'
epi = """DESCRIPTION:

Output:
	1. Optional, graph of spectrum vs frequency per cycle ( $PREFIX_[organism]_freq.png )
	2. Table (.CSV format) of final results (see raspir repo README) ( $PREFIX_final_stats.csv )
"""
parser = argparse.ArgumentParser(description=desc, epilog=epi,
                                 formatter_class=CustomFormatter)
argparse.ArgumentDefaultsHelpFormatter
parser.add_argument('csv_file', metavar='input_file', type=str,
                    help='input file (.CSV format)')
parser.add_argument('out_prefix', metavar='output_prefix', type=str,
                    help='output file (prefix only)')
parser.add_argument('-d', '--display_images', action='store_true',
                     help='display frequency plots per taxa')
parser.add_argument('-e', '--set_error', type=float, default=0.01,
                    help='std-error cutoff parameter')
parser.add_argument('-a', '--set_alpha', type=float, default=0.05,
                    help='alpha parameter')
parser.add_argument('--version', action='version', version='1.0.2')


# Global parameters
norm_cpm = 1000000

# functions
def frange(start, stop, step):
    """
    Define range of the reference's read distance
    """
    return takewhile(lambda x: x < stop, count(start, step))

def read_count(x, min_reads=4):
    """
    Read count table
    """
    items_diff = [abs(j - i) for i, j in zip(x['Position'], x['Position'][1:])]
    a = [0]
    items_diff = a + items_diff
    x['items_diff'] = items_diff
    x['gap'] = np.where(x['items_diff'] == 1, 0, 1)
    x['readCount'] = np.cumsum(x['gap'])
    x = x.drop_duplicates(subset='readCount', keep='first')
    x = x.drop(['items_diff', 'gap'], 1)
    
    # remove all Organisms with less than `min_read` reads
    x['indexNames'] = x['readCount'].sum()
    store_index = x[x['indexNames'] < int(min_reads)].index
    x.drop(store_index, inplace=True)
    return x
 
def normalise_genome_position(x):
    """
    Normalise position (circular genome)
    """
    x['PositionNorm0'] = np.where(x['Position'] > (x['GenomeLength'] / 2),
                                  (x['GenomeLength'] - x['Position']),
                                  x['Position'])
    x['PositionNorm'] = x['PositionNorm0']**(1/2)
    
    # Reference position
    n_reads = x['readCount'].max()
    start_position_ref = int(1)
    end_position_ref = x['GenomeLength'].iloc[0]
    end_position_ref = end_position_ref + n_reads
    increase_by = (end_position_ref / n_reads)
    x['ref_Position'] = list(frange(start_position_ref, end_position_ref,
                                    increase_by))
    x['ref_Position'] = x['ref_Position'].astype(int)
    x['PositionNorm_ref0'] = np.where(x['ref_Position'] > (x['GenomeLength'] / 2),
                                      (x['GenomeLength'] - x['ref_Position']),
                                      x['ref_Position'])
    x['PositionNorm_ref'] = x['PositionNorm_ref0'].astype(int)
    return x


def make_time_domain(x):
    """
    Time domain signal
    """
    # check if read count of organism matches with minimum requirement
    mean_depth = int(x['Depth'].mean())
    n_reads = int(x['readCount'].max())
    species_name = x['Organism'].iloc[0]
    x['PositionNorm'] = x['PositionNorm'] * x['Depth']
    x['PositionNorm_ref'] = x['PositionNorm_ref'] * mean_depth

    reference_combinations_distances_sort = []
    real_combinations_distances_sort = []

    # calculate the biological distance
    real_read_positions = sorted(x['PositionNorm'])
    reference_read_positions = sorted(x['PositionNorm_ref'])

    if n_reads > int(1000):
        random.seed(222)
        real_select_random = random.sample(real_read_positions, 400)
        real_read_combinations_sub = list(itertools.combinations(real_select_random, 2))
        real_combinations_distances_sub = [abs(i - j) for i, j in real_read_combinations_sub]
        real_combinations_distances_sort_sub = sorted(real_combinations_distances_sub)
        real_combinations_distances_sort.append(real_combinations_distances_sort_sub)

        # calculate the reference distance
        reference_select_random = random.sample(reference_read_positions, 400)
        reference_read_combinations_sub = list(itertools.combinations(reference_select_random, 2))
        reference_combinations_distances_sub = [abs(i - j) for i, j in reference_read_combinations_sub]
        reference_combinations_distances_sort_sub = sorted(reference_combinations_distances_sub)
        reference_combinations_distances_sort.append(reference_combinations_distances_sort_sub)
    else:
        real_read_combinations = list(itertools.combinations(real_read_positions, 2))
        real_combinations_distances = [abs(i - j) for i, j in real_read_combinations]
        real_combinations_distances_sort1 = sorted(real_combinations_distances)
        real_combinations_distances_sort.append(real_combinations_distances_sort1)

        # calculate the reference distance
        reference_read_combinations = list(itertools.combinations(reference_read_positions, 2))
        reference_combinations_distances = [abs(i - j) for i, j in reference_read_combinations]
        reference_combinations_distances_sort1 = sorted(reference_combinations_distances)
        reference_combinations_distances_sort.append(reference_combinations_distances_sort1)

    # create output data frame
    df = pd.DataFrame(list(zip(reference_combinations_distances_sort, real_combinations_distances_sort)),
                      columns=['Reference', 'Real'])
    df['Organism'] = species_name
    df2 = df.apply(lambda i: i.explode() if i.name in ['Reference', 'Real'] else i)
    return df2


# frequency domain signal (fds)
def fourier_trans(x):
    species_name = x['Organism'].iloc[0]
    sep = '_'
    # check for separators and use try except to avoid errors
    if sep in species_name:
        try:
            stripped_name = species_name.split(sep)
            stripped_name2 = stripped_name[3] + ' ' + stripped_name[4]
        except:
            logging.info('Warning  Name could not be parsed correctly using "_" splits: {}', str(species_name))
            stripped_name = species_name
            stripped_name2 = species_name
            
    else:
        stripped_name = species_name
        stripped_name2 = species_name

    x['fft_ref1'] = np.fft.fft(x['Reference'])
    x['fft_bio1'] = np.fft.fft(x['Real'])

    x['fft_ref'] = [complex(np.around(items2.real), np.around(items2.imag)) for items2 in x['fft_ref1']]
    x['fft_bio'] = [complex(np.around(items2.real), np.around(items2.imag)) for items2 in x['fft_bio1']]

    x['fft_abs_ref'] = np.abs(x['fft_ref'])
    x['fft_abs_bio'] = np.abs(x['fft_bio'])
    x['fft_abs_ref_sqrt'] = np.around(x['fft_abs_ref'] / norm_cpm, 2)
    x['fft_abs_bio_sqrt'] = np.around(x['fft_abs_bio'] / norm_cpm, 2)

    # Pearson correlation
    if (sum(x['fft_abs_ref_sqrt']) > 0) & (sum(x['fft_abs_bio_sqrt']) > 0):
        pearson_corr = linregress(x['fft_abs_ref_sqrt'], x['fft_abs_bio_sqrt'])
        pearson_standard_error0 = pearson_corr[4]
        pearson_corr_r0, pearson_corr_p0 = stats.pearsonr(x['fft_abs_ref_sqrt'], x['fft_abs_bio_sqrt'])
        pearson_corr_r = round(pearson_corr_r0, 4)
        pearson_corr_p = round(pearson_corr_p0, 10)
        euclidean_dist_0 = np.linalg.norm(x['fft_abs_ref_sqrt']-x['fft_abs_bio_sqrt'])
        euclidean_dist = round(euclidean_dist_0, 1)
        pearson_standard_error = round(pearson_standard_error0, 5)
        return species_name, pearson_corr_r, pearson_corr_p, pearson_standard_error, euclidean_dist
    else:
        pearson_corr_r2 = 0
        pearson_corr_p2 = 0
        pearson_standard_error2 = 1
        euclidean_dist2 = 1
        return species_name, pearson_corr_r2, pearson_corr_p2, pearson_standard_error2, euclidean_dist2


def make_freq_images(x, set_images):
    logging.basicConfig(level=logging.ERROR)
    if set_images is False:
        pass
    elif set_images is True:
        species_name = x['Organism'].iloc[0]
        path_real = x['PathName'].iloc[0]
		
        x['fft_ref1'] = np.fft.fft(x['Reference'])
        x['fft_bio1'] = np.fft.fft(x['Real'])

        x['fft_ref'] = [complex(np.around(items2.real), np.around(items2.imag)) for items2 in x['fft_ref1']]
        x['fft_bio'] = [complex(np.around(items2.real), np.around(items2.imag)) for items2 in x['fft_bio1']]

        x['fft_abs_ref'] = np.abs(x['fft_ref'])
        x['fft_abs_bio'] = np.abs(x['fft_bio'])
        x['fft_abs_ref_sqrt'] = np.around(x['fft_abs_ref'] / norm_cpm, 2)
        x['fft_abs_bio_sqrt'] = np.around(x['fft_abs_bio'] / norm_cpm, 2)

        # plot frequency signal
        val_bio = x['Real']
        val_ref = x['Reference']
        x_bio0 = x['fft_abs_bio_sqrt']
        x_reference0 = x['fft_abs_ref_sqrt']

        x_bio1 = x_bio0.sort_values()
        x_bio2 = pd.concat([x_bio1[::2], x_bio1[len(x_bio1)-2:0:-2]])
        x_bio = x_bio2.tolist()
        x_reference1 = x_reference0.sort_values()
        x_reference2 = pd.concat([x_reference1[::2], x_reference1[len(x_reference1)-2:0:-2]])
        x_reference = x_reference2.tolist()
        freqs_bio0 = fftpack.fftfreq(len(val_bio))
        freqs_bio = sorted(freqs_bio0)

        freqs_ref0 = fftpack.fftfreq(len(val_ref))
        freqs_ref = sorted(freqs_ref0)
        a = (len(freqs_ref) - len(x_reference))
        b = (len(freqs_bio) - len(x_bio))
        x_reference += [0]*a
        x_bio += [0]*b
        x_reference3 = np.sqrt(x_reference)
        sep = '_'
        # add error handling in case separator not present for some taxa. eg chrY etc from mouse
        if sep in species_name:
            try:
                stripped_name = species_name.split(sep)
                stripped_name2 = stripped_name[3] + ' ' + stripped_name[4]
                stripped_name3 = stripped_name[3] + '_' + stripped_name[4]
            except:
                logging.info('Warning  Name could not be parsed correctly using "_" splits: {}', str(species_name))
                # don't try to strip the names for chromosomes without separators
                # we don't need to attribute stripped_name2 since it is not used elsewhere
                stripped_name = species_name
                stripped_name2 = species_name
                stripped_name3 = species_name
        else:
            logging.info('Warning  Name could not be parsed correctly using "_" splits: {}', str(species_name))
            stripped_name = species_name
            stripped_name2 = species_name
            stripped_name3 = species_name

        fig, ax1 = plt.subplots(1, 1, figsize=(2.5, 2))
        fig.suptitle(stripped_name2, style='italic', fontsize=4)
        ax1.plot(freqs_ref, x_reference3, "black", linewidth=1, linestyle='--', label="Reference", alpha=0.4)
        ax1.plot(freqs_bio, x_bio, 'blue', linewidth=1, linestyle='--', label="Sample")
        ax1.legend(framealpha=1, loc='upper right', fontsize=3)
        ax1.fill_between(freqs_ref, x_reference3, x_bio, facecolor='pink', alpha=0.2, interpolate=True)
        fig.text(0.5, 0.025, "Frequency per cycle",  ha='center', va='center', fontsize=4)
        fig.text(0.010, 0.5, "Spectrum", ha='center', va='center', rotation='vertical', fontsize=4)
        plt.xticks(fontsize=3)
        plt.yticks(fontsize=3)
        outfile = '_'.join([path_real, stripped_name3, 'freq.png'])
        plt.savefig(outfile, dpi=600)
        logging.basicConfig(level=logging.DEBUG)
        logging.info('  File written: {}'.format(outfile))
        plt.close()

def final_table(x, set_error, set_alpha):
    """
    Create file data table
    """
    a0 = pd.DataFrame(x, columns=['Pearson'])
    a = a0.dropna(axis=0, how='all')
    a = a[['Species', 'r_value', 'p_value', 'stError', 'euclideanR0']] = pd.DataFrame(a['Pearson'].tolist())
    a.columns = ['a', 'b', 'c', 'd', 'e',
                 'Species', 'r_value', 'p_value', 'stError', 'euclideanR0']
    a.reset_index(drop=True, inplace=True)
    a = a.drop(a.columns[[0, 1]], axis=1)
    a = a[['Species', 'r_value', 'p_value', 'stError', 'euclideanR0']]
    a['euclidean'] = np.around((1 / a['euclideanR0']) * 1000, 3)
    a['distribution'] = np.where(
        (a['p_value'] < set_alpha) & (a['r_value'] > 0.5) & (a['stError'] < set_error) & (a['euclidean'] < 0.5),
        'uniform', 'nonuniform')
    b1 = a[['Species', 'r_value', 'p_value', 'stError', 'euclidean', 'distribution']]
    b2 = b1[b1.distribution == 'uniform']
    return b2


def process_csv(file_name, out_prefix, args):
    with open(file_name, newline='') as inF:
        df = pd.read_csv(inF, delimiter=',')

        # filtering reads attrib to human chromosomes
        pattern_del = '1_1_1_'
        filter_approach = df['Organism'].str.contains(pattern_del, na=False)
        df = df[~filter_approach]
        logging.info('1a) Human reads have been removed')
        
        # filtering  reads attributed to organisms starting with chr eg some human, mouse chrs
        pattern_del = 'chr'
        filter_approach = df['Organism'].str.startswith(pattern_del, na=False)
        df = df[~filter_approach]
        logging.info('1b) Human chr reads have been removed')

        # counting reads per organism
        df = df.dropna(subset=['GenomeLength'])
        df_filter1 = df.groupby('Organism').apply(read_count)
        logging.info('2) Continue with the first position of each read')

        if df_filter1.empty is True:
            logging.warning('Note: Dataset is empty')
            outfile = '{}_filtered.csv'.format(out_prefix)
            df_filter1.to_csv(outfile, index=False)
            logging.info('  File written: {}'.format(outfile))
            return None
        else:
            df_filter = df_filter1.reset_index(drop=True)
            df_filter2 = df_filter.groupby('Organism').apply(normalise_genome_position)
            logging.info('3) Genome position normalised')

            position_domain0 = df_filter2.reset_index(drop=True)
            position_domain = position_domain0.groupby("Organism").apply(make_time_domain)
            logging.info('4) Time-domain signal built')

            frequency_domain0 = position_domain.reset_index(drop=True)
            frequency_domain0['PathName'] = os.path.join(out_prefix)
            frequency_domain = frequency_domain0.groupby('Organism').apply(fourier_trans)
            logging.info('5) Frequency-domain signal generated')

            frequency_domain0.groupby('Organism').apply(make_freq_images, set_images=args.display_images)
            logging.info('6) Frequency plots produced (optional)')

            stat_table = final_table(frequency_domain,
                                     set_error=args.set_error,
                                     set_alpha=args.set_alpha)
            logging.info('7) Output table has been generated')

            outfile = '{}_final_stats.csv'.format(out_prefix)
            stat_table.to_csv(outfile, index=False)
            logging.info('  File written: {}'.format(outfile))
            logging.info('8) Run successful')

def main(args):
    """
    Main interface
    """
    # output directory
    outdir = os.path.split(args.out_prefix)[0]
    if outdir != '' and not os.path.isdir(outdir):
        os.makedirs(args.outdir)
    # processing each file
    process_csv(args.csv_file, out_prefix=args.out_prefix, args=args)

if __name__ == "__main__":
    args = parser.parse_args()
    main(args)
