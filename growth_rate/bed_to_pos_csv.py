# Script to make pos.csv (inputfile for growth determiner) out of bed files
# Needs bed file as input. Associated bam.txt file needs to be in the same directory.
# Author: Sophia Poertner, 2021

import pandas as pd
import matplotlib.pylab as plt
import click
import random
import os
import sys


def read_bed(bed, path):
    df_input_file = pd.read_csv("{}/{}".format(path, bed), delimiter="\t", header=None, index_col=0)
    return df_input_file


def get_pos(df_input_file, filename, path, plot_samples, min_read_count):
    for organism in set(df_input_file.index.values):  # for every organism in the file
        organism_df = df_input_file.loc[[organism]]  # make a df containing only that organism
        if len(organism_df.index) >= min_read_count:
            print("INFO: working on organism {}".format(organism))
            positions = organism_df[1].tolist()
            genome_length = get_genome_length(filename, path, organism)
            normalised_position = norm_shuffle(positions, genome_length)
            save_as_csv(filename, path, normalised_position, organism)
            if plot_samples:
                plot_reads(normalised_position, filename, path, organism)
#    return positions


# normalise the read position so they are between 0 and 1. Shuffle the read positions
def norm_shuffle(position, genome_length):
    normalised_position = []
    for pos in position:
        normalised_position.append(pos / genome_length)
    random.shuffle(normalised_position)
    return normalised_position


# get genome length from bam.txt file
def get_genome_length(file, path, organism):
    if not os.path.isfile("{}/{}".format(path, file.replace("bed", "bam.txt").replace(".filt", ""))):
        print("Error: Could not find expected file: "+ file.replace("bed", "bam.txt").replace(".filt", "") + " at path: " + path)
        sys.exit("Error in bed_to_pos_csv.py: no bam.txt file found for file {}  - every bed file needs its associated bam.txt file!".format(file))
    else:
        bam_txt = pd.read_csv("{}/{}".format(path, file.replace("bed", "bam.txt").replace(".filt", "")), delimiter="\t", header=None,
                              index_col=0)
        genome_length = bam_txt[1][organism]
    return genome_length


def plot_reads(normalised_position, file, path, organism):
    round_list = [round(elem, 3) for elem in normalised_position]
    f = plt.figure()
    plt.hist(round_list, bins=100)
    plt.title(organism)
    f.savefig(path + "/" + file.replace(".bed", "_" + str(organism) + "_plot.pdf"), bbox_inches='tight')


def save_as_csv(file, path, normalised_position, organism):
    pd.DataFrame([float(elem) for elem in normalised_position]).to_csv(
        path + '/' + file.replace(".bed", "_" + organism + "_pos.csv"),
        header=None,
        index=None)


@click.command()
@click.option('--input_file', '-i', help='Name of the input file', required=True)
@click.option('--input_path', '-p', help='Name of the input path', required=True)
@click.option('--plot_samples', '-s', default=False, help='Should the samples be plotted? Default False')
@click.option('--min_read_count', '-min', default=1000, help='Minimum amount of reads per sample. Default 1000')
def main(input_file, input_path, plot_samples, min_read_count):
    df_input_file = read_bed(input_file, input_path)
    get_pos(df_input_file, input_file, input_path, plot_samples, min_read_count)


if __name__ == "__main__":
    main()
