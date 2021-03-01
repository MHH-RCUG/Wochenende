"""
A reporting script for the results of the Wochenende Pipeline
Author: Tobias Scheithauer, August 2018
Author: Erik Wessels, October 2019
Author: Fabian Charly Friedrich, February 2020
Author: Colin Davenport
Author: Sophia Poertner

TODOs

Changelog
1.5.4 new references added
1.5.3 docs and comments improved
1.5.2 lint with tool black
1.5.1 add debug switch for printing, otherwise in quiet mode.
1.5 comment out pysam as can cause miniconda cryptic error and currently unused (used in slow BAM methods only), bugfix RPMM
1.4 add version, reduce output file name length. reporting -> rep, sorted -> s etc
1.3 restructured the script
1.2 bugfix reads_per_human_cell
1.1 add reads_per_human_cell for metagenomes
1.0 implement the script
"""

import sys
import os
from Bio import SeqIO, SeqUtils

# import pysam
import pandas as pd
import click


version = "1.5.4"

###################################
# functions
##################################
def print_start(input_file, reference, sequencer):
    click.echo(f"starting {sequencer} reporting")
    click.echo(f'Started {"bam" if input_file[-4:] == ".bam" else "txt"} mode.')
    click.echo(f"Using {input_file} as alignment file")
    click.echo(f"Using {reference} as reference")
    click.echo()


def check_arguments(input_file, reference, sequencer):
    # check if input file exist and if it is the correct file format
    if os.path.isfile(input_file):
        if not input_file.endswith(".bam.txt") and not input_file.endswith(".bam"):
            click.echo(
                f"The input file {input_file} \t has a incorrect file format.\nType --help for help"
            )
            sys.exit()
    else:
        click.echo(
            f"The input file {input_file} \t does not exist\nType --help for help"
        )
        sys.exit()

    # check if refseq file exist and if it is the correct file format
    if os.path.isfile(reference):
        if (
            not reference.endswith(".fasta")
            and not reference.endswith(".fa")
            and not reference.endswith(".fna")
        ):
            click.echo(
                f"The input file {reference} \t has a incorrect file format.\nType --help for help"
            )
            sys.exit()
    else:
        click.echo(
            f"The input file {reference} \t does not exist\nType --help for help"
        )
        sys.exit()
    # check the sequencer
    if not sequencer == "illumina" and not sequencer == "solid":
        click.echo(
            f"The sequencer {sequencer} \t is not supported\nType --help for help"
        )
        sys.exit()
    click.echo("Arguments checked successfully")


def create_res_df(input_file, reference):
    if input_file[-4:] == ".bam":
        return create_res_df_from_bam(input_file, reference)
    else:
        return create_res_df_from_txt(input_file, reference)


def create_res_df_from_txt(input_file, refseq_file):
    res_df = pd.read_csv(
        input_file,
        sep="\t",
        header=None,
        names=["species", "chr_length", "read_count"],
        usecols=[0, 1, 2],
    )[:-1]
    # get gc content of ref sequences
    gc_ref_dict = {}
    for seq_record in SeqIO.parse(refseq_file, "fasta"):
        gc_ref_dict[seq_record.name] = SeqUtils.GC(str(seq_record.seq).replace("N", ""))
    res_df["gc_ref"] = [gc_ref_dict[s] for s in res_df["species"]]
    return res_df


def create_res_df_from_bam(input_file, reference):
    (
        species_list,
        chr_length_list,
        read_count_list,
        basecount_list,
        gc_ref_list,
        gc_reads_list,
    ) = ([], [], [], [], [], [])

    for seq_record in SeqIO.parse(reference, "fasta"):
        # joining all reads
        joined_reads = "".join(
            [
                read.query_sequence
                for read in pysam.AlignmentFile(input_file, "rb").fetch(
                    contig=seq_record.name
                )
            ]
        )

        # appending to all Lists
        species_list.append(seq_record.name)
        chr_length_list.append(len(seq_record.seq))
        read_count_list.append(
            pysam.AlignmentFile(input_file, "rb").count(contig=seq_record.name)
        )
        gc_ref_list.append(SeqUtils.GC(seq_record.seq))
        gc_reads_list.append(SeqUtils.GC(joined_reads))
        basecount_list.append(sum([len(joined_reads)]))

    # create and return dataframe
    return pd.DataFrame(
        data={
            "species": species_list,
            "chr_length": chr_length_list,
            "gc_ref": gc_ref_list,
            "gc_reads": gc_reads_list,
            "read_count": read_count_list,
            "basecount": basecount_list,
        }
    )


def solid_gc_normalization(gc):
    # TODO implement a linear scale model for the solid reads
    return 1


def solid_normalization(res_df):
    res_df["norm_factor"] = [solid_gc_normalization(res_df["gc_ref"])]
    res_df["read_count"] = res_df["read_count"] * res_df["norm_factor"]
    res_df["basecount"] = res_df["basecount"] * res_df["norm_factor"]
    res_df["norm_factor"] = None
    return res_df


def bacteria_per_human_cell(res_df):
    # calculating bacteria per human cell
    # check for human_refs to be correct!
    # the mitochondrial reads have NOT been added to the sum of human reads, as the bacteria/human ratio  would have been extremly small.
    # this needs to be further discussed
    human_refs = [
        "1_1_1_1",
        "1_1_1_2",
        "1_1_1_3",
        "1_1_1_4",
        "1_1_1_5",
        "1_1_1_6",
        "1_1_1_7",
        "1_1_1_8",
        "1_1_1_9",
        "1_1_1_10",
        "1_1_1_11",
        "1_1_1_12",
        "1_1_1_13",
        "1_1_1_14",
        "1_1_1_15",
        "1_1_1_16",
        "1_1_1_17",
        "1_1_1_18",
        "1_1_1_19",
        "1_1_1_20",
        "1_1_1_21",
        "1_1_1_22",
        "1_1_1_X",
        "1_1_1_Y",
    ]
    human_cov = res_df[res_df["species"].isin(human_refs)]["read_count"].sum()
    res_df["bacteria_per_human_cell"] = (
        6191.39 * res_df["reads_per_million_ref_bases"]
    ) / human_cov
    return res_df


def compute_res_df(res_df, sequencer):

    # special solid normalization
    # if sequencer == 'solid': res_df = solid_normalization(res_df)

    # calculating reads_per_million_ref_bases
    res_df["reads_per_million_ref_bases"] = res_df["read_count"] / (
        res_df["chr_length"] / 1000000
    ).round(5)

    # calculating reads_per_million_reads_in_experiment
    res_df["reads_per_million_reads_in_experiment"] = res_df["read_count"] / (
        res_df["read_count"].sum() / 1000000
    ).round(5)

    # total normalization RPMM. Corrected in autumn 2020 after bug found.
    res_df["RPMM"] = res_df["read_count"] / (
        (res_df["chr_length"] / 1000000) * (res_df["read_count"].sum() / 1000000)
    ).round(5)

    # calculating bacteria per human cell
    res_df = bacteria_per_human_cell(res_df)

    return res_df


def export_res(res_df, output):
    res_df.to_csv(f"{output}.rep.us.csv", sep=",", float_format="%.5f", index=False)
    res_df_filtered_and_sorted = res_df.loc[res_df["read_count"] >= 20].sort_values(
        by="RPMM", ascending=False
    )
    res_df_filtered_and_sorted.to_csv(
        f"{output}.rep.s.csv", sep=",", float_format="%.5f", index=False
    )


###################################
# main
##################################
@click.command()
@click.option(
    "--input_file",
    "-i",
    help="File in .bam.txt (recommended) or .bam format (slow) from the Wochenende pipeline output",
)
@click.option(
    "--reference",
    "-r",
    help="File in .fasta format. This must be the reference used by the Wochenende pipeline",
)
@click.option(
    "--sequencer",
    "-s",
    default="illumina",
    help="Sequencer technology used. Only illumina is supported (prev. SOLiD), default: illumina ",
)
@click.option(
    "--output_name",
    "-o",
    default="report",
    help="Name for the output file(sample name), default report",
)
def main(input_file, reference, sequencer, output_name):
    """
    This Python3.6 script can be used to report the results of the Wochenende pipeline.
    The .bam.txt file as input is recommended (fast).
    The .bam file will take longer and generate more information.

    The column reads_per_human_cell is only for metagenomes from human hosts.

    Reports for solid sequencing data are not supported, a special GC normalisation model has to be implemented first.
    """
    # Set debug = True to print debugging info such as input file, reference, sequencer
    debug = False
    if debug:
        print(version)

    check_arguments(input_file, reference, sequencer)
    if debug:
        print_start(input_file, reference, sequencer)
    res_df = create_res_df(input_file, reference)
    res_df = compute_res_df(res_df, sequencer)
    export_res(res_df, output_name)


if __name__ == "__main__":
    main()
