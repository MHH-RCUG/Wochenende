#!/usr/bin/env python
# coding: utf-8


# Author: Konstantinos Sifakis, Sophia Poertner, Colin Davenport
# Script written: Jan 2020-Feb 2020, adapted DEC 2020-JAN 2020
# Script created for: MHH - RCUG
# # Plot coverage data after the runbatch_sambamba_depth.sh sambamba coverage step of the Wochenende pipeline.
# Create some png files,one txt that includes the organisms that seem related
# and a tsv file with information from tsv file that was given.
# Usage: python3 wochenende_plot.py dup.calmd_cov_window.txt.filt.csv


# imports:

import argparse
import pandas as pd
import matplotlib
import os
import statistics
import warnings
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Changelog

# 0.2.2 change perhaps present to low_med_score and prob_present to high_score
# 0.2.1 default mean coverage 0.2# 0.2 reduce --minMeanCov default from 1.0 to 0.2 after testing, add version.
# 0.1 initial working script

version = "0.2.2"

# use_line_collection should be True

warnings.filterwarnings("ignore", category=UserWarning)

# Needed to read file  [ --filename ]  and then [ file directory ] so that file would open as a tsv file

parser = argparse.ArgumentParser(
    description="Plot coverage data after the sambamba (or wochenende_postprocess.sh) coverage step of the Wochenende "
                "pipeline. Creates png files of absolute and mean coverage per taxon, and a statistical overview file."
)
parser.add_argument(
    "filename1",
    default="check_string_for_empty",
    help="Usage: python3 wochenende_plot.py dup.calmd_cov_window.txt.filt.csv",
)
parser.add_argument(
    "--minMeanCov",
    default=0.2,
    help="""The minimum mean coverage. Default: 0.2
                    This does not affect the score. After the score has been calculated, this number will be the 
                    threshold for both low scored and highly scoring taxa.\n
                    ## Case1: Poor scoring that have a maximum score above this minMeanCov value -> Images go to folder 
                    low_med_score.\n
                    ## Case2: Poor scoring will have a maximum score below or equal to this minMeanCov value -> Images 
                    will not be printed (unless --createAllPngs is given).\n
                    ## Case3: Highly scored taxa that have a minimum score above this minMeanCov value -> Images will 
                    go to folder high_score\n
                    ## Case4: Highly scored taxa that have a minimum score below or equal to this minMeanCov value -> 
                    Images will go to folder low_med_score.\n
                    """,
)
parser.add_argument(
    "--createAllPngs",
    default=0,
    help="Choose if you want to create figures for the poorly scored taxa which also have their max "
         "(mean_coverage) <= cov input argument(--minMeanCov)\nIf you wish to create the figures you "
         "should write --createAllPngs 1\n Default: 0, do not create poor scoring figures",
)
parser.add_argument(
    "--sclim",
    default=0.0,
    help="Choose the score limit to filter whether a taxon is going to be well or poorly "
         "scored\nDefault: 0.0, negative scores are termed bad, positive values are highly scored.",
)
parser.add_argument(
    "--minWindows",
    default=5,
    help="Choose the least number of (default - set in runbatch_sambamba_depth.sh) 100kbp windows per organism that should be taken into account.\n "
         "Generally, more then 5 is suggested, otherwise quartiles q1 or Q3 will be useless. Since we "
         "are focussing on bacteria and bins are 100,000 bp length (50,000bp overlap), there should be "
         "no point on having less then 5 window bins accepted. Default: 5.",
)
args = parser.parse_args()

# _____________________________________________________

# transfer arguments to variables

filename = args.filename1
plot_all = int(args.createAllPngs)
Scorelimitforminormax = float(args.minMeanCov)
Scorelimit = float(args.sclim)
minWindowsNeeded = int(args.minWindows)

# CONSTANTS
MY_DPI = 96
MULTIPDPI = 3


# AnotherScorelimit=2.0
# defaultorgnames = 2  # not used??


# functions ____________________________________________


def shorten_filenames():
    if filename == "":
        print("#### No argument was given")
    else:
        print("#### Filename entered correctly ####")
        # Shorten filename, otherwise it causes errors on windows systems
        # typical filename: KGCF14D_S2_R1.ndp.lc.trm.s.mq30.01mm.dup.bam.txt.rep.s.csv
        index_dot = filename.find(".")
        outfile = filename[:index_dot]

        return outfile


def create_directories(dirpath, outfile):
    # make new directories
    images_path = "{}/images".format(dirpath)
    sample_outpath = "{}/images/{}".format(dirpath, outfile)
    perhaps_path = "{}/images/{}/low_med_score".format(dirpath, outfile)
    prob_path = "{}/images/{}/high_score".format(dirpath, outfile)

    for path in (images_path, sample_outpath, perhaps_path, prob_path):
        try:
            os.mkdir(path)
        except OSError:
            print("Creation of the directory %s failed" % path)
            print(
                "###__Maybe this directory already exists__###__OR you are not allowed to create it__###"
            )
        else:
            print("Successfully created the directory %s " % path)
    return images_path, sample_outpath, perhaps_path, prob_path


def del_create_output_files(prob_path):
    # Delete the old .txt files that keeps the scores and GenBankID (from a previous run)
    if os.path.isfile("{}/Scores.csv".format(prob_path)):
        os.remove("{}/Scores.csv".format(prob_path))
    if os.path.isfile("{}/GenBankID.txt".format(prob_path)):
        os.remove("{}/GenBankID.txt".format(prob_path))

    # Make a new .txt file for Scores and Yes or No column
    f_score = open("{}/Scores.csv".format(prob_path), "a")
    f_score.write(
        "GenBank_ID\tGenus\tSpecies\tStrain\tStrain_info\tscore\tScorelimit>={}\tMight_be_related_to_smthg".format(
            Scorelimit
        )
    )
    genbankid = open("{}/GenBankID.txt".format(prob_path), "a")
    return genbankid, f_score


def shorten_organism_names(organism):
    # Shorten names of the reference sequence, otherwise it causes errors on windows systems
    # typical organism name: 1_AE004091_2_Pseudomonas_aeruginosa_PAO1__complete_genome_BAC
    tmpfile = organism
    tmpfile = tmpfile.replace("complete_genome", "")
    tmpfile = tmpfile.replace("complete", "")
    tmpfile = tmpfile.replace("sequence", "")
    tmpfile = tmpfile.replace("chromosome", "")
    tmpfile = tmpfile.replace(".ndp.lc", "")
    tmpfile = tmpfile.replace(".trm.s", "")
    tmpfile = tmpfile.replace(".mq30.01mm", "")
    tmpfile = tmpfile.replace(".dup.bam", "")
    tmpfile = tmpfile.replace(".rep.s.csv", "")
    tmpfile = tmpfile.replace(".dup_cov_window", "")
    tmpfile = tmpfile.replace(".txt.filt", "")
    tmpfile = tmpfile.replace(".csv", "")
    tmpfile = tmpfile.replace("__", "_")
    organism = tmpfile
    return organism


def plot_figures(
        df_organism,
        score,
        scoreforlowmean,
        organism,
        passed_q1,
        q1,
        passed_ms,
        mean_2_stdv,
        max_chr_start,
        simplified_max_mean_coverage,
):
    # Temporary rc parameters in effect
    fig, ax1 = plt.subplots()
    # plt.figure(figsize=(48,36), dpi= 70)
    plt.figure(figsize=(8000 / MY_DPI, 8000 / MY_DPI), dpi=MY_DPI)
    # Used to change name in Figure-Title depending on results. And title will include info about Windows and
    # the Mean Coverage each one has.
    # score=round(score,2)
    #            fig.suptitle('score= {}'.format(scoreforlowmean), fontsize=12,fontweight='bold')
    # Scorelimit is adjustable some lines above(at the top)
    if score >= Scorelimit and df_organism["meanCoverage"].shape[0] >= minWindowsNeeded:
        fig.suptitle(
            "score= {}".format(scoreforlowmean), fontsize=12, fontweight="bold"
        )
        ax1.set_title(
            "\n{} \n\n {}% bins had mean Coverage above {} \n {}% bins had mean Coverage above {}".format(
                organism, passed_q1, q1, passed_ms, mean_2_stdv
            ),
            fontdict={"fontsize": 8, "fontweight": "medium"},
        )

    elif (
            score >= Scorelimit and df_organism["meanCoverage"].shape[0] <= minWindowsNeeded
    ):  # Scorelimit is adjustable some lines above(at the top)

        fig.suptitle(
            "score= {} ,however data are not enough".format(scoreforlowmean),
            fontsize=12,
            fontweight="bold",
        )
        ax1.set_title(
            "\n{} \n\n {}% bins had mean Coverage above {} \n {}% bins had mean Coverage below {}".format(
                organism, passed_q1, q1, passed_ms, mean_2_stdv
            ),
            fontdict={"fontsize": 8, "fontweight": "medium"},
        )

    elif (
            score < Scorelimit and df_organism["meanCoverage"].shape[0] < minWindowsNeeded
    ):  # Scorelimit is adjustable some lines above(at the top)
        if q1 == mean_2_stdv:
            fig.suptitle(
                "low score= {}".format(scoreforlowmean), fontsize=12, fontweight="bold"
            )
            ax1.set_title(
                "\n{}\n\n{}% bins had mean Coverage equal to {}".format(
                    organism, passed_ms, mean_2_stdv
                ),
                fontdict={"fontsize": 8, "fontweight": "medium"},
            )

        else:
            fig.suptitle(
                "low score= {}".format(scoreforlowmean), fontsize=12, fontweight="bold"
            )
            ax1.set_title(
                "\n{}\n\n{}% bins had mean Coverage below or equal to {} \n {}% bins had mean Coverage above or equal "
                "to{}".format(
                    organism, passed_ms, mean_2_stdv, passed_q1, q1
                ),
                fontdict={"fontsize": 8, "fontweight": "medium"},
            )

    else:
        if q1 == mean_2_stdv:
            fig.suptitle(
                "low score= {}".format(scoreforlowmean), fontsize=12, fontweight="bold"
            )
            ax1.set_title(
                "\n{}\n\n{}% bins had mean Coverage above or equal to {}".format(
                    organism, passed_ms, mean_2_stdv
                ),
                fontdict={"fontsize": 8, "fontweight": "medium"},
            )

        else:
            fig.suptitle(
                "low score= {}".format(scoreforlowmean), fontsize=12, fontweight="bold"
            )
            ax1.set_title(
                "\n{}\n\n{}% bins had mean Coverage above or equal to {} \n {}% bins had mean Coverage above or equal "
                "to{}".format(
                    organism, passed_ms, mean_2_stdv, passed_q1, q1
                ),
                fontdict={"fontsize": 8, "fontweight": "medium"},
            )
    ax1.set(ylabel="Read Counts", ylim=(0, 30))
    ax1.set(xlabel="Position", xlim=(0, max_chr_start))
    ax1.yaxis.label.set_color("blue")
    # ax1.stem(df_organism["chromStart"],df_organism["readCount"], markerfmt='bo',linefmt="tab:green",
    # use_line_collection=True)
    markers, stems, base = ax1.stem(
        df_organism["chromStart"], df_organism["readCount"], markerfmt="bo", linefmt="tab:cyan"
    )
    # TODO stem not iterable any more matplotlib 3.3.4
    # https://matplotlib.org/stable/gallery/shapes_and_collections/line_collection.html
    #stems.set_linewidth(10)
    #for stem in stems:
        # use_line_collection = True
        # stem.set_linewidth(10)
        # stem.set_use_line_collection=True
    ax1.tick_params(axis="y", labelcolor="blue")
    ax1.grid()
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
    # print("MAX mean coverage:",np.ceil(max(B1["meanCoverage"])))
    # print("MAX read counts:", np.ceil(max (B1["readCount"])))
    # assign("last.warning"=NULL,envir=baseenv())
    ax2.set(ylabel="Mean Coverage", ylim=(0, simplified_max_mean_coverage))
    ax2.yaxis.label.set_color("red")
    ax2.stem(
        df_organism["chromStart"],
        df_organism["meanCoverage"],
        markerfmt="rd",
        linefmt="white",
        use_line_collection=True,
    )  # bottom 1 could be used for those with a good score and bad Mean Coverage
    ax2.tick_params(axis="y", labelcolor="tab:red")
    fig.tight_layout()
    return fig


def split_organism_names(organism):
    if len(organism.split("_")) >= 7:
        organism0 = organism.split("_")[1]
        organism1 = organism.split("_")[3]
        organism2 = organism.split("_")[4]
        organism3 = organism.split("_")[5]
        organism4 = organism.split("_")[6]
    else:
        print(
            "Warning:\n  Ref.-Seq.",
            organism,
            "has the wrong name format (only affects the score file)"
        )
        organism0 = organism
        organism1 = "_"
        organism2 = "_"
        organism3 = "_"
        organism4 = "_"
    for organism_part in (organism1, organism2, organism3, organism4):
        if organism_part == "":
            organism1 = "_"  # but why?
    return organism0, organism1, organism2, organism3, organism4


def write_in_files(scoreforlowmean,
                   min_mean_coverage,
                   df_organism,
                   max_mean_coverage,
                   genbankid,
                   organism0,
                   organism1,
                   organism2,
                   organism3,
                   organism4,
                   f_score,
                   scoreforlowmeantxt):
    # NEW COLUMN FOR MORE SCORES
    # depending on the scoreforlowmean, the min_mean_coverage and the len(df_organism["mean_coverage"]) different things
    # are written in the score file for the organisms. Either "yes" or "no" for the columns "Scorelimit>=0.0" and
    # "Might_be_related_to_smthg" ( YES NO/ YES YES/ NO NO/ NO YES )
    # Only if all three are greater than the limit, the GenbankID of the organism is written in the GenBankID.txt.
    # These are also the plots in the high_cov directory
    if (
            scoreforlowmean >= Scorelimit
            and min_mean_coverage >= Scorelimitforminormax
            and len(df_organism["meanCoverage"] >= minWindowsNeeded)
    ):
        genbankid.write("\n%s" % organism0)
        # round_max_mean_coverage
        #                if scoreforlowmean>AnotherScorelimit:
        # f_score.write('%s\t%s\t---------->>>>>>YES<<<<<<----------\n' % (organism, scoreforlowmeantxt))
        f_score.write(
            "\n%s\t%s\t%s\t%s\t%s\t%s\t---------->>>>>>YES<<<<<<----------\t---------->>>>>>YES<<<<<<----------"
            % (
                organism0,
                organism1,
                organism2,
                organism3,
                organism4,
                scoreforlowmeantxt,
            )
        )
    elif (
            scoreforlowmean >= Scorelimit
            and min_mean_coverage < Scorelimitforminormax
            and len(df_organism["meanCoverage"] >= minWindowsNeeded)
    ):
        f_score.write(
            "\n%s\t%s\t%s\t%s\t%s\t%s\t---------->>>>>>YES<<<<<<----------\t---------->>>>>>YES<<<<<<----------"
            % (
                organism0,
                organism1,
                organism2,
                organism3,
                organism4,
                scoreforlowmeantxt,
            )
        )
    elif (
            scoreforlowmean < Scorelimit
            and max_mean_coverage >= Scorelimitforminormax
            and len(df_organism["meanCoverage"] >= minWindowsNeeded)
    ):
        f_score.write(
            "\n%s\t%s\t%s\t%s\t%s\t%s\tNO\t---------->>>>>>YES<<<<<<----------"
            % (
                organism0,
                organism1,
                organism2,
                organism3,
                organism4,
                scoreforlowmeantxt,
            )
        )
    else:
        f_score.write(
            "\n%s\t%s\t%s\t%s\t%s\t%s\tNO\tNO"
            % (
                organism0,
                organism1,
                organism2,
                organism3,
                organism4,
                scoreforlowmeantxt,
            )
        )


def save_figures(scoreforlowmean,
                 min_mean_coverage,
                 df_organism,
                 fig,
                 prob_path,
                 organism,
                 round_score,
                 perhaps_path,
                 round_max_mean_coverage,
                 max_mean_coverage,
                 sample_outpath,
                 prob,
                 perhaps,
                 other
                 ):
    # RENAMING the figures depending on their score
    #                         &
    # Making the csv of the Scores  \t (tab) separated
    if (
            scoreforlowmean >= Scorelimit
            and min_mean_coverage >= Scorelimitforminormax
            and len(df_organism["meanCoverage"] >= minWindowsNeeded)
    ):
        prob += 1
        # fig.savefig('{}/{}_MC_RCSCORED{}.png'.format(prob_path,organism,round_score),dpi=MY_DPI*MULTIPDPI)
        fig.savefig("{}/{}_{}.png".format(prob_path, organism, round_score), dpi=MY_DPI * MULTIPDPI)
    elif (
            scoreforlowmean >= Scorelimit
            and min_mean_coverage < Scorelimitforminormax
            and len(df_organism["meanCoverage"] >= minWindowsNeeded)
    ):
        perhaps += 1
        # fig.savefig('{}/{}_MC_RC_plot_all_Might{}.png'.format(perhaps_path,organism,round_max_mean_coverage),dpi=MY_DPI*MULTIPDPI)
        fig.savefig("{}/{}_{}.png".format(perhaps_path, organism, round_max_mean_coverage), dpi=MY_DPI * MULTIPDPI)
    elif (
            scoreforlowmean < Scorelimit
            and max_mean_coverage >= Scorelimitforminormax
            and len(df_organism["meanCoverage"] >= minWindowsNeeded)
    ):
        perhaps += 1
        # fig.savefig('{}/{}_MC_RC_plot_all_Might{}.png'.format(perhaps_path,organism,round_max_mean_coverage),dpi=MY_DPI*MULTIPDPI)
        fig.savefig("{}/{}_{}.png".format(perhaps_path, organism, round_max_mean_coverage), dpi=MY_DPI * MULTIPDPI)
    else:
        other += 1
        # No reason anymore to create so many figures of none usefull data-organisms
        if plot_all == 1:
            fig.savefig(
                "{}/{}_MC_RC_plot_all.png".format(sample_outpath, organism), dpi=MY_DPI * MULTIPDPI
            )
        else:
            pass
        # This can be adjusted. Would have been an alternative second score.
    return prob, perhaps, other


def generate_output(
        df_organism,
        score,
        scoreforlowmean,
        organism,
        passed_q1,
        q1,
        passed_ms,
        mean_2_stdv,
        max_chr_start,
        simplified_max_mean_coverage,
        max_mean_coverage,
        min_mean_coverage,
        prob_path,
        genbankid,
        f_score,
        scoreforlowmeantxt,
        perhaps_path,
        sample_outpath,
        prob,
        perhaps,
        other
):
    with plt.rc_context(
            {
                "axes.edgecolor": "orange",
                "xtick.color": "black",
                "ytick.color": "green",
                "figure.facecolor": "white",
            }
    ):
        fig = plot_figures(
            df_organism,
            score,
            scoreforlowmean,
            organism,
            passed_q1,
            q1,
            passed_ms,
            mean_2_stdv,
            max_chr_start,
            simplified_max_mean_coverage,
        )
        round_score = round(score)
        round_max_mean_coverage = round(max_mean_coverage)
        organism0, organism1, organism2, organism3, organism4 = split_organism_names(organism)
        prob, perhaps, other = save_figures(scoreforlowmean,
                                            min_mean_coverage,
                                            df_organism,
                                            fig,
                                            prob_path,
                                            organism,
                                            round_score,
                                            perhaps_path,
                                            round_max_mean_coverage,
                                            max_mean_coverage,
                                            sample_outpath,
                                            prob,
                                            perhaps,
                                            other)
        write_in_files(scoreforlowmean,
                       min_mean_coverage,
                       df_organism,
                       max_mean_coverage,
                       genbankid,
                       organism0,
                       organism1,
                       organism2,
                       organism3,
                       organism4,
                       f_score,
                       scoreforlowmeantxt)
        plt.close("all")
    return prob, perhaps, other


def calculate_scores(
        df_input_file,
        prob_path,
        genbankid,
        f_score,
        perhaps_path,
        sample_outpath,
):
    # counting how many organisms are high_score, low_med_score (images created) and probably not present (no image created)
    prob = 0
    perhaps = 0
    other = 0
    for organism in set(df_input_file.index.values):

        # Read Counts
        # __________________________________________________________________________________#

        # df for each organism, find min&max readcount&chromstart for each chromsome
        df_organism = df_input_file.loc[[organism]]
        # B1maxread = max(df_organism["readCount"])
        # B1minread = min(df_organism["readCount"])
        max_chr_start = max(df_organism["chromStart"])
        # B1minStart = min(df_organism["chromStart"])

        organism = shorten_organism_names(organism)
        # MEAN COVERAGE
        # ____________________________________________________________________________________#

        # max_mean_coverage will show us, in which directory to put the plots.
        # The later simplified verison ist for the axis
        max_mean_coverage = max(df_organism["meanCoverage"])

        min_mean_coverage = min(df_organism["meanCoverage"])
        bin_size = (max(df_organism["chromEnd"]) - max(df_organism["chromStart"])) / 2
        # /2 because every window is "covered" 2 times
        # Bin_size_Multi_howmany = bin_size * (
        #         len(df_organism["meanCoverage"]) + 1
        # )  # +1 Correction for overlaping windows ----- not used?
        len_empty = max(df_organism["chromEnd"]) - (len(df_organism["meanCoverage"]) * bin_size)
        howmany_empty_bins = len_empty / bin_size
        howmany_empty_bins_cor = (
                howmany_empty_bins - 1
        )  # Correction for overlaping windows
        #        print(howmany_empty_bins,"------",organism)
        #        print(bin_size,"bin_size",organism)
        #        print(Bin_size_Multi_howmany,"Bin_size_Multi_howmany")
        #        print(howmany_empty_bins_cor,"Howmany Bins?")
        #        print(max(df_organism["chromEnd"]),"______________",organism)
        j = max_mean_coverage
        # simplified_max_mean_coverage = k
        if 0 <= j <= 1:
            k = 1
        elif 1 <= j <= 2:
            k = 2
        elif 1 <= j <= 5:
            k = 5
        elif 5 <= j <= 10:
            k = 10
        elif j / 10 >= 1 and j <= 20:
            k = 20
        elif j / 20 >= 1 and j <= 50:
            k = 50
        elif 1 <= 50 < j <= 100:
            k = 100
        else:
            k = 200
        simplified_max_mean_coverage = k
        # Here i am doing the scaling for axis_ not the best way to do it tho............

        score = statistics.mean(df_organism["meanCoverage"])  # mean meanCoverage

        # We will name the average mean.coverage of our samples =score and this will help us remember that (average)
        # should ( ~ ) the final (score)

        # NEW SCORE:
        if len(df_organism["meanCoverage"]) > 10:
            # dividedby=statistics.stdev(df_organism["meanCoverage"])            

            norm_cov_list = []
            for meanCove in df_organism["meanCoverage"]:
                # mean_Cove - min_mean_coverage -> the lowest point is now 0
                # / (simplified_max_mean_coverage - min_mean_coverage) -> the highest point is now 1
                # normalize the Data so it is between 0 and 1
                normed_coverage = (meanCove - min_mean_coverage) / (simplified_max_mean_coverage - min_mean_coverage)
                norm_cov_list.append(normed_coverage)
            norm_mean_cov = statistics.mean(norm_cov_list)  # normed mean of mean

            dividedby = norm_mean_cov
            score2 = (
                             (score * len(df_organism["meanCoverage"])) - (howmany_empty_bins_cor * score)
                     ) / ((dividedby + 1) * len(df_organism["meanCoverage"]))
            # /dividedby should not be 0 ever + if it's less then 1 then it automatically becomes a sort of
            # multiplication
            # So we should use Divided by +1 in order to be sure that this will be working as it supposed to be
            # working... Reusing our score every time depending on
            # how much the NORMALIZED Standard Deviation is. +1 in order to really have a good result serving our goal
            # to reduce score

            # we also divide by the len(df_organism["meanCoverage"]) which is basically our data so we divide by the
            # number of the data so that we manipulate score,
            # and it should not get too high if there are too many but low mean Coverage windows.

            score = score2  # something could be done for more score calculations

        elif 2 <= len(df_organism["meanCoverage"]) <= 10:
            norm_cov_list = []
            for meanCove in df_organism["meanCoverage"]:
                normed_coverage = (meanCove - min_mean_coverage) / (simplified_max_mean_coverage - min_mean_coverage)
                norm_cov_list.append(normed_coverage)
            norm_mean_cov = statistics.mean(norm_cov_list)

            dividedby = norm_mean_cov
            score2 = (
                             (
                                     (score * len(df_organism["meanCoverage"]))
                                     - (howmany_empty_bins_cor * score)
                             )
                             / ((dividedby + 1) * len(df_organism["meanCoverage"]))
                     ) - (
                             2 * min_mean_coverage
                     )  # penalty
            # for having less then 10 windows score - 2 times the min of those windows.
            score = score2  # MORE to be done? maybe ?.
            # NEW SCORE
            # print(dividedby,"__Divided_by___",organism)
            # print(score,"SCORED2")
            # That should happen only if not good mean and not good length of the df_organism data.frame
        else:
            score = -100  # Typical not good score
            # GIVE number of windows that are above mean-2*sd and those that are above q1 of the df_organism data.frame
            # minWindowsNeeded=5       # more then 5 is suggested otherwise q1 and Q3 are useless.....

        # Those are more or less describing a panda.describe() function input
        mean_describe = 1
        std_describe = 2
        min_describe = 3
        quartile_1_describe = 4
        # Q2_describe = 5  # Not used
        quartile_3_describe = 6
        max_describe = 7

        mean_2_stdv = df_organism.describe()["meanCoverage"][mean_describe] - (
                2 * df_organism.describe()["meanCoverage"][std_describe]
        )  # mean - 2 * std

        # Checking if mean is more then 0 and if df_organism has more then 5 inputs.
        if df_organism["meanCoverage"].shape[0] >= minWindowsNeeded and mean_2_stdv >= 0:
            q1 = df_organism.describe()["meanCoverage"][quartile_1_describe]
            passed_q1_f = df_organism[df_organism["meanCoverage"] >= q1].shape[0]
            passed_q1 = (passed_q1_f * 100) / df_organism["meanCoverage"].shape[0]
            # Not needed , we have calculated it before if.
            # mean_2_stdv=df_organism.describe()["meanCoverage"][mean_describe]-(2*df_organism.describe()["meanCoverage"][std_describe])
            passed_msf = df_organism[df_organism["meanCoverage"] >= mean_2_stdv].shape[0]
            passed_ms = (passed_msf * 100) / df_organism["meanCoverage"].shape[0]

            # Many figures and their titles describing the bins are showing that normal destribution is not followed
            # always but Poisson seems more accurate distribution
            # in order to describe our data.

            # a new column into score.txt file would be added in regards to which distribution is followed in the future

            if passed_ms >= 100:
                # print("Normal Distribution??")
                mean_2_stdv = df_organism.describe()["meanCoverage"][min_describe]
                passed_msf = df_organism[df_organism["meanCoverage"] >= mean_2_stdv].shape[0]
                passed_ms = (passed_msf * 100) / df_organism["meanCoverage"].shape[0]

        # WORKING.................
        # IF LOW SCORE SHOW THE MAX and everything else is above that max in the figure title.
        # NOT DONE

        elif df_organism["meanCoverage"].shape[0] >= minWindowsNeeded and mean_2_stdv < 0:
            q1 = df_organism.describe()["meanCoverage"][quartile_1_describe]
            passed_q1_f = df_organism[df_organism["meanCoverage"] >= q1].shape[0]
            passed_q1 = (passed_q1_f * 100) / df_organism["meanCoverage"].shape[0]
            # mean not usefull lets use Q3 this time
            mean_2_stdv = df_organism.describe()["meanCoverage"][
                quartile_3_describe
            ]  # this is Q3 with the name of mean_2_stdv
            passed_msf = df_organism[df_organism["meanCoverage"] >= mean_2_stdv].shape[0]
            passed_ms = (passed_msf * 100) / df_organism["meanCoverage"].shape[0]

        # This is part of an old version and should not be called at all with then new changes made( the_least windowss
        # that we need is >=5)

        else:  # We will use min and max cause not enough data to get a proper q1 or mean or Q3.
            q1 = df_organism.describe()["meanCoverage"][
                min_describe
            ]  # this would be our min with a different name
            passed_q1_f = df_organism[df_organism["meanCoverage"] >= q1].shape[0]
            passed_q1 = (passed_q1_f * 100) / df_organism["meanCoverage"].shape[0]

            mean_2_stdv = df_organism.describe()["meanCoverage"][
                max_describe
            ]  # this would be out max with different name
            passed_msf = df_organism[df_organism["meanCoverage"] <= mean_2_stdv].shape[0]
            passed_ms = (passed_msf * 100) / df_organism["meanCoverage"].shape[0]

        passed_ms = round(passed_ms, 2)
        passed_q1 = round(passed_q1, 2)
        mean_2_stdv = round(mean_2_stdv, 5)
        q1 = round(q1, 5)
        scoreforlowmeantxt = score
        scoreforlowmean = round(score, 4)
        if scoreforlowmean == 0.0:
            scoreforlowmean = round(score, 10)

        # plot Figures, write in Scores file
        prob, perhaps, other = generate_output(
            df_organism,
            score,
            scoreforlowmean,
            organism,
            passed_q1,
            q1,
            passed_ms,
            mean_2_stdv,
            max_chr_start,
            simplified_max_mean_coverage,
            max_mean_coverage,
            min_mean_coverage,
            prob_path,
            genbankid,
            f_score,
            scoreforlowmeantxt,
            perhaps_path,
            sample_outpath,
            prob,
            perhaps,
            other
        )
    return prob, perhaps, other


def main():
    # info and get directory path + outfile
    print("\n###############################")
    print("Wochenende plotting, Version: %s" % version)
    print("###############################\n")
    outfile = shorten_filenames()
    dirpath = os.getcwd()
    print("INFO: Current directory is : " + dirpath)
    print("INFO: Current file path is : " + dirpath + "/" + filename)
    print("INFO: Input filename: " + str(filename))
    print("INFO: Output filename: %s" % outfile)

    df_input_file = pd.read_csv("{}/{}".format(dirpath, filename), delimiter="\t", index_col="# chrom")

    # create output directories
    images_path, sample_outpath, perhaps_path, prob_path = create_directories(
        dirpath, outfile
    )

    genbankid, f_score = del_create_output_files(prob_path)

    # This function starts a big loop over each chromosome found in the input file.
    prob, perhaps, other = calculate_scores(
        df_input_file,
        prob_path,
        genbankid,
        f_score,
        perhaps_path,
        sample_outpath
    )
    print("organisms in file:", len(set(df_input_file.index.values)))
    print("High coverage high scoring organisms:", prob)
    print("Low medium coverage medium scoring organisms:", perhaps)
    print("other organisms - no images generated:", other)
    print("########### File completed #######")
    f_score.close()
    print("Wochenende plot script completed")


if __name__ == "__main__":
    main()
