#!/usr/bin/python3

"""
A whole genome/metagenome analysis pipeline in Python3
Author: Tobias Scheithauer
Author: Dr. Colin Davenport
Author: Fabian Friedrich

TODOs:
- Also trial with far more comprehensive bbmap adapters
- handle TruSeq and NEB Next adapters
  and test this vs alternatives to Trimmomatic, eg

Changelog
1.6 add ngmlr aligner, --longreads now omits Picard remove_dups by default (fails)
1.5.1 improve SOLiD adapter removal with fastp - configure var adapter_fastp
1.5 restructure wochenende_reporting, requires Python3.6+
1.4 add wochenende_plot.py file plotting
1.3 add samtools flagstat to get per cent reads aligned
1.3 add --testWochenende tests to test pipeline functionality and report success on a small reference
1.2 add --no-prinseq option to cancel prinseq read exclusion
1.2 add generation of unaligned reads (function runGetUnmappedReads) in FASTQ format
1.0 add reporting

"""

import sys
import os
import subprocess
import shutil
import argparse
import time

version = "1.6 - Mar 2020"

##############################
# CONFIGURATION
##############################

## Paths to commands - please edit as appropriate. If it is in your PATH just type the command. We recommend conda.
path_fastqc      = 'fastqc'
path_afterqc     = '/mnt/ngsnfs/tools/afterQC/AfterQC-0.9.6/after.py'
path_fastp       = 'fastp'
path_prinseq     = 'prinseq-lite.pl'
path_perl        = 'perl'
path_perldup     = '/mnt/ngsnfs/tools/Wochenende/dependencies/remove_pcr_duplicates.pl'
path_fastuniq    = 'fastuniq'
path_trimmomatic = 'trimmomatic'
path_fastq_mcf   = 'fastq_mcf'
path_bwa         = 'bwa'
path_samtools    = '/mnt/ngsnfs/tools/miniconda3/envs/wochenende/bin/samtools'
path_bamtools    = '/mnt/ngsnfs/tools/miniconda3/envs/wochenende/bin/bamtools'
path_sambamba    = 'sambamba'
path_java        = 'java'
path_abra_jar    = '/mnt/ngsnfs/tools/abra2/abra2_latest.jar'
path_minimap2    = 'minimap2'
path_ngmlr       = 'ngmlr'
## Paths to reference seqs. Edit as appropriate!
path_refseq_dict = {
    "2016_06_1p_genus" : "/working2/tuem/metagen/refs/2016/bwa/2016_06_PPKC_metagenome_test_1p_genus.fa",
    "2016_06_1p_spec_corrected" : "/lager2/rcug/seqres/metagenref/bwa/2016_06_PPKC_metagenome_test_1p_spec_change_cln.fa",
    "2016_06_1p_spec" : "/working2/tuem/metagen/refs/2016/bwa/2016_06_PPKC_metagenome_test_1p_spec_change.fa",
    "2019_01_meta" : "/lager2/rcug/seqres/metagenref/bwa/all_kingdoms_refseq_2019_Jan_final.fasta",
    "2019_10_meta_human" : "/lager2/rcug/seqres/metagenref/bwa/refSeqs_allKingdoms_201910_3.fasta",
    "2019_01_meta_mouse" : "/lager2/rcug/seqres/metagenref/bwa/all_kingdoms_refseq_2019_Jan_final_mm10_no_human.fasta",
    "2019_01_meta_mouse_ASF_OMM" : "/lager2/rcug/seqres/metagenref/bwa/mm10_plus_ASF_OMM.fasta",
    "2019_01_meta_mouse_ASF" : "/lager2/rcug/seqres/metagenref/bwa/mm10_plus_ASF.fasta",
    "2019_01_meta_mouse_OMM" : "/lager2/rcug/seqres/metagenref/bwa/mm10_plus_OMM.fasta",
    "hg19": "/lager2/rcug/seqres/HS/bwa/hg19.fa",
    "GRCh37": "/lager2/rcug/seqres/HS/bwa/GRCh37.fa",
    "GRCh38-45GB": "/lager2/rcug/seqres/HS/bwa/Homo_sapiens.GRCh38.dna.toplevel.fa",
    "GRCh38-noalt": "/lager2/rcug/seqres/HS/bwa/GRCh38_no_alt.fa",
    "GRCh38-mito" : "/lager2/rcug/seqres/HS/bwa/Homo_sapiens.GRCh38.dna.chromosome.MT.fa",
    "mm10": "/lager2/rcug/seqres/MM/bwa/mm10.fa",
    "rn6": "/lager2/rcug/seqres/RN/bwa/Rattus_norvegicus.Rnor_6.0.dna.toplevel.fa",
    "rat_1AR1_ont": "/lager2/rcug/seqres/RN/bwa/1AR1_2019_ONT_final.fasta",
    "zf10": "/lager2/rcug/seqres/DR/bwa/GRCz10.fa",
    "ss11": "/lager2/rcug/seqres/SS/bwa/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa",
    "PA14": "/lager2/rcug/seqres/PA/bwa/NC_008463.fna",
    "nci_viruses": "/lager2/rcug/seqres/metagenref/bwa/nci_viruses.fa",
    "testdb": "testdb/ref.fa",
}
# Adapters - edit as appropriate
ea_adapter_fasta = '/lager2/rcug/seqres/contaminants/2020_02/adapters/adapters.fa'
adapter_fasta = '/mnt/ngsnfs/tools/miniconda3/envs/wochenende/share/trimmomatic-0.38-0/adapters/TruSeq3-PE.fa'
adapter_fastp = '/lager2/rcug/seqres/contaminants/2020_02/adapters/adapters_solid.fa'
## Other
path_tmpdir = '/ngsssd1/rcug/tmp/'

##############################
# INITIALIZATION AND ORGANIZATIONAL FUNCTIONS
##############################

stage_outfile = ""
stage_infile = ""
fileList = []
global IOthreadsConstant
IOthreadsConstant='8'
global args
os.makedirs(path_tmpdir, exist_ok=True)

print('Wochenende - Whole Genome/Metagenome Sequencing Alignment Pipeline')
print('Wochenende was created by Dr. Colin Davenport and Tobias Scheithauer')
print('version: ' + version)
print()


def check_arguments(args):
    # Check argument cobination
    if(args.aligner == "minimap2" and not args.longread):
        args.longrad = True
        print("WARNING: Usage of minimap2 optimized for ONT data only. Added --longread flag.")

    if(args.readType == "PE" and args.aligner == "minimap2"):
        print("ERROR: Usage of minimap2 optimized for ONT data only. Combination of '--readType PE' and '--aligner minimap2' is not allowed.")
        sys.exit(1)

    if(args.readType == "PE" and args.longread):
        print("ERROR: Combination of '--readType PE' and '--longread' is not allowed.")
        sys.exit(1)

    if(args.fastp and args.aligner == "minimap2"):
        print("ERROR: Combination of '--fastp' and '--aligner minimap2' is not allowed.")
        sys.exit(1)

    if(args.fastp and args.longread):
        print("ERROR: Combination of '--fastp' and '--longread' is not allowed.")
        sys.exit(1)
    return args


def createProgressFile(args):
    # Read or create progress file
    with open(progress_file, mode='a+') as f:
        f.seek(0)
        progress = f.readlines()
    if progress == [] or progress[1].replace("\n", "") == "<current file>" or args.force_restart:
        with open(progress_file, mode='w') as f:
            f.writelines(["# PROGRESS FILE FOR Wochenende\n", "<current file>\n"])
        return None
    else:
        print("Found progress file x.tmp, attempting to resume after last completed stage. If not desired, use --force_restart or delete the .tmp progress files.")
        return progress[1].replace("\n", "")


def addToProgress(func_name, c_file):
    # Add run functions to progress file
    with open(progress_file, mode='r') as f:
        progress_lines = f.readlines()
        progress_lines[1] = c_file + "\n"
        if func_name+"\n" not in progress_lines:
            progress_lines.append(func_name + "\n")
    with open(progress_file, mode='w') as f:
        f.writelines(progress_lines)
    return progress_lines[1].replace("\n", "")


def runFunc(func_name, func, cF, newCurrentFile, *extraArgs):
    # Run function and add it to the progress file
    with open(progress_file, mode='r') as f:
        done = func_name in ''.join(f.readlines())
    if not done:
        if newCurrentFile:
            cF = func(cF, *extraArgs)
        else:
            func(cF, *extraArgs)
    return addToProgress(func_name, cF)


def runStage(stage, programCommand):
    # Run a stage of this Pipeline
    print("######  "+ stage + "  ######")
    try:
        #print(programCommand)
        process = subprocess.Popen(programCommand, stdout=subprocess.PIPE)
        output, error = process.communicate()
    except OSError as e:
        print(programCommand)
        print("Execution failed:", e, file=sys.stderr)
        sys.exit(1)


def deriveRead2Name(seRead):
    # Get name for paired end read based on single end
    read1 = seRead
    if "fastq" in seRead:
        if "_R1" in seRead:
            read2 = seRead.replace("_R1","_R2")
        else:
            print("seRead: "+str(seRead))
            raise NameError("Invalid format for Paired-End-Reads 1")
    elif "fq" in seRead:
        if "_R1" in seRead:
            read2 = seRead.replace("_R1","_R2")
        else:
            raise NameError("Invalid format for Paired-End-Reads 2")
    else:
        raise NameError("Invalid format for Paired-End-Reads 3")
    print("Read1 was: " + read1 + ", read2 derived as: " + read2)
    return read2


def rejigFiles(stage, stage_infile, stage_outfile):
    # Record, then prepare files for next stage
    fileList.append(stage)
    fileList.append(stage_infile)
    fileList.append(stage_outfile)
    stage_infile = stage_outfile


##############################
# TOOL FUNCTIONS
##############################

def runFastQC(stage_infile):
    # quality control
    stage = "FastQC"
    fastqc_out = stage_infile + '_fastqc_out'
    try:
        os.mkdir(fastqc_out)
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)
    fastQCcmd = [path_fastqc, '-t','4','-quiet', '-o', fastqc_out, stage_infile]
    runStage(stage, fastQCcmd)
    stage_outfile = stage_infile
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def runAfterQC(stage_infile):
    # automatic filtering, trimming, error removing and quality control
    stage = "AfterQC"
    print("######  "+ stage + "  ######")
    afterQCcmd = ['python', path_afterqc]
    runStage(stage, afterQCcmd)
    stage_outfile = ""
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def runFastpSE(stage_infile, noThreads):
    # all-in-one FASTQ-preprocessor - single end reads
    stage = "fastp - SE"
    print("######  "+ stage + "  ######")
    prefix = stage_infile.replace(".fastq","")
    stage_outfile = prefix + '.trm.fastq'
    fastpcmd = [path_fastp, '--in1=' + stage_infile, '--out1=' + stage_outfile,
                '--disable_quality_filtering', '--disable_length_filtering',
                '--length_required=40',
                '--adapter_fasta=' + adapter_fastp,
                '--cut_front', '--cut_window_size=5',
                '--cut_mean_quality=15', '--html='+ prefix + '.html',
                '--json='+ prefix + '.json', '--thread='+ noThreads]
    runStage(stage, fastpcmd)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def runFastpPE(stage_infile_1, stage_infile_2, noThreads):
    # all-in-one FASTQ-preprocessor - paired end reads
    stage = "fastp - PE"
    print("######  "+ stage + "  ######")
    prefix = stage_infile_1.replace(".fastq","")
    stage_outfile = prefix + '.fastp.fastq'
    fastpcmd = [path_fastp, '--in1='+ stage_infile_1, '--out1='+ stage_outfile,
                '--in2='+ stage_infile_2, '--out2='+ deriveRead2Name(stage_outfile),
                '--disable_quality_filtering', '--disable_length_filtering',
                '--length_required=40',
                '--adapter_fasta=' + adapter_fastp,
                '--cut_front', '--cut_window_size=5',
                '--cut_mean_quality=15', '--html='+ prefix + '.html',
                '--json='+ prefix + '.json', '--thread='+ noThreads]
    runStage(stage, fastpcmd)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def runPrinseq(stage_infile):
    # low complexity reads removal - single end reads
    stage = "Remove low complexity reads with Prinseq"
    stage_outfile = stage_infile
    prefix = stage_outfile.replace(".fastq","")
    stage_outfile = prefix + '.lc.fastq'
    prinseqCmd = [path_prinseq, '-fastq', stage_infile, '-lc_method', 'dust',
                  '-lc_threshold', '3', '-out_good', stage_outfile,
                  '-out_bad', prefix + '.lc_seqs.fq']
    runStage(stage, prinseqCmd)
    # prinseq adds extra .fastq by itself. Remove this by moving the file to
    # the filename expected by downstream apps
    prinseqOutfile= stage_outfile + ".fastq"
    try:
        shutil.move(prinseqOutfile, stage_outfile)
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def runPrinseqPE(stage_infile_1, stage_infile_2):
    # low complexity reads removal - paired end reads
    stage = "Remove low complexity reads with Prinseq"
    stage_outfile = stage_infile_1
    prefix = stage_outfile.replace(".fastq", "")
    stage_outfile = prefix + '.lc.fastq'
    prinseqCmd = [path_prinseq, '-fastq', stage_infile_1,
                  '-fastq2', stage_infile_2, '-lc_method', 'dust',
                  '-lc_threshold', '3', '-out_good', stage_outfile,
                  '-out_bad', prefix + '.lc_seqs.fq']
    runStage(stage, prinseqCmd)
    # prinseq adds extra .fastq by itself. Remove this by moving the file to
    # the filename expected by downstream apps
    try:
        shutil.move(stage_outfile + "_1.fastq", stage_outfile)
        shutil.move(stage_outfile + "_2.fastq", deriveRead2Name(stage_outfile))
        os.remove(stage_outfile + "_1_singletons.fastq")
        os.remove(stage_outfile + "_2_singletons.fastq")
    except OSError as e:
        print("Execution failed:", e, file=sys.stderr)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def runPerlDup(stage_infile):
    # duplicate reads removal
    stage = "Remove non-unique reads with Perldup"
    stage_outfile = stage_infile
    prefix = stage_outfile.replace(".fastq", '')
    stage_outfile = prefix + ".ndp.fastq"
    runPerlDupCmd = [path_perl, path_perldup, stage_infile, stage_outfile]
    runStage(stage, runPerlDupCmd)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def runFastUniq(stage_infile):
    # duplicate reads removal
    stage = "Remove non-unique reads with FastUniq"
    stage_outfile = stage_infile
    prefix = stage_outfile.replace(".fastq", '')
    stage_outfile = prefix + ".ndp.fastq"
    with open('readlist.tmp', 'a') as readlist:
        readlist.write(stage_infile.replace(os.getcwd()+'/', '') + '\n')
        #readlist.write('\n')
        readlist.write(deriveRead2Name(stage_infile).replace(os.getcwd()+'/', '') + '\n')
    runFastuniqCmd = [path_fastuniq, '-i', 'readlist.tmp', '-t', 'q',
                      '-o', stage_outfile, '-p', deriveRead2Name(stage_outfile),
                      '-c', '0']
    runStage(stage, runFastuniqCmd)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def runTMTrimming(stage_infile):
    # adapter and quality trimming - single end
    stage = "Trimming with Trimmomatic - SE"
    stage_outfile = stage_infile
    prefix = stage_infile.replace(".fastq","")
    stage_outfile = prefix + '.trm.fastq'
    # Also trial with far more comprehensive bbmap adapters:
    # /mnt/ngsnfs/tools/miniconda2/pkgs/bbmap-37.17-0/opt/bbmap-37.17/resources/adapters.fa
    trimCmd = [path_trimmomatic, 'SE', '-threads', IOthreadsConstant, '-phred33',
               stage_infile, stage_outfile,
               'ILLUMINACLIP:'+adapter_fasta+':2:30:10', 'LEADING:3',
               'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:36']
    runStage(stage, trimCmd)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def runTMTrimmingPE(stage_infile):
    # adapter and quality trimming - paired end
    stage = "Trimming with Trimmomatic - PE"
    stage_outfile = stage_infile
    prefix = stage_infile.replace(".fastq","")
    stage_outfile = prefix + '.trm.fastq'
    tmpfile1 = prefix + '1.tmp'
    tmpfile2 = prefix + '2.tmp'
    trimCmd = [path_trimmomatic, 'PE', '-threads', IOthreadsConstant, '-phred33',
               stage_infile, deriveRead2Name(stage_infile), stage_outfile,
               tmpfile1, deriveRead2Name(stage_outfile), tmpfile2,
               'ILLUMINACLIP:'+adapter_fasta+':2:30:10', 'LEADING:3',
               'TRAILING:3', 'SLIDINGWINDOW:4:15', 'MINLEN:36']
    runStage(stage, trimCmd)
    os.remove(tmpfile1)
    os.remove(tmpfile2)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile



def runEATrimming(stage_infile):
    # adapter trimming
    stage = "Trimming with EA-utils"
    prefix = stage_infile.replace(".fastq","")
    stage_outfile = prefix + '.tre.fastq'
    trimCmd = [path_fastq_mcf, '-f', '-o', stage_outfile, ea_adqapter_fasta,
               stage_infile]
    runStage(stage, trimCmd)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def runAligner(stage_infile, aligner, index, noThreads, readType):
    # Alignment - Short-read single and paired end using bwa-mem. minimap2 or ngmlr for long reads

    ngmlrMinIdentity = 0.85 # Aligner ngmlr only: minimum identity (fraction) of read to reference

    stage= "Alignment"
    print("######  "+ stage + "  ######")

    prefix = stage_infile.replace(".fastq","")
    stage_outfile = prefix + '.bam'
    global inputFastq
    readGroup = os.path.basename(inputFastq.replace(".fastq",""))

    alignerCmd = ""
    if "minimap2" in aligner:
        alignerCmd = [path_minimap2, '-x', 'map-ont', '-a', '-t', str(noThreads), str(index), stage_infile]
    elif "ngmlr" in aligner:
        alignerCmd = [path_ngmlr, '-x', 'ont', '-i', str(ngmlrMinIdentity), '-t', str(noThreads), '-r', str(index), '-q', stage_infile]
    elif "PE" in readType:
        stage_infile2 = deriveRead2Name(stage_infile)
        alignerCmd = [path_bwa, 'mem', '-t', str(noThreads), '-R',
                  '"@RG\\tID:' + readGroup + '_001\\tSM:' + readGroup + '"',
                  str(index), stage_infile, stage_infile2]
    elif "SE" in readType:
        alignerCmd = [path_bwa, 'mem', '-t', str(noThreads), '-R',
                  '"@RG\\tID:' + readGroup + '_001\\tSM:' + readGroup + '"',
                  str(index), stage_infile]
    else:
        print("Read type not defined")
        system.exit(1)

    samtoolsCmd = ['|', path_samtools, 'view', '-@', IOthreadsConstant, '-bhS',
                   '>', stage_outfile]
    wholeCmd = alignerCmd + samtoolsCmd
    print(' '.join(wholeCmd))
    wholeCmdString=' '.join(wholeCmd)

    try:
        # could not get subprocess.run, .call etc to work with pipes and redirect '>'
        os.system(wholeCmdString)
    except:
        sys.exit(1)

    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def runBAMsort(stage_infile):
    # Ministage runStage(stage, bwaCmd)
    stage = "Sort BAM"
    prefix = stage_infile.replace(".bam","")
    stage_outfile = prefix + '.s.bam'
    samtoolsSortCmd = [path_samtools, 'sort', '-@', IOthreadsConstant,
                       stage_infile, '-o', stage_outfile]
    runStage(stage, samtoolsSortCmd)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def runBAMindex(stage_infile):
    # Stage output not used further in flow
    stage = "Index BAM"
    samtoolsIndexCmd = [path_samtools, 'index', stage_infile]
    runStage(stage, samtoolsIndexCmd)
    # No rejigfiles needed as dead end
    return 0

def runSamtoolsFlagstat(stage_infile):
    # Stage output not used further in flow
    stage = "samtools flagstat"
    print("######  "+ stage + "  ######")
    flagstatOutput = [stage_infile, '.flagstat.txt']
    flagstatOutputStr = ""
    flagstatOutputStr = ''.join(flagstatOutput)
    #print(flagstatOutputStr)
    samtoolsFlagstatCmd = [path_samtools, 'flagstat', stage_infile, '>', flagstatOutputStr]
    samtoolsFlagstatCmdStr = ""
    samtoolsFlagstatCmdStr = ' '.join(samtoolsFlagstatCmd)
    print(samtoolsFlagstatCmdStr)

    try:
        # could not get subprocess.run, .call etc to work with pipes and redirect '>'
        os.system(samtoolsFlagstatCmdStr)
    except:
        print ("########################################")
        print ("##### Error with flagstat")
        print ("########################################")
        sys.exit(1)
    # No rejigfiles needed as dead end
    return 0

def runGetUnmappedReads(stage_infile):
    # Stage output not used further in flow
    stage = "Get unmapped reads from BAM"
    print("######  "+ stage + "  ######")
    unmappedFastq = [stage_infile, '.unmapped.fastq']
    unmappedFastqString = ''
    unmappedFastqString =''.join(unmappedFastq)
    samtoolsGetUnmappedCmd = [path_samtools, 'view', '-f', '4', stage_infile, '|', path_samtools, 'bam2fq', '-', '>', unmappedFastqString ]
    print(' '.join(samtoolsGetUnmappedCmd))
    samtoolsGetUnmappedCmdString=' '.join(samtoolsGetUnmappedCmd)
    print(samtoolsGetUnmappedCmdString)

    try:
        # could not get subprocess.run, .call etc to work with pipes and redirect '>'
        os.system(samtoolsGetUnmappedCmdString)
    except:
        print ("########################################")
        print ("##### Error getting unmapped reads with samtools")
        print ("########################################")
        sys.exit(1)

    # No rejigfiles needed as dead end
    return 0


def runMQ30(stage_infile):
    # Remove reads with less than MQ30
    stage = "Remove MQ30 reads"
    prefix = stage_infile.replace(".bam","")
    stage_outfile = prefix + '.mq30.bam'
    samtoolsMQ30Cmd = [path_samtools, 'view','-@', IOthreadsConstant, '-b', '-q', '30', stage_infile, '-o', stage_outfile]
    runStage(stage, samtoolsMQ30Cmd)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile

def runBamtools(stage_infile):
    # Keep only reads with 0 or 1 mismatch
    stage = "Keep only reads with 0 or 1 mismatch - for metagenomics"
    prefix = stage_infile.replace(".bam","")
    stage_outfile = prefix + '.01mm.bam'
    tmpfile0 = prefix + '.0mm.tmp.bam'
    tmpfile1 = prefix + '.1mm.tmp.bam'
    keep0mm = [path_bamtools, 'filter', '-in', stage_infile, '-out', tmpfile0, '-tag', 'NM:0']
    keep1mm = [path_bamtools, 'filter', '-in', stage_infile, '-out', tmpfile1, '-tag', 'NM:1']
    bam_merge = [path_samtools, 'merge', '-@', IOthreadsConstant, stage_outfile, tmpfile0, tmpfile1]
    bamtools_cmd1 = ' ' .join(keep0mm)
    bamtools_cmd2 = ' ' .join(keep1mm)
    merge = ' ' .join(bam_merge)


    try:
        # could not get subprocess.run, .call etc to work with "&&"
        #print(bamtools_cmd)
        #os.system("path_bamtools filter -in stage_infile -out tmpfile0 -tag NM:0 && keep1mm = path_bamtools filter -in stage_infile -out tmpfile1 -tag NM:1 && bam_merge = path_samtools merge -@ IOthreadsConstant stage_outfile tmpfile0 tmpfile1")
        os.system(bamtools_cmd1)
        os.system(bamtools_cmd2)
        os.system(merge)
        os.remove(tmpfile0)
        os.remove(tmpfile1)

    except:
        print ("########################################")
        print ("############ Error with BAMtools commands")
        print ("########################################")
        sys.exit(1)

    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def markDups(stage_infile):
    # duplicate removal in bam
    stage = "Sambamba mark duplicates"
    prefix = stage_infile.replace(".bam","")
    stage_outfile = prefix + '.dup.bam'
    markDupsCmd = [path_sambamba, 'markdup', '--remove-duplicates',
                   '-t', IOthreadsConstant, '--sort-buffer-size=4096',
                   '--hash-table-size=512288', '--overflow-list-size=200000',
                   '--tmpdir=tmp', stage_infile, stage_outfile]
    runStage(stage, markDupsCmd)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def runIDXstats(stage_infile):
    # simple alignment statistics
    stage = "Samtools index stats - chromosome counts"
    print("######  "+ stage + "  ######")

    prefix = stage_infile.replace(".bam",".bam.txt")
    stage_outfile = prefix

    samtoolsidxCmd = [path_samtools, 'idxstats', stage_infile, '>', stage_outfile]
    wholeCmdString=' '.join(samtoolsidxCmd)
    try:
        # could not get subprocess.run, .call etc to work with pipes and redirect '>'
        os.system(wholeCmdString)
    except:
        sys.exit(1)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def calmd(stage_infile, fasta):
    # MD tag generation
    stage = "Samtools calmd"
    prefix = stage_infile.replace(".bam","")
    stage_outfile = prefix + '.calmd.bam'
    #samtools calmd -@ 8 -b $filename $ref > $outfile
    calmdCmd = [path_samtools, 'calmd', '-Q', '-@', IOthreadsConstant,
                '-b', stage_infile, fasta, '>', stage_outfile]
    wholeCmdString=' '.join(calmdCmd)
    try:
        # could not get subprocess.run, .call etc to work with pipes and redirect '>'
        os.system(wholeCmdString)
    except:
        sys.exit(1)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile


def abra(stage_infile, fasta, threads):
    stage = "Abra read realignment"
    prefix = stage_infile.replace(".bam","")
    stage_outfile = prefix + '.abra.bam'
    #java -Xmx16G -jar /mnt/ngsnfs/tools/abra2/abra2_latest.jar --in $bam --out $bam.abra.bam --ref $ref --threads 14 --dist 1000 --tmpdir /data/tmp/ > abra.log
    abra_tmpdir = os.path.join(path_tmpdir, 'abra_' +  str(int(time.time())))
    os.makedirs(abra_tmpdir, exist_ok=True)
    abraCmd = [path_java, '-Xmx16G', '-jar', path_abra_jar, '--in', stage_infile,
               '--out', stage_outfile, '--ref', fasta, '--threads', threads,
               '--dist', '1000', '--tmpdir', abra_tmpdir]
    runStage(stage, abraCmd)
    rejigFiles(stage, stage_infile, stage_outfile)
    return stage_outfile

def runTests(stage_infile):
    stage = "Running internal tests"
    # use sbatch script to start tests
    # This section should check the output
    # Test output will be printed on std out
    print("\n\n")
    print("####################################################################")
    print("######################### Starting tests ###########################")
    print("####################################################################")

    failedCount=0

    print("\n\nTest tempfile length")
    tempFile=""
    tempFile="testdb/reads_R1.fastqprogress.tmp"
    try:
        with open(tempFile, mode='r') as f:
            f.seek(0)
            testList = f.readlines()
    except:
        print("Test FAILED - could not open !" + str(tempFile))

    print(len(testList))
    if (len(testList) == 18):
        print("Test tempfile length 18 lines  ...  passed!")
    else:
        print("Test FAILED!")
        failedCount=failedCount+1


    print("\n\nTest unmapped file length")
    tempFile=""
    tempFile="testdb/reads_R1.ndp.lc.trm.s.bam.unmapped.fastq"
    try:
        with open(tempFile, mode='r') as f:
            f.seek(0)
            testList = f.readlines()
            print(len(testList))
    except:
        print("Test FAILED, could not open file!")

    if (len(testList) == 44):
        print("Test unmapped length 44 lines ...  passed!")
    else:
        print("Test FAILED!")
        failedCount=failedCount+1


    print("\n\nTest bam.txt file contents")
    tempFile=""
    tempFile="testdb/reads_R1.ndp.lc.trm.s.mq30.01mm.dup.bam.txt"
    with open(tempFile, mode='r') as f:
        f.seek(0)
        testList = f.readlines()
    print(testList[0])
    if (testList[0] == "1	599940	411	0\n"):
        print("Test stats contents  ...  passed!")
    else:
        print("Test FAILED!")
        failedCount=failedCount+1

    print(str(failedCount) + " tests failed!\n")

    print("\nTesting completed, cleaning up testdb directory")
    os.system("bash wochenende_test_cleanup.sh")




##############################
# MAIN FUNCTION (PIPELINE DEFINITION)
##############################


def main(args, sys_argv):
    args = check_arguments(args)
    global progress_file
    progress_file = args.fastq + "progress.tmp"
    currentFile = createProgressFile(args)
    threads = args.threads
    global inputFastq
    inputFastq = args.fastq
    if currentFile is None:
        currentFile = inputFastq

    print ("Meta/genome selected: " + args.metagenome)

    # Single ended input reads
    if args.readType == "SE":
        if not args.longread and not args.no_fastqc:
            currentFile = runFunc("runFastQC", runFastQC, currentFile, False)
        if not args.longread and not args.no_duplicate_removal:
            currentFile = runFunc("runPerlDup", runPerlDup, currentFile, True)
        if not args.longread and not args.no_prinseq:
            currentFile = runFunc("runPrinseq", runPrinseq, currentFile, True)
        if args.fastp and not args.longread:
            currentFile= runFunc("runFastpSE", runFastpSE, currentFile, True, args.threads)
        if not args.longread and not args.fastp:
            currentFile = runFunc("runTMTrimming", runTMTrimming, currentFile, True)
        #if not args.longread:
            #currentFile = runFunc("runEATrimming", runEATrimming, currentFile, True)
        currentFile = runFunc("runAligner", runAligner, currentFile, True,
                                 args.aligner, path_refseq_dict.get(args.metagenome),
                                 args.threads, args.readType)
        currentFile = runFunc("runBAMsort", runBAMsort, currentFile, True)
        currentFile = runFunc("runBAMindex", runBAMindex, currentFile, False)
        currentFile = runFunc("runSamtoolsFlagstat", runSamtoolsFlagstat, currentFile, False)
        currentFile = runFunc("runGetUnmappedReads", runGetUnmappedReads, currentFile, False)
        if args.mq30:
            currentFile = runFunc("runMQ30", runMQ30, currentFile, True)
        if args.remove_mismatching:
            currentFile = runFunc("runBamtools", runBamtools, currentFile, True)
        if not args.no_duplicate_removal and not args.longread:
            currentFile = runFunc("markDups", markDups, currentFile, True)
        currentFile = runFunc("runIDXstats", runIDXstats, currentFile, False)
        if not args.no_abra and not args.longread:
            currentFile = runFunc("abra", abra, currentFile, True,
                                 path_refseq_dict.get(args.metagenome), threads)
        currentFile = runFunc("calmd", calmd, currentFile, True,
                                 path_refseq_dict.get(args.metagenome))
        currentFile = runFunc("runBAMindex2", runBAMindex, currentFile, False)
        currentFile = runFunc("runIDXstats", runIDXstats, currentFile, False)



    # Paired end input reads. Long reads cannot be paired end (true 2020).
    elif args.readType == "PE":
        print('Input File 1 : ' + currentFile)
        print('Input File 2 : ' + deriveRead2Name(currentFile))
        if not args.longread and not args.no_fastqc:
             runFunc("runFastQC1", runFastQC, currentFile, False)
             runFunc("runFastQC2", runFastQC, deriveRead2Name(currentFile), False)
        # PerlDup does not work for PE data
        # read1 = runFunc("runPerlDup1", runPerlDup, read1, True)
        # read2 = runFunc("runPerlDup2", runPerlDup, read2, True)
        # currentFile = runFunc("runFastUniq", runFastUniq, currentFile, True)
        # Prinseq not tested for PE data
        # currentFile = runFunc("runPrinseq", runPrinseqPE, currentFile, True, deriveRead2Name(currentFile))
        if args.fastp:
            currentFile = runFunc("runFastpPE", runFastpPE, currentFile, True,
            deriveRead2Name(currentFile), args.threads)
        currentFile = runFunc("runTMTrimmingPE", runTMTrimmingPE, currentFile, True)
        # currentFile = runFunc("runEATrimming", runEATrimming, currentFile, True)
        currentFile = runFunc("runAligner", runAligner, currentFile, True,
                                 args.aligner, path_refseq_dict.get(args.metagenome),
                                 args.threads, args.readType)
        currentFile = runFunc("runBAMsort", runBAMsort, currentFile, True)
        currentFile = runFunc("runBAMindex", runBAMindex, currentFile, False)
        currentFile = runFunc("runSamtoolsFlagstat", runSamtoolsFlagstat, currentFile, False)
        currentFile = runFunc("runGetUnmappedReads", runGetUnmappedReads, currentFile, False)
        if args.mq30:
            currentFile = runFunc("runMQ30", runMQ30, currentFile, True)
        if args.remove_mismatching:
            currentFile = runFunc("runBamtools", runBamtools, currentFile, True)
        if not args.no_duplicate_removal:
            currentFile = runFunc("markDups", markDups, currentFile, True)
        currentFile = runFunc("runIDXstats", runIDXstats, currentFile, False)
        if not args.no_abra:
            currentFile = runFunc("abra", abra, currentFile, True,
                                 path_refseq_dict.get(args.metagenome), threads)
        currentFile = runFunc("calmd", calmd, currentFile, True,
                                 path_refseq_dict.get(args.metagenome))
        currentFile = runFunc("runBAMindex2", runBAMindex, currentFile, False)
        currentFile = runFunc("runIDXstats", runIDXstats, currentFile, False)

    else:
        print(" --readType must be set to either SE or PE (meaning single ended or paired-end)")




    # Report all files
    if args.debug:
        i=0
        for word in fileList:
            print("Filelist item: " + fileList[i])
            i = i + 1

    # Report all percentage mapped


    if args.testWochenende:

        # Run internal tests
        currentFile = runFunc("runTests", runTests, currentFile, False)




##############################
# COMMAND LINE ARGUMENTS DEFINITION
##############################

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        epilog="We recommend using bioconda for the installation of the tools. Remember to run 'conda activate <environment name>' before you start if you are using bioconda. Details about the installation are available on https://github.com/MHH-RCUG/Wochenende#installation")

    parser.add_argument("fastq", help="_R1.fastq Input read1 fastq file",
                        type=lambda x: (os.path.abspath(os.path.expanduser(x))))

    parser.add_argument("--aligner", help="Aligner to use, either bwamem, ngmlr minimap2. Usage of minimap2 and ngmlr currently optimized for nanopore data only.",
                        action="store", choices=["bwamem","minimap2", "ngmlr"], default="bwamem")

    parser.add_argument("--readType", help="Single end or paired end data",
                        action="store", choices=["PE", "SE"], default="SE")

    parser.add_argument("--metagenome", help="Meta/genome reference to use",
                        action="store", choices=list(path_refseq_dict))

    parser.add_argument("--threads", help="Number of cores, default = 16",
                        action="store", default="16")

    parser.add_argument("--fastp", help="Use fastp instead of fastqc and trimmomatic", action="store_true")

    parser.add_argument("--debug", help="Report all files", action="store_true")

    parser.add_argument("--longread", help="Only do steps relevant for long PacBio/ONT reads eg. no trimming, alignment & bam conversion", action="store_true")

    parser.add_argument("--no_duplicate_removal", help="Skips steps for duplicate removal. Recommended for amplicon sequencing.", action="store_true")

    parser.add_argument("--no_prinseq", help="Skips prinseq step for low_complexity sequence removal.", action="store_true")

    parser.add_argument("--no_fastqc", help="Skips FastQC quality control step.", action="store_true")

    parser.add_argument("--no_abra", help="Skips steps for Abra realignment. Recommended for metagenome and amplicon analysis.", action="store_true")

    parser.add_argument("--mq30", help="Remove reads with mapping quality less than 30. Recommended for metagenome and amplicon analysis.", action="store_true")

    parser.add_argument("--remove_mismatching", help="Remove reads with 2 or more mismatches (via the NM bam tag)", action="store_true")

    parser.add_argument("--force_restart", help="Force restart, without regard to existing progress", action="store_true")

    parser.add_argument("--testWochenende", help="Run pipeline tests vs testdb, needs the subdirectory testdb, default false", action="store_true")

    if len(sys.argv) == 1:
        parser.print_help()
    else:
        main(parser.parse_args(), sys.argv[1:])
