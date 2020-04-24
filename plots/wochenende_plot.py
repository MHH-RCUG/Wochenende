#!/usr/bin/env python
# coding: utf-8


#Author: Konstantinos Sifakis
#Script made: Jan 2020-Feb 2020
#Script made for: MHH - RCUG
# Plot coverage data after the wochenende_posthoc_filter.sh sambamba coverage step of the Wochenende pipeline.
# Create some png files,one txt that includes the organisms that seem related
# and a tsv file with information from tsv file that was given.
# Usage: python3 wochenende_plot.py dup.calmd_cov_window.txt.filt.csv


#________CHECK what is needed:



import argparse
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import os
import math
import statistics


#Changelog

#0.2 reduce --minMeanCov default from 1.0 to 0.2 after testing, add version.
#0.1 initial working script


# use_line_collection should be True
import warnings
warnings.filterwarnings("ignore",category=UserWarning)


#Needed to read file  [ --filename ]  and then [ file directory ] so that file would open as a tsv file

parser = argparse.ArgumentParser(description='Plot coverage data after the wochenende_posthoc_filter.sh sambamba coverage step of the Wochenende pipeline. Creates png files of absolute and mean coverage per taxon, and a statistical overview file.')
parser.add_argument('filename1',default='check_string_for_empty',
                    help='Usage: python3 wochenende_plot.py dup.calmd_cov_window.txt.filt.csv')
parser.add_argument('--minMeanCov', default=1,
                     help='''The minimum mean coverage. Default: 1 
			This does not affect the score. After the score has been calculated, this number will be the threshold for both bad scored and good scored taxa.\n
			## Case1: Poor scoring that have a maximum score above this minMeanCov value -> Images go to folder potentially_present.\n
			## Case2: Poor scoring will have a maximum score below or equal to this minMeanCov value -> Images will not be printed (unless --createAllPngs is given).\n
			## Case3: Highly scored taxa that have a minimum score above this minMeanCov value -> Images will go to folder probably_present\n
			## Case4: Highly scored taxa that have a minimum score below or equal to this minMeanCov value -> Images will go to folder potentially present.\n				 
			''')
parser.add_argument('--createAllPngs', default=0,
                     help='Choose if you want to create figures for the poorly scored taxa which also have their max (mean_coverage) <= cov input argument(--minMeanCov)\nIf you wish to create the figures you should write --createAllPngs 1\n Default: 0, do not create poor scoring figures')
parser.add_argument('--sclim', default=0.0 , help= 'Choose the Score limit to filter whether a taxon is going to be well or poorly scored\nDefault: 0.0, negative scores are termed bad, positive values are highly scored.')
parser.add_argument('--minWindows', default=5 , help='Choose the least number of 100kbp windows per organism that should be taken into account.\n Generally, more then 5 is suggested, otherwise quartiles Q1 or Q3 will be useless. Since we are focussing on bacteria and bins are 100,000 bp length (50,000bp overlap), there should be no point on having less then 5 window bins accepted. Default: 5.')
args = parser.parse_args()


#_____________________________________________________


version = 0.2
print("\n###############################")
print("Wochenende plotting, Version: %s" % version)
print("###############################\n")

filename=args.filename1
BDSC=int(args.createAllPngs)
Scorelimitforminormax=float(args.minMeanCov)
Scorelimit=float(args.sclim)
minWindowsNeeded=int(args.minWindows)



#############################SCORE LIMITS############################

#AnotherScorelimit=2.0


defaultorgnames=2

#filename=str(args)
#filename=filename.split("(")[1].split(")")[0].split("=")[1].replace("'","") # this is the worst way to do it"


if filename == 'check_string_for_empty':
    print ('#### No argument was given')
else :
    print ('#### Filename entered correctly ####')

    # Shorten filename, otherwise it causes errors on windows systems
    # typical filename: KGCF14D_S2_R1.ndp.lc.trm.s.mq30.01mm.dup.bam.txt.rep.s.csv
    tmpfile=filename
    tmpfile=tmpfile.replace("_complete_genome","")
    tmpfile=tmpfile.replace(".ndp.lc","")
    tmpfile=tmpfile.replace(".trm.s","")
    tmpfile=tmpfile.replace(".mq30.01mm","")
    tmpfile=tmpfile.replace(".dup.bam","")
    tmpfile=tmpfile.replace(".rep.s.csv","")
    tmpfile=tmpfile.replace(".dup_cov_window","")
    tmpfile=tmpfile.replace(".txt.filt","")
    tmpfile=tmpfile.replace(".csv","")
    outfile=tmpfile

    #DO I HAVE A TSV FILE? #####IS IT A TSV?????############
    dirpath = os.getcwd()
    print("INFO: Current directory is : " + dirpath)
    print("INFO: Current file path is : " + dirpath +"/" + filename)
    #Panda_input_file=pd.read_csv("{}/{}".format(dirpath,filename),delimiter="\t")
    Panda_input_file=pd.read_csv("{}/{}".format(dirpath,filename),delimiter="\t")
    #fileitself=outfile.split("/")[-1]
    #FN.head(1)

    # Files summary
    print("INFO: Input filename: "+ str(filename))
    print("INFO: Output filename: %s" % outfile)
    #print("Output filename2: %s" % fileitself)


    #############MAKE NEW DIRECTORIES
    wocdirpath="{}/images".format(dirpath)
    newdirpath="{}/images/{}".format(dirpath,outfile)
    mightdirpath="{}/images/{}/perhaps_present".format(dirpath,outfile)
    newdirpathSC="{}/images/{}/prob_present".format(dirpath,outfile)


    try:
        os.mkdir(wocdirpath)
    except OSError:
        print ("Creation of the directory %s failed" % wocdirpath)
        print( '###__Maybe this directory already exists__###__OR you are not allowed to create it__###')
    else:
        print ("Successfully created the directory %s " % wocdirpath)

    try:
        os.mkdir(newdirpath)
    except OSError:
        print ("Creation of the directory %s failed" % newdirpath)
        print( '###__Maybe this directory already exists__###__OR you are not allowed to create it__###')
    else:
        print ("Successfully created the directory %s " % newdirpath)

    try:
        os.mkdir(newdirpathSC)
    except OSError:
        print ("Creation of the directory %s failed" % newdirpathSC)
        print( '###__Maybe this directory already exists__###__OR you are not allowed to create it__###')
    else:
        print ("Successfully created the directory %s " % newdirpathSC)

    try:
        os.mkdir(mightdirpath)
    except OSError:
        print ("Creation of the directory %s failed" % mightdirpath)
        print( '###__Maybe this directory already exists__###__OR you are not allowed to create it__###')
    else:
        print ("Successfully created the directory %s " % mightdirpath)

    RC="readCount"
    MC="meanCoverage"

    #set(Blut_1_S7_R1["# chrom"])
    B1=Panda_input_file.set_index('# chrom')

##########Delete .txt file that keeps the scores

    try:
        os.remove('{}/Scores.csv'.format(newdirpathSC))
    except:
        pass

    try:
        os.remove('{}/GenBankID.txt'.format(newdirpathSC))
    except:
        pass



#############################SCORE LIMITS############################

#    Scorelimit=0.0
#    AnotherScorelimit=2.0


#Make a new .txt file for Scores and Yes or No column

    f_score = open ('{}/Scores.csv'.format(newdirpathSC),'a' )
    f_score.write("GenBank_ID\tGenus\tSpecies\tStrain\tStrain_info\tScore\tScorelimit>={}\tMight_be_related_to_smthg".format(Scorelimit))
    GENBANKID = open ('{}/GENBANKID.txt'.format(newdirpathSC),'a')
#   GENBANKID.write("GenBank_ID\tGenus\tSpecies\tStrain\tStrain_info\tScore\tYes_or_No")



    #This section starts a big loop over each chromosome found in the file of input windows.

    for i in set(Panda_input_file["# chrom"]):


#########################################Read Counts################################
#__________________________________________________________________________________#

        small=B1.loc[[i]]
        B1maxread=max(small["readCount"])
        B1minread=min(small["readCount"])
        B1maxStart=max(small["chromStart"])
        B1minStart=min(small["chromStart"])


        # Shorten filename, otherwise it causes errors on windows systems
        # typical filename: KGCF14D_S2_R1.ndp.lc.trm.s.mq30.01mm.dup.bam.txt.rep.s.csv
        tmpfile=i
        tmpfile=tmpfile.replace("complete_genome","")
        tmpfile=tmpfile.replace("complete","")
        tmpfile=tmpfile.replace("sequence","")
        tmpfile=tmpfile.replace("chromosome","")
        tmpfile=tmpfile.replace(".ndp.lc","")
        tmpfile=tmpfile.replace(".trm.s","")
        tmpfile=tmpfile.replace(".mq30.01mm","")
        tmpfile=tmpfile.replace(".dup.bam","")
        tmpfile=tmpfile.replace(".rep.s.csv","")
        tmpfile=tmpfile.replace(".dup_cov_window","")
        tmpfile=tmpfile.replace(".txt.filt","")
        tmpfile=tmpfile.replace(".csv","")
        tmpfile=tmpfile.replace("__","_")
        i=tmpfile


        #########################################MEAN COVERAGE################################
        #____________________________________________________________________________________#

#______________        #Those are used mostly in order to create figures and their axis man and min.
        B1maxcove=max(small["meanCoverage"])                   #This one is for axis
        B1maxcove2=max(small["meanCoverage"])                  #This one is for real use.. (Max of this will show us where to put it, in which directory)

        B1mincove=min(small["meanCoverage"])
        Bin_size=(max(small["chromEnd"])-max(small["chromStart"]))/2 # /2 because every window is "covered" 2 times
        Bin_size_Multi_howmany=Bin_size*(len(small["meanCoverage"])+1)    #+1 Correction for overlaping windows
        Len_empty=max(small["chromEnd"])-(len(small["meanCoverage"])*Bin_size)
        Howmany_empty_Bins=Len_empty/Bin_size
        Howmany_empty_Bins_cor=Howmany_empty_Bins -1 # Correction for overlaping windows
#        print(Howmany_empty_Bins,"------",i)
#        print(Bin_size,"Bin_size",i)
#        print(Bin_size_Multi_howmany,"Bin_size_Multi_howmany")
#        print(Howmany_empty_Bins_cor,"Howmany Bins?")
#        print(max(small["chromEnd"]),"______________",i)
        j=B1maxcove
        if j>=0 and j<=1:
            k=1
        elif j>=1 and j<=2:
            k=2
        elif j>=1 and j<=5:
            k=5
        elif j>=5 and j<=10:
            k=10
        elif j/10>=1 and j<=20:
            k=20
        elif j/20>=1 and j<=50:
            k=50
        elif j>50>=1 and j<=100:
            k=100
        else:
            k=200
        B1maxcove=k
#Here i am doing tha scaling for axis_ not the best way to do it tho............

        Score=statistics.mean(small["meanCoverage"])

# We will name the average mean.coverage of our samples =Score and this will help us remember that (average) should ( ~ ) the final (Score)

        #print(Score,"_______________________" ,i)
#Checking if mean is more then 1.?? and if small dataframe has more then 2 inputs..
        if len(small["meanCoverage"])>10:
           # Dividedby=statistics.stdev(small["meanCoverage"])

            #Score=((Score)**3)/Dividedby                       ###########################################################################
#NEW SCORE:

            lst=["What?"]
            for meanCove in small["meanCoverage"]:
                Normalized_SD=(meanCove-B1mincove)/(B1maxcove-B1mincove)
                lst.append(Normalized_SD)
            Norm_SD=statistics.mean(lst[1:])

            Dividedby=Norm_SD
            Score2=((Score*len(small["meanCoverage"]))-(Howmany_empty_Bins_cor*Score))/((Dividedby+1)*len(small["meanCoverage"])) 
            #/Dividedby should not be 0 ever + if it's less then 1 then it automatically becomes a sort of multiplication
            # So we should use Divided by +1 in order to be sure that this will be working as it supposed to be working... Reusing our Score every time depending on
            # how much the NORMALIZED Standard Deviation is. +1 in order to really have a good result serving our goal to reduce Score

            # we also divide by the len(small["meanCoverage"]) which is basically our data so we divide by the #number of the data so that we manipulate Score,
            # and it should not get too high if there are too many but low mean Coverage windows.

            Score=Score2            #something could be done for more score calculations

        elif len(small["meanCoverage"])>=2 and len(small["meanCoverage"])<=10:
            lst=["What?"]
            for meanCove in small["meanCoverage"]:
                Normalized_SD=(meanCove-B1mincove)/(B1maxcove-B1mincove)
                lst.append(Normalized_SD)
            Norm_SD=statistics.mean(lst[1:])

            Dividedby=Norm_SD
            Score2=(((Score*len(small["meanCoverage"]))-(Howmany_empty_Bins_cor*Score))/((Dividedby+1)*len(small["meanCoverage"]))) -(2*B1mincove) #penalty 
            # for having less then 10 windows Score - 2 times the min of those windows.
            Score=Score2            # MORE to be done? maybe ?.
            #NEW SCORE
            # print(Dividedby,"__Divided_by___",i)
            # print(Score,"SCORED2")
            # That should happen only if not good mean and not good length of the small data.frame
        else:
            Score=-100 #  Typical not good score
            # GIVE number of windows that are above mean-2 sd and those that are above Q1 of the small data.frame
            # minWindowsNeeded=5       # more then 5 is suggested otherwise Q1 and Q3 are useless.....


        mean_describe=1
        std_describe=2
        min_describe=3
        Q1_describe=4             #Those are more or less describing a panda.describe() function input
        Q2_describe=5             #Not used
        Q3_describe=6
        max_describe=7

        Mean_2SD=small.describe()["meanCoverage"][mean_describe]-(2*small.describe()["meanCoverage"][std_describe])

        if small["meanCoverage"].shape[0]>=minWindowsNeeded and Mean_2SD>=0:
            Q1=small.describe()["meanCoverage"][Q1_describe]
            passedQ1F=small[small["meanCoverage"]>=Q1].shape[0]
            passedQ1=(passedQ1F*100)/small["meanCoverage"].shape[0]
            #Not needed , we have calculated it before if.
            #Mean_2SD=small.describe()["meanCoverage"][mean_describe]-(2*small.describe()["meanCoverage"][std_describe])
            passedMSF=small[small["meanCoverage"]>=Mean_2SD].shape[0]
            passedMS=(passedMSF*100)/small["meanCoverage"].shape[0]

            #Many figures and their titles describing the bins are showing that normal destribution is not followed always but Poisson seems more accurate distribution
            # in order to describe our data.

            # a new column into Score.txt file would be added in regards to which distribution is followed in the future.

            if passedMS>=100:
               # print("Normal Distribution??")
                Mean_2SD=small.describe()["meanCoverage"][min_describe]
                passedMSF=small[small["meanCoverage"]>=Mean_2SD].shape[0]
                passedMS=(passedMSF*100)/small["meanCoverage"].shape[0]

	#######################WORKING.................
	#IF BAD SCORE SHOW THE MAX and everything else is above that max in the figure title.
	#########################NOT DONE

        elif small["meanCoverage"].shape[0]>=minWindowsNeeded and Mean_2SD<0:
            Q1=small.describe()["meanCoverage"][Q1_describe]
            passedQ1F=small[small["meanCoverage"]>=Q1].shape[0]
            passedQ1=(passedQ1F*100)/small["meanCoverage"].shape[0]
            ############################ mean not usefull lets use Q3 this time
            Mean_2SD=small.describe()["meanCoverage"][Q3_describe]     #this is Q3 with the name of Mean_2SD
            passedMSF=small[small["meanCoverage"]>=Mean_2SD].shape[0]
            passedMS=(passedMSF*100)/small["meanCoverage"].shape[0]


	#This is part of an old version and should not be called at all with then new changes made( the_least windowss that we need is >=5)


        else: # We will use min and max cause not enough data to get a proper Q1 or mean or Q3.
            Q1=small.describe()["meanCoverage"][min_describe] # this would be our min with a different name
            passedQ1F=small[small["meanCoverage"]>=Q1].shape[0]
            passedQ1=(passedQ1F*100)/small["meanCoverage"].shape[0]

            Mean_2SD=small.describe()["meanCoverage"][max_describe] # this would be out max with different name
            passedMSF=small[small["meanCoverage"]<=Mean_2SD].shape[0]
            passedMS=(passedMSF*100)/small["meanCoverage"].shape[0]

        passedMS=round(passedMS,2)
        passedQ1=round(passedQ1,2)
        Mean_2SD=round(Mean_2SD,5)
        Q1=round(Q1,5)
        Scoreforlowmeantxt=Score
        Scoreforlowmean=round(Score,4)
        Scoreforhighmean=round(Score,2)
        if Scoreforlowmean==0.0:
            Scoreforlowmean=round(Score,10)

###################################Figures _____MC & RC##########################

        with plt.rc_context({'axes.edgecolor':'orange', 'xtick.color':'black', 'ytick.color':'green', 'figure.facecolor':'white'}):
            # Temporary rc parameters in effect
            fig, ax1 = plt.subplots()
           # plt.figure(figsize=(48,36), dpi= 70)
            my_dpi=96
            multipdpi=3
            plt.figure(figsize=(8000/my_dpi, 8000/my_dpi), dpi=my_dpi)
#Used to change name in Figure-Title depending on results. And title will include info about Windows and the Mean Coverage each one has.
            #Score=round(Score,2)
#            fig.suptitle('Score= {}'.format(Scoreforlowmean), fontsize=12,fontweight='bold')
            if   Score>=Scorelimit and small["meanCoverage"].shape[0]>=minWindowsNeeded:          #Scorelimit is adjustable some lines above(at the top)

                fig.suptitle('Score= {}'.format(Scoreforlowmean), fontsize=12,fontweight='bold')
                ax1.set_title('\n{} \n\n {}% bins had mean Coverage above {} \n {}% bins had mean Coverage above {}'.format(i,passedQ1,Q1,passedMS,Mean_2SD), fontdict={'fontsize': 8, 'fontweight': 'medium'})

            elif Score>=Scorelimit and small["meanCoverage"].shape[0]<=minWindowsNeeded:          #Scorelimit is adjustable some lines above(at the top)

                fig.suptitle('Score= {} ,however data are not enough'.format(Scoreforlowmean), fontsize=12,fontweight='bold')
                ax1.set_title('\n{} \n\n {}% bins had mean Coverage above {} \n {}% bins had mean Coverage below {}'.format(i,passedQ1,Q1,passedMS,Mean_2SD), fontdict={'fontsize': 8, 'fontweight': 'medium'})

            elif Score<Scorelimit and small["meanCoverage"].shape[0]<minWindowsNeeded:          #Scorelimit is adjustable some lines above(at the top)
                if Q1==Mean_2SD:
                    fig.suptitle('Bad Score= {}'.format(Scoreforlowmean), fontsize=12, fontweight='bold')
                    ax1.set_title('\n{}\n\n{}% bins had mean Coverage equal to {}'.format(i,passedMS,Mean_2SD), fontdict={'fontsize': 8, 'fontweight': 'medium'})

                else:
                    fig.suptitle('Bad Score= {}'.format(Scoreforlowmean), fontsize=12, fontweight='bold')
                    ax1.set_title('\n{}\n\n{}% bins had mean Coverage below or equal to {} \n {}% bins had mean Coverage above or equal to{}'.format(i,passedMS,Mean_2SD,passedQ1,Q1), fontdict={'fontsize': 8, 'fontweight': 'medium'})

            else:
                if Q1==Mean_2SD:
                    fig.suptitle('Bad Score= {}'.format(Scoreforlowmean), fontsize=12, fontweight='bold')
                    ax1.set_title('\n{}\n\n{}% bins had mean Coverage above or equal to {}'.format(i,passedMS,Mean_2SD), fontdict={'fontsize': 8, 'fontweight': 'medium'})

                else:
                    fig.suptitle('Bad Score= {}'.format(Scoreforlowmean), fontsize=12, fontweight='bold')
                    ax1.set_title('\n{}\n\n{}% bins had mean Coverage above or equal to {} \n {}% bins had mean Coverage above or equal to{}'.format(i,passedMS,Mean_2SD,passedQ1,Q1), fontdict={'fontsize': 8, 'fontweight': 'medium'})

            ax1.set(ylabel='Read Counts',ylim=(0,30))
            ax1.set(xlabel='Position',xlim=(0,B1maxStart))
            ax1.yaxis.label.set_color('blue')
            #ax1.stem(small["chromStart"],small["readCount"], markerfmt='bo',linefmt="tab:green",use_line_collection=True)
            markers,stems, base =ax1.stem(small["chromStart"],small["readCount"], markerfmt='bo',linefmt="tab:cyan")
            for stem in stems:
                use_line_collection=True
                stem.set_linewidth(10)
                #stem.set_use_line_collection=True
            ax1.tick_params(axis='y',labelcolor="blue")
            ax1.grid()
            ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis
            #print("MAX mean coverage:",np.ceil(max(B1["meanCoverage"])))
            #print("MAX read counts:", np.ceil(max (B1["readCount"])))
            #assign("last.warning"=NULL,envir=baseenv())
            ax2.set( ylabel='Mean Coverage',ylim=(0,B1maxcove))
            ax2.yaxis.label.set_color('red')
            ax2.stem(small["chromStart"],small["meanCoverage"], markerfmt='rd',linefmt="white",use_line_collection=True)  # bottom 1 could be used for those with a good Score and bad Mean Coverage
            ax2.tick_params(axis='y',labelcolor="tab:red")
            fig.tight_layout()
            RS=round(Score)
            #organism1=i.split("_")[3]+"_"+i.split("_")[4]+"_"+i.split("_")[5]+"_"+i.split("_")[6]
            #organism1=i.split("_")[3]+"_"+i.split("_")[4]+"_"+i.split("_")[5]+"_"+i.split("_")[6]
            organism0=i.split("_")[1]
            organism1=i.split("_")[3]
            organism2=i.split("_")[4]
            organism3=i.split("_")[5]
            organism4=i.split("_")[6]
            if organism1=="":
                organism1="_"
               # print("O")
            if organism2=="":
                organism1="_"
               # print("O")
            if organism3=="":
                organism1="_"
               # print("O")
            if organism4=="":
                organism1="_"
               # print("O")

#NEW COLUMN FOR MORE SCORES
#YES NO/ YES YES/ NO NO/ NO YES

#RENAMING the figures depending on their score
#                           &
#Making the csv of the Scores  \t (tab) seperated

            RMMC=round(B1maxcove2)

            if   Scoreforlowmean>=Scorelimit and B1mincove>=Scorelimitforminormax and len(small['meanCoverage']>=minWindowsNeeded):
                #fig.savefig('{}/{}_MC_RCSCORED{}.png'.format(newdirpathSC,i,RS),dpi=my_dpi*multipdpi)
                fig.savefig('{}/{}_{}.png'.format(newdirpathSC,i,RS),dpi=my_dpi*multipdpi)
                GENBANKID.write('\n%s'%(organism0))

#                if Scoreforlowmean>AnotherScorelimit:
                    #f_score.write('%s\t%s\t---------->>>>>>YES<<<<<<----------\n' % (organism, Scoreforlowmeantxt))
                f_score.write('\n%s\t%s\t%s\t%s\t%s\t%s\t---------->>>>>>YES<<<<<<----------\t---------->>>>>>YES<<<<<<----------' % (organism0,organism1,organism2,organism3,organism4, Scoreforlowmeantxt))
#                else:
#                    f_score.write('\n%s\t%s\t%s\t%s\t%s\t%s\t---------->>>>>>YES<<<<<<----------\tNO' % (organism0,organism1,organism2,organism3,organism4, Scoreforlowmeantxt))

            #needs names.+ FOLDER NAMES******************************
            elif Scoreforlowmean>=Scorelimit and B1mincove<Scorelimitforminormax and len(small['meanCoverage']>=minWindowsNeeded):
                #fig.savefig('{}/{}_MC_RC_BDSC_Might{}.png'.format(mightdirpath,i,RMMC),dpi=my_dpi*multipdpi)
                fig.savefig('{}/{}_{}.png'.format(mightdirpath,i,RMMC),dpi=my_dpi*multipdpi)
#                if Scoreforlowmean>=AnotherScorelimit:
                    #f_score.write('%s\t%s\t---------->>>>>>YES<<<<<<----------\n' % (organism, Scoreforlowmeantxt))
                f_score.write('\n%s\t%s\t%s\t%s\t%s\t%s\t---------->>>>>>YES<<<<<<----------\t---------->>>>>>YES<<<<<<----------' % (organism0,organism1,organism2,organism3,organism4, Scoreforlowmeantxt))
#                else:
#                    f_score.write('\n%s\t%s\t%s\t%s\t%s\t%s\tNO\tNO' % (organism0,organism1,organism2,organism3,organism4, Scoreforlowmeantxt))

            elif Scoreforlowmean<Scorelimit and B1maxcove2>=Scorelimitforminormax and len(small['meanCoverage']>=minWindowsNeeded):
                #fig.savefig('{}/{}_MC_RC_BDSC_Might{}.png'.format(mightdirpath,i,RMMC),dpi=my_dpi*multipdpi)
                fig.savefig('{}/{}_{}.png'.format(mightdirpath,i,RMMC),dpi=my_dpi*multipdpi)
#                if Scoreforlowmean>=AnotherScorelimit:
                    #f_score.write('%s\t%s\t---------->>>>>>YES<<<<<<----------\n' % (organism, Scoreforlowmeantxt))
                f_score.write('\n%s\t%s\t%s\t%s\t%s\t%s\tNO\t---------->>>>>>YES<<<<<<----------' % (organism0,organism1,organism2,organism3,organism4, Scoreforlowmeantxt))
#                else:
#                    f_score.write('\n%s\t%s\t%s\t%s\t%s\t%s\tNO\tNO' % (organism0,organism1,organism2,organism3,organism4, Scoreforlowmeantxt))

            else:
		#No reason anymore to create so many figures of none usefull data-organisms
                if BDSC==1:
                    fig.savefig('{}/{}_MC_RC_BDSC.png'.format(newdirpath,i),dpi=my_dpi*multipdpi)
                else:
                    pass
		################################################ This can be adjusted. Would have been an alternative second score.

#                if Scoreforlowmean>=AnotherScorelimit:
                    #f_score.write('%s\t%s\t---------->>>>>>YES<<<<<<----------\n' % (organism, Scoreforlowmeantxt))
                f_score.write('\n%s\t%s\t%s\t%s\t%s\t%s\tNO\tNO' % (organism0,organism1,organism2,organism3,organism4, Scoreforlowmeantxt))
#                else:
#                    f_score.write('\n%s\t%s\t%s\t%s\t%s\t%s\tNO\tNO' % (organism0,organism1,organism2,organism3,organism4, Scoreforlowmeantxt))



        #    print(small)
            plt.close("all")


    print("########### File completed #######")
f_score.close()
print("Wochenende plot script completed  \t*\t*\t*\t*\t* \n\n")
#print(i[11:])


