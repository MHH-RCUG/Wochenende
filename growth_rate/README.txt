Requirements:
python>=3.4
matplotlib
numpy
pandas
tqdm
lmfit	(https://lmfit.github.io/lmfit-py/installation.html)


# Automated - preferred
Run wochenende_postprocess.sh after the Wochenende pipeline. Data preprocessing and growth rate calculations will be made as part of the pipeline.

# Manual
Copy all the bam files from the wanted samples and the associated bam.txt files in the growth_rate_code directory.
Activate the wochenende enviroment.
First run the runbatch_bed_to_csv.sh script.
Then run the run_reproduction_determiner.sh


In general, the program takes as input the path of a folder containing one or more csv files. Each csv file contains the read positions of an organism.
It produces a single result csv file named 'results_summary.csv'. Additionally, a png file is saved for each input csv file visualizing the fitting result, if 'save_plots' is set to true.
The results are saved in fit_results/experiment_name/data_folder_name/
'experiment_name' can be set in main.py

The result csv file consists of 7 columns and a row for each organism.
In case the growth rate cannot be computed, it is set to 0.
When the fitting error cannot be computed, it is set to -1.

Error codes explanation:
[] : no error
-1 : insufficient median coverage
-2 : not enough reads left after filtering
-3 : not enough bins left after filtering
-4 : growth rate way too high, probably corrupt fitting
-5 : fit error too large


