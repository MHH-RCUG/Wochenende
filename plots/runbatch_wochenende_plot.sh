#!/bin/bash
# Colin Davenport, Jan 2020
# Run Wochenende plotting as batch

for i in $(ls *cov_window.txt.filt.csv)
        do
	# filter out human reads and sort by total reads per window
	python wochenende_plot.py "$i" > output.log

done

