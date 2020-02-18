#!/bin/bash
# Colin Davenport, Jan 2020

# Run posthoc filters on a directory of output from the Wochenende pipeline after all jobs finish

set -euo pipefail

# 1 filter and sort bam.txt files
bash runbatch_metagen_awk_filter.sh

# 2 calculate depth of coverage for final bam files in the directory
bash runbatch_sambamba_depth.sh

# 3 filter out human genome depth coverage and sort the bacterial regions by read coverage
bash runbatch_metagen_window_filter.sh


# Check exit codes
if [ $? -eq 0 ]; then
    echo All steps OK
else
    echo Script FAILED
fi

