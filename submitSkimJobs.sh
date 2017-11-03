#!/bin/bash

for era in "B" "C" "D" "E" "F"; do
    echo "Submitting job for era ${era}"
    bsub -q 2nd submitSkimJobsHelper.sh ${era}
done
