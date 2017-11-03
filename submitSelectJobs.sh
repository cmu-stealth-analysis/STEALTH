#!/bin/bash

for era in "B" "C" "D" "E" "F"; do
    echo "Submitting job to select events for era ${era}"
    bsub -q 1nd submitSelectJobsHelper.sh ${era}
done
