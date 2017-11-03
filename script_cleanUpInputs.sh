#!/bin/bash

for era in "B" "C" "D" "E" "F" "G" "H"; do
    echo "Cleaning up input files for era = ${era}"
    cat inputFileList_JetHT_era${era}.txt | grep -v "failed" | sed "s|^/eos/uscms/|root://cmsxrootd.fnal.gov//|"  > tmp.txt
    mv tmp.txt inputFileList_JetHT_era${era}.txt
done
