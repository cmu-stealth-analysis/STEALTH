#!/bin/bash

rm -r publicationPlots/*
./plotSTDistributionsWithErrors.py --outputFileName "formatted"
./plotSTDistributionComparisons.py --inputFilePath "${EOSPREFIX}/store/user/lpcsusystealth/selections/combinedControl/data_DoubleEG_201*.root" --outputFileName "control_STDistributions"
./plotLimits.py --outputSuffix "fullChain_r_0_1" --contour_signalStrength 0.1
./plotLimits.py
