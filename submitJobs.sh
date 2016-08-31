#!/bin/bash

WORKDIR=/afs/cern.ch/user/m/mandrews/work/PHOTONID

#bsub -J Ntuplizer -o $WORKDIR/logs/stdout.log -e $WORKDIR/logs/stderr.log -q 8nh < $WORKDIR/runNtuplizer.sh
bsub -q 1nd < $WORKDIR/runNtuplizer.sh
