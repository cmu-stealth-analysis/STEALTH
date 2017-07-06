#!/bin/bash

WORKDIR=/afs/cern.ch/user/m/mandrews/work/STEALTH
#ggFILE=SinglePhoton_2016B_evtSt.root

python $WORKDIR/getST.py

cp * $WORKDIR/ggNTUPLES/

