#!/bin/bash

out_file='st_scaling_params.dat'
if [ -f $out_file ]; then
  rm st_scaling_params.dat 
fi

#nJtMax=5 #actual number of jets
nJtMax=6 #actual number of jets
for ((ijt=0;ijt<=$((nJtMax-2));ijt++)); do
  python getBkg_Tweak.py -j $ijt
	#mv bkgfit.png BKG/GJet_selA/bkgFit_$((ijt+2))jet.png 
	#mv bkgfit.png BKG/JetHT_SepRereco/bkgFit_$((ijt+2))jet.png 
	mv bkgfit.png BKG/JetHT_ReminiAOD/bkgFit_$((ijt+2))jet.png 
	#mv bkgfit.png BKG/DiPhoton_selA/bkgFit_$((ijt+2))jet.png 
done
