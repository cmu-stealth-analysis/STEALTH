# STEALTH
README for 2016 Stealth Photon Analysis

(1) Given a set of ggNtuples, one can first run a skim to keep only evts with >= 1 loose photon:

$ python skimPhoLoose_Data.py

with the input and output files defined therein.

(2) Then, using the skims, run the selection criteria on them:

$ python selectEvents.py -s A -e G -H 60

where:

  -s selection scheme A(2 photons) or B(1 photon)
  
  -e run era: C,D,...
  
  -H evt HT cut

This will additionally store the ST and nJets for the event in the branches b_evtST and b_nJets, respectively.

(3) To create ST plots of the selection output, do:

$ python plotST.py -i file1.root file2.root -l 1000. -r 3500. -b 5

where:

  -i is the space-delimited list of input files
  
  -l lower ST limit
  
  -r upper ST limit
  
  -b number of bins to plot over

This additionally outputs the ff. which are used for bkg estimation: 

hSTs.root -> raw histograms of the nJet distributions

normRatios.txt -> scaling ratios between control and signal jet distributions

For running batch jobs of the ST plots, you may prefer to use the runPlot.py script:

$ ./runPlot.py

##############################
!! Still undergoing refinement: !!

(4) To run the bkg estimation, do: 

$ python getBkg.py

which will use the hSTs.root and normRatios.txt from the last execution of plotST.py
