# Stealth SUSY Analysis Instructions
## Carnegie Mellon Group
### (M. Andrews, N. Bower, T. Mudholkar, M. Paulini, M. Sun, M. Weinberg)

#### For full details about the compulsory and optional arguments to pass to any executable (including the C++ executables!), please do: `executableName -h`; for example, `./mergeFiles.py -h`

## Creating N-tuples
### Running ntuplizer
Please see: [ggNtuplizer with modifications](https://github.com/tanmaymudholkar/ggAnalysis/tree/94X)

The crab utility is documented [here](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCrab).

### Event selection
To run selection: `./submitEventSelectionJobs.py --optionalIdentifier example_identifier`

By default this won't submit the jobs, only create the jdl files which you can examine. To submit the jobs, pass the flag `--isProductionRun`.

Outputs will be stored in the folder `/store/user/lpcsusystealth/selections/DoublePhoton<_optional_identifier_if_set>`

### Merging output files
To merge all the files created by the selection script: `./runSelectionMerge.py --optionalIdentifier example_identifier`

## Full analysis chain
To run full analysis chain: `./runAnalysis.py --selectionSuffix <example_selection_identifier> --optionalIdentifier <example_analysis_identifier>`

To run a specific chain, do: `./runAnalysis.sh --chain "type"`.

If the flag `--runUnblinded` is passed, the chain runs with the signal region unblinded: the observed data for the signal sample is now plotted as well, and observed limit contours are drawn in addition to the expected limit contours.

* Chain `data`: This chain examines ST distributions in the control and signal region **data**. This outputs the uncertainty estimates including the ST scaling estimate from the control data, the expected NEvents distributions, and also the number of events recorded (outside the blinded region).
* Chain `MC`: This chain examines ST distributions in the signal region **MC samples**. This chain runs over the full signal MC n-tuples and generates histograms of the expected number of events in each (m_gluino, m_neutralino) bin, as well as the expected number of events with the various shifts. In addition, this chain generates MC systematics (luminosity, statistical uncertainty on number of passing events, JEC, unclustered MET, JER, prefiring, photon scale factor) using the distributions created in the `data` chain.
* Chain `combine`: Using the expected and observed background number of events from the `data` chain, the expected number of signal events from the `MC` chain, and all systematic uncertainties from both chains, this chain creates three data cards per (m_gluino, m_neutralino) bin: with the nominal cross section and with it shifted up and down by the theoretical uncertainty. In addition, this chain submits the Condor jobs to run the combine tool on these datacards.
* Chain `signalContamination`: This chain generates histograms of the potential signal contamination throughout the (m_gluino, m_neutralino) phase-space.
* Chain `ancillaryPlots`: This chain generates publication-quality histograms of the expected number of events (for the signal sample) and the observed number of events (for the control sample) with a few useful distributions from MC overlaid.
* Chain `limits`: This chain generates the 95% expected limit plots with contours at signal_strength = 1.
