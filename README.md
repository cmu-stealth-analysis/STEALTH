# Stealth SUSY Analysis Instructions
## Carnegie Mellon Group
### (M. Andrews, N. Bower, T. Mudholkar, M. Paulini, M. Sun, M. Weinberg)

#### For full details about the compulsory and optional arguments to pass to any executable (including the C++ executables!), please do: `executableName -h`; for example, `./mergeFiles.py -h`

## Creating N-tuples
### Running ntuplizer
Please see: [ggNtuplizer with modifications](https://github.com/tanmaymudholkar/ggAnalysis/tree/94X)

The crab utility is documented [here](https://twiki.cern.ch/twiki/bin/view/CMSPublic/SWGuideCrab).

### Event selection
To run complete selection for 2017: `./submitEventSelectionJobs.py --year 2017`

By default this won't submit the jobs, only create the jdl files which you can examine. To submit the jobs, pass the flag `--isProductionRun`.

If the jobs have been submitted at least once before, the total number of events in the input files is likely stored in a cache file, which saves a lot of time. By default the script does not read the number of events from the cache. To enable the cache, pass flag `--enable_cache`.

### Merging output files
To merge all the files created by the selection script: `./runSelectionMerge.sh`

## Full analysis chain
To run full analysis chain: `./runAnalysis.sh`

To run a specific step N in the chain, do: `./runAnalysis.sh specificStep=N`.

* Step 1: Examine ST distributions in the control region **data**. This outputs the uncertainty estimate due to a deviation from perfect ST scaling.
* Step 2: Examine ST distributions in the signal region **data**. This outputs the uncertainty estimates due to an error on ST shape, due to the statistical uncertainty on the number of events in the normalization region, and due to the uncertainty on the parameter $\rho$. In addition this outputs the expected and observed number of events in each signal region.
* Step 3: Examine ST distributions in the signal region **MC samples**. This step runs over the full signal MC n-tuples and generates histograms of the expected number of events in each (m_gluino, m_neutralino) bin, as well as the expected number of events with the various shifts.
* Step 4: Generate MC systematics (luminosity, statistical uncertainty on number of passing events, JEC, unclustered MET, JER, prefiring, photon scale factor) using the distributions created in step 3.
* Step 5: Using the expected and observed number of events created in steps 2 and 4, and using the fractional uncertainties created in steps 1, 2, and 4, create three data cards per (m_gluino, m_neutralino) bin: with the nominal cross section and with it shifted up and down by the theoretical uncertainty.
* Step 6: Submit Condor jobs to run the combine tool on the datacards created in step 5.

## Make plots for paper

These scripts create plots to be used in the paper:
* `./plotLimits.py`: make 2D distribution of limit plots with proper TDR-style formatting
* `./plotSTDistributionsWithErrors.py`: plot expected and observed ST distributions with correct errors and proper TDR-style formatting
