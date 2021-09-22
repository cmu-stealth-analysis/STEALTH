#include "../include/makeHistogramsLowST.h"

void fill_distributions_from_file(const std::vector<std::string> & inputPaths, TH2D & dist2D, std::map<int, TH1D> & distST, TH1D & distNJets, const double & minST_nJetsDistributions) {
  std::cout << "Filling distributions from inputs..." << std::endl;
  TChain inputChain("ggNtuplizer/EventTree");
  inputChain.SetMaxTreeSize(100000000000); // 1 TB

  for (const std::string & inputPath: inputPaths) {
    std::cout << "Adding events from file: " << inputPath << std::endl;
    int read_status = inputChain.Add(inputPath.c_str(), 0);
    assert(read_status == 1);
  }

  inputChain.SetBranchStatus("*", 0); // so that only the needed branches, explicitly activated below, are read in per event
  // TTreeReader inputTreeReader(inputChain);
  // TTreeReaderValue<float> evt_ST(inputTreeReader, "b_evtST");
  // TTreeReaderValue<int> evt_nJets(inputTreeReader, "b_nJetsDR");
  // TTreeReaderValue<double> evt_MCXSecWeight(inputTreeReader, "b_MCXSecWeight");
  // TTreeReaderValue<float> evt_prefiringWeight(inputTreeReader, "b_evtPrefiringWeight");
  // TTreeReaderValue<float> evt_photonMCScaleFactor(inputTreeReader, "b_evtphotonMCScaleFactor");
  float evt_ST = -1.;
  inputChain.SetBranchStatus("b_evtST", 1);
  inputChain.SetBranchAddress("b_evtST", &evt_ST);
  int evt_nJetsDR = -1.;
  inputChain.SetBranchStatus("b_nJetsDR", 1);
  inputChain.SetBranchAddress("b_nJetsDR", &evt_nJetsDR);
  double evt_MCXSecWeight = -1.;
  inputChain.SetBranchStatus("b_MCXSecWeight", 1);
  inputChain.SetBranchAddress("b_MCXSecWeight", &evt_MCXSecWeight);
  float evt_prefiringWeight = -1.;
  inputChain.SetBranchStatus("b_evtPrefiringWeight", 1);
  inputChain.SetBranchAddress("b_evtPrefiringWeight", &evt_prefiringWeight);
  float evt_photonMCScaleFactor = -1.;
  inputChain.SetBranchStatus("b_evtphotonMCScaleFactor", 1);
  inputChain.SetBranchAddress("b_evtphotonMCScaleFactor", &evt_photonMCScaleFactor);

  Long64_t nEntries = inputChain.GetEntries();
  std::cout << "Total number of available events: " << nEntries << std::endl;

  // long entryIndex = 0;
  tmProgressBar progressBar = tmProgressBar(static_cast<long>(nEntries));
  long progressBarUpdatePeriod = ((nEntries < 50) ? 1 : static_cast<long>(0.5 + 1.0*(nEntries/50)));
  progressBar.initialize();
  for (Long64_t entryIndex = 0; entryIndex < nEntries; ++entryIndex) {
    // while (inputTreeReader.Next()) {
    // ++entryIndex;
    if ((entryIndex == 0) ||
	(entryIndex % progressBarUpdatePeriod == 0) ||
	(entryIndex == (nEntries-1))) progressBar.updateBar(1.0*entryIndex/nEntries, entryIndex);

    Long64_t loadStatus = inputChain.LoadTree(entryIndex); assert(loadStatus >= 0);
    // if (loadStatus < 0) {
    //   std::cout << "Warning: loadStatus < 0 for entry index: " << entryIndex << "; load status = " << loadStatus << "; fileName: " << (inputChain.GetFile())->GetName() << std::endl;
    //   ++nProblematicEntries;
    //   assert(nProblematicEntries <= N_PROBLEMATIC_ENTRIES_THRESHOLD);
    //   continue;
    // }
    int nBytesRead = inputChain.GetEntry(entryIndex, 0); // Get only the required branches
    assert(nBytesRead > 0);
    // if (nBytesRead <= 0) {
    //   std::cout << "Warning: failed to read SOME information from entry at index: " << entryIndex << "; nBytesRead = " << nBytesRead << "; fileName: " << (inputChain.GetFile())->GetName() << std::endl;
    //   ++nProblematicEntries;
    //   assert(nProblematicEntries <= N_PROBLEMATIC_ENTRIES_THRESHOLD);
    //   continue;
    // }

    int nJetsBin = evt_nJetsDR;
    if (nJetsBin > 6) nJetsBin = 6;

    if (evt_ST > 1300.) continue;

    double STBinWidth = dist2D.GetXaxis()->GetBinWidth(evt_ST);
    double eventWeightNoBinWidth = ((evt_MCXSecWeight)*(evt_prefiringWeight)*(evt_photonMCScaleFactor));

    dist2D.Fill(evt_ST, nJetsBin, eventWeightNoBinWidth/STBinWidth);
    (distST.at(nJetsBin)).Fill(evt_ST, eventWeightNoBinWidth/STBinWidth);
    bool fill_distNJets = true;
    if (minST_nJetsDistributions > 0.) fill_distNJets = (evt_ST >= minST_nJetsDistributions);
    if (fill_distNJets) distNJets.Fill(nJetsBin, eventWeightNoBinWidth);
  }
  progressBar.terminate();
  std::cout << "Done." << std::endl;
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  TH1::AddDirectory(kFALSE);
  tmArgumentParser argumentParser = tmArgumentParser("Save ST distributions at low ST.");
  argumentParser.addArgument("inputFilePaths", "", true, "Path to files with n-tuplized events (if more than one, separate by a semicolon).");
  argumentParser.addArgument("outputFolder", "", true, "Output folder.");
  argumentParser.addArgument("outputFileName", "", true, "Output file name.");
  argumentParser.addArgument("nJetsDistributionsMinST", "-1.0", false, "Min value of ST for events to contribute to the 1D nJets distributions.");
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);

  std::cout << "Options passed:" << std::endl << options << std::endl;

  double STEdges[14] = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1100., 1200., 1300.};

  TH2D distribution2D("dist2D", "2D distribution (weighted);ST;nJets;weighted events/GeV", 13, STEdges, 7, -0.5, 6.5);
  distribution2D.Sumw2();
  std::map<int, TH1D> STDistributions;
  for (int nJetsBin = 0; nJetsBin <= 6; ++nJetsBin) {
    STDistributions[nJetsBin] = TH1D(("distST_" + std::to_string(nJetsBin) + "JetsBin").c_str(), "ST distribution (weighted);ST;weighted events/GeV", 13, STEdges);
    (STDistributions.at(nJetsBin)).Sumw2();
  }
  TH1D nJetsDistribution = TH1D("dist_nJets", "nJets distribution (weighted);nJets;weighted events", 7, -0.5, 6.5);
  nJetsDistribution.Sumw2();
  fill_distributions_from_file(options.inputFilePaths, distribution2D, STDistributions, nJetsDistribution, options.nJetsDistributionsMinST);

  TFile *outputFile = TFile::Open((options.outputFolder + "/" + options.outputFileName).c_str(), "RECREATE");
  outputFile->WriteTObject(&distribution2D);
  for (int nJetsBin = 0; nJetsBin <= 6; ++nJetsBin) {
    outputFile->WriteTObject(&(STDistributions.at(nJetsBin)));
  }
  outputFile->WriteTObject(&nJetsDistribution);
  outputFile->Close();

  return EXIT_SUCCESS;
}
