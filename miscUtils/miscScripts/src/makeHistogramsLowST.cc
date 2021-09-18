#include "../include/makeHistogramsLowST.h"

void fill_distributions_from_file(const std::string & inputPath, TH2D & dist2D, std::map<int, TH1D> & distST, TH1D & distNJets) {
  std::cout << "Filling distributions from file at path: " << inputPath << std::endl;
  TChain* inputChain = new TChain("ggNtuplizer/EventTree");
  inputChain->SetMaxTreeSize(100000000000); // 1 TB

  // for (auto&& inputFileName: inputFileNames) {
  //   std::cout << "Adding events from file: " << inputFileName << std::endl;
  //   int read_status = inputChain->Add(inputFileName.c_str(), 0);
  //   assert(read_status == 1);
  // }

  int read_status = inputChain->Add(inputPath.c_str(), 0);
  assert(read_status == 1);

  Long64_t nEntries = inputChain->GetEntries();
  std::cout << "Total number of available events: " << nEntries << std::endl;

  TTreeReader inputTreeReader(inputChain);
  TTreeReaderValue<float> evt_ST(inputTreeReader, "b_evtST");
  TTreeReaderValue<int> evt_nJets(inputTreeReader, "b_nJetsDR");
  TTreeReaderValue<double> evt_MCCustomWeight(inputTreeReader, "b_MCCustomWeight");
  TTreeReaderValue<float> evt_prefiringWeight(inputTreeReader, "b_evtPrefiringWeight");
  TTreeReaderValue<float> evt_photonMCScaleFactor(inputTreeReader, "b_evtphotonMCScaleFactor");
  long entryIndex = 0;
  tmProgressBar progressBar = tmProgressBar(static_cast<long>(nEntries));
  long progressBarUpdatePeriod = ((nEntries < 50) ? 1 : static_cast<long>(0.5 + 1.0*(nEntries/50)));
  progressBar.initialize();
  while (inputTreeReader.Next()) {
    ++entryIndex;
    if ((entryIndex == 1) ||
	(entryIndex % progressBarUpdatePeriod == 0) ||
	(entryIndex == nEntries)) progressBar.updateBar(1.0*entryIndex/nEntries, entryIndex);

    int nJetsBin = *evt_nJets;
    if (nJetsBin > 6) nJetsBin = 6;

    double STBinWidth = dist2D.GetXaxis()->GetBinWidth(*evt_ST);
    double eventWeightNoBinWidth = (((*evt_MCCustomWeight)/1000.0)*(*evt_prefiringWeight)*(*evt_photonMCScaleFactor));

    dist2D.Fill(*evt_ST, nJetsBin, eventWeightNoBinWidth/STBinWidth);
    (distST.at(nJetsBin)).Fill(*evt_ST, eventWeightNoBinWidth/STBinWidth);
    if (*evt_ST >= 800.) distNJets.Fill(nJetsBin, eventWeightNoBinWidth);
  }
  progressBar.terminate();
  std::cout << "Done." << std::endl;
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  TH1::AddDirectory(kFALSE);
  tmArgumentParser argumentParser = tmArgumentParser("Save ST distributions at low ST.");
  argumentParser.addArgument("inputFilePath", "", true, "Path to file with n-tuplized events.");
  argumentParser.addArgument("outputFolder", "", true, "Output folder.");
  argumentParser.addArgument("outputFileName", "", true, "Output file name.");
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);

  std::cout << "Options passed:" << std::endl << options << std::endl;

  // std::vector<std::string> inputFileNames;
  // std::ifstream fileWithInputFilesList((options.inputFilesList).c_str());
  // if (!fileWithInputFilesList.is_open()) {
  //   std::cout << "ERROR: Failed to open file with path: " << options.inputFilesList << std::endl;
  //   std::exit(EXIT_FAILURE);
  // }
  // while (!fileWithInputFilesList.eof()) {
  //   std::string inputFileName;
  //   fileWithInputFilesList >> inputFileName;
  //   if (!(inputFileName.empty())) inputFileNames.push_back(inputFileName);
  // }
  // fileWithInputFilesList.close();

  double STEdges[13] = {0., 100., 200., 300., 400., 500., 600., 700., 800., 900., 1000., 1200., 1500.};

  TH2D distribution2D("dist2D", "2D distribution (weighted);ST;nJets;weighted events/GeV", 12, STEdges, 7, -0.5, 6.5);
  distribution2D.Sumw2();
  std::map<int, TH1D> STDistributions;
  for (int nJetsBin = 0; nJetsBin <= 6; ++nJetsBin) {
    STDistributions[nJetsBin] = TH1D(("distST_" + std::to_string(nJetsBin) + "JetsBin").c_str(), "ST distribution (weighted);ST;weighted events/GeV", 12, STEdges);
    (STDistributions.at(nJetsBin)).Sumw2();
  }
  TH1D nJetsDistribution = TH1D("dist_nJets", "nJets distribution (weighted);nJets;weighted events", 7, -0.5, 6.5);
  nJetsDistribution.Sumw2();
  fill_distributions_from_file(options.inputFilePath, distribution2D, STDistributions, nJetsDistribution);

  TFile *outputFile = TFile::Open((options.outputFolder + "/" + options.outputFileName).c_str(), "RECREATE");
  outputFile->WriteTObject(&distribution2D);
  for (int nJetsBin = 0; nJetsBin <= 6; ++nJetsBin) {
    outputFile->WriteTObject(&(STDistributions.at(nJetsBin)));
  }
  outputFile->WriteTObject(&nJetsDistribution);
  outputFile->Close();

  return EXIT_SUCCESS;
}
