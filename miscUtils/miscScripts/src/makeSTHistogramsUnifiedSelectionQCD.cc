#include "../include/makeSTHistogramsUnifiedSelectionQCD.h"

float getEventWeight(float genHT, const std::vector<float>& eventWeights) {
  // hack for now
  if (genHT < 300.) return 0.;
  if (genHT < 500.) return eventWeights[0];
  if (genHT < 700.) return eventWeights[1];
  if (genHT < 1000.) return eventWeights[2];
  if (genHT < 1500.) return eventWeights[3];
  if (genHT < 2000.) return eventWeights[4];
  return eventWeights[5];
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  tmArgumentParser argumentParser = tmArgumentParser("Merge outputs of event selection script into a single file.");
  argumentParser.addArgument("inputFilesList", "", true, "Path to file containing list of paths with n-tuplized events.");
  argumentParser.addArgument("outputFolder", "", true, "Output folder.");
  argumentParser.addArgument("outputFileName", "", true, "Name of output file.");
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);

  std::cout << "Options passed:" << std::endl << options << std::endl;

  std::vector<float> eventWeights;
  eventWeights.push_back(41.5273/0.003089);
  eventWeights.push_back(41.5273/0.03316);
  eventWeights.push_back(41.5273/0.1566);
  eventWeights.push_back(41.5273/0.9067);
  eventWeights.push_back(41.5273/9.867);
  eventWeights.push_back(41.5273/47.5);

  std::vector<std::string> inputFileNames;
  std::ifstream fileWithInputFilesList((options.inputFilesList).c_str());
  if (!fileWithInputFilesList.is_open()) {
    std::cout << "ERROR: Failed to open file with path: " << options.inputFilesList << std::endl;
    std::exit(EXIT_FAILURE);
  }
  while (!fileWithInputFilesList.eof()) {
    std::string inputFileName;
    fileWithInputFilesList >> inputFileName;
    if (!(inputFileName.empty())) inputFileNames.push_back(inputFileName);
  }
  fileWithInputFilesList.close();

  TChain* inputChain = new TChain("ggNtuplizer/EventTree");
  inputChain->SetMaxTreeSize(100000000000); // 1 TB

  for (auto&& inputFileName: inputFileNames) {
    std::cout << "Adding events from file: " << inputFileName << std::endl;
    int read_status = inputChain->Add(inputFileName.c_str(), 0);
    assert(read_status == 1);
  }

  Long64_t nEntries = inputChain->GetEntries();
  std::cout << "Total number of available events: " << nEntries << std::endl;

  TFile *outputFile = TFile::Open((options.outputFolder + "/" + options.outputFileName).c_str(), "RECREATE");
  std::map<int, TH1F*> STHistograms;
  for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
    std::string histogramName = "STDistribution_" + std::to_string(nJetsBin) + "JetsBin";
    std::string histogramTitle = "STDistribution, " + std::to_string(nJetsBin) + " Jets Bin";
    STHistograms[nJetsBin] = new TH1F(histogramName.c_str(), histogramTitle.c_str(), 51, 0., 5100.);
  }

  TTreeReader inputTreeReader(inputChain);
  TTreeReaderValue<float> evt_ST(inputTreeReader, "b_evtST");
  TTreeReaderValue<int> evt_nJets(inputTreeReader, "b_nJetsDR");
  TTreeReaderValue<float> evt_genHT(inputTreeReader, "genHT");

  long entryIndex = 0;
  tmProgressBar progressBar = tmProgressBar(static_cast<long>(nEntries));
  long progressBarUpdatePeriod = ((nEntries < 50) ? 1 : static_cast<long>(0.5 + 1.0*(nEntries/50)));
  progressBar.initialize();
  while (inputTreeReader.Next()) {
    if (entryIndex % progressBarUpdatePeriod == 0) progressBar.updateBar(1.0*entryIndex/nEntries, entryIndex);
    ++entryIndex;

    int nJetsBin = *evt_nJets;
    if (nJetsBin < 2) continue;
    if (nJetsBin > 6) nJetsBin = 6;

    STHistograms[nJetsBin]->Fill(*evt_ST, getEventWeight(*evt_genHT, eventWeights));
  }

  for (int nJetsBin = 2; nJetsBin <= 6; ++nJetsBin) {
    outputFile->WriteTObject(STHistograms[nJetsBin]);
  }

  outputFile->Close();

  return EXIT_SUCCESS;
}
