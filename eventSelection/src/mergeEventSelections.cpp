#include "../include/mergeEventSelections.h"

bool file_has_zero_events(const std::string & file_path) {
  TFile *test_file = TFile::Open(file_path.c_str(), "READ");
  assert((test_file->IsOpen()) && (!(test_file->IsZombie())));
  TTree *eventTree = (TTree*)(test_file->Get("ggNtuplizer/EventTree"));
  return (eventTree == nullptr);
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  TH1::AddDirectory(kFALSE);

  tmArgumentParser argumentParser = tmArgumentParser("Merge outputs of event selection script into a single file.");
  argumentParser.addArgument("inputFilesList", "", true, "Path to file containing list of paths with n-tuplized events.");
  argumentParser.addArgument("outputFolder", "root://cmseos.fnal.gov//store/user/lpcsusystealth/selections/combined_DoublePhoton", false, "Output folder.");
  argumentParser.addArgument("outputFileName", "", true, "Name of output file.");
  // argumentParser.addArgument("addWeightBranch", "", false, "If this argument is set to w, then a new branch named \"b_MCCustomWeight\" is created for which every event is set to the value w. This is intended to be used in the workflow to allow merged datasets with different event weights.");
  // argumentParser.addArgument("isMC", "false", true, "Takes value \"true\" if there are additional plots relevant for MC samples only.");
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);
  std::cout << "Beginning to merge event selections..." << std::endl;

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
  TH1D * counts_TH1 = nullptr;

  for (auto&& inputFileName: inputFileNames) {
    std::cout << "Adding events from file: " << inputFileName << std::endl;
    if (file_has_zero_events(inputFileName)) std::cout << "WARNING: File with path \"" << inputFileName << "\" has 0 events!" << std::endl;
    else {
      int read_status = inputChain->Add(inputFileName.c_str(), 0);
      assert(read_status == 1);
      TH1D * counts_TH1_input = nullptr;
      TFile * input_file_handle = TFile::Open(inputFileName.c_str(), "READ");
      assert((input_file_handle->IsOpen() && !(input_file_handle->IsZombie())));
      input_file_handle->GetObject("counters_partial", counts_TH1_input);
      if (counts_TH1 == nullptr) {
        counts_TH1 = static_cast<TH1D*>(counts_TH1_input->Clone("counters"));
        assert(counts_TH1 != nullptr);
      }
      else {
        assert(counts_TH1->Add(counts_TH1_input));
      }
      input_file_handle->Close();
    }
  }
  counts_TH1->GetXaxis()->SetBinLabel(1, "analyzed");
  counts_TH1->GetXaxis()->SetBinLabel(2, "selected");
  for (int criterionIndex = eventSelectionCriterionFirst; criterionIndex != static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria); ++criterionIndex) {
    eventSelectionCriterion criterion = static_cast<eventSelectionCriterion>(criterionIndex);
    int offset = 0; std::string label = "";
    offset = 3; label = "pass " + eventSelectionCriterionNames.at(criterion);
    counts_TH1->GetXaxis()->SetBinLabel(offset+criterionIndex, label.c_str());
    offset = 3+(static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria)); label = "pass up to " + eventSelectionCriterionNames.at(criterion);
    counts_TH1->GetXaxis()->SetBinLabel(offset+criterionIndex, label.c_str());
    offset = 3+(2*(static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria))); label = "pass besides " + eventSelectionCriterionNames.at(criterion);
    counts_TH1->GetXaxis()->SetBinLabel(offset+criterionIndex, label.c_str());
  }

  Long64_t nEntries = inputChain->GetEntries();
  std::cout << "Total number of available events: " << nEntries << std::endl;

  TFile *outputFile = TFile::Open(("~/cmslpc_scratch/merged/" + options.outputFileName).c_str(), "RECREATE");
  TDirectory* outputDirectory = outputFile->mkdir("ggNtuplizer");
  outputDirectory->cd();

  if (nEntries == 0) {
    std::cout << "WARNING: no events available to merge!" << std::endl;
    TTree *dummyTree = new TTree("EventTree", ""); // tree with 0 entries
    // double evtWeight = options.eventWeight;
    // if (evtWeight > 0.) dummyTree->Branch("b_MCCustomWeight", &evtWeight, "b_MCCustomWeight/D"); // create a dummy branch
    // TTree *dummyTree = inputChain->CloneTree(0);
    // outputFile->WriteTObject(dummyTree);
    outputFile->WriteTObject(dummyTree);
    outputFile->Close();
  }
  else {
    TTree *outputTree = inputChain->CloneTree(0);
    // double evtWeight = options.eventWeight;
    // if (evtWeight > 0.) outputTree->Branch("b_MCCustomWeight", &evtWeight, "b_MCCustomWeight/D"); // only create a branch for event weight if needed
    tmProgressBar progressBar = tmProgressBar(static_cast<int>(nEntries));
    int progressBarUpdatePeriod = ((nEntries < 50) ? 1 : static_cast<int>(0.5 + 1.0*(nEntries/50)));
    progressBar.initialize();
    for (Long64_t eventIndex = 0; eventIndex < nEntries; ++eventIndex) {
      Long64_t treeStatus = inputChain->LoadTree(eventIndex);
      if (treeStatus < 0) {
        std::cout << "ERROR in LoadTree!" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      Int_t eventStatus = inputChain->GetEntry(eventIndex, 1);
      if (eventStatus == 0) {
        std::cout << "ERROR in GetEntry!" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      if (eventIndex > 0 && ((static_cast<int>(eventIndex) % progressBarUpdatePeriod == 0) || eventIndex == static_cast<int>(nEntries-1))) progressBar.updateBar(static_cast<double>(1.0*eventIndex/nEntries), eventIndex);
      outputTree->Fill();
    }
    progressBar.terminate();

    std::cout << "Checking number of entries in the output file..." << std::endl;
    Long64_t nEvtsOutput = outputTree->GetEntries();
    if (!(nEvtsOutput == nEntries)) {
      std::cout << "ERROR: Unexpected number of events in output tree!" << std::endl;
      std::exit(EXIT_FAILURE);
    }
    // outputFile->Write();
    outputDirectory->WriteTObject(outputTree);
    outputFile->WriteTObject(counts_TH1);
    outputFile->Close();
  }

  int xrdcp_return_status = system(("set -x && xrdcp --silent --nopbar --force --path --streams 15 ~/cmslpc_scratch/merged/" + options.outputFileName + " " + options.outputFolder + "/" + options.outputFileName + " && rm -f ~/cmslpc_scratch/merged/" + options.outputFileName + " && set +x").c_str());
  if (xrdcp_return_status != 0) {
    std::cout << "ERROR: xrdcp likely failed with status "<< xrdcp_return_status << std::endl;
  }

  std::cout << "Merge finished!" << std::endl;
  return 0;
}
