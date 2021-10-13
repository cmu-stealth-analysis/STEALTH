#include "../include/getSumMCWeights.h"

TChain * get_chain_from_input_paths_files(const std::vector<std::string> & inputPathsFiles) {
  TChain * inputChain = new TChain("ggNtuplizer/EventTree");
  inputChain->SetMaxTreeSize(10000000000000LL); // 10 TB
  for (const std::string & inputPathsFile : inputPathsFiles) {
    std::cout << "Adding paths from file: " << inputPathsFile << std::endl;
    std::ifstream inputPathsFileStream;
    inputPathsFileStream.open(inputPathsFile.c_str());
    assert(inputPathsFileStream.is_open());
    while (!(inputPathsFileStream.eof())) {
      std::string inputPath;
      inputPathsFileStream >> inputPath;
      if (!(inputPath.empty())) {
	if (!(file_has_zero_events(inputPath))) inputChain->Add(inputPath.c_str());
      }
    }
    inputPathsFileStream.close();
  }
  return inputChain;
}

void setup_chain(TChain * inputChain, float & event_weight) {
  inputChain->SetBranchStatus("*", 0);
  inputChain->SetBranchStatus("genWeight", 1);
  inputChain->SetBranchAddress("genWeight", &event_weight);
}

genInfoStruct loop_over_chain_events(TChain * inputChain, float & event_weight) {
  long nEntries = inputChain->GetEntries();
  genInfoStruct info;
  info.totalNEvts = 0;
  info.sumWeights = 0.;
  tmProgressBar progressBar(nEntries);
  int tmp = static_cast<int>(0.5 + 1.0*nEntries/20);
  int progressBarUpdatePeriod = tmp > 1 ? tmp : 1;
  progressBar.initialize();
  for (Long64_t entryIndex = 0; entryIndex < nEntries; ++entryIndex) {
    Long64_t loadStatus = inputChain->LoadTree(entryIndex);
    assert(loadStatus >= 0);
    int nBytesRead = inputChain->GetEntry(entryIndex, 0); // Get only the required branches
    assert(nBytesRead > 0);
    if ((entryIndex == 0) ||
        (entryIndex == (nEntries-1)) ||
        ((entryIndex % static_cast<Long64_t>(progressBarUpdatePeriod)) == 0)) progressBar.updateBar(static_cast<double>(1.0*entryIndex/nEntries), entryIndex);
    (info.sumWeights) += event_weight;
    ++(info.totalNEvts);
  }
  progressBar.terminate();
  return info;
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  TH1::AddDirectory(kFALSE);
  std::cout << "Saving sum of weights into output JSON..." << std::endl;
  argumentsStruct arguments = get_command_line_arguments(argc, argv);
  TChain * inputChain = get_chain_from_input_paths_files(arguments.inputPathsFiles);
  float event_weight = -1.;
  setup_chain(inputChain, event_weight);
  genInfoStruct info = loop_over_chain_events(inputChain, event_weight);
  write_output_to_file(arguments.outputFolder + std::string("/") + arguments.outputFileName, info);
  std::cout << "All done." << std::endl;
  return EXIT_SUCCESS;
}
