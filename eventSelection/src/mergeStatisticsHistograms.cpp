#include "../include/mergeStatisticsHistograms.h"

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  tmArgumentParser argumentParser = tmArgumentParser("Calculate systematics due to uncertainty on jet energy corrections.");
  argumentParser.addArgument("inputFilesList", "", true, "Path to file containing list of paths with histograms.");
  argumentParser.addArgument("outputFolder", "root://cmseos.fnal.gov//store/user/lpcsusystealth/statistics/merged/", false, "Output folder.");
  argumentParser.addArgument("outputFileName", "statistics.root", false, "Name of output file.");
  argumentParser.addArgument("isMC", "false", true, "Takes value \"true\" if there are additional plots relevant for MC samples only.");
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);
  std::cout << "Beginning to merge statistics histograms..." << std::endl;
  statisticsHistograms histogramsList = statisticsHistograms(options.isMC);

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
  
  for (auto&& inputFileName: inputFileNames) {
    std::cout << "Adding histograms from file: " << inputFileName << std::endl;
    TFile *inputFile = TFile::Open(inputFileName.c_str(), "READ");
    if (!(inputFile->IsOpen()) || inputFile->IsZombie()) {
      std::cout << "ERROR: Unable to open input file to read from. Attempted to read from file with path: " << inputFileName << std::endl;
    }
    for (auto&& statsElement: (histogramsList.stats)) {
      TH1F *inputHistogram;
      auto& histogramName = statsElement.first;
      inputFile->GetObject(histogramName.c_str(), inputHistogram);
      if (inputHistogram) {
        (statsElement.second).Add(inputHistogram);
      }
      else {
        std::cout << "ERROR: Unable to find histogram named \"" << histogramName << "\"." << std::endl;
      }
    }
    inputFile->Close();
  }

  histogramsList.writeToFile(options.outputFileName);
  return 0;
}
