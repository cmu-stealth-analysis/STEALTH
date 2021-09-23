#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include "tmArgumentParser.h"
#include "tmProgressBar.h"
#include "tmMiscellaneous.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"

struct optionsStruct {
  std::string outputFolder, outputFileName;
  std::vector<std::string> inputFilePaths;
  double nJetsDistributionsMinST;
  bool useMCWeights;

  friend std::ostream& operator<< (std::ostream& out, const optionsStruct& options) {
    std::string tmp = std::string("(");
    for (const std::string & inputFilePath : (options.inputFilePaths)) tmp += (inputFilePath + std::string(", "));
    size_t tmp_length = tmp.length();
    tmp.erase(tmp_length-2, 2); /* remove trailing ", " */
    tmp += std::string(")");
    out << "inputFilePaths: " << tmp << std::endl
	<< "outputFolder: " << options.outputFolder << std::endl
	<< "outputFileName: " << options.outputFileName << std::endl
	<< "nJetsDistributionsMinST: " << options.nJetsDistributionsMinST << std::endl
	<< "useMCWeights: " << (options.useMCWeights ? "true" : "false") << std::endl;
    return out;
  }
};

optionsStruct getOptionsFromParser(tmArgumentParser& argumentParser) {
  optionsStruct options = optionsStruct();
  std::string inputFilePathsRaw = argumentParser.getArgumentString("inputFilePaths");
  (options.inputFilePaths).clear();
  options.inputFilePaths = tmMiscUtils::getSplitString(inputFilePathsRaw, std::string(";"));
  options.outputFolder = argumentParser.getArgumentString("outputFolder");
  options.outputFileName = argumentParser.getArgumentString("outputFileName");
  options.nJetsDistributionsMinST = std::stod(argumentParser.getArgumentString("nJetsDistributionsMinST"));
  std::string useMCWeightsRaw = argumentParser.getArgumentString("useMCWeights");
  if (useMCWeightsRaw == "true") options.useMCWeights = true;
  else if (useMCWeightsRaw == "false") options.useMCWeights = false;
  else {
    std::cout << "ERROR: argument \"useMCWeights\" must be set to \"true\" or \"false\". Currently it is set to: " << useMCWeightsRaw << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return options;
}
