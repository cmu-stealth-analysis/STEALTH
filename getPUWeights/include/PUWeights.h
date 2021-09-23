#include <cstdlib>
#include <algorithm>
#include <initializer_list>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cassert>

#include "tmArgumentParser.h"
#include "tmProgressBar.h"
#include "tmMiscellaneous.h"
#include "tmROOTSaverUtils.h"
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TFile.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH1D.h"
#include "TObjArray.h"

struct argumentsStruct {
  std::string inputDataPath, outputFolder, outputFileName;
  std::vector<std::string> inputMCPaths;
  bool addMCXSecWeight;
};

argumentsStruct getArgumentsFromParser(tmArgumentParser& argumentParser) {
  argumentsStruct arguments = argumentsStruct();
  arguments.inputDataPath = argumentParser.getArgumentString("inputDataPath");
  std::string inputMCPathsRaw = argumentParser.getArgumentString("inputMCPaths");
  (arguments.inputMCPaths).clear();
  arguments.inputMCPaths = tmMiscUtils::getSplitString(inputMCPathsRaw, std::string(","));
  assert((arguments.inputMCPaths).size() >= 1);
  arguments.outputFolder = argumentParser.getArgumentString("outputFolder");
  arguments.outputFileName = argumentParser.getArgumentString("outputFileName");
  std::string addMCXSecWeightRaw = argumentParser.getArgumentString("addMCXSecWeight");
  if (addMCXSecWeightRaw == "true") arguments.addMCXSecWeight = true;
  else if (addMCXSecWeightRaw == "false") arguments.addMCXSecWeight = false;
  else {
    std::cout << "ERROR: unrecognized value for argument addMCXSecWeight, needs to be \"true\" or \"false\". Currently, value: " << addMCXSecWeightRaw << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  return arguments;
}
