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
  std::string inputDataPath, inputMCPath, outputFolder, outputFileName;
  bool addRelativeMCCustomWeight;
};

argumentsStruct getArgumentsFromParser(tmArgumentParser& argumentParser) {
  argumentsStruct arguments = argumentsStruct();
  arguments.inputDataPath = argumentParser.getArgumentString("inputDataPath");
  arguments.inputMCPath = argumentParser.getArgumentString("inputMCPath");
  arguments.outputFolder = argumentParser.getArgumentString("outputFolder");
  arguments.outputFileName = argumentParser.getArgumentString("outputFileName");
  std::string addRelativeMCCustomWeightRaw = argumentParser.getArgumentString("addRelativeMCCustomWeight");
  if (addRelativeMCCustomWeightRaw == "true") arguments.addRelativeMCCustomWeight = true;
  else if (addRelativeMCCustomWeightRaw == "false") arguments.addRelativeMCCustomWeight = false;
  else {
    std::cout << "ERROR: unrecognized value for argument addRelativeMCCustomWeight, needs to be \"true\" or \"false\". Currently, value: " << addRelativeMCCustomWeightRaw << std::endl;
    std::exit(EXIT_FAILURE);
  }
  
  return arguments;
}
