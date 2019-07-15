#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include "tmArgumentParser.h"
#include "TROOT.h"
#include "TFile.h"
#include "TObject.h"
#include "TNamed.h"

#include "selectionCriteria.h"
#include "TH2I.h"
#include "TH2F.h"
#include "TFile.h"
#include "parameters.h"

struct optionsStruct {
  std::string inputFileName, outputFolder, outputFileName;
  bool isMC;
};

optionsStruct getOptionsFromParser(tmArgumentParser& argumentParser) {
  optionsStruct options = optionsStruct();
  options.inputFileName = argumentParser.getArgumentString("inputFileName");
  options.outputFolder = argumentParser.getArgumentString("outputFolder");
  options.outputFileName = argumentParser.getArgumentString("outputFileName");
  return options;
}
