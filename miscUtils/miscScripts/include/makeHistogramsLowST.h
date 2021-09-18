#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include "tmArgumentParser.h"
#include "tmProgressBar.h"
#include "TROOT.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

struct optionsStruct {
  std::string inputFilePath, outputFolder, outputFileName;

  friend std::ostream& operator<< (std::ostream& out, const optionsStruct& options) {
    out << "inputFilePath: " << options.inputFilePath << std::endl
	<< "outputFolder: " << options.outputFolder << std::endl
	<< "outputFileName: " << options.outputFileName << std::endl;
    return out;
  }
};

optionsStruct getOptionsFromParser(tmArgumentParser& argumentParser) {
  optionsStruct options = optionsStruct();
  options.inputFilePath = argumentParser.getArgumentString("inputFilePath");
  options.outputFolder = argumentParser.getArgumentString("outputFolder");
  options.outputFileName = argumentParser.getArgumentString("outputFileName");
  return options;
}
