#include <cstdlib>
#include <cmath>
#include <iostream>
#include <vector>
#include <map>
#include <cassert>

#include "TROOT.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"

#include "tmArgumentParser.h"
#include "tmProgressBar.h"

#include "../../../eventSelection/include/STRegionsStruct.h"

struct optionsStruct {
  std::string sourceFilePath, outputFolder, selection, yearString;
  double STNormTarget, STNormMax;
  STRegionsStruct STRegions;

  friend std::ostream& operator<< (std::ostream& out, const optionsStruct& options) {
    out << "sourceFilePath: " << options.sourceFilePath << std::endl
	<< "outputFolder: " << options.outputFolder << std::endl
        << "selection: " << options.selection << std::endl
        << "yearString: " << options.yearString << std::endl
        << "STRegions: " << options.STRegions << std::endl
        << "STNormTarget: " << options.STNormTarget << std::endl
        << "STNormMax: " << options.STNormMax << std::endl;
    return out;
  }
};

optionsStruct getOptionsFromParser(tmArgumentParser& argumentParser) {
  optionsStruct options = optionsStruct();
  options.sourceFilePath = argumentParser.getArgumentString("sourceFilePath");
  options.outputFolder = argumentParser.getArgumentString("outputFolder");
  options.selection = argumentParser.getArgumentString("selection");
  options.yearString = argumentParser.getArgumentString("yearString");
  std::string STBoundariesSourceFile = argumentParser.getArgumentString("STBoundariesSourceFile");
  options.STRegions = STRegionsStruct(STBoundariesSourceFile, 3500.0);
  options.STNormTarget = std::stod(argumentParser.getArgumentString("STNormTarget"));
  options.STNormMax = std::stod(argumentParser.getArgumentString("STNormMax"));
  return options;
}
