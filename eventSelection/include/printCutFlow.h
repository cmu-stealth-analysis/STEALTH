#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <cassert>

#include "TROOT.h"
#include "Rtypes.h"
#include "TFile.h"
#include "TH1.h"
#include "TH1D.h"

#include "selectionCriteria.h"

typedef int iyear;

struct optionsStruct {
  std::map<iyear, std::string> inputFilePaths;
  std::string outputFilePath;

  optionsStruct() {
    inputFilePaths[2016] = "";
    inputFilePaths[2017] = "";
    inputFilePaths[2018] = "";
    outputFilePath = "";
  }

  optionsStruct(int argc, char* argv[]) {
    assert(argc == 5); // first argument is name of executable
    inputFilePaths[2016] = std::string(argv[1]);
    inputFilePaths[2017] = std::string(argv[2]);
    inputFilePaths[2018] = std::string(argv[3]);
    outputFilePath = std::string(argv[4]);
  }

  friend std::ostream& operator<< (std::ostream& out, const optionsStruct& options) {
    out << "inputFilePaths: {2016: " << options.inputFilePaths.at(2016) << ", 2017: " << options.inputFilePaths.at(2017) << ", 2018: " << options.inputFilePaths.at(2018) << "}, outputFilePath: " << options.outputFilePath;
    return out;
  }
};
