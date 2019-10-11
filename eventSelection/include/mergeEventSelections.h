#include <cstdlib>
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
#include "TROOT.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"

struct optionsStruct {
  std::string inputFilesList, outputFolder, outputFileName;
  /* bool isMC; */
};

optionsStruct getOptionsFromParser(tmArgumentParser& argumentParser) {
  optionsStruct options = optionsStruct();
  options.inputFilesList = argumentParser.getArgumentString("inputFilesList");
  options.outputFolder = argumentParser.getArgumentString("outputFolder");
  options.outputFileName = argumentParser.getArgumentString("outputFileName");
  /* std::string MCString = argumentParser.getArgumentString("isMC"); */
  /* if (MCString == "true") { */
  /*   options.isMC = true; */
  /* } */
  /* else if (MCString == "false") { */
  /*   options.isMC = false; */
  /* } */
  /* else { */
  /*   std::cout << "ERROR: argument \"isMC\" can be either the string \"true\" or the string \"false\"; current value: " << MCString << std::endl; */
  /*   std::exit(EXIT_FAILURE); */
  /* } */
  return options;
}
