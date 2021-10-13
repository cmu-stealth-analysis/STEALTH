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

#include "tmProgressBar.h"
#include "tmArgumentParser.h"
#include "tmMiscellaneous.h"

#include "Rtypes.h"
#include "TROOT.h"
#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH1D.h"

#define DIPH_STRINGIFY(x) #x
#define DIPH_TH1_NAME DIPH_STRINGIFY(nJets_in_normST)
#define DIPH_STNORM_MIN 1200.0
#define DIPH_STNORM_MAX 1300.0
#define DIPH_TOLERANCE_SANITY_CHECK 0.0001

struct argumentsStruct {
  std::string outputFolder, outputFileName;
  std::vector<std::string> inputPathsFiles;
};

struct genInfoStruct {
  long long totalNEvts = 0;
  double sumWeights = 0.;
};

argumentsStruct get_command_line_arguments(int argc, char** argv) {
  tmArgumentParser argumentParser = tmArgumentParser("Save sum of weights from an input list of files to an output JSON file.");
  argumentParser.addArgument("inputPathsFiles", "", true, "Comma-separated list of paths to files containing newline-separated paths to input files with ntuples.");
  argumentParser.addArgument("outputFolder", "", true, "Folder in which to store JSON output.");
  argumentParser.addArgument("outputFileName", "", true, "Name of output file.");
  argumentParser.setPassedStringValues(argc, argv);

  argumentsStruct arguments = argumentsStruct();
  std::string inputPathsFilesRaw = argumentParser.getArgumentString("inputPathsFiles");
  (arguments.inputPathsFiles).clear();
  arguments.inputPathsFiles = tmMiscUtils::getSplitString(inputPathsFilesRaw, std::string(","));
  assert((arguments.inputPathsFiles).size() >= 1);
  arguments.outputFolder = argumentParser.getArgumentString("outputFolder");
  arguments.outputFileName = argumentParser.getArgumentString("outputFileName");
  return arguments;
}


bool file_has_zero_events(const std::string & file_path) {
  TFile *test_file = TFile::Open(file_path.c_str(), "READ");
  assert((test_file->IsOpen()) && (!(test_file->IsZombie())));
  TTree *eventTree = (TTree*)(test_file->Get("ggNtuplizer/EventTree"));
  bool has_zero_events = (eventTree == nullptr);
  test_file->Close();
  return has_zero_events;
}

void write_output_to_file(const std::string & out_file_path, const genInfoStruct & info) {
  std::ofstream output_file_handle;
  output_file_handle.open(out_file_path.c_str());
  output_file_handle << "{" << std::endl;
  output_file_handle << "    \"total_nevts_raw\": " << info.totalNEvts << "," << std::endl;
  output_file_handle << "    \"total_nevts_mc_weighted\": " << std::fixed << std::setprecision(2) << info.sumWeights << std::endl;
  output_file_handle << "}" << std::endl;
}
