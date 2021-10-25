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

#define MAX_N_EVENTS_WITHOUT_PU_INFO 10

struct argumentsStruct {
  std::string dataPUSourceWithXRDPrefix, outputFileNameWeights, outputFileNamePU;
  std::vector<std::string> inputPathsFiles;
};

struct outputInfoStruct {
  long long totalNEvts = 0;
  double sumWeights = 0.;
  TH1D pu_MC = TH1D("pileup_MC", "pileup_MC", 100, 0., 100.);
};

struct eventInfoStruct {
  float MCWeight = 0.;
  std::vector<int> * BX_for_PU = nullptr;
  std::vector<float> * PUTrue = nullptr;
};

argumentsStruct get_command_line_arguments(int argc, char** argv) {
  tmArgumentParser argumentParser = tmArgumentParser("Save sum of weights from an input list of files to an output JSON file.");
  argumentParser.addArgument("inputPathsFiles", "", true, "Comma-separated list of paths to files containing newline-separated paths to input files with ntuples.");
  argumentParser.addArgument("dataPUSourceWithXRDPrefix", "", true, "Path to ROOT file containing PU distribution of target data.");
  argumentParser.addArgument("outputFileNameWeights", "", true, "Name of output json file in which to store sum of MC weights.");
  argumentParser.addArgument("outputFileNamePU", "", true, "Name of output root file in which to store PU weights.");
  argumentParser.setPassedStringValues(argc, argv);

  argumentsStruct arguments = argumentsStruct();
  std::string inputPathsFilesRaw = argumentParser.getArgumentString("inputPathsFiles");
  (arguments.inputPathsFiles).clear();
  arguments.inputPathsFiles = tmMiscUtils::getSplitString(inputPathsFilesRaw, std::string(","));
  assert((arguments.inputPathsFiles).size() >= 1);
  arguments.dataPUSourceWithXRDPrefix = argumentParser.getArgumentString("dataPUSourceWithXRDPrefix");
  arguments.outputFileNameWeights = argumentParser.getArgumentString("outputFileNameWeights");
  arguments.outputFileNamePU = argumentParser.getArgumentString("outputFileNamePU");
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

void write_weight_outputs_to_json_file(const std::string & out_file_path, outputInfoStruct & output_info) {
  std::ofstream output_file_handle;
  output_file_handle.open(out_file_path.c_str());
  output_file_handle << "{" << std::endl;
  output_file_handle << "    \"total_nevts_raw\": " << output_info.totalNEvts << "," << std::endl;
  output_file_handle << "    \"total_nevts_mc_weighted\": " << std::fixed << std::setprecision(2) << output_info.sumWeights << std::endl;
  output_file_handle << "}" << std::endl;
}

void write_histogram_outputs_to_root_file(const std::string & out_file_path, TH1D * pileup_MC, TH1D * PUWeightsHistogram) {
  TFile *outputFileHandle = TFile::Open(out_file_path.c_str(), "RECREATE");
  assert((outputFileHandle->IsOpen()) && (!(outputFileHandle->IsZombie())));
  outputFileHandle->WriteTObject(pileup_MC);
  outputFileHandle->WriteTObject(PUWeightsHistogram);
  outputFileHandle->Close();
}
