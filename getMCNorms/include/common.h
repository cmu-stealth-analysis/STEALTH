#ifndef GETMCNORMS_INCLUDE_COMMON_H
#define GETMCNORMS_INCLUDE_COMMON_H

#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cassert>
#include <map>
#include <vector>

#include "TH1D.h"
#include "TTree.h"
#include "TChain.h"
#include "TFile.h"

#include "tmArgumentParser.h"
#include "tmMiscellaneous.h"

#include "../../eventSelection/include/STRegionsStruct.h"

#define ST_MAX_RANGE 3500.0

namespace common {
  struct argumentsStruct {
    std::string outputFolder, outputFileName;
    std::vector<std::string> inputPathsFiles;
    bool addMCWeights;
    STRegionsStruct STRegions;
  };

  argumentsStruct get_command_line_arguments(int, char**);
  bool file_has_zero_events(const std::string &);
  TChain * get_chain_from_input_paths_files(const std::vector<std::string> &);
  void write_output_th1s_to_file(const std::string &, std::map<std::string, TH1D> &);
  void move_via_xrdcp(const std::string &, const std::string &);
}

#endif
