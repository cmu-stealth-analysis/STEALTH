#include "../include/common.h"

common::argumentsStruct
common::get_command_line_arguments(int argc, char** argv) {
  tmArgumentParser argumentParser = tmArgumentParser("Save various distributions from input files.");
  argumentParser.addArgument("inputPathsFiles", "", true, "Comma-separated list of paths to files containing newline-separated paths to input files with ntuples.");
  argumentParser.addArgument("outputFolder", "root://cmseos.fnal.gov//store/user/lpcsusystealth/analysisEOSAreas/analysis", false, "Output folder.");
  argumentParser.addArgument("outputFileName", "", true, "Name of output file.");
  argumentParser.addArgument("addMCWeights", "false", true, "If this argument is set, then relative weights are read in from an additional branch, used for GJet MC samples.");
  argumentParser.setPassedStringValues(argc, argv);

  argumentsStruct arguments = argumentsStruct();
  std::string inputPathsFilesRaw = argumentParser.getArgumentString("inputPathsFiles");
  (arguments.inputPathsFiles).clear();
  arguments.inputPathsFiles = tmMiscUtils::getSplitString(inputPathsFilesRaw, std::string(","));
  assert((arguments.inputPathsFiles).size() >= 1);
  arguments.outputFolder = argumentParser.getArgumentString("outputFolder");
  arguments.outputFileName = argumentParser.getArgumentString("outputFileName");
  std::string addMCWeightsRaw = argumentParser.getArgumentString("addMCWeights");
  if (addMCWeightsRaw == "true") arguments.addMCWeights = true;
  else if (addMCWeightsRaw == "false") arguments.addMCWeights = false;
  else {
    std::cout << "ERROR: unrecognized value for argument addMCWeights, needs to be \"true\" or \"false\". Currently, value: " << addMCWeightsRaw << std::endl;
    std::exit(EXIT_FAILURE);
  }

  return arguments;
}

bool
common::file_has_zero_events(const std::string & file_path) {
  TFile *test_file = TFile::Open(file_path.c_str(), "READ");
  assert((test_file->IsOpen()) && (!(test_file->IsZombie())));
  TTree *eventTree = (TTree*)(test_file->Get("ggNtuplizer/EventTree"));
  bool has_zero_events = (eventTree == nullptr);
  test_file->Close();
  return has_zero_events;
}

TChain *
common::get_chain_from_input_paths_files(const std::vector<std::string> & inputPathsFiles) {
  TChain * inputChain = new TChain("ggNtuplizer/EventTree");
  inputChain->SetMaxTreeSize(10000000000000LL); // 10 TB
  for (const std::string & inputPathsFile : inputPathsFiles) {
    std::cout << "Adding paths from file: " << inputPathsFile << std::endl;
    std::ifstream inputPathsFileStream;
    inputPathsFileStream.open(inputPathsFile.c_str());
    assert(inputPathsFileStream.is_open());
    while (!(inputPathsFileStream.eof())) {
      std::string inputPath;
      inputPathsFileStream >> inputPath;
      if (!(inputPath.empty())) {
	if (!(common::file_has_zero_events(inputPath))) inputChain->Add(inputPath.c_str());
      }
    }
    inputPathsFileStream.close();
  }
  return inputChain;
}

void
common::write_output_th1s_to_file(const std::string & out_file_path, std::map<std::string, TH1D> & output_th1s) {
  TFile * outputFile = TFile::Open(out_file_path.c_str(), "RECREATE");
  assert((outputFile->IsOpen()) && (!(outputFile->IsZombie())));
  for (auto & name_histogram_pair : output_th1s) {
    outputFile->WriteTObject(&(name_histogram_pair.second));
  }
  outputFile->Close();
}

void
common::move_via_xrdcp(const std::string & source, const std::string & target) {
  int xrdcp_return_status = system(("set -x && xrdcp --nopbar --silent --force --path --streams 15 " + source + " " + target + " && rm -f " + source + " && set +x").c_str());
  assert (xrdcp_return_status == 0);
}
