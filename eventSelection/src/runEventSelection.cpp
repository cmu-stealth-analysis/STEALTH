#include "../include/runEventSelection.h"

optionsStruct getOptionsFromParser(tmArgumentParser& argumentParser) {
  optionsStruct options = optionsStruct();
  // options.inputMCPath = argumentParser.getArgumentString("inputMCPath");
  options.inputFilePath = argumentParser.getArgumentString("inputFilePath");
  options.outputFilePath = argumentParser.getArgumentString("outputFilePath");
  options.counterStartInclusive = std::stol(argumentParser.getArgumentString("counterStartInclusive"));
  options.counterEndInclusive = std::stol(argumentParser.getArgumentString("counterEndInclusive"));
  options.photonSelectionType = argumentParser.getArgumentString("photonSelectionType");
  options.year = std::stoi(argumentParser.getArgumentString("year"));
  options.JECUncertainty = std::stoi(argumentParser.getArgumentString("JECUncertainty"));
  return options;
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  parametersStruct parameters = parametersStruct();
  tmArgumentParser argumentParser = tmArgumentParser("Calculate systematics due to uncertainty on jet energy corrections.");
  argumentParser.addArgument("inputFilePath", "", true, "Path to input file.");
  argumentParser.addArgument("outputFilePath", "", true, "Path to output file.");
  argumentParser.addArgument("counterStartInclusive", "", true, "Event number from input file from which to start. The event with this index is included in the processing.");
  argumentParser.addArgument("counterEndInclusive", "", true, "Event number from input file at which to end. The event with this index is included in the processing.");
  argumentParser.addArgument("photonSelectionType", "fake", true, "Photon selection type: can be any one of: \"fake\", \"medium\", \"mediumfake\", \"fakeMC\", \"mediumMC\", \"mediumfakeMC\"");
  argumentParser.addArgument("year", "-1", false, "Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger, and the photon ID cuts which are based on year-dependent recommendations. Default year: -1, which is used for MC and means the trigger is disabled and the 2017 photon ID recommendations are implemented.");
  argumentParser.addArgument("JECUncertainty", "0", false, "Apply a uniform upward or downward jet energy uncertainty correction to jet pt. Default: 0, i.e. do not apply any other correction. +/-1 are allowed as well, shifting all jet pt up or down respectively by the relevant jet energy correction.");
  argumentParser.setPassedStringValues(argc, argv);
  optionsStruct options = getOptionsFromParser(argumentParser);
  parameters.setHLTPhotonBitForYear(options.year);
  return 0;
}
