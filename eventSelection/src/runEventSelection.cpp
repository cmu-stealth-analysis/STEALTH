#include "../include/runEventSelection.h"

optionsStruct getOptionsFromParser(tmArgumentParser& argumentParser) {
  optionsStruct options = optionsStruct();
  options.inputFilesList = argumentParser.getArgumentString("inputFilesList");
  options.outputFilePath = argumentParser.getArgumentString("outputFilePath");
  std::string MCString = argumentParser.getArgumentString("isMC");
  if (MCString == "true") {
    options.isMC = true;
  }
  else if (MCString == "false") {
    options.isMC = false;
  }
  else {
    std::cout << "ERROR: argument \"isMC\" can be either the string \"true\" or the string \"false\"; current value: " << MCString << std::endl;
    std::exit(EXIT_FAILURE);
  }
  options.counterStartInclusive = std::stol(argumentParser.getArgumentString("counterStartInclusive"));
  options.counterEndInclusive = std::stol(argumentParser.getArgumentString("counterEndInclusive"));
  std::string photonSelectionTypeString = argumentParser.getArgumentString("photonSelectionType");
  if (photonSelectionTypeString == "medium") {
    options.photonSelectionType = PhotonSelectionType::medium;
  }
  else if (photonSelectionTypeString == "fake") {
    options.photonSelectionType = PhotonSelectionType::fake;
  }
  else if (photonSelectionTypeString == "mediumfake") {
    options.photonSelectionType = PhotonSelectionType::mediumfake;
  }
  else {
    std::cout << "ERROR: argument \"photonSelectionType\" can be one of \"medium\", \"fake\", or \"mediumfake\"; current value: " << photonSelectionTypeString << std::endl;
    std::exit(EXIT_FAILURE);
  }
  options.year = std::stoi(argumentParser.getArgumentString("year"));
  if (!(options.year == 2016 || options.year == 2017 || options.year == -1)) {
    std::cout << "ERROR: argument \"year\" can be one of 2016, 2017, or -1; current value: " << options.year << std::endl;
    std::exit(EXIT_FAILURE);
  }
  options.JECUncertainty = std::stoi(argumentParser.getArgumentString("JECUncertainty"));
  if (std::abs(options.JECUncertainty) > 1) {
    std::cout << "ERROR: argument \"JECUncertainty\" can be one of -1, 0, or +1; current value: " << options.JECUncertainty << std::endl;
    std::exit(EXIT_FAILURE);
  }
  return options;
}

std::string getNDashes(const int& n) {
  std::stringstream dashes;
  for (int counter = 0; counter < n; ++counter) dashes << "-";
  return dashes.str();
}

void initializeCounters(countersStruct &counters) {
  for (int categoryIndex = photonFailureCategoryFirst; categoryIndex != static_cast<int>(photonFailureCategory::nPhotonFailureCategories); ++categoryIndex) {
    photonFailureCategory category = static_cast<photonFailureCategory>(categoryIndex);
    counters.photonFailureCounters[category] = 0l;
  }

  for (int categoryIndex = jetFailureCategoryFirst; categoryIndex != static_cast<int>(jetFailureCategory::nJetFailureCategories); ++categoryIndex) {
    jetFailureCategory category = static_cast<jetFailureCategory>(categoryIndex);
    counters.jetFailureCounters[category] = 0l;
  }

  for (int categoryIndex = eventFailureCategoryFirst; categoryIndex != static_cast<int>(eventFailureCategory::nEventFailureCategories); ++categoryIndex) {
    eventFailureCategory category = static_cast<eventFailureCategory>(categoryIndex);
    counters.eventFailureCounters[category] = 0l;
  }

  for (int miscCounterIndex = miscCounterFirst; miscCounterIndex != static_cast<int>(miscCounter::nMiscCounters); ++miscCounterIndex) {
    miscCounter miscCounterEnumIndex = static_cast<miscCounter>(miscCounterIndex);
    counters.miscCounters[miscCounterEnumIndex] = 0l;
  }
}

void printCounters(countersStruct &counters) {
  std::cout << "Photon counters: " << std::endl;
  for (const auto& counterValuePair : counters.photonFailureCounters) {
    std::cout << photonFailureCategoryNames[counterValuePair.first] << " : " << counterValuePair.second << std::endl;
  }
  std::cout << getNDashes(100) << std::endl;

  std::cout << "Jet counters: " << std::endl;
  for (const auto& counterValuePair : counters.jetFailureCounters) {
    std::cout << jetFailureCategoryNames[counterValuePair.first] << " : " << counterValuePair.second << std::endl;
  }
  std::cout << getNDashes(100) << std::endl;

  std::cout << "Event counters: " << std::endl;
  for (const auto& counterValuePair : counters.eventFailureCounters) {
    std::cout << eventFailureCategoryNames[counterValuePair.first] << " : " << counterValuePair.second << std::endl;
  }
  std::cout << getNDashes(100) << std::endl;

  std::cout << "Miscellaneous counters: " << std::endl;
  for (const auto& counterValuePair : counters.miscCounters) {
    std::cout << miscCounterNames[counterValuePair.first] << " : " << counterValuePair.second << std::endl;
  }
}

void incrementCounters(photonFailureCategory photonCategory, countersStruct& counters) {
  ++((counters.photonFailureCounters)[photonCategory]);
}

void incrementCounters(jetFailureCategory jetCategory, countersStruct& counters) {
  ++((counters.jetFailureCounters)[jetCategory]);
}

void incrementCounters(eventFailureCategory eventCategory, countersStruct& counters) {
  ++((counters.eventFailureCounters)[eventCategory]);
}

void incrementCounters(miscCounter miscCounterEnumIndex, countersStruct& counters) {
  ++((counters.miscCounters)[miscCounterEnumIndex]);
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  tmArgumentParser argumentParser = tmArgumentParser("Calculate systematics due to uncertainty on jet energy corrections.");
  argumentParser.addArgument("inputFilesList", "", true, "Path to file containing list of input files.");
  argumentParser.addArgument("outputFilePath", "", true, "Path to output file.");
  argumentParser.addArgument("isMC", "false", false, "Input file is a MC sample -- disable HLT photon trigger and enable additional MC selection.");
  argumentParser.addArgument("counterStartInclusive", "", true, "Event number from input file from which to start. The event with this index is included in the processing.");
  argumentParser.addArgument("counterEndInclusive", "", true, "Event number from input file at which to end. The event with this index is included in the processing.");
  argumentParser.addArgument("photonSelectionType", "fake", true, "Photon selection type: can be any one of: \"fake\", \"medium\", \"mediumfake\", \"fakeMC\", \"mediumMC\", \"mediumfakeMC\"");
  argumentParser.addArgument("year", "2017", false, "Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.");
  argumentParser.addArgument("JECUncertainty", "0", false, "Apply a uniform upward or downward jet energy uncertainty correction to jet pt. Default: 0, i.e. do not apply any other correction. +/-1 are allowed as well, shifting all jet pt up or down respectively by 1.0 times the uncertainty on the jet energy correction.");
  argumentParser.setPassedStringValues(argc, argv);

  optionsStruct options = getOptionsFromParser(argumentParser);

  parametersStruct parameters = parametersStruct();
  parameters.tuneParametersForYear(options.year);

  countersStruct counters = countersStruct();
  initializeCounters(counters);

  std::stringstream optionsStringstream;
  optionsStringstream << options;
  TNamed *optionsObject = new TNamed("optionsString", optionsStringstream.str().c_str());
  std::stringstream parametersStringstream;
  parametersStringstream << parameters;
  TNamed *parametersObject = new TNamed("parametersString", parametersStringstream.str().c_str());
  TFile *outputFile = TFile::Open(options.outputFilePath.c_str(), "RECREATE");
  if (!(outputFile->IsOpen()) || outputFile->IsZombie()) {
    std::cout << "ERROR: Unable to open output file to write. File path: " << options.outputFilePath << std::endl;
  }

  outputFile->WriteTObject(parametersObject);
  outputFile->WriteTObject(optionsObject);
  outputFile->Close();

  std::cout << getNDashes(100) << std::endl;
  // For testing: first print counters, should have all zeros
  std::cout << "Empty counters:" << std::endl;
  printCounters(counters);
  incrementCounters(photonFailureCategory::pT, counters);
  incrementCounters(photonFailureCategory::pT, counters);
  incrementCounters(jetFailureCategory::eta, counters);
  incrementCounters(eventFailureCategory::MCGenInformation, counters);
  incrementCounters(eventFailureCategory::MCGenInformation, counters);
  incrementCounters(eventFailureCategory::MCGenInformation, counters);
  incrementCounters(miscCounter::acceptedEvents, counters);
  // Should now have 1 in these categories
  printCounters(counters);
  std::cout << getNDashes(100) << std::endl
            << "Options:" << std::endl
            << optionsStringstream.str() << std::endl
            << getNDashes(100) << std::endl;
  std::cout << getNDashes(100) << std::endl
            << "Parameters:" << std::endl
            << parametersStringstream.str() << std::endl
            << getNDashes(100) << std::endl;
  return 0;
}
