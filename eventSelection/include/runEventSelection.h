#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cassert>
#include <regex>

#include "tmArgumentParser.h"
#include "tmProgressBar.h"
#include "tmMiscellaneous.h"
#include "TROOT.h"
#include "TMath.h"
#include "TNamed.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH2I.h"

#include "constants.h"
#include "miscDataStructures.h"
#include "parameters.h"
#include "collections.h"
#include "shiftedObservablesStruct.h"
#include "STRegionsStruct.h"
#include "selectionCriteria.h"
#include "objectProperties.h"
#include "statisticsHistograms.h"
#include "MCRegions.h"
#include "selectionRegionUtils.h"
#include "hltEmulation.h"

#define N_PROBLEMATIC_ENTRIES_THRESHOLD 3

struct optionsStruct {
  std::string inputPathsFile, selectionType;
  std::vector<std::string> inputPaths;
  bool disableJetSelection, invertElectronVeto;
  int lineNumberStartInclusive, lineNumberEndInclusive, year;

  /* Not read from the command line, but instead inferred */
  bool isMC;
  std::string MC_eventProgenitor;

  friend std::ostream& operator<< (std::ostream& out, const optionsStruct& options) {
    out << "inputPathsFile: " << options.inputPathsFile << std::endl
	<< "selectionType: " << options.selectionType << std::endl
        << "disableJetSelection: " << (options.disableJetSelection? "true": "false") << std::endl
        << "Line range (for looping over paths from input file): [" << options.lineNumberStartInclusive << ", " << options.lineNumberEndInclusive << "]" << std::endl
        << "year: " << options.year << std::endl
	<< "invertElectronVeto: " << (options.invertElectronVeto? "true": "false") << std::endl
        << "isMC: " << (options.isMC? "true": "false") << std::endl
	<< "eventProgenitor: " << options.MC_eventProgenitor << std::endl;
    return out;
  }
};

optionsStruct getOptionsFromParser(tmArgumentParser& argumentParser) {
  optionsStruct options = optionsStruct();
  options.inputPathsFile = argumentParser.getArgumentString("inputPathsFile");
  std::string selectionTypeString = argumentParser.getArgumentString("selectionType");
  if (selectionTypeString == "MC_stealth_t5") {
    options.isMC = true;
    options.MC_eventProgenitor = "gluino";
  }
  else if (selectionTypeString == "MC_stealth_t6") {
    options.isMC = true;
    options.MC_eventProgenitor = "squark";
  }
  else if (selectionTypeString == "MC_EMEnrichedQCD") {
    options.isMC = false; // hack
    options.MC_eventProgenitor = "";
  }
  // else if (selectionTypeString == "MC_GJet") {
  else if (std::regex_match(selectionTypeString, std::regex("^MC_GJet[0-9]*$"))) {
    options.isMC = false; // hack
    options.MC_eventProgenitor = "";
  }
  // else if (selectionTypeString == "MC_GJet_singlephoton") {
  else if (std::regex_match(selectionTypeString, std::regex("^MC_GJet_singlephoton[0-9]*$"))) {
    options.isMC = false; // hack
    options.MC_eventProgenitor = "";
  }
  // else if (selectionTypeString == "MC_QCD") {
  else if (std::regex_match(selectionTypeString, std::regex("^MC_QCD[0-9]*$"))) {
    options.isMC = false; // hack
    options.MC_eventProgenitor = "";
  }
  // else if (selectionTypeString == "MC_QCD_singlephoton") {
  else if (std::regex_match(selectionTypeString, std::regex("^MC_QCD_singlephoton[0-9]*$"))) {
    options.isMC = false; // hack
    options.MC_eventProgenitor = "";
  }
  else if (selectionTypeString == "MC_hgg") {
    options.isMC = false; // hack
    options.MC_eventProgenitor = "";
  }
  else if (selectionTypeString == "data") {
    options.isMC = false;
    options.MC_eventProgenitor = "";
  }
  else if (selectionTypeString == "data_singlephoton") {
    options.isMC = false;
    options.MC_eventProgenitor = "";
  }
  else if (selectionTypeString == "data_jetHT") {
    options.isMC = false;
    options.MC_eventProgenitor = "";
  }
  else {
    std::cout << "ERROR: argument \"selectionType\" can only be any one of \"data\", \"data_singlephoton\", \"data_jetHT\", \"MC_stealth_t5\", \"MC_stealth_t6\", \"MC_EMEnrichedQCD\", \"MC_GJet[N]\", \"MC_GJet_singlephoton[N]\", \"MC_QCD[N]\", \"MC_QCD_singlephoton[N]\", or \"MC_hgg\"; current value: " << selectionTypeString << std::endl;
    std::exit(EXIT_FAILURE);
  }
  options.selectionType = selectionTypeString;

  std::string disableJetString = argumentParser.getArgumentString("disableJetSelection");
  if (disableJetString == "true") {
    options.disableJetSelection = true;
  }
  else if (disableJetString == "false") {
    options.disableJetSelection = false;
  }
  else {
    std::cout << "ERROR: argument \"disableJetSelection\" can be either the string \"true\" or the string \"false\"; current value: " << disableJetString << std::endl;
    std::exit(EXIT_FAILURE);
  }

  options.lineNumberStartInclusive = std::stoi(argumentParser.getArgumentString("lineNumberStartInclusive"));
  options.lineNumberEndInclusive = std::stoi(argumentParser.getArgumentString("lineNumberEndInclusive"));
  assert (options.lineNumberEndInclusive >= options.lineNumberStartInclusive);
  std::ifstream inputPathsFileStream;
  inputPathsFileStream.open(options.inputPathsFile);
  assert(inputPathsFileStream.is_open());
  int currentLineIndex = 0;
  while (!(inputPathsFileStream.eof())) {
    ++currentLineIndex;
    if (currentLineIndex > options.lineNumberEndInclusive) break;
    std::string currentLineContents;
    inputPathsFileStream >> currentLineContents;
    if (currentLineIndex >= options.lineNumberStartInclusive) {
      assert(!(currentLineContents.empty()));
      (options.inputPaths).push_back(currentLineContents);
    }
  }
  inputPathsFileStream.close();
  assert(static_cast<int>((options.inputPaths).size()) == (1 + options.lineNumberEndInclusive - options.lineNumberStartInclusive));

  options.year = std::stoi(argumentParser.getArgumentString("year"));
  if (!((options.year == 2016) || (options.year == 2017) || (options.year == 2018))) {
    std::cout << "ERROR: argument \"year\" can be one of 2016, 2017, or 2018; current value: " << options.year << std::endl;
    std::exit(EXIT_FAILURE);
  }

  std::string invertElectronVetoString = argumentParser.getArgumentString("invertElectronVeto");
  if (invertElectronVetoString == "true") {
    options.invertElectronVeto = true;
  }
  else if (invertElectronVetoString == "false") {
    options.invertElectronVeto = false;
  }
  else {
    std::cout << "ERROR: argument \"invertElectronVeto\" can be either the string \"true\" or the string \"false\"; current value: " << invertElectronVetoString << std::endl;
    std::exit(EXIT_FAILURE);
  }

  return options;
}

std::string getNDashes(const int& n) {
  std::stringstream dashes;
  for (int counter = 0; counter < n; ++counter) dashes << "-";
  return dashes.str();
}

struct MCExaminationResultsStruct{
  bool isPhotonWithDesiredMom = false;
  bool isJetCandidateFromEventProgenitor = false;
  bool isJetCandidateFromSinglet = false;
  bool isJetCandidateFromStealthSource = false;
  float eventProgenitorMass = -1.;
  float neutralinoMass = -1.;
  truthPhotonProperties truth_photon_properties;
  truthJetCandidateProperties truth_jetCandidate_properties;
};

struct photonExaminationResultsStruct{
  photonType photon_type = photonType::nPhotonTypes;
  bool isMarginallyUnselected = false;
  mediumPhotonCriterion marginallyUnselectedMediumCriterion = mediumPhotonCriterion::nMediumPhotonCriteria;
  vetoedPhotonCriterion marginallyUnselectedVetoedCriterion = vetoedPhotonCriterion::nVetoedPhotonCriteria;
  fakePhotonCriterion marginallyUnselectedFakeCriterion = fakePhotonCriterion::nFakePhotonCriteria;
  bool contributesToMisc2DHistograms = false;
  photonProperties pho_properties;
  eventWeightsStruct MCScaleFactors;
  float energy = 0.;
};

struct jetExaminationResultsStruct{
  bool passesSelectionJECNominal = false;
  bool passesSelectionJECDown = false;
  bool passesSelectionJECUp = false;
  bool passesSelectionMissingHEMDown = false;
  bool passesSelectionMissingHEMUp = false; // just for completeness...
  float missing_HEM_adjustment_pT = -1.0;
  bool contributesToHT = false; // passes all selection criteria except deltaR from nearest photon
  bool isMarginallyUnselected = false;
  bool isAwayFromCaloPhoton = false;
  bool hasGenVariablesSet = false;
  bool hasEventProgenitorPartonMom = false;
  bool hasSingletPartonMom = false;
  genJetProperties gen_jet_properties;
  jetCriterion marginallyUnselectedCriterion = jetCriterion::nJetCriteria;
  jetProperties jet_properties;
  float jecFractionalUncertainty = 0.;
  eventWeightsStruct prefireWeights;
};

std::map<shiftType, float> get_empty_STMap() {
  std::map<shiftType, float> outputMap;
  for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
    shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
    outputMap[typeIndex] = 0.;
  }
  return outputMap;
}

std::map<shiftType, int> get_empty_NJetsMap() {
  std::map<shiftType, int> outputMap;
  for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
    shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
    outputMap[typeIndex] = 0;
  }
  return outputMap;
}

struct eventExaminationResultsStruct{
  Long64_t eventIndex = 0;
  selectionRegion evt_region = selectionRegion::nSelectionRegions;
  bool isInterestingEvent = false;
  int evt_nJetsDR = 0;
  float evt_ST_electromagnetic = 0.;
  float evt_ST_hadronic = 0.;
  float evt_ST_MET = 0.;
  float evt_ST = 0.;
  float evt_photonPT_leading = 0.;
  float evt_photonPT_subLeading = 0.;
  float evt_photonEta_leading = 0.;
  float evt_photonEta_subLeading = 0.;
  float evt_photonMVA_leading = 0.;
  float evt_photonMVA_subLeading = 0.;
  eventWeightsStruct evt_prefireWeights;
  eventWeightsStruct evt_photonMCScaleFactors;
  std::map<shiftType, float> evt_shifted_ST = get_empty_STMap();
  std::map<shiftType, int> evt_shifted_nJetsDR = get_empty_NJetsMap();
};

bool passesBitMask(const UShort_t& bitCollection, const UShort_t& bitMask) {
  return ((bitCollection&bitMask) == bitMask);
}

int getMaxNJets(std::map<shiftType, int>& evt_shifted_nJetsDR) {
  int maxNJets = -1;
  for (auto&& evt_shifted_nJetsDRElement : evt_shifted_nJetsDR) {
    int evt_shifted_nJets = evt_shifted_nJetsDRElement.second;
    if ((maxNJets < 0) || (evt_shifted_nJets > maxNJets)) maxNJets = evt_shifted_nJets;
  }
  return maxNJets;
}

void addShiftedEToSTMap(const float& E, std::map<shiftType,float>& evt_shifted_ST, shiftType shift_type) {
  evt_shifted_ST[shift_type] += E;
}

void incrementNJetsMap(std::map<shiftType, int>& evt_shifted_nJetsDR, shiftType shift_type) {
  evt_shifted_nJetsDR[shift_type] += 1;
}

bool checkHLTBit(const ULong64_t& inputHLTBits, const int& indexOfBitToCheck) {
  return (((inputHLTBits>>indexOfBitToCheck)&1) == 1);
}

template<typename criterion>
int getNFalseBits(std::map<criterion, bool>& bits) {
  int nFalseBits = 0;
  for (auto&& bitsElement: bits) {
    if (!(bitsElement.second)) ++nFalseBits;
  }
  return nFalseBits;
}

template<typename criterion>
criterion getFirstFalseCriterion(std::map<criterion, bool>& bits) {
  for (auto&& bitsElement: bits) {
    if (!(bitsElement.second)) return (bitsElement.first);
  }
  // Control shouldn't reach here
  std::cout << "ERROR: getFirstFalseCriterion called with a collection of bits of which none is false." << std::endl;
  std::exit(EXIT_FAILURE);

  // Formality just to get code to compile
  return (bits.begin())->first;
}
