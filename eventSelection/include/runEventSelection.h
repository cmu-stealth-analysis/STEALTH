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

struct optionsStruct {
  std::string inputFilesList/*, outputFilePrefix , inputFile_STRegionBoundaries */;
  bool isMC, disableJetSelection;
  /* PhotonSelectionType photonSelectionType; */
  Long64_t counterStartInclusive, counterEndInclusive;
  int year;
  /* int nGluinoMassBins, nNeutralinoMassBins; */
  /* double minGluinoMass, maxGluinoMass, minNeutralinoMass, maxNeutralinoMass; */

  friend std::ostream& operator<< (std::ostream& out, const optionsStruct& options) {
    out << "inputFilesList: " << options.inputFilesList << std::endl
        /* << "outputFilePrefix: " << options.outputFilePrefix << std::endl */
        /* << "inputFile_STRegionBoundaries: " << options.inputFile_STRegionBoundaries << std::endl */
        << "isMC: " << (options.isMC? "true": "false") << std::endl
        << "disableJetSelection: " << (options.disableJetSelection? "true": "false") << std::endl
        << "Event range: [" << options.counterStartInclusive << ", " << options.counterEndInclusive << "]" << std::endl
        << "year: " << options.year << std::endl;
        /* << "nGluinoMassBins: " << options.nGluinoMassBins << std::endl */
        /* << "nNeutralinoMassBins: " << options.nNeutralinoMassBins << std::endl */
        /* << "(minGluinoMass, maxGluinoMass): (" << options.minGluinoMass << ", " << options.maxGluinoMass << ")" << std::endl */
        /* << "(minNeutralinoMass, maxNeutralinoMass): (" << options.minNeutralinoMass << ", " << options.maxNeutralinoMass << ")" << std::endl; */
    return out;
  }
};

optionsStruct getOptionsFromParser(tmArgumentParser& argumentParser) {
  optionsStruct options = optionsStruct();
  options.inputFilesList = argumentParser.getArgumentString("inputFilesList");
  /* options.outputFilePrefix = argumentParser.getArgumentString("outputFilePrefix"); */
  // options.inputFile_STRegionBoundaries = argumentParser.getArgumentString("inputFile_STRegionBoundaries");
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

  options.counterStartInclusive = std::stol(argumentParser.getArgumentString("counterStartInclusive"));
  options.counterEndInclusive = std::stol(argumentParser.getArgumentString("counterEndInclusive"));
  options.year = std::stoi(argumentParser.getArgumentString("year"));
  if (!(options.year == 2016 || options.year == 2017)) {
    std::cout << "ERROR: argument \"year\" can be one of 2016 or 2017; current value: " << options.year << std::endl;
    std::exit(EXIT_FAILURE);
  }
  /* options.nGluinoMassBins = std::stoi(argumentParser.getArgumentString("nGluinoMassBins")); */
  /* options.minGluinoMass = std::stod(argumentParser.getArgumentString("minGluinoMass")); */
  /* options.maxGluinoMass = std::stod(argumentParser.getArgumentString("maxGluinoMass")); */
  /* options.nNeutralinoMassBins = std::stoi(argumentParser.getArgumentString("nNeutralinoMassBins")); */
  /* options.minNeutralinoMass = std::stod(argumentParser.getArgumentString("minNeutralinoMass")); */
  /* options.maxNeutralinoMass = std::stod(argumentParser.getArgumentString("maxNeutralinoMass")); */
  return options;
}

std::string getNDashes(const int& n) {
  std::stringstream dashes;
  for (int counter = 0; counter < n; ++counter) dashes << "-";
  return dashes.str();
}

struct MCExaminationResultsStruct{
  bool isPhotonWithDesiredMom = false;
  bool isJetCandidateFromGluino = false;
  bool isJetCandidateFromSinglet = false;
  bool isJetCandidateFromStealthSource = false;
  float gluinoMass = -1.;
  float neutralinoMass = -1.;
  truthPhotonProperties truth_photon_properties;
  truthJetCandidateProperties truth_jetCandidate_properties;
};

struct photonExaminationResultsStruct{
  bool isSelectedFake = false;
  bool isSelectedMedium = false;
  bool isMarginallyUnselectedFake = false;
  bool isMarginallyUnselectedMedium = false;
  mediumPhotonCriterion marginallyUnselectedMediumCriterion = mediumPhotonCriterion::nMediumPhotonCriteria;
  fakePhotonCriterion marginallyUnselectedFakeCriterion = fakePhotonCriterion::nFakePhotonCriteria;
  photonProperties pho_properties;
  eventWeightsStruct MCScaleFactors;
  float energy = 0.;
};

struct jetExaminationResultsStruct{
  bool passesSelectionJECNominal = false;
  bool passesSelectionJECDown = false;
  bool passesSelectionJECUp = false;
  bool contributesToHT = false; // passes all selection criteria except deltaR from nearest photon
  bool isMarginallyUnselected = false;
  bool isAwayFromCaloPhoton = false;
  bool hasGenVariablesSet = false;
  bool hasGluinoPartonMom = false;
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
  float evt_ST = 0.;
  float evt_photonPT_leading = 0.;
  float evt_photonPT_subLeading = 0.;
  float evt_photonEta_leading = 0.;
  float evt_photonEta_subLeading = 0.;
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
