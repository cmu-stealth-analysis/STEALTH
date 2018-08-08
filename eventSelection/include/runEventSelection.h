#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

#include "tmArgumentParser.h"
#include "tmProgressBar.h"
#include "tmMiscellaneous.h"
/* #include "tmROOTSaverUtils.h" */
#include "TROOT.h"
#include "TTree.h"
/* #include "TTreeReader.h" */
/* #include "TTreeReaderValue.h" */
/* #include "TTreeReaderArray.h" */
#include "TFile.h"
#include "TNamed.h"

struct PIDsStruct {
  const int photon = 22;
  const int gluino = 1000021;
  const int neutralino = 1000022;
  friend std::ostream& operator<< (std::ostream& out, const PIDsStruct& PIDs) {
    out << "photon --> " << PIDs.photon << ", "
        << "gluino --> " << PIDs.gluino << ", "
        << "neutralino --> " << PIDs.neutralino;
    return out;
  }
};

struct rangeStruct{
  float rangeLower, rangeUpper;
rangeStruct(float rangeLower_, float rangeUpper_) : rangeLower(rangeLower_),
    rangeUpper(rangeUpper_){}
  bool isInside(float candidate) {
    return (candidate >= rangeLower && candidate < rangeUpper);
  }
  friend std::ostream& operator<< (std::ostream& out, const rangeStruct& range) {
    out << "[" << range.rangeLower << ", " << range.rangeUpper << ")";
    return out;
  }
};

struct quadraticPolynomialStruct{
  float constCoefficient, linearCoefficient, squareCoefficient;
quadraticPolynomialStruct(float constCoefficient_, float linearCoefficient_, float squareCoefficient_) : constCoefficient(constCoefficient_),
    linearCoefficient(linearCoefficient_),
    squareCoefficient(squareCoefficient_) {}
  float getPolynomialValue(float pT) {
    return (constCoefficient + linearCoefficient*pT + squareCoefficient*pT*pT);
  }
  friend std::ostream& operator<< (std::ostream& out, const quadraticPolynomialStruct& polynomial) {
    out << "const --> " << polynomial.constCoefficient << ", "
        << "linear--> " << polynomial.linearCoefficient << ", "
        << "square --> " << polynomial.squareCoefficient;
    return out;
  }
};

struct EAValuesStruct{
  float regionUpperBound, chargedHadronsEA, neutralHadronsEA, photonsEA;
EAValuesStruct(float regionUpperBound_, float chargedHadronsEA_, float neutralHadronsEA_, float photonsEA_) : regionUpperBound(regionUpperBound_),
    chargedHadronsEA(chargedHadronsEA_),
    neutralHadronsEA(neutralHadronsEA_),
    photonsEA(photonsEA_) {}

  friend std::ostream& operator<< (std::ostream& out, const EAValuesStruct& EAValues) {
    out << "region upper eta bound --> " << EAValues.regionUpperBound << ", effective areas: "
        << "charged hadrons --> " << EAValues.chargedHadronsEA << ", "
        << "neutral hadrons --> " << EAValues.neutralHadronsEA << ", "
        << "photons --> " << EAValues.photonsEA;
    return out;
  }
};

struct parametersStruct {
  const float pTCutSubLeading = 25.0;
  const float pTCutLeading = 35.0;
  const float photonEtaCut = 1.442;
  const float R9Cut = 1.0;
  const int nSubLeadingPhotons = 2;
  const int nLeadingPhotons = 1;
  const float jetEtaCut = 2.4;
  const float jetpTCut = 30.;
  const float jetPUIDThreshold = 0.61;
  const int jetSelectionID = 6;
  const float minDeltaRCut = 0.4;
  const int nJetsCut = 2;
  const float HTCut = 60.;
  const float electronPtCut = 15.;
  const float electronEtaCut = 2.5;
  const float electronDzCut = 0.1;
  const float electronPFPUIsoCut = 0.1;
  const int nElectronsCut = 0;
  const float muonPtCut = 15.0;
  const float muonPFPUIsoCut = 0.12;
  const int nMuonsCut = 0;
  const float METThreshold = 15.0;
  const EAValuesStruct region1EAs = EAValuesStruct(1.0, 0.0385, 0.0636, 0.124);
  const EAValuesStruct region2EAs = EAValuesStruct(1.479, 0.0468, 0.1103, 0.1093);
  const PIDsStruct PIDs;
  const int nPhotonsWithNeutralinoMom = 2;
  const float towerHOverECut = 0.035;
  const rangeStruct sigmaietaietaRange = rangeStruct(0.0103, 0.015);
  const rangeStruct chargedIsolationRange = rangeStruct(1.416, 15.0);
  const quadraticPolynomialStruct neutralIsolationCut = quadraticPolynomialStruct(2.491, 0.0126, 0.000026);
  const quadraticPolynomialStruct photonIsolationCut = quadraticPolynomialStruct(2.952, 0.004, 0.0);
  int HLTPhotonBit = -1;
  const int HLTPhotonBit2016 = 14;
  const int HLTPhotonBit2017 = 36;
  void tuneParametersForYear(int year) {
    if (year == 2017) HLTPhotonBit = HLTPhotonBit2017;
    else if (year == 2016) HLTPhotonBit = HLTPhotonBit2016;
  }
  friend std::ostream& operator<< (std::ostream& out, const parametersStruct& parameters) {
    out << "PIDs: " << parameters.PIDs << std::endl;

    out << "Photon cuts:" << std::endl
        << "pT_SubLeading: " << parameters.pTCutSubLeading << ", "
        << "pT_Leading: " << parameters.pTCutLeading << ", "
        << "eta: " << parameters.photonEtaCut << ", "
        << "R9: " << parameters.R9Cut << ", "
        << "HOverE: " << parameters.towerHOverECut << ", "
        << "sigmaietaietaRange: " << parameters.sigmaietaietaRange << ", "
        << "chargedIsolationRange: " << parameters.chargedIsolationRange << ", "
        << "neutral isolation cut coefficients: " << parameters.neutralIsolationCut << ", "
        << "photon isolation cut coefficients: " << parameters.photonIsolationCut << std::endl;

    out << "Region 1 effective areas: " << parameters.region1EAs << std::endl
        << "Region 2 effective areas: " << parameters.region2EAs << std::endl;

    out << "Jet cuts:" << std::endl
        << "pT: " << parameters.jetpTCut << ", "
        << "eta: " << parameters.jetEtaCut << ", "
        << "PUID: " << parameters.jetPUIDThreshold << ", "
        << "jet ID: " << parameters.jetSelectionID << ", "
        << "minDeltaR: " << parameters.minDeltaRCut << std::endl;

    out << "Electron cuts:" << std::endl
        << "pT: " << parameters.electronPtCut << ", "
        << "eta: " << parameters.electronEtaCut << ", "
        << "Dz: " << parameters.electronDzCut << ", "
        << "PFPUIso: " << parameters.electronPFPUIsoCut << std::endl;

    out << "Muon cuts:" << std::endl
        << "pT: " << parameters.muonPtCut << ", "
        << "PFPUIso: " << parameters.muonPFPUIsoCut << std::endl;

    out << "MET cuts:" << std::endl
        << "MET: " << parameters.METThreshold << std::endl;

    out << "Event cuts:" << std::endl
        << "nSubLeading photons: " << parameters.nSubLeadingPhotons << ", "
        << "nLeading photons min: " << parameters.nLeadingPhotons << ", "
        << "photon HLT bit index: " << parameters.HLTPhotonBit << ", "
        << "nJets min: " << parameters.nJetsCut << ", "
        << "HT Cut: " << parameters.HTCut << ", "
        << "nElectrons: " << parameters.nElectronsCut << ", "
        << "nMuons: " << parameters.nMuonsCut << ", "
        << "nPhotonsWithNeutralinoMom: " << parameters.nPhotonsWithNeutralinoMom;
    return out;
  }
};

enum class PhotonSelectionType{medium, fake, mediumfake};
std::string getPhotonSelectionTypeString(PhotonSelectionType type) {
  std::string outputString;
  switch(type) {
  case (PhotonSelectionType::medium) :
    outputString = "PhotonSelectionType::medium";
    break;
  case (PhotonSelectionType::fake) :
    outputString = "PhotonSelectionType::fake";
    break;
  case (PhotonSelectionType::mediumfake) :
    outputString = "PhotonSelectionType::mediumfake";
    break;
  default:
    std::cout << "ERROR: Unknown photon selection type!"<< std::endl;
    std::exit(EXIT_FAILURE);
  }
  return outputString;
}

struct optionsStruct {
  std::string inputFilesList, outputFilePath;
  bool isMC;
  PhotonSelectionType photonSelectionType;
  long counterStartInclusive, counterEndInclusive;
  int year, JECUncertainty;

  friend std::ostream& operator<< (std::ostream& out, const optionsStruct& options) {
    out << "inputFilesList: " << options.inputFilesList << std::endl
        << "outputFilePath: " << options.outputFilePath << std::endl
        << "isMC: " << (options.isMC? "true": "false") << ", "
        << "photonSelectionType: " << getPhotonSelectionTypeString(options.photonSelectionType) << std::endl
        << "Event range: [" << options.counterStartInclusive << ", " << options.counterEndInclusive << "]" << std::endl
        << "year: " << options.year << ", "
        << "JEC Uncertainty scale: " << options.JECUncertainty;
    return out;
  }
};

enum class photonFailureCategory{eta=0, pT, hOverE, neutralIsolation, photonIsolation, conversionSafeElectronVeto, R9, sigmaietaiataXORchargedIsolation, mediumIDCut, nPhotonFailureCategories};
int photonFailureCategoryFirst = static_cast<int>(photonFailureCategory::eta);
std::map<photonFailureCategory, std::string> photonFailureCategoryNames = {
  {photonFailureCategory::eta, "eta"},
  {photonFailureCategory::pT, "pT"},
  {photonFailureCategory::hOverE, "hOverE"},
  {photonFailureCategory::neutralIsolation, "neutralIsolation"},
  {photonFailureCategory::photonIsolation, "photonIsolation"},
  {photonFailureCategory::conversionSafeElectronVeto, "conversionSafeElectronVeto"},
  {photonFailureCategory::R9, "R9"},
  {photonFailureCategory::sigmaietaiataXORchargedIsolation, "sigmaietaiataXORchargedIsolation"},
  {photonFailureCategory::mediumIDCut, "mediumIDCut"}
};

enum class jetFailureCategory{eta=0, pT, PFLooseID, puID, jetID, nJetFailureCategories};
int jetFailureCategoryFirst = static_cast<int>(jetFailureCategory::eta);
std::map<jetFailureCategory, std::string> jetFailureCategoryNames = {
  {jetFailureCategory::eta, "eta"},
  {jetFailureCategory::pT, "pT"},
  {jetFailureCategory::PFLooseID, "PFLooseID"},
  {jetFailureCategory::puID, "puID"},
  {jetFailureCategory::jetID, "jetID"}
};

enum class eventFailureCategory{HLTPhoton=0, wrongNMediumOrFakePhotons, wrongNPhotons, HLTJet, wrongNJets, hTCut, electronVeto, muonVeto, MCGenInformation, nEventFailureCategories};
int eventFailureCategoryFirst = static_cast<int>(eventFailureCategory::HLTPhoton);
std::map<eventFailureCategory, std::string> eventFailureCategoryNames = {
  {eventFailureCategory::HLTPhoton, "HLTPhoton"},
  {eventFailureCategory::wrongNMediumOrFakePhotons, "wrongNMediumOrFakePhotons"},
  {eventFailureCategory::wrongNPhotons, "wrongNPhotons"},
  {eventFailureCategory::HLTJet, "HLTJet"},
  {eventFailureCategory::wrongNJets, "wrongNJets"},
  {eventFailureCategory::hTCut, "hTCut"},
  {eventFailureCategory::electronVeto, "electronVeto"},
  {eventFailureCategory::muonVeto, "muonVeto"},
  {eventFailureCategory::MCGenInformation, "MCGenInformation"}
};

enum class miscCounter{failingPhotons=0, failingJets, failingEvents, acceptedEvents, nMiscCounters};
int miscCounterFirst = static_cast<int>(miscCounter::failingPhotons);
std::map<miscCounter, std::string> miscCounterNames = {
  {miscCounter::failingPhotons, "failingPhotons"},
  {miscCounter::failingJets, "failingJets"},
  {miscCounter::failingEvents, "failingEvents"},
  {miscCounter::acceptedEvents, "acceptedEvents"}
};

struct countersStruct{
  std::map<photonFailureCategory, long> photonFailureCounters;
  std::map<jetFailureCategory, long> jetFailureCounters;
  std::map<eventFailureCategory, long> eventFailureCounters;
  std::map<miscCounter, long> miscCounters;
};
