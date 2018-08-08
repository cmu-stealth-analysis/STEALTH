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
#include "TROOT.h"
#include "TNamed.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"

const int TRUETOINTT = ((Int_t)(true)); // readability

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
  rangeStruct (float rangeLower_, float rangeUpper_) : rangeLower(rangeLower_),
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
  quadraticPolynomialStruct (float constCoefficient_, float linearCoefficient_, float squareCoefficient_) : constCoefficient(constCoefficient_),
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

enum class PFTypesForEA{chargedHadron=0, neutralHadron, photon};
struct EAValuesStruct{
  float regionUpperBound, chargedHadronsEA, neutralHadronsEA, photonsEA;
  
  EAValuesStruct (float regionUpperBound_, float chargedHadronsEA_, float neutralHadronsEA_, float photonsEA_) : regionUpperBound(regionUpperBound_),
    chargedHadronsEA(chargedHadronsEA_),
    neutralHadronsEA(neutralHadronsEA_),
    photonsEA(photonsEA_) {}

  float getEffectiveArea(const PFTypesForEA& PFType) const{
    float effectiveArea = 0.0;
    switch(PFType) {
    case (PFTypesForEA::chargedHadron) :
      effectiveArea = chargedHadronsEA;
      break;
    case (PFTypesForEA::neutralHadron) :
      effectiveArea = neutralHadronsEA;
      break;
    case (PFTypesForEA::photon) :
      effectiveArea = photonsEA;
      break;
    default :
      std::cout << "ERROR: Unknown PF type for EA!"<< std::endl;
      std::exit(EXIT_FAILURE);
    }
    return effectiveArea;
  }
  
  friend std::ostream& operator<< (std::ostream& out, const EAValuesStruct& EAValues) {
    out << "region upper eta bound --> " << EAValues.regionUpperBound << ", effective areas: "
        << "charged hadrons --> " << EAValues.chargedHadronsEA << ", "
        << "neutral hadrons --> " << EAValues.neutralHadronsEA << ", "
        << "photons --> " << EAValues.photonsEA;
    return out;
  }
};

enum class PhotonSelectionType{medium=0, fake, mediumfake};
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

struct parametersStruct {
  const float pTCutSubLeading = 25.0;
  const float pTCutLeading = 35.0;
  const float photonEtaCut = 1.442;
  const float R9Cut = 1.0;
  const float jetEtaCut = 2.4;
  const float jetpTCut = 30.;
  const float jetPUIDThreshold = 0.61;
  const float minDeltaRCut = 0.4;
  const float HTCut = 60.;
  const float electronPtCut = 15.;
  const float electronEtaCut = 2.5;
  const float electronDzCut = 0.1;
  const float electronPFPUIsolationCut = 0.1;
  const float muonPtCut = 15.0;
  const float muonPFPUIsolationCut = 0.12;
  const float METThreshold = 15.0;
  const PIDsStruct PIDs;

  // Default values for 2-medium-photon selection
  int nMediumPhotonsRequired = 2;
  int nFakePhotonsRequired = 0;
  void tuneParametersForPhotonSelectionType(const PhotonSelectionType& selectionType) {
    if (selectionType == PhotonSelectionType::fake) {
      nMediumPhotonsRequired = 0;
      nFakePhotonsRequired = 2;
    }
    else if (selectionType == PhotonSelectionType::mediumfake) {
      nMediumPhotonsRequired = 1;
      nFakePhotonsRequired = 1;
    }
  }

  // Default values for MC (year = -1): use 2017 parameters
  float towerHOverECut = 0.035;
  rangeStruct sigmaietaietaRange = rangeStruct(0.0103, 0.015);
  rangeStruct chargedIsolationRange = rangeStruct(1.416, 15.0);
  quadraticPolynomialStruct neutralIsolationCut = quadraticPolynomialStruct(2.491, 0.0126, 0.000026);
  quadraticPolynomialStruct photonIsolationCut = quadraticPolynomialStruct(2.952, 0.004, 0.0);
  EAValuesStruct region1EAs = EAValuesStruct(1.0, 0.0385, 0.0636, 0.124);
  EAValuesStruct region2EAs = EAValuesStruct(1.479, 0.0468, 0.1103, 0.1093);
  int HLTPhotonBit = -1;
  void tuneParametersForYear(const int& year) {
    if (year == 2017) { // only need to change HLTPhotonBit
      HLTPhotonBit = 36;
    }
    else if (year == 2016) {
      towerHOverECut = 0.0396;
      sigmaietaietaRange = rangeStruct(0.01022, 0.015);
      chargedIsolationRange = rangeStruct(0.441, 15.0);
      neutralIsolationCut = quadraticPolynomialStruct(2.725, 0.0148, 0.000017);
      photonIsolationCut = quadraticPolynomialStruct(2.571, 0.0047, 0.0);
      region1EAs = EAValuesStruct(1.0, 0.036, 0.0597, 0.121);
      region2EAs = EAValuesStruct(1.479, 0.0377, 0.0807, 0.1107);
      HLTPhotonBit = 14;
    }
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
        << "minDeltaR: " << parameters.minDeltaRCut << std::endl;

    out << "Electron cuts:" << std::endl
        << "pT: " << parameters.electronPtCut << ", "
        << "eta: " << parameters.electronEtaCut << ", "
        << "Dz: " << parameters.electronDzCut << ", "
        << "PFPUIso: " << parameters.electronPFPUIsolationCut << std::endl;

    out << "Muon cuts:" << std::endl
        << "pT: " << parameters.muonPtCut << ", "
        << "PFPUIso: " << parameters.muonPFPUIsolationCut << std::endl;

    out << "MET min: " << parameters.METThreshold << std::endl;

    out << "Event cuts:" << std::endl
        << "photon HLT bit index: " << parameters.HLTPhotonBit << ", "
        << "HT Cut: " << parameters.HTCut;
    return out;
  }
};

struct optionsStruct {
  std::string inputFilesList, outputFilePath;
  bool isMC;
  PhotonSelectionType photonSelectionType;
  Long64_t counterStartInclusive, counterEndInclusive;
  int year, JECUncertainty;

  friend std::ostream& operator<< (std::ostream& out, const optionsStruct& options) {
    out << "inputFilesList: " << options.inputFilesList << std::endl
        << "outputFilePath: " << options.outputFilePath << std::endl
        << "isMC: " << (options.isMC? "true": "false") << std::endl
        << "photonSelectionType: " << getPhotonSelectionTypeString(options.photonSelectionType) << std::endl
        << "Event range: [" << options.counterStartInclusive << ", " << options.counterEndInclusive << "]" << std::endl
        << "year: " << options.year << std::endl
        << "JEC Uncertainty scale: " << options.JECUncertainty;
    return out;
  }
};

enum class photonFailureCategory{eta=0, pT, hOverE, neutralIsolation, photonIsolation, conversionSafeElectronVeto, R9, sigmaietaiataORchargedIsolation, mediumIDCut, nPhotonFailureCategories};
int photonFailureCategoryFirst = static_cast<int>(photonFailureCategory::eta);
std::map<photonFailureCategory, std::string> photonFailureCategoryNames = {
  {photonFailureCategory::eta, "eta"},
  {photonFailureCategory::pT, "pT"},
  {photonFailureCategory::hOverE, "hOverE"},
  {photonFailureCategory::neutralIsolation, "neutralIsolation"},
  {photonFailureCategory::photonIsolation, "photonIsolation"},
  {photonFailureCategory::conversionSafeElectronVeto, "conversionSafeElectronVeto"},
  {photonFailureCategory::R9, "R9"},
  {photonFailureCategory::sigmaietaiataORchargedIsolation, "sigmaietaiataORchargedIsolation"},
  {photonFailureCategory::mediumIDCut, "mediumIDCut"}
};

enum class jetFailureCategory{eta=0, pT, PFLooseID, puID, jetID, deltaR, nJetFailureCategories};
int jetFailureCategoryFirst = static_cast<int>(jetFailureCategory::eta);
std::map<jetFailureCategory, std::string> jetFailureCategoryNames = {
  {jetFailureCategory::eta, "eta"},
  {jetFailureCategory::pT, "pT"},
  {jetFailureCategory::PFLooseID, "PFLooseID"},
  {jetFailureCategory::puID, "puID"},
  {jetFailureCategory::jetID, "jetID"},
  {jetFailureCategory::deltaR, "deltaR"}
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

enum class miscCounter{failingPhotons=0, passingPhotons, failingJets, passingJets, failingEvents, acceptedEvents, nMiscCounters};
int miscCounterFirst = static_cast<int>(miscCounter::failingPhotons);
std::map<miscCounter, std::string> miscCounterNames = {
  {miscCounter::failingPhotons, "failingPhotons"},
  {miscCounter::passingPhotons, "passingPhotons"},
  {miscCounter::failingJets, "failingJets"},
  {miscCounter::passingJets, "passingJets"},
  {miscCounter::failingEvents, "failingEvents"},
  {miscCounter::acceptedEvents, "acceptedEvents"}
};

enum class counterType{differential=0, global, nCounterTypes};
int counterTypeFirst = static_cast<int>(counterType::differential);
std::map<counterType, std::string> counterTypeNames = {
  {counterType::differential, "differential"},
  {counterType::global, "global"}
};
std::map<std::string, counterType> counterTypes = {
  {"differential", counterType::differential},
  {"global", counterType::global}
};

struct countersStruct{
  std::map<counterType, std::map<photonFailureCategory, Long64_t> > photonFailureCounters;
  std::map<counterType, std::map<jetFailureCategory, Long64_t> > jetFailureCounters;
  std::map<counterType, std::map<eventFailureCategory, Long64_t> > eventFailureCounters;
  std::map<miscCounter, Long64_t> miscCounters;
};

struct eventDetailsStruct{
  ULong64_t HLTPhotonBits;
  float eventRho;
  Int_t nPhotons;
  Int_t nJets;
  Int_t nElectrons;
  Int_t nMuons;
  float PFMET;
  Int_t nMCParticles;
};

struct MCCollectionStruct{
  std::vector<int> MCPIDs;
  std::vector<int> MCMomPIDs;
};

struct photonsCollectionStruct{
  std::vector<float> pT;
  std::vector<float> eta;
  std::vector<float> phi;
  std::vector<float> HOverE;
  std::vector<float> sigmaIEtaIEta;
  std::vector<float> PFChargedIsolationUncorrected;
  std::vector<float> PFNeutralIsolationUncorrected;
  std::vector<float> PFPhotonIsolationUncorrected;
  std::vector<UShort_t> ID;
  std::vector<int> electronVeto;
};

struct jetsCollectionStruct{
  std::vector<float> pT;
  std::vector<float> eta;
  std::vector<float> phi;
  std::vector<float> JECUncertainty;
  std::vector<float> PUID;
  std::vector<bool> looseID;
  std::vector<int> ID;
};

struct electronsCollectionStruct{
  std::vector<float> pT;
  std::vector<float> eta;
  std::vector<float> dz;
  std::vector<float> PFPUIsolation;
  std::vector<UShort_t> ID;
};

struct muonsCollectionStruct{
  std::vector<float> pT;
  std::vector<float> PFPUIsolation;
  std::vector<UShort_t> ID;
};

struct photonExaminationResultsStruct{
  bool passesSelectionAsMedium, passesSelectionAsFake, passesLeadingpTCut;
  float eta, phi, pT;

  photonExaminationResultsStruct (bool passesSelectionAsMedium_, bool passesSelectionAsFake_, bool passesLeadingpTCut_, float eta_, float phi_, float pT_) : passesSelectionAsMedium(passesSelectionAsMedium_), passesSelectionAsFake(passesSelectionAsFake_),passesLeadingpTCut(passesLeadingpTCut_), eta(eta_), phi(phi_), pT(pT_){}
};

struct jetExaminationResultsStruct{
  bool passesSelection;
  float eta, phi, pT;
  jetExaminationResultsStruct (bool passesSelection_, float eta_, float phi_, float pT_) : passesSelection(passesSelection_),
    eta(eta_),
    phi(phi_),
    pT(pT_) {}
};

struct angularVariablesStruct{
  float eta, phi;
  angularVariablesStruct (float eta_, float phi_) : eta(eta_), phi(phi_) {}
};
