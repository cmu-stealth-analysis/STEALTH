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
#include "TLorentzVector.h"

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
  const float pTCutSubLeading = 25.0f;
  const float pTCutLeading = 35.0f;
  const float photonEtaCut = 1.442f;
  const float invariantMassCut = 60.0f;
  const float jetEtaCut = 2.5f;
  const float jetpTCut = 30.f;
  const float jetPUIDThreshold = 0.61f;
  const float minDeltaRCut = 0.4f;
  const float HTCut = 60.0f;
  const float electronPtCut = 15.f;
  const float electronEtaCut = 2.5f;
  const float electronDzCut = 0.1f;
  const float electronPFPUIsolationCut = 0.1f;
  const float muonPtCut = 15.0f;
  const float muonPFPUIsolationCut = 0.12f;
  const PIDsStruct PIDs;
  const UShort_t MCStatusFlagConstraint = static_cast<UShort_t>(7);

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
  float towerHOverECut = 0.035f;
  rangeStruct sigmaietaietaRange = rangeStruct(0.0103f, 0.02f);
  rangeStruct chargedIsolationRange = rangeStruct(1.416f, 6.0f);
  quadraticPolynomialStruct neutralIsolationCut = quadraticPolynomialStruct(2.491f, 0.0126f, 0.000026f);
  quadraticPolynomialStruct photonIsolationCut = quadraticPolynomialStruct(2.952f, 0.004f, 0.0f);
  EAValuesStruct region1EAs = EAValuesStruct(1.0f, 0.0385f, 0.0636f, 0.124f);
  EAValuesStruct region2EAs = EAValuesStruct(1.479f, 0.0468f, 0.1103f, 0.1093f);
  int HLTPhotonBit = -1;
  void tuneParametersForYear(const int& year) {
    if (year == 2017) { // only need to change HLTPhotonBit
      HLTPhotonBit = 36;
    }
    else if (year == 2016) {
      // Disabling everything except the HLT photon bit at the moment
      /* towerHOverECut = 0.0396; */
      /* sigmaietaietaRange = rangeStruct(0.01022, 0.015); */
      /* chargedIsolationRange = rangeStruct(0.441, 15.0); */
      /* neutralIsolationCut = quadraticPolynomialStruct(2.725, 0.0148, 0.000017); */
      /* photonIsolationCut = quadraticPolynomialStruct(2.571, 0.0047, 0.0); */
      /* region1EAs = EAValuesStruct(1.0, 0.036, 0.0597, 0.121); */
      /* region2EAs = EAValuesStruct(1.479, 0.0377, 0.0807, 0.1107); */
      HLTPhotonBit = 14;
    }
  }
  friend std::ostream& operator<< (std::ostream& out, const parametersStruct& parameters) {
    out << "PIDs: " << parameters.PIDs << std::endl;

    out << "Photon cuts:" << std::endl
        << "pT_SubLeading: " << parameters.pTCutSubLeading << ", "
        << "pT_Leading: " << parameters.pTCutLeading << ", "
        << "eta: " << parameters.photonEtaCut << ", "
        << "HOverE: " << parameters.towerHOverECut << ", "
        << "sigmaietaietaRange: " << parameters.sigmaietaietaRange << ", "
        << "chargedIsolationRange: " << parameters.chargedIsolationRange << ", "
        << "neutral isolation cut coefficients: " << parameters.neutralIsolationCut << ", "
        << "photon isolation cut coefficients: " << parameters.photonIsolationCut << std::endl;

    out << "Region 1 effective areas: " << parameters.region1EAs << std::endl
        << "Region 2 effective areas: " << parameters.region2EAs << std::endl;

    out << "Invariant mass cut: " << parameters.invariantMassCut << std::endl;

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

enum class photonFailureCategory{eta=0, pT, hOverE, neutralIsolation, photonIsolation, conversionSafeElectronVeto, sigmaietaiataORchargedIsolation, mediumIDCut, nPhotonFailureCategories};
int photonFailureCategoryFirst = static_cast<int>(photonFailureCategory::eta);
std::map<photonFailureCategory, std::string> photonFailureCategoryNames = {
  {photonFailureCategory::eta, "eta"},
  {photonFailureCategory::pT, "pT"},
  {photonFailureCategory::hOverE, "hOverE"},
  {photonFailureCategory::neutralIsolation, "neutralIsolation"},
  {photonFailureCategory::photonIsolation, "photonIsolation"},
  {photonFailureCategory::conversionSafeElectronVeto, "conversionSafeElectronVeto"},
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

enum class eventFailureCategory{HLTPhoton=0, wrongNMediumOrFakePhotons, wrongNPhotons, lowInvariantMass, HLTJet, wrongNJets, hTCut, electronVeto, muonVeto, MCGenInformation, nEventFailureCategories};
int eventFailureCategoryFirst = static_cast<int>(eventFailureCategory::HLTPhoton);
std::map<eventFailureCategory, std::string> eventFailureCategoryNames = {
  {eventFailureCategory::HLTPhoton, "HLTPhoton"},
  {eventFailureCategory::wrongNMediumOrFakePhotons, "wrongNMediumOrFakePhotons"},
  {eventFailureCategory::wrongNPhotons, "wrongNPhotons"},
  {eventFailureCategory::lowInvariantMass, "lowInvariantMass"},
  {eventFailureCategory::HLTJet, "HLTJet"},
  {eventFailureCategory::wrongNJets, "wrongNJets"},
  {eventFailureCategory::hTCut, "hTCut"},
  {eventFailureCategory::electronVeto, "electronVeto"},
  {eventFailureCategory::muonVeto, "muonVeto"},
  {eventFailureCategory::MCGenInformation, "MCGenInformation"}
};

enum class miscCounter{failingPhotons=0, passingPhotons, totalPhotons, failingJets, passingJets, totalJets, failingEvents, acceptedEvents, totalEvents, nMiscCounters};
int miscCounterFirst = static_cast<int>(miscCounter::failingPhotons);
std::map<miscCounter, std::string> miscCounterNames = {
  {miscCounter::failingPhotons, "failingPhotons"},
  {miscCounter::passingPhotons, "passingPhotons"},
  {miscCounter::totalPhotons, "totalPhotons"},
  {miscCounter::failingJets, "failingJets"},
  {miscCounter::passingJets, "passingJets"},
  {miscCounter::totalJets, "totalJets"},
  {miscCounter::failingEvents, "failingEvents"},
  {miscCounter::acceptedEvents, "acceptedEvents"},
  {miscCounter::totalEvents, "totalEvents"}
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

  eventDetailsStruct(TChain &inputChain, const bool& isMC) {
    inputChain.SetBranchAddress("HLTPho", &(HLTPhotonBits));
    inputChain.SetBranchStatus("HLTPho", 1);
    inputChain.SetBranchAddress("rho", &(eventRho));
    inputChain.SetBranchStatus("rho", 1);
    inputChain.SetBranchAddress("nPho", &(nPhotons));
    inputChain.SetBranchStatus("nPho", 1);
    inputChain.SetBranchAddress("nJet", &(nJets));
    inputChain.SetBranchStatus("nJet", 1);
    inputChain.SetBranchAddress("nEle", &(nElectrons));
    inputChain.SetBranchStatus("nEle", 1);
    inputChain.SetBranchAddress("nMu", &(nMuons));
    inputChain.SetBranchStatus("nMu", 1);
    inputChain.SetBranchAddress("pfMET", &(PFMET));
    inputChain.SetBranchStatus("pfMET", 1);
    if (isMC) {
      inputChain.SetBranchAddress("nMC", &(nMCParticles));
      inputChain.SetBranchStatus("nMC", 1);
    }
  }
};

struct MCCollectionStruct{
  std::vector<int> * MCPIDs = nullptr;
  std::vector<int> * MCMomPIDs = nullptr;
  std::vector<UShort_t> * MCStatusFlags = nullptr;

  MCCollectionStruct(TChain &inputChain, const bool& isMC) {
    if (isMC) {
      inputChain.SetBranchAddress("mcPID", &(MCPIDs));
      inputChain.SetBranchStatus("mcPID", 1);
      inputChain.SetBranchAddress("mcMomPID", &(MCMomPIDs));
      inputChain.SetBranchStatus("mcMomPID", 1);
      inputChain.SetBranchAddress("mcStatusFlag", &(MCStatusFlags));
      inputChain.SetBranchStatus("mcStatusFlag", 1);
    }
  }
};

struct photonsCollectionStruct{
  std::vector<float> * pT = nullptr;
  std::vector<float> * eta = nullptr;
  std::vector<float> * phi = nullptr;
  std::vector<float> * HOverE = nullptr;
  std::vector<float> * sigmaIEtaIEta = nullptr;
  std::vector<float> * PFChargedIsolationUncorrected = nullptr;
  std::vector<float> * PFNeutralIsolationUncorrected = nullptr;
  std::vector<float> * PFPhotonIsolationUncorrected = nullptr;
  std::vector<UShort_t> * ID = nullptr;
  std::vector<int> * electronVeto = nullptr;
  std::vector<float> * energy = nullptr;

  photonsCollectionStruct(TChain &inputChain) {
    inputChain.SetBranchAddress("phoEt", &(pT));
    inputChain.SetBranchStatus("phoEt", 1);
    inputChain.SetBranchAddress("phoEta", &(eta));
    inputChain.SetBranchStatus("phoEta", 1);
    inputChain.SetBranchAddress("phoPhi", &(phi));
    inputChain.SetBranchStatus("phoPhi", 1);
    inputChain.SetBranchAddress("phoHoverE", &(HOverE));
    inputChain.SetBranchStatus("phoHoverE", 1);
    inputChain.SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &(sigmaIEtaIEta));
    inputChain.SetBranchStatus("phoSigmaIEtaIEtaFull5x5", 1);
    inputChain.SetBranchAddress("phoPFChIso", &(PFChargedIsolationUncorrected));
    inputChain.SetBranchStatus("phoPFChIso", 1);
    inputChain.SetBranchAddress("phoPFNeuIso", &(PFNeutralIsolationUncorrected));
    inputChain.SetBranchStatus("phoPFNeuIso", 1);
    inputChain.SetBranchAddress("phoPFPhoIso", &(PFPhotonIsolationUncorrected));
    inputChain.SetBranchStatus("phoPFPhoIso", 1);
    inputChain.SetBranchAddress("phoIDbit", &(ID));
    inputChain.SetBranchStatus("phoIDbit", 1);
    inputChain.SetBranchAddress("phoEleVeto", &(electronVeto));
    inputChain.SetBranchStatus("phoEleVeto", 1);
    inputChain.SetBranchAddress("phoE", &(energy));
    inputChain.SetBranchStatus("phoE", 1);
  }
};

struct jetsCollectionStruct{
  std::vector<float> * pT = nullptr;
  std::vector<float> * eta = nullptr;
  std::vector<float> * phi = nullptr;
  std::vector<float> * JECUncertainty = nullptr;
  std::vector<float> * PUID = nullptr;
  std::vector<bool> * looseID = nullptr;
  std::vector<int> * ID = nullptr;

  jetsCollectionStruct(TChain &inputChain, const bool& isMC) {
    inputChain.SetBranchAddress("jetPt", &(pT));
    inputChain.SetBranchStatus("jetPt", 1);
    inputChain.SetBranchAddress("jetEta", &(eta));
    inputChain.SetBranchStatus("jetEta", 1);
    inputChain.SetBranchAddress("jetPhi", &(phi));
    inputChain.SetBranchStatus("jetPhi", 1);
    inputChain.SetBranchAddress("jetPUID", &(PUID));
    inputChain.SetBranchStatus("jetPUID", 1);
    inputChain.SetBranchAddress("jetPFLooseId", &(looseID));
    inputChain.SetBranchStatus("jetPFLooseId", 1);
    inputChain.SetBranchAddress("jetID", &(ID));
    inputChain.SetBranchStatus("jetID", 1);
    if (isMC) {
      inputChain.SetBranchAddress("jetJECUnc", &(JECUncertainty));
      inputChain.SetBranchStatus("jetJECUnc", 1);
    }
  }
};

struct electronsCollectionStruct{
  std::vector<float> * pT = nullptr;
  std::vector<float> * eta = nullptr;
  std::vector<float> * dz = nullptr;
  std::vector<float> * PFPUIsolation = nullptr;
  std::vector<UShort_t> * ID = nullptr;

  electronsCollectionStruct(TChain &inputChain) {
    inputChain.SetBranchAddress("elePt", &(pT));
    inputChain.SetBranchStatus("elePt", 1);
    inputChain.SetBranchAddress("eleEta", &(eta));
    inputChain.SetBranchStatus("eleEta", 1);
    inputChain.SetBranchAddress("eleDz", &(dz));
    inputChain.SetBranchStatus("eleDz", 1);
    inputChain.SetBranchAddress("elePFPUIso", &(PFPUIsolation));
    inputChain.SetBranchStatus("elePFPUIso", 1);
    inputChain.SetBranchAddress("eleIDbit", &(ID));
    inputChain.SetBranchStatus("eleIDbit", 1);
  }
};

struct muonsCollectionStruct{
  std::vector<float> * pT = nullptr;
  std::vector<float> * PFPUIsolation = nullptr;
  std::vector<UShort_t> * ID = nullptr;

  muonsCollectionStruct(TChain &inputChain) {
    inputChain.SetBranchAddress("muPt", &(pT));
    inputChain.SetBranchStatus("muPt", 1);
    inputChain.SetBranchAddress("muPFPUIso", &(PFPUIsolation));
    inputChain.SetBranchStatus("muPFPUIso", 1);
    inputChain.SetBranchAddress("muIDbit", &(ID));
    inputChain.SetBranchStatus("muIDbit", 1);
  }
};

struct photonExaminationResultsStruct{
  bool passesSelectionAsMedium, passesSelectionAsFake, passesLeadingpTCut;
  float eta, phi, pT, energy;

  photonExaminationResultsStruct (bool passesSelectionAsMedium_, bool passesSelectionAsFake_, bool passesLeadingpTCut_, float eta_, float phi_, float pT_, float energy_) : passesSelectionAsMedium(passesSelectionAsMedium_), passesSelectionAsFake(passesSelectionAsFake_),passesLeadingpTCut(passesLeadingpTCut_), eta(eta_), phi(phi_), pT(pT_), energy(energy_){}
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

struct extraEventInfoStruct{
  Long64_t eventIndex;
  int evt_nJetsDR;
  float evt_ST;

  extraEventInfoStruct (Long64_t eventIndex_, int evt_nJetsDR_, float evt_ST_) : eventIndex(eventIndex_),
    evt_nJetsDR(evt_nJetsDR_),
    evt_ST(evt_ST_) {}
};

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
  if (options.isMC && !(options.year == -1)) {
    std::cout << "ERROR: Expected argument year=-1 with isMC=true; current values: isMC: " << (options.isMC ? "true" : "false") << ", year: " << options.year << std::endl;
  }
  return options;
}

std::string getNDashes(const int& n) {
  std::stringstream dashes;
  for (int counter = 0; counter < n; ++counter) dashes << "-";
  return dashes.str();
}

void initializeCounters(countersStruct &counters) {
  for (int counterIndex = counterTypeFirst; counterIndex != static_cast<int>(counterType::nCounterTypes); ++counterIndex) {
    counterType typeIndex = static_cast<counterType>(counterIndex);
    for (int categoryIndex = photonFailureCategoryFirst; categoryIndex != static_cast<int>(photonFailureCategory::nPhotonFailureCategories); ++categoryIndex) {
      photonFailureCategory category = static_cast<photonFailureCategory>(categoryIndex);
      counters.photonFailureCounters[typeIndex][category] = 0l;
    }

    for (int categoryIndex = jetFailureCategoryFirst; categoryIndex != static_cast<int>(jetFailureCategory::nJetFailureCategories); ++categoryIndex) {
      jetFailureCategory category = static_cast<jetFailureCategory>(categoryIndex);
      counters.jetFailureCounters[typeIndex][category] = 0l;
    }

    for (int categoryIndex = eventFailureCategoryFirst; categoryIndex != static_cast<int>(eventFailureCategory::nEventFailureCategories); ++categoryIndex) {
      eventFailureCategory category = static_cast<eventFailureCategory>(categoryIndex);
      counters.eventFailureCounters[typeIndex][category] = 0l;
    }
  }

  for (int miscCounterIndex = miscCounterFirst; miscCounterIndex != static_cast<int>(miscCounter::nMiscCounters); ++miscCounterIndex) {
    miscCounter miscCounterEnumIndex = static_cast<miscCounter>(miscCounterIndex);
    counters.miscCounters[miscCounterEnumIndex] = 0l;
  }
}

void printCounters(countersStruct &counters) {
  for (const auto& counterTypeMapElement : counterTypes) {
    std::string counterTypeString = counterTypeMapElement.first;
    std::cout << counterTypeString << " photon failure counters: " << std::endl;
    for (const auto& counterValuePair : counters.photonFailureCounters[counterTypeMapElement.second]) {
      std::cout << photonFailureCategoryNames[counterValuePair.first] << " : " << counterValuePair.second << " = " << std::setprecision(3) << 100.0*(static_cast<double>(counterValuePair.second)/((counters.miscCounters)[miscCounter::failingPhotons])) << " %" << std::endl;
    }
    std::cout << getNDashes(100) << std::endl;

    std::cout << counterTypeString << " jet failure counters: " << std::endl;
    for (const auto& counterValuePair : counters.jetFailureCounters[counterTypeMapElement.second]) {
      std::cout << jetFailureCategoryNames[counterValuePair.first] << " : " << counterValuePair.second << " = " << std::setprecision(3) << 100.0*(static_cast<double>(counterValuePair.second)/((counters.miscCounters)[miscCounter::failingJets])) << " %" << std::endl;
    }
    std::cout << getNDashes(100) << std::endl;

    std::cout << counterTypeString << " event failure counters: " << std::endl;
    for (const auto& counterValuePair : counters.eventFailureCounters[counterTypeMapElement.second]) {
      std::cout << eventFailureCategoryNames[counterValuePair.first] << " : " << counterValuePair.second << " = " << std::setprecision(3) << 100.0*(static_cast<double>(counterValuePair.second)/((counters.miscCounters)[miscCounter::failingEvents])) << " %" << std::endl;
    }
    std::cout << getNDashes(100) << std::endl;

  }
  
  std::cout << "Miscellaneous counters: " << std::endl;
  for (const auto& counterValuePair : counters.miscCounters) {
    std::cout << miscCounterNames[counterValuePair.first] << " : " << counterValuePair.second << std::endl;
  }

  std::cout << "Accepted events: " << ((counters.miscCounters)[miscCounter::acceptedEvents]) << "/" << ((counters.miscCounters)[miscCounter::totalEvents]) << " = " << std::setprecision(4) << (100.0*(static_cast<double>((counters.miscCounters)[miscCounter::acceptedEvents]))/((counters.miscCounters)[miscCounter::totalEvents]))<< " %" << std::endl;
}

void incrementCounters(const photonFailureCategory& photonCategory, const counterType& counterTypeIndex, countersStruct& counters) {
  ++(((counters.photonFailureCounters)[counterTypeIndex])[photonCategory]);
}

void incrementCounters(const jetFailureCategory& jetCategory, const counterType& counterTypeIndex, countersStruct& counters) {
  ++(((counters.jetFailureCounters)[counterTypeIndex])[jetCategory]);
}

void incrementCounters(const eventFailureCategory& eventCategory, const counterType& counterTypeIndex, countersStruct& counters) {
  ++(((counters.eventFailureCounters)[counterTypeIndex])[eventCategory]);
}

void incrementCounters(const miscCounter& miscCounterEnumIndex, countersStruct& counters) {
  ++((counters.miscCounters)[miscCounterEnumIndex]);
}

template<typename failureCategory>
void applyCondition(countersStruct &counters, const failureCategory& category, bool& passesAllCriteriaSoFar, const bool& testCondition) {
  if (!testCondition) {
    incrementCounters(category, counterType::global, counters);
    if (passesAllCriteriaSoFar) {
      incrementCounters(category, counterType::differential, counters);
      passesAllCriteriaSoFar = false;
    }
  }
}
