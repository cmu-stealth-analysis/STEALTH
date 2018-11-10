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
#include "TMath.h"
#include "TNamed.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TChain.h"
#include "TLorentzVector.h"
#include "TH2F.h"
#include "TEfficiency.h"

namespace constants{ // for readability
  const int TRUETOINTT = ((Int_t)(true));
  const float VALUEOFTWOPI = static_cast<float>(2.0*TMath::Pi());
}

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

struct quadraticPolynomialStruct{
  float constCoefficient, linearCoefficient, squareCoefficient;
  quadraticPolynomialStruct () : constCoefficient(0.),
    linearCoefficient(0.),
    squareCoefficient(0.) {}

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

  EAValuesStruct () : regionUpperBound(0.),
    chargedHadronsEA(0.),
    neutralHadronsEA(0.),
    photonsEA(0.) {}
  
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
  const float jetEtaCut = 2.5f;
  const float jetpTCut = 30.f;
  const float jetPUIDThreshold = 0.61f;
  const float minDeltaRCut = 0.4f;
  const float HTCut = 60.0f;
  /* const float electronPtCut = 15.0f; */
  /* const float electronEtaCut = 2.5f; */
  /* const float electronDzCut = 0.1f; */
  /* const float electronPFPUIsolationCut = 0.1f; */
  /* const float muonPtCut = 15.0f; */
  /* const float muonPFPUIsolationCut = 0.12f; */
  const PIDsStruct PIDs;
  const UShort_t MCStatusFlagBitMask = static_cast<UShort_t>(7u);

  int nMediumPhotonsRequired, nFakePhotonsRequired;
  void tuneParametersForPhotonSelectionType(const PhotonSelectionType& selectionType) {
    if (selectionType == PhotonSelectionType::medium) {
      nMediumPhotonsRequired = 2;
      nFakePhotonsRequired = 0;
    }
    else if (selectionType == PhotonSelectionType::fake) {
      nMediumPhotonsRequired = 0;
      nFakePhotonsRequired = 2;
    }
    else if (selectionType == PhotonSelectionType::mediumfake) {
      nMediumPhotonsRequired = 1;
      nFakePhotonsRequired = 1;
    }
  }

  int HLTPhotonBit;
  float invariantMassCut, towerHOverECut, sigmaIEtaIEtaCut, sigmaIEtaIEtaCutLoose, chargedIsolationCut, chargedIsolationCutLoose;
  quadraticPolynomialStruct neutralIsolationCut, photonIsolationCut;
  EAValuesStruct region1EAs, region2EAs;
  TFile* sourceFile_efficiencyMap_2016;
  TEfficiency* efficiencyMap_2016;
  TFile* sourceFile_efficiencyMap_2017;
  TH2F* efficiencyMap_2017;
  void tuneParametersForYear(const int& year, const bool& isMC) {
    if (year == 2017) {
      if (isMC) HLTPhotonBit = -1;
      else HLTPhotonBit = 37;
      invariantMassCut = 60.0f;

      towerHOverECut = 0.02197f;
      sigmaIEtaIEtaCut = 0.01015f;
      sigmaIEtaIEtaCutLoose = 0.02f;
      chargedIsolationCut = 1.141f;
      chargedIsolationCutLoose = 6.0f;
      neutralIsolationCut = quadraticPolynomialStruct(1.189f, 0.01512f, 0.00002259f);
      photonIsolationCut = quadraticPolynomialStruct(2.08f, 0.004017f, 0.0f);
      region1EAs = EAValuesStruct(1.0f, 0.0112f, 0.0668f, 0.1113f);
      region2EAs = EAValuesStruct(1.479f, 0.0108f, 0.1054f, 0.0953f);

      sourceFile_efficiencyMap_2017 = TFile::Open("eventSelection/data/L1prefiring_jet_2017BtoF.root", "READ");
      if (!(sourceFile_efficiencyMap_2017->IsOpen()) || sourceFile_efficiencyMap_2017->IsZombie()) {
        std::cout << "ERROR: Unable to open file with path: eventSelection/data/L1prefiring_jet_2017BtoF.root" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      sourceFile_efficiencyMap_2017->GetObject("L1prefiring_jet_2017BtoF", efficiencyMap_2017);
      if (efficiencyMap_2017) std::cout << "Opened prefiring efficiency map for 2017" << std::endl;
      else {
        std::cout << "ERROR: Unable to open histogram with path: L1prefiring_jet2017BtoF" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
    else if (year == 2016) {
      if (isMC) HLTPhotonBit = -1;
      else HLTPhotonBit = 16;
      invariantMassCut = 60.0f;

      towerHOverECut = 0.0396f;
      sigmaIEtaIEtaCut = 0.01022f;
      sigmaIEtaIEtaCutLoose = 0.02f;
      chargedIsolationCut = 0.441f;
      chargedIsolationCutLoose = 6.0f;
      neutralIsolationCut = quadraticPolynomialStruct(2.725f, 0.0148f, 0.000017f);
      photonIsolationCut = quadraticPolynomialStruct(2.571f, 0.0047f, 0.0f);
      region1EAs = EAValuesStruct(1.0f, 0.036f, 0.0597f, 0.121f);
      region2EAs = EAValuesStruct(1.479f, 0.0377f, 0.0807f, 0.1107f);

      sourceFile_efficiencyMap_2016 = TFile::Open("eventSelection/data/Map_Jet_L1IsoEG30eff_bxm1_looseJet_SingleMuon_Run2016B-H.root", "READ");
      if (!(sourceFile_efficiencyMap_2016->IsOpen()) || sourceFile_efficiencyMap_2016->IsZombie()) {
        std::cout << "ERROR: Unable to open file with path: eventSelection/data/Map_Jet_L1IsoEG30eff_bxm1_looseJet_SingleMuon_Run2016B-H.root" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      sourceFile_efficiencyMap_2016->GetObject("prefireEfficiencyMap", efficiencyMap_2016);
      if (efficiencyMap_2016) std::cout << "Opened prefiring efficiency map for 2016" << std::endl;
      else {
        std::cout << "ERROR: Unable to open histogram with path: prefireEfficiencyMap" << std::endl;
        std::exit(EXIT_FAILURE);
      }
    }
  }
  friend std::ostream& operator<< (std::ostream& out, const parametersStruct& parameters) {
    out << "PIDs: " << parameters.PIDs << std::endl;

    out << "Photon cuts:" << std::endl
        << "pT_SubLeading: " << parameters.pTCutSubLeading << ", "
        << "pT_Leading: " << parameters.pTCutLeading << ", "
        << "eta: " << parameters.photonEtaCut << ", "
        << "HOverE: " << parameters.towerHOverECut << ", "
        << "sigmaietaieta cut: " << parameters.sigmaIEtaIEtaCut << ", "
        << "sigmaietaieta cut (loose): " << parameters.sigmaIEtaIEtaCutLoose << ", "
        << "charged isolation cut: " << parameters.chargedIsolationCut << ", "
        << "charged isolation cut (loose): " << parameters.chargedIsolationCutLoose << ", "
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

    /* out << "Electron cuts:" << std::endl */
    /*     << "pT: " << parameters.electronPtCut << ", " */
    /*     << "eta: " << parameters.electronEtaCut << ", " */
    /*     << "Dz: " << parameters.electronDzCut << ", " */
    /*     << "PFPUIso: " << parameters.electronPFPUIsolationCut << std::endl; */

    /* out << "Muon cuts:" << std::endl */
    /*     << "pT: " << parameters.muonPtCut << ", " */
    /*     << "PFPUIso: " << parameters.muonPFPUIsolationCut << std::endl; */

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

enum class photonFailureCategory{eta=0, pT, hOverE, neutralIsolation, photonIsolation, conversionSafeElectronVeto, sigmaietaiataANDchargedIsolation, sigmaietaiataANDchargedIsolationLoose, nPhotonFailureCategories};
int photonFailureCategoryFirst = static_cast<int>(photonFailureCategory::eta);
std::map<photonFailureCategory, std::string> photonFailureCategoryNames = {
  {photonFailureCategory::eta, "eta"},
  {photonFailureCategory::pT, "pT"},
  {photonFailureCategory::hOverE, "hOverE"},
  {photonFailureCategory::neutralIsolation, "neutralIsolation"},
  {photonFailureCategory::photonIsolation, "photonIsolation"},
  {photonFailureCategory::conversionSafeElectronVeto, "conversionSafeElectronVeto"},
  {photonFailureCategory::sigmaietaiataANDchargedIsolation, "sigmaietaiataANDchargedIsolation"},
  {photonFailureCategory::sigmaietaiataANDchargedIsolationLoose, "sigmaietaiataANDchargedIsolationLoose"}
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

enum class eventFailureCategory{HLTPhoton=0, wrongNSelectedPhotons, incompatiblePhotonSelectionType, lowInvariantMass, HLTJet, wrongNJets, hTCut, electronVeto, muonVeto, MCGenInformation, nEventFailureCategories};
int eventFailureCategoryFirst = static_cast<int>(eventFailureCategory::HLTPhoton);
std::map<eventFailureCategory, std::string> eventFailureCategoryNames = {
  {eventFailureCategory::HLTPhoton, "HLTPhoton"},
  {eventFailureCategory::wrongNSelectedPhotons, "wrongNSelectedPhotons"},
  {eventFailureCategory::incompatiblePhotonSelectionType, "incompatiblePhotonSelectionType"},
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

/* struct electronsCollectionStruct{ */
/*   std::vector<float> * pT = nullptr; */
/*   std::vector<float> * eta = nullptr; */
/*   std::vector<float> * dz = nullptr; */
/*   std::vector<float> * PFPUIsolation = nullptr; */
/*   std::vector<UShort_t> * ID = nullptr; */

/*   electronsCollectionStruct(TChain &inputChain) { */
/*     inputChain.SetBranchAddress("elePt", &(pT)); */
/*     inputChain.SetBranchStatus("elePt", 1); */
/*     inputChain.SetBranchAddress("eleEta", &(eta)); */
/*     inputChain.SetBranchStatus("eleEta", 1); */
/*     inputChain.SetBranchAddress("eleDz", &(dz)); */
/*     inputChain.SetBranchStatus("eleDz", 1); */
/*     inputChain.SetBranchAddress("elePFPUIso", &(PFPUIsolation)); */
/*     inputChain.SetBranchStatus("elePFPUIso", 1); */
/*     inputChain.SetBranchAddress("eleIDbit", &(ID)); */
/*     inputChain.SetBranchStatus("eleIDbit", 1); */
/*   } */
/* }; */

/* struct muonsCollectionStruct{ */
/*   std::vector<float> * pT = nullptr; */
/*   std::vector<float> * PFPUIsolation = nullptr; */
/*   std::vector<UShort_t> * ID = nullptr; */

/*   muonsCollectionStruct(TChain &inputChain) { */
/*     inputChain.SetBranchAddress("muPt", &(pT)); */
/*     inputChain.SetBranchStatus("muPt", 1); */
/*     inputChain.SetBranchAddress("muPFPUIso", &(PFPUIsolation)); */
/*     inputChain.SetBranchStatus("muPFPUIso", 1); */
/*     inputChain.SetBranchAddress("muIDbit", &(ID)); */
/*     inputChain.SetBranchStatus("muIDbit", 1); */
/*   } */
/* }; */

struct photonExaminationResultsStruct{
  bool passesSelectionAsMedium, passesSelectionAsFake, passesLeadingpTCut;
  float eta, phi, pT, energy;

  photonExaminationResultsStruct (bool passesSelectionAsMedium_, bool passesSelectionAsFake_, bool passesLeadingpTCut_, float eta_, float phi_, float pT_, float energy_) : passesSelectionAsMedium(passesSelectionAsMedium_), passesSelectionAsFake(passesSelectionAsFake_),passesLeadingpTCut(passesLeadingpTCut_), eta(eta_), phi(phi_), pT(pT_), energy(energy_){}
};

struct jetExaminationResultsStruct{
  bool passesSelection;
  float eta, phi, pT, prefireWeight;
  jetExaminationResultsStruct (bool passesSelection_, float eta_, float phi_, float pT_, float prefireWeight_) : passesSelection(passesSelection_),
    eta(eta_),
    phi(phi_),
    pT(pT_),
    prefireWeight(prefireWeight_) {}
};

struct angularVariablesStruct{
  float eta, phi;
  angularVariablesStruct (float eta_, float phi_) : eta(eta_), phi(phi_) {}
};

struct extraEventInfoStruct{
  Long64_t eventIndex;
  int evt_nJetsDR;
  float evt_ST;
  float scaleFactor;

  extraEventInfoStruct (Long64_t eventIndex_, int evt_nJetsDR_, float evt_ST_, float scaleFactor_) : eventIndex(eventIndex_),
    evt_nJetsDR(evt_nJetsDR_),
    evt_ST(evt_ST_),
    scaleFactor(scaleFactor_) {}
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
  if (!(options.year == 2016 || options.year == 2017)) {
    std::cout << "ERROR: argument \"year\" can be one of 2016 or 2017; current value: " << options.year << std::endl;
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

bool passesBitMask(const UShort_t& bitCollection, const UShort_t& bitMask) {
  return ((bitCollection&bitMask) == bitMask);
}
