#ifndef H_COLLECTIONS
#define H_COLLECTIONS

#include <vector>

#include "TROOT.h"
#include "TChain.h"

struct eventDetailsStruct{
  ULong64_t HLTPhotonBits;
  float eventRho;
  Int_t nPhotons;
  Int_t nJets;
  Int_t nElectrons;
  Int_t nMuons;
  float PFMET;
  float PFMET_UnclusteredDown;
  float PFMET_UnclusteredUp;
  float PFMET_JERDown;
  float PFMET_JERUp;
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
      inputChain.SetBranchAddress("pfMET_T1UESDo", &(PFMET_UnclusteredDown));
      inputChain.SetBranchStatus("pfMET_T1UESDo", 1);
      inputChain.SetBranchAddress("pfMET_T1UESUp", &(PFMET_UnclusteredUp));
      inputChain.SetBranchStatus("pfMET_T1UESUp", 1);
      inputChain.SetBranchAddress("pfMET_T1JERDo", &(PFMET_JERDown));
      inputChain.SetBranchStatus("pfMET_T1JERDo", 1);
      inputChain.SetBranchAddress("pfMET_T1JERUp", &(PFMET_JERUp));
      inputChain.SetBranchStatus("pfMET_T1JERUp", 1);
      inputChain.SetBranchAddress("nMC", &(nMCParticles));
      inputChain.SetBranchStatus("nMC", 1);
    }
  }
};

struct MCCollectionStruct{
  std::vector<int> * MCPIDs = nullptr;
  std::vector<int> * MCMomPIDs = nullptr;
  std::vector<UShort_t> * MCStatusFlags = nullptr;
  std::vector<float> * MCMasses = nullptr;
  std::vector<float> * MCMomMasses = nullptr;
  std::vector<float> * MCEts = nullptr;
  std::vector<float> * MCEtas = nullptr;
  std::vector<float> * MCPhis = nullptr;

  MCCollectionStruct(TChain &inputChain, const bool& isMC) {
    if (isMC) {
      inputChain.SetBranchAddress("mcPID", &(MCPIDs));
      inputChain.SetBranchStatus("mcPID", 1);
      inputChain.SetBranchAddress("mcMomPID", &(MCMomPIDs));
      inputChain.SetBranchStatus("mcMomPID", 1);
      inputChain.SetBranchAddress("mcStatusFlag", &(MCStatusFlags));
      inputChain.SetBranchStatus("mcStatusFlag", 1);
      inputChain.SetBranchAddress("mcMass", &(MCMasses));
      inputChain.SetBranchStatus("mcMass", 1);
      inputChain.SetBranchAddress("mcMomMass", &(MCMomMasses));
      inputChain.SetBranchStatus("mcMomMass", 1);
      inputChain.SetBranchAddress("mcEt", &(MCEts));
      inputChain.SetBranchStatus("mcEt", 1);
      inputChain.SetBranchAddress("mcEta", &(MCEtas));
      inputChain.SetBranchStatus("mcEta", 1);
      inputChain.SetBranchAddress("mcPhi", &(MCPhis));
      inputChain.SetBranchStatus("mcPhi", 1);
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

#endif
