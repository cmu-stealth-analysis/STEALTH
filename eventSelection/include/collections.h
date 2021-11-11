#ifndef H_COLLECTIONS
#define H_COLLECTIONS

#include <vector>

#include "TROOT.h"
#include "TChain.h"

struct eventDetailsStruct{
  ULong64_t HLTPhotonBits;
  ULong64_t HLTJetBits;
  float eventRho;
  Int_t nPhotons;
  Int_t nJets;
  Int_t nElectrons;
  Int_t nMuons;
  float PFMET;
  float PFMET_phi;
  float PFMET_UnclusteredDown;
  float PFMET_UnclusteredUp;
  float PFMET_JERDown;
  float PFMET_JERUp;
  Int_t nMCParticles;

  std::vector<int> * event_BX_for_PU = nullptr;
  std::vector<float> * event_PU = nullptr;

  eventDetailsStruct(TChain &inputChain, const bool& readMCCollections, const bool& calculateShiftedDistributions, const bool& savePUWeights) {
    if (savePUWeights) {
      inputChain.SetBranchStatus("puBX", 1);
      inputChain.SetBranchAddress("puBX", &(event_BX_for_PU));
      inputChain.SetBranchStatus("puTrue", 1);
      inputChain.SetBranchAddress("puTrue", &(event_PU));
    }
    inputChain.SetBranchStatus("HLTPho", 1);
    inputChain.SetBranchAddress("HLTPho", &(HLTPhotonBits));
    inputChain.SetBranchStatus("HLTJet", 1);
    inputChain.SetBranchAddress("HLTJet", &(HLTJetBits));
    inputChain.SetBranchStatus("rho", 1);
    inputChain.SetBranchAddress("rho", &(eventRho));
    inputChain.SetBranchStatus("nPho", 1);
    inputChain.SetBranchAddress("nPho", &(nPhotons));
    inputChain.SetBranchStatus("nJet", 1);
    inputChain.SetBranchAddress("nJet", &(nJets));
    inputChain.SetBranchStatus("nEle", 1);
    inputChain.SetBranchAddress("nEle", &(nElectrons));
    inputChain.SetBranchStatus("nMu", 1);
    inputChain.SetBranchAddress("nMu", &(nMuons));
    inputChain.SetBranchStatus("pfMET", 1);
    inputChain.SetBranchAddress("pfMET", &(PFMET));
    if (calculateShiftedDistributions) {
      inputChain.SetBranchStatus("pfMETPhi", 1);
      inputChain.SetBranchAddress("pfMETPhi", &(PFMET_phi));
      inputChain.SetBranchStatus("pfMET_T1UESDo", 1);
      inputChain.SetBranchAddress("pfMET_T1UESDo", &(PFMET_UnclusteredDown));
      inputChain.SetBranchStatus("pfMET_T1UESUp", 1);
      inputChain.SetBranchAddress("pfMET_T1UESUp", &(PFMET_UnclusteredUp));
      inputChain.SetBranchStatus("pfMET_T1JERDo", 1);
      inputChain.SetBranchAddress("pfMET_T1JERDo", &(PFMET_JERDown));
      inputChain.SetBranchStatus("pfMET_T1JERUp", 1);
      inputChain.SetBranchAddress("pfMET_T1JERUp", &(PFMET_JERUp));
    }
    if (readMCCollections) {
      inputChain.SetBranchStatus("nMC", 1);
      inputChain.SetBranchAddress("nMC", &(nMCParticles));
    }
  }
};

struct MCCollectionStruct{
  std::vector<int> * MCPIDs = nullptr;
  std::vector<int> * MCMomPIDs = nullptr;
  std::vector<int> * MCStatuses = nullptr;
  std::vector<UShort_t> * MCStatusFlags = nullptr;
  std::vector<float> * MCMasses = nullptr;
  std::vector<float> * MCMomMasses = nullptr;
  std::vector<float> * MCEts = nullptr;
  std::vector<float> * MCEtas = nullptr;
  std::vector<float> * MCPhis = nullptr;

  MCCollectionStruct(TChain &inputChain, const bool& readMCCollections) {
    if (readMCCollections) {
      inputChain.SetBranchStatus("mcPID", 1);
      inputChain.SetBranchAddress("mcPID", &(MCPIDs));
      inputChain.SetBranchStatus("mcMomPID", 1);
      inputChain.SetBranchAddress("mcMomPID", &(MCMomPIDs));
      inputChain.SetBranchStatus("mcStatus", 1);
      inputChain.SetBranchAddress("mcStatus", &(MCStatuses));
      inputChain.SetBranchStatus("mcStatusFlag", 1);
      inputChain.SetBranchAddress("mcStatusFlag", &(MCStatusFlags));
      inputChain.SetBranchStatus("mcMass", 1);
      inputChain.SetBranchAddress("mcMass", &(MCMasses));
      inputChain.SetBranchStatus("mcMomMass", 1);
      inputChain.SetBranchAddress("mcMomMass", &(MCMomMasses));
      inputChain.SetBranchStatus("mcEt", 1);
      inputChain.SetBranchAddress("mcEt", &(MCEts));
      inputChain.SetBranchStatus("mcEta", 1);
      inputChain.SetBranchAddress("mcEta", &(MCEtas));
      inputChain.SetBranchStatus("mcPhi", 1);
      inputChain.SetBranchAddress("mcPhi", &(MCPhis));
    }
  }
};

struct photonsCollectionStruct{
  std::vector<float> * pT = nullptr;
  std::vector<float> * eta = nullptr;
  std::vector<float> * mva = nullptr;
  std::vector<float> * phi = nullptr;
  std::vector<float> * HOverE = nullptr;
  std::vector<float> * sigmaIEtaIEta = nullptr;
  std::vector<float> * PFChargedIsolationUncorrected = nullptr;
  std::vector<float> * PFNeutralIsolationUncorrected = nullptr;
  std::vector<float> * PFPhotonIsolationUncorrected = nullptr;
  std::vector<UShort_t> * ID = nullptr;
  std::vector<int> * electronVeto = nullptr;
  std::vector<int> * hasPixelSeed = nullptr;
  std::vector<float> * energy = nullptr;
  std::vector<float> * R9 = nullptr;
  std::vector<float> * ecalClusIso = nullptr;
  std::vector<float> * trkIso = nullptr;

  photonsCollectionStruct(TChain &inputChain) {
    inputChain.SetBranchStatus("phoEt", 1);
    inputChain.SetBranchAddress("phoEt", &(pT));
    inputChain.SetBranchStatus("phoEta", 1);
    inputChain.SetBranchAddress("phoEta", &(eta));
    inputChain.SetBranchStatus("phoIDMVA", 1);
    inputChain.SetBranchAddress("phoIDMVA", &(mva));
    inputChain.SetBranchStatus("phoPhi", 1);
    inputChain.SetBranchAddress("phoPhi", &(phi));
    inputChain.SetBranchStatus("phoHoverE", 1);
    inputChain.SetBranchAddress("phoHoverE", &(HOverE));
    inputChain.SetBranchStatus("phoSigmaIEtaIEtaFull5x5", 1);
    inputChain.SetBranchAddress("phoSigmaIEtaIEtaFull5x5", &(sigmaIEtaIEta));
    inputChain.SetBranchStatus("phoPFChIso", 1);
    inputChain.SetBranchAddress("phoPFChIso", &(PFChargedIsolationUncorrected));
    inputChain.SetBranchStatus("phoPFNeuIso", 1);
    inputChain.SetBranchAddress("phoPFNeuIso", &(PFNeutralIsolationUncorrected));
    inputChain.SetBranchStatus("phoPFPhoIso", 1);
    inputChain.SetBranchAddress("phoPFPhoIso", &(PFPhotonIsolationUncorrected));
    inputChain.SetBranchStatus("phoIDbit", 1);
    inputChain.SetBranchAddress("phoIDbit", &(ID));
    inputChain.SetBranchStatus("phoEleVeto", 1);
    inputChain.SetBranchAddress("phoEleVeto", &(electronVeto));
    inputChain.SetBranchStatus("phohasPixelSeed", 1);
    inputChain.SetBranchAddress("phohasPixelSeed", &(hasPixelSeed));
    inputChain.SetBranchStatus("phoE", 1);
    inputChain.SetBranchAddress("phoE", &(energy));
    inputChain.SetBranchStatus("phoR9Full5x5", 1);
    inputChain.SetBranchAddress("phoR9Full5x5", &(R9));
    inputChain.SetBranchStatus("phoPFEcalClusIso", 1);
    inputChain.SetBranchAddress("phoPFEcalClusIso", &(ecalClusIso));
    inputChain.SetBranchStatus("phoTrkIso", 1);
    inputChain.SetBranchAddress("phoTrkIso", &(trkIso));
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
  std::vector<float> * jetGenPT = nullptr;
  std::vector<float> * jetGenEta = nullptr;
  std::vector<float> * jetGenPhi = nullptr;
  std::vector<int> * jetGenPartonID = nullptr;
  std::vector<int> * jetGenPartonMomID = nullptr;

  jetsCollectionStruct(TChain &inputChain, const bool& saveMCObjects, const bool& calculateShiftedDistributions) {
    inputChain.SetBranchStatus("jetPt", 1);
    inputChain.SetBranchAddress("jetPt", &(pT));
    inputChain.SetBranchStatus("jetEta", 1);
    inputChain.SetBranchAddress("jetEta", &(eta));
    inputChain.SetBranchStatus("jetPhi", 1);
    inputChain.SetBranchAddress("jetPhi", &(phi));
    inputChain.SetBranchStatus("jetPUID", 1);
    inputChain.SetBranchAddress("jetPUID", &(PUID));
    inputChain.SetBranchStatus("jetPFLooseId", 1);
    inputChain.SetBranchAddress("jetPFLooseId", &(looseID));
    inputChain.SetBranchStatus("jetID", 1);
    inputChain.SetBranchAddress("jetID", &(ID));
    if (saveMCObjects) {
      inputChain.SetBranchStatus("jetGenJetPt", 1);
      inputChain.SetBranchAddress("jetGenJetPt", &(jetGenPT));
      inputChain.SetBranchStatus("jetGenJetEta", 1);
      inputChain.SetBranchAddress("jetGenJetEta", &(jetGenEta));
      inputChain.SetBranchStatus("jetGenJetPhi", 1);
      inputChain.SetBranchAddress("jetGenJetPhi", &(jetGenPhi));
      inputChain.SetBranchStatus("jetGenPartonID", 1);
      inputChain.SetBranchAddress("jetGenPartonID", &(jetGenPartonID));
      inputChain.SetBranchStatus("jetGenPartonMomID", 1);
      inputChain.SetBranchAddress("jetGenPartonMomID", &(jetGenPartonMomID));
    }
    if (calculateShiftedDistributions) {
      inputChain.SetBranchStatus("jetJECUnc", 1);
      inputChain.SetBranchAddress("jetJECUnc", &(JECUncertainty));
    }
  }
};

#endif
