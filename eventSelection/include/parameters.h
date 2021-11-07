#ifndef H_PARAMETERS
#define H_PARAMETERS

#include <vector>
#include <map>
#include <iostream>
#include <sstream>

#include "miscDataStructures.h"
#include "selectionCriteria.h"
#include "constants.h"

#include "TROOT.h"
#include "TFile.h"
#include "TH2F.h"

struct parametersStruct {
  const float photonBarrelEtaCut = 1.442f;
  const float photonEndcapEtaLow = 1.52f;
  const float photonEndcapEtaHigh = 2.4f;
  const float jetEtaCut = 2.4f;
  const float jetpTCut = 30.f;
  const float jetPUIDThreshold = 0.61f;
  const float deltaRScale_truthMatching = 0.05f;
  const float deltaRScale_jetPhotonDistance = 0.4f;
  const float HTCut = 60.0f;
  const float preNormalizationBuffer = 200.0f;
  const int jetCandidateStatusConstraint = 23;

  /*
    bit convention:
    bit 0: fromHardProcessFinalState
    bit 1: isPromptFinalState
    bit 2: isHardProcess

    7 = 111, satisfies all three conditions
    4 = 100, isHardProcess (but no constraint on the other conditions)
    2 = 010, isPromptFinalState (but no constraints on other conditions)
  */
  const UShort_t MCStatusFlagBitMask = static_cast<UShort_t>(7u);
  const UShort_t MCStatusFlagBitMask_promptOnly = static_cast<UShort_t>(2u);

  triggerType HLT_triggerType;
  int HLTBit_photon;
  int HLTBit_jet;

  float pTCutSubLeading;
  float pTCutLeading;
  float invariantMassCut;
  photonQualityCutsStruct photonQualityCutsBarrel;
  photonQualityCutsStruct photonQualityCutsEndcap;
  EAValuesStruct effectiveAreas[7];
  bool calculatePrefiringWeights;
  TFile* sourceFile_prefiringEfficiencyMap;
  TH2F* prefiringEfficiencyMap;
  TFile* sourceFile_photonMCScaleFactorsMap_medium;
  TFile* sourceFile_photonMCScaleFactorsMap_loose;
  std::map<photonType, TH2F*> photonMCScaleFactorsMaps;
  TFile* sourceFile_PUWeights_handle;
  TH1D* PUWeights;
  void tuneParameters(const int& year, const bool& calculateMCScaleFactorWeights, const bool& savePUWeights, const std::string& PUWeightsPathWithXRDPrefix, const std::string& selectionType) {
    if (savePUWeights) {
      assert(PUWeightsPathWithXRDPrefix != "/dev/null");
      sourceFile_PUWeights_handle = TFile::Open(PUWeightsPathWithXRDPrefix.c_str(), "READ");
      assert((sourceFile_PUWeights_handle->IsOpen()) && (!(sourceFile_PUWeights_handle->IsZombie())));
      sourceFile_PUWeights_handle->GetObject("pileupWeights", PUWeights);
      assert(PUWeights != nullptr);
      std::cout << "Read in pileup weights from file: " << PUWeightsPathWithXRDPrefix << std::endl;
    }
    if (year == 2018) { // very similar to 2017. Differences: no ECAL prefiring in 2018, and different scale factors.
      /* "interesting" photon bits: */
      /* 16: HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v */
      /* 37: HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v */
      /* 22: HLT_DoublePhoton70_v */
      /* 10: HLT_Photon200_v */

      /* "interesting" jet bits: */
      /* 37: HLT_PFHT1050_v */

      HLT_triggerType = triggerType::nTriggerTypes;
      HLTBit_photon = -1;
      HLTBit_jet = -1;
      if ((selectionType == "data") ||
          (selectionType == "MC_hgg") ||
          (std::regex_match(selectionType, std::regex("^MC_DiPhotonJets$"))) ||
	  (std::regex_match(selectionType, std::regex("^MC_EMEnrichedGJetPt[0-9]*_[0-9]*$"))) ||
	  (std::regex_match(selectionType, std::regex("^MC_HighHTQCD[0-9]*_[0-9]*$"))) ||
	  (std::regex_match(selectionType, std::regex("^MC_GJetHT[0-9]*_[0-9]*$")))) {
        HLT_triggerType = triggerType::photon;
        HLTBit_photon = 37;
      }
      else if ((selectionType == "data_singlephoton") ||
               (std::regex_match(selectionType, std::regex("^MC_DiPhotonJets_singlephoton$"))) ||
	       (std::regex_match(selectionType, std::regex("^MC_EMEnrichedGJetPt[0-9]*_singlephoton_[0-9]*$"))) ||
	       (std::regex_match(selectionType, std::regex("^MC_HighHTQCD[0-9]*_singlephoton_[0-9]*$"))) ||
	       (std::regex_match(selectionType, std::regex("^MC_GJetHT[0-9]*_singlephoton_[0-9]*$")))) {
        HLT_triggerType = triggerType::photon;
        HLTBit_photon = 10;
      }
      pTCutSubLeading = 25.0f;
      pTCutLeading = 35.0f;
      invariantMassCut = 90.0f;

      photonQualityCutsBarrel = photonQualityCutsStruct(0.02197f, 0.04596f, 0.01015f, 0.0106f, 1.141f, 1.694f, 6.0f, 1.189f, 0.01512f, 0.00002259f, 24.032f, 0.01512f, 0.00002259f, 2.080f, 0.004017f, 2.876f, 0.004017f);
      photonQualityCutsEndcap = photonQualityCutsStruct(0.03260f, 0.05900f, 0.02720f, 0.0272f, 1.051f, 2.089f, 6.0f, 2.718f, 0.01170f, 0.00002300f, 19.722f, 0.01170f, 0.00002300f, 3.867f, 0.003700f, 4.162f, 0.003700f);

      effectiveAreas[0] = EAValuesStruct(1.0f, 0.0112f, 0.0668f, 0.1113f);
      effectiveAreas[1] = EAValuesStruct(1.479f, 0.0108f, 0.1054f, 0.0953f);
      effectiveAreas[2] = EAValuesStruct(2.0f, 0.0106f, 0.0786f, 0.0619f);
      effectiveAreas[3] = EAValuesStruct(2.2f, 0.01002f, 0.0233f, 0.0837f);
      effectiveAreas[4] = EAValuesStruct(2.3f, 0.0098f, 0.0078f, 0.1070f);
      effectiveAreas[5] = EAValuesStruct(2.4f, 0.0089f, 0.0028f, 0.1212f);
      effectiveAreas[6] = EAValuesStruct(10.0f, 0.0087f, 0.0137f, 0.1466f);

      calculatePrefiringWeights = false;

      if (calculateMCScaleFactorWeights) {
        sourceFile_photonMCScaleFactorsMap_medium = TFile::Open("eventSelection/data/2018_PhotonsMedium.root", "READ");
        if (!(sourceFile_photonMCScaleFactorsMap_medium->IsOpen()) || sourceFile_photonMCScaleFactorsMap_medium->IsZombie()) {
          std::cout << "ERROR: Unable to open file with path: eventSelection/data/2018_PhotonsMedium.root" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        TH2F *photonMCScaleFactorsMapMedium = new TH2F();
        sourceFile_photonMCScaleFactorsMap_medium->GetObject("EGamma_SF2D", photonMCScaleFactorsMapMedium);
        assert(photonMCScaleFactorsMapMedium != nullptr);
        photonMCScaleFactorsMaps[photonType::medium] = photonMCScaleFactorsMapMedium;
        std::cout << "Opened medium photon MC scale factors map for 2018" << std::endl;

        sourceFile_photonMCScaleFactorsMap_loose = TFile::Open("eventSelection/data/2018_PhotonsLoose.root", "READ");
        if (!(sourceFile_photonMCScaleFactorsMap_loose->IsOpen()) || sourceFile_photonMCScaleFactorsMap_loose->IsZombie()) {
          std::cout << "ERROR: Unable to open file with path: eventSelection/data/2018_PhotonsLoose.root" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        TH2F *photonMCScaleFactorsMapLoose = new TH2F();
        sourceFile_photonMCScaleFactorsMap_loose->GetObject("EGamma_SF2D", photonMCScaleFactorsMapLoose);
        assert(photonMCScaleFactorsMapLoose != nullptr);
        photonMCScaleFactorsMaps[photonType::vetoed] = photonMCScaleFactorsMapLoose;
        std::cout << "Opened loose photon MC scale factors map for 2018" << std::endl;
      }
    }
    if (year == 2017) {
      /* "interesting" photon bits: */
      /* 16: HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v */
      /* 37: HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_PixelVeto_Mass55_v */
      /* 22: HLT_DoublePhoton70_v */
      /* 10: HLT_Photon200_v */
      /* 2: HLT_Photon33_v */

      /* "interesting" jet bits: */
      /* 37: HLT_PFHT1050_v */

      HLT_triggerType = triggerType::nTriggerTypes;
      HLTBit_photon = -1;
      HLTBit_jet = -1;
      if ((selectionType == "data") ||
          (selectionType == "MC_hgg") ||
          (std::regex_match(selectionType, std::regex("^MC_DiPhotonJets$"))) ||
	  (std::regex_match(selectionType, std::regex("^MC_EMEnrichedGJetPt[0-9]*_[0-9]*$"))) ||
	  (std::regex_match(selectionType, std::regex("^MC_HighHTQCD[0-9]*_[0-9]*$"))) ||
	  (std::regex_match(selectionType, std::regex("^MC_GJetHT[0-9]*_[0-9]*$")))) {
        HLT_triggerType = triggerType::photon;
        HLTBit_photon = 37;
      }
      else if ((selectionType == "data_singlephoton") ||
               (std::regex_match(selectionType, std::regex("^MC_DiPhotonJets_singlephoton$"))) ||
	       (std::regex_match(selectionType, std::regex("^MC_EMEnrichedGJetPt[0-9]*_singlephoton_[0-9]*$"))) ||
	       (std::regex_match(selectionType, std::regex("^MC_HighHTQCD[0-9]*_singlephoton_[0-9]*$"))) ||
	       (std::regex_match(selectionType, std::regex("^MC_GJetHT[0-9]*_singlephoton_[0-9]*$")))) {
        HLT_triggerType = triggerType::photon;
        HLTBit_photon = 10;
      }
      pTCutSubLeading = 25.0f;
      pTCutLeading = 35.0f;
      invariantMassCut = 90.0f;

      photonQualityCutsBarrel = photonQualityCutsStruct(0.02197f, 0.04596f, 0.01015f, 0.0106f, 1.141f, 1.694f, 6.0f, 1.189f, 0.01512f, 0.00002259f, 24.032f, 0.01512f, 0.00002259f, 2.080f, 0.004017f, 2.876f, 0.004017f);
      photonQualityCutsEndcap = photonQualityCutsStruct(0.03260f, 0.05900f, 0.02720f, 0.0272f, 1.051f, 2.089f, 6.0f, 2.718f, 0.01170f, 0.00002300f, 19.722f, 0.01170f, 0.00002300f, 3.867f, 0.003700f, 4.162f, 0.003700f);

      effectiveAreas[0] = EAValuesStruct(1.0f, 0.0112f, 0.0668f, 0.1113f);
      effectiveAreas[1] = EAValuesStruct(1.479f, 0.0108f, 0.1054f, 0.0953f);
      effectiveAreas[2] = EAValuesStruct(2.0f, 0.0106f, 0.0786f, 0.0619f);
      effectiveAreas[3] = EAValuesStruct(2.2f, 0.01002f, 0.0233f, 0.0837f);
      effectiveAreas[4] = EAValuesStruct(2.3f, 0.0098f, 0.0078f, 0.1070f);
      effectiveAreas[5] = EAValuesStruct(2.4f, 0.0089f, 0.0028f, 0.1212f);
      effectiveAreas[6] = EAValuesStruct(10.0f, 0.0087f, 0.0137f, 0.1466f);

      calculatePrefiringWeights = true;
      sourceFile_prefiringEfficiencyMap = TFile::Open("eventSelection/data/L1prefiring_jetpt_2017BtoF.root", "READ");
      if (!(sourceFile_prefiringEfficiencyMap->IsOpen()) || sourceFile_prefiringEfficiencyMap->IsZombie()) {
        std::cout << "ERROR: Unable to open file with path: eventSelection/data/L1prefiring_jetpt_2017BtoF.root" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      sourceFile_prefiringEfficiencyMap->GetObject("L1prefiring_jetpt_2017BtoF", prefiringEfficiencyMap);
      if (prefiringEfficiencyMap) std::cout << "Opened prefiring efficiency map for 2017" << std::endl;
      else {
        std::cout << "ERROR: Unable to open histogram with path: L1prefiring_jetpt_2017BtoF" << std::endl;
        std::exit(EXIT_FAILURE);
      }

      if (calculateMCScaleFactorWeights) {
        sourceFile_photonMCScaleFactorsMap_medium = TFile::Open("eventSelection/data/2017_PhotonsMedium.root", "READ");
        if (!(sourceFile_photonMCScaleFactorsMap_medium->IsOpen()) || sourceFile_photonMCScaleFactorsMap_medium->IsZombie()) {
          std::cout << "ERROR: Unable to open file with path: eventSelection/data/2017_PhotonsMedium.root" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        TH2F *photonMCScaleFactorsMapMedium = new TH2F();
        sourceFile_photonMCScaleFactorsMap_medium->GetObject("EGamma_SF2D", photonMCScaleFactorsMapMedium);
        assert(photonMCScaleFactorsMapMedium != nullptr);
        photonMCScaleFactorsMaps[photonType::medium] = photonMCScaleFactorsMapMedium;
        std::cout << "Opened medium photon MC scale factors map for 2017" << std::endl;

        sourceFile_photonMCScaleFactorsMap_loose = TFile::Open("eventSelection/data/2017_PhotonsLoose.root", "READ");
        if (!(sourceFile_photonMCScaleFactorsMap_loose->IsOpen()) || sourceFile_photonMCScaleFactorsMap_loose->IsZombie()) {
          std::cout << "ERROR: Unable to open file with path: eventSelection/data/2017_PhotonsLoose.root" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        TH2F *photonMCScaleFactorsMapLoose = new TH2F();
        sourceFile_photonMCScaleFactorsMap_loose->GetObject("EGamma_SF2D", photonMCScaleFactorsMapLoose);
        assert(photonMCScaleFactorsMapLoose != nullptr);
        photonMCScaleFactorsMaps[photonType::vetoed] = photonMCScaleFactorsMapLoose;
        std::cout << "Opened loose photon MC scale factors map for 2017" << std::endl;
      }
    }
    else if (year == 2016) {
      /* "interesting" photon bits: */
      /* 16: HLT_Diphoton30PV_18PV_R9Id_AND_IsoCaloId_AND_HE_R9Id_DoublePixelVeto_Mass55_v */
      /* 22: HLT_DoublePhoton60_v */
      /* 7: HLT_Photon175_v */
      /* 1: HLT_Photon30_v */
      /* 12: HLT_Photon165_HE10_v */

      /* "interesting" jet bits: */
      /* 33: HLT_PFHT900_v */

      HLT_triggerType = triggerType::nTriggerTypes;
      HLTBit_photon = -1;
      HLTBit_jet = -1;
      if ((selectionType == "data") ||
          (selectionType == "MC_hgg") ||
          (std::regex_match(selectionType, std::regex("^MC_DiPhotonJets$"))) ||
	  (std::regex_match(selectionType, std::regex("^MC_EMEnrichedGJetPt[0-9]*_[0-9]*$"))) ||
	  (std::regex_match(selectionType, std::regex("^MC_HighHTQCD[0-9]*_[0-9]*$"))) ||
	  (std::regex_match(selectionType, std::regex("^MC_GJetHT[0-9]*_[0-9]*$")))) {
        HLT_triggerType = triggerType::photon;
        HLTBit_photon = 16;
      }
      else if ((selectionType == "data_singlephoton") ||
               (std::regex_match(selectionType, std::regex("^MC_DiPhotonJets_singlephoton$"))) ||
	       (std::regex_match(selectionType, std::regex("^MC_EMEnrichedGJetPt[0-9]*_singlephoton_[0-9]*$"))) ||
	       (std::regex_match(selectionType, std::regex("^MC_HighHTQCD[0-9]*_singlephoton_[0-9]*$"))) ||
	       (std::regex_match(selectionType, std::regex("^MC_GJetHT[0-9]*_singlephoton_[0-9]*$")))) {
        HLT_triggerType = triggerType::photon;
        HLTBit_photon = 7;
	/* HLT_triggerType = triggerType::jet; */
      }
      pTCutSubLeading = 25.0f;
      pTCutLeading = 35.0f;
      invariantMassCut = 90.0f;

      photonQualityCutsBarrel = photonQualityCutsStruct(0.0396f, 0.0597f, 0.01022f, 0.01031f, 0.441f, 1.295f, 6.0f, 2.725f, 0.0148f, 0.000017f, 10.91f, 0.0148f, 0.000017f, 2.571f, 0.0047f, 3.630f, 0.0047f);
      photonQualityCutsEndcap = photonQualityCutsStruct(0.0219f, 0.0481f, 0.03001f, 0.03013f, 0.442f, 1.011f, 6.0f, 1.715f, 0.0163f, 0.000014f, 5.931f, 0.0163f, 0.000014f, 3.863f, 0.0034f, 6.641f, 0.0034f);

      effectiveAreas[0] = EAValuesStruct(1.0f, 0.036f, 0.0597f, 0.121f);
      effectiveAreas[1] = EAValuesStruct(1.479f, 0.0377f, 0.0807f, 0.1107f);
      effectiveAreas[2] = EAValuesStruct(2.0f, 0.0306f, 0.0629f, 0.0699f);
      effectiveAreas[3] = EAValuesStruct(2.2f, 0.0283f, 0.0197f, 0.1056f);
      effectiveAreas[4] = EAValuesStruct(2.3f, 0.0254f, 0.0184f, 0.1457f);
      effectiveAreas[5] = EAValuesStruct(2.4f, 0.0217f, 0.0284f, 0.1719f);
      effectiveAreas[6] = EAValuesStruct(10.0f, 0.0167f, 0.0591f, 0.1998f);

      calculatePrefiringWeights = true;
      sourceFile_prefiringEfficiencyMap = TFile::Open("eventSelection/data/L1prefiring_jetpt_2016BtoH.root", "READ");
      if (!(sourceFile_prefiringEfficiencyMap->IsOpen()) || sourceFile_prefiringEfficiencyMap->IsZombie()) {
        std::cout << "ERROR: Unable to open file with path: eventSelection/data/L1prefiring_jetpt_2016BtoH.root" << std::endl;
        std::exit(EXIT_FAILURE);
      }
      sourceFile_prefiringEfficiencyMap->GetObject("L1prefiring_jetpt_2016BtoH", prefiringEfficiencyMap);
      if (prefiringEfficiencyMap) std::cout << "Opened prefiring efficiency map for 2016" << std::endl;
      else {
        std::cout << "ERROR: Unable to open histogram with path: L1prefiring_jetpt_2016BtoH" << std::endl;
        std::exit(EXIT_FAILURE);
      }

      if (calculateMCScaleFactorWeights) {
        sourceFile_photonMCScaleFactorsMap_medium = TFile::Open("eventSelection/data/80X_2016_Medium_photons.root", "READ");
        if (!(sourceFile_photonMCScaleFactorsMap_medium->IsOpen()) || sourceFile_photonMCScaleFactorsMap_medium->IsZombie()) {
          std::cout << "ERROR: Unable to open file with path: eventSelection/data/80X_2016_Medium_photons.root" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        TH2F *photonMCScaleFactorsMapMedium = new TH2F();
        sourceFile_photonMCScaleFactorsMap_medium->GetObject("EGamma_SF2D", photonMCScaleFactorsMapMedium);
        assert(photonMCScaleFactorsMapMedium != nullptr);
        photonMCScaleFactorsMaps[photonType::medium] = photonMCScaleFactorsMapMedium;
        std::cout << "Opened medium photon MC scale factors map for 2016" << std::endl;

        sourceFile_photonMCScaleFactorsMap_loose = TFile::Open("eventSelection/data/80X_2016_Loose_photons.root", "READ");
        if (!(sourceFile_photonMCScaleFactorsMap_loose->IsOpen()) || sourceFile_photonMCScaleFactorsMap_loose->IsZombie()) {
          std::cout << "ERROR: Unable to open file with path: eventSelection/data/80X_2016_Loose_photons.root" << std::endl;
          std::exit(EXIT_FAILURE);
        }
        TH2F *photonMCScaleFactorsMapLoose = new TH2F();
        sourceFile_photonMCScaleFactorsMap_loose->GetObject("EGamma_SF2D", photonMCScaleFactorsMapLoose);
        assert(photonMCScaleFactorsMapLoose != nullptr);
        photonMCScaleFactorsMaps[photonType::vetoed] = photonMCScaleFactorsMapLoose;
        std::cout << "Opened loose photon MC scale factors map for 2016" << std::endl;
      }
    }
  }
  friend std::ostream& operator<< (std::ostream& out, const parametersStruct& parameters) {
    out << "Photon cuts:" << std::endl
        << "pT_SubLeading: " << parameters.pTCutSubLeading << ", "
        << "pT_Leading: " << parameters.pTCutLeading << ", "
        << "eta cut, barrel: " << parameters.photonBarrelEtaCut << ", "
        << "eta endcap, low: " << parameters.photonEndcapEtaLow << ", "
        << "eta endcap, high: " << parameters.photonEndcapEtaHigh << ", "
        << "Quality cuts (barrel): " << parameters.photonQualityCutsBarrel << ", "
        << "Quality cuts (endcap): " << parameters.photonQualityCutsEndcap << std::endl;

    out << "Effective Areas: " << std::endl;
    for (unsigned int areaCounter = 0; areaCounter < 7; ++areaCounter) {
      out << parameters.effectiveAreas[areaCounter] << std::endl;
    }

    out << "Invariant mass cut: " << parameters.invariantMassCut << std::endl;

    out << "Jet cuts:" << std::endl
        << "pT: " << parameters.jetpTCut << ", "
        << "eta: " << parameters.jetEtaCut << ", "
        << "PUID: " << parameters.jetPUIDThreshold << ", "
        << "deltaRScale_truthMatching: " << parameters.deltaRScale_truthMatching << ", "
        << "deltaRScale_jetPhotonDistance: " << parameters.deltaRScale_jetPhotonDistance << std::endl;

    out << "Event cuts:" << std::endl
        << "HLT bit index, photon: " << parameters.HLTBit_photon << ", "
        << "HLT bit index, jet: " << parameters.HLTBit_jet << ", "
        << "HT Cut: " << parameters.HTCut;
    return out;
  }
};
#endif
