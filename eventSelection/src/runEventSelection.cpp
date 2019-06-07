#include "../include/runEventSelection.h"

float getRhoCorrectedIsolation(const float& uncorrectedIsolation, const PFTypesForEA& PFType, const float& absEta, const float& rho, EAValuesStruct (&EAValues)[7]) {
  float effectiveArea = 0.;
  unsigned int areaCounter = 0;
  while (areaCounter < 7) {
    if (absEta <= (EAValues[areaCounter]).regionUpperBound) {
      effectiveArea = (EAValues[areaCounter]).getEffectiveArea(PFType);
      break;
    }
    ++areaCounter;
  }
  if (areaCounter == 7) {
    std::cout << "ERROR: Area counter value = " << areaCounter << "; getRhoCorrectedIsolation called for eta = " << absEta << ", which is above all available upper bounds." << std::endl;
    std::exit(EXIT_FAILURE);
  }
  float correctedIsolationTest = uncorrectedIsolation - rho*effectiveArea;
  return ((correctedIsolationTest > 0.0) ? correctedIsolationTest : 0.0);
}

eventWeightsStruct findMCScaleFactors(const float& eta, const float& pT, TH2F* MCScaleFactorsMap) {
  eventWeightsStruct MCScaleFactors;
  float binContent = MCScaleFactorsMap->GetBinContent(MCScaleFactorsMap->FindFixBin(eta, pT));
  if (binContent > 0.) {
    float binError = MCScaleFactorsMap->GetBinError(MCScaleFactorsMap->FindFixBin(eta, pT));
    float scaleFactorNominal = binContent;
    float scaleFactorUp = binContent + binError;
    float scaleFactorDown = std::max(binContent - binError, 0.0f);
    if (scaleFactorDown > scaleFactorUp) {// sanity check
      std::cout << "ERROR: Garbage scale factors." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    MCScaleFactors = eventWeightsStruct(scaleFactorNominal, scaleFactorDown, scaleFactorUp);
  }
  else {
    MCScaleFactors = eventWeightsStruct(1.0f, 1.0f, 1.0f);
  }
  return MCScaleFactors;
}

photonExaminationResultsStruct examinePhoton(optionsStruct &options, parametersStruct &parameters, // countersStruct &counters, 
                                             const float& rho, const photonsCollectionStruct& photonsCollection, const int& photonIndex, std::vector<angularVariablesStruct> &selectedTruePhotonAngles) {
  photonExaminationResultsStruct results;
  photonProperties& properties = results.pho_properties;
  eventWeightsStruct& scaleFactors = results.MCScaleFactors;

  std::map<mediumPhotonCriterion, bool> medium_bits;
  std::map<fakePhotonCriterion, bool> fake_bits;

  // Kinematic cuts
  properties[photonProperty::eta] = (photonsCollection.eta)->at(photonIndex);
  properties[photonProperty::phi] = (photonsCollection.phi)->at(photonIndex);
  angularVariablesStruct photonAngle(properties[photonProperty::eta], properties[photonProperty::phi]);
  properties[photonProperty::deltaR_nearestTruePhoton] = photonAngle.getMinDeltaR(selectedTruePhotonAngles);
  float absEta = std::fabs(properties[photonProperty::eta]);
  bool passesEta = (absEta < parameters.photonBarrelEtaCut);
  medium_bits[mediumPhotonCriterion::eta] = passesEta;
  fake_bits[fakePhotonCriterion::eta] = passesEta;

  properties[photonProperty::pT] = std::fabs((photonsCollection.pT)->at(photonIndex));
  bool passesPT = (properties[photonProperty::pT] > parameters.pTCutSubLeading);
  medium_bits[mediumPhotonCriterion::pT] = passesPT;
  fake_bits[fakePhotonCriterion::pT] = passesPT;

  // Electron veto
  bool passesConvSafeVeto = (((photonsCollection.electronVeto)->at(photonIndex)) == (Int_t)(true));
  medium_bits[mediumPhotonCriterion::conversionSafeElectronVeto] = passesConvSafeVeto;
  fake_bits[fakePhotonCriterion::conversionSafeElectronVeto] = passesConvSafeVeto;

  // Quality cuts
  photonQualityCutsStruct* qualityCuts = &(parameters.photonQualityCutsBarrel);
  if (absEta > parameters.photonBarrelEtaCut) qualityCuts = &(parameters.photonQualityCutsEndcap);

  properties[photonProperty::hOverE] = (photonsCollection.HOverE)->at(photonIndex);
  bool passesHOverE = (properties[photonProperty::hOverE] < qualityCuts->towerHOverE);
  medium_bits[mediumPhotonCriterion::hOverE] = passesHOverE;
  fake_bits[fakePhotonCriterion::hOverE] = passesHOverE;

  float pTDependentNeutralIsolationCut = (qualityCuts->neutralIsolation).getPolynomialValue(properties[photonProperty::pT]);
  properties[photonProperty::rhoCorrectedNeutralIsolation] = getRhoCorrectedIsolation(((photonsCollection.PFNeutralIsolationUncorrected)->at(photonIndex)), PFTypesForEA::neutralHadron, absEta, rho, parameters.effectiveAreas);
  bool passesNeutralIsolation = (properties[photonProperty::rhoCorrectedNeutralIsolation] < pTDependentNeutralIsolationCut);
  medium_bits[mediumPhotonCriterion::neutralIsolation] = passesNeutralIsolation;
  fake_bits[fakePhotonCriterion::neutralIsolation] = passesNeutralIsolation;

  float pTDependentPhotonIsolationCut = (qualityCuts->photonIsolation).getPolynomialValue(properties[photonProperty::pT]);
  properties[photonProperty::rhoCorrectedPhotonIsolation] = getRhoCorrectedIsolation(((photonsCollection.PFPhotonIsolationUncorrected)->at(photonIndex)), PFTypesForEA::photon, absEta, rho, parameters.effectiveAreas);
  bool passesPhotonIsolation = (properties[photonProperty::rhoCorrectedPhotonIsolation] < pTDependentPhotonIsolationCut);
  medium_bits[mediumPhotonCriterion::photonIsolation] = passesPhotonIsolation;
  fake_bits[fakePhotonCriterion::photonIsolation] = passesPhotonIsolation;

  properties[photonProperty::rhoCorrectedChargedIsolation] = getRhoCorrectedIsolation(((photonsCollection.PFChargedIsolationUncorrected)->at(photonIndex)), PFTypesForEA::chargedHadron, absEta, rho, parameters.effectiveAreas);
  bool passesMedium_chargedIsolationCut = (properties[photonProperty::rhoCorrectedChargedIsolation] < qualityCuts->chargedIsolation);
  bool passesLoose_chargedIsolationCut = (properties[photonProperty::rhoCorrectedChargedIsolation] < qualityCuts->chargedIsolationLoose);
  medium_bits[mediumPhotonCriterion::chargedIsolation] = passesMedium_chargedIsolationCut;
  fake_bits[fakePhotonCriterion::chargedIsolationLoose] = passesLoose_chargedIsolationCut;

  properties[photonProperty::sigmaIEtaIEta] = ((photonsCollection.sigmaIEtaIEta)->at(photonIndex));
  bool passesMedium_sigmaIEtaIEtaCut = (properties[photonProperty::sigmaIEtaIEta] < qualityCuts->sigmaIEtaIEta);
  bool passesLoose_sigmaIEtaIEtaCut = (properties[photonProperty::sigmaIEtaIEta] < qualityCuts->sigmaIEtaIEtaLoose);
  medium_bits[mediumPhotonCriterion::sigmaIEtaIEta] = passesMedium_sigmaIEtaIEtaCut;
  fake_bits[fakePhotonCriterion::sigmaIEtaIEtaLoose] = passesLoose_sigmaIEtaIEtaCut;

  fake_bits[fakePhotonCriterion::notTightChIsoAndSigmaIEtaIEta] = !(passesMedium_sigmaIEtaIEtaCut && passesMedium_chargedIsolationCut);

  assert(static_cast<int>(fake_bits.size()) == static_cast<int>(fakePhotonCriterion::nFakePhotonCriteria));
  int nFalseBits_fake = getNFalseBits(fake_bits);
  assert(static_cast<int>(medium_bits.size()) == static_cast<int>(mediumPhotonCriterion::nMediumPhotonCriteria));
  int nFalseBits_medium = getNFalseBits(medium_bits);
  results.isSelectedFake = (nFalseBits_fake == 0);
  results.isSelectedMedium = (nFalseBits_medium == 0);
  results.isMarginallyUnselectedFake = (nFalseBits_fake == 1);
  if (results.isMarginallyUnselectedFake) {
    results.marginallyUnselectedFakeCriterion = getFirstFalseCriterion(fake_bits);
  }
  results.isMarginallyUnselectedMedium = (nFalseBits_medium == 1);
  if (results.isMarginallyUnselectedMedium) {
    results.marginallyUnselectedMediumCriterion = getFirstFalseCriterion(medium_bits);
  }

  if (options.isMC && (results.isSelectedFake || results.isSelectedMedium)) {
    scaleFactors = findMCScaleFactors(((photonsCollection.eta)->at(photonIndex)), ((photonsCollection.pT)->at(photonIndex)), parameters.photonMCScaleFactorsMap);
  }

  results.energy = (photonsCollection.energy)->at(photonIndex);

  if ((results.isMarginallyUnselectedMedium || results.isMarginallyUnselectedFake) || (results.isSelectedMedium || results.isSelectedFake)) assert(static_cast<int>((results.pho_properties).size()) == static_cast<int>(photonProperty::nPhotonProperties));

  return results;
}

MCExaminationResultsStruct examineMCParticle(parametersStruct &parameters, // countersStruct &counters, 
                                             // const int& nMCParticles, 
                                             const MCCollectionStruct& MCCollection, const int& MCIndex) {
  MCExaminationResultsStruct MCExaminationResults;
  truthPhotonProperties& pho_properties = MCExaminationResults.truth_photon_properties;

  int particle_mcPID = (MCCollection.MCPIDs)->at(MCIndex);
  int particle_mcMomPID = (MCCollection.MCMomPIDs)->at(MCIndex);

  if ((particle_mcPID == parameters.PIDs.photon) && (particle_mcMomPID == parameters.PIDs.neutralino)) {
    if (passesBitMask((MCCollection.MCStatusFlags)->at(MCIndex), parameters.MCStatusFlagBitMask)) {
      MCExaminationResults.isPhotonWithNeutralinoMom = true;
      pho_properties[truthPhotonProperty::eta] = (MCCollection.MCEtas)->at(MCIndex);
      pho_properties[truthPhotonProperty::phi] = (MCCollection.MCPhis)->at(MCIndex);
      pho_properties[truthPhotonProperty::pT] = (MCCollection.MCEts)->at(MCIndex);
    }
  }
  if (particle_mcPID == parameters.PIDs.gluino) MCExaminationResults.gluinoMass = (MCCollection.MCMasses)->at(MCIndex);
  if (particle_mcMomPID == parameters.PIDs.neutralino) MCExaminationResults.neutralinoMass = (MCCollection.MCMomMasses)->at(MCIndex);

  if (MCExaminationResults.isPhotonWithNeutralinoMom) assert(static_cast<int>((MCExaminationResults.truth_photon_properties).size()) == static_cast<int>(truthPhotonProperty::nTruthPhotonProperties));
  return MCExaminationResults;
}

float getDiphotonInvariantMass(const std::vector<TLorentzVector>& selectedPhotonFourMomentaList) {
  TLorentzVector eventSum;
  for (const auto& selectedPhotonFourMomentum : selectedPhotonFourMomentaList) {
    eventSum += selectedPhotonFourMomentum;
  }
  return (eventSum.M());
}

eventWeightsStruct findPrefireWeights(const float& eta, const float& pT, TH2F* efficiencyMap) {
  eventWeightsStruct prefireWeights;
  float binContent = efficiencyMap->GetBinContent(efficiencyMap->FindFixBin(eta, pT));
  if (binContent > 0.) {
    float binError = efficiencyMap->GetBinError(efficiencyMap->FindFixBin(eta, pT));
    float prefiringFractionNominal = binContent;
    float prefiringFractionUp = std::min(binContent + binError, 1.0f);
    float prefiringFractionDown = std::max(binContent - binError, 0.0f);
    if (prefiringFractionDown > prefiringFractionUp) {// sanity check
      std::cout << "ERROR: Garbage prefiring fractions." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    prefireWeights = eventWeightsStruct(1.0f - prefiringFractionNominal, 1.0f - prefiringFractionUp, 1.0f - prefiringFractionDown); // probability down = 1 - prefiringFractionUp
  }
  else {
    prefireWeights = eventWeightsStruct(1.0f, 1.0f, 1.0f);
  }
  return prefireWeights;
}

jetExaminationResultsStruct examineJet(optionsStruct &options, parametersStruct &parameters, // countersStruct &counters, 
                                       const jetsCollectionStruct& jetsCollection, const int& jetIndex, std::vector<angularVariablesStruct> &selectedPhotonAngles, std::vector<angularVariablesStruct> &selectedTruePhotonAngles) {
  jetExaminationResultsStruct results;
  jetProperties& properties = results.jet_properties;

  std::map<jetCriterion, bool> bits;
  bool passesPT_nominal = false;
  bool passesPT_JECDown = false;
  bool passesPT_JECUp = false;

  //Kinematic cuts: eta, pT
  properties[jetProperty::eta] = (jetsCollection.eta)->at(jetIndex);
  properties[jetProperty::phi] = (jetsCollection.phi)->at(jetIndex);
  angularVariablesStruct jetAngle = angularVariablesStruct(properties[jetProperty::eta], properties[jetProperty::phi]);
  float minDeltaR = jetAngle.getMinDeltaR(selectedPhotonAngles);
  properties[jetProperty::deltaR_nearestCaloPhoton] = minDeltaR;
  bits[jetCriterion::deltaR_photon] = ((minDeltaR > parameters.minDeltaRCut) || (minDeltaR < 0.));
  results.isAwayFromCaloPhoton = bits[jetCriterion::deltaR_photon];
  properties[jetProperty::deltaR_nearestTruePhoton] = jetAngle.getMinDeltaR(selectedTruePhotonAngles);

  float absEta = std::fabs(properties[jetProperty::eta]);
  bits[jetCriterion::eta] = (absEta < parameters.jetEtaCut);

  float jet_pT = ((jetsCollection.pT)->at(jetIndex));
  properties[jetProperty::pT] = jet_pT;
  passesPT_nominal = (jet_pT > parameters.jetpTCut);

  float jecFractionalUncertainty = 0.;
  if (options.isMC) {
    jecFractionalUncertainty = (jetsCollection.JECUncertainty)->at(jetIndex);
    float jet_pT_JECDown = (1.0 - jecFractionalUncertainty)*jet_pT;
    float jet_pT_JECUp = (1.0 + jecFractionalUncertainty)*jet_pT;
    passesPT_JECDown = (jet_pT_JECDown > parameters.jetpTCut);
    passesPT_JECUp = (jet_pT_JECUp > parameters.jetpTCut);
  }

  results.prefireWeights = findPrefireWeights((jetsCollection.eta)->at(jetIndex), jet_pT, parameters.prefiringEfficiencyMap);

  // ID cuts: PUID, jetID
  properties[jetProperty::PUID] = (jetsCollection.PUID)->at(jetIndex);
  bits[jetCriterion::puID] = (properties[jetProperty::PUID] > parameters.jetPUIDThreshold);

  properties[jetProperty::jetID] = (jetsCollection.ID)->at(jetIndex);
  bits[jetCriterion::jetID] = (properties[jetProperty::jetID] == 6);

  int nNonPTFalseBits = getNFalseBits(bits);
  bool passesNonPTCriteria = (nNonPTFalseBits == 0);
  results.passesSelectionJECDown = (passesNonPTCriteria && passesPT_JECDown);
  results.passesSelectionJECUp = (passesNonPTCriteria && passesPT_JECUp);

  bits[jetCriterion::pT] = passesPT_nominal;
  assert(static_cast<int>(bits.size()) == static_cast<int>(jetCriterion::nJetCriteria));
  int nFalseBits = getNFalseBits(bits);
  results.passesSelectionJECNominal = (nFalseBits == 0);
  results.isMarginallyUnselected = (nFalseBits == 1);
  if (results.isMarginallyUnselected) results.marginallyUnselectedCriterion = getFirstFalseCriterion(bits);

  results.contributesToHT = false;
  if (nFalseBits == 0) results.contributesToHT = true;
  else if ((nFalseBits == 1) && results.marginallyUnselectedCriterion == jetCriterion::deltaR_photon) results.contributesToHT = true;

  if ((results.isMarginallyUnselected || results.passesSelectionJECNominal) || (results.passesSelectionJECUp || results.passesSelectionJECDown)) assert(static_cast<int>((results.jet_properties).size()) == static_cast<int>(jetProperty::nJetProperties));
  return results;
}

int getMCBinIndex(const float& generated_gluinoMass, const float& generated_neutralinoMass) {
  if ((generated_gluinoMass > 1675.) && (generated_gluinoMass < 1725.)) {
    if ((generated_neutralinoMass > 850.) && (generated_neutralinoMass < 950.)) {
      return 1;
    }
  }
  else if ((generated_gluinoMass > 1275.) && (generated_gluinoMass < 1325.)) {
    if ((generated_neutralinoMass > 1280.) && (generated_neutralinoMass < 1300.)) {
      return 2;
    }
  }
  else if ((generated_gluinoMass > 1275.) && (generated_gluinoMass < 1325.)) {
    if ((generated_neutralinoMass > 100.) && (generated_neutralinoMass < 118.75)) {
      return 3;
    }
  }
  return 0;
}

eventExaminationResultsStruct examineEvent(optionsStruct &options, parametersStruct &parameters, // countersStruct &counters, 
                                           Long64_t& entryIndex, eventDetailsStruct& eventDetails, MCCollectionStruct &MCCollection, photonsCollectionStruct &photonsCollection, jetsCollectionStruct &jetsCollection, // const STRegionsStruct& STRegions, 
                                           statisticsHistograms& statistics) {
  eventExaminationResultsStruct eventResult;

  eventResult.eventIndex = entryIndex;
  selectionRegion& region = eventResult.evt_region;
  std::map<eventSelectionCriterion, bool> selectionBits;
  eventProperties event_properties;
  float& event_ST = eventResult.evt_ST;
  int& n_jetsDR = eventResult.evt_nJetsDR;
  std::map<shiftType, float>& shifted_ST = eventResult.evt_shifted_ST;
  std::map<shiftType, int>& shifted_nJetsDR = eventResult.evt_shifted_nJetsDR;

  selectionBits[eventSelectionCriterion::HLTPhoton] = true;
  if (!(options.isMC) && parameters.HLTPhotonBit >= 0) { // Apply HLT photon selection iff input is not MC and HLTBit is set to a positive integer
    selectionBits[eventSelectionCriterion::HLTPhoton] = checkHLTBit(eventDetails.HLTPhotonBits, parameters.HLTPhotonBit);
  }

  // bool passesNonHLTEventSelection = true;
  // Additional selection, only for MC
  float generated_gluinoMass = 0.;
  float generated_neutralinoMass = 0.;
  selectionBits[eventSelectionCriterion::MCGenInformation] = true;
  truthPhotonPropertiesCollection selectedTruePhotonProperties;
  std::vector<angularVariablesStruct> selectedTruePhotonAngles;
  int nPhotonsWithNeutralinoMom = 0;
  int MCBinIndex = 0;
  if (options.isMC) {
    bool gluinoMassIsSet = false;
    bool neutralinoMassIsSet = false;
    for (int MCIndex = 0; MCIndex < eventDetails.nMCParticles; ++MCIndex) {
      MCExaminationResultsStruct MCExaminationResults = examineMCParticle(parameters, // counters, 
                                                                        // (eventDetails.nMCParticles), 
                                                                          MCCollection, MCIndex);
      if (MCExaminationResults.isPhotonWithNeutralinoMom){
        ++nPhotonsWithNeutralinoMom;
        selectedTruePhotonProperties.push_back(MCExaminationResults.truth_photon_properties);
        selectedTruePhotonAngles.push_back(angularVariablesStruct((MCExaminationResults.truth_photon_properties)[truthPhotonProperty::eta], (MCExaminationResults.truth_photon_properties)[truthPhotonProperty::phi]));
      }
      if ((MCExaminationResults.gluinoMass > 0.) && !(gluinoMassIsSet)) {
        generated_gluinoMass = MCExaminationResults.gluinoMass;
        gluinoMassIsSet = true;
      }
      if ((MCExaminationResults.neutralinoMass > 0.) && !(neutralinoMassIsSet)) {
        generated_neutralinoMass = MCExaminationResults.neutralinoMass;
        neutralinoMassIsSet = true;
      }
    }
    selectionBits[eventSelectionCriterion::MCGenInformation] = (nPhotonsWithNeutralinoMom == 2);
    if (selectionBits[eventSelectionCriterion::MCGenInformation] && (!(gluinoMassIsSet && neutralinoMassIsSet))) {
      std::cout << "ERROR: Unable to find gluino or neutralino mass in an event that passes MC selection." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    MCBinIndex = getMCBinIndex(generated_gluinoMass, generated_neutralinoMass);
  }
  event_properties[eventProperty::MC_nPhotonsWithNeutralinoMom] = nPhotonsWithNeutralinoMom;

  // Photon selection
  std::vector<angularVariablesStruct> selectedPhotonAngles;
  std::vector<TLorentzVector> selectedPhotonFourMomenta;
  photonPropertiesCollection selectedMediumPhotonProperties;
  photonPropertiesCollection selectedMediumPhotonProperties_closeToTruePhoton;
  photonPropertiesCollection selectedMediumPhotonProperties_awayFromTruePhoton;
  photonPropertiesCollection selectedFakePhotonProperties;
  photonPropertiesCollection selectedFakePhotonProperties_closeToTruePhoton;
  photonPropertiesCollection selectedFakePhotonProperties_awayFromTruePhoton;

  unselectedFakePhotonPropertiesCollection unselected_fake_pho_properties;
  unselectedFakePhotonPropertiesCollection unselected_fake_pho_properties_closeToTruePhoton;
  unselectedFakePhotonPropertiesCollection unselected_fake_pho_properties_awayFromTruePhoton;
  unselectedMediumPhotonPropertiesCollection unselected_medium_pho_properties;
  unselectedMediumPhotonPropertiesCollection unselected_medium_pho_properties_closeToTruePhoton;
  unselectedMediumPhotonPropertiesCollection unselected_medium_pho_properties_awayFromTruePhoton;

  int n_selectedPhotonsPassingSubLeadingpTCut = 0;
  int n_selectedPhotonsPassingLeadingpTCut = 0;
  int n_mediumPhotons = 0;
  int n_fakePhotons = 0;
  // int nVetoPhotons = 0;
  for (Int_t photonIndex = 0; photonIndex < (eventDetails.nPhotons); ++photonIndex) {
    photonExaminationResultsStruct photonExaminationResults = examinePhoton(options, parameters, // counters, 
                                                                            (eventDetails.eventRho), photonsCollection, photonIndex, selectedTruePhotonAngles);
    if (photonExaminationResults.isSelectedFake || photonExaminationResults.isSelectedMedium) {
      ++n_selectedPhotonsPassingSubLeadingpTCut;
      float photon_ET = (photonExaminationResults.pho_properties)[photonProperty::pT];
      if (photon_ET > parameters.pTCutLeading) ++n_selectedPhotonsPassingLeadingpTCut;
      event_ST += photon_ET;
      if (options.isMC) {
        for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
          shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
          addShiftedEToSTMap(photon_ET, shifted_ST, typeIndex); // no effect on photon's contribution to ST due to any of the shifts
        }
        (eventResult.evt_photonMCScaleFactors).nominal *= (photonExaminationResults.MCScaleFactors).nominal;
        (eventResult.evt_photonMCScaleFactors).down *= (photonExaminationResults.MCScaleFactors).down;
        (eventResult.evt_photonMCScaleFactors).up *= (photonExaminationResults.MCScaleFactors).up;
      }
      selectedPhotonAngles.push_back(angularVariablesStruct((photonExaminationResults.pho_properties)[photonProperty::eta], (photonExaminationResults.pho_properties)[photonProperty::phi]));
      TLorentzVector photonFourMomentum;
      photonFourMomentum.SetPtEtaPhiE(photon_ET, (photonExaminationResults.pho_properties)[photonProperty::eta], (photonExaminationResults.pho_properties)[photonProperty::phi], photonExaminationResults.energy);
      selectedPhotonFourMomenta.push_back(photonFourMomentum);
    }
    if (photonExaminationResults.isSelectedMedium) {
      ++n_mediumPhotons;
      selectedMediumPhotonProperties.push_back(photonExaminationResults.pho_properties);
      if (options.isMC) {
        float nearestTruePhotonDeltaR = (photonExaminationResults.pho_properties)[photonProperty::deltaR_nearestTruePhoton];
        if (nearestTruePhotonDeltaR >= parameters.minDeltaRCut) selectedMediumPhotonProperties_awayFromTruePhoton.push_back(photonExaminationResults.pho_properties);
        else if (nearestTruePhotonDeltaR > 0.) selectedMediumPhotonProperties_closeToTruePhoton.push_back(photonExaminationResults.pho_properties);
      }
    }
    else if (photonExaminationResults.isSelectedFake) {
      ++n_fakePhotons;
      selectedFakePhotonProperties.push_back(photonExaminationResults.pho_properties);
      if (options.isMC) {
        float nearestTruePhotonDeltaR = (photonExaminationResults.pho_properties)[photonProperty::deltaR_nearestTruePhoton];
        if (nearestTruePhotonDeltaR >= parameters.minDeltaRCut) selectedFakePhotonProperties_awayFromTruePhoton.push_back(photonExaminationResults.pho_properties);
        else if (nearestTruePhotonDeltaR > 0.) selectedFakePhotonProperties_closeToTruePhoton.push_back(photonExaminationResults.pho_properties);
      }
    }
    if (photonExaminationResults.isMarginallyUnselectedFake) {
      unselected_fake_pho_properties.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedFakeCriterion, photonExaminationResults.pho_properties));
      if (options.isMC) {
        float nearestTruePhotonDeltaR = (photonExaminationResults.pho_properties)[photonProperty::deltaR_nearestTruePhoton];
        if (nearestTruePhotonDeltaR >= parameters.minDeltaRCut) unselected_fake_pho_properties_awayFromTruePhoton.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedFakeCriterion, photonExaminationResults.pho_properties));
        else if (nearestTruePhotonDeltaR > 0.) unselected_fake_pho_properties_closeToTruePhoton.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedFakeCriterion, photonExaminationResults.pho_properties));
      }
    }
    if (photonExaminationResults.isMarginallyUnselectedMedium) {
      unselected_medium_pho_properties.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedMediumCriterion, photonExaminationResults.pho_properties));
      if (options.isMC) {
        float nearestTruePhotonDeltaR = (photonExaminationResults.pho_properties)[photonProperty::deltaR_nearestTruePhoton];
        if (nearestTruePhotonDeltaR >= parameters.minDeltaRCut) unselected_medium_pho_properties_awayFromTruePhoton.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedMediumCriterion, photonExaminationResults.pho_properties));
        else if (nearestTruePhotonDeltaR > 0.) unselected_medium_pho_properties_closeToTruePhoton.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedMediumCriterion, photonExaminationResults.pho_properties));
      }
    }
  }

  selectionBits[eventSelectionCriterion::photonEnergy] = ((n_selectedPhotonsPassingSubLeadingpTCut >= 2) && (n_selectedPhotonsPassingLeadingpTCut >= 1));

  selectionBits[eventSelectionCriterion::photonQuality] = false;
  if (n_mediumPhotons == 2) {
    selectionBits[eventSelectionCriterion::photonQuality] = true;
    region = selectionRegion::signal;
  }
  else if ((n_mediumPhotons == 1) && (n_fakePhotons >= 1)) {
    selectionBits[eventSelectionCriterion::photonQuality] = true;
    region = selectionRegion::control_mediumfake;
  }
  else if ((n_mediumPhotons == 0) && (n_fakePhotons >= 2)) {
    selectionBits[eventSelectionCriterion::photonQuality] = true;
    region = selectionRegion::control_fakefake;
  }

  float evt_invariantMass = -1.0;
  selectionBits[eventSelectionCriterion::invariantMass] = true;
  if ((n_mediumPhotons + n_fakePhotons) >= 2) {
    evt_invariantMass = getDiphotonInvariantMass(selectedPhotonFourMomenta);
    selectionBits[eventSelectionCriterion::invariantMass] = (evt_invariantMass >= parameters.invariantMassCut);
  }

  event_properties[eventProperty::nMediumPhotons] = 1.0*n_mediumPhotons;
  event_properties[eventProperty::nFakePhotons] = 1.0*n_fakePhotons;
  event_properties[eventProperty::nSelectedPhotonsPassingSubLeadingpTCut] = n_selectedPhotonsPassingSubLeadingpTCut;
  event_properties[eventProperty::nSelectedPhotonsPassingLeadingpTCut] = n_selectedPhotonsPassingLeadingpTCut;
  event_properties[eventProperty::invariantMass] = evt_invariantMass;

  // bool eventPassesNonHLTPhotonSelection = passesNonHLTEventSelection;

  // Jet selection
  float evt_hT = 0;
  jetPropertiesCollection selectedJetProperties;
  jetPropertiesCollection selectedJetProperties_closeToTruePhoton;
  jetPropertiesCollection selectedJetProperties_awayFromTruePhoton;
  unselectedJetPropertiesCollection unselected_jet_properties;
  unselectedJetPropertiesCollection unselected_jet_properties_closeToTruePhoton;
  unselectedJetPropertiesCollection unselected_jet_properties_awayFromTruePhoton;

  for (Int_t jetIndex = 0; jetIndex < (eventDetails.nJets); ++jetIndex) {
    jetExaminationResultsStruct jetExaminationResults = examineJet(options, parameters, // counters, 
                                                                   jetsCollection, jetIndex, selectedPhotonAngles, selectedTruePhotonAngles);
    // if (options.isMC) counters.jetTotalCountersMCMap->Fill(generated_gluinoMass, generated_neutralinoMass);
    (eventResult.evt_prefireWeights).nominal *= (jetExaminationResults.prefireWeights).nominal; // All jets, whether or not they pass any of the cuts, contribute to the prefiring weight
    (eventResult.evt_prefireWeights).down *= (jetExaminationResults.prefireWeights).down;
    (eventResult.evt_prefireWeights).up *= (jetExaminationResults.prefireWeights).up;
    if (jetExaminationResults.passesSelectionJECNominal) selectedJetProperties.push_back(jetExaminationResults.jet_properties);
    else if (jetExaminationResults.isMarginallyUnselected) {
      unselected_jet_properties.push_back(std::make_pair(jetExaminationResults.marginallyUnselectedCriterion, jetExaminationResults.jet_properties));
      if (options.isMC) {
        float nearestTruePhotonDeltaR = (jetExaminationResults.jet_properties)[jetProperty::deltaR_nearestTruePhoton];
        if (nearestTruePhotonDeltaR >= parameters.minDeltaRCut) unselected_jet_properties_awayFromTruePhoton.push_back(std::make_pair(jetExaminationResults.marginallyUnselectedCriterion, jetExaminationResults.jet_properties));
        else if (nearestTruePhotonDeltaR > 0.) unselected_jet_properties_closeToTruePhoton.push_back(std::make_pair(jetExaminationResults.marginallyUnselectedCriterion, jetExaminationResults.jet_properties));
      }
    }
    if (jetExaminationResults.contributesToHT) {
      evt_hT += jetExaminationResults.jet_properties[jetProperty::pT]; // Add to hT whether or not jet passes deltaR check
      if (jetExaminationResults.passesSelectionJECNominal) {
        event_ST += jetExaminationResults.jet_properties[jetProperty::pT]; // Add to sT only if jet passes deltaR check, to avoid double-counting
        ++n_jetsDR; // Count only those jets that are sufficiently away from a photon
        selectedJetProperties.push_back(jetExaminationResults.jet_properties);
        if (options.isMC) {
          float nearestTruePhotonDeltaR = (jetExaminationResults.jet_properties)[jetProperty::deltaR_nearestTruePhoton];
          if (nearestTruePhotonDeltaR >= parameters.minDeltaRCut) selectedJetProperties_awayFromTruePhoton.push_back(jetExaminationResults.jet_properties);
          else if (nearestTruePhotonDeltaR > 0.) selectedJetProperties_closeToTruePhoton.push_back(jetExaminationResults.jet_properties);
        }
      }
    }
    if (options.isMC && ((jetExaminationResults.passesSelectionJECDown || jetExaminationResults.passesSelectionJECDown) || jetExaminationResults.passesSelectionJECNominal)) { // Actually we just need to check JECDown
      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
        shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);

        bool passes_selection = jetExaminationResults.passesSelectionJECNominal;
        float shifted_contribution = (jetExaminationResults.jet_properties[jetProperty::pT]);
        if (typeIndex == shiftType::JECDown) {
          passes_selection = jetExaminationResults.passesSelectionJECDown;
          shifted_contribution = (1.0 - (jetExaminationResults.jecFractionalUncertainty))*(jetExaminationResults.jet_properties[jetProperty::pT]);
        }
        else if (typeIndex == shiftType::JECUp) {
          passes_selection = jetExaminationResults.passesSelectionJECUp;
          shifted_contribution = (1.0 + (jetExaminationResults.jecFractionalUncertainty))*(jetExaminationResults.jet_properties[jetProperty::pT]);
        }

        if (passes_selection) {
          addShiftedEToSTMap(shifted_contribution, shifted_ST, typeIndex);
          incrementNJetsMap(shifted_nJetsDR, typeIndex);
        }
      }
    }
  }
  event_properties[eventProperty::hT] = evt_hT;
  event_properties[eventProperty::nJetsDR] = n_jetsDR;
  int max_nJets = n_jetsDR;
  if (options.isMC) { // this makes sure that the nJets used to make the decision whether or not to save the event is the maximum nJets accounting for all the shifts
    int maxNJetsShifted = getMaxNJets(shifted_nJetsDR);
    if (maxNJetsShifted > max_nJets) max_nJets = maxNJetsShifted;
  }

  selectionBits[eventSelectionCriterion::NJets] = (max_nJets >= 2);
  // Add MET to ST
  event_ST += eventDetails.PFMET;
  event_properties[eventProperty::ST] = event_ST;

  if (options.isMC) {
    // Add shifted energies
    addShiftedEToSTMap(eventDetails.PFMET, shifted_ST, shiftType::JECDown);
    addShiftedEToSTMap(eventDetails.PFMET, shifted_ST, shiftType::JECUp);
    addShiftedEToSTMap(eventDetails.PFMET_UnclusteredDown, shifted_ST, shiftType::UnclusteredMETDown);
    addShiftedEToSTMap(eventDetails.PFMET_UnclusteredUp, shifted_ST, shiftType::UnclusteredMETUp);
    addShiftedEToSTMap(eventDetails.PFMET_JERDown, shifted_ST, shiftType::JERMETDown);
    addShiftedEToSTMap(eventDetails.PFMET_JERUp, shifted_ST, shiftType::JERMETUp);

    // // Fill acceptance maps
    // if ((n_jetsDR >= 2) && (event_ST >= STRegions.STNormRangeMin)) {
    //   int STRegionIndex = (STRegions.STAxis).FindFixBin(event_ST);
    //   // if (passesMCTruth) {
    //   //   counters.acceptanceMCMap_eventPassesTruth[STRegionIndex]->Fill(generated_gluinoMass, generated_neutralinoMass);
    //   //   if (eventPassesNonHLTPhotonSelection) {
    //   //     counters.acceptanceMCMap_eventPassesSelection[STRegionIndex]->Fill(generated_gluinoMass, generated_neutralinoMass);
    //   //   }
    //   // }
    // }
  }
  // else if ((parameters.HLTPhotonBit >= 0) && passesNonHLTEventSelection) {
  //   if ((n_jetsDR >= 2) && (event_ST >= STRegions.STNormRangeMin)) {
  //     int nJetsBin = (n_jetsDR > 6 ? 6 : n_jetsDR);
  //     int STRegionIndex = (STRegions.STAxis).FindFixBin(event_ST);
  //     (counters.nTriggeredEvents_cuts)->Fill(1.0*STRegionIndex, 1.0*nJetsBin);
  //     if (selectionBits[eventSelectionCriterion::HLTPhoton]) {
  //       (counters.nTriggeredEvents_cutsANDtrigger)->Fill(1.0*STRegionIndex, 1.0*nJetsBin);
  //     }
  //   }
  // }

  assert(static_cast<int>(selectionBits.size()) == static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria));
  int nEventFalseBits = getNFalseBits(selectionBits);
  eventProperties temp1 = initialize_eventProperties_with_defaults(); // temp1 and temp2 are dummies -- they won't contribute to the histograms
  if (nEventFalseBits == 0) {
    unselectedEventProperties temp2 = std::make_pair(eventSelectionCriterion::nEventSelectionCriteria, temp1);
    statistics.fillStatisticsHistograms(event_properties, false, temp2,
                             selectedTruePhotonProperties,
                             selectedMediumPhotonProperties,
                             selectedMediumPhotonProperties_closeToTruePhoton,
                             selectedMediumPhotonProperties_awayFromTruePhoton,
                             unselected_medium_pho_properties,
                             unselected_medium_pho_properties_closeToTruePhoton,
                             unselected_medium_pho_properties_awayFromTruePhoton,
                             selectedFakePhotonProperties,
                             selectedFakePhotonProperties_closeToTruePhoton,
                             selectedFakePhotonProperties_awayFromTruePhoton,
                             unselected_fake_pho_properties,
                             unselected_fake_pho_properties_closeToTruePhoton,
                             unselected_fake_pho_properties_awayFromTruePhoton,
                             selectedJetProperties,
                             selectedJetProperties_closeToTruePhoton,
                             selectedJetProperties_awayFromTruePhoton,
                             unselected_jet_properties,
                             unselected_jet_properties_closeToTruePhoton,
                             unselected_jet_properties_awayFromTruePhoton,
                             region, options.isMC, MCBinIndex);
  }
  else if (nEventFalseBits == 1) {
    eventSelectionCriterion marginallyUnselectedEventCriterion = getFirstFalseCriterion(selectionBits);
    unselectedEventProperties unselected_event_properties = std::make_pair(marginallyUnselectedEventCriterion, event_properties);
    statistics.fillStatisticsHistograms(temp1, true, unselected_event_properties,
                             selectedTruePhotonProperties,
                             selectedMediumPhotonProperties,
                             selectedMediumPhotonProperties_closeToTruePhoton,
                             selectedMediumPhotonProperties_awayFromTruePhoton,
                             unselected_medium_pho_properties,
                             unselected_medium_pho_properties_closeToTruePhoton,
                             unselected_medium_pho_properties_awayFromTruePhoton,
                             selectedFakePhotonProperties,
                             selectedFakePhotonProperties_closeToTruePhoton,
                             selectedFakePhotonProperties_awayFromTruePhoton,
                             unselected_fake_pho_properties,
                             unselected_fake_pho_properties_closeToTruePhoton,
                             unselected_fake_pho_properties_awayFromTruePhoton,
                             selectedJetProperties,
                             selectedJetProperties_closeToTruePhoton,
                             selectedJetProperties_awayFromTruePhoton,
                             unselected_jet_properties,
                             unselected_jet_properties_closeToTruePhoton,
                             unselected_jet_properties_awayFromTruePhoton,
                             region, options.isMC, MCBinIndex);
  }

  if (nEventFalseBits <= 1) assert(static_cast<int>(event_properties.size()) == static_cast<int>(eventProperty::nEventProperties));
  return eventResult;
}

void loopOverEvents(optionsStruct &options, parametersStruct &parameters, // countersStruct &counters, const STRegionsStruct& STRegions,
                    std::vector<eventExaminationResultsStruct>& selectedEventsInfo, statisticsHistograms& statistics) {
  // std::vector<eventExaminationResultsStruct> selectedEventsInfo;

  std::ifstream fileWithInputFilesList(options.inputFilesList);
  if (!fileWithInputFilesList.is_open()) {
    std::cout << "ERROR: Failed to open file with path: " << options.inputFilesList << std::endl;
    std::exit(EXIT_FAILURE);
  }
  TChain inputChain("ggNtuplizer/EventTree");
  std::cout << "Starting to add files to chain..." << std::endl;
  while (!fileWithInputFilesList.eof()) {
    std::string inputFileName;
    fileWithInputFilesList >> inputFileName;
    if (!inputFileName.empty()) {
      // std::cout << "Adding... " << inputFileName << std::endl;
      inputChain.Add(inputFileName.c_str()); // Add files to TChain
    }
  }
  std::cout << "Finished adding files to chain!" << std::endl;
  fileWithInputFilesList.close();

  Long64_t nEvts = inputChain.GetEntries();
  Long64_t nEntriesToProcess = 1 + options.counterEndInclusive - options.counterStartInclusive;
  std::cout << "Number of available events: " << nEvts << std::endl;

  inputChain.SetBranchStatus("*", 0); // so that only the needed branches, explicitly activated below, are read in per event

  eventDetailsStruct eventDetails = eventDetailsStruct(inputChain, options.isMC);
  photonsCollectionStruct photonsCollection = photonsCollectionStruct(inputChain);
  jetsCollectionStruct jetsCollection = jetsCollectionStruct(inputChain, options.isMC);
  MCCollectionStruct MCCollection = MCCollectionStruct(inputChain, options.isMC);

  tmProgressBar progressBar = tmProgressBar(static_cast<int>(nEntriesToProcess));
  int progressBarUpdatePeriod = ((nEntriesToProcess < 1000) ? 1 : static_cast<int>(0.5 + 1.0*(nEntriesToProcess/1000)));
  progressBar.initialize();
  for (Long64_t entryIndex = options.counterStartInclusive; entryIndex <= options.counterEndInclusive; ++entryIndex) {
    if (entryIndex > nEvts) {
      std::cout << "ERROR: Entry index falls outside event range. Index: " << entryIndex << ", available number of events: " << nEvts << std::endl;
      std::exit(EXIT_FAILURE);
    }
    Long64_t loadStatus = inputChain.LoadTree(entryIndex);
    if (loadStatus < 0) {
      std::cout << "ERROR in loading tree for entry index: " << entryIndex << "; load status = " << loadStatus << std::endl;
      std::exit(EXIT_FAILURE);
    }
    int nBytesRead = inputChain.GetEntry(entryIndex, 0); // Get only the required branches
    if (nBytesRead <= 0) {
      std::cout << "ERROR: Failed to read SOME information from entry at index: " << entryIndex << "; nBytesRead = " << nBytesRead << std::endl;
      std::exit(EXIT_FAILURE);
    }

    int entryProcessing = static_cast<int>(entryIndex - options.counterStartInclusive);
    if (entryProcessing > 0 && ((static_cast<int>(entryProcessing) % progressBarUpdatePeriod == 0) || entryProcessing == static_cast<int>(nEntriesToProcess-1))) progressBar.updateBar(static_cast<double>(1.0*entryProcessing/nEntriesToProcess), entryProcessing);

    eventExaminationResultsStruct eventExaminationResults = examineEvent(options, parameters, // counters, 
                                                                         entryIndex, eventDetails, MCCollection, photonsCollection, jetsCollection, // STRegions, 
                                                                         statistics);
    bool passesEventSelection = (eventExaminationResults.evt_region != selectionRegion::nSelectionRegions);
    // incrementCounters(miscCounter::totalEvents, counters, false, 0., 0.);
    if (!(passesEventSelection)) {
      // incrementCounters(miscCounter::failingEvents, counters, false, 0., 0.);
      continue;
    }
    // incrementCounters(miscCounter::acceptedEvents, counters, false, 0., 0.);
    selectedEventsInfo.push_back(eventExaminationResults);
  }
  progressBar.terminate();
}

void writeSelectionToFile(optionsStruct &options, TFile *outputFile, const std::vector<eventExaminationResultsStruct>& selectedEventsInfo, selectionRegion& region) {
  std::string regionName = selectionRegionNames[region];
  std::cout << "Beginning to write selected events to file for selection type: " <<  regionName << std::endl;
  std::ifstream fileWithInputFilesList(options.inputFilesList);
  if (!fileWithInputFilesList.is_open()) {
    std::cout << "ERROR: Failed to open file with path: " << options.inputFilesList << std::endl;
    std::exit(EXIT_FAILURE);
  }
  TChain inputChain("ggNtuplizer/EventTree");
  std::cout << "Starting to add files to chain..." << std::endl;
  while (!fileWithInputFilesList.eof()) {
    std::string inputFileName;
    fileWithInputFilesList >> inputFileName;
    if (!inputFileName.empty()) {
      // std::cout << "Adding... " << inputFileName << std::endl;
      inputChain.Add(inputFileName.c_str()); // Add files to TChain
    }
  }
  std::cout << "Finished adding files to chain!" << std::endl;
  fileWithInputFilesList.close();

  TDirectory *outDir = outputFile->mkdir("ggNtuplizer");
  outDir->cd();
  TTree *outputTree = inputChain.CloneTree(0);

  int nJetsDR; // stores number of jets in event passing deltaR cut
  outputTree->Branch("b_nJets", &nJetsDR, "b_nJets/I");
  float ST; // stores event sT
  outputTree->Branch("b_evtST", &ST, "b_evtST/F");
  eventWeightsStruct prefireWeights = eventWeightsStruct(1.0f, 1.0f, 1.0f); // stores prefiring weights and errors
  outputTree->Branch("b_evtPrefiringWeight", &(prefireWeights.nominal), "b_evtPrefiringWeight/F");
  outputTree->Branch("b_evtPrefiringWeightDown", &(prefireWeights.down), "b_evtPrefiringWeightDown/F");
  outputTree->Branch("b_evtPrefiringWeightUp", &(prefireWeights.up), "b_evtPrefiringWeightUp/F");
  eventWeightsStruct photonMCScaleFactors = eventWeightsStruct(1.0f, 1.0f, 1.0f); // stores scale factors for MC samples

  std::map<shiftType, float> shifted_ST;
  std::map<shiftType, int> shifted_nJetsDR;
  if (options.isMC) {
    std::string branchPrefix_ST = "b_evtST_shifted_";
    std::string branchPrefix_nJets = "b_nJets_shifted_";
    for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
      shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
      shifted_ST[typeIndex] = 0.;
      outputTree->Branch((getShiftedVariableBranchName(typeIndex, "evtST")).c_str(), &(shifted_ST[typeIndex]), (getShiftedVariableBranchName(typeIndex, "evtST") + "/F").c_str());
      shifted_nJetsDR[typeIndex] = 0;
      outputTree->Branch((getShiftedVariableBranchName(typeIndex, "nJets")).c_str(), &(shifted_nJetsDR[typeIndex]), (getShiftedVariableBranchName(typeIndex, "nJets") + "/I").c_str());
    }
    outputTree->Branch("b_evtphotonMCScaleFactor", &(photonMCScaleFactors.nominal), "b_evtphotonMCScaleFactor/F");
    outputTree->Branch("b_evtphotonMCScaleFactorDown", &(photonMCScaleFactors.down), "b_evtphotonMCScaleFactorDown/F");
    outputTree->Branch("b_evtphotonMCScaleFactorUp", &(photonMCScaleFactors.up), "b_evtphotonMCScaleFactorUp/F");
  }

  int nSelectedEvents = static_cast<int>(0.5 + selectedEventsInfo.size());
  tmProgressBar progressBar = tmProgressBar(nSelectedEvents);
  int progressBarUpdatePeriod = ((nSelectedEvents < 1000) ? 1 : static_cast<int>(0.5 + 1.0*(nSelectedEvents/1000)));
  int processingIndex = -1;
  progressBar.initialize();

  for (auto& selectedEventInfo : selectedEventsInfo) {
    ++processingIndex;

    Long64_t index = selectedEventInfo.eventIndex;
    nJetsDR = selectedEventInfo.evt_nJetsDR;
    ST = selectedEventInfo.evt_ST;
    prefireWeights = eventWeightsStruct((selectedEventInfo.evt_prefireWeights).nominal, (selectedEventInfo.evt_prefireWeights).down, (selectedEventInfo.evt_prefireWeights).up);
    if (options.isMC) {
      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
        shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
        shifted_ST[typeIndex] = (selectedEventInfo.evt_shifted_ST).at(typeIndex);
        shifted_nJetsDR[typeIndex] = (selectedEventInfo.evt_shifted_nJetsDR).at(typeIndex);
      }
      photonMCScaleFactors = eventWeightsStruct((selectedEventInfo.evt_photonMCScaleFactors).nominal, (selectedEventInfo.evt_photonMCScaleFactors).down, (selectedEventInfo.evt_photonMCScaleFactors).up);
    }

    Long64_t loadStatus = inputChain.LoadTree(index);
    if (loadStatus < 0) {
      std::cout << "ERROR in loading tree for entry index: " << index << "; load status = " << loadStatus << std::endl;
      std::exit(EXIT_FAILURE);
    }
    int nBytesRead = inputChain.GetEntry(index, 1); // Get all branches before filling output tree
    if (nBytesRead <= 0) {
      std::cout << "ERROR: For selected event, failed to read ALL information from entry at index: " << index << "; nBytesRead = " << nBytesRead << std::endl;
      std::exit(EXIT_FAILURE);
    }
    if (processingIndex > 0 && ((processingIndex%progressBarUpdatePeriod == 0) || processingIndex == static_cast<int>(nSelectedEvents-1))) progressBar.updateBar(static_cast<double>(1.0*processingIndex/nSelectedEvents), processingIndex);
    if (selectedEventInfo.evt_region != region) continue;

    outputTree->Fill();
  }
  progressBar.terminate();
  outputFile->Write();
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  do_sanity_checks_objectProperties();
  do_sanity_checks_selectionCriteria();
  tmArgumentParser argumentParser = tmArgumentParser("Run the event selection.");
  argumentParser.addArgument("inputFilesList", "", true, "Path to file containing list of input files.");
  argumentParser.addArgument("isMC", "false", false, "Input file is a MC sample -- disable HLT photon trigger and enable additional MC selection.");
  argumentParser.addArgument("counterStartInclusive", "", true, "Event number from input file from which to start. The event with this index is included in the processing.");
  argumentParser.addArgument("counterEndInclusive", "", true, "Event number from input file at which to end. The event with this index is included in the processing.");
  // argumentParser.addArgument("photonSelectionType", "fake", true, "Photon selection type: can be any one of: \"fake\", \"medium\", \"mediumfake\", or \"singlemedium\".");
  argumentParser.addArgument("year", "2017", false, "Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.");
  // argumentParser.addArgument("inputFile_STRegionBoundaries", "STRegionBoundaries.dat", false, "Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity."); // for trigger efficiency studies
  // all remaining arguments are only used in MC samples to construct the histograms that help in diagnosing efficiency issues.
  // argumentParser.addArgument("nGluinoMassBins", "20", false, "nBins on the gluino mass axis"); // (800 - 25) GeV --> (1750 + 25) GeV in steps of 50 GeV
  // argumentParser.addArgument("minGluinoMass", "775.0", false, "Min gluino mass for the 2D plots.");
  // argumentParser.addArgument("maxGluinoMass", "1775.0", false, "Max gluino mass for the 2D plots.");
  // argumentParser.addArgument("nNeutralinoMassBins", "133", false, "nBins on the neutralino mass axis.");
  // argumentParser.addArgument("minNeutralinoMass", "93.75", false, "Min neutralino mass for the 2D plots.");
  // argumentParser.addArgument("maxNeutralinoMass", "1756.25", false, "Max neutralino mass for the 2D plots."); // (100 - 6.25) GeV --> (1750 + 6.25) GeV in steps of 12.5 GeV
  argumentParser.setPassedStringValues(argc, argv);

  optionsStruct options = getOptionsFromParser(argumentParser);

  // STRegionsStruct STRegions(options.inputFile_STRegionBoundaries);

  parametersStruct parameters = parametersStruct();
  parameters.tuneParametersForYear(options.year, options.isMC);

  // countersStruct counters = countersStruct();
  // initializeCounters(counters, options, STRegions.nSTSignalBins);

  std::stringstream optionsStringstream;
  optionsStringstream << options;
  // TNamed *optionsObject = new TNamed("optionsString", optionsStringstream.str().c_str());
  std::stringstream parametersStringstream;
  parametersStringstream << parameters;
  // TNamed *parametersObject = new TNamed("parametersString", parametersStringstream.str().c_str());

  std::vector<eventExaminationResultsStruct> selectedEventsInfo;

  statisticsHistograms statistics = statisticsHistograms(options.isMC, parameters.MCBinNames);

  loopOverEvents(options, parameters, // counters, STRegions,
                 selectedEventsInfo, statistics);

  statistics.writeToFile("statisticsHistograms.root");

  for (int selectionRegionIndex = shiftTypeFirst; selectionRegionIndex != static_cast<int>(selectionRegion::nSelectionRegions); ++selectionRegionIndex) {
    selectionRegion region = static_cast<selectionRegion>(selectionRegionIndex);
    std::string outputFilePath = std::string("selection_") + selectionRegionNames[region] + std::string(".root");
    TFile *outputFile = TFile::Open(outputFilePath.c_str(), "RECREATE");
    if (!(outputFile->IsOpen()) || outputFile->IsZombie()) {
      std::cout << "ERROR: Unable to open output file to write. Attempted to create file with path: " << outputFilePath << std::endl;
    }
    writeSelectionToFile(options, outputFile, selectedEventsInfo, region);
    outputFile->Close();
  }

  std::cout << getNDashes(100) << std::endl;

  std::cout << getNDashes(100) << std::endl
            << "Options:" << std::endl
            << optionsStringstream.str() << std::endl
            << getNDashes(100) << std::endl;
  std::cout << getNDashes(100) << std::endl
            << "Parameters:" << std::endl
            << parametersStringstream.str() << std::endl
            << getNDashes(100) << std::endl;
  std::cout << "All done!" << std::endl;
  return 0;
}
