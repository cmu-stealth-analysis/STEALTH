#include "../include/runEventSelection.h"

bool checkHLTBit(const ULong64_t& inputHLTBits, const int& indexOfBitToCheck) {
  return (((inputHLTBits>>indexOfBitToCheck)&1) == 1);
}

float getRhoCorrectedIsolation(const float& uncorrectedIsolation, const PFTypesForEA& PFType, const float& absEta, const float& rho, const EAValuesStruct EAValues[6]) {
  float effectiveArea = 0.;
  unsigned int areaCounter = 0;
  while (areaCounter < 6) {
    if (absEta <= (EAValues[areaCounter]).regionUpperBound) {
      effectiveArea = (EAValues[areaCounter]).getEffectiveArea(PFType);
      break;
    }
    ++areaCounter;
  }
  if (areaCounter >= 6) {
    std::cout << "ERROR: Area counter value = " << areaCounter << ";getRhoCorrectedIsolation called for eta = " << absEta << ", which is above all available upper bounds." << std::endl;
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

photonExaminationResultsStruct examinePhoton(optionsStruct &options, parametersStruct &parameters, countersStruct &counters, const float& rho, const photonsCollectionStruct& photonsCollection, const int& photonIndex, const float& generated_gluinoMass, const float& generated_neutralinoMass) {
  bool isBarrelMedium = false;
  bool isBarrelFake = false;
  bool isEndcapMedium = false;
  bool isEndcapFake = false;
  bool passesCommonCuts = true;
  eventWeightsStruct MCScaleFactors = eventWeightsStruct(1.0f, 1.0f, 1.0f);

  // Kinematic cuts
  float absEta = std::fabs((photonsCollection.eta)->at(photonIndex));
  bool isBarrelPhoton = (absEta < parameters.photonBarrelEtaCut);
  bool isEndcapPhoton = ((absEta > parameters.photonEndcapEtaLow) && (absEta < parameters.photonEndcapEtaHigh));
  applyCondition(counters, photonFailureCategory::eta, passesCommonCuts, (isBarrelPhoton || isEndcapPhoton), options.isMC, generated_gluinoMass, generated_neutralinoMass);
  float photon_ET = std::fabs((photonsCollection.pT)->at(photonIndex));
  bool passesLeadingpTCut = (photon_ET > parameters.pTCutLeading);
  applyCondition(counters, photonFailureCategory::pT, passesCommonCuts, (photon_ET > parameters.pTCutSubLeading), options.isMC, generated_gluinoMass, generated_neutralinoMass);

  if (!(isBarrelPhoton || isEndcapPhoton)) {
    photonExaminationResultsStruct photonExaminationResults = photonExaminationResultsStruct(false, false, false, false, passesLeadingpTCut, ((photonsCollection.eta)->at(photonIndex)), ((photonsCollection.phi)->at(photonIndex)), ((photonsCollection.pT)->at(photonIndex)), ((photonsCollection.energy)->at(photonIndex)), MCScaleFactors);
    return photonExaminationResults;
  }

  // Electron veto
  applyCondition(counters, photonFailureCategory::conversionSafeElectronVeto, passesCommonCuts, (((photonsCollection.electronVeto)->at(photonIndex)) == constants::TRUETOINTT), options.isMC, generated_gluinoMass, generated_neutralinoMass);

  // Quality cuts
  photonQualityCutsStruct* qualityCuts = &(parameters.photonQualityCutsBarrel);
  if (isEndcapPhoton) qualityCuts = &(parameters.photonQualityCutsEndcap);

  applyCondition(counters, photonFailureCategory::hOverE, passesCommonCuts, (((photonsCollection.HOverE)->at(photonIndex)) < qualityCuts->towerHOverE), options.isMC, generated_gluinoMass, generated_neutralinoMass); // H-over-E criterion: same for medium and fake selection

  float pTDependentNeutralIsolationCut = (qualityCuts->neutralIsolation).getPolynomialValue(photon_ET);
  float rhoCorrectedNeutralIsolation = getRhoCorrectedIsolation(((photonsCollection.PFNeutralIsolationUncorrected)->at(photonIndex)), PFTypesForEA::neutralHadron, absEta, rho, parameters.effectiveAreas);
  applyCondition(counters, photonFailureCategory::neutralIsolation, passesCommonCuts, (rhoCorrectedNeutralIsolation < pTDependentNeutralIsolationCut), options.isMC, generated_gluinoMass, generated_neutralinoMass); // neutral isolation criterion: same for medium and fake selection

  float pTDependentPhotonIsolationCut = (qualityCuts->photonIsolation).getPolynomialValue(photon_ET);
  float rhoCorrectedPhotonIsolation = getRhoCorrectedIsolation(((photonsCollection.PFPhotonIsolationUncorrected)->at(photonIndex)), PFTypesForEA::photon, absEta, rho, parameters.effectiveAreas);
  applyCondition(counters, photonFailureCategory::photonIsolation, passesCommonCuts, (rhoCorrectedPhotonIsolation < pTDependentPhotonIsolationCut), options.isMC, generated_gluinoMass, generated_neutralinoMass); // photon isolation criterion: same for medium and fake selection

  float photon_sigmaIEtaIEta = ((photonsCollection.sigmaIEtaIEta)->at(photonIndex));
  bool passesMedium_sigmaIEtaIEtaCut = (photon_sigmaIEtaIEta < qualityCuts->sigmaIEtaIEta);
  bool passesLoose_sigmaIEtaIEtaCut = (photon_sigmaIEtaIEta < qualityCuts->sigmaIEtaIEtaLoose);

  float photon_chargedIsolation = getRhoCorrectedIsolation(((photonsCollection.PFChargedIsolationUncorrected)->at(photonIndex)), PFTypesForEA::chargedHadron, absEta, rho, parameters.effectiveAreas);
  bool passesMedium_chargedIsolationCut = (photon_chargedIsolation < qualityCuts->chargedIsolation);
  bool passesLoose_chargedIsolationCut = (photon_chargedIsolation < qualityCuts->chargedIsolationLoose);

  bool passesMedium_sigmaIEtaIEtaANDChargedIsolationCuts = (passesMedium_sigmaIEtaIEtaCut && passesMedium_chargedIsolationCut);
  bool passesFake_sigmaIEtaIEtaANDChargedIsolationCuts = ((!(passesMedium_sigmaIEtaIEtaANDChargedIsolationCuts)) && (passesLoose_sigmaIEtaIEtaCut && passesLoose_chargedIsolationCut)); // if either sigma-ieta-ieta or charged isolation fail the medium cut, then check if both pass the loose cuts

  bool passesSelectionAsMedium = passesCommonCuts;
  applyCondition(counters, photonFailureCategory::sigmaietaiataANDchargedIsolation, passesSelectionAsMedium, passesMedium_sigmaIEtaIEtaANDChargedIsolationCuts, options.isMC, generated_gluinoMass, generated_neutralinoMass);
  bool passesSelectionAsFake = passesCommonCuts;
  applyCondition(counters, photonFailureCategory::sigmaietaiataANDchargedIsolationLoose, passesSelectionAsFake, passesFake_sigmaIEtaIEtaANDChargedIsolationCuts, options.isMC, generated_gluinoMass, generated_neutralinoMass);

  if (passesSelectionAsFake || passesSelectionAsMedium) incrementCounters(miscCounter::passingPhotons, counters, options.isMC, generated_gluinoMass, generated_neutralinoMass);
  else incrementCounters(miscCounter::failingPhotons, counters, options.isMC, generated_gluinoMass, generated_neutralinoMass);
  incrementCounters(miscCounter::totalPhotons, counters, options.isMC, generated_gluinoMass, generated_neutralinoMass);

  if (options.isMC && (passesSelectionAsFake || passesSelectionAsMedium)) {
    MCScaleFactors = findMCScaleFactors(((photonsCollection.eta)->at(photonIndex)), ((photonsCollection.pT)->at(photonIndex)), parameters.photonMCScaleFactorsMap);
  }

  if (isBarrelPhoton) {
    isBarrelMedium = passesSelectionAsMedium;
    isBarrelFake = passesSelectionAsFake;
  }
  else if (isEndcapPhoton) {
    isEndcapMedium = passesSelectionAsMedium;
    isEndcapFake = passesSelectionAsFake;
  }

  photonExaminationResultsStruct photonExaminationResults = photonExaminationResultsStruct(isBarrelMedium, isBarrelFake, isEndcapMedium, isEndcapFake, passesLeadingpTCut, ((photonsCollection.eta)->at(photonIndex)), ((photonsCollection.phi)->at(photonIndex)), ((photonsCollection.pT)->at(photonIndex)), ((photonsCollection.energy)->at(photonIndex)), MCScaleFactors);
  return photonExaminationResults;
}

MCExaminationResultsStruct examineMC(parametersStruct &parameters, const int& nMCParticles, const MCCollectionStruct& MCCollection) {
  int nPhotonsWithNeutralinoMom = 0;
  float generated_gluinoMass = 0.;
  bool gluinoMassIsSet = false;
  float generated_neutralinoMass = 0.;
  bool neutralinoMassIsSet = false;
  for (int MCIndex = 0; MCIndex < nMCParticles; ++MCIndex) {
    int particle_mcPID = (MCCollection.MCPIDs)->at(MCIndex);
    int particle_mcMomPID = (MCCollection.MCMomPIDs)->at(MCIndex);
    if ((particle_mcPID == parameters.PIDs.photon) && (particle_mcMomPID == parameters.PIDs.neutralino)) {
      if (passesBitMask((MCCollection.MCStatusFlags)->at(MCIndex), parameters.MCStatusFlagBitMask)) {
        ++nPhotonsWithNeutralinoMom;
      }
    }
    if (!(gluinoMassIsSet) && particle_mcPID == parameters.PIDs.gluino) {
      generated_gluinoMass = (MCCollection.MCMasses)->at(MCIndex);
      gluinoMassIsSet = true;
    }
    if (!(neutralinoMassIsSet) && particle_mcMomPID == parameters.PIDs.neutralino) {
      generated_neutralinoMass = (MCCollection.MCMomMasses)->at(MCIndex);
      neutralinoMassIsSet = true;
    }
  }
  bool passesMCSelection = (nPhotonsWithNeutralinoMom == 2);
  if (passesMCSelection && (!(gluinoMassIsSet && neutralinoMassIsSet))) {
    std::cout << "ERROR: Unable to find gluino or neutralino mass in an event that passes MC selection." << std::endl;
    std::exit(EXIT_FAILURE);
  }

  MCExaminationResultsStruct MCExaminationResults = MCExaminationResultsStruct(passesMCSelection, generated_gluinoMass, generated_neutralinoMass);
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

jetExaminationResultsStruct examineJet(optionsStruct &options, parametersStruct &parameters, countersStruct &counters, const jetsCollectionStruct& jetsCollection, const int& jetIndex, const float& generated_gluinoMass, const float& generated_neutralinoMass) {
  bool passesJetSelection = true;
  bool passesJetSelectionJECDown = false;
  bool passesJetSelectionJECUp = false;

  //Kinematic cuts: eta, pT
  float absEta = std::fabs((jetsCollection.eta)->at(jetIndex));
  applyCondition(counters, jetFailureCategory::eta, passesJetSelection, (absEta < parameters.jetEtaCut), options.isMC, generated_gluinoMass, generated_neutralinoMass);
  float jet_pT = ((jetsCollection.pT)->at(jetIndex));
  float jecFractionalUncertainty = 0.;
  if (options.isMC) {
    jecFractionalUncertainty = (jetsCollection.JECUncertainty)->at(jetIndex);
    float jet_pT_JECDown = (1.0 - jecFractionalUncertainty)*jet_pT;
    float jet_pT_JECUp = (1.0 + jecFractionalUncertainty)*jet_pT;
    passesJetSelectionJECDown = (jet_pT_JECDown > parameters.jetpTCut);
    passesJetSelectionJECUp = (jet_pT_JECUp > parameters.jetpTCut);
  }
  applyCondition(counters, jetFailureCategory::pT, passesJetSelection, (jet_pT > parameters.jetpTCut), options.isMC, generated_gluinoMass, generated_neutralinoMass);

  eventWeightsStruct prefireWeights = findPrefireWeights((jetsCollection.eta)->at(jetIndex), jet_pT, parameters.prefiringEfficiencyMap);

  // ID cuts: loose ID, PUID, jetID
  // applyCondition(counters, jetFailureCategory::PFLooseID, passesJetSelection, (((jetsCollection.looseID)->at(jetIndex))), options.isMC, generated_gluinoMass, generated_neutralinoMass); // disable until better understanding
  applyCondition(counters, jetFailureCategory::puID, passesJetSelection, (((jetsCollection.PUID)->at(jetIndex)) > parameters.jetPUIDThreshold), options.isMC, generated_gluinoMass, generated_neutralinoMass);
  applyCondition(counters, jetFailureCategory::jetID, passesJetSelection, (((jetsCollection.ID)->at(jetIndex)) == 6), options.isMC, generated_gluinoMass, generated_neutralinoMass);

  if (passesJetSelection) incrementCounters(miscCounter::passingJets, counters, options.isMC, generated_gluinoMass, generated_neutralinoMass);
  else incrementCounters(miscCounter::failingJets, counters, options.isMC, generated_gluinoMass, generated_neutralinoMass);
  incrementCounters(miscCounter::totalJets, counters, options.isMC, generated_gluinoMass, generated_neutralinoMass);

  jetExaminationResultsStruct result = jetExaminationResultsStruct(passesJetSelection, passesJetSelectionJECDown, passesJetSelectionJECUp, ((jetsCollection.eta)->at(jetIndex)), ((jetsCollection.phi)->at(jetIndex)), jet_pT, jecFractionalUncertainty, prefireWeights);
  return result;
}

float getMinDeltaR(const float& jetEta, const float& jetPhi, const std::vector<angularVariablesStruct>& selectedPhotonAnglesList) {
  float min_dR = -1.0;
  for (const auto& selectedPhotonAngles : selectedPhotonAnglesList) {
    float deltaEta = selectedPhotonAngles.eta - jetEta;
    float smallerPhi, largerPhi;
    if (selectedPhotonAngles.phi > jetPhi) {
      largerPhi = selectedPhotonAngles.phi;
      smallerPhi = jetPhi;
    }
    else {
      smallerPhi = selectedPhotonAngles.phi;
      largerPhi = jetPhi;
    }
    float deltaPhi_direction1 = largerPhi - smallerPhi;
    float deltaPhi_direction2 = constants::VALUEOFTWOPI - deltaPhi_direction1;
    float deltaPhi = (deltaPhi_direction1 < deltaPhi_direction2) ? deltaPhi_direction1 : deltaPhi_direction2;
    float dR = std::sqrt(deltaEta*deltaEta + deltaPhi*deltaPhi);
    if ((min_dR < 0) || (dR < min_dR)) min_dR = dR;
  }
  return min_dR;
}

eventExaminationResultsStruct examineEvent(optionsStruct &options, parametersStruct &parameters, countersStruct &counters, Long64_t& entryIndex, eventDetailsStruct& eventDetails, MCCollectionStruct &MCCollection, photonsCollectionStruct &photonsCollection, jetsCollectionStruct &jetsCollection) {
  int evt_nJetsDR = 0;
  float evt_ST = 0.0;
  eventWeightsStruct evt_prefireWeights(1.0f, 1.0f, 1.0f);
  eventWeightsStruct evt_photonMCScaleFactors(1.0f, 1.0f, 1.0f);
  std::map<shiftType, float> shifted_ST = empty_STMap();
  std::map<shiftType, int> shifted_nJetsDR = empty_NJetsMap();
  bool passesEventSelection = true;

  // Additional selection, only for MC
  float generated_gluinoMass = 0.;
  float generated_neutralinoMass = 0.;
  if (options.isMC) {
    MCExaminationResultsStruct MCExaminationResults = examineMC(parameters, (eventDetails.nMCParticles), MCCollection);
    generated_gluinoMass = MCExaminationResults.gluinoMass;
    generated_neutralinoMass = MCExaminationResults.neutralinoMass;
    applyCondition(counters, eventFailureCategory::MCGenInformation, passesEventSelection, MCExaminationResults.passesMCSelection, options.isMC, generated_gluinoMass, generated_neutralinoMass);
  }

  if (!(options.isMC) && parameters.HLTPhotonBit >= 0) { // Apply HLT photon selection iff input is not MC and HLTBit is set to a positive integer
    applyCondition(counters, eventFailureCategory::HLTPhoton, passesEventSelection, checkHLTBit(eventDetails.HLTPhotonBits, parameters.HLTPhotonBit), options.isMC, generated_gluinoMass, generated_neutralinoMass);
  }

  // Photon selection
  std::vector<angularVariablesStruct> selectedPhotonAnglesList;
  std::vector<TLorentzVector> selectedPhotonFourMomentaList;
  int nSelectedPhotonsPassingSubLeadingpTCut = 0;
  int nSelectedPhotonsPassingLeadingpTCut = 0;
  int nMediumPhotons = 0;
  int nFakePhotons = 0;
  int nVetoPhotons = 0;
  for (Int_t photonIndex = 0; photonIndex < (eventDetails.nPhotons); ++photonIndex) {
    photonExaminationResultsStruct photonExaminationResults = examinePhoton(options, parameters, counters, (eventDetails.eventRho), photonsCollection, photonIndex, generated_gluinoMass, generated_neutralinoMass);
    if (photonExaminationResults.isEndcapMedium) ++nVetoPhotons;
    else if (photonExaminationResults.isBarrelMedium || photonExaminationResults.isBarrelFake) {
      ++nSelectedPhotonsPassingSubLeadingpTCut;
      if (photonExaminationResults.passesLeadingpTCut) ++nSelectedPhotonsPassingLeadingpTCut;
      evt_ST += photonExaminationResults.pT;
      if (options.isMC) {
        for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
          shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
          addShiftedEToSTMap(photonExaminationResults.pT, shifted_ST, typeIndex); // no effect on photon's contribution to ST due to any of the shifts
        }
        evt_photonMCScaleFactors.nominal *= (photonExaminationResults.MCScaleFactors).nominal;
        evt_photonMCScaleFactors.down *= (photonExaminationResults.MCScaleFactors).down;
        evt_photonMCScaleFactors.up *= (photonExaminationResults.MCScaleFactors).up;
      }
      angularVariablesStruct angularVariables = angularVariablesStruct(photonExaminationResults.eta, photonExaminationResults.phi);
      selectedPhotonAnglesList.push_back(angularVariables);
      TLorentzVector photonFourMomentum;
      photonFourMomentum.SetPtEtaPhiE(photonExaminationResults.pT, photonExaminationResults.eta, photonExaminationResults.phi, photonExaminationResults.energy);
      selectedPhotonFourMomentaList.push_back(photonFourMomentum);
    }
    if (photonExaminationResults.isBarrelMedium) ++nMediumPhotons;
    else if (photonExaminationResults.isBarrelFake) ++nFakePhotons;
  }

  if (options.photonSelectionType == PhotonSelectionType::singlemedium) {
    applyCondition(counters, eventFailureCategory::lowEnergyPhotons, passesEventSelection, (nSelectedPhotonsPassingLeadingpTCut >= 1), options.isMC, generated_gluinoMass, generated_neutralinoMass);
    applyCondition(counters, eventFailureCategory::wrongNMediumPhotons, passesEventSelection, (nMediumPhotons == 1), options.isMC, generated_gluinoMass, generated_neutralinoMass);
    applyCondition(counters, eventFailureCategory::endcapPhotonVeto, passesEventSelection, (nVetoPhotons == 0), options.isMC, generated_gluinoMass, generated_neutralinoMass);
  }

  else {
    applyCondition(counters, eventFailureCategory::lowEnergyPhotons, passesEventSelection, ((nSelectedPhotonsPassingSubLeadingpTCut >= 2) && (nSelectedPhotonsPassingLeadingpTCut >= 1)), options.isMC, generated_gluinoMass, generated_neutralinoMass);
    if (options.photonSelectionType == PhotonSelectionType::medium) {
      applyCondition(counters, eventFailureCategory::wrongNMediumPhotons, passesEventSelection, (nMediumPhotons == 2), options.isMC, generated_gluinoMass, generated_neutralinoMass);
    }
    else if (options.photonSelectionType == PhotonSelectionType::mediumfake) { // nMediumPhotons != 2
      applyCondition(counters, eventFailureCategory::wrongNMediumPhotons, passesEventSelection, (nMediumPhotons == 1) && ((nMediumPhotons + nFakePhotons) >= 2), options.isMC, generated_gluinoMass, generated_neutralinoMass);
      applyCondition(counters, eventFailureCategory::endcapPhotonVeto, passesEventSelection, (nVetoPhotons == 0), options.isMC, generated_gluinoMass, generated_neutralinoMass);
    }
    else if (options.photonSelectionType == PhotonSelectionType::fake) { // nMediumPhotons != 2 and (nMediumPhotons != 1 or (nMediumPhotons + nFakePhotons) < 2)
      applyCondition(counters, eventFailureCategory::wrongNMediumPhotons, passesEventSelection, (nMediumPhotons == 0) && (nFakePhotons >= 2), options.isMC, generated_gluinoMass, generated_neutralinoMass);
      applyCondition(counters, eventFailureCategory::endcapPhotonVeto, passesEventSelection, (nVetoPhotons == 0), options.isMC, generated_gluinoMass, generated_neutralinoMass);
    }
  }

  // Apply invariant mass cut only in the double photon selections
  float evt_invariantMass = -1.0;
  if ((nSelectedPhotonsPassingSubLeadingpTCut >= 2) && (options.photonSelectionType != PhotonSelectionType::singlemedium)) {
    evt_invariantMass = getDiphotonInvariantMass(selectedPhotonFourMomentaList);
    applyCondition(counters, eventFailureCategory::lowInvariantMass, passesEventSelection, (evt_invariantMass >= parameters.invariantMassCut), options.isMC, generated_gluinoMass, generated_neutralinoMass);
  }

  // Jet selection
  float evt_HT = 0;
  for (Int_t jetIndex = 0; jetIndex < (eventDetails.nJets); ++jetIndex) {
    jetExaminationResultsStruct jetExaminationResults = examineJet(options, parameters, counters, jetsCollection, jetIndex, generated_gluinoMass, generated_neutralinoMass);
    evt_prefireWeights.nominal *= (jetExaminationResults.prefireWeights).nominal; // All jets, whether or not they pass any of the cuts, contribute to the prefiring weight
    evt_prefireWeights.down *= (jetExaminationResults.prefireWeights).down;
    evt_prefireWeights.up *= (jetExaminationResults.prefireWeights).up;
    float min_dR = getMinDeltaR(jetExaminationResults.eta, jetExaminationResults.phi, selectedPhotonAnglesList);
    bool passesDeltaRCut = ((min_dR > parameters.minDeltaRCut) || (min_dR < 0.0));
    if (!passesDeltaRCut) {
      incrementCounters(jetFailureCategory::deltaR, counterType::global, counters, options.isMC, generated_gluinoMass, generated_neutralinoMass);
      if (jetExaminationResults.passesSelectionJECNominal) incrementCounters(jetFailureCategory::deltaR, counterType::differential, counters, options.isMC, generated_gluinoMass, generated_neutralinoMass);
    }
    if (jetExaminationResults.passesSelectionJECNominal) {
      evt_HT += jetExaminationResults.pT; // Add hT whether or not jet passes deltaR check
      if (passesDeltaRCut) {
        evt_ST += jetExaminationResults.pT; // Add to sT only if jet passes deltaR check, to avoid double-counting
        ++evt_nJetsDR; // Count only those jets that are sufficiently away from a photon
      }
    }
    if (options.isMC && passesDeltaRCut) {
      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
        shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
        float shifted_contribution = (jetExaminationResults.pT);
        bool passes_selection = jetExaminationResults.passesSelectionJECNominal;
        if (typeIndex == shiftType::JECDown) {
          passes_selection = jetExaminationResults.passesSelectionJECDown;
          shifted_contribution = (1.0 - (jetExaminationResults.jecFractionalUncertainty))*(jetExaminationResults.pT);
        }
        else if (typeIndex == shiftType::JECUp) {
          passes_selection = jetExaminationResults.passesSelectionJECUp;
          shifted_contribution = (1.0 + (jetExaminationResults.jecFractionalUncertainty))*(jetExaminationResults.pT);
        }
        if (passes_selection) {
          addShiftedEToSTMap(shifted_contribution, shifted_ST, typeIndex);
          incrementNJetsMap(shifted_nJetsDR, typeIndex);
        }
      }
    }
  }
  int max_nJets = evt_nJetsDR;
  if (options.isMC) { // this makes sure that the nJets used to make the decision whether or not to save the event is the maximum nJets accounting for all the shifts
    int maxNJetsShifted = getMaxNJets(shifted_nJetsDR);
    if (maxNJetsShifted > max_nJets) max_nJets = maxNJetsShifted;
  }
  applyCondition(counters, eventFailureCategory::wrongNJets, passesEventSelection, (max_nJets >= 2), options.isMC, generated_gluinoMass, generated_neutralinoMass);
  applyCondition(counters, eventFailureCategory::hTCut, passesEventSelection, (evt_HT >= parameters.HTCut), options.isMC, generated_gluinoMass, generated_neutralinoMass);

  // Add MET to ST
  evt_ST += eventDetails.PFMET;

  // Add shifted energies
  if (options.isMC) {
    addShiftedEToSTMap(eventDetails.PFMET, shifted_ST, shiftType::JECDown);
    addShiftedEToSTMap(eventDetails.PFMET, shifted_ST, shiftType::JECUp);
    addShiftedEToSTMap(eventDetails.PFMET_UnclusteredDown, shifted_ST, shiftType::UnclusteredMETDown);
    addShiftedEToSTMap(eventDetails.PFMET_UnclusteredUp, shifted_ST, shiftType::UnclusteredMETUp);
    addShiftedEToSTMap(eventDetails.PFMET_JERDown, shifted_ST, shiftType::JERMETDown);
    addShiftedEToSTMap(eventDetails.PFMET_JERUp, shifted_ST, shiftType::JERMETUp);
  }

  eventExaminationResultsStruct eventResult = eventExaminationResultsStruct(entryIndex, passesEventSelection, evt_ST, evt_nJetsDR, evt_prefireWeights, evt_photonMCScaleFactors, shifted_ST, shifted_nJetsDR);
  return eventResult;
}

std::vector<eventExaminationResultsStruct> getSelectedEventsWithInfo(optionsStruct &options, parametersStruct &parameters, countersStruct &counters) {
  std::vector<eventExaminationResultsStruct> selectedEventsInfo;

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

    eventExaminationResultsStruct eventExaminationResults = examineEvent(options, parameters, counters, entryIndex, eventDetails, MCCollection, photonsCollection, jetsCollection);
    bool passesEventSelection = eventExaminationResults.passesSelection;
    incrementCounters(miscCounter::totalEvents, counters, false, 0., 0.);
    if (!(passesEventSelection)) {
      incrementCounters(miscCounter::failingEvents, counters, false, 0., 0.);
      continue;
    }
    incrementCounters(miscCounter::acceptedEvents, counters, false, 0., 0.);
    selectedEventsInfo.push_back(eventExaminationResults);
  }
  progressBar.terminate();
  return selectedEventsInfo;
}

void writeSelectedEventsToFile(optionsStruct &options, TFile *outputFile, const std::vector<eventExaminationResultsStruct>& selectedEventsInfo) {
  std::cout << "Beginning to write selected events to file..." << std::endl;
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

    outputTree->Fill();
  }
  progressBar.terminate();
  outputFile->Write();
}

int main(int argc, char* argv[]) {
  gROOT->SetBatch();
  tmArgumentParser argumentParser = tmArgumentParser("Calculate systematics due to uncertainty on jet energy corrections.");
  argumentParser.addArgument("inputFilesList", "", true, "Path to file containing list of input files.");
  argumentParser.addArgument("outputFilePath", "", true, "Path to output file.");
  argumentParser.addArgument("isMC", "false", false, "Input file is a MC sample -- disable HLT photon trigger and enable additional MC selection.");
  argumentParser.addArgument("counterStartInclusive", "", true, "Event number from input file from which to start. The event with this index is included in the processing.");
  argumentParser.addArgument("counterEndInclusive", "", true, "Event number from input file at which to end. The event with this index is included in the processing.");
  argumentParser.addArgument("photonSelectionType", "fake", true, "Photon selection type: can be any one of: \"fake\", \"medium\", \"mediumfake\", or \"singlemedium\".");
  argumentParser.addArgument("year", "2017", false, "Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.");
  argumentParser.addArgument("inputFile_STRegionBoundaries", "STRegionBoundaries.dat", false, "Path to file with ST region boundaries. First bin is the normalization bin, and the last bin is the last boundary to infinity."); // for trigger efficiency studies
  // all remaining arguments are only used in MC samples to construct the histograms that help in diagnosing efficiency issues.
  argumentParser.addArgument("nGluinoMassBins", "20", false, "nBins on the gluino mass axis"); // (800 - 25) GeV --> (1750 + 25) GeV in steps of 50 GeV
  argumentParser.addArgument("minGluinoMass", "775.0", false, "Min gluino mass for the 2D plots.");
  argumentParser.addArgument("maxGluinoMass", "1775.0", false, "Max gluino mass for the 2D plots.");
  argumentParser.addArgument("nNeutralinoMassBins", "133", false, "nBins on the neutralino mass axis.");
  argumentParser.addArgument("minNeutralinoMass", "93.75", false, "Min neutralino mass for the 2D plots.");
  argumentParser.addArgument("maxNeutralinoMass", "1756.25", false, "Max neutralino mass for the 2D plots."); // (100 - 6.25) GeV --> (1750 + 6.25) GeV in steps of 12.5 GeV
  argumentParser.setPassedStringValues(argc, argv);

  optionsStruct options = getOptionsFromParser(argumentParser);

  STRegionsStruct STRegions(options.inputFile_STRegionBoundaries);

  parametersStruct parameters = parametersStruct();
  parameters.tuneParametersForYear(options.year, options.isMC);

  countersStruct counters = countersStruct();
  initializeCounters(counters, options, STRegions.nSTSignalBins);

  std::stringstream optionsStringstream;
  optionsStringstream << options;
  TNamed *optionsObject = new TNamed("optionsString", optionsStringstream.str().c_str());
  std::stringstream parametersStringstream;
  parametersStringstream << parameters;
  TNamed *parametersObject = new TNamed("parametersString", parametersStringstream.str().c_str());
  TFile *outputFile = TFile::Open(options.outputFilePath.c_str(), "RECREATE");
  if (!(outputFile->IsOpen()) || outputFile->IsZombie()) {
    std::cout << "ERROR: Unable to open output file to write. File path: " << options.outputFilePath << std::endl;
  }

  std::vector<eventExaminationResultsStruct> selectedEventsInfo = getSelectedEventsWithInfo(options, parameters, counters);

  writeSelectedEventsToFile(options, outputFile, selectedEventsInfo);

  outputFile->WriteTObject(parametersObject);
  outputFile->WriteTObject(optionsObject);
  outputFile->Close();

  std::cout << getNDashes(100) << std::endl;
  printAndSaveCounters(counters, options.isMC, "MCStatisticsDetails.root", "triggerEfficiencyRawEventCounters.root");

  std::cout << getNDashes(100) << std::endl
            << "Options:" << std::endl
            << optionsStringstream.str() << std::endl
            << getNDashes(100) << std::endl;
  std::cout << getNDashes(100) << std::endl
            << "Parameters:" << std::endl
            << parametersStringstream.str() << std::endl
            << getNDashes(100) << std::endl;
  return 0;
}
