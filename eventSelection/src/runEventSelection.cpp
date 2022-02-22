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
    std::cout << "ERROR: Area counter value = " << areaCounter << "; getRhoCorrectedIsolation called for type = " << PFTypesForEANames[PFType] << ", eta = " << absEta << ", which is above all available upper bounds." << std::endl;
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

photonExaminationResultsStruct examinePhoton(optionsStruct &options, parametersStruct &parameters, const float& rho, const photonsCollectionStruct& photonsCollection, const int& photonIndex, std::vector<angularVariablesStruct> &selectedTruePhotonAngles) {
  photonExaminationResultsStruct results;
  photonProperties& properties = results.pho_properties;
  eventWeightsStruct& scaleFactors = results.MCScaleFactors;

  std::map<mediumPhotonCriterion, bool> medium_bits;
  std::map<vetoedPhotonCriterion, bool> vetoed_bits;
  std::map<fakePhotonCriterion, bool> fake_bits;

  // mva
  properties[photonProperty::mva] = (photonsCollection.mva)->at(photonIndex);

  // Kinematic cuts
  properties[photonProperty::eta] = (photonsCollection.eta)->at(photonIndex);
  properties[photonProperty::phi] = (photonsCollection.phi)->at(photonIndex);
  angularVariablesStruct photonAngle(properties[photonProperty::eta], properties[photonProperty::phi]);
  properties[photonProperty::deltaR_nearestTruePhoton] = photonAngle.getMinDeltaR(selectedTruePhotonAngles);
  float absEta = std::fabs(properties[photonProperty::eta]);
  bool passesEta = (absEta < parameters.photonBarrelEtaCut);
  medium_bits[mediumPhotonCriterion::eta] = passesEta;
  vetoed_bits[vetoedPhotonCriterion::eta] = passesEta;
  fake_bits[fakePhotonCriterion::eta] = passesEta;

  properties[photonProperty::pT] = std::fabs((photonsCollection.pT)->at(photonIndex));
  bool passesPT = (properties[photonProperty::pT] > parameters.pTCutSubLeading);
  medium_bits[mediumPhotonCriterion::pT] = passesPT;
  vetoed_bits[vetoedPhotonCriterion::pT] = passesPT;
  fake_bits[fakePhotonCriterion::pT] = passesPT;

  // // Electron veto (old, conversion-safe)
  // bool passesConvSafeVetoRaw = (((photonsCollection.electronVeto)->at(photonIndex)) == (Int_t)(true));
  // bool passesConvSafeVeto = passesConvSafeVetoRaw;
  // if (options.invertElectronVeto) passesConvSafeVeto = !(passesConvSafeVetoRaw); // this is easier than renaming criteria and changing the names in a dozen places
  // medium_bits[mediumPhotonCriterion::conversionSafeElectronVeto] = passesConvSafeVeto;
  // vetoed_bits[vetoedPhotonCriterion::conversionSafeElectronVeto] = passesConvSafeVeto;
  // fake_bits[fakePhotonCriterion::conversionSafeElectronVeto] = passesConvSafeVeto;

  // Pixel veto
  bool passesPixelVetoRaw = (((photonsCollection.hasPixelSeed)->at(photonIndex)) == (Int_t)(false));
  bool passesPixelVeto = passesPixelVetoRaw;
  if (options.invertElectronVeto) passesPixelVeto = !(passesPixelVetoRaw); // this is easier than renaming criteria and changing the names in a dozen places
  medium_bits[mediumPhotonCriterion::conversionSafeElectronVeto] = passesPixelVeto;
  vetoed_bits[vetoedPhotonCriterion::conversionSafeElectronVeto] = passesPixelVeto;
  fake_bits[fakePhotonCriterion::conversionSafeElectronVeto] = passesPixelVeto;

  // Quality cuts
  photonQualityCutsStruct* qualityCuts = &(parameters.photonQualityCutsBarrel);
  if (absEta > parameters.photonBarrelEtaCut) qualityCuts = &(parameters.photonQualityCutsEndcap);

  properties[photonProperty::hOverE] = (photonsCollection.HOverE)->at(photonIndex);
  bool passesHOverE = (properties[photonProperty::hOverE] < qualityCuts->towerHOverE);
  medium_bits[mediumPhotonCriterion::hOverE] = passesHOverE;
  bool passesHOverELoose = (properties[photonProperty::hOverE] < qualityCuts->towerHOverELoose);
  vetoed_bits[vetoedPhotonCriterion::hOverE] = passesHOverELoose;

  float pTDependentNeutralIsolationCut = (qualityCuts->neutralIsolation).getPolynomialValue(properties[photonProperty::pT]);
  properties[photonProperty::rhoCorrectedNeutralIsolation] = getRhoCorrectedIsolation(((photonsCollection.PFNeutralIsolationUncorrected)->at(photonIndex)), PFTypesForEA::neutralHadron, absEta, rho, parameters.effectiveAreas);
  bool passesNeutralIsolation = (properties[photonProperty::rhoCorrectedNeutralIsolation] < pTDependentNeutralIsolationCut);
  medium_bits[mediumPhotonCriterion::neutralIsolation] = passesNeutralIsolation;
  float pTDependentNeutralIsolationCutLoose = (qualityCuts->neutralIsolationLoose).getPolynomialValue(properties[photonProperty::pT]);
  bool passesNeutralIsolationLoose = (properties[photonProperty::rhoCorrectedNeutralIsolation] < pTDependentNeutralIsolationCutLoose);
  vetoed_bits[vetoedPhotonCriterion::neutralIsolation] = passesNeutralIsolationLoose;

  float pTDependentPhotonIsolationCut = (qualityCuts->photonIsolation).getPolynomialValue(properties[photonProperty::pT]);
  properties[photonProperty::rhoCorrectedPhotonIsolation] = getRhoCorrectedIsolation(((photonsCollection.PFPhotonIsolationUncorrected)->at(photonIndex)), PFTypesForEA::photon, absEta, rho, parameters.effectiveAreas);
  bool passesPhotonIsolation = (properties[photonProperty::rhoCorrectedPhotonIsolation] < pTDependentPhotonIsolationCut);
  medium_bits[mediumPhotonCriterion::photonIsolation] = passesPhotonIsolation;
  float pTDependentPhotonIsolationCutLoose = (qualityCuts->photonIsolationLoose).getPolynomialValue(properties[photonProperty::pT]);
  bool passesPhotonIsolationLoose = (properties[photonProperty::rhoCorrectedPhotonIsolation] < pTDependentPhotonIsolationCutLoose);
  vetoed_bits[vetoedPhotonCriterion::photonIsolation] = passesPhotonIsolationLoose;

  // vetoed_bits[vetoedPhotonCriterion::passesNeutIsoAndPhoIsoLooseCriteria] = (passesNeutralIsolationLoose && passesPhotonIsolationLoose);
  // fake_bits[fakePhotonCriterion::passesNeutIsoAndPhoIsoLooseCriteria] = (passesNeutralIsolationLoose && passesPhotonIsolationLoose);

  properties[photonProperty::rawChargedIsolation] = (photonsCollection.PFChargedIsolationUncorrected)->at(photonIndex);
  properties[photonProperty::rhoCorrectedChargedIsolation] = getRhoCorrectedIsolation(((photonsCollection.PFChargedIsolationUncorrected)->at(photonIndex)), PFTypesForEA::chargedHadron, absEta, rho, parameters.effectiveAreas);
  bool passesChargedIsolation = (properties[photonProperty::rhoCorrectedChargedIsolation] < qualityCuts->chargedIsolation);
  medium_bits[mediumPhotonCriterion::chargedIsolation] = passesChargedIsolation;
  // vetoed_bits[vetoedPhotonCriterion::chIsoBetweenMedAndLoose] = ((!(passesChargedIsolation)) && (properties[photonProperty::rhoCorrectedChargedIsolation] < qualityCuts->chargedIsolationLoose));
  // fake_bits[fakePhotonCriterion::chIsoBetweenLooseAndExtraLoose] = ((properties[photonProperty::rhoCorrectedChargedIsolation] >= qualityCuts->chargedIsolationLoose) && (properties[photonProperty::rhoCorrectedChargedIsolation] <= qualityCuts->chargedIsolationExtraLoose));
  bool passesChargedIsolationLoose = (properties[photonProperty::rhoCorrectedChargedIsolation] < qualityCuts->chargedIsolationLoose);
  vetoed_bits[vetoedPhotonCriterion::chargedIsolation] = passesChargedIsolationLoose;
  fake_bits[fakePhotonCriterion::passesInvertedChIso] = (properties[photonProperty::rhoCorrectedChargedIsolation] >= qualityCuts->chargedIsolationLoose);

  properties[photonProperty::sigmaIEtaIEta] = ((photonsCollection.sigmaIEtaIEta)->at(photonIndex));
  bool passesSigmaIEtaIEta = (properties[photonProperty::sigmaIEtaIEta] < qualityCuts->sigmaIEtaIEta);
  medium_bits[mediumPhotonCriterion::sigmaIEtaIEta] = passesSigmaIEtaIEta;
  bool passesSigmaIEtaIEtaLoose = (properties[photonProperty::sigmaIEtaIEta] < qualityCuts->sigmaIEtaIEtaLoose);
  vetoed_bits[vetoedPhotonCriterion::sigmaIEtaIEta] = passesSigmaIEtaIEtaLoose;

  // fake_bits[fakePhotonCriterion::passesOtherLooseCuts] = (passesHOverELoose || passesSigmaIEtaIEtaLoose || passesNeutralIsolationLoose);

  bool passes_showerShapeMedIDCuts = ((properties[photonProperty::sigmaIEtaIEta] < qualityCuts->sigmaIEtaIEta) && (properties[photonProperty::hOverE] < qualityCuts->towerHOverE));
  // bool passes_showerShapeLooseIDCuts = ((properties[photonProperty::sigmaIEtaIEta] < qualityCuts->sigmaIEtaIEtaLoose) && (properties[photonProperty::hOverE] < qualityCuts->towerHOverELoose));
  // vetoed_bits[vetoedPhotonCriterion::passesShowerShapeLooseIDCuts] = passes_showerShapeLooseIDCuts;
  fake_bits[fakePhotonCriterion::passesShowerShapeMedIDCuts] = passes_showerShapeMedIDCuts;

  properties[photonProperty::R9] = ((photonsCollection.R9)->at(photonIndex));
  properties[photonProperty::ecalClusIso] = ((photonsCollection.ecalClusIso)->at(photonIndex));
  properties[photonProperty::trkIso] = ((photonsCollection.trkIso)->at(photonIndex));
  properties[photonProperty::energy] = (photonsCollection.energy)->at(photonIndex);

  int nFalseBits_medium = getNFalseBits(medium_bits);
  bool fails_mediumID = (nFalseBits_medium >= 1);
  vetoed_bits[vetoedPhotonCriterion::failsMediumID] = fails_mediumID;
  fake_bits[fakePhotonCriterion::failsMediumID] = fails_mediumID;

  assert(static_cast<int>(medium_bits.size()) == static_cast<int>(mediumPhotonCriterion::nMediumPhotonCriteria));
  assert(static_cast<int>(vetoed_bits.size()) == static_cast<int>(vetoedPhotonCriterion::nVetoedPhotonCriteria));
  assert(static_cast<int>(fake_bits.size()) == static_cast<int>(fakePhotonCriterion::nFakePhotonCriteria));

  int nFalseBits_vetoed = getNFalseBits(vetoed_bits);
  int nFalseBits_fake = getNFalseBits(fake_bits);
  if (nFalseBits_medium == 0) {
    assert((nFalseBits_vetoed != 0) && (nFalseBits_fake != 0));
    results.photon_type = photonType::medium;
  }
  else if (nFalseBits_vetoed == 0) {
    assert(nFalseBits_medium != 0);
    results.photon_type = photonType::vetoed;
  }
  else if (nFalseBits_fake == 0) {
    assert((nFalseBits_medium != 0) && (nFalseBits_vetoed != 0));
    results.photon_type = photonType::fake;
  }
  results.isMarginallyUnselected = ((nFalseBits_medium == 1) ||
				    (nFalseBits_vetoed == 1) ||
				    (nFalseBits_fake == 1));
  if (results.isMarginallyUnselected) {
    if (nFalseBits_medium == 1) {
      results.marginallyUnselectedMediumCriterion = getFirstFalseCriterion(medium_bits);
    }
    if (nFalseBits_vetoed == 1) {
      results.marginallyUnselectedVetoedCriterion = getFirstFalseCriterion(vetoed_bits);
    }
    if (nFalseBits_fake == 1) {
      results.marginallyUnselectedFakeCriterion = getFirstFalseCriterion(fake_bits);
    }
  }

  if (options.calculateMCScaleFactorWeights && (results.photon_type != photonType::nPhotonTypes)) {
    if (results.photon_type == photonType::fake) {
      scaleFactors = findMCScaleFactors(((photonsCollection.eta)->at(photonIndex)), ((photonsCollection.pT)->at(photonIndex)), (parameters.photonMCScaleFactorsMaps).at(photonType::vetoed));
    }
    else {
      scaleFactors = findMCScaleFactors(((photonsCollection.eta)->at(photonIndex)), ((photonsCollection.pT)->at(photonIndex)), (parameters.photonMCScaleFactorsMaps).at(results.photon_type));
    }
  }

  results.energy = (photonsCollection.energy)->at(photonIndex);

  // results.contributesToMisc2DHistograms = (passesEta && passesPT && passesConvSafeVeto && fails_mediumID && fake_bits[fakePhotonCriterion::passesShowerShapeMedIDCuts]);
  results.contributesToMisc2DHistograms = (passesEta && passesPT && passesPixelVeto && fails_mediumID && fake_bits[fakePhotonCriterion::passesShowerShapeMedIDCuts]);

  return results;
}

MCExaminationResultsStruct examineMCParticle(optionsStruct &options, parametersStruct &parameters, const MCCollectionStruct& MCCollection, const int& MCIndex, const bool& doPromptOnlyBitMask) {
  MCExaminationResultsStruct MCExaminationResults;
  int particle_mcPID = (MCCollection.MCPIDs)->at(MCIndex);
  // int custom_particle_ID = PIDUtils::getCustomParticleID(particle_mcPID);
  int particle_mcMomPID = (MCCollection.MCMomPIDs)->at(MCIndex);
  int custom_mom_ID = PIDUtils::getCustomParticleID(particle_mcMomPID);
  UShort_t particle_statusFlag = static_cast<UShort_t>(((MCCollection.MCStatusFlags)->at(MCIndex))&(static_cast<UShort_t>(7u))); // picks out only first 3 bits
  bool hasRequiredMom = (PIDUtils::isNeutralinoPID(particle_mcMomPID) || PIDUtils::isHiggsPID(particle_mcMomPID));
  if (PIDUtils::isPhotonPID(particle_mcPID) && hasRequiredMom) {
    bool passes_bit_mask = passesBitMask(particle_statusFlag, parameters.MCStatusFlagBitMask);
    if (doPromptOnlyBitMask) {
      passes_bit_mask = passesBitMask(particle_statusFlag, parameters.MCStatusFlagBitMask_promptOnly);
    }
    if (passes_bit_mask) {
      MCExaminationResults.isPhotonWithDesiredMom = true;
      truthPhotonProperties& pho_properties = MCExaminationResults.truth_photon_properties;
      pho_properties[truthPhotonProperty::eta] = (MCCollection.MCEtas)->at(MCIndex);
      pho_properties[truthPhotonProperty::phi] = (MCCollection.MCPhis)->at(MCIndex);
      pho_properties[truthPhotonProperty::pT] = (MCCollection.MCEts)->at(MCIndex);
      pho_properties[truthPhotonProperty::status] = (MCCollection.MCStatuses)->at(MCIndex);
      // assert(static_cast<int>((MCExaminationResults.truth_photon_properties).size()) == static_cast<int>(truthPhotonProperty::nTruthPhotonProperties)); // distance to nearest truth jet candidate needs to be set later, do this check then
    }
  }
  if (!((options.selectionType == "MC_stealth_t5") || (options.selectionType == "MC_stealth_t6"))) { // fill in "fake" eventProgenitor and neutralino masses for non-stealth events
    MCExaminationResults.eventProgenitorMass = 1500.;
    MCExaminationResults.neutralinoMass = 800.;
  }

  if (PIDUtils::isJetCandidatePID(particle_mcPID) && ((MCCollection.MCStatuses)->at(MCIndex) == parameters.jetCandidateStatusConstraint)) {
    if (options.MC_eventProgenitor == "gluino") {
      if (PIDUtils::isGluinoPID(particle_mcMomPID)) MCExaminationResults.isJetCandidateFromEventProgenitor = true;
    }
    else if (options.MC_eventProgenitor == "squark") {
      if (PIDUtils::isSquarkPID(particle_mcMomPID)) MCExaminationResults.isJetCandidateFromEventProgenitor = true;
    }
    if (PIDUtils::isSingletPID(particle_mcMomPID)) MCExaminationResults.isJetCandidateFromSinglet = true;
    if (MCExaminationResults.isJetCandidateFromEventProgenitor || MCExaminationResults.isJetCandidateFromSinglet) {
      MCExaminationResults.isJetCandidateFromStealthSource = true;
      truthJetCandidateProperties& jetCandidate_properties = MCExaminationResults.truth_jetCandidate_properties;
      jetCandidate_properties[truthJetCandidateProperty::eta] = (MCCollection.MCEtas)->at(MCIndex);
      jetCandidate_properties[truthJetCandidateProperty::phi] = (MCCollection.MCPhis)->at(MCIndex);
      jetCandidate_properties[truthJetCandidateProperty::pT] = (MCCollection.MCEts)->at(MCIndex);
      jetCandidate_properties[truthJetCandidateProperty::momID] = custom_mom_ID;
      jetCandidate_properties[truthJetCandidateProperty::status] = (MCCollection.MCStatuses)->at(MCIndex);
      jetCandidate_properties[truthJetCandidateProperty::statusFlag] = particle_statusFlag;
      // assert(static_cast<int>((MCExaminationResults.truth_jetCandidate_properties).size()) == static_cast<int>(truthJetCandidateProperty::nTruthJetCandidateProperties)); // distance to nearest true photon needs to be set later, do this check then
    }
  }

  if (options.MC_eventProgenitor == "gluino") {
    if (PIDUtils::isGluinoPID(particle_mcPID)) MCExaminationResults.eventProgenitorMass = (MCCollection.MCMasses)->at(MCIndex);
  }
  else if (options.MC_eventProgenitor == "squark") {
    if (PIDUtils::isSquarkPID(particle_mcPID)) MCExaminationResults.eventProgenitorMass = (MCCollection.MCMasses)->at(MCIndex);
  }
  if (PIDUtils::isNeutralinoPID(particle_mcMomPID)) MCExaminationResults.neutralinoMass = (MCCollection.MCMomMasses)->at(MCIndex);

  return MCExaminationResults;
}

float getDiphotonInvariantMass(const std::vector<TLorentzVector>& list_selectedPhotonFourMomenta) {
  TLorentzVector eventSum;
  for (const auto& selectedPhotonFourMomentum : list_selectedPhotonFourMomenta) {
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

jetExaminationResultsStruct examineJet(optionsStruct &options, parametersStruct &parameters, const jetsCollectionStruct& jetsCollection, const int& jetIndex, std::vector<angularVariablesStruct> &list_selectedPhotonAngles, std::vector<angularVariablesStruct> &selectedTruePhotonAngles, std::vector<angularVariablesStruct> &selectedTrueJetCandidateAngles_all) {
  jetExaminationResultsStruct results;

  if (options.saveMCObjects) {
    genJetProperties& gen_properties = results.gen_jet_properties;
    gen_properties[genJetProperty::pT] = (jetsCollection.jetGenPT)->at(jetIndex);
    gen_properties[genJetProperty::eta] = (jetsCollection.jetGenEta)->at(jetIndex);
    gen_properties[genJetProperty::phi] = (jetsCollection.jetGenPhi)->at(jetIndex);
    gen_properties[genJetProperty::partonID] = 1.0*PIDUtils::getCustomParticleID((jetsCollection.jetGenPartonID)->at(jetIndex));
    gen_properties[genJetProperty::partonMomID] = 1.0*PIDUtils::getCustomParticleID((jetsCollection.jetGenPartonMomID)->at(jetIndex));
    assert(static_cast<int>(gen_properties.size()) == static_cast<int>(genJetProperty::nGenJetProperties));

    results.hasGenVariablesSet = ((gen_properties[genJetProperty::pT] > -100.) &&
                                  (gen_properties[genJetProperty::eta] > -100.) &&
                                  (gen_properties[genJetProperty::phi] > -100.)); // -999 is the default value set by the n-tuplizer
    if (results.hasGenVariablesSet) {
      if (options.MC_eventProgenitor == "gluino") {
	results.hasEventProgenitorPartonMom = PIDUtils::isGluinoPID((jetsCollection.jetGenPartonMomID)->at(jetIndex));
      }
      else if (options.MC_eventProgenitor == "squark") {
	results.hasEventProgenitorPartonMom = PIDUtils::isSquarkPID((jetsCollection.jetGenPartonMomID)->at(jetIndex));
      }
      results.hasSingletPartonMom = PIDUtils::isSingletPID((jetsCollection.jetGenPartonMomID)->at(jetIndex));
    }
  }

  jetProperties& properties = results.jet_properties;

  std::map<jetCriterion, bool> bits;
  bool passesPT_nominal = false;
  bool passesPT_JECDown = false;
  bool passesPT_JECUp = false;
  bool passesPT_missingHEMDown = false;
  bool passesPT_missingHEMUp = false; // for syntactic consistency

  //Kinematic cuts: eta, pT
  properties[jetProperty::eta] = (jetsCollection.eta)->at(jetIndex);
  properties[jetProperty::phi] = (jetsCollection.phi)->at(jetIndex);

  angularVariablesStruct jetAngle = angularVariablesStruct(properties[jetProperty::eta], properties[jetProperty::phi]);
  float minDeltaR = jetAngle.getMinDeltaR(list_selectedPhotonAngles);
  properties[jetProperty::deltaR_nearestCaloPhoton] = minDeltaR;
  bits[jetCriterion::deltaR_photon] = ((minDeltaR > parameters.deltaRScale_jetPhotonDistance) || (minDeltaR < 0.));
  results.isAwayFromCaloPhoton = bits[jetCriterion::deltaR_photon];
  properties[jetProperty::deltaR_nearestTruePhoton] = -0.005;
  properties[jetProperty::deltaR_nearestTrueJetCandidate] = -0.005;
  properties[jetProperty::truthPTRatio] = -1.0;
  properties[jetProperty::deltaR_genJet] = -0.005;
  if (options.saveMCObjects) {
    properties[jetProperty::deltaR_nearestTruePhoton] = jetAngle.getMinDeltaR(selectedTruePhotonAngles);
    results.isCloseToTruePhoton = (((properties[jetProperty::deltaR_nearestTruePhoton]) > 0.) && (properties[jetProperty::deltaR_nearestTruePhoton] <= parameters.deltaRScale_jetPhotonDistance));
    properties[jetProperty::deltaR_nearestTrueJetCandidate] = jetAngle.getMinDeltaR(selectedTrueJetCandidateAngles_all);
    if (results.hasGenVariablesSet) {
      properties[jetProperty::truthPTRatio] = ((jetsCollection.pT)->at(jetIndex))/((results.gen_jet_properties)[genJetProperty::pT]);
      angularVariablesStruct gen_angle = angularVariablesStruct((results.gen_jet_properties)[genJetProperty::eta], (results.gen_jet_properties)[genJetProperty::phi]);
      properties[jetProperty::deltaR_genJet] = jetAngle.get_deltaR(gen_angle);
    }
  }

  float absEta = std::fabs(properties[jetProperty::eta]);
  bool passes_eta = (absEta < parameters.jetEtaCut);
  bits[jetCriterion::eta] = passes_eta;

  float jet_pT = ((jetsCollection.pT)->at(jetIndex));
  properties[jetProperty::pT] = jet_pT;
  passesPT_nominal = (jet_pT > parameters.jetpTCut);

  float jecFractionalUncertainty = 0.;
  if (options.calculateShiftedDistributions) {
    jecFractionalUncertainty = (jetsCollection.JECUncertainty)->at(jetIndex);
    float jet_pT_JECDown = (1.0 - jecFractionalUncertainty)*jet_pT;
    float jet_pT_JECUp = (1.0 + jecFractionalUncertainty)*jet_pT;
    passesPT_JECDown = (jet_pT_JECDown > parameters.jetpTCut);
    passesPT_JECUp = (jet_pT_JECUp > parameters.jetpTCut);
    float jet_pT_missingHEMDown = jet_pT;
    float jet_pT_missingHEMUp = jet_pT;
    results.missing_HEM_adjustment_pT = 0.;
    if ((properties[jetProperty::phi] > -1.57) && (properties[jetProperty::phi] < -0.87)) {
      if ((properties[jetProperty::eta] > -2.5) && (properties[jetProperty::eta] <= -1.3)) {
	jet_pT_missingHEMDown = 0.8*jet_pT;
	jet_pT_missingHEMUp = 1.0*jet_pT;
	results.missing_HEM_adjustment_pT = 0.2*jet_pT;
      }
      else if ((properties[jetProperty::eta] > -3.0) && (properties[jetProperty::eta] <= -2.5)) {
	jet_pT_missingHEMDown = 0.65*jet_pT;
	jet_pT_missingHEMUp = 1.0*jet_pT;
	results.missing_HEM_adjustment_pT = 0.35*jet_pT;
      }
    }
    passesPT_missingHEMDown = (jet_pT_missingHEMDown > parameters.jetpTCut);
    passesPT_missingHEMUp = (jet_pT_missingHEMUp > parameters.jetpTCut);
  }
  results.jecFractionalUncertainty = jecFractionalUncertainty;

  if (parameters.calculatePrefiringWeights) {
    results.prefireWeights = findPrefireWeights((jetsCollection.eta)->at(jetIndex), jet_pT, parameters.prefiringEfficiencyMap);
  }
  else {
    results.prefireWeights = eventWeightsStruct(1.0f, 1.0f, 1.0f);
  }

  // ID cuts: PUID, jetID
  properties[jetProperty::PUID] = (jetsCollection.PUID)->at(jetIndex);
  bits[jetCriterion::puID] = (properties[jetProperty::PUID] > parameters.jetPUIDThreshold);

  properties[jetProperty::jetID] = (jetsCollection.ID)->at(jetIndex);
  bits[jetCriterion::jetID] = ((jetsCollection.ID)->at(jetIndex) == 6);

  int nNonPTFalseBits = getNFalseBits(bits);
  bool passesNonPTCriteria = (nNonPTFalseBits == 0);
  results.passesSelectionDRJECDown = (passesNonPTCriteria && passesPT_JECDown);
  results.passesSelectionDRJECUp = (passesNonPTCriteria && passesPT_JECUp);
  results.passesSelectionDRMissingHEMDown = (passesNonPTCriteria && passesPT_missingHEMDown);
  results.passesSelectionDRMissingHEMUp = (passesNonPTCriteria && passesPT_missingHEMUp);

  bool passesNonPTNonDRCriteria = ((nNonPTFalseBits == 0) || ((nNonPTFalseBits == 1) && (minDeltaR <= parameters.deltaRScale_jetPhotonDistance)));
  results.passesSelectionAllJECDown = (passesNonPTNonDRCriteria && passesPT_JECDown);
  results.passesSelectionAllJECUp = (passesNonPTNonDRCriteria && passesPT_JECUp);
  results.passesSelectionAllMissingHEMDown = (passesNonPTNonDRCriteria && passesPT_missingHEMDown);
  results.passesSelectionAllMissingHEMUp = (passesNonPTNonDRCriteria && passesPT_missingHEMUp);

  bits[jetCriterion::pT] = passesPT_nominal;
  assert(static_cast<int>(bits.size()) == static_cast<int>(jetCriterion::nJetCriteria));
  int nFalseBits = getNFalseBits(bits);
  results.passesSelectionDRJECNominal = (nFalseBits == 0);
  results.passesSelectionAllJECNominal = ((nFalseBits == 0) || ((nFalseBits == 1) && (minDeltaR <= parameters.deltaRScale_jetPhotonDistance)));
  results.isMarginallyUnselected = (nFalseBits == 1);
  if (results.isMarginallyUnselected) results.marginallyUnselectedCriterion = getFirstFalseCriterion(bits);

  results.contributesToHT = false;
  if (nFalseBits == 0) results.contributesToHT = true;
  else if ((nFalseBits == 1) && results.marginallyUnselectedCriterion == jetCriterion::deltaR_photon) results.contributesToHT = true;

  if (((results.isMarginallyUnselected || results.passesSelectionDRJECNominal) || (results.passesSelectionDRJECUp || results.passesSelectionDRJECDown)) || (results.passesSelectionDRMissingHEMDown || results.passesSelectionDRMissingHEMUp)) assert(static_cast<int>((results.jet_properties).size()) == static_cast<int>(jetProperty::nJetProperties));
  return results;
}

void setSelectedPhotonClosestJet(photonPropertiesCollection& photon_properties_collection, std::vector<angularVariablesStruct>& selectedJetAngles, std::vector<angularVariablesStruct>& genJetAngles, std::vector<angularVariablesStruct>& eventProgenitorMomGenJetAngles, std::vector<angularVariablesStruct>& singletMomGenJetAngles) {
  for (auto&& photon_properties: photon_properties_collection) {
    angularVariablesStruct photonAngle = angularVariablesStruct(photon_properties.at(photonProperty::eta), photon_properties.at(photonProperty::phi));
    photon_properties[photonProperty::deltaR_nearestSelectedJet] = photonAngle.getMinDeltaR(selectedJetAngles);
    photon_properties[photonProperty::deltaR_nearestGenJet] = photonAngle.getMinDeltaR(genJetAngles);
    photon_properties[photonProperty::deltaR_nearestEventProgenitorMomGenJet] = photonAngle.getMinDeltaR(eventProgenitorMomGenJetAngles);
    photon_properties[photonProperty::deltaR_nearestSingletMomGenJet] = photonAngle.getMinDeltaR(singletMomGenJetAngles);
    assert(static_cast<int>((photon_properties).size()) == static_cast<int>(photonProperty::nPhotonProperties));
  }
}

void setUnselectedMediumPhotonClosestJet(unselectedMediumPhotonPropertiesCollection& unselected_photon_properties_collection, std::vector<angularVariablesStruct>& selectedJetAngles, std::vector<angularVariablesStruct>& genJetAngles, std::vector<angularVariablesStruct>& eventProgenitorMomGenJetAngles, std::vector<angularVariablesStruct>& singletMomGenJetAngles) {
  for (auto&& unselected_photon_properties_pair: unselected_photon_properties_collection) {
    photonProperties& photon_properties = unselected_photon_properties_pair.second;
    angularVariablesStruct photonAngle = angularVariablesStruct(photon_properties.at(photonProperty::eta), photon_properties.at(photonProperty::phi));
    photon_properties[photonProperty::deltaR_nearestSelectedJet] = photonAngle.getMinDeltaR(selectedJetAngles);
    photon_properties[photonProperty::deltaR_nearestGenJet] = photonAngle.getMinDeltaR(genJetAngles);
    photon_properties[photonProperty::deltaR_nearestEventProgenitorMomGenJet] = photonAngle.getMinDeltaR(eventProgenitorMomGenJetAngles);
    photon_properties[photonProperty::deltaR_nearestSingletMomGenJet] = photonAngle.getMinDeltaR(singletMomGenJetAngles);
    assert(static_cast<int>((photon_properties).size()) == static_cast<int>(photonProperty::nPhotonProperties));
  }
}

void setUnselectedVetoedPhotonClosestJet(unselectedVetoedPhotonPropertiesCollection& unselected_photon_properties_collection, std::vector<angularVariablesStruct>& selectedJetAngles, std::vector<angularVariablesStruct>& genJetAngles, std::vector<angularVariablesStruct>& eventProgenitorMomGenJetAngles, std::vector<angularVariablesStruct>& singletMomGenJetAngles) {
  for (auto&& unselected_photon_properties_pair: unselected_photon_properties_collection) {
    photonProperties& photon_properties = unselected_photon_properties_pair.second;
    angularVariablesStruct photonAngle = angularVariablesStruct(photon_properties.at(photonProperty::eta), photon_properties.at(photonProperty::phi));
    photon_properties[photonProperty::deltaR_nearestSelectedJet] = photonAngle.getMinDeltaR(selectedJetAngles);
    photon_properties[photonProperty::deltaR_nearestGenJet] = photonAngle.getMinDeltaR(genJetAngles);
    photon_properties[photonProperty::deltaR_nearestEventProgenitorMomGenJet] = photonAngle.getMinDeltaR(eventProgenitorMomGenJetAngles);
    photon_properties[photonProperty::deltaR_nearestSingletMomGenJet] = photonAngle.getMinDeltaR(singletMomGenJetAngles);
    assert(static_cast<int>((photon_properties).size()) == static_cast<int>(photonProperty::nPhotonProperties));
  }
}

void setUnselectedFakePhotonClosestJet(unselectedFakePhotonPropertiesCollection& unselected_photon_properties_collection, std::vector<angularVariablesStruct>& selectedJetAngles, std::vector<angularVariablesStruct>& genJetAngles, std::vector<angularVariablesStruct>& eventProgenitorMomGenJetAngles, std::vector<angularVariablesStruct>& singletMomGenJetAngles) {
  for (auto&& unselected_photon_properties_pair: unselected_photon_properties_collection) {
    photonProperties& photon_properties = unselected_photon_properties_pair.second;
    angularVariablesStruct photonAngle = angularVariablesStruct(photon_properties.at(photonProperty::eta), photon_properties.at(photonProperty::phi));
    photon_properties[photonProperty::deltaR_nearestSelectedJet] = photonAngle.getMinDeltaR(selectedJetAngles);
    photon_properties[photonProperty::deltaR_nearestGenJet] = photonAngle.getMinDeltaR(genJetAngles);
    photon_properties[photonProperty::deltaR_nearestEventProgenitorMomGenJet] = photonAngle.getMinDeltaR(eventProgenitorMomGenJetAngles);
    photon_properties[photonProperty::deltaR_nearestSingletMomGenJet] = photonAngle.getMinDeltaR(singletMomGenJetAngles);
    assert(static_cast<int>((photon_properties).size()) == static_cast<int>(photonProperty::nPhotonProperties));
  }
}

bool checkPhotonsInBarrel(const std::vector<angularVariablesStruct> & selectedTruePhotonAngles, const float & eta_cut) {
  assert(selectedTruePhotonAngles.size() == 2);
  for (const angularVariablesStruct & photon_angle_info : selectedTruePhotonAngles) {
    if (std::fabs(photon_angle_info.eta) >= eta_cut) return false;
  }
  return true;
}

eventExaminationResultsStruct examineEvent(optionsStruct &options, parametersStruct &parameters, Long64_t& entryIndex, int & n_events_without_pu_info, // const int& year,
                                           eventDetailsStruct& eventDetails, MCCollectionStruct &MCCollection, photonsCollectionStruct &photonsCollection, jetsCollectionStruct &jetsCollection, statisticsHistograms& statistics, STRegionsStruct& STRegions) {
  eventExaminationResultsStruct eventResult;

  eventResult.eventIndex = entryIndex;
  selectionRegion& region = eventResult.evt_region;
  std::map<eventSelectionCriterion, bool> selectionBits;
  eventProperties event_properties;
  float& event_ST_electromagnetic = eventResult.evt_ST_electromagnetic;
  float& event_ST_hadronic = eventResult.evt_ST_hadronic;
  float& event_ST_MET = eventResult.evt_ST_MET;
  float& event_ST = eventResult.evt_ST;
  int n_goodJetsCloseToSelectedPhoton = 0;
  int& n_jetsDR = eventResult.evt_nJetsDR;
  int& n_jetsAll = eventResult.evt_nJetsAll;
  float& invariant_mass = eventResult.evt_invariantMass;
  std::map<shiftType, float>& shifted_ST = eventResult.evt_shifted_ST;
  std::map<shiftType, int>& shifted_nJetsDR = eventResult.evt_shifted_nJetsDR;
  std::map<shiftType, int>& shifted_nJetsAll = eventResult.evt_shifted_nJetsAll;

  // First save PU weights if needed
  eventResult.evt_PUWeight = -1.0;
  bool eventPUIsSensible = true;
  if (options.savePUWeights) {
    float event_PU_true = -1.;
    for (unsigned int BXCounter = 0; BXCounter < static_cast<unsigned int>((eventDetails.event_BX_for_PU)->size()); ++BXCounter) {
      int bx = (eventDetails.event_BX_for_PU)->at(BXCounter);
      if (bx == 0) {
	event_PU_true = (eventDetails.event_PU)->at(BXCounter);
	break;
      }
    }
    if ((event_PU_true < PU_MINVAL) || (event_PU_true > PU_MAXVAL)) {
      ++n_events_without_pu_info;
      if (static_cast<double>(1.0*entryIndex) > static_cast<double>(1.0/MAX_FRAC_EVENTS_WITHOUT_PU_INFO)) {
	assert(static_cast<double>(1.0*n_events_without_pu_info/entryIndex) < static_cast<double>(MAX_FRAC_EVENTS_WITHOUT_PU_INFO));
      }
      eventPUIsSensible = false;
    }
    else {
      eventResult.evt_PUWeight = static_cast<double>((parameters.PUWeights)->GetBinContent((parameters.PUWeights)->GetXaxis()->FindFixBin(event_PU_true)));
    }
  }

  // Additional selection, only for MC
  bool passesExtendedMCSelection = false;
  float generated_eventProgenitorMass = 0.;
  float generated_neutralinoMass = 0.;
  selectionBits[eventSelectionCriterion::MCGenInformation] = true;
  int nPhotonsWithDesiredMom = 0;
  truthPhotonPropertiesCollection selectedTruePhotonProperties;
  std::vector<angularVariablesStruct> selectedTruePhotonAngles;
  int nJetCandidatesWithStealthMom = 0;
  int nJetCandidatesWithEventProgenitorMom = 0;
  int nJetCandidatesWithSingletMom = 0;
  int nStealthJetsCloseToTruePhoton = 0;
  truthJetCandidatePropertiesCollection selectedTrueJetCandidateProperties_all;
  std::vector<angularVariablesStruct> selectedTrueJetCandidateAngles_all;
  truthJetCandidatePropertiesCollection selectedTrueJetCandidateProperties_fromEventProgenitor;
  std::vector<angularVariablesStruct> selectedTrueJetCandidateAngles_fromEventProgenitor;
  truthJetCandidatePropertiesCollection selectedTrueJetCandidateProperties_fromSinglet;
  std::vector<angularVariablesStruct> selectedTrueJetCandidateAngles_fromSinglet; // wasteful, fix later...
  int MCRegionIndex = 0;
  if ((options.enableMCEventFilter) && (!(options.MC_eventProgenitor == ""))) {
    bool eventProgenitorMassIsSet = false;
    bool neutralinoMassIsSet = false;
    for (int MCIndex = 0; MCIndex < eventDetails.nMCParticles; ++MCIndex) {
      // MCExaminationResultsStruct MCExaminationResults = examineMCParticle(options, parameters, MCCollection, MCIndex, doPromptOnlyBitMask);
      MCExaminationResultsStruct MCExaminationResults = examineMCParticle(options, parameters, MCCollection, MCIndex, false);
      if (MCExaminationResults.isPhotonWithDesiredMom) {
        ++nPhotonsWithDesiredMom;
        // min deltaR will be filled just outside the MC loop
        selectedTruePhotonProperties.push_back(MCExaminationResults.truth_photon_properties);
        selectedTruePhotonAngles.push_back(angularVariablesStruct((MCExaminationResults.truth_photon_properties)[truthPhotonProperty::eta], (MCExaminationResults.truth_photon_properties)[truthPhotonProperty::phi]));
      }
      if (MCExaminationResults.isJetCandidateFromStealthSource) {
        ++nJetCandidatesWithStealthMom;
        // min deltaR will be filled just outside the MC loop
        selectedTrueJetCandidateProperties_all.push_back(MCExaminationResults.truth_jetCandidate_properties);
        selectedTrueJetCandidateAngles_all.push_back(angularVariablesStruct((MCExaminationResults.truth_jetCandidate_properties)[truthJetCandidateProperty::eta], (MCExaminationResults.truth_jetCandidate_properties)[truthJetCandidateProperty::phi]));
        if (MCExaminationResults.isJetCandidateFromEventProgenitor) {
          ++nJetCandidatesWithEventProgenitorMom;
          selectedTrueJetCandidateProperties_fromEventProgenitor.push_back(MCExaminationResults.truth_jetCandidate_properties);
          selectedTrueJetCandidateAngles_fromEventProgenitor.push_back(angularVariablesStruct((MCExaminationResults.truth_jetCandidate_properties)[truthJetCandidateProperty::eta], (MCExaminationResults.truth_jetCandidate_properties)[truthJetCandidateProperty::phi]));
        }
        else if (MCExaminationResults.isJetCandidateFromSinglet) {
          ++nJetCandidatesWithSingletMom;
          selectedTrueJetCandidateProperties_fromSinglet.push_back(MCExaminationResults.truth_jetCandidate_properties);
          selectedTrueJetCandidateAngles_fromSinglet.push_back(angularVariablesStruct((MCExaminationResults.truth_jetCandidate_properties)[truthJetCandidateProperty::eta], (MCExaminationResults.truth_jetCandidate_properties)[truthJetCandidateProperty::phi]));
        }
      }
      if ((MCExaminationResults.eventProgenitorMass > 0.) && !(eventProgenitorMassIsSet)) {
        generated_eventProgenitorMass = MCExaminationResults.eventProgenitorMass;
        eventProgenitorMassIsSet = true;
      }
      if ((MCExaminationResults.neutralinoMass > 0.) && !(neutralinoMassIsSet)) {
        generated_neutralinoMass = MCExaminationResults.neutralinoMass;
        neutralinoMassIsSet = true;
      }
    } // ends loop over MC particles
    bool MCCriterion = ((eventProgenitorMassIsSet) && (nPhotonsWithDesiredMom == 2));
    selectionBits[eventSelectionCriterion::MCGenInformation] = MCCriterion;
    if (selectionBits[eventSelectionCriterion::MCGenInformation] && (!(eventProgenitorMassIsSet && neutralinoMassIsSet))) {
      std::cout << "ERROR: Unable to find eventProgenitor or neutralino mass in an event that passes MC selection." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    MCRegionIndex = MCRegions::getRegionIndex(generated_eventProgenitorMass, generated_neutralinoMass);
    passesExtendedMCSelection = MCCriterion && (checkPhotonsInBarrel(selectedTruePhotonAngles, parameters.photonBarrelEtaCut));
  }
  event_properties[eventProperty::MC_nPhotonsWithDesiredMom] = nPhotonsWithDesiredMom;
  event_properties[eventProperty::MC_nJetCandidatesWithStealthMom] = nJetCandidatesWithStealthMom;
  event_properties[eventProperty::MC_nJetCandidatesWithEventProgenitorMom] = nJetCandidatesWithEventProgenitorMom;
  event_properties[eventProperty::MC_nJetCandidatesWithSingletMom] = nJetCandidatesWithSingletMom;
  event_properties[eventProperty::MC_nStealthJetsCloseToTruePhoton] = nStealthJetsCloseToTruePhoton;

  // Photon selection maps. Index: photonIndex from the n-tuples
  std::map<int, photonType> selectedPhotonTypes;
  std::map<int, float> selectedPhotonPTs;
  std::map<int, float> selectedPhotonEtas;
  std::map<int, float> selectedPhotonMVAs;
  std::map<int, eventWeightsStruct> selectedPhotonScaleFactors;
  std::map<int, photonProperties> selectedPhotonProperties;
  std::map<int, angularVariablesStruct> selectedPhotonAngles;
  std::map<int, TLorentzVector> selectedPhotonFourMomenta;

  // Index maps are from the index within a selection to the photonIndex from the n-tuples
  // e.g. selectedPhotonPTs.at(selectedMediumPhotonIndices.at(i)) == (selectedMediumPhotonProperties.at(i)).at(photonProperty::pT) would be true for all i
  photonPropertiesCollection selectedMediumPhotonProperties;
  std::map<int, int> selectedMediumPhotonIndices;
  photonPropertiesCollection selectedMediumPhotonProperties_closeToTruePhoton;
  photonPropertiesCollection selectedMediumPhotonProperties_awayFromTruePhoton;
  photonPropertiesCollection selectedVetoedPhotonProperties;
  std::map<int, int> selectedVetoedPhotonIndices;
  photonPropertiesCollection selectedVetoedPhotonProperties_closeToTruePhoton;
  photonPropertiesCollection selectedVetoedPhotonProperties_awayFromTruePhoton;
  photonPropertiesCollection selectedFakePhotonProperties;
  std::map<int, int> selectedFakePhotonIndices;
  photonPropertiesCollection selectedFakePhotonProperties_closeToTruePhoton;
  photonPropertiesCollection selectedFakePhotonProperties_awayFromTruePhoton;

  unselectedMediumPhotonPropertiesCollection unselected_medium_pho_properties;
  unselectedMediumPhotonPropertiesCollection unselected_medium_pho_properties_closeToTruePhoton;
  unselectedMediumPhotonPropertiesCollection unselected_medium_pho_properties_awayFromTruePhoton;
  unselectedVetoedPhotonPropertiesCollection unselected_vetoed_pho_properties;
  unselectedVetoedPhotonPropertiesCollection unselected_vetoed_pho_properties_closeToTruePhoton;
  unselectedVetoedPhotonPropertiesCollection unselected_vetoed_pho_properties_awayFromTruePhoton;
  unselectedFakePhotonPropertiesCollection unselected_fake_pho_properties;
  unselectedFakePhotonPropertiesCollection unselected_fake_pho_properties_closeToTruePhoton;
  unselectedFakePhotonPropertiesCollection unselected_fake_pho_properties_awayFromTruePhoton;

  int n_mediumPhotons = 0;
  int n_mediumPhotonsPassingLeadingPTCut = 0;
  int n_truthMatchedMediumPhotons = 0;
  int n_vetoedPhotons = 0;
  int n_vetoedPhotonsPassingLeadingPTCut = 0;
  int n_truthMatchedVetoedPhotons = 0;
  int n_fakePhotons = 0;
  int n_fakePhotonsPassingLeadingPTCut = 0;
  int n_truthMatchedFakePhotons = 0;
  for (Int_t photonIndex = 0; photonIndex < (eventDetails.nPhotons); ++photonIndex) {
    photonExaminationResultsStruct photonExaminationResults = examinePhoton(options, parameters, (eventDetails.eventRho), photonsCollection, photonIndex, selectedTruePhotonAngles);
    if (photonExaminationResults.contributesToMisc2DHistograms) {
      statistics.fillMisc2DHistograms((photonExaminationResults.pho_properties)[photonProperty::rhoCorrectedChargedIsolation], (photonExaminationResults.pho_properties)[photonProperty::rhoCorrectedNeutralIsolation], (photonExaminationResults.pho_properties)[photonProperty::rhoCorrectedPhotonIsolation]);
    }
    if (photonExaminationResults.photon_type != photonType::nPhotonTypes) {
      selectedPhotonTypes[photonIndex] = photonExaminationResults.photon_type;
      // Make sure the photons are ordered in PT
      if (n_mediumPhotons > 0) assert(((photonExaminationResults.pho_properties)[photonProperty::pT]) <= selectedPhotonPTs.at(selectedMediumPhotonIndices.at(n_mediumPhotons-1)));
      if (n_vetoedPhotons > 0) assert(((photonExaminationResults.pho_properties)[photonProperty::pT]) <= selectedPhotonPTs.at(selectedVetoedPhotonIndices.at(n_vetoedPhotons-1)));
      if (n_fakePhotons > 0) assert(((photonExaminationResults.pho_properties)[photonProperty::pT]) <= selectedPhotonPTs.at(selectedFakePhotonIndices.at(n_fakePhotons-1)));
      selectedPhotonPTs[photonIndex] = (photonExaminationResults.pho_properties)[photonProperty::pT];
      selectedPhotonEtas[photonIndex] = (photonExaminationResults.pho_properties)[photonProperty::eta];
      selectedPhotonMVAs[photonIndex] = (photonExaminationResults.pho_properties)[photonProperty::mva];
      selectedPhotonScaleFactors[photonIndex] = eventWeightsStruct((photonExaminationResults.MCScaleFactors).nominal, (photonExaminationResults.MCScaleFactors).down, (photonExaminationResults.MCScaleFactors).up);
      selectedPhotonProperties[photonIndex] = photonExaminationResults.pho_properties;
      selectedPhotonAngles[photonIndex] = angularVariablesStruct((photonExaminationResults.pho_properties)[photonProperty::eta], (photonExaminationResults.pho_properties)[photonProperty::phi]);
      TLorentzVector photonFourMomentum;
      photonFourMomentum.SetPtEtaPhiE((photonExaminationResults.pho_properties).at(photonProperty::pT), (photonExaminationResults.pho_properties)[photonProperty::eta], (photonExaminationResults.pho_properties)[photonProperty::phi], photonExaminationResults.energy);
      selectedPhotonFourMomenta[photonIndex] = photonFourMomentum;
    }

    if (photonExaminationResults.photon_type == photonType::medium) {
      selectedMediumPhotonIndices[n_mediumPhotons] = photonIndex;
      ++n_mediumPhotons; // serves as index of selectedMediumPhoton
      if ((photonExaminationResults.pho_properties).at(photonProperty::pT) > parameters.pTCutLeading) ++n_mediumPhotonsPassingLeadingPTCut;
      selectedMediumPhotonProperties.push_back(photonExaminationResults.pho_properties);
      if (options.saveMCObjects) {
        float nearestTruePhotonDeltaR = (photonExaminationResults.pho_properties)[photonProperty::deltaR_nearestTruePhoton];
        if (nearestTruePhotonDeltaR >= parameters.deltaRScale_truthMatching) selectedMediumPhotonProperties_awayFromTruePhoton.push_back(photonExaminationResults.pho_properties);
        else if (nearestTruePhotonDeltaR > 0.) {
          ++n_truthMatchedMediumPhotons;
          selectedMediumPhotonProperties_closeToTruePhoton.push_back(photonExaminationResults.pho_properties);
        }
      }
    }
    else if (photonExaminationResults.photon_type == photonType::vetoed) {
      selectedVetoedPhotonIndices[n_vetoedPhotons] = photonIndex;
      ++n_vetoedPhotons;
      if ((photonExaminationResults.pho_properties).at(photonProperty::pT) > parameters.pTCutLeading) ++n_vetoedPhotonsPassingLeadingPTCut;
      selectedVetoedPhotonProperties.push_back(photonExaminationResults.pho_properties);
      if (options.saveMCObjects) {
        float nearestTruePhotonDeltaR = (photonExaminationResults.pho_properties)[photonProperty::deltaR_nearestTruePhoton];
        if (nearestTruePhotonDeltaR >= parameters.deltaRScale_truthMatching) selectedVetoedPhotonProperties_awayFromTruePhoton.push_back(photonExaminationResults.pho_properties);
        else if (nearestTruePhotonDeltaR > 0.) {
          ++n_truthMatchedVetoedPhotons;
          selectedVetoedPhotonProperties_closeToTruePhoton.push_back(photonExaminationResults.pho_properties);
        }
      }
    }
    else if (photonExaminationResults.photon_type == photonType::fake) {
      selectedFakePhotonIndices[n_fakePhotons] = photonIndex;
      ++n_fakePhotons;
      if ((photonExaminationResults.pho_properties).at(photonProperty::pT) > parameters.pTCutLeading) ++n_fakePhotonsPassingLeadingPTCut;
      selectedFakePhotonProperties.push_back(photonExaminationResults.pho_properties);
      if (options.saveMCObjects) {
        float nearestTruePhotonDeltaR = (photonExaminationResults.pho_properties)[photonProperty::deltaR_nearestTruePhoton];
        if (nearestTruePhotonDeltaR >= parameters.deltaRScale_truthMatching) selectedFakePhotonProperties_awayFromTruePhoton.push_back(photonExaminationResults.pho_properties);
        else if (nearestTruePhotonDeltaR > 0.) {
          ++n_truthMatchedFakePhotons;
          selectedFakePhotonProperties_closeToTruePhoton.push_back(photonExaminationResults.pho_properties);
        }
      }
    }

    if (photonExaminationResults.isMarginallyUnselected) {
      if (photonExaminationResults.marginallyUnselectedMediumCriterion != mediumPhotonCriterion::nMediumPhotonCriteria) {
	unselected_medium_pho_properties.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedMediumCriterion, photonExaminationResults.pho_properties));
	if (options.saveMCObjects) {
	  float nearestTruePhotonDeltaR = (photonExaminationResults.pho_properties)[photonProperty::deltaR_nearestTruePhoton];
	  if (nearestTruePhotonDeltaR >= parameters.deltaRScale_truthMatching) unselected_medium_pho_properties_awayFromTruePhoton.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedMediumCriterion, photonExaminationResults.pho_properties));
	  else if (nearestTruePhotonDeltaR > 0.) unselected_medium_pho_properties_closeToTruePhoton.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedMediumCriterion, photonExaminationResults.pho_properties));
	}
      }
      if (photonExaminationResults.marginallyUnselectedVetoedCriterion != vetoedPhotonCriterion::nVetoedPhotonCriteria) {
	unselected_vetoed_pho_properties.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedVetoedCriterion, photonExaminationResults.pho_properties));
	if (options.saveMCObjects) {
	  float nearestTruePhotonDeltaR = (photonExaminationResults.pho_properties)[photonProperty::deltaR_nearestTruePhoton];
	  if (nearestTruePhotonDeltaR >= parameters.deltaRScale_truthMatching) unselected_vetoed_pho_properties_awayFromTruePhoton.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedVetoedCriterion, photonExaminationResults.pho_properties));
	  else if (nearestTruePhotonDeltaR > 0.) unselected_vetoed_pho_properties_closeToTruePhoton.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedVetoedCriterion, photonExaminationResults.pho_properties));
	}
      }
      if (photonExaminationResults.marginallyUnselectedFakeCriterion != fakePhotonCriterion::nFakePhotonCriteria) {
	unselected_fake_pho_properties.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedFakeCriterion, photonExaminationResults.pho_properties));
	if (options.saveMCObjects) {
	  float nearestTruePhotonDeltaR = (photonExaminationResults.pho_properties)[photonProperty::deltaR_nearestTruePhoton];
	  if (nearestTruePhotonDeltaR >= parameters.deltaRScale_truthMatching) unselected_fake_pho_properties_awayFromTruePhoton.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedFakeCriterion, photonExaminationResults.pho_properties));
	  else if (nearestTruePhotonDeltaR > 0.) unselected_fake_pho_properties_closeToTruePhoton.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedFakeCriterion, photonExaminationResults.pho_properties));
	}
      }
    }
  }
  assert(static_cast<int>(selectedMediumPhotonProperties.size()) == n_mediumPhotons);
  assert(static_cast<int>(selectedVetoedPhotonProperties.size()) == n_vetoedPhotons);
  assert(static_cast<int>(selectedFakePhotonProperties.size()) == n_fakePhotons);

  event_properties[eventProperty::MC_nTruthMatchedMediumPhotons] = n_truthMatchedMediumPhotons;
  event_properties[eventProperty::MC_nTruthMatchedVetoedPhotons] = n_truthMatchedVetoedPhotons;
  event_properties[eventProperty::MC_nTruthMatchedFakePhotons] = n_truthMatchedFakePhotons;

  selectionBits[eventSelectionCriterion::doublePhoton] = false;
  selectionRegionDetailsStruct selection_region_details = selectionRegionUtils::getSelectionRegion(options.doSinglePhotonSelection, n_mediumPhotons, n_mediumPhotonsPassingLeadingPTCut, selectedMediumPhotonIndices, n_vetoedPhotons, n_vetoedPhotonsPassingLeadingPTCut, selectedVetoedPhotonIndices, n_fakePhotons, n_fakePhotonsPassingLeadingPTCut, selectedFakePhotonIndices, selectedPhotonPTs);
  int index_leadingPhoton = -1;
  int index_subLeadingPhoton = -1;
  if (selection_region_details.selection_region != selectionRegion::nSelectionRegions) {
    selectionBits[eventSelectionCriterion::doublePhoton] = true;
    index_leadingPhoton = selection_region_details.indexLeadingPhoton;
    index_subLeadingPhoton = selection_region_details.indexSubLeadingPhoton;
    region = selection_region_details.selection_region;
  }
  if (options.disablePhotonSelection) selectionBits[eventSelectionCriterion::doublePhoton] = true;

  int type_leadingPhoton = -1;
  float pT_leadingPhoton = -1.;
  float eta_leadingPhoton = -1000.;
  float mva_leadingPhoton = -2.;
  eventWeightsStruct scaleFactors_leadingPhoton;
  photonProperties properties_leadingPhoton;
  angularVariablesStruct photonAngle_leadingPhoton = angularVariablesStruct(-999., -1.);
  TLorentzVector photonFourMomentum_leadingPhoton;
  int type_subLeadingPhoton = -1;
  float pT_subLeadingPhoton = -1.;
  float eta_subLeadingPhoton = -1000.;
  float mva_subLeadingPhoton = -2.;
  eventWeightsStruct scaleFactors_subLeadingPhoton;
  photonProperties properties_subLeadingPhoton;
  angularVariablesStruct photonAngle_subLeadingPhoton = angularVariablesStruct(-999., -1.);
  TLorentzVector photonFourMomentum_subLeadingPhoton;
  std::vector<angularVariablesStruct> list_selectedPhotonAngles;
  std::vector<float> list_selectedPhotonPTs;
  std::vector<angularVariablesStruct> list_selectedFakePhotonAngles;
  std::vector<angularVariablesStruct> list_selectedLoosePhotonAngles;
  std::vector<angularVariablesStruct> list_selectedNominalPhotonAngles;
  std::vector<TLorentzVector> list_selectedPhotonFourMomenta;
  if (region != selectionRegion::nSelectionRegions) {
    assert(index_leadingPhoton >= 0);
    assert(((options.doSinglePhotonSelection) && (index_subLeadingPhoton == -1)) || (index_subLeadingPhoton >= 0));
    assert(index_leadingPhoton != index_subLeadingPhoton);
    type_leadingPhoton = static_cast<int>(selectedPhotonTypes.at(index_leadingPhoton));
    pT_leadingPhoton = selectedPhotonPTs.at(index_leadingPhoton);
    eta_leadingPhoton = selectedPhotonEtas.at(index_leadingPhoton);
    mva_leadingPhoton = selectedPhotonMVAs.at(index_leadingPhoton);
    scaleFactors_leadingPhoton = eventWeightsStruct((selectedPhotonScaleFactors.at(index_leadingPhoton)).nominal, (selectedPhotonScaleFactors.at(index_leadingPhoton)).down, (selectedPhotonScaleFactors.at(index_leadingPhoton)).up);
    properties_leadingPhoton = selectedPhotonProperties.at(index_leadingPhoton);
    photonAngle_leadingPhoton = angularVariablesStruct((selectedPhotonAngles.at(index_leadingPhoton)).eta, (selectedPhotonAngles.at(index_leadingPhoton)).phi);
    list_selectedPhotonAngles.push_back(photonAngle_leadingPhoton);
    list_selectedPhotonPTs.push_back(pT_leadingPhoton);
    if (type_leadingPhoton == static_cast<int>(photonType::fake)) list_selectedFakePhotonAngles.push_back(photonAngle_leadingPhoton);
    else if (type_leadingPhoton == static_cast<int>(photonType::vetoed)) list_selectedLoosePhotonAngles.push_back(photonAngle_leadingPhoton);
    else if (type_leadingPhoton == static_cast<int>(photonType::medium)) list_selectedNominalPhotonAngles.push_back(photonAngle_leadingPhoton);
    list_selectedPhotonFourMomenta.push_back(selectedPhotonFourMomenta.at(index_leadingPhoton));
    if (!(options.doSinglePhotonSelection)) {
      type_subLeadingPhoton = static_cast<int>(selectedPhotonTypes.at(index_subLeadingPhoton));
      pT_subLeadingPhoton = selectedPhotonPTs.at(index_subLeadingPhoton);
      eta_subLeadingPhoton = selectedPhotonEtas.at(index_subLeadingPhoton);
      mva_subLeadingPhoton = selectedPhotonMVAs.at(index_subLeadingPhoton);
      scaleFactors_subLeadingPhoton = eventWeightsStruct((selectedPhotonScaleFactors.at(index_subLeadingPhoton)).nominal, (selectedPhotonScaleFactors.at(index_subLeadingPhoton)).down, (selectedPhotonScaleFactors.at(index_subLeadingPhoton)).up);
      properties_subLeadingPhoton = selectedPhotonProperties.at(index_subLeadingPhoton);
      photonAngle_subLeadingPhoton = angularVariablesStruct((selectedPhotonAngles.at(index_subLeadingPhoton)).eta, (selectedPhotonAngles.at(index_subLeadingPhoton)).phi);
      list_selectedPhotonAngles.push_back(photonAngle_subLeadingPhoton);
      list_selectedPhotonPTs.push_back(pT_subLeadingPhoton);
      if (type_leadingPhoton == static_cast<int>(photonType::fake)) list_selectedFakePhotonAngles.push_back(photonAngle_subLeadingPhoton);
      else if (type_leadingPhoton == static_cast<int>(photonType::vetoed)) list_selectedLoosePhotonAngles.push_back(photonAngle_subLeadingPhoton);
      else if (type_leadingPhoton == static_cast<int>(photonType::medium)) list_selectedNominalPhotonAngles.push_back(photonAngle_subLeadingPhoton);
      list_selectedPhotonFourMomenta.push_back(selectedPhotonFourMomenta.at(index_subLeadingPhoton));
    }
    assert(pT_leadingPhoton >= pT_subLeadingPhoton);

    event_ST_electromagnetic += pT_leadingPhoton;
    if (!(options.doSinglePhotonSelection)) event_ST_electromagnetic += pT_subLeadingPhoton;
    event_ST += event_ST_electromagnetic;
    if (options.calculateShiftedDistributions) {
      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
	shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
	// no effect on photon's contribution to ST due to any of the shifts
	addShiftedEToSTMap(pT_leadingPhoton, shifted_ST, typeIndex);
	addShiftedEToSTMap(pT_subLeadingPhoton, shifted_ST, typeIndex);
      }
    }
    if (options.calculateMCScaleFactorWeights) {
      (eventResult.evt_photonMCScaleFactors).nominal *= ((scaleFactors_leadingPhoton.nominal)*(scaleFactors_subLeadingPhoton.nominal));
      (eventResult.evt_photonMCScaleFactors).down *= ((scaleFactors_leadingPhoton.down)*(scaleFactors_subLeadingPhoton.down));
      (eventResult.evt_photonMCScaleFactors).up *= ((scaleFactors_leadingPhoton.up)*(scaleFactors_subLeadingPhoton.up));
    }
    if (list_selectedPhotonAngles.size() == 2) {
      eventResult.evt_deltaR_photons = (list_selectedPhotonAngles.at(0)).get_deltaR(list_selectedPhotonAngles.at(1));
    }
  }

  selectionBits[eventSelectionCriterion::invariantMass] = true;
  if ((region != selectionRegion::nSelectionRegions) && (!(options.doSinglePhotonSelection))) {
    invariant_mass = getDiphotonInvariantMass(list_selectedPhotonFourMomenta);
    selectionBits[eventSelectionCriterion::invariantMass] = (invariant_mass >= parameters.invariantMassCut);
  }
  if (options.disablePhotonSelection) selectionBits[eventSelectionCriterion::invariantMass] = true;

  event_properties[eventProperty::nMediumPhotons] = 1.0*n_mediumPhotons;
  event_properties[eventProperty::nVetoedPhotons] = 1.0*n_vetoedPhotons;
  event_properties[eventProperty::nFakePhotons] = 1.0*n_fakePhotons;
  event_properties[eventProperty::leadingPhotonPT] = pT_leadingPhoton;
  event_properties[eventProperty::subLeadingPhotonPT] = pT_subLeadingPhoton;
  event_properties[eventProperty::leadingPhotonEta] = eta_leadingPhoton;
  event_properties[eventProperty::subLeadingPhotonEta] = eta_subLeadingPhoton;
  event_properties[eventProperty::leadingPhotonType] = type_leadingPhoton;
  event_properties[eventProperty::subLeadingPhotonType] = type_subLeadingPhoton;
  event_properties[eventProperty::invariantMass] = invariant_mass;

  // The following reads in MC collections but requires photons to already have been reconstructed
  GenLevelEventInfoStruct& gen_level_info = eventResult.evt_gen_level_info;
  gen_level_info.nKinematicPhotons = 0;
  gen_level_info.nRecoPhotonsMatchedToGenPhotons = 0;
  gen_level_info.deltaR_genPhoton_mom_matchingLeadingPhoton = -0.5;
  gen_level_info.deltaR_genPhoton_mom_matchingSubLeadingPhoton = -0.5;
  std::vector<angularVariablesStruct> finalStateGenLevelKinematicPhotonAngles;
  std::vector<float> finalStateGenLevelKinematicPhotonPTs;
  std::vector<float> finalStateGenLevelKinematicPhotonDeltaRWRTMoms;
  std::vector<angularVariablesStruct> promptFSGenLevelPhotonAngles;
  if (options.saveMCGenLevelInfo) {
    for (int MCIndex = 0; MCIndex < eventDetails.nMCParticles; ++MCIndex) {
      int particle_mcPID = (MCCollection.MCPIDs)->at(MCIndex);
      if (PIDUtils::isPhotonPID(particle_mcPID)) {
	UShort_t particle_statusFlag = static_cast<UShort_t>(((MCCollection.MCStatusFlags)->at(MCIndex))&(static_cast<UShort_t>(7u))); // picks out only first 3 bits
	bool passes_bit_mask = passesBitMask(particle_statusFlag, parameters.MCStatusFlagBitMask_promptOnly);
	float gen_photon_eta = (MCCollection.MCEtas)->at(MCIndex);
	float gen_photon_phi = (MCCollection.MCPhis)->at(MCIndex);
	float gen_photon_et = (MCCollection.MCEts)->at(MCIndex);
	if ((std::fabs(gen_photon_eta) < parameters.photonBarrelEtaCut) &&
	    (std::fabs(gen_photon_et) > parameters.pTCutSubLeading) &&
	    passes_bit_mask) {
	  ++(gen_level_info.nKinematicPhotons);
	  angularVariablesStruct gen_photon_angle = angularVariablesStruct(gen_photon_eta, gen_photon_phi);
	  int mc_mom_pid = -1;
	  if (options.saveMCMomInfo) mc_mom_pid = (MCCollection.MCMomPIDs)->at(MCIndex);
	  float deltaR_wrt_mom = -0.5;
	  if (mc_mom_pid > 0) {
	    angularVariablesStruct gen_photon_mom_angle = angularVariablesStruct((MCCollection.MCMomEtas)->at(MCIndex), (MCCollection.MCMomPhis)->at(MCIndex));
	    deltaR_wrt_mom = gen_photon_angle.get_deltaR(gen_photon_mom_angle);
	  }
	  finalStateGenLevelKinematicPhotonAngles.push_back(gen_photon_angle);
	  finalStateGenLevelKinematicPhotonPTs.push_back(gen_photon_et);
	  finalStateGenLevelKinematicPhotonDeltaRWRTMoms.push_back(deltaR_wrt_mom);
	  if (options.doOverlapRemoval) {
	    promptFSGenLevelPhotonAngles.push_back(gen_photon_angle);
	  }
	}
      }
    }
    assert(finalStateGenLevelKinematicPhotonAngles.size() == finalStateGenLevelKinematicPhotonPTs.size());
    assert(finalStateGenLevelKinematicPhotonDeltaRWRTMoms.size() == finalStateGenLevelKinematicPhotonPTs.size());
    if (list_selectedPhotonAngles.size() > 0) {
      assert(list_selectedPhotonAngles.size() == list_selectedPhotonPTs.size());
      bool leadingMatchedGenMomDeltaRIsSet = false;
      bool subLeadingMatchedGenMomDeltaRIsSet = false;
      for (size_t selectedPhotonIndex = 0; selectedPhotonIndex < list_selectedPhotonAngles.size(); ++selectedPhotonIndex) {
	std::pair<int, float> index_and_min_dr = (list_selectedPhotonAngles.at(selectedPhotonIndex)).getIndexAndMinDeltaR(finalStateGenLevelKinematicPhotonAngles);
	int min_dr_index = index_and_min_dr.first;
	float min_dr = index_and_min_dr.second;
	if ((min_dr > 0.) && (min_dr <= parameters.deltaRScale_truthMatching)) {
	  assert(min_dr_index >= 0);
	  float matched_gen_pt = finalStateGenLevelKinematicPhotonPTs.at(min_dr_index);
	  if (std::fabs((matched_gen_pt/(list_selectedPhotonPTs.at(selectedPhotonIndex))) - 1.0) < parameters.pTRatioMinusOneThreshold_truthMatching) {
	    ++(gen_level_info.nRecoPhotonsMatchedToGenPhotons);
	    if (!(leadingMatchedGenMomDeltaRIsSet)) {
	      (gen_level_info.deltaR_genPhoton_mom_matchingLeadingPhoton) = finalStateGenLevelKinematicPhotonDeltaRWRTMoms.at(min_dr_index);
	      leadingMatchedGenMomDeltaRIsSet = true;
	    }
	    else if (!(subLeadingMatchedGenMomDeltaRIsSet)) {
	      (gen_level_info.deltaR_genPhoton_mom_matchingSubLeadingPhoton) = finalStateGenLevelKinematicPhotonDeltaRWRTMoms.at(min_dr_index);
	      subLeadingMatchedGenMomDeltaRIsSet = true;
	    }
	  }
	}
      }
    }
    if (options.doOverlapRemoval) {
      eventResult.evt_nPhotonsMatchedToGenPromptFS = 0;
      for (angularVariablesStruct & selected_photon_angle : list_selectedPhotonAngles) {
	float min_deltaR_wrt_promptFSGenLevelPhotons = selected_photon_angle.getMinDeltaR(promptFSGenLevelPhotonAngles);
	if ((min_deltaR_wrt_promptFSGenLevelPhotons > 0.) && (min_deltaR_wrt_promptFSGenLevelPhotons < parameters.deltaRScale_truthMatching)) {
	  ++(eventResult.evt_nPhotonsMatchedToGenPromptFS);
	}
      }
    }
  }

  selectionBits[eventSelectionCriterion::overlap] = true;
  if (options.doOverlapRemoval) {
    selectionBits[eventSelectionCriterion::overlap] = (eventResult.evt_nPhotonsMatchedToGenPromptFS <= options.overlapRemoval_maxNPromptPhotons);
  }

  // Jet selection
  float evt_hT = 0;
  float pT_leadingJet = -1.;
  jetPropertiesCollection selectedJetProperties;
  jetPropertiesCollection selectedJetProperties_closeToTruePhoton;
  jetPropertiesCollection selectedJetProperties_awayFromTruePhoton;
  unselectedJetPropertiesCollection unselected_jet_properties;
  unselectedJetPropertiesCollection unselected_jet_properties_closeToTruePhoton;
  unselectedJetPropertiesCollection unselected_jet_properties_awayFromTruePhoton;
  genJetPropertiesCollection gen_jet_properties_collection;
  genJetPropertiesCollection eventProgenitor_mom_gen_jet_properties_collection;
  genJetPropertiesCollection singlet_mom_gen_jet_properties_collection;
  std::vector<angularVariablesStruct> selectedJetAngles;
  std::vector<angularVariablesStruct> genJetAngles;
  std::vector<angularVariablesStruct> eventProgenitorMomGenJetAngles;
  std::vector<angularVariablesStruct> singletMomGenJetAngles;
  int n_genJets = 0;
  int n_eventProgenitorMomGenJets = 0;
  int n_singletMomGenJets = 0;
  float MET_HEMAdjustmentX = 0.;
  float MET_HEMAdjustmentY = 0.;

  for (Int_t jetIndex = 0; jetIndex < (eventDetails.nJets); ++jetIndex) {
    jetExaminationResultsStruct jetExaminationResults = examineJet(options, parameters, jetsCollection, jetIndex, list_selectedPhotonAngles, selectedTruePhotonAngles, selectedTrueJetCandidateAngles_all);
    if (options.saveMCObjects && jetExaminationResults.hasGenVariablesSet) {
      if ((jetExaminationResults.gen_jet_properties)[genJetProperty::pT] > parameters.pTCutSubLeading) {
        ++n_genJets;
        gen_jet_properties_collection.push_back(jetExaminationResults.gen_jet_properties);
        genJetAngles.push_back(angularVariablesStruct((jetExaminationResults.gen_jet_properties)[genJetProperty::eta], (jetExaminationResults.gen_jet_properties)[genJetProperty::phi]));
        if (jetExaminationResults.hasEventProgenitorPartonMom) {
          ++n_eventProgenitorMomGenJets;
          eventProgenitor_mom_gen_jet_properties_collection.push_back(jetExaminationResults.gen_jet_properties);
          eventProgenitorMomGenJetAngles.push_back(angularVariablesStruct((jetExaminationResults.gen_jet_properties)[genJetProperty::eta], (jetExaminationResults.gen_jet_properties)[genJetProperty::phi]));
        }
        if (jetExaminationResults.hasSingletPartonMom) {
          ++n_singletMomGenJets;
          singlet_mom_gen_jet_properties_collection.push_back(jetExaminationResults.gen_jet_properties);
          singletMomGenJetAngles.push_back(angularVariablesStruct((jetExaminationResults.gen_jet_properties)[genJetProperty::eta], (jetExaminationResults.gen_jet_properties)[genJetProperty::phi]));
        }
      }
    }
    (eventResult.evt_prefireWeights).nominal *= (jetExaminationResults.prefireWeights).nominal; // All jets, whether or not they pass any of the cuts, contribute to the prefiring weight
    (eventResult.evt_prefireWeights).down *= (jetExaminationResults.prefireWeights).down;
    (eventResult.evt_prefireWeights).up *= (jetExaminationResults.prefireWeights).up;
    if (jetExaminationResults.passesSelectionDRJECNominal) selectedJetProperties.push_back(jetExaminationResults.jet_properties);
    else if (jetExaminationResults.isMarginallyUnselected) {
      unselected_jet_properties.push_back(std::make_pair(jetExaminationResults.marginallyUnselectedCriterion, jetExaminationResults.jet_properties));
      if (options.saveMCObjects) {
        float nearestTruePhotonDeltaR = (jetExaminationResults.jet_properties)[jetProperty::deltaR_nearestTruePhoton];
        if (nearestTruePhotonDeltaR >= parameters.deltaRScale_truthMatching) unselected_jet_properties_awayFromTruePhoton.push_back(std::make_pair(jetExaminationResults.marginallyUnselectedCriterion, jetExaminationResults.jet_properties));
        else if (nearestTruePhotonDeltaR > 0.) unselected_jet_properties_closeToTruePhoton.push_back(std::make_pair(jetExaminationResults.marginallyUnselectedCriterion, jetExaminationResults.jet_properties));
      }
    }
    if (jetExaminationResults.contributesToHT) {
      evt_hT += jetExaminationResults.jet_properties[jetProperty::pT]; // Add to hT whether or not jet passes deltaR check
      ++n_jetsAll; // n_jetsAll counts all jets that pass the other selections, whether or not they are close to a photon
      selectedJetAngles.push_back(angularVariablesStruct((jetExaminationResults.jet_properties)[jetProperty::eta], (jetExaminationResults.jet_properties)[jetProperty::phi]));
      angularVariablesStruct jetAngle = angularVariablesStruct((jetExaminationResults.jet_properties)[jetProperty::eta], (jetExaminationResults.jet_properties)[jetProperty::phi]);
      float nearestFakePhotonDeltaR = jetAngle.getMinDeltaR(list_selectedFakePhotonAngles);
      float nearestLoosePhotonDeltaR = jetAngle.getMinDeltaR(list_selectedLoosePhotonAngles);
      float nearestNominalPhotonDeltaR = jetAngle.getMinDeltaR(list_selectedNominalPhotonAngles);
      statistics.fillRecoEfficiencyStatisticsHistograms(jetExaminationResults.isCloseToTruePhoton, jetExaminationResults.jet_properties[jetProperty::pT],
                                             ((nearestFakePhotonDeltaR > 0.) && (nearestFakePhotonDeltaR <= parameters.deltaRScale_jetPhotonDistance)),
                                             ((nearestLoosePhotonDeltaR > 0.) && (nearestLoosePhotonDeltaR <= parameters.deltaRScale_jetPhotonDistance)),
                                             ((nearestNominalPhotonDeltaR > 0.) && (nearestNominalPhotonDeltaR <= parameters.deltaRScale_jetPhotonDistance)));

      if (jetExaminationResults.passesSelectionDRJECNominal) {
	event_ST_hadronic += jetExaminationResults.jet_properties[jetProperty::pT];
        event_ST += jetExaminationResults.jet_properties[jetProperty::pT]; // Add to sT only if jet passes deltaR check, to avoid double-counting
        ++n_jetsDR; // n_jetsDR counts only those jets that are sufficiently away from a photon
        if ((pT_leadingJet < 0.) || (pT_leadingJet < jetExaminationResults.jet_properties[jetProperty::pT])) pT_leadingJet = jetExaminationResults.jet_properties[jetProperty::pT];
        selectedJetProperties.push_back(jetExaminationResults.jet_properties);
        if (options.saveMCObjects) {
          float nearestTruePhotonDeltaR = (jetExaminationResults.jet_properties)[jetProperty::deltaR_nearestTruePhoton];
          if (nearestTruePhotonDeltaR >= parameters.deltaRScale_truthMatching) selectedJetProperties_awayFromTruePhoton.push_back(jetExaminationResults.jet_properties);
          else if (nearestTruePhotonDeltaR > 0.) selectedJetProperties_closeToTruePhoton.push_back(jetExaminationResults.jet_properties);
        }
      }
      else {
        ++n_goodJetsCloseToSelectedPhoton;
      }
    }

    if (options.calculateShiftedDistributions && (((jetExaminationResults.passesSelectionDRJECDown || jetExaminationResults.passesSelectionDRJECUp) || jetExaminationResults.passesSelectionDRJECNominal) || (jetExaminationResults.passesSelectionDRMissingHEMDown || jetExaminationResults.passesSelectionDRMissingHEMUp))) { // Actually we just need to check JECDown and MissingHEMDown...
      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
        shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);

        bool passes_selection = jetExaminationResults.passesSelectionDRJECNominal;
        float shifted_contribution = (jetExaminationResults.jet_properties[jetProperty::pT]);
        if (typeIndex == shiftType::JECDown) {
          passes_selection = jetExaminationResults.passesSelectionDRJECDown;
          shifted_contribution = (1.0 - (jetExaminationResults.jecFractionalUncertainty))*(jetExaminationResults.jet_properties[jetProperty::pT]);
        }
        else if (typeIndex == shiftType::JECUp) {
          passes_selection = jetExaminationResults.passesSelectionDRJECUp;
          shifted_contribution = (1.0 + (jetExaminationResults.jecFractionalUncertainty))*(jetExaminationResults.jet_properties[jetProperty::pT]);
        }
	else if (typeIndex == shiftType::missingHEMDown) {
          passes_selection = jetExaminationResults.passesSelectionDRMissingHEMDown;
          shifted_contribution = jetExaminationResults.jet_properties[jetProperty::pT] - jetExaminationResults.missing_HEM_adjustment_pT;
        }
	else if (typeIndex == shiftType::missingHEMUp) {
          passes_selection = jetExaminationResults.passesSelectionDRMissingHEMUp;
          shifted_contribution = 1.0*(jetExaminationResults.jet_properties[jetProperty::pT]);
        }

        if (passes_selection) {
          addShiftedEToSTMap(shifted_contribution, shifted_ST, typeIndex);
          incrementNJetsMap(shifted_nJetsDR, typeIndex);
        }
      }
      if (jetExaminationResults.missing_HEM_adjustment_pT > 0) {
        MET_HEMAdjustmentX += (-1.0)*(jetExaminationResults.missing_HEM_adjustment_pT)*std::cos(jetExaminationResults.jet_properties[jetProperty::phi]);
        MET_HEMAdjustmentY += (-1.0)*(jetExaminationResults.missing_HEM_adjustment_pT)*std::sin(jetExaminationResults.jet_properties[jetProperty::phi]);
      }
    }
    if (options.calculateShiftedDistributions && (((jetExaminationResults.passesSelectionAllJECDown || jetExaminationResults.passesSelectionAllJECUp) || jetExaminationResults.passesSelectionAllJECNominal) || (jetExaminationResults.passesSelectionAllMissingHEMDown || jetExaminationResults.passesSelectionAllMissingHEMUp))) { // Actually we just need to check JECDown and MissingHEMDown...
      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
        shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);

        bool passes_selection = jetExaminationResults.passesSelectionAllJECNominal;
        if (typeIndex == shiftType::JECDown) {
          passes_selection = jetExaminationResults.passesSelectionAllJECDown;
        }
        else if (typeIndex == shiftType::JECUp) {
          passes_selection = jetExaminationResults.passesSelectionAllJECUp;
        }
	else if (typeIndex == shiftType::missingHEMDown) {
          passes_selection = jetExaminationResults.passesSelectionAllMissingHEMDown;
        }
	else if (typeIndex == shiftType::missingHEMUp) {
          passes_selection = jetExaminationResults.passesSelectionAllMissingHEMUp;
        }

        if (passes_selection) {
          // addShiftedEToSTMap(shifted_contribution, shifted_ST, typeIndex); // only jets passing deltaR criterion contribute to shifted_ST
          incrementNJetsMap(shifted_nJetsAll, typeIndex);
        }
      }
      // if (jetExaminationResults.missing_HEM_adjustment_pT > 0) {
      //   MET_HEMAdjustmentX += (-1.0)*(jetExaminationResults.missing_HEM_adjustment_pT)*std::cos(jetExaminationResults.jet_properties[jetProperty::phi]);
      //   MET_HEMAdjustmentY += (-1.0)*(jetExaminationResults.missing_HEM_adjustment_pT)*std::sin(jetExaminationResults.jet_properties[jetProperty::phi]);
      // }
    }
  }
  event_properties[eventProperty::hT] = evt_hT;
  event_properties[eventProperty::nGoodJetsCloseToSelectedPhoton] = n_goodJetsCloseToSelectedPhoton;
  event_properties[eventProperty::nJetsDR] = n_jetsDR;
  event_properties[eventProperty::nJetsAll] = n_jetsAll;
  int max_nJets = n_jetsDR;
  if (options.calculateShiftedDistributions) { // this makes sure that the nJets used to make the decision whether or not to save the event is the maximum nJets accounting for all the shifts
    int maxNJetsShifted = getMaxNJets(shifted_nJetsDR);
    if (maxNJetsShifted > max_nJets) max_nJets = maxNJetsShifted;
  }

  selectionBits[eventSelectionCriterion::HLTSelection] = true;
  if (options.doHLTSelection) {
    assert((parameters.HLTBit_photon >= 0) || (parameters.HLTBit_jet >= 0));
    if (parameters.HLT_triggerType == triggerType::jet) {
      selectionBits[eventSelectionCriterion::HLTSelection] = checkHLTBit(eventDetails.HLTJetBits, parameters.HLTBit_jet);
    }
    else if (parameters.HLT_triggerType == triggerType::photon) {
      selectionBits[eventSelectionCriterion::HLTSelection] = checkHLTBit(eventDetails.HLTPhotonBits, parameters.HLTBit_photon);
    }
    else {
      std::cout << "ERROR: parameter \"HLT_triggerType\" is neither \"photon\" nor \"jet\"." << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }

  event_properties[eventProperty::MC_nGenJets] = n_genJets;
  event_properties[eventProperty::MC_nEventProgenitorMomGenJets] = n_eventProgenitorMomGenJets;
  event_properties[eventProperty::MC_nSingletMomGenJets] = n_singletMomGenJets;
  if (options.saveMCObjects) {
    assert(static_cast<int>(selectedTruePhotonProperties.size()) == nPhotonsWithDesiredMom);
    assert(static_cast<int>(selectedTruePhotonAngles.size()) == nPhotonsWithDesiredMom);
    for (int truePhotonCounter = 0; truePhotonCounter < nPhotonsWithDesiredMom; ++truePhotonCounter) {
      (selectedTruePhotonProperties[truePhotonCounter])[truthPhotonProperty::deltaR_nearestTruthJetCandidate] = (selectedTruePhotonAngles[truePhotonCounter]).getMinDeltaR(selectedTrueJetCandidateAngles_all);
      (selectedTruePhotonProperties[truePhotonCounter])[truthPhotonProperty::deltaR_nearestGenJet] = (selectedTruePhotonAngles[truePhotonCounter]).getMinDeltaR(genJetAngles);
      (selectedTruePhotonProperties[truePhotonCounter])[truthPhotonProperty::deltaR_nearestEventProgenitorMomGenJet] = (selectedTruePhotonAngles[truePhotonCounter]).getMinDeltaR(eventProgenitorMomGenJetAngles);
      (selectedTruePhotonProperties[truePhotonCounter])[truthPhotonProperty::deltaR_nearestSingletMomGenJet] = (selectedTruePhotonAngles[truePhotonCounter]).getMinDeltaR(singletMomGenJetAngles);
      assert(static_cast<int>((selectedTruePhotonProperties[truePhotonCounter]).size()) == static_cast<int>(truthPhotonProperty::nTruthPhotonProperties));
    }

    assert(static_cast<int>(selectedTrueJetCandidateProperties_all.size()) == nJetCandidatesWithStealthMom);
    assert(static_cast<int>(selectedTrueJetCandidateAngles_all.size()) == nJetCandidatesWithStealthMom);
    for (int trueJetCandidateCounter = 0; trueJetCandidateCounter < nJetCandidatesWithStealthMom; ++trueJetCandidateCounter) {
      float min_deltaR = (selectedTrueJetCandidateAngles_all[trueJetCandidateCounter]).getMinDeltaR(selectedTruePhotonAngles);
      if (min_deltaR < parameters.deltaRScale_jetPhotonDistance) ++nStealthJetsCloseToTruePhoton;
      (selectedTrueJetCandidateProperties_all[trueJetCandidateCounter])[truthJetCandidateProperty::deltaR_nearestTruePhoton] = min_deltaR;
      assert(static_cast<int>((selectedTrueJetCandidateProperties_all[trueJetCandidateCounter]).size()) == static_cast<int>(truthJetCandidateProperty::nTruthJetCandidateProperties));
    }
    assert(static_cast<int>(selectedTrueJetCandidateProperties_fromEventProgenitor.size()) == nJetCandidatesWithEventProgenitorMom);
    assert(static_cast<int>(selectedTrueJetCandidateAngles_fromEventProgenitor.size()) == nJetCandidatesWithEventProgenitorMom);
    for (int trueJetCandidateCounter = 0; trueJetCandidateCounter < nJetCandidatesWithEventProgenitorMom; ++trueJetCandidateCounter) {
      (selectedTrueJetCandidateProperties_fromEventProgenitor[trueJetCandidateCounter])[truthJetCandidateProperty::deltaR_nearestTruePhoton] = (selectedTrueJetCandidateAngles_fromEventProgenitor[trueJetCandidateCounter]).getMinDeltaR(selectedTruePhotonAngles);
      assert(static_cast<int>((selectedTrueJetCandidateProperties_fromEventProgenitor[trueJetCandidateCounter]).size()) == static_cast<int>(truthJetCandidateProperty::nTruthJetCandidateProperties));
    }
    assert(static_cast<int>(selectedTrueJetCandidateProperties_fromSinglet.size()) == nJetCandidatesWithSingletMom);
    assert(static_cast<int>(selectedTrueJetCandidateAngles_fromSinglet.size()) == nJetCandidatesWithSingletMom);
    for (int trueJetCandidateCounter = 0; trueJetCandidateCounter < nJetCandidatesWithSingletMom; ++trueJetCandidateCounter) {
      (selectedTrueJetCandidateProperties_fromSinglet[trueJetCandidateCounter])[truthJetCandidateProperty::deltaR_nearestTruePhoton] = (selectedTrueJetCandidateAngles_fromSinglet[trueJetCandidateCounter]).getMinDeltaR(selectedTruePhotonAngles);
      assert(static_cast<int>((selectedTrueJetCandidateProperties_fromSinglet[trueJetCandidateCounter]).size()) == static_cast<int>(truthJetCandidateProperty::nTruthJetCandidateProperties));
    }
  }

  if (options.saveMCObjects) {
    setSelectedPhotonClosestJet(selectedMediumPhotonProperties, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setSelectedPhotonClosestJet(selectedMediumPhotonProperties_closeToTruePhoton, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setSelectedPhotonClosestJet(selectedMediumPhotonProperties_awayFromTruePhoton, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setSelectedPhotonClosestJet(selectedVetoedPhotonProperties, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setSelectedPhotonClosestJet(selectedVetoedPhotonProperties_closeToTruePhoton, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setSelectedPhotonClosestJet(selectedVetoedPhotonProperties_awayFromTruePhoton, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setSelectedPhotonClosestJet(selectedFakePhotonProperties, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setSelectedPhotonClosestJet(selectedFakePhotonProperties_closeToTruePhoton, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setSelectedPhotonClosestJet(selectedFakePhotonProperties_awayFromTruePhoton, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);

    setUnselectedMediumPhotonClosestJet(unselected_medium_pho_properties, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setUnselectedMediumPhotonClosestJet(unselected_medium_pho_properties_closeToTruePhoton, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setUnselectedMediumPhotonClosestJet(unselected_medium_pho_properties_awayFromTruePhoton, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setUnselectedVetoedPhotonClosestJet(unselected_vetoed_pho_properties, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setUnselectedVetoedPhotonClosestJet(unselected_vetoed_pho_properties_closeToTruePhoton, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setUnselectedVetoedPhotonClosestJet(unselected_vetoed_pho_properties_awayFromTruePhoton, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setUnselectedFakePhotonClosestJet(unselected_fake_pho_properties, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setUnselectedFakePhotonClosestJet(unselected_fake_pho_properties_closeToTruePhoton, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
    setUnselectedFakePhotonClosestJet(unselected_fake_pho_properties_awayFromTruePhoton, selectedJetAngles, genJetAngles, eventProgenitorMomGenJetAngles, singletMomGenJetAngles);
  }

  selectionBits[eventSelectionCriterion::NJets] = (max_nJets >= 2);
  if (options.disableJetSelection) selectionBits[eventSelectionCriterion::NJets] = true;
  // Add MET to ST
  event_ST_MET += eventDetails.PFMET;
  event_ST += eventDetails.PFMET;
  event_properties[eventProperty::ST_electromagnetic] = event_ST_electromagnetic;
  event_properties[eventProperty::ST_hadronic] = event_ST_hadronic;
  event_properties[eventProperty::ST_MET] = event_ST_MET;
  event_properties[eventProperty::ST] = event_ST;
  event_properties[eventProperty::selectionRegionIndex] = 1.0*static_cast<int>(region);

  if (options.calculateShiftedDistributions) {
    // Add shifted energies
    float HEMAdjustedMETX = (eventDetails.PFMET)*std::cos(eventDetails.PFMET_phi) + MET_HEMAdjustmentX;
    float HEMAdjustedMETY = (eventDetails.PFMET)*std::sin(eventDetails.PFMET_phi) + MET_HEMAdjustmentY;
    float HEMAdjustedMETMagnitude = std::sqrt(std::pow(HEMAdjustedMETX, 2) + std::pow(HEMAdjustedMETY, 2));

    for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
      shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);

      // bool passes_selection = jetExaminationResults.passesSelectionDRJECNominal;
      float shifted_METcontribution = eventDetails.PFMET;
      if (typeIndex == shiftType::UnclusteredMETDown) {
	shifted_METcontribution = eventDetails.PFMET_UnclusteredDown;
      }
      else if (typeIndex == shiftType::UnclusteredMETUp) {
	shifted_METcontribution = eventDetails.PFMET_UnclusteredUp;
      }
      else if (typeIndex == shiftType::JERMETDown) {
	shifted_METcontribution = eventDetails.PFMET_JERDown;
      }
      else if (typeIndex == shiftType::JERMETUp) {
	shifted_METcontribution = eventDetails.PFMET_JERUp;
      }
      else if (typeIndex == shiftType::missingHEMDown) {
	shifted_METcontribution = HEMAdjustedMETMagnitude;
      }
      else if (typeIndex == shiftType::missingHEMUp) {
	shifted_METcontribution = eventDetails.PFMET;
      }
      addShiftedEToSTMap(shifted_METcontribution, shifted_ST, typeIndex);
    }
    // addShiftedEToSTMap(eventDetails.PFMET, shifted_ST, shiftType::JECDown);
    // addShiftedEToSTMap(eventDetails.PFMET, shifted_ST, shiftType::JECUp);
    // addShiftedEToSTMap(eventDetails.PFMET_UnclusteredDown, shifted_ST, shiftType::UnclusteredMETDown);
    // addShiftedEToSTMap(eventDetails.PFMET_UnclusteredUp, shifted_ST, shiftType::UnclusteredMETUp);
    // addShiftedEToSTMap(eventDetails.PFMET_JERDown, shifted_ST, shiftType::JERMETDown);
    // addShiftedEToSTMap(eventDetails.PFMET_JERUp, shifted_ST, shiftType::JERMETUp);
    // addShiftedEToSTMap(HEMAdjustedMETMagnitude, shifted_ST, shiftType::missingHEMDown);
    // addShiftedEToSTMap(eventDetails.PFMET, shifted_ST, shiftType::missingHEMUp);
  }

  assert(static_cast<int>(selectionBits.size()) == static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria));

  int nEventFalseBits = getNFalseBits(selectionBits);
  statistics.fillIDEfficiencyTimesAcceptanceStatisticsHistograms(event_ST, n_jetsDR, (nEventFalseBits == 0), region, MCRegionIndex);
  statistics.fillIDEfficiencyStatisticsHistograms(event_ST, n_jetsDR, (nEventFalseBits == 0), passesExtendedMCSelection, region, MCRegionIndex);

  eventProperties temp1 = initialize_eventProperties_with_defaults(); // temp1 and temp2 are dummies -- they won't contribute to the histograms
  if (nEventFalseBits == 0) {
    unselectedEventProperties temp2 = std::make_pair(eventSelectionCriterion::nEventSelectionCriteria, temp1);
    statistics.fill1DStatisticsHistograms(event_properties, false, temp2,
                                          selectedTruePhotonProperties,
                                          selectedTrueJetCandidateProperties_all,
                                          selectedTrueJetCandidateProperties_fromEventProgenitor,
                                          selectedTrueJetCandidateProperties_fromSinglet,
                                          selectedMediumPhotonProperties,
                                          selectedMediumPhotonProperties_closeToTruePhoton,
                                          selectedMediumPhotonProperties_awayFromTruePhoton,
                                          unselected_medium_pho_properties,
                                          unselected_medium_pho_properties_closeToTruePhoton,
                                          unselected_medium_pho_properties_awayFromTruePhoton,
					  selectedVetoedPhotonProperties,
                                          selectedVetoedPhotonProperties_closeToTruePhoton,
                                          selectedVetoedPhotonProperties_awayFromTruePhoton,
                                          unselected_vetoed_pho_properties,
                                          unselected_vetoed_pho_properties_closeToTruePhoton,
                                          unselected_vetoed_pho_properties_awayFromTruePhoton,
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
                                          gen_jet_properties_collection,
                                          eventProgenitor_mom_gen_jet_properties_collection,
                                          singlet_mom_gen_jet_properties_collection,
                                          region, MCRegionIndex);
  }
  else if (nEventFalseBits == 1) {
    eventSelectionCriterion marginallyUnselectedEventCriterion = getFirstFalseCriterion(selectionBits);
    unselectedEventProperties unselected_event_properties = std::make_pair(marginallyUnselectedEventCriterion, event_properties);
    statistics.fill1DStatisticsHistograms(temp1, true, unselected_event_properties,
                                          selectedTruePhotonProperties,
                                          selectedTrueJetCandidateProperties_all,
                                          selectedTrueJetCandidateProperties_fromEventProgenitor,
                                          selectedTrueJetCandidateProperties_fromSinglet,
                                          selectedMediumPhotonProperties,
                                          selectedMediumPhotonProperties_closeToTruePhoton,
                                          selectedMediumPhotonProperties_awayFromTruePhoton,
                                          unselected_medium_pho_properties,
                                          unselected_medium_pho_properties_closeToTruePhoton,
                                          unselected_medium_pho_properties_awayFromTruePhoton,
                                          selectedVetoedPhotonProperties,
                                          selectedVetoedPhotonProperties_closeToTruePhoton,
                                          selectedVetoedPhotonProperties_awayFromTruePhoton,
                                          unselected_vetoed_pho_properties,
                                          unselected_vetoed_pho_properties_closeToTruePhoton,
                                          unselected_vetoed_pho_properties_awayFromTruePhoton,
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
                                          gen_jet_properties_collection,
                                          eventProgenitor_mom_gen_jet_properties_collection,
                                          singlet_mom_gen_jet_properties_collection,
                                          region, MCRegionIndex);
  }

  if ((((selectionBits.at(eventSelectionCriterion::doublePhoton)) && (selectionBits.at(eventSelectionCriterion::invariantMass))) &&
       (selectionBits.at(eventSelectionCriterion::MCGenInformation))) && (!(options.disablePhotonSelection))) {
    statistics.fillHLTEfficiencyStatisticsHistograms(eta_leadingPhoton, pT_leadingPhoton,
                                                     eta_subLeadingPhoton, pT_subLeadingPhoton,
                                                     selectionBits.at(eventSelectionCriterion::HLTSelection),
                                                     region);
  }

  eventResult.isInterestingEvent = (((nEventFalseBits == 0) && (event_ST >= (STRegions.STNormRangeMin - parameters.preNormalizationBuffer))) && eventPUIsSensible);
  if (options.disablePhotonSelection) { // option "disablePhotonSelection" is meant to be used to produce pure MC samples.
    eventResult.isInterestingEvent = ((nEventFalseBits == 0) && eventPUIsSensible); // will only select events with at least 2 jets, with no criterion on ST
  }
  if (options.disableJetSelection) {
    eventResult.isInterestingEvent = ((nEventFalseBits == 0) && eventPUIsSensible); // if jet selection is disabled, we want to save events regardless of ST
  }

  eventResult.evt_photonIndex_leading = index_leadingPhoton;
  eventResult.evt_photonIndex_subLeading = index_subLeadingPhoton;
  eventResult.evt_photonPT_leading = pT_leadingPhoton;
  eventResult.evt_photonPT_subLeading = pT_subLeadingPhoton;
  eventResult.evt_photonEta_leading = eta_leadingPhoton;
  eventResult.evt_photonEta_subLeading = eta_subLeadingPhoton;
  eventResult.evt_photonMVA_leading = mva_leadingPhoton;
  eventResult.evt_photonMVA_subLeading = mva_subLeadingPhoton;
  eventResult.evt_jetPT_leading = pT_leadingJet;

  if (nEventFalseBits <= 1) assert(static_cast<int>(event_properties.size()) == static_cast<int>(eventProperty::nEventProperties));
  return eventResult;
}

void loopOverEvents(optionsStruct &options, parametersStruct &parameters, // const int& year,
                    std::vector<eventExaminationResultsStruct>& selectedEventsInfo, statisticsHistograms& statistics, STRegionsStruct& STRegions) {
  TChain inputChain("ggNtuplizer/EventTree");
  std::cout << "Starting to add files to chain..." << std::endl;
  for (const std::string& inputPath: (options.inputPaths)) {
    std::cout << "Adding... " << inputPath << std::endl;
    int read_status = inputChain.Add(inputPath.c_str());
    assert(read_status == 1);
  }
  std::cout << "Finished adding files to chain!" << std::endl;

  inputChain.SetBranchStatus("*", 0); // so that only the needed branches, explicitly activated below, are read in per event

  if (options.doOverlapRemoval) {
    assert(options.overlapRemoval_maxNPromptPhotons >= 0);
    assert(options.saveMCGenLevelInfo);
  }
  if (options.enableMCEventFilter) assert(options.saveMCMomInfo);
  eventDetailsStruct eventDetails = eventDetailsStruct(inputChain, ((options.saveMCGenLevelInfo) || (options.enableMCEventFilter) || (options.doOverlapRemoval)), options.calculateShiftedDistributions, options.savePUWeights);
  photonsCollectionStruct photonsCollection = photonsCollectionStruct(inputChain);
  jetsCollectionStruct jetsCollection = jetsCollectionStruct(inputChain, options.saveMCObjects, options.calculateShiftedDistributions);
  MCCollectionStruct MCCollection = MCCollectionStruct(inputChain, ((options.saveMCGenLevelInfo) || (options.enableMCEventFilter) || (options.doOverlapRemoval)));

  Long64_t nEvts = inputChain.GetEntries();
  std::cout << "Number of events to process: " << nEvts << std::endl;

  tmProgressBar progressBar = tmProgressBar(static_cast<int>(nEvts));
  int progressBarUpdatePeriod = ((nEvts < 50) ? 1 : static_cast<int>(0.5 + 1.0*(nEvts/50)));
  progressBar.initialize();
  int nProblematicEntries = 0;
  int n_events_without_pu_info = 0;
  for (Long64_t entryIndex = 0; entryIndex < nEvts; ++entryIndex) {
    Long64_t loadStatus = inputChain.LoadTree(entryIndex);
    if (loadStatus < 0) {
      std::cout << "Warning: loadStatus < 0 for entry index: " << entryIndex << "; load status = " << loadStatus << "; fileName: " << (inputChain.GetFile())->GetName() << std::endl;
      ++nProblematicEntries;
      assert(nProblematicEntries <= N_PROBLEMATIC_ENTRIES_THRESHOLD);
      continue;
    }
    int nBytesRead = inputChain.GetEntry(entryIndex, 0); // Get only the required branches
    if (nBytesRead <= 0) {
      std::cout << "Warning: failed to read SOME information from entry at index: " << entryIndex << "; nBytesRead = " << nBytesRead << "; fileName: " << (inputChain.GetFile())->GetName() << std::endl;
      ++nProblematicEntries;
      assert(nProblematicEntries <= N_PROBLEMATIC_ENTRIES_THRESHOLD);
      continue;
    }

    int entryProcessing = static_cast<int>(entryIndex);
    if (entryProcessing > 0 && ((static_cast<int>(entryProcessing) % progressBarUpdatePeriod == 0) || entryProcessing == static_cast<int>(nEvts-1))) progressBar.updateBar(static_cast<double>(1.0*entryProcessing/nEvts), entryProcessing);

    eventExaminationResultsStruct eventExaminationResults = examineEvent(options, parameters, // counters, 
                                                                         entryIndex, n_events_without_pu_info, // year,
                                                                         eventDetails, MCCollection,
									 photonsCollection, jetsCollection,
									 statistics, STRegions);
    bool eventIsToBeStored = eventExaminationResults.isInterestingEvent;
    // incrementCounters(miscCounter::totalEvents, counters, false, 0., 0.);
    if (!(eventIsToBeStored)) {
      // incrementCounters(miscCounter::failingEvents, counters, false, 0., 0.);
      continue;
    }
    // incrementCounters(miscCounter::acceptedEvents, counters, false, 0., 0.);
    selectedEventsInfo.push_back(eventExaminationResults);
  }
  progressBar.terminate();
  if (n_events_without_pu_info > 0) {
    std::cout << "WARNING: " << n_events_without_pu_info << "/" << nEvts << " events (" << ((100.0*n_events_without_pu_info)/(1.0*nEvts)) << "%) without PU info." << std::endl;
  }
}

void writeSelectionToFile(optionsStruct &options, TFile *outputFile, const std::vector<eventExaminationResultsStruct>& selectedEventsInfo, const selectionRegion& region, const bool& restrictToRegion) {
  std::string regionName = "";
  if (restrictToRegion) {
    regionName = selectionRegionNames.at(region);
    std::cout << "Beginning to write selected events to file for selection type: " <<  regionName << std::endl;
  }
  else {
    std::cout << "Beginning to write unified selected events to file..." <<  regionName << std::endl;
  }

  TChain inputChain("ggNtuplizer/EventTree");
  std::cout << "Starting to add files to chain..." << std::endl;
  for (const std::string& inputPath: (options.inputPaths)) {
    // std::cout << "Adding... " << inputPath << std::endl;
    int read_status = inputChain.Add(inputPath.c_str());
    assert(read_status == 1);
  }
  std::cout << "Finished adding files to chain!" << std::endl;
  inputChain.SetBranchStatus("*", 1); // so that all branches are read in

  TDirectory *outDir = outputFile->mkdir("ggNtuplizer");
  outDir->cd();
  TTree *outputTree = inputChain.CloneTree(0);

  double PUWeight; // stores PU weight, might be used for PU reweighting
  if (options.savePUWeights) {
    outputTree->Branch("b_PUWeightNoSelection", &PUWeight, "b_PUWeightNoSelection/D");
  }
  int nKinematicMCPhotons; // stores number of MC photons in the barrel passing subleading photon ET cut
  outputTree->Branch("b_nKinematicMCPhotons", &nKinematicMCPhotons, "b_nKinematicMCPhotons/I");
  int nRecoPhotonsMatchedToMCGenPhotons; // stores number of MC photons in the barrel passing subleading photon ET cut
  outputTree->Branch("b_nRecoPhotonsMatchedToMCGenPhotons", &nRecoPhotonsMatchedToMCGenPhotons, "b_nRecoPhotonsMatchedToMCGenPhotons/I");
  float deltaR_genPhoton_MCMom_matchingLeadingPhoton;
  outputTree->Branch("b_deltaR_genPhoton_MCMom_matchingLeadingPhoton", &deltaR_genPhoton_MCMom_matchingLeadingPhoton, "b_deltaR_genPhoton_MCMom_matchingLeadingPhoton/F");
  float deltaR_genPhoton_MCMom_matchingSubLeadingPhoton;
  outputTree->Branch("b_deltaR_genPhoton_MCMom_matchingSubLeadingPhoton", &deltaR_genPhoton_MCMom_matchingSubLeadingPhoton, "b_deltaR_genPhoton_MCMom_matchingSubLeadingPhoton/F");
  float event_deltaR_photons; // stores deltaR between two reco-level photons
  outputTree->Branch("b_deltaR_photons", &event_deltaR_photons, "b_deltaR_photons/F");
  int nJetsDR; // stores number of jets in event passing deltaR cut
  outputTree->Branch("b_nJetsDR", &nJetsDR, "b_nJetsDR/I");
  int nJetsAll; // stores total number of jets in event whether or not they pass deltaR cut
  outputTree->Branch("b_nJetsAll", &nJetsAll, "b_nJetsAll/I");
  float ST; // stores event sT
  outputTree->Branch("b_evtST", &ST, "b_evtST/F");
  float invMass; // stores invariant mass of two selected photons in the event
  outputTree->Branch("b_invMass", &invMass, "b_invMass/F");
  float ST_electromagnetic; // stores event sT
  outputTree->Branch("b_evtST_electromagnetic", &ST_electromagnetic, "b_evtST_electromagnetic/F");
  float ST_hadronic; // stores event sT
  outputTree->Branch("b_evtST_hadronic", &ST_hadronic, "b_evtST_hadronic/F");
  float ST_MET; // stores event sT
  outputTree->Branch("b_evtST_MET", &ST_MET, "b_evtST_MET/F");
  int photonIndex_leading;
  outputTree->Branch("b_photonIndex_leading", &photonIndex_leading, "b_photonIndex_leading/I");
  int photonIndex_subLeading;
  outputTree->Branch("b_photonIndex_subLeading", &photonIndex_subLeading, "b_photonIndex_subLeading/I");
  float photonPT_leading; // stores PT of leading photon, useful for HLT efficiency
  outputTree->Branch("b_photonPT_leading", &photonPT_leading);
  float photonPT_subLeading; // stores PT of subleading photon, useful for HLT efficiency
  outputTree->Branch("b_photonPT_subLeading", &photonPT_subLeading);
  float photonEta_leading; // stores eta of leading photon, useful for HLT efficiency
  outputTree->Branch("b_photonEta_leading", &photonEta_leading);
  float photonEta_subLeading; // stores eta of subleading photon, useful for HLT efficiency
  outputTree->Branch("b_photonEta_subLeading", &photonEta_subLeading);
  float photonMVA_leading; // stores MVA ID of leading photon
  outputTree->Branch("b_photonMVA_leading", &photonMVA_leading);
  float photonMVA_subLeading; // stores MVA ID of subleading photon
  outputTree->Branch("b_photonMVA_subLeading", &photonMVA_subLeading);
  int nPhotonsMatchedToGenPromptFS; // stores the number of photons matched to gen-level photons with the "prompt final state" flag set
  outputTree->Branch("b_nPhotonsMatchedToGenPromptFS", &nPhotonsMatchedToGenPromptFS, "b_nPhotonsMatchedToGenPromptFS/I");
  float jetPT_leading; // stores pT of leading jet
  outputTree->Branch("b_jetPT_leading", &jetPT_leading);
  eventWeightsStruct prefireWeights = eventWeightsStruct(1.0f, 1.0f, 1.0f); // stores prefiring weights and errors
  outputTree->Branch("b_evtPrefiringWeight", &(prefireWeights.nominal), "b_evtPrefiringWeight/F");
  outputTree->Branch("b_evtPrefiringWeightDown", &(prefireWeights.down), "b_evtPrefiringWeightDown/F");
  outputTree->Branch("b_evtPrefiringWeightUp", &(prefireWeights.up), "b_evtPrefiringWeightUp/F");
  eventWeightsStruct photonMCScaleFactors = eventWeightsStruct(1.0f, 1.0f, 1.0f); // stores scale factors for MC samples

  std::map<shiftType, float> shifted_ST;
  std::map<shiftType, int> shifted_nJetsDR;
  std::map<shiftType, int> shifted_nJetsAll;
  if (options.calculateShiftedDistributions) {
    for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
      shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
      shifted_ST[typeIndex] = 0.;
      outputTree->Branch((getShiftedVariableBranchName(typeIndex, "evtST")).c_str(), &(shifted_ST[typeIndex]), (getShiftedVariableBranchName(typeIndex, "evtST") + "/F").c_str());
      shifted_nJetsDR[typeIndex] = 0;
      outputTree->Branch((getShiftedVariableBranchName(typeIndex, "nJetsDR")).c_str(), &(shifted_nJetsDR[typeIndex]), (getShiftedVariableBranchName(typeIndex, "nJetsDR") + "/I").c_str());
      shifted_nJetsAll[typeIndex] = 0;
      outputTree->Branch((getShiftedVariableBranchName(typeIndex, "nJetsAll")).c_str(), &(shifted_nJetsAll[typeIndex]), (getShiftedVariableBranchName(typeIndex, "nJetsAll") + "/I").c_str());
    }
  }
  if (options.calculateMCScaleFactorWeights) {
    outputTree->Branch("b_evtphotonMCScaleFactor", &(photonMCScaleFactors.nominal), "b_evtphotonMCScaleFactor/F");
    outputTree->Branch("b_evtphotonMCScaleFactorDown", &(photonMCScaleFactors.down), "b_evtphotonMCScaleFactorDown/F");
    outputTree->Branch("b_evtphotonMCScaleFactorUp", &(photonMCScaleFactors.up), "b_evtphotonMCScaleFactorUp/F");
  }

  double MCXSecWeight = options.MCBackgroundWeight;
  if (options.saveMCBackgroundWeight) {
    outputTree->Branch("b_MCXSecWeight", &MCXSecWeight, "b_MCXSecWeight/D");
  }

  int nSelectedEvents = static_cast<int>(0.5 + selectedEventsInfo.size());
  tmProgressBar progressBar = tmProgressBar(nSelectedEvents);
  int progressBarUpdatePeriod = ((nSelectedEvents < 50) ? 1 : static_cast<int>(0.5 + 1.0*(nSelectedEvents/50)));
  int processingIndex = -1;
  progressBar.initialize();

  for (auto& selectedEventInfo : selectedEventsInfo) {
    ++processingIndex;

    Long64_t index = selectedEventInfo.eventIndex;
    PUWeight = selectedEventInfo.evt_PUWeight;
    nKinematicMCPhotons = (selectedEventInfo.evt_gen_level_info).nKinematicPhotons;
    nRecoPhotonsMatchedToMCGenPhotons = (selectedEventInfo.evt_gen_level_info).nRecoPhotonsMatchedToGenPhotons;
    deltaR_genPhoton_MCMom_matchingLeadingPhoton = (selectedEventInfo.evt_gen_level_info).deltaR_genPhoton_mom_matchingLeadingPhoton;
    deltaR_genPhoton_MCMom_matchingSubLeadingPhoton = (selectedEventInfo.evt_gen_level_info).deltaR_genPhoton_mom_matchingSubLeadingPhoton;
    event_deltaR_photons = selectedEventInfo.evt_deltaR_photons;
    nJetsDR = selectedEventInfo.evt_nJetsDR;
    nJetsAll = selectedEventInfo.evt_nJetsAll;
    invMass = selectedEventInfo.evt_invariantMass;
    ST_electromagnetic = selectedEventInfo.evt_ST_electromagnetic;
    ST_hadronic = selectedEventInfo.evt_ST_hadronic;
    ST_MET = selectedEventInfo.evt_ST_MET;
    ST = selectedEventInfo.evt_ST;
    photonIndex_leading = selectedEventInfo.evt_photonIndex_leading;
    photonIndex_subLeading = selectedEventInfo.evt_photonIndex_subLeading;
    photonPT_leading = selectedEventInfo.evt_photonPT_leading;
    photonPT_subLeading = selectedEventInfo.evt_photonPT_subLeading;
    photonEta_leading = selectedEventInfo.evt_photonEta_leading;
    photonEta_subLeading = selectedEventInfo.evt_photonEta_subLeading;
    photonMVA_leading = selectedEventInfo.evt_photonMVA_leading;
    photonMVA_subLeading = selectedEventInfo.evt_photonMVA_subLeading;
    nPhotonsMatchedToGenPromptFS = selectedEventInfo.evt_nPhotonsMatchedToGenPromptFS;
    jetPT_leading = selectedEventInfo.evt_jetPT_leading;
    prefireWeights = eventWeightsStruct((selectedEventInfo.evt_prefireWeights).nominal, (selectedEventInfo.evt_prefireWeights).down, (selectedEventInfo.evt_prefireWeights).up);
    if (options.calculateShiftedDistributions) {
      for (int shiftTypeIndex = shiftTypeFirst; shiftTypeIndex != static_cast<int>(shiftType::nShiftTypes); ++shiftTypeIndex) {
        shiftType typeIndex = static_cast<shiftType>(shiftTypeIndex);
        shifted_ST[typeIndex] = (selectedEventInfo.evt_shifted_ST).at(typeIndex);
        shifted_nJetsDR[typeIndex] = (selectedEventInfo.evt_shifted_nJetsDR).at(typeIndex);
        shifted_nJetsAll[typeIndex] = (selectedEventInfo.evt_shifted_nJetsAll).at(typeIndex);
      }
    }
    if (options.calculateMCScaleFactorWeights) {
      photonMCScaleFactors = eventWeightsStruct((selectedEventInfo.evt_photonMCScaleFactors).nominal, (selectedEventInfo.evt_photonMCScaleFactors).down, (selectedEventInfo.evt_photonMCScaleFactors).up);
    }
    // std::cout << "Nominal ST: " << ST << std::endl;
    // std::cout << "Nominal nJets: " << nJetsDR << std::endl;
    // std::cout << "shifted_nJets: " << std::endl;
    // printShiftedVariablesMap(shifted_nJetsDR, "shifted_nJetsDR");
    // std::cout << "shifted_ST: " << std::endl;
    // printShiftedVariablesMap(shifted_ST, "shifted_ST");

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
    if (restrictToRegion) {
      if (selectedEventInfo.evt_region != region) continue;
    }

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
  argumentParser.addArgument("inputPathsFile", "", true, "Path to file containing list of input files.");
  argumentParser.addArgument("selectionType", "default", true, "Selection type. Currently only allowed to be \"data\", \"data_singlephoton\", \"data_jetHT\", \"MC_stealth_t5\", \"MC_stealth_t6\", \"MC_DiPhotonJets(|Box)\", \"MC_EMEnrichedGJetPt16_[N]\", \"MC_EMEnrichedGJetPt17_[N]\", \"MC_EMEnrichedGJetPt18_[N]\", \"MC_HighHTQCD16_[N]\", \"MC_HighHTQCD17_[N]\", \"MC_HighHTQCD18_[N]\", \"MC_GJetHT16_[N]\", \"MC_GJetHT17_[N]\", \"MC_GJetHT18_[N]\", or \"MC_hgg\".");
  argumentParser.addArgument("disablePhotonSelection", "default", true, "Do not filter on photon criteria. Used to build \"pure MC\" samples.");
  argumentParser.addArgument("disableJetSelection", "default", true, "Do not filter on nJets.");
  argumentParser.addArgument("lineNumberStartInclusive", "", true, "Line number from input file from which to start. The file with this index is included in the processing.");
  argumentParser.addArgument("lineNumberEndInclusive", "", true, "Line number from input file at which to end. The file with this index is included in the processing.");
  argumentParser.addArgument("year", "", true, "Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.");
  argumentParser.addArgument("invertElectronVeto", "default", true, "Invert the electron veto condition on selected photons; meant to be used to estimate trigger efficiency.");
  argumentParser.addArgument("MCBackgroundWeight", "-1.0", true, "(meant for background MC samples) create a branch to save MC event weights given the cross section of the process.");
  argumentParser.addArgument("PUWeightsPathWithXRDPrefix", "/dev/null", true, "If event-dependent PU weights are to be saved, then this option should be set to the path to the source file from which to obtain the histogram with PU-dependent weights.");
  argumentParser.setPassedStringValues(argc, argv);

  optionsStruct options = getOptionsFromParser(argumentParser);

  parametersStruct parameters = parametersStruct();
  parameters.tuneParameters(options.year, options.calculateMCScaleFactorWeights, options.savePUWeights, options.PUWeightsPathWithXRDPrefix, options.selectionType);

  std::stringstream optionsStringstream;
  optionsStringstream << options;
  std::stringstream parametersStringstream;
  parametersStringstream << parameters;

  std::vector<eventExaminationResultsStruct> selectedEventsInfo;

  STRegionsStruct STRegions = STRegionsStruct("STRegionBoundaries.dat", 20000.0);

  statisticsHistograms statistics = statisticsHistograms(options.saveMCObjects, false, HLTEmulation::etaBinEdges, HLTEmulation::pTBinEdges, STRegions.STBoundaries);

  loopOverEvents(options, parameters, // options.year,
                 selectedEventsInfo, statistics, STRegions);

  statistics.writeToFile("statisticsHistograms.root");

  for (const selectionRegion& region: options.selectionsToWrite) {
    // for (int selectionRegionIndex = selectionRegionFirst; selectionRegionIndex != static_cast<int>(selectionRegion::nSelectionRegions); ++selectionRegionIndex) {
    // selectionRegion region = static_cast<selectionRegion>(selectionRegionIndex);
    // bool write_selection = true;
    // if ((region == selectionRegion::control_singlemedium) || (region == selectionRegion::control_singleloose) || (region == selectionRegion::control_singlefake)) {
    //   write_selection = false;
    //   if (((options.selectionType == "data_singlephoton") ||
    //        (std::regex_match(options.selectionType, std::regex("^MC_GJet16_singlephoton[0-9]*$"))) ||
    //        (std::regex_match(options.selectionType, std::regex("^MC_GJet17_singlephoton[0-9]*$"))) ||
    //        (std::regex_match(options.selectionType, std::regex("^MC_GJet18_singlephoton[0-9]*$"))) ||
    //        (std::regex_match(options.selectionType, std::regex("^MC_QCD16_singlephoton[0-9]*$"))) ||
    //        (std::regex_match(options.selectionType, std::regex("^MC_QCD17_singlephoton[0-9]*$"))) ||
    //        (std::regex_match(options.selectionType, std::regex("^MC_QCD18_singlephoton[0-9]*$")))) &&
    //       (!(options.isMC))) write_selection = true;
    // }
    // else {
    //   if ((options.selectionType == "data_singlephoton") ||
    //       (std::regex_match(options.selectionType, std::regex("^MC_GJet16_singlephoton[0-9]*$"))) ||
    //       (std::regex_match(options.selectionType, std::regex("^MC_GJet17_singlephoton[0-9]*$"))) ||
    //       (std::regex_match(options.selectionType, std::regex("^MC_GJet18_singlephoton[0-9]*$"))) ||
    //       (std::regex_match(options.selectionType, std::regex("^MC_QCD16_singlephoton[0-9]*$"))) ||
    //       (std::regex_match(options.selectionType, std::regex("^MC_QCD17_singlephoton[0-9]*$"))) ||
    //       (std::regex_match(options.selectionType, std::regex("^MC_QCD18_singlephoton[0-9]*$")))) write_selection = false;
    // }
    // if (options.selectionType == "data_jetHT") write_selection = false; // we are only interested in the statistics histograms for JetHT data
    // if (options.disablePhotonSelection) write_selection = false; // special case covered below
    // if (!(write_selection)) continue;
    if (options.disablePhotonSelection) continue; // special case covered below
    std::string outputFilePath = std::string("selection_") + selectionRegionNames.at(region) + std::string(".root");
    TFile *outputFile = TFile::Open(outputFilePath.c_str(), "RECREATE");
    if (!(outputFile->IsOpen()) || outputFile->IsZombie()) {
      std::cout << "ERROR: Unable to open output file to write. Attempted to create file with path: " << outputFilePath << std::endl;
    }
    writeSelectionToFile(options, outputFile, selectedEventsInfo, region, true);
    outputFile->Close();
  }

  if (options.disablePhotonSelection) {
    std::string outputFilePath = std::string("selection_unified.root");
    TFile *outputFile = TFile::Open(outputFilePath.c_str(), "RECREATE");
    if (!(outputFile->IsOpen()) || outputFile->IsZombie()) {
      std::cout << "ERROR: Unable to open output file to write. Attempted to create file with path: " << outputFilePath << std::endl;
    }
    selectionRegion region = selectionRegion::nSelectionRegions;
    writeSelectionToFile(options, outputFile, selectedEventsInfo, region, false);
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
