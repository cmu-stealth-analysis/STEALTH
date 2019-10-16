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

bool passesHLTEmulation(const int& year, photonProperties& properties_leadingPhoton, photonProperties& properties_subLeadingPhoton, const int& triggerBit) {
  bool passesEmulation = false;
  switch(year) {
  case 2016:
    {
      if (triggerBit == 16) {
	bool leadingPassesEmulation = false;
	float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
	float& leading_photon_energy = properties_leadingPhoton[photonProperty::energy];
	float& leading_photon_R9 = properties_leadingPhoton[photonProperty::R9];
	float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];
	float& leading_photon_sigmaIEtaIEta = properties_leadingPhoton[photonProperty::sigmaIEtaIEta];
	float& leading_photon_ecalClusIso = properties_leadingPhoton[photonProperty::ecalClusIso];

	if (leading_photon_R9 >= 0.85) {
	  leadingPassesEmulation = ((leading_photon_pT > 30.) &&
				    (leading_photon_HOverE <= 0.1*leading_photon_energy));
	}
	else if ((leading_photon_R9 >= 0.5) && (leading_photon_R9 < 0.85)) {
	  leadingPassesEmulation = ((leading_photon_pT > 30.) &&
				    (leading_photon_HOverE <= 0.1*leading_photon_energy) &&
				    (leading_photon_sigmaIEtaIEta <= 0.015) &&
				    (leading_photon_ecalClusIso <= (6.0 + 0.012*leading_photon_pT)));
	}

	bool subLeadingPassesEmulation = false;
	float& subLeading_photon_pT = properties_subLeadingPhoton[photonProperty::pT];
	float& subLeading_photon_energy = properties_subLeadingPhoton[photonProperty::energy];
	float& subLeading_photon_R9 = properties_subLeadingPhoton[photonProperty::R9];
	float& subLeading_photon_HOverE = properties_subLeadingPhoton[photonProperty::hOverE];
	float& subLeading_photon_sigmaIEtaIEta = properties_subLeadingPhoton[photonProperty::sigmaIEtaIEta];
	float& subLeading_photon_ecalClusIso = properties_subLeadingPhoton[photonProperty::ecalClusIso];
	float& subLeading_photon_trkIso = properties_subLeadingPhoton[photonProperty::trkIso];

	if (subLeading_photon_R9 >= 0.85) {
	  subLeadingPassesEmulation = ((subLeading_photon_pT > 18.) &&
				       (subLeading_photon_HOverE <= 0.1*subLeading_photon_energy));
	}
	else if ((subLeading_photon_R9 >= 0.5) && (subLeading_photon_R9 < 0.85)) {
	  subLeadingPassesEmulation = ((subLeading_photon_pT > 18.) &&
				       (subLeading_photon_HOverE <= 0.1*subLeading_photon_energy) &&
				       (subLeading_photon_sigmaIEtaIEta <= 0.015) &&
				       (subLeading_photon_ecalClusIso <= (6.0 + 0.012*subLeading_photon_pT)) &&
				       (subLeading_photon_trkIso <= (6.0 + 0.002*subLeading_photon_pT)));
	}

	passesEmulation = leadingPassesEmulation && subLeadingPassesEmulation;
      }
      else if (triggerBit == 22) {
	bool leadingPassesEmulation = false;
	float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
	float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];

	leadingPassesEmulation = ((leading_photon_pT > 60.) &&
				  (leading_photon_HOverE <= 0.15));

	bool subLeadingPassesEmulation = false;
	float& subLeading_photon_pT = properties_subLeadingPhoton[photonProperty::pT];
	float& subLeading_photon_HOverE = properties_subLeadingPhoton[photonProperty::hOverE];

	subLeadingPassesEmulation = ((subLeading_photon_pT > 60.) &&
				     (subLeading_photon_HOverE) <= 0.15);

	passesEmulation = leadingPassesEmulation && subLeadingPassesEmulation;
      }
      else {
	std::cout << "ERROR: unsupported trigger bit: " << triggerBit << std::endl;
	std::exit(EXIT_FAILURE);
      }
      break;
    }
  case 2017:
  case 2018:
    {
      if (triggerBit == 37) {
	bool leadingPassesEmulation = false;
	float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
	float& leading_photon_energy = properties_leadingPhoton[photonProperty::energy];
	float& leading_photon_R9 = properties_leadingPhoton[photonProperty::R9];
	float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];
	float& leading_photon_sigmaIEtaIEta = properties_leadingPhoton[photonProperty::sigmaIEtaIEta];
	float& leading_photon_ecalClusIso = properties_leadingPhoton[photonProperty::ecalClusIso];

	if (leading_photon_R9 >= 0.85) {
	  leadingPassesEmulation = ((leading_photon_pT > 30.) &&
				    (leading_photon_HOverE <= 0.1*leading_photon_energy));
	}
	else if ((leading_photon_R9 >= 0.5) && (leading_photon_R9 < 0.85)) {
	  leadingPassesEmulation = ((leading_photon_pT > 30.) &&
				    (leading_photon_HOverE <= 0.1*leading_photon_energy) &&
				    (leading_photon_sigmaIEtaIEta <= 0.015) &&
				    (leading_photon_ecalClusIso <= (6.0 + 0.012*leading_photon_pT)));
	}

	bool subLeadingPassesEmulation = false;
	float& subLeading_photon_pT = properties_subLeadingPhoton[photonProperty::pT];
	float& subLeading_photon_energy = properties_subLeadingPhoton[photonProperty::energy];
	float& subLeading_photon_R9 = properties_subLeadingPhoton[photonProperty::R9];
	float& subLeading_photon_HOverE = properties_subLeadingPhoton[photonProperty::hOverE];
	float& subLeading_photon_sigmaIEtaIEta = properties_subLeadingPhoton[photonProperty::sigmaIEtaIEta];
	float& subLeading_photon_ecalClusIso = properties_subLeadingPhoton[photonProperty::ecalClusIso];
	float& subLeading_photon_trkIso = properties_subLeadingPhoton[photonProperty::trkIso];

	if (subLeading_photon_R9 >= 0.85) {
	  subLeadingPassesEmulation = ((subLeading_photon_pT > 18.) &&
				       (subLeading_photon_HOverE <= 0.1*subLeading_photon_energy));
	}
	else if ((subLeading_photon_R9 >= 0.5) && (subLeading_photon_R9 < 0.85)) {
	  subLeadingPassesEmulation = ((subLeading_photon_pT > 18.) &&
				       (subLeading_photon_HOverE <= 0.1*subLeading_photon_energy) &&
				       (subLeading_photon_sigmaIEtaIEta <= 0.015) &&
				       (subLeading_photon_ecalClusIso <= (6.0 + 0.012*subLeading_photon_pT)) &&
				       (subLeading_photon_trkIso <= (6.0 + 0.002*subLeading_photon_pT)));
	}

	passesEmulation = leadingPassesEmulation && subLeadingPassesEmulation;
      }
      else if (triggerBit == 22) {
	bool leadingPassesEmulation = false;
	float& leading_photon_pT = properties_leadingPhoton[photonProperty::pT];
	float& leading_photon_HOverE = properties_leadingPhoton[photonProperty::hOverE];

	leadingPassesEmulation = ((leading_photon_pT > 70.) &&
				  (leading_photon_HOverE <= 0.15));

	bool subLeadingPassesEmulation = false;
	float& subLeading_photon_pT = properties_subLeadingPhoton[photonProperty::pT];
	float& subLeading_photon_HOverE = properties_subLeadingPhoton[photonProperty::hOverE];

	subLeadingPassesEmulation = ((subLeading_photon_pT > 70.) &&
				     (subLeading_photon_HOverE) <= 0.15);

	passesEmulation = leadingPassesEmulation && subLeadingPassesEmulation;
      }
      else {
	std::cout << "ERROR: unsupported trigger bit: " << triggerBit << std::endl;
	std::exit(EXIT_FAILURE);
      }
      break;
    }
  default:
    {
      std::cout << "ERROR: unsupported year: " << year << std::endl;
      std::exit(EXIT_FAILURE);
    }
  }
  return passesEmulation;
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

  // Quality cuts
  photonQualityCutsStruct* qualityCuts = &(parameters.photonQualityCutsBarrel);
  if (absEta > parameters.photonBarrelEtaCut) qualityCuts = &(parameters.photonQualityCutsEndcap);

  properties[photonProperty::hOverE] = (photonsCollection.HOverE)->at(photonIndex);
  bool passesHOverE = (properties[photonProperty::hOverE] < qualityCuts->towerHOverE);
  medium_bits[mediumPhotonCriterion::hOverE] = passesHOverE;
  // bool passesHOverELoose = (properties[photonProperty::hOverE] < qualityCuts->towerHOverELoose);

  float pTDependentNeutralIsolationCut = (qualityCuts->neutralIsolation).getPolynomialValue(properties[photonProperty::pT]);
  properties[photonProperty::rhoCorrectedNeutralIsolation] = getRhoCorrectedIsolation(((photonsCollection.PFNeutralIsolationUncorrected)->at(photonIndex)), PFTypesForEA::neutralHadron, absEta, rho, parameters.effectiveAreas);
  bool passesNeutralIsolation = (properties[photonProperty::rhoCorrectedNeutralIsolation] < pTDependentNeutralIsolationCut);
  medium_bits[mediumPhotonCriterion::neutralIsolation] = passesNeutralIsolation;
  // float pTDependentNeutralIsolationCutLoose = (qualityCuts->neutralIsolationLoose).getPolynomialValue(properties[photonProperty::pT]);
  // bool passesNeutralIsolationLoose = (properties[photonProperty::rhoCorrectedNeutralIsolation] < pTDependentNeutralIsolationCutLoose);

  float pTDependentPhotonIsolationCut = (qualityCuts->photonIsolation).getPolynomialValue(properties[photonProperty::pT]);
  properties[photonProperty::rhoCorrectedPhotonIsolation] = getRhoCorrectedIsolation(((photonsCollection.PFPhotonIsolationUncorrected)->at(photonIndex)), PFTypesForEA::photon, absEta, rho, parameters.effectiveAreas);
  bool passesPhotonIsolation = (properties[photonProperty::rhoCorrectedPhotonIsolation] < pTDependentPhotonIsolationCut);
  medium_bits[mediumPhotonCriterion::photonIsolation] = passesPhotonIsolation;
  // float pTDependentPhotonIsolationCutLoose = (qualityCuts->photonIsolationLoose).getPolynomialValue(properties[photonProperty::pT]);
  // fake_bits[fakePhotonCriterion::passesPhoIsoVeto] = (properties[photonProperty::rhoCorrectedPhotonIsolation] >= (parameters.looseIDFactor)*(pTDependentPhotonIsolationCutLoose));

  properties[photonProperty::rawChargedIsolation] = (photonsCollection.PFChargedIsolationUncorrected)->at(photonIndex);
  properties[photonProperty::rhoCorrectedChargedIsolation] = getRhoCorrectedIsolation(((photonsCollection.PFChargedIsolationUncorrected)->at(photonIndex)), PFTypesForEA::chargedHadron, absEta, rho, parameters.effectiveAreas);
  bool passesChargedIsolation = (properties[photonProperty::rhoCorrectedChargedIsolation] < qualityCuts->chargedIsolation);
  medium_bits[mediumPhotonCriterion::chargedIsolation] = passesChargedIsolation;
  fake_bits[fakePhotonCriterion::passesChIsoVeto] = (properties[photonProperty::rhoCorrectedChargedIsolation] >= qualityCuts->chargedIsolationLoose);

  properties[photonProperty::sigmaIEtaIEta] = ((photonsCollection.sigmaIEtaIEta)->at(photonIndex));
  bool passesSigmaIEtaIEta = (properties[photonProperty::sigmaIEtaIEta] < qualityCuts->sigmaIEtaIEta);
  medium_bits[mediumPhotonCriterion::sigmaIEtaIEta] = passesSigmaIEtaIEta;
  // bool passesSigmaIEtaIEtaLoose = (properties[photonProperty::sigmaIEtaIEta] < qualityCuts->sigmaIEtaIEtaLoose);

  // fake_bits[fakePhotonCriterion::passesOtherLooseCuts] = (passesHOverELoose || passesSigmaIEtaIEtaLoose || passesNeutralIsolationLoose);

  fake_bits[fakePhotonCriterion::passesShowerShapeMedIDCuts] = ((properties[photonProperty::sigmaIEtaIEta] < qualityCuts->sigmaIEtaIEta) && (properties[photonProperty::hOverE] < qualityCuts->towerHOverE));

  properties[photonProperty::R9] = ((photonsCollection.R9)->at(photonIndex));
  properties[photonProperty::ecalClusIso] = ((photonsCollection.ecalClusIso)->at(photonIndex));
  properties[photonProperty::trkIso] = ((photonsCollection.trkIso)->at(photonIndex));
  properties[photonProperty::energy] = (photonsCollection.energy)->at(photonIndex);

  assert(static_cast<int>(medium_bits.size()) == static_cast<int>(mediumPhotonCriterion::nMediumPhotonCriteria));
  int nFalseBits_medium = getNFalseBits(medium_bits);
  fake_bits[fakePhotonCriterion::failsMediumID] = (nFalseBits_medium >= 1);

  assert(static_cast<int>(fake_bits.size()) == static_cast<int>(fakePhotonCriterion::nFakePhotonCriteria));
  int nFalseBits_fake = getNFalseBits(fake_bits);
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

  // if ((results.isMarginallyUnselectedMedium || results.isMarginallyUnselectedFake) || (results.isSelectedMedium || results.isSelectedFake)) assert(static_cast<int>((results.pho_properties).size()) == static_cast<int>(photonProperty::nPhotonProperties)); // delay until later, when the angle to the closest jet is filled

  return results;
}

MCExaminationResultsStruct examineMCParticle(parametersStruct &parameters, const MCCollectionStruct& MCCollection, const int& MCIndex) {
  MCExaminationResultsStruct MCExaminationResults;
  int particle_mcPID = (MCCollection.MCPIDs)->at(MCIndex);
  // int custom_particle_ID = PIDUtils::getCustomParticleID(particle_mcPID);
  int particle_mcMomPID = (MCCollection.MCMomPIDs)->at(MCIndex);
  int custom_mom_ID = PIDUtils::getCustomParticleID(particle_mcMomPID);
  UShort_t particle_statusFlag = static_cast<UShort_t>(((MCCollection.MCStatusFlags)->at(MCIndex))&(static_cast<UShort_t>(7u))); // picks out only first 3 bits

  if (PIDUtils::isPhotonPID(particle_mcPID) && (PIDUtils::isNeutralinoPID(particle_mcMomPID) || PIDUtils::isHiggsPID(particle_mcMomPID))) {
    if (passesBitMask(particle_statusFlag, parameters.MCStatusFlagBitMask)) {
      MCExaminationResults.isPhotonWithDesiredMom = true;
      truthPhotonProperties& pho_properties = MCExaminationResults.truth_photon_properties;
      pho_properties[truthPhotonProperty::eta] = (MCCollection.MCEtas)->at(MCIndex);
      pho_properties[truthPhotonProperty::phi] = (MCCollection.MCPhis)->at(MCIndex);
      pho_properties[truthPhotonProperty::pT] = (MCCollection.MCEts)->at(MCIndex);
      pho_properties[truthPhotonProperty::status] = (MCCollection.MCStatuses)->at(MCIndex);
      // assert(static_cast<int>((MCExaminationResults.truth_photon_properties).size()) == static_cast<int>(truthPhotonProperty::nTruthPhotonProperties)); // distance to nearest truth jet candidate needs to be set later, do this check then
      if (PIDUtils::isHiggsPID(particle_mcMomPID)) { // fill in "fake" gluino and neutralino masses for Hgg events
        MCExaminationResults.gluinoMass = 1500.;
        MCExaminationResults.neutralinoMass = 800.;
      }
    }
  }
  // if (PIDUtils::isNeutralinoPID(particle_mcMomPID)) std::cout << std::endl << "Found MC particle with neutralino mom: ID = " << particle_mcPID << std::endl;
  // if (PIDUtils::isSinglinoPID(particle_mcMomPID)) std::cout << std::endl << "Found MC particle with singlino mom. Its ID = " << particle_mcPID << std::endl;
  // if (PIDUtils::isSingletPID(particle_mcMomPID)) {
  //   std::cout << std::endl << "Found unexpected MC particle with singlet mom. Its ID = " << particle_mcPID << std::endl;
  //   std::exit(EXIT_FAILURE);
  // }
  // if (PIDUtils::isGluinoPID(particle_mcMomPID)) std::cout << std::endl << "Found MC particle with gluino mom: ID = " << particle_mcPID << std::endl;
  // if (PIDUtils::isNeutralinoPID(particle_mcPID)) std::cout << std::endl << "Found MC neutralino. Its mom has: ID = " << particle_mcMomPID << std::endl;
  if (PIDUtils::isJetCandidatePID(particle_mcPID) && ((MCCollection.MCStatuses)->at(MCIndex) == parameters.jetCandidateStatusConstraint)) {
    if (PIDUtils::isGluinoPID(particle_mcMomPID)) MCExaminationResults.isJetCandidateFromGluino = true;
    if (PIDUtils::isSingletPID(particle_mcMomPID)) MCExaminationResults.isJetCandidateFromSinglet = true;
    if (MCExaminationResults.isJetCandidateFromGluino || MCExaminationResults.isJetCandidateFromSinglet) {
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
  if (PIDUtils::isGluinoPID(particle_mcPID)) MCExaminationResults.gluinoMass = (MCCollection.MCMasses)->at(MCIndex);
  if (PIDUtils::isNeutralinoPID(particle_mcMomPID)) MCExaminationResults.neutralinoMass = (MCCollection.MCMomMasses)->at(MCIndex);

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

jetExaminationResultsStruct examineJet(optionsStruct &options, parametersStruct &parameters, const jetsCollectionStruct& jetsCollection, const int& jetIndex, std::vector<angularVariablesStruct> &selectedPhotonAngles, std::vector<angularVariablesStruct> &selectedTruePhotonAngles, std::vector<angularVariablesStruct> &selectedTrueJetCandidateAngles_all) {
  jetExaminationResultsStruct results;

  if (options.isMC) {
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
      results.hasGluinoPartonMom = PIDUtils::isGluinoPID((jetsCollection.jetGenPartonMomID)->at(jetIndex));
      results.hasSingletPartonMom = PIDUtils::isSingletPID((jetsCollection.jetGenPartonMomID)->at(jetIndex));
    }
  }

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
  bits[jetCriterion::deltaR_photon] = ((minDeltaR > parameters.deltaRScale_jetPhotonDistance) || (minDeltaR < 0.));
  results.isAwayFromCaloPhoton = bits[jetCriterion::deltaR_photon];
  properties[jetProperty::deltaR_nearestTruePhoton] = -0.005;
  properties[jetProperty::deltaR_nearestTrueJetCandidate] = -0.005;
  properties[jetProperty::truthPTRatio] = -1.0;
  properties[jetProperty::deltaR_genJet] = -0.005;
  if (options.isMC) {
    properties[jetProperty::deltaR_nearestTruePhoton] = jetAngle.getMinDeltaR(selectedTruePhotonAngles);
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
  if (options.isMC) {
    jecFractionalUncertainty = (jetsCollection.JECUncertainty)->at(jetIndex);
    float jet_pT_JECDown = (1.0 - jecFractionalUncertainty)*jet_pT;
    float jet_pT_JECUp = (1.0 + jecFractionalUncertainty)*jet_pT;
    passesPT_JECDown = (jet_pT_JECDown > parameters.jetpTCut);
    passesPT_JECUp = (jet_pT_JECUp > parameters.jetpTCut);
  }

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

void setSelectedPhotonClosestJet(photonPropertiesCollection& photon_properties_collection, std::vector<angularVariablesStruct>& genJetAngles, std::vector<angularVariablesStruct>& gluinoMomGenJetAngles, std::vector<angularVariablesStruct>& singletMomGenJetAngles) {
  for (auto&& photon_properties: photon_properties_collection) {
    angularVariablesStruct photonAngle = angularVariablesStruct(photon_properties.at(photonProperty::eta), photon_properties.at(photonProperty::phi));
    photon_properties[photonProperty::deltaR_nearestGenJet] = photonAngle.getMinDeltaR(genJetAngles);
    photon_properties[photonProperty::deltaR_nearestGluinoMomGenJet] = photonAngle.getMinDeltaR(gluinoMomGenJetAngles);
    photon_properties[photonProperty::deltaR_nearestSingletMomGenJet] = photonAngle.getMinDeltaR(singletMomGenJetAngles);
    assert(static_cast<int>((photon_properties).size()) == static_cast<int>(photonProperty::nPhotonProperties));
  }
}

void setUnselectedFakePhotonClosestJet(unselectedFakePhotonPropertiesCollection& unselected_photon_properties_collection, std::vector<angularVariablesStruct>& genJetAngles, std::vector<angularVariablesStruct>& gluinoMomGenJetAngles, std::vector<angularVariablesStruct>& singletMomGenJetAngles) {
  for (auto&& unselected_photon_properties_pair: unselected_photon_properties_collection) {
    photonProperties& photon_properties = unselected_photon_properties_pair.second;
    angularVariablesStruct photonAngle = angularVariablesStruct(photon_properties.at(photonProperty::eta), photon_properties.at(photonProperty::phi));
    photon_properties[photonProperty::deltaR_nearestGenJet] = photonAngle.getMinDeltaR(genJetAngles);
    photon_properties[photonProperty::deltaR_nearestGluinoMomGenJet] = photonAngle.getMinDeltaR(gluinoMomGenJetAngles);
    photon_properties[photonProperty::deltaR_nearestSingletMomGenJet] = photonAngle.getMinDeltaR(singletMomGenJetAngles);
    assert(static_cast<int>((photon_properties).size()) == static_cast<int>(photonProperty::nPhotonProperties));
  }
}

void setUnselectedMediumPhotonClosestJet(unselectedMediumPhotonPropertiesCollection& unselected_photon_properties_collection, std::vector<angularVariablesStruct>& genJetAngles, std::vector<angularVariablesStruct>& gluinoMomGenJetAngles, std::vector<angularVariablesStruct>& singletMomGenJetAngles) {
  for (auto&& unselected_photon_properties_pair: unselected_photon_properties_collection) {
    photonProperties& photon_properties = unselected_photon_properties_pair.second;
    angularVariablesStruct photonAngle = angularVariablesStruct(photon_properties.at(photonProperty::eta), photon_properties.at(photonProperty::phi));
    photon_properties[photonProperty::deltaR_nearestGenJet] = photonAngle.getMinDeltaR(genJetAngles);
    photon_properties[photonProperty::deltaR_nearestGluinoMomGenJet] = photonAngle.getMinDeltaR(gluinoMomGenJetAngles);
    photon_properties[photonProperty::deltaR_nearestSingletMomGenJet] = photonAngle.getMinDeltaR(singletMomGenJetAngles);
    assert(static_cast<int>((photon_properties).size()) == static_cast<int>(photonProperty::nPhotonProperties));
  }
}

eventExaminationResultsStruct examineEvent(optionsStruct &options, parametersStruct &parameters, Long64_t& entryIndex, const int& year, eventDetailsStruct& eventDetails, MCCollectionStruct &MCCollection, photonsCollectionStruct &photonsCollection, jetsCollectionStruct &jetsCollection, statisticsHistograms& statistics, STRegionsStruct& STRegions) {
  eventExaminationResultsStruct eventResult;

  eventResult.eventIndex = entryIndex;
  selectionRegion& region = eventResult.evt_region;
  std::map<eventSelectionCriterion, bool> selectionBits;
  eventProperties event_properties;
  float& event_ST = eventResult.evt_ST;
  int n_goodJetsCloseToSelectedPhoton = 0;
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
  int nPhotonsWithDesiredMom = 0;
  truthPhotonPropertiesCollection selectedTruePhotonProperties;
  std::vector<angularVariablesStruct> selectedTruePhotonAngles;
  int nJetCandidatesWithStealthMom = 0;
  int nJetCandidatesWithGluinoMom = 0;
  int nJetCandidatesWithSingletMom = 0;
  int nStealthJetsCloseToTruePhoton = 0;
  truthJetCandidatePropertiesCollection selectedTrueJetCandidateProperties_all;
  std::vector<angularVariablesStruct> selectedTrueJetCandidateAngles_all;
  truthJetCandidatePropertiesCollection selectedTrueJetCandidateProperties_fromGluino;
  std::vector<angularVariablesStruct> selectedTrueJetCandidateAngles_fromGluino;
  truthJetCandidatePropertiesCollection selectedTrueJetCandidateProperties_fromSinglet;
  std::vector<angularVariablesStruct> selectedTrueJetCandidateAngles_fromSinglet; // wasteful, fix later...
  int MCRegionIndex = 0;
  if (options.isMC) {
    bool gluinoMassIsSet = false;
    bool neutralinoMassIsSet = false;
    for (int MCIndex = 0; MCIndex < eventDetails.nMCParticles; ++MCIndex) {
      MCExaminationResultsStruct MCExaminationResults = examineMCParticle(parameters, MCCollection, MCIndex);
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
        if (MCExaminationResults.isJetCandidateFromGluino) {
          ++nJetCandidatesWithGluinoMom;
          selectedTrueJetCandidateProperties_fromGluino.push_back(MCExaminationResults.truth_jetCandidate_properties);
          selectedTrueJetCandidateAngles_fromGluino.push_back(angularVariablesStruct((MCExaminationResults.truth_jetCandidate_properties)[truthJetCandidateProperty::eta], (MCExaminationResults.truth_jetCandidate_properties)[truthJetCandidateProperty::phi]));
        }
        else if (MCExaminationResults.isJetCandidateFromSinglet) {
          ++nJetCandidatesWithSingletMom;
          selectedTrueJetCandidateProperties_fromSinglet.push_back(MCExaminationResults.truth_jetCandidate_properties);
          selectedTrueJetCandidateAngles_fromSinglet.push_back(angularVariablesStruct((MCExaminationResults.truth_jetCandidate_properties)[truthJetCandidateProperty::eta], (MCExaminationResults.truth_jetCandidate_properties)[truthJetCandidateProperty::phi]));
        }
      }
      if ((MCExaminationResults.gluinoMass > 0.) && !(gluinoMassIsSet)) {
        generated_gluinoMass = MCExaminationResults.gluinoMass;
        gluinoMassIsSet = true;
      }
      if ((MCExaminationResults.neutralinoMass > 0.) && !(neutralinoMassIsSet)) {
        generated_neutralinoMass = MCExaminationResults.neutralinoMass;
        neutralinoMassIsSet = true;
      }
    } // ends loop over MC particles
    selectionBits[eventSelectionCriterion::MCGenInformation] = (nPhotonsWithDesiredMom == 2);
    if (selectionBits[eventSelectionCriterion::MCGenInformation] && (!(gluinoMassIsSet && neutralinoMassIsSet))) {
      std::cout << "ERROR: Unable to find gluino or neutralino mass in an event that passes MC selection." << std::endl;
      std::exit(EXIT_FAILURE);
    }
    MCRegionIndex = MCRegions::getRegionIndex(generated_gluinoMass, generated_neutralinoMass);
  }
  event_properties[eventProperty::MC_nPhotonsWithDesiredMom] = nPhotonsWithDesiredMom;
  event_properties[eventProperty::MC_nJetCandidatesWithStealthMom] = nJetCandidatesWithStealthMom;
  event_properties[eventProperty::MC_nJetCandidatesWithGluinoMom] = nJetCandidatesWithGluinoMom;
  event_properties[eventProperty::MC_nJetCandidatesWithSingletMom] = nJetCandidatesWithSingletMom;
  event_properties[eventProperty::MC_nStealthJetsCloseToTruePhoton] = nStealthJetsCloseToTruePhoton;

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
  int n_truthMatchedMediumPhotons = 0;
  int n_fakePhotons = 0;
  int n_truthMatchedFakePhotons = 0;
  float eta_leadingPhoton = -1000.;
  float pT_leadingPhoton = -1.;
  float eta_subLeadingPhoton = -1000.;
  float pT_subLeadingPhoton = -1.;
  photonProperties properties_leadingPhoton;
  photonProperties properties_subLeadingPhoton;
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
      if ((eta_leadingPhoton < -999.) && (pT_leadingPhoton < 0.)) {
        eta_leadingPhoton = (photonExaminationResults.pho_properties)[photonProperty::eta];
        pT_leadingPhoton = (photonExaminationResults.pho_properties)[photonProperty::pT];
        properties_leadingPhoton = (photonExaminationResults.pho_properties);
      }
      else {
        if (((photonExaminationResults.pho_properties)[photonProperty::pT]) > pT_leadingPhoton) {
          pT_subLeadingPhoton = pT_leadingPhoton;
          eta_subLeadingPhoton = eta_leadingPhoton;
          properties_subLeadingPhoton = properties_leadingPhoton;
          pT_leadingPhoton = (photonExaminationResults.pho_properties)[photonProperty::pT];
          eta_leadingPhoton = (photonExaminationResults.pho_properties)[photonProperty::eta];
          properties_leadingPhoton = (photonExaminationResults.pho_properties);
        }
        else {
          pT_subLeadingPhoton = (photonExaminationResults.pho_properties)[photonProperty::pT];
          eta_subLeadingPhoton = (photonExaminationResults.pho_properties)[photonProperty::eta];
          properties_subLeadingPhoton = (photonExaminationResults.pho_properties);
        }
      }
    }
    if (photonExaminationResults.isSelectedMedium) {
      ++n_mediumPhotons;
      selectedMediumPhotonProperties.push_back(photonExaminationResults.pho_properties);
      if (options.isMC) {
        float nearestTruePhotonDeltaR = (photonExaminationResults.pho_properties)[photonProperty::deltaR_nearestTruePhoton];
        if (nearestTruePhotonDeltaR >= parameters.deltaRScale_truthMatching) selectedMediumPhotonProperties_awayFromTruePhoton.push_back(photonExaminationResults.pho_properties);
        else if (nearestTruePhotonDeltaR > 0.) {
          ++n_truthMatchedMediumPhotons;
          selectedMediumPhotonProperties_closeToTruePhoton.push_back(photonExaminationResults.pho_properties);
        }
      }
    }
    else if (photonExaminationResults.isSelectedFake) {
      ++n_fakePhotons;
      selectedFakePhotonProperties.push_back(photonExaminationResults.pho_properties);
      if (options.isMC) {
        float nearestTruePhotonDeltaR = (photonExaminationResults.pho_properties)[photonProperty::deltaR_nearestTruePhoton];
        if (nearestTruePhotonDeltaR >= parameters.deltaRScale_truthMatching) selectedFakePhotonProperties_awayFromTruePhoton.push_back(photonExaminationResults.pho_properties);
        else if (nearestTruePhotonDeltaR > 0.) {
          ++n_truthMatchedFakePhotons;
          selectedFakePhotonProperties_closeToTruePhoton.push_back(photonExaminationResults.pho_properties);
        }
      }
    }
    if (photonExaminationResults.isMarginallyUnselectedFake) {
      unselected_fake_pho_properties.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedFakeCriterion, photonExaminationResults.pho_properties));
      if (options.isMC) {
        float nearestTruePhotonDeltaR = (photonExaminationResults.pho_properties)[photonProperty::deltaR_nearestTruePhoton];
        if (nearestTruePhotonDeltaR >= parameters.deltaRScale_truthMatching) unselected_fake_pho_properties_awayFromTruePhoton.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedFakeCriterion, photonExaminationResults.pho_properties));
        else if (nearestTruePhotonDeltaR > 0.) unselected_fake_pho_properties_closeToTruePhoton.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedFakeCriterion, photonExaminationResults.pho_properties));
      }
    }
    if (photonExaminationResults.isMarginallyUnselectedMedium) {
      unselected_medium_pho_properties.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedMediumCriterion, photonExaminationResults.pho_properties));
      if (options.isMC) {
        float nearestTruePhotonDeltaR = (photonExaminationResults.pho_properties)[photonProperty::deltaR_nearestTruePhoton];
        if (nearestTruePhotonDeltaR >= parameters.deltaRScale_truthMatching) unselected_medium_pho_properties_awayFromTruePhoton.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedMediumCriterion, photonExaminationResults.pho_properties));
        else if (nearestTruePhotonDeltaR > 0.) unselected_medium_pho_properties_closeToTruePhoton.push_back(std::make_pair(photonExaminationResults.marginallyUnselectedMediumCriterion, photonExaminationResults.pho_properties));
      }
    }
  }

  bool passes_HLTEmulation = passesHLTEmulation(year, properties_leadingPhoton, properties_subLeadingPhoton, parameters.HLTPhotonBit);
  event_properties[eventProperty::MC_nTruthMatchedMediumPhotons] = n_truthMatchedMediumPhotons;
  event_properties[eventProperty::MC_nTruthMatchedFakePhotons] = n_truthMatchedFakePhotons;

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
  genJetPropertiesCollection gen_jet_properties_collection;
  genJetPropertiesCollection gluino_mom_gen_jet_properties_collection;
  genJetPropertiesCollection singlet_mom_gen_jet_properties_collection;
  std::vector<angularVariablesStruct> genJetAngles;
  std::vector<angularVariablesStruct> gluinoMomGenJetAngles;
  std::vector<angularVariablesStruct> singletMomGenJetAngles;
  int n_genJets = 0;
  int n_gluinoMomGenJets = 0;
  int n_singletMomGenJets = 0;

  for (Int_t jetIndex = 0; jetIndex < (eventDetails.nJets); ++jetIndex) {
    jetExaminationResultsStruct jetExaminationResults = examineJet(options, parameters, jetsCollection, jetIndex, selectedPhotonAngles, selectedTruePhotonAngles, selectedTrueJetCandidateAngles_all);
    // if (options.isMC) counters.jetTotalCountersMCMap->Fill(generated_gluinoMass, generated_neutralinoMass);
    if (options.isMC && jetExaminationResults.hasGenVariablesSet) {
      if ((jetExaminationResults.gen_jet_properties)[genJetProperty::pT] > parameters.pTCutSubLeading) {
        ++n_genJets;
        gen_jet_properties_collection.push_back(jetExaminationResults.gen_jet_properties);
        genJetAngles.push_back(angularVariablesStruct((jetExaminationResults.gen_jet_properties)[genJetProperty::eta], (jetExaminationResults.gen_jet_properties)[genJetProperty::phi]));
        if (jetExaminationResults.hasGluinoPartonMom) {
          ++n_gluinoMomGenJets;
          gluino_mom_gen_jet_properties_collection.push_back(jetExaminationResults.gen_jet_properties);
          gluinoMomGenJetAngles.push_back(angularVariablesStruct((jetExaminationResults.gen_jet_properties)[genJetProperty::eta], (jetExaminationResults.gen_jet_properties)[genJetProperty::phi]));
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
    if (jetExaminationResults.passesSelectionJECNominal) selectedJetProperties.push_back(jetExaminationResults.jet_properties);
    else if (jetExaminationResults.isMarginallyUnselected) {
      unselected_jet_properties.push_back(std::make_pair(jetExaminationResults.marginallyUnselectedCriterion, jetExaminationResults.jet_properties));
      if (options.isMC) {
        float nearestTruePhotonDeltaR = (jetExaminationResults.jet_properties)[jetProperty::deltaR_nearestTruePhoton];
        if (nearestTruePhotonDeltaR >= parameters.deltaRScale_truthMatching) unselected_jet_properties_awayFromTruePhoton.push_back(std::make_pair(jetExaminationResults.marginallyUnselectedCriterion, jetExaminationResults.jet_properties));
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
          if (nearestTruePhotonDeltaR >= parameters.deltaRScale_truthMatching) selectedJetProperties_awayFromTruePhoton.push_back(jetExaminationResults.jet_properties);
          else if (nearestTruePhotonDeltaR > 0.) selectedJetProperties_closeToTruePhoton.push_back(jetExaminationResults.jet_properties);
        }
      }
      else {
        ++n_goodJetsCloseToSelectedPhoton;
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
  event_properties[eventProperty::nGoodJetsCloseToSelectedPhoton] = n_goodJetsCloseToSelectedPhoton;
  event_properties[eventProperty::nJetsDR] = n_jetsDR;
  int max_nJets = n_jetsDR;
  if (options.isMC) { // this makes sure that the nJets used to make the decision whether or not to save the event is the maximum nJets accounting for all the shifts
    int maxNJetsShifted = getMaxNJets(shifted_nJetsDR);
    if (maxNJetsShifted > max_nJets) max_nJets = maxNJetsShifted;
  }

  event_properties[eventProperty::MC_nGenJets] = n_genJets;
  event_properties[eventProperty::MC_nGluinoMomGenJets] = n_gluinoMomGenJets;
  event_properties[eventProperty::MC_nSingletMomGenJets] = n_singletMomGenJets;
  if (options.isMC) {
    assert(static_cast<int>(selectedTruePhotonProperties.size()) == nPhotonsWithDesiredMom);
    assert(static_cast<int>(selectedTruePhotonAngles.size()) == nPhotonsWithDesiredMom);
    for (int truePhotonCounter = 0; truePhotonCounter < nPhotonsWithDesiredMom; ++truePhotonCounter) {
      (selectedTruePhotonProperties[truePhotonCounter])[truthPhotonProperty::deltaR_nearestTruthJetCandidate] = (selectedTruePhotonAngles[truePhotonCounter]).getMinDeltaR(selectedTrueJetCandidateAngles_all);
      (selectedTruePhotonProperties[truePhotonCounter])[truthPhotonProperty::deltaR_nearestGenJet] = (selectedTruePhotonAngles[truePhotonCounter]).getMinDeltaR(genJetAngles);
      (selectedTruePhotonProperties[truePhotonCounter])[truthPhotonProperty::deltaR_nearestGluinoMomGenJet] = (selectedTruePhotonAngles[truePhotonCounter]).getMinDeltaR(gluinoMomGenJetAngles);
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
    assert(static_cast<int>(selectedTrueJetCandidateProperties_fromGluino.size()) == nJetCandidatesWithGluinoMom);
    assert(static_cast<int>(selectedTrueJetCandidateAngles_fromGluino.size()) == nJetCandidatesWithGluinoMom);
    for (int trueJetCandidateCounter = 0; trueJetCandidateCounter < nJetCandidatesWithGluinoMom; ++trueJetCandidateCounter) {
      (selectedTrueJetCandidateProperties_fromGluino[trueJetCandidateCounter])[truthJetCandidateProperty::deltaR_nearestTruePhoton] = (selectedTrueJetCandidateAngles_fromGluino[trueJetCandidateCounter]).getMinDeltaR(selectedTruePhotonAngles);
      assert(static_cast<int>((selectedTrueJetCandidateProperties_fromGluino[trueJetCandidateCounter]).size()) == static_cast<int>(truthJetCandidateProperty::nTruthJetCandidateProperties));
    }
    assert(static_cast<int>(selectedTrueJetCandidateProperties_fromSinglet.size()) == nJetCandidatesWithSingletMom);
    assert(static_cast<int>(selectedTrueJetCandidateAngles_fromSinglet.size()) == nJetCandidatesWithSingletMom);
    for (int trueJetCandidateCounter = 0; trueJetCandidateCounter < nJetCandidatesWithSingletMom; ++trueJetCandidateCounter) {
      (selectedTrueJetCandidateProperties_fromSinglet[trueJetCandidateCounter])[truthJetCandidateProperty::deltaR_nearestTruePhoton] = (selectedTrueJetCandidateAngles_fromSinglet[trueJetCandidateCounter]).getMinDeltaR(selectedTruePhotonAngles);
      assert(static_cast<int>((selectedTrueJetCandidateProperties_fromSinglet[trueJetCandidateCounter]).size()) == static_cast<int>(truthJetCandidateProperty::nTruthJetCandidateProperties));
    }
    if (n_genJets > 0) {
      setSelectedPhotonClosestJet(selectedMediumPhotonProperties, genJetAngles, gluinoMomGenJetAngles, singletMomGenJetAngles);
      setSelectedPhotonClosestJet(selectedMediumPhotonProperties_closeToTruePhoton, genJetAngles, gluinoMomGenJetAngles, singletMomGenJetAngles);
      setSelectedPhotonClosestJet(selectedMediumPhotonProperties_awayFromTruePhoton, genJetAngles, gluinoMomGenJetAngles, singletMomGenJetAngles);
      setSelectedPhotonClosestJet(selectedFakePhotonProperties, genJetAngles, gluinoMomGenJetAngles, singletMomGenJetAngles);
      setSelectedPhotonClosestJet(selectedFakePhotonProperties_closeToTruePhoton, genJetAngles, gluinoMomGenJetAngles, singletMomGenJetAngles);
      setSelectedPhotonClosestJet(selectedFakePhotonProperties_awayFromTruePhoton, genJetAngles, gluinoMomGenJetAngles, singletMomGenJetAngles);

      setUnselectedFakePhotonClosestJet(unselected_fake_pho_properties, genJetAngles, gluinoMomGenJetAngles, singletMomGenJetAngles);
      setUnselectedFakePhotonClosestJet(unselected_fake_pho_properties_closeToTruePhoton, genJetAngles, gluinoMomGenJetAngles, singletMomGenJetAngles);
      setUnselectedFakePhotonClosestJet(unselected_fake_pho_properties_awayFromTruePhoton, genJetAngles, gluinoMomGenJetAngles, singletMomGenJetAngles);
      setUnselectedMediumPhotonClosestJet(unselected_medium_pho_properties, genJetAngles, gluinoMomGenJetAngles, singletMomGenJetAngles);
      setUnselectedMediumPhotonClosestJet(unselected_medium_pho_properties_closeToTruePhoton, genJetAngles, gluinoMomGenJetAngles, singletMomGenJetAngles);
      setUnselectedMediumPhotonClosestJet(unselected_medium_pho_properties_awayFromTruePhoton, genJetAngles, gluinoMomGenJetAngles, singletMomGenJetAngles);
    }
  }

  selectionBits[eventSelectionCriterion::NJets] = (max_nJets >= 2);
  if (options.disableJetSelection) selectionBits[eventSelectionCriterion::NJets] = true;
  // Add MET to ST
  event_ST += eventDetails.PFMET;
  event_properties[eventProperty::ST] = event_ST;
  event_properties[eventProperty::selectionRegionIndex] = 1.0*static_cast<int>(region);

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

  statistics.fillIDEfficiencyStatisticsHistograms(event_ST, region);

  assert(static_cast<int>(selectionBits.size()) == static_cast<int>(eventSelectionCriterion::nEventSelectionCriteria));
  int nEventFalseBits = getNFalseBits(selectionBits);
  eventProperties temp1 = initialize_eventProperties_with_defaults(); // temp1 and temp2 are dummies -- they won't contribute to the histograms
  if (nEventFalseBits == 0) {
    unselectedEventProperties temp2 = std::make_pair(eventSelectionCriterion::nEventSelectionCriteria, temp1);
    statistics.fill1DStatisticsHistograms(event_properties, false, temp2,
                                          selectedTruePhotonProperties,
                                          selectedTrueJetCandidateProperties_all,
                                          selectedTrueJetCandidateProperties_fromGluino,
                                          selectedTrueJetCandidateProperties_fromSinglet,
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
                                          gen_jet_properties_collection,
                                          gluino_mom_gen_jet_properties_collection,
                                          singlet_mom_gen_jet_properties_collection,
                                          region, options.isMC, MCRegionIndex);
    statistics.fillHLTEmulationStatisticsHistograms(eta_leadingPhoton, pT_leadingPhoton,
                                                    eta_subLeadingPhoton, pT_subLeadingPhoton,
                                                    passes_HLTEmulation, region);
  }
  else if (nEventFalseBits == 1) {
    eventSelectionCriterion marginallyUnselectedEventCriterion = getFirstFalseCriterion(selectionBits);
    unselectedEventProperties unselected_event_properties = std::make_pair(marginallyUnselectedEventCriterion, event_properties);
    statistics.fill1DStatisticsHistograms(temp1, true, unselected_event_properties,
                                          selectedTruePhotonProperties,
                                          selectedTrueJetCandidateProperties_all,
                                          selectedTrueJetCandidateProperties_fromGluino,
                                          selectedTrueJetCandidateProperties_fromSinglet,
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
                                          gen_jet_properties_collection,
                                          gluino_mom_gen_jet_properties_collection,
                                          singlet_mom_gen_jet_properties_collection,
                                          region, options.isMC, MCRegionIndex);
  }

  eventResult.isInterestingEvent = ((nEventFalseBits == 0) && (event_ST >= (STRegions.STNormRangeMin - parameters.preNormalizationBuffer)));

  eventResult.evt_photonPT_leading = pT_leadingPhoton;
  eventResult.evt_photonPT_subLeading = pT_subLeadingPhoton;
  eventResult.evt_photonEta_leading = eta_leadingPhoton;
  eventResult.evt_photonEta_subLeading = eta_subLeadingPhoton;

  if (nEventFalseBits <= 1) assert(static_cast<int>(event_properties.size()) == static_cast<int>(eventProperty::nEventProperties));
  return eventResult;
}

void loopOverEvents(optionsStruct &options, parametersStruct &parameters, const int& year, std::vector<eventExaminationResultsStruct>& selectedEventsInfo, statisticsHistograms& statistics, STRegionsStruct& STRegions) {
  TChain inputChain("ggNtuplizer/EventTree");
  std::cout << "Starting to add files to chain..." << std::endl;
  for (const std::string& inputPath: (options.inputPaths)) {
    std::cout << "Adding... " << inputPath << std::endl;
    int read_status = inputChain.Add(inputPath.c_str());
    assert(read_status == 1);
  }
  std::cout << "Finished adding files to chain!" << std::endl;

  inputChain.SetBranchStatus("*", 0); // so that only the needed branches, explicitly activated below, are read in per event

  eventDetailsStruct eventDetails = eventDetailsStruct(inputChain, options.isMC);
  photonsCollectionStruct photonsCollection = photonsCollectionStruct(inputChain);
  jetsCollectionStruct jetsCollection = jetsCollectionStruct(inputChain, options.isMC);
  MCCollectionStruct MCCollection = MCCollectionStruct(inputChain, options.isMC);

  Long64_t nEvts = inputChain.GetEntries();
  std::cout << "Number of events to process: " << nEvts << std::endl;

  tmProgressBar progressBar = tmProgressBar(static_cast<int>(nEvts));
  int progressBarUpdatePeriod = ((nEvts < 1000) ? 1 : static_cast<int>(0.5 + 1.0*(nEvts/1000)));
  progressBar.initialize();
  for (Long64_t entryIndex = 0; entryIndex < nEvts; ++entryIndex) {
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

    int entryProcessing = static_cast<int>(entryIndex);
    if (entryProcessing > 0 && ((static_cast<int>(entryProcessing) % progressBarUpdatePeriod == 0) || entryProcessing == static_cast<int>(nEvts-1))) progressBar.updateBar(static_cast<double>(1.0*entryProcessing/nEvts), entryProcessing);

    eventExaminationResultsStruct eventExaminationResults = examineEvent(options, parameters, // counters, 
                                                                         entryIndex,
                                                                         year, eventDetails, MCCollection,
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
}

void writeSelectionToFile(optionsStruct &options, TFile *outputFile, const std::vector<eventExaminationResultsStruct>& selectedEventsInfo, selectionRegion& region) {
  std::string regionName = selectionRegionNames[region];
  std::cout << "Beginning to write selected events to file for selection type: " <<  regionName << std::endl;

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

  int nJetsDR; // stores number of jets in event passing deltaR cut
  outputTree->Branch("b_nJets", &nJetsDR, "b_nJets/I");
  float ST; // stores event sT
  outputTree->Branch("b_evtST", &ST, "b_evtST/F");
  float photonPT_leading; // stores PT of leading photon, useful for HLT efficiency
  outputTree->Branch("b_photonPT_leading", &photonPT_leading);
  float photonPT_subLeading; // stores PT of subleading photon, useful for HLT efficiency
  outputTree->Branch("b_photonPT_subLeading", &photonPT_subLeading);
  float photonEta_leading; // stores eta of leading photon, useful for HLT efficiency
  outputTree->Branch("b_photonEta_leading", &photonEta_leading);
  float photonEta_subLeading; // stores eta of subleading photon, useful for HLT efficiency
  outputTree->Branch("b_photonEta_subLeading", &photonEta_subLeading);
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
    photonPT_leading = selectedEventInfo.evt_photonPT_leading;
    photonPT_subLeading = selectedEventInfo.evt_photonPT_subLeading;
    photonEta_leading = selectedEventInfo.evt_photonEta_leading;
    photonEta_subLeading = selectedEventInfo.evt_photonEta_subLeading;
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
  argumentParser.addArgument("inputPathsFile", "", true, "Path to file containing list of input files.");
  argumentParser.addArgument("isMC", "default", true, "Input file is a MC sample -- disable HLT photon trigger and enable additional MC selection.");
  argumentParser.addArgument("disableJetSelection", "default", true, "Do not filter on nJets.");
  argumentParser.addArgument("lineNumberStartInclusive", "", true, "Line number from input file from which to start. The file with this index is included in the processing.");
  argumentParser.addArgument("lineNumberEndInclusive", "", true, "Line number from input file at which to end. The file with this index is included in the processing.");
  argumentParser.addArgument("year", "", true, "Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.");
  argumentParser.setPassedStringValues(argc, argv);

  optionsStruct options = getOptionsFromParser(argumentParser);

  parametersStruct parameters = parametersStruct();
  parameters.tuneParametersForYear(options.year, options.isMC);

  std::stringstream optionsStringstream;
  optionsStringstream << options;
  std::stringstream parametersStringstream;
  parametersStringstream << parameters;

  std::vector<eventExaminationResultsStruct> selectedEventsInfo;

  STRegionsStruct STRegions = STRegionsStruct("STRegionBoundaries.dat");

  statisticsHistograms statistics = statisticsHistograms(options.isMC, HLTEmulation::etaBinEdges, HLTEmulation::pTBinEdges, STRegions.STBoundaries);

  loopOverEvents(options, parameters, options.year, selectedEventsInfo, statistics, STRegions);

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
