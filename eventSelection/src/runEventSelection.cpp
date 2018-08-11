#include "../include/runEventSelection.h"

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
    std::cout << counterTypeString << " photon counters: " << std::endl;
    for (const auto& counterValuePair : counters.photonFailureCounters[counterTypeMapElement.second]) {
      std::cout << photonFailureCategoryNames[counterValuePair.first] << " : " << counterValuePair.second << std::endl;
    }
    std::cout << getNDashes(100) << std::endl;

    std::cout << counterTypeString << " jet counters: " << std::endl;
    for (const auto& counterValuePair : counters.jetFailureCounters[counterTypeMapElement.second]) {
      std::cout << jetFailureCategoryNames[counterValuePair.first] << " : " << counterValuePair.second << std::endl;
    }
    std::cout << getNDashes(100) << std::endl;

    std::cout << counterTypeString << " event counters: " << std::endl;
    for (const auto& counterValuePair : counters.eventFailureCounters[counterTypeMapElement.second]) {
      std::cout << eventFailureCategoryNames[counterValuePair.first] << " : " << counterValuePair.second << std::endl;
    }
    std::cout << getNDashes(100) << std::endl;

  }
  
  std::cout << "Miscellaneous counters: " << std::endl;
  for (const auto& counterValuePair : counters.miscCounters) {
    std::cout << miscCounterNames[counterValuePair.first] << " : " << counterValuePair.second << std::endl;
  }
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

float getRhoCorrectedIsolation(const float& uncorrectedIsolation, const PFTypesForEA& PFType, const float& absEta, const float& rho, const EAValuesStruct& region1EAs, const EAValuesStruct& region2EAs) {
  float effectiveArea;
  if (absEta <= region1EAs.regionUpperBound) effectiveArea = region1EAs.getEffectiveArea(PFType);
  else if (absEta <= region2EAs.regionUpperBound) effectiveArea = region2EAs.getEffectiveArea(PFType);
  else effectiveArea = 0.0; // photons with eta > 1.479 are going to fail the kinematic cut anyway
  float correctedIsolationTest = uncorrectedIsolation - rho*effectiveArea;
  return ((correctedIsolationTest > 0.0) ? correctedIsolationTest : 0.0);
}

photonExaminationResultsStruct examinePhoton(parametersStruct &parameters, countersStruct &counters, const float& rho, const photonsCollectionStruct& photonsCollection, const int& photonIndex) {
  bool passesCommonCuts = true;
  // Kinematic cuts
  applyCondition(counters, photonFailureCategory::eta, passesCommonCuts, (std::fabs((photonsCollection.eta)[photonIndex]) < parameters.photonEtaCut));
  applyCondition(counters, photonFailureCategory::pT, passesCommonCuts, (std::fabs((photonsCollection.pT)[photonIndex]) > parameters.pTCutSubLeading));

  // Electron veto
  applyCondition(counters, photonFailureCategory::conversionSafeElectronVeto, passesCommonCuts, ((photonsCollection.electronVeto)[photonIndex] == TRUETOINTT));

  bool passesLeadingpTCut = ((photonsCollection.pT)[photonIndex] > parameters.pTCutLeading);

  bool passesSelectionAsMedium = passesCommonCuts;
  bool passesSelectionAsFake = passesCommonCuts;

  // Medium ID for medium photon selection
  bool passesMediumID = (((((photonsCollection.ID)[photonIndex])>>1)&1) == 1);
  applyCondition(counters, photonFailureCategory::mediumIDCut, passesSelectionAsMedium, passesMediumID);

  if (passesMediumID) {
    passesSelectionAsFake = false; // don't bother applying the fake selections and set passesFake = False in return struct
  }
  else { // only apply fake selections if the photon is not medium
    applyCondition(counters, photonFailureCategory::hOverE, passesSelectionAsFake, ((photonsCollection.HOverE)[photonIndex] < parameters.towerHOverECut)); // HOverE <-- same as medium selection

    bool passesSigmaIEtaIEtaCut = (parameters.sigmaietaietaRange).isInside((photonsCollection.sigmaIEtaIEta)[photonIndex]); // sigmaietaieta <-- INVERTED from medium selection
    bool passesChargedIsolationCut = (parameters.chargedIsolationRange).isInside(getRhoCorrectedIsolation((photonsCollection.PFChargedIsolationUncorrected)[photonIndex], PFTypesForEA::chargedHadron, std::fabs((photonsCollection.eta)[photonIndex]), rho, parameters.region1EAs, parameters.region2EAs)); // Rho-corrected charged isolation <-- INVERTED from medium selection
    applyCondition(counters, photonFailureCategory::sigmaietaiataORchargedIsolation, passesSelectionAsFake, (passesSigmaIEtaIEtaCut || passesChargedIsolationCut)); // n.b. OR, not XOR

    float pTDependentNeutralIsolationCut = (parameters.neutralIsolationCut).getPolynomialValue(std::fabs((photonsCollection.pT)[photonIndex]));
    applyCondition(counters, photonFailureCategory::neutralIsolation, passesSelectionAsFake, (getRhoCorrectedIsolation((photonsCollection.PFNeutralIsolationUncorrected)[photonIndex], PFTypesForEA::neutralHadron, std::fabs((photonsCollection.eta)[photonIndex]), rho, parameters.region1EAs, parameters.region2EAs) < pTDependentNeutralIsolationCut)); // Neutral isolation <-- same as medium selection

    float pTDependentPhotonIsolationCut = (parameters.photonIsolationCut).getPolynomialValue(std::fabs((photonsCollection.pT)[photonIndex]));
    applyCondition(counters, photonFailureCategory::photonIsolation, passesSelectionAsFake, (getRhoCorrectedIsolation((photonsCollection.PFPhotonIsolationUncorrected)[photonIndex], PFTypesForEA::photon, std::fabs((photonsCollection.eta)[photonIndex]), rho, parameters.region1EAs, parameters.region2EAs) < pTDependentPhotonIsolationCut)); // Neutral isolation <-- same as medium selection
  }

  if (passesSelectionAsFake || passesSelectionAsMedium) incrementCounters(miscCounter::passingPhotons, counters);
  else incrementCounters(miscCounter::failingPhotons, counters);

  photonExaminationResultsStruct photonExaminationResults = photonExaminationResultsStruct(passesSelectionAsMedium, passesSelectionAsFake, passesLeadingpTCut, (photonsCollection.eta)[photonIndex], (photonsCollection.phi)[photonIndex], (photonsCollection.pT)[photonIndex]);
  return photonExaminationResults;
}

bool passesMCSelection(parametersStruct &parameters, const int& nMCParticles, const MCCollectionStruct& MCCollection) {
  int nPhotonsWithNeutralinoMom = 0;
  for (int MCIndex = 0; MCIndex < nMCParticles; ++MCIndex) {
    if ((MCCollection.MCPIDs[MCIndex] == parameters.PIDs.photon) && (MCCollection.MCMomPIDs[MCIndex] == parameters.PIDs.neutralino)) ++nPhotonsWithNeutralinoMom;
  }
  return (nPhotonsWithNeutralinoMom == 2);
}

jetExaminationResultsStruct examineJet(optionsStruct &options, parametersStruct &parameters, countersStruct &counters, const jetsCollectionStruct& jetsCollection, const int& jetIndex) {
  bool passesJetSelection = true;

  //Kinematic cuts: eta, pT
  applyCondition(counters, jetFailureCategory::eta, passesJetSelection, ((jetsCollection.eta)[jetIndex] < parameters.jetEtaCut));
  float jet_pT = (jetsCollection.pT)[jetIndex];
  if (options.isMC) jet_pT += (((jetsCollection.JECUncertainty)[jetIndex])*((jetsCollection.pT)[jetIndex])*options.JECUncertainty);
  applyCondition(counters, jetFailureCategory::pT, passesJetSelection, (jet_pT > parameters.jetpTCut));

  // ID cuts: loose ID, PUID, jetID
  applyCondition(counters, jetFailureCategory::PFLooseID, passesJetSelection, ((jetsCollection.looseID)[jetIndex]));
  applyCondition(counters, jetFailureCategory::puID, passesJetSelection, ((jetsCollection.PUID)[jetIndex] > parameters.jetPUIDThreshold));
  applyCondition(counters, jetFailureCategory::jetID, passesJetSelection, ((jetsCollection.ID)[jetIndex] == 6));

  if (passesJetSelection) incrementCounters(miscCounter::passingJets, counters);
  else incrementCounters(miscCounter::failingJets, counters);

  jetExaminationResultsStruct result = jetExaminationResultsStruct(passesJetSelection, (jetsCollection.eta)[jetIndex], (jetsCollection.phi)[jetIndex], jet_pT);
  return result;
}

bool examineElectron(parametersStruct &parameters, const electronsCollectionStruct& electronsCollection, const int& electronIndex) {
  bool passesElectronSelection = (((electronsCollection.pT)[electronIndex] > parameters.electronPtCut) &&
                                  ((electronsCollection.eta)[electronIndex] < parameters.electronEtaCut) &&
                                  (((((electronsCollection.ID)[electronIndex])>>3)&1) == 1) && // tight electron
                                  ((electronsCollection.dz)[electronIndex] < parameters.electronDzCut) &&
                                  ((electronsCollection.PFPUIsolation)[electronIndex] < parameters.electronPFPUIsolationCut));
  return passesElectronSelection;
}

bool examineMuon(parametersStruct &parameters, const muonsCollectionStruct& muonsCollection, const int& muonIndex) {
  bool passesMuonSelection = (((muonsCollection.pT)[muonIndex] > parameters.muonPtCut) &&
                              ((muonsCollection.PFPUIsolation)[muonIndex] < parameters.muonPFPUIsolationCut) &&
                              (((((muonsCollection.ID)[muonIndex])>>2)&1) == 1)); // tight muon
  return passesMuonSelection;
}

bool examineEvent(optionsStruct &options, parametersStruct &parameters, countersStruct &counters, int& evt_nJetsDR, float& evt_ST, eventDetailsStruct& eventDetails, MCCollectionStruct &MCCollection, photonsCollectionStruct &photonsCollection, jetsCollectionStruct &jetsCollection, electronsCollectionStruct &electronsCollection, muonsCollectionStruct &muonsCollection) {
  evt_nJetsDR = 0;
  evt_ST = 0.0;
  bool passesEventSelection = true;
  if (!(options.isMC) && parameters.HLTPhotonBit >= 0) { // Apply HLT photon selection iff input is not MC and HLTBit is set to a positive integer
    applyCondition(counters, eventFailureCategory::HLTPhoton, passesEventSelection, (((((*(eventDetails.HLTPhotonBits)))>>(parameters.HLTPhotonBit))&1) == 1));
  }

  // Photon selection
  std::vector<angularVariablesStruct> selectedPhotonAnglesList;
  int nPhotonsPassingSubLeadingpTCut = 0;
  int nPhotonsPassingLeadingpTCut = 0;
  int nMediumPhotons = 0;
  int nFakePhotons = 0;
  for (Int_t photonIndex = 0; photonIndex < (*(eventDetails.nPhotons)); ++photonIndex) {
    photonExaminationResultsStruct photonExaminationResults = examinePhoton(parameters, counters, (*(eventDetails.eventRho)), photonsCollection, photonIndex);
    if (photonExaminationResults.passesSelectionAsMedium || photonExaminationResults.passesSelectionAsFake) {
      ++nPhotonsPassingSubLeadingpTCut;
      evt_ST += photonExaminationResults.pT;
      angularVariablesStruct angularVariables = angularVariablesStruct(photonExaminationResults.eta, photonExaminationResults.phi);
      selectedPhotonAnglesList.push_back(angularVariables);
    }
    if (photonExaminationResults.passesLeadingpTCut) ++nPhotonsPassingLeadingpTCut;
    if (photonExaminationResults.passesSelectionAsMedium) ++nMediumPhotons;
    else if (photonExaminationResults.passesSelectionAsFake) ++nFakePhotons;
  }

  applyCondition(counters, eventFailureCategory::wrongNMediumOrFakePhotons, passesEventSelection, ((nMediumPhotons == parameters.nMediumPhotonsRequired) && (nFakePhotons == parameters.nFakePhotonsRequired)));
  applyCondition(counters, eventFailureCategory::wrongNPhotons, passesEventSelection, ((nPhotonsPassingSubLeadingpTCut == 2) && (nPhotonsPassingLeadingpTCut >= 1)));

  // Additional photon selection, only for MC
  if (options.isMC) applyCondition(counters, eventFailureCategory::MCGenInformation, passesEventSelection, (passesMCSelection(parameters, (*(eventDetails.nMCParticles)), MCCollection)));

  // Jet selection
  float evt_HT = 0;
  int nJetsPassingSelection = 0;
  for (Int_t jetIndex = 0; jetIndex < (*(eventDetails.nJets)); ++jetIndex) {
    jetExaminationResultsStruct jetExaminationResults = examineJet(options, parameters, counters, jetsCollection, jetIndex);
    if (jetExaminationResults.passesSelection) {
      ++nJetsPassingSelection;
      float min_dR = -1.0;
      for (const auto& selectedPhotonAngles : selectedPhotonAnglesList) {
        float dR = std::sqrt(std::pow((selectedPhotonAngles.eta - jetExaminationResults.eta), 2.0f) + std::pow((selectedPhotonAngles.phi - jetExaminationResults.phi), 2.0f));
        if ((min_dR < 0) || (dR < min_dR)) min_dR = dR;
      }
      evt_HT += jetExaminationResults.pT; // Add hT whether or not jet passes deltaR check
      applyCondition(counters, jetFailureCategory::deltaR, jetExaminationResults.passesSelection, ((min_dR > parameters.minDeltaRCut) || (min_dR < 0.0)));
    }
    if (jetExaminationResults.passesSelection) {// Could have changed while applying the delta-R condition, so need to check again
      evt_ST += jetExaminationResults.pT; // Add sT only if jet passes deltaR check, to avoid double-counting
      ++evt_nJetsDR; // Count only those jets that are sufficiently away from a photon
    }
  }
  applyCondition(counters, eventFailureCategory::wrongNJets, passesEventSelection, (nJetsPassingSelection >= 2));
  applyCondition(counters, eventFailureCategory::hTCut, passesEventSelection, (evt_HT >= parameters.HTCut));

  // Electron veto
  int nTightElectrons = 0;
  for (Int_t electronIndex = 0; electronIndex < (*(eventDetails.nElectrons)); ++electronIndex) {
    bool passesElectronSelection = examineElectron(parameters, electronsCollection, electronIndex);
    if (passesElectronSelection) ++nTightElectrons;
  }
  applyCondition(counters, eventFailureCategory::electronVeto, passesEventSelection, (nTightElectrons == 0));

  // Muon veto
  int nTightMuons = 0;
  for (Int_t muonIndex = 0; muonIndex < (*(eventDetails.nMuons)); ++muonIndex) {
    bool passesMuonSelection = examineMuon(parameters, muonsCollection, muonIndex);
    if (passesMuonSelection) ++nTightMuons;
  }
  applyCondition(counters, eventFailureCategory::muonVeto, passesEventSelection, (nTightMuons == 0));

  // Add MET to ST only if it clears threshold
  if ((*(eventDetails.PFMET)) > parameters.METThreshold) evt_ST += (*(eventDetails.PFMET));

  return passesEventSelection;
}

void runSelection(optionsStruct &options, parametersStruct &parameters, countersStruct &counters, TFile *outputFile) {
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
      std::cout << "Adding... " << inputFileName << std::endl;
      inputChain.Add(inputFileName.c_str()); // Add files to TChain
    }
  }
  std::cout << "Finished adding files to chain!" << std::endl;

  // Initialize tree in output file as empty clone
  TDirectory *outDir = outputFile->mkdir("ggNtuplizer");
  outDir->cd();
  TTree *outputTree = inputChain.CloneTree(0);

  // Initialize output branches
  int evt_nJetsDR; // stores number of jets in event passing deltaR cut
  float evt_ST; // stores event sT
  outputTree->Branch("b_nJets", &evt_nJetsDR, "b_nJets/I");
  outputTree->Branch("b_evtST", &evt_ST, "b_evtST/F");

  Long64_t nEvts = inputChain.GetEntries();
  Long64_t nEntriesToProcess = 1 + options.counterEndInclusive - options.counterStartInclusive;
  std::cout << "Number of available events: " << nEvts << std::endl;

  // inputChain.SetBranchStatus("*", 0); // so that only the needed branches, explicitly activated below, are read in per event

  TTreeReader chainReader(&inputChain);

  eventDetailsStruct eventDetails = eventDetailsStruct(chainReader, options.isMC);
  photonsCollectionStruct photonsCollection = photonsCollectionStruct(chainReader);
  jetsCollectionStruct jetsCollection = jetsCollectionStruct(chainReader, options.isMC);
  electronsCollectionStruct electronsCollection = electronsCollectionStruct(chainReader);
  muonsCollectionStruct muonsCollection = muonsCollectionStruct(chainReader);
  MCCollectionStruct MCCollection = MCCollectionStruct();
  if (options.isMC) MCCollection = MCCollectionStruct(chainReader);

  tmProgressBar progressBar = tmProgressBar(static_cast<int>(nEntriesToProcess));
  int progressBarUpdatePeriod = static_cast<int>(nEntriesToProcess) < 1000 ? 1 : static_cast<int>(0.5 + 1.0*nEntriesToProcess/1000);
  Long64_t currentEntryIndex = -1l + options.counterStartInclusive; // so that first call to Next changes it to options.counterStartInclusive
  progressBar.initialize();

  // for (Long64_t entryIndex = options.counterStartInclusive; entryIndex <= options.counterEndInclusive; ++entryIndex) {
    // if (entryIndex > nEvts) {
    //   std::cout << "ERROR: Entry index falls outside event range. Index: " << entryIndex << ", available number of events: " << nEvts << std::endl;
    //   std::exit(EXIT_FAILURE);
    // }
    // Long64_t loadStatus = inputChain.LoadTree(entryIndex);
    // if (loadStatus < 0) {
    //   std::cout << "ERROR in loading tree for entry index: " << entryIndex << "; load status = " << loadStatus << std::endl;
    //   std::exit(EXIT_FAILURE);
    // }
    // int nBytesRead = inputChain.GetEntry(entryIndex, 0); // Get only the required branches
    // if (nBytesRead <= 0) {
    //   std::cout << "ERROR: Failed to read any information from entry at index: " << entryIndex << "; nBytesRead = " << nBytesRead << std::endl;
    //   std::exit(EXIT_FAILURE);
    // }

    // int entryProcessing = static_cast<int>(entryIndex - options.counterStartInclusive);

  chainReader.SetEntriesRange(options.counterStartInclusive, 1l + options.counterEndInclusive);
  while (chainReader.Next()) {
    ++currentEntryIndex;
    Long64_t TReaderEntryIndex = chainReader.GetCurrentEntry();
    if (!(TReaderEntryIndex == currentEntryIndex)) {
      std::cout << "ERROR: TReaderEntryIndex = " << TReaderEntryIndex << " is not the same as currentEntryIndex = " << currentEntryIndex << std::endl;
      std::exit(EXIT_FAILURE);
    }
    int entryProcessing = static_cast<int>(currentEntryIndex - options.counterStartInclusive);
    if (entryProcessing > 0 && ((static_cast<int>(entryProcessing) % progressBarUpdatePeriod == 0) || entryProcessing == static_cast<int>(nEntriesToProcess-1))) progressBar.updateBar(static_cast<double>(1.0*entryProcessing/nEntriesToProcess), entryProcessing);
    
    evt_nJetsDR = 0;
    evt_ST = 0.0;
    bool passesEventSelection = examineEvent(options, parameters, counters, evt_nJetsDR, evt_ST, eventDetails, MCCollection, photonsCollection, jetsCollection, electronsCollection, muonsCollection);
    if (!(passesEventSelection)) {
      incrementCounters(miscCounter::failingEvents, counters);
      continue;
    }
    nBytesRead = inputChain.GetEntry(currentEntryIndex, 1); // Get all branches before filling output tree
    if (nBytesRead <= 0) {
      std::cout << "ERROR: For selected event, failed to read any information from entry at index: " << currentEntryIndex << "; nBytesRead = " << nBytesRead << std::endl;
      std::exit(EXIT_FAILURE);
    }
    outputTree->Fill();
    incrementCounters(miscCounter::acceptedEvents, counters);
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
  argumentParser.addArgument("photonSelectionType", "fake", true, "Photon selection type: can be any one of: \"fake\", \"medium\", \"mediumfake\", \"fakeMC\", \"mediumMC\", \"mediumfakeMC\"");
  argumentParser.addArgument("year", "2017", false, "Year of data-taking. Affects the HLT photon Bit index in the format of the n-tuplizer on which to trigger (unless sample is MC), and the photon ID cuts which are based on year-dependent recommendations.");
  argumentParser.addArgument("JECUncertainty", "0", false, "Apply a uniform upward or downward jet energy uncertainty correction to jet pt. Default: 0, i.e. do not apply any other correction. +/-1 are allowed as well, shifting all jet pt up or down respectively by 1.0 times the uncertainty on the jet energy correction.");
  argumentParser.setPassedStringValues(argc, argv);

  optionsStruct options = getOptionsFromParser(argumentParser);

  parametersStruct parameters = parametersStruct();
  parameters.tuneParametersForYear(options.year);
  parameters.tuneParametersForPhotonSelectionType(options.photonSelectionType);

  countersStruct counters = countersStruct();
  initializeCounters(counters);

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

  runSelection(options, parameters, counters, outputFile);

  outputFile->WriteTObject(parametersObject);
  outputFile->WriteTObject(optionsObject);
  outputFile->Close();

  std::cout << getNDashes(100) << std::endl;
  printCounters(counters);

  // std::cout << getNDashes(100) << std::endl;
  // // For testing: first print counters, should have all zeros
  // std::cout << "Empty counters:" << std::endl;
  // printCounters(counters);
  // incrementCounters(photonFailureCategory::pT, counterType::differential, counters);
  // incrementCounters(photonFailureCategory::pT, counterType::global, counters);
  // incrementCounters(jetFailureCategory::eta, counterType::global, counters);
  // incrementCounters(eventFailureCategory::MCGenInformation, counterType::differential, counters);
  // incrementCounters(eventFailureCategory::MCGenInformation, counterType::differential, counters);
  // incrementCounters(eventFailureCategory::MCGenInformation, counterType::differential, counters);
  // incrementCounters(miscCounter::acceptedEvents, counters);
  // // Should now have 1 in these categories
  // printCounters(counters);

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
